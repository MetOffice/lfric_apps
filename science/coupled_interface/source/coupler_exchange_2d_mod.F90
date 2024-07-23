!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Provides functionality that performs coupling sends and receives
!> @details  Extends the abstract external field class to provide the
!> "copy_from_lfric" and "copy_to_lfric" functions (renamed in this module to
!> cpl_field_send and cpl_field_receive) that perform coupling send and receive
!> operations on two-dimensional horizontal field data

module coupler_exchange_2d_mod

#ifdef MCT
  use mod_oasis,                only: oasis_get,                       &
                                      oasis_get_ncpl, oasis_get_freqs, &
                                      oasis_put, oasis_put_inquire,    &
                                      oasis_sent, oasis_sentout,       &
                                      oasis_recvd, oasis_recvout, oasis_out
#endif
  use abstract_external_field_mod, only: abstract_external_field_type
  use constants_mod,            only: i_def, r_def, l_def, str_def, imdi
  use field_mod,                only: field_type, field_proxy_type
  use field_collection_mod,     only: field_collection_type
  use function_space_collection_mod,  &
                                only: function_space_collection
  use function_space_mod,       only: function_space_type
  use log_mod,                  only: log_event,       &
                                      LOG_LEVEL_DEBUG, &
                                      LOG_LEVEL_INFO,  &
                                      LOG_LEVEL_ERROR, &
                                      log_scratch_space
  use model_clock_mod,          only: model_clock_type
  use mpi_mod,                  only: global_mpi

  implicit none

  private

type, extends(abstract_external_field_type), public :: coupler_exchange_2d_type
  private
  ! Time of coupling in seconds from start of the run
  integer(i_def) :: coupling_time
  !> Size  of the coupling send and receive fields
  integer(i_def) :: coupling_size
  !> Index used to sort data passed through the coupler
  integer(i_def), allocatable :: sorting_index(:)
contains
  private
  !> Initialises the object
  procedure, public :: initialise
  !> Copy data from the LFRic field and pass to the coupler
  procedure, public :: copy_from_lfric => coupler_send_2d
  !> Copy data from the coupler into the LFRic field
  procedure, public :: copy_to_lfric => coupler_receive_2d
  !> Sets the time ready for coupling
  procedure, public :: set_time
  !> Checks if the currently set time is scheduled for a coupling operation
  procedure, public :: is_coupling_time
  !> Tidies up on destruction
  final             :: finalise
end type  coupler_exchange_2d_type

  contains

  !> @brief Initialises the external field used for coupling
  !> @param [in] lfric_field_ptr Pointer to an lfric field
  !> @param [in] sorting_index   Index to sort data for coupling
  !
  subroutine initialise( self, lfric_field_ptr, sorting_index )
  implicit none

  class(coupler_exchange_2d_type), intent(inout) :: self
  type(field_type), pointer,       intent(in)    :: lfric_field_ptr
  integer(i_def),                  intent(in)    :: sorting_index(:)


  call self%abstract_external_field_initialiser(lfric_field_ptr)

  self%coupling_time = 0
  self%coupling_size = size(sorting_index)
  allocate(self%sorting_index, source=sorting_index)

  end subroutine initialise


  !>@brief Sends field data through the coupler to another model component
  !>@param return_code Optional return code from the copy_from procedure
  !
  subroutine coupler_send_2d(self, return_code)

  implicit none

  class(coupler_exchange_2d_type), intent(inout) :: self
  integer(i_def), optional,        intent(out)   :: return_code

#ifdef MCT
  ! Field from which the data will be sent - maybe a multi-data field
  type(field_type), pointer                 :: field
  ! Proxy of the multidata field
  type( field_proxy_type )                  :: field_proxy
  ! Number of multi-data levels in field
  integer(i_def)                            :: ndata
  ! Processed and sorted data ready to be passed to Oasis
  real(r_def)                               :: sorted_data(self%coupling_size)
  ! Oasis id for variable or data level being sent
  integer(i_def)                            :: var_id
  ! Name of the field being sent
  character(str_def)                        :: name
  ! OASIS error code
  integer(i_def)                            :: kinfo
  ! Error code to return from all receives
  integer(i_def)                            :: ierror
  ! Looping index over the field data
  integer(i_def)                            :: i
  ! Looping index over the multi-data levels
  integer(i_def)                            :: nmulti

  field       => self%get_lfric_field_ptr()
  name        =  trim(adjustl(field%get_name()))
  field_proxy =  field%get_proxy()
  ndata       =  field_proxy%vspace%get_ndata()

  ierror = 0

  do nmulti = 1, ndata
    var_id = field%get_cpl_id(nmulti)
    if (var_id /= imdi) then
      ! Reorder the LFRic data into the outgoing Oasis send buffer
      ! Get just the single category we are sending
      do i = 1, self%coupling_size
        sorted_data(i) = &
            field_proxy%data((self%sorting_index(i)-1)*ndata+nmulti)
      enddo

      ! Send the data to the coupler
      call oasis_put(var_id, self%coupling_time, sorted_data(:), kinfo)

      write(log_scratch_space, '(3A, 2E12.3)' ) "cpl_field_send: field ", &
                  trim(name)," sent with min,max = ",                     &
                  minval(sorted_data), maxval(sorted_data)
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

    else
      ierror = 1
      write(log_scratch_space, '(3A)' ) "PROBLEM cpl_field_send: field ", &
                                        trim(name), " cpl_id NOT set"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    endif
  enddo

  if( present(return_code) ) return_code = ierror

#else
  if( present(return_code) ) return_code = 1
  write(log_scratch_space, '(A)' ) &
                 "cpl_field_send: to use OASIS cpp directive MCT must be set"
  call log_event( log_scratch_space, LOG_LEVEL_ERROR )

#endif

  end subroutine coupler_send_2d


  !> @brief Receives field data through the coupler from another model component
  !> @param return_code The return code from the copy_to procedure
  !
  subroutine coupler_receive_2d( self, return_code )
  implicit none

  class(coupler_exchange_2d_type), intent(inout) :: self
  integer(i_def), optional,        intent(out)   :: return_code

#ifdef MCT
  ! Field into which the data will be received - maybe a multi-data field
  type(field_type), pointer                 :: field
  ! Proxy of the receiving field
  type(field_proxy_type)                    :: field_proxy
  ! Number of multi-data levels in field
  integer(i_def)                            :: ndata
  ! Data received from OASIS
  real(r_def)                               :: sorted_data(self%coupling_size)
  ! Oasis id for varialble or data level
  integer(i_def)                            :: var_id
  ! Name of the variable being received
  character(str_def)                        :: name
  ! OASIS error code
  integer(i_def)                            :: kinfo
  ! Error code to return from all receives
  integer(i_def)                            :: ierror
  ! Looping index over the field data
  integer(i_def)                            :: i
  ! Looping index over the multi-data levels
  integer(i_def)                            :: nmulti


  field       => self%get_lfric_field_ptr()
  name        =  trim(adjustl(field%get_name()))
  field_proxy =  field%get_proxy()
  ndata       =  field_proxy%vspace%get_ndata()

  ierror = 0

  do nmulti = 1, ndata
    var_id = field%get_cpl_id(nmulti)
    if (var_id /= imdi) then

      call oasis_get(var_id, self%coupling_time, sorted_data(:), kinfo)

      if (kinfo == oasis_recvd .or. kinfo == oasis_recvout) then
        do i = 1, self%coupling_size
          field_proxy%data((self%sorting_index(i)-1)*ndata+nmulti) = &
                                                                 sorted_data(i)
        enddo
        write(log_scratch_space, '(3A)' ) "cpl_field_receive: field ", &
                           trim(name), " received"
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      else
        ierror = 1
        write(log_scratch_space, '(3A)' ) "cpl_field_receive: field ", &
                           trim(name), " NOT exchanged on this timestep"
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      endif
    else
      write(log_scratch_space, '(3A)' ) "PROBLEM cpl_field_receive: field ", &
                                         trim(name), " cpl_id NOT set"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    endif
  enddo

  call field_proxy%set_dirty()

  if( present(return_code) ) return_code = ierror

#else
  if( present(return_code) ) return_code = 1
  write(log_scratch_space, '(A)' ) &
               "cpl_field_receive: to use OASIS cpp directive MCT must be set"
  call log_event( log_scratch_space, LOG_LEVEL_ERROR )

#endif
  end subroutine coupler_receive_2d


  !> @brief Sets the time ready for coupling
  !> @param [in] model_clock  Time within the model.
  !
  subroutine set_time(self, model_clock)
  implicit none
  class(coupler_exchange_2d_type), intent(inout) :: self
  class(model_clock_type),         intent(in)    :: model_clock

  ! Store time of coupling in seconds from start of the run
  self%coupling_time = &
       int( model_clock%seconds_from_steps(model_clock%get_step())          &
            - model_clock%seconds_from_steps(model_clock%get_first_step()), &
            i_def)

  end subroutine set_time

  !> @brief   Checks if the currently set time is scheduled for a coupling
  !>          operation
  !> @return  Whether the currently set time is scheduled for coupling
  !
  function is_coupling_time(self) result(is_coupling)
  implicit none
  class(coupler_exchange_2d_type), intent(inout) :: self

  logical(l_def)            :: is_coupling
  ! Field from which the data will be sent
  type(field_type), pointer :: field
  ! Oasis id for variable or data level being sent
  integer(i_def)            :: var_id
  ! Oasis error code
  integer(i_def)            :: kinfo

  field => self%get_lfric_field_ptr()
  var_id = field%get_cpl_id(1)

  is_coupling = .false.
#ifdef MCT
  call oasis_put_inquire(var_id, self%coupling_time, kinfo)
  if (kinfo == oasis_sent .or. kinfo == oasis_sentout) then
    is_coupling = .true.
  end if
#endif

  end function is_coupling_time

  ! Finaliser
  !
  subroutine finalise(self)
  implicit none
  type(coupler_exchange_2d_type), intent(inout) :: self

  if (allocated(self%sorting_index)) deallocate(self%sorting_index)

  end subroutine finalise

end module coupler_exchange_2d_mod
