!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Provides functionality that performs coupling of a scalar quantity

module coupler_exchange_0d_mod

#ifdef MCT
  use mod_oasis,                only: oasis_put, oasis_put_inquire,    &
                                      oasis_sent, oasis_sentout
#endif
  use constants_mod,            only: i_def, r_def, str_def, imdi
  use field_mod,                only: field_type
  use field_collection_mod,     only: field_collection_type
  use log_mod,                  only: log_event,       &
                                      LOG_LEVEL_DEBUG, &
                                      LOG_LEVEL_INFO,  &
                                      LOG_LEVEL_ERROR, &
                                      log_scratch_space
  use model_clock_mod,          only: model_clock_type
  use mpi_mod,                  only: global_mpi

  implicit none

  private

  public coupler_send_0d

  contains

  !>@brief Sends a 0d variable (i.e. a scalar) to another component
  !> @param [in] scalar      The value to be passed to the other component
  !> @param [in] name        Name of the field used to generate the scalar
  !> @param [in] var_id      Oasis id for scalar being sent
  !> @param [in] model_clock Time within the model.
  !>
  subroutine coupler_send_0d( scalar, name, var_id, model_clock )

  implicit none

  real(r_def),             intent(in) :: scalar
  character(str_def),      intent(in) :: name
  integer(i_def),          intent(in) :: var_id
  class(model_clock_type), intent(in) :: model_clock

#ifdef MCT
  ! Oasis requires that scalars are sent as a one-element array
  real(r_def)                         :: tiny_array(1)
  ! Returned Oasis error code
  integer(i_def)                      :: ierror
  ! Rank number of current PE
  integer(i_def)                      :: local_rank
  ! Time of coupling in seconds from start of the run
  integer(i_def)                      :: coupling_time

  coupling_time = &
       int( model_clock%seconds_from_steps(model_clock%get_step())          &
            - model_clock%seconds_from_steps(model_clock%get_first_step()), &
            i_def)

  if (var_id /= imdi) then
    ! Send the data to the coupler
    call oasis_put_inquire(var_id, coupling_time, ierror)
    if (ierror == oasis_sent .or. ierror == oasis_sentout) then
      ! Oasis_put expects an array
      tiny_array(:)=scalar
      ! 0D coupling - so only pass from pe0
      local_rank  = global_mpi%get_comm_rank()
      if(local_rank == 0) then
         call oasis_put(var_id, coupling_time , tiny_array, ierror)

         write(log_scratch_space, '(3A, 2E12.3)' ) &
                            "coupler_send_0d: field ", &
                            trim(name), &
                            " sent with value = ", scalar
         call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      endif
    else
      write(log_scratch_space, '(3A)' ) "coupler_send_0d: field ", &
                     trim(name), " NOT exchanged on this timestep"
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
    endif
  else
    write(log_scratch_space, '(3A)' ) "PROBLEM coupler_send_0d: field ", &
                                      trim(name), " cpl_id NOT set"
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

#else
  write(log_scratch_space, '(A)' ) &
                 "coupler_send_0d: to use OASIS cpp directive MCT must be set"
  call log_event( log_scratch_space, LOG_LEVEL_ERROR )

#endif

  end subroutine coupler_send_0d

  ! Scalars are not currently received by LFRic,
  ! so "coupler_receive_0d" has not been implemented

end module coupler_exchange_0d_mod
