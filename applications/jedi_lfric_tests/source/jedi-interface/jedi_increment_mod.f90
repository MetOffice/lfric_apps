!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the JEDI Increment emulator class.
!>
!> @details This module holds a JEDI Increment emulator class that includes
!>          only the functionality required by the LFRic-JEDI model interface
!>          on JEDI-LFRIC. This includes i) read/write to a defined set of
!>          fields and ii) ability to interoperate between LFRic fields and
!>          JEDI fields (Atlas fields here).
module jedi_increment_mod

  use, intrinsic :: iso_fortran_env, only : real64
  use atlas_field_emulator_mod,      only : atlas_field_emulator_type
  use atlas_field_interface_mod,     only : atlas_field_interface_type
  use io_context_mod,                only : io_context_type
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use jedi_increment_config_mod,     only : jedi_increment_config_type
  use jedi_lfric_field_meta_mod,     only : jedi_lfric_field_meta_type
  use field_collection_mod,          only : field_collection_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_INFO,     &
                                            LOG_LEVEL_ERROR
  use constants_mod,                 only : i_def, l_def, str_def
  use model_clock_mod,               only : model_clock_type

  implicit none

  private

type, public :: jedi_increment_type
  private

  !> These fields emulate a set of Atlas fields are and used purely for testing
  type ( atlas_field_emulator_type ), allocatable :: fields(:)

  !> An object that stores the field meta-data associated with the fields
  type( jedi_lfric_field_meta_type )              :: field_meta_data

  !> Field collection to perform IO
  type( field_collection_type ), public           :: io_collection

  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formatting for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'
  type( jedi_datetime_type )                      :: inc_time

  !> The jedi_geometry object
  type( jedi_geometry_type ), pointer, public     :: geometry => null()

contains

  !> Jedi increment initialiser.
  procedure :: initialise => increment_initialiser_zero
  procedure :: increment_initialiser

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the fields in a LFRic field collection
  procedure, private :: setup_interface_to_field_collection

  !> Set the internal Atlas field emulators from a input field_collection
  procedure, public :: set_from_field_collection

  !> Get via the internal Atlas field emulators via a copy to a field_collection
  procedure, public :: get_to_field_collection

  !> Return the curent time
  procedure, public :: valid_time

  !> Read model fields from file into fields
  procedure, public :: read_file

  !> @todo Write model fields to file from the fields
  !> procedure, public :: write_file

  !> Zero the model fields
  procedure, public :: zero

  !> Update the curent time
  procedure, public :: update_time

  !> Print field
  procedure, public :: print_field

  !> Set the internal Atlas field emulators from the model prognostics fields
  procedure, public :: set_from_model_prognostics

  !> Get the model prognostics fields from the internal Atlas field emulators
  procedure, public :: get_to_model_prognostics

  !> Finalizer
  final             :: jedi_increment_destructor

end type jedi_increment_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! Methods required by the JEDI (OOPS) model interface
!-------------------------------------------------------------------------------

!> @brief    Initialiser via read for jedi_increment_type
!>
!> @param [in] geometry The geometry object required to construct the increment
!> @param [in] config   The configuration object including the required
!>                      information to construct a increment and read a file to
!>                      initialise the fields
subroutine increment_initialiser_zero( self, geometry, config )

  implicit none

  class( jedi_increment_type ),       intent(inout) :: self
  type( jedi_geometry_type ), target,    intent(in) :: geometry
  type( jedi_increment_config_type ), intent(inout) :: config

  ! Create
  call self%increment_initialiser( geometry, config )
  ! init fields
  call self%zero()

end subroutine increment_initialiser_zero

!> @brief    Initialiser for jedi_increment_type
!>
!> @param [in] geometry The geometry object required to construct the increment
!> @param [in] config   A configuration object including the required
!>                      information to construct a increment
subroutine increment_initialiser( self, geometry, config )

  use fs_continuity_mod,     only : W3, Wtheta

  implicit none

  class( jedi_increment_type ),       intent(inout) :: self
  type( jedi_geometry_type ), target,    intent(in) :: geometry
  type( jedi_increment_config_type ), intent(inout) :: config

  ! Local
  integer(i_def) :: n_horizontal
  integer(i_def) :: n_levels
  integer(i_def) :: n_layers
  integer(i_def) :: ivar
  integer(i_def) :: n_variables
  integer(i_def) :: fs_id
  logical(l_def) :: twod_field

  ! Setup
  self%field_meta_data = config%field_meta_data
  self%inc_time = config%inc_time
  self%geometry => geometry
  n_variables = self%field_meta_data%get_n_variables()

  allocate( self%fields(n_variables) )

  n_horizontal=geometry%get_n_horizontal()
  n_layers=geometry%get_n_layers()

  do ivar=1,n_variables

    fs_id      = self%field_meta_data%get_variable_function_space(ivar)
    twod_field = self%field_meta_data%get_variable_is_2d(ivar)

    select case (fs_id)
    case (W3)
      if (twod_field) then
        n_levels = 1
      else
        n_levels = n_layers
      end if

    case (Wtheta)
      n_levels = n_layers + 1

    case default
      log_scratch_space = 'The requested LFRic function space is not supported.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

    call self%fields(ivar)%initialise(        &
                                n_levels,     &
                                n_horizontal, &
                                self%field_meta_data%get_variable_name(ivar) )
  end do

  ! Setup the io_collection (empty ready for read/write operations)
  call self%io_collection%initialise(name = 'io_collection', table_len=100)

end subroutine increment_initialiser

!> @brief    Returns the current time of the increment
!>
!> @return   time   The current time of the state
function valid_time( self ) result(time)

  implicit none

  class( jedi_increment_type ), intent(inout) :: self
  type( jedi_datetime_type )                  :: time

  time = self%inc_time

end function valid_time

!> @brief    A method to update the internal Atlas field emulators
!>
!> @param [in] read_time   The data datetime to be read
!> @param [in] file_prefix Character array that specifies the file to read from
subroutine read_file( self, read_time, file_prefix )

  use jedi_lfric_io_update_mod,      only : update_io_field_collection
  use lfric_xios_read_mod,           only : read_state

  implicit none

  class( jedi_increment_type ), intent(inout) :: self
  type( jedi_datetime_type ),      intent(in) :: read_time
  character(len=*),                intent(in) :: file_prefix

  ! Local
  character( len=str_def ), allocatable :: variable_names(:)
  class( io_context_type ),     pointer :: context_ptr

  ! Ensure the JEDI-IO context is set to current
  context_ptr => self%geometry%get_io_context()
  call context_ptr%set_current()

  ! Set the clock to the desired read time
  call set_clock(self, read_time)

  ! Ensure the io_collection contains the variables defined in the list
  ! stored in the type
  call update_io_field_collection( self%io_collection,            &
                                   self%geometry%get_mesh(),      &
                                   self%geometry%get_twod_mesh(), &
                                   self%field_meta_data )

  ! Read the increment into the io_collection
  call read_state( self%io_collection, prefix=file_prefix )

  ! Get the list of internal variables
  call self%field_meta_data%get_variable_names( variable_names )

  ! Set from the io_collection populated by the read to the Atlas field
  ! emulators
  call self%set_from_field_collection( variable_names, self%io_collection )

end subroutine read_file

!> @brief    A method to set all internal Atlas field emulators to zero
!>
subroutine zero( self )

  implicit none

  class( jedi_increment_type ),    intent(inout) :: self

  ! Local
  integer :: n_variables
  integer :: ivar

  ! Zero the Atlas fields
  n_variables = self%field_meta_data%get_n_variables()
  do ivar=1,n_variables
    call self%fields(ivar)%zero()
  enddo

end subroutine zero

!------------------------------------------------------------------------------
! Local methods to support LFRic-JEDI implementation
!------------------------------------------------------------------------------

!> @brief    Set the Atlas field emulators from the model prognostics stored in
!>           the model prognostics field_collection
!>
!> @param [inout] model_prognostics The model prognostic field collection
!>
subroutine set_from_model_prognostics( self, model_prognostics )

  use transform_winds_mod, only : wind_vector_to_scalar

  implicit none

  class( jedi_increment_type ),  intent(inout) :: self
  type( field_collection_type ), intent(inout) :: model_prognostics

  ! Local
  character( len=str_def ), allocatable :: variable_names(:)

  !> @todo This is a placeholder for code to come
  !>       as part of #3734
  !>
  !>       In this code:
  !>       The cell-centered winds (W3/Wtheta) stored in the input prognostics
  !>       field collection are updated by interpolation from the W2 wind
  !>       field. This is performed in wind_vector_to_scalar via
  !>       map_physics_winds and split_wind_alg. The updated cell centred winds
  !>       are then updated via the copy: set_from_field_collection

  call wind_vector_to_scalar( model_prognostics )

  ! Get a list of variables to copy
  call self%field_meta_data%get_variable_names( variable_names )

  ! Set model_prognostics to the Atlas field emulators
  call self%set_from_field_collection( variable_names, model_prognostics )

end subroutine set_from_model_prognostics

!> @brief   Get the model prognostics fields from the internal Atlas field
!>          emulators
!>
!> @param [inout] model_prognostics The model prognostic field collection
!>
subroutine get_to_model_prognostics( self, model_prognostics )

  use transform_winds_mod, only : wind_scalar_to_vector

  implicit none

  class( jedi_increment_type ),  intent(inout) :: self
  type( field_collection_type ), intent(inout) :: model_prognostics

  ! Local
  character( len=str_def ), allocatable :: variable_names(:)

  ! Get a list of variables to copy
  call self%field_meta_data%get_variable_names( variable_names )

  ! Get Atlas field emulators to the model_prognostics
  call self%get_to_field_collection( variable_names, model_prognostics )

  !> @todo This is a placeholder for code to come
  !>       as part of #3734
  !>
  !>       In this code:
  !>       The cell-centered winds stored in Atlas have been copied via
  !>       get_to_field_collection to the cell centered (W3/Wtheta) model
  !>       prognostic field collection. The following performs the
  !>       interpolation of the winds from the cell centered (W3/Wtheta) to W2
  !>       using the set_wind method.
  !>
  call wind_scalar_to_vector( model_prognostics )

end subroutine get_to_model_prognostics

!> @brief    Setup atlas_lfric_interface_fields that enables copying
!>           between Atlas field emulators and the LFRic fields in
!>           io_collection
!>
!> @param [inout] atlas_lfric_interface_fields The atlas lfric interface bundle
!> @param [in]    variable_names               The list of variables to setup
!> @param [inout] field_collection             The fields to link to
!>
subroutine setup_interface_to_field_collection( self, &
                                                atlas_lfric_interface_fields, &
                                                variable_names, &
                                                field_collection )

  use jedi_lfric_utils_mod, only : get_model_field
  use field_mod,            only : field_type

  implicit none

  class( jedi_increment_type ),       intent(inout) :: self
  type( atlas_field_interface_type ), intent(inout) :: atlas_lfric_interface_fields(:)
  character( len=str_def ),           intent(in)    :: variable_names(:)
  type( field_collection_type ),      intent(inout) :: field_collection

  ! Local
  integer(i_def)            :: ivar
  type(field_type), pointer :: lfric_field_ptr
  real(real64),     pointer :: atlas_data_ptr(:,:)
  integer(i_def),   pointer :: horizontal_map_ptr(:)
  integer(i_def)            :: n_variables
  logical                   :: all_variables_exists

  ! Check that the increment conatins all the required fields
  all_variables_exists = self%field_meta_data%check_variables_exist(variable_names)
  if (.not. all_variables_exists) then
    log_scratch_space = 'The Atlas increment does not conatin all the required fields.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  ! Link the Atlas emulator fields with lfric fields
  call self%geometry%get_horizontal_map(horizontal_map_ptr)
  n_variables=size(variable_names)
  do ivar = 1, n_variables
    ! Get the required data
    !! field_meta_data and field_collection
    call get_model_field( variable_names(ivar), &
                          field_collection, lfric_field_ptr )
    atlas_data_ptr => self%fields(ivar)%get_data()
    call atlas_lfric_interface_fields(ivar)%initialise( atlas_data_ptr,     &
                                                        horizontal_map_ptr, &
                                                        lfric_field_ptr )
  end do

end subroutine setup_interface_to_field_collection


!> @brief    Set the internal Atlas field emulators by copying the fields
!>           stored in a field_collection
!>
!> @param [in]    variable_names    The list of variables to copy
!> @param [inout] field_collection  The fields to copy from
!>
subroutine set_from_field_collection( self, variable_names, field_collection )

  implicit none

  class( jedi_increment_type ),  intent(inout) :: self
  character( len=str_def ),      intent(in)    :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer(i_def) :: ivar
  integer(i_def) :: n_variables
  type( atlas_field_interface_type ), allocatable :: atlas_lfric_interface_fields(:)

  ! Allocate space for the interface fields
  n_variables = size(variable_names)
  allocate( atlas_lfric_interface_fields(n_variables) )

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_field_collection( atlas_lfric_interface_fields, &
                                                 variable_names, &
                                                 field_collection )

  ! Copy the LFRic fields in the field_collection to the Atlas field emulators
  do ivar = 1, n_variables
    call atlas_lfric_interface_fields(ivar)%copy_from_lfric()
  end do

end subroutine set_from_field_collection

!> @brief    Get a copy of the internal Atlas field emulators by copying into a
!>           field_collection
!>
!> @param [in]    variable_names    The list of variables to copy
!> @param [inout] field_collection  The fields to copy to
!>
subroutine get_to_field_collection( self, variable_names, field_collection )

  implicit none

  class( jedi_increment_type ),  intent(inout) :: self
  character( len=str_def ),         intent(in) :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer(i_def) :: ivar
  integer(i_def) :: n_variables
  type( atlas_field_interface_type ), allocatable :: atlas_lfric_interface_fields(:)

  ! Allocate space for the interface fields
  n_variables = size(variable_names)
  allocate( atlas_lfric_interface_fields(n_variables) )

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_field_collection( atlas_lfric_interface_fields, &
                                                 variable_names, &
                                                 field_collection )

  ! Copy the LFRic fields in the field_collection to the Atlas field emulators
  do ivar = 1, n_variables
    call atlas_lfric_interface_fields(ivar)%copy_to_lfric()
  end do

end subroutine get_to_field_collection

!> @brief    Update the inc_time by a single time-step
!>
!> @param [in] time_step  Update the inc_time by the time_step
!>                        duration
subroutine update_time( self, time_step )

  implicit none

  class( jedi_increment_type ), intent(inout) :: self
  type( jedi_duration_type ),      intent(in) :: time_step

  self%inc_time = self%inc_time + time_step

end subroutine update_time

!> @brief    Set the IO clock to the time specified by the input datetime
!>
!> @param [in] new_datetime  The datetime to be used to update the LFRic clock
subroutine set_clock( self, new_time )

  implicit none

  class( jedi_increment_type ), intent(inout) :: self
  type( jedi_datetime_type ),      intent(in) :: new_time

  type( jedi_duration_type )        :: time_difference
  type( jedi_duration_type )        :: time_step
  type( model_clock_type ), pointer :: io_clock
  logical( l_def )                  :: clock_running

  io_clock => self%geometry%get_clock()
  call time_step%init( int( io_clock%get_seconds_per_step(), kind=i_def ) )

  time_difference = new_time - self%inc_time

  if ( time_difference == 0_i_def ) then
    write ( log_scratch_space, '(A)' ) &
      "New time is the same as the current time"
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  else if ( time_difference < 0_i_def ) then
    write ( log_scratch_space, '(A)' ) &
      "The xios clock can not go backwards."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Tick the clock to the required time.
  clock_running = .true.
  do while ( new_time%is_ahead( self%inc_time ) )
    clock_running = io_clock%tick()
    self%inc_time = self%inc_time + time_step
  end do

  ! Check the clock is still running
  if ( .not. clock_running ) then
    write ( log_scratch_space, '(A)' ) &
      "State::set_clock::The LFRic IO clock has stopped."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

end subroutine set_clock

!> @brief    The jedi_increment_type finalizer
!>
subroutine jedi_increment_destructor( self )

  implicit none

  type( jedi_increment_type ), intent(inout) :: self

  nullify(self%geometry)
  if ( allocated(self%fields ) ) then
    deallocate(self%fields )
  end if

end subroutine jedi_increment_destructor

!! For testing

!> Print the 1st point in each field
subroutine print_field( self )

  implicit none

  class( jedi_increment_type ), intent(inout) :: self

  ! Local
  real(real64), pointer :: atlas_data_ptr(:,:)
  integer(i_def)        :: ivar
  character(str_def)    :: iso_datetime

  ! Printing data
  call log_event( "Increment print ----", LOG_LEVEL_INFO )
  call self%inc_time%to_string( iso_datetime )
  write ( log_scratch_space, '(2A)' ) 'Time: ', iso_datetime
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  do ivar = 1, self%field_meta_data%get_n_variables()
    atlas_data_ptr => self%fields(ivar)%get_data()
    write ( log_scratch_space, '(2A,F20.10)' ) &
      trim(self%field_meta_data%get_variable_name(ivar)), &
      ", atlas_data_ptr(1,1) = ", atlas_data_ptr(1,1)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do

end subroutine print_field

end module jedi_increment_mod
