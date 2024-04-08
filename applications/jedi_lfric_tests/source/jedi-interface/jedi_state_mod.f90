!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the JEDI State emulator class.
!>
!> @details This module holds a JEDI State emulator class that includes only
!>          the functionality required by the LFRic-JEDI model interface on the
!>          LFRic-API. This includes i) read/write to a defined set of fields,
!>          ii) ability to interoperate between LFRic fields and JEDI fields
!>          (Atlas fields here), and iii) storage of model_data instance to be
!>          used for IO and model time stepping.
module jedi_state_mod

  use, intrinsic :: iso_fortran_env, only : real64
  use atlas_field_emulator_mod,      only : atlas_field_emulator_type
  use atlas_field_interface_mod,     only : atlas_field_interface_type
  use calendar_mod,                  only : calendar_type
  use constants_mod,                 only : i_def, l_def, str_def
  use driver_model_data_mod,         only : model_data_type
  use driver_time_mod,               only : init_time, final_time
  use field_collection_mod,          only : field_collection_type
  use io_context_mod,                only : io_context_type
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use jedi_state_config_mod,         only : jedi_state_config_type
  use jedi_lfric_field_meta_mod,     only : jedi_lfric_field_meta_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_INFO,     &
                                            LOG_LEVEL_ERROR
  use model_clock_mod,               only : model_clock_type

  implicit none

  private

type, public :: jedi_state_type
  private

  !> These fields emulate a set of Atlas fields are and used purely for testing
  type ( atlas_field_emulator_type ), allocatable :: fields(:)

  !> An object that stores the field meta-data associated with the fields
  type( jedi_lfric_field_meta_type )              :: field_meta_data

  !> Interface field linking the Atlas emulator fields and LFRic fields in the
  !> model data (to do field copies)
  type( atlas_field_interface_type ), allocatable :: fields_to_model_data(:)

  !> Model data that stores the fields to propagate
  !> (will be subsumed into modelDB)
  type( model_data_type ),  public                :: model_data

  !> Model clock associated with the model_data (will be subsumed into modelDB)
  type( model_clock_type ), allocatable, public   :: model_clock

  !> Model calendar associated with the model_data (will be subsumed into modelDB)
  class( calendar_type ), allocatable             :: calendar

  !> Field collection to perform IO
  type( field_collection_type ), public           :: io_collection

  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formatting for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'
  type( jedi_datetime_type )                      :: state_time

  !> The jedi_geometry object
  type( jedi_geometry_type ), pointer, public     :: geometry => null()

contains

  !> Jedi state initialiser.
  procedure :: initialise => state_initialiser_read
  procedure :: state_initialiser

  !> Initialise the model_data
  procedure, private :: initialise_model_data

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the LFRic fields in the model_data
  procedure, private :: setup_interface_to_model_data

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the fields in a LFRic field collection
  procedure, private :: setup_interface_to_field_collection

  !> Get Atlas field from the bundle using the fieldname
  procedure, private :: get_field_data

  !> Set the internal Atlas field emulators from a input field_collection
  procedure, public :: set_from_field_collection

  !> Get via the internal Atlas field emulators via a copy to a field_collection
  procedure, public :: get_to_field_collection

  !> Return the curent time
  procedure, public :: valid_time

  !> Read model fields from file into fields
  procedure, public :: read_file

  !> Write model fields to file
  procedure, public :: write_file

  !> Create the model_data
  procedure, public :: create_model_data

  !> Update the curent time
  procedure, public :: update_time

  !> Print field
  procedure, public :: print_field

  !> Copy the data in the LFRic fields stored in the model_data to the internal
  !> Atlas field emulators
  procedure, public :: from_model_data

  !> Finalizer
  final             :: jedi_state_destructor

end type jedi_state_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! Methods required by the JEDI (OOPS) model interface
!-------------------------------------------------------------------------------

!> @brief    Initialiser via read for jedi_state_type
!>
!> @param [in] geometry The geometry object required to construct the state
!> @param [in] config   The configuration object including the required
!>                      information to construct a state and read a file to
!>                      initialise the fields
subroutine state_initialiser_read( self, geometry, config )

  implicit none

  class( jedi_state_type ),           intent(inout) :: self
  type( jedi_geometry_type ), target, intent(in)    :: geometry
  type( jedi_state_config_type ),     intent(inout) :: config

  call self%state_initialiser( geometry, config )

  ! Initialise the Atlas field emulators via the model_data or the
  ! io_collection
  if ( .not. config%use_pseudo_model ) then
    ! This calls the models initialise method that populates model_data and
    ! does a copy from model_data to the fields
    call self%initialise_model_data()
  else
    ! We are not running the non-linear model so read the file directly and
    ! do a copy from io_collection to the fields
    call self%read_file( config%state_time, config%read_file_prefix )
  end if

end subroutine state_initialiser_read

!> @brief    Initialiser for jedi_state_type
!>
!> @param [in] geometry The geometry object required to construct the state
!> @param [in] config   A configuration object including the required
!>                      information to construct a state
subroutine state_initialiser( self, geometry, config )

  use fs_continuity_mod,     only : W3, Wtheta

  implicit none

  class( jedi_state_type ),        intent(inout) :: self
  type( jedi_geometry_type ), target, intent(in) :: geometry
  type( jedi_state_config_type ),  intent(inout) :: config

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
  self%state_time = config%state_time
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

  ! If running the model, create model data and link to fields .
  if ( .not. config%use_pseudo_model ) then
    call init_time( self%model_clock, self%calendar )
    call self%create_model_data()
  end if

end subroutine state_initialiser

!> @brief    A method to initialise the model_data and copy into the Atlas
!>           field emulators
!>
subroutine initialise_model_data( self )

  use jedi_lfric_fake_nl_mod, only : initialise_fake_nl_model_data

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Read the state into the LFRic model data
  call initialise_fake_nl_model_data( self%model_data )

  ! Copy model_data fields to the Atlas field emulators
  call self%from_model_data()

end subroutine initialise_model_data

!> @brief    Returns the current time of the increment
!>
!> @return   time   The current time of the state
function valid_time( self ) result(time)

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  type( jedi_datetime_type )              :: time

  time = self%state_time

end function valid_time

!> @brief    A method to update the internal Atlas field emulators
!>
!> @param [in] read_time   The datetime to be read from the file
!> @param [in] file_prefix Character array that specifies the file to read from
subroutine read_file( self, read_time, file_prefix )

  use jedi_lfric_io_update_mod, only: update_io_field_collection
  use lfric_xios_read_mod,      only: read_state

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  type( jedi_datetime_type ),  intent(in) :: read_time
  character(len=*),            intent(in) :: file_prefix

  ! Local
  character( len=str_def ), allocatable :: variable_names(:)
  character( len=str_def )              :: current_datetime
  class( io_context_type ),     pointer :: context_ptr

  ! Ensure the JEDI-IO context is set to current
  context_ptr => self%geometry%get_io_context()
  call context_ptr%set_current()

  ! Set the clock to the desired read time
  call set_clock(self, read_time)

  call read_time%to_string( current_datetime )
  call log_event( "jedi_state_mod: reading file at &
                & " // current_datetime, LOG_LEVEL_INFO )

  ! Ensure the io_collection contains the variables defined in the list
  ! stored in the type
  call update_io_field_collection( self%io_collection,            &
                                   self%geometry%get_mesh(),      &
                                   self%geometry%get_twod_mesh(), &
                                   self%field_meta_data )

  ! Read the state into the io_collection
  call read_state( self%io_collection, prefix=file_prefix )

  ! Copy model_data fields to the Atlas field emulators
  call self%field_meta_data%get_variable_names(variable_names)

  call self%set_from_field_collection(variable_names, self%io_collection)

end subroutine read_file

!> Write model fields to file
!>
!> @param [in] write_time  The datetime to be write to the file
!> @param [in] file_prefix Character array that specifies the file to write to
subroutine write_file( self, write_time, file_prefix )

  use jedi_lfric_io_update_mod, only : update_io_field_collection
  use lfric_xios_write_mod,     only : write_state

  implicit none

  class( jedi_state_type ),   intent(inout) :: self
  type( jedi_datetime_type ),    intent(in) :: write_time
  character(len=*),              intent(in) :: file_prefix

  ! Local
  character( len=str_def ), allocatable :: variable_names(:)
  character( len=str_def )              :: current_datetime
  class( io_context_type ),     pointer :: context_ptr

  ! Ensure the JEDI-IO context is set to current
  context_ptr => self%geometry%get_io_context()
  call context_ptr%set_current()

  ! Set the clock to the desired write time
  call set_clock( self, write_time )
  call write_time%to_string( current_datetime )
  call log_event( "jedi_state_mod: writing file at &
                & " // current_datetime, LOG_LEVEL_INFO )

  ! Ensure the io_collection contains only the required data for writing,
  ! e.g. remove redundant fields.
  call update_io_field_collection( self%io_collection,            &
                                   self%geometry%get_mesh(),      &
                                   self%geometry%get_twod_mesh(), &
                                   self%field_meta_data )
  ! Copy the Atlas field emulators to the io_collection fields
  call self%field_meta_data%get_variable_names(variable_names)

  call self%get_to_field_collection(variable_names, self%io_collection)

  ! Write the state into the io_collection
  call write_state( self%io_collection, prefix=file_prefix )

  call log_event( "jedi_state_mod: write_file: finished", LOG_LEVEL_INFO )

end subroutine write_file

!> @brief    Create an instance of the model_data for jedi_state_type
!>
subroutine create_model_data( self )

  use jedi_lfric_fake_nl_mod,   only : create_fake_nl_model_data

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Create model data and then link the Atlas fields to them
  call create_fake_nl_model_data( self%geometry%get_mesh(),      &
                                  self%geometry%get_twod_mesh(), &
                                  self%model_data )
  call self%setup_interface_to_model_data()

end subroutine create_model_data

!> @brief    Setup fields_to_model_data variable that enables copying between
!>           Atlas field emulators and the LFRic fields in the model_data
!>
subroutine setup_interface_to_model_data( self )

  use jedi_lfric_utils_mod,  only: get_model_field
  use field_mod,             only: field_type

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def)            :: ivar
  type(field_type), pointer :: lfric_field_ptr
  real(real64),     pointer :: atlas_data_ptr(:,:)
  integer(i_def),   pointer :: horizontal_map_ptr(:)
  integer(i_def)            :: n_variables

  type( field_collection_type ), pointer :: depository

  nullify(depository)
  depository => self%model_data%get_field_collection("depository")

  n_variables = self%field_meta_data%get_n_variables()

  ! Allocate space for the interface fields
  if ( allocated( self%fields_to_model_data ) ) then
    deallocate( self%fields_to_model_data )
  endif
  allocate( self%fields_to_model_data( n_variables ) )

  ! Link the Atlas emulator fields with lfric fields
  call self%geometry%get_horizontal_map(horizontal_map_ptr)
  do ivar=1, n_variables

    ! Get the required data
    call get_model_field( self%field_meta_data%get_variable_name(ivar), &
                          depository, lfric_field_ptr )

    atlas_data_ptr => self%fields(ivar)%get_data()

    call self%fields_to_model_data(ivar)%initialise( atlas_data_ptr,     &
                                                     horizontal_map_ptr, &
                                                     lfric_field_ptr )

  end do

end subroutine setup_interface_to_model_data

!> @brief    Copy from model_data to the Atlas field emulators
!>
subroutine from_model_data( self )

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def) :: ivar

  !> @todo Will need some sort of transform for winds and
  !>       possibly other higher order elements.
  !>
  !>       call transform_winds(model_data)

  ! copy to the Atlas emulator fields
  do ivar = 1, size(self%fields_to_model_data)
    call self%fields_to_model_data(ivar)%copy_from_lfric()
  end do

end subroutine from_model_data

!> @brief    Setup Atlas-LFRic interface_fields that enables copying
!>           between Atlas field emulators and the LFRic fields in
!>           io_collection
!>
!> @param [inout] interface_fields  The Atlas-LFRic interafce object
!> @param [in]    variable_names    The name of the fields to setup
!> @param [inout] field_collection  The field collection to link the Atlas fields to
subroutine setup_interface_to_field_collection( self,             &
                                                interface_fields, &
                                                variable_names,   &
                                                field_collection )

  use jedi_lfric_utils_mod,  only: get_model_field
  use field_mod,             only: field_type

  implicit none

  class( jedi_state_type ),           intent(inout) :: self
  type( atlas_field_interface_type ), intent(inout) :: interface_fields(:)
  character( len=str_def ),              intent(in) :: variable_names(:)
  type( field_collection_type ),      intent(inout) :: field_collection

  ! Local
  integer(i_def)            :: ivar
  type(field_type), pointer :: lfric_field_ptr
  real(real64),     pointer :: atlas_data_ptr(:,:)
  integer(i_def),   pointer :: horizontal_map_ptr(:)
  integer(i_def)            :: n_variables
  logical                   :: all_variables_exists

  ! Check that the state contains all the required fields
  all_variables_exists = &
        self%field_meta_data%check_variables_exist( variable_names )
  if ( .not. all_variables_exists ) then
    log_scratch_space = &
        'The Atlas state does not conatin all the required fields.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  ! Link the Atlas emulator fields with lfric fields
  call self%geometry%get_horizontal_map( horizontal_map_ptr )
  n_variables=size( variable_names )
  do ivar = 1, n_variables

    ! Get the required data and setup field interface
    call get_model_field( variable_names(ivar), &
                          field_collection, lfric_field_ptr )
    call self%get_field_data(variable_names(ivar), atlas_data_ptr)
    call interface_fields(ivar)%initialise( atlas_data_ptr,     &
                                            horizontal_map_ptr, &
                                            lfric_field_ptr )
  end do

end subroutine setup_interface_to_field_collection

!> @brief  Method to get a pointer to the field for a given variable name
!>
!> @param [in]  variable_names  The name of the field to retrieve
!> @param [out] atlas_data_ptr  The field pointer to the atlas_emulator array
subroutine get_field_data( self, variable_name, atlas_data_ptr )

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  character( len=str_def ),    intent(in) :: variable_name
  real( real64 ), pointer,    intent(out) :: atlas_data_ptr(:,:)

  ! Local
  integer( i_def ) :: ivar
  integer( i_def ) :: n_variables
  logical( l_def ) :: field_found

  ! Find the field requested
  n_variables = size(self%fields)

  field_found = .false.
  do ivar = 1, n_variables
    if (variable_name==self%fields(ivar)%get_field_name()) then
      atlas_data_ptr => self%fields(ivar)%get_data()
      field_found = .true.
      return
    endif
  end do

  if (.not. field_found) then
    nullify(atlas_data_ptr)
    write ( log_scratch_space, '(3A)' ) &
      "The Atlas field with the name ", trim(variable_name), &
      ", could not be found."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

end subroutine get_field_data

!> @brief    Set the internal Atlas field emulators by copying the fields
!>           stored in a field_collection
!>
!> @param [in]    variable_names    The name of the fields to setup
!> @param [inout] field_collection  The field collection to copy from
subroutine set_from_field_collection( self, &
                                      variable_names, &
                                      field_collection )

  implicit none

  class( jedi_state_type ),      intent(inout) :: self
  character( len=str_def ),         intent(in) :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer(i_def) :: ivar
  integer(i_def) :: n_variables
  type( atlas_field_interface_type ), allocatable :: interface_fields(:)

  ! Allocate space for the interface fields
  n_variables = size(variable_names)
  allocate( interface_fields(n_variables) )

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_field_collection( interface_fields, &
                                                 variable_names, &
                                                 field_collection )

  ! Copy the LFRic fields in the field_collection to the Atlas field emulators
  do ivar = 1, n_variables
    call interface_fields(ivar)%copy_from_lfric()
  end do

end subroutine set_from_field_collection

!> @brief    Get a copy of the internal Atlas field emulators by copying into a
!>           field_collection
!>
!> @param [in]    variable_names    The name of the fields to setup
!> @param [inout] field_collection  The field collection to copy to
subroutine get_to_field_collection( self, variable_names, field_collection )

  implicit none

  class( jedi_state_type ),      intent(inout) :: self
  character( len=str_def ),         intent(in) :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer(i_def) :: ivar
  integer(i_def) :: n_variables
  type( atlas_field_interface_type ), allocatable :: interface_fields(:)

  ! Allocate space for the interface fields
  n_variables = size( variable_names )
  allocate( interface_fields(n_variables) )

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_field_collection( interface_fields, &
                                                 variable_names, &
                                                 field_collection )

  ! Copy the LFRic fields in the field_collection to the Atlas field emulators
  do ivar = 1, n_variables
    call interface_fields(ivar)%copy_to_lfric()
  end do

end subroutine get_to_field_collection

!> @brief    Update the state time by a single time-step
!>
!> @param [in] time_step  Update the state time by the time_step
!>                        duration
subroutine update_time( self, time_step )

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  type( jedi_duration_type ),  intent(in) :: time_step

  self%state_time = self%state_time + time_step

end subroutine update_time

!> @brief    Set the IO clock to the time specified by the input time
!>
!> @param [in] new_time  The time to be used to update the LFRic clock
subroutine set_clock( self, new_time )

  implicit none

  class( jedi_state_type ),   intent(inout) :: self
  type( jedi_datetime_type ), intent(in)    :: new_time

  ! Local
  type( jedi_duration_type )        :: time_difference
  type( jedi_duration_type )        :: time_step
  type( model_clock_type ), pointer :: io_clock
  logical( kind=l_def )             :: clock_running

  nullify(io_clock)
  io_clock => self%geometry%get_clock()

  call time_step%init( int( io_clock%get_seconds_per_step(), kind=i_def ) )
  time_difference = new_time - self%state_time

  if ( time_difference == 0_i_def ) then
    write ( log_scratch_space, '(A)' ) &
      "New time is the same as the current state time"
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  else if ( time_difference < 0_i_def ) then
    write ( log_scratch_space, '(A)' ) &
      "The xios clock can not go backwards."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Tick the clock to the required time.
  clock_running = .true.
  do while ( new_time%is_ahead( self%state_time ) )
    clock_running = io_clock%tick()
    self%state_time = self%state_time + time_step
  end do

  ! Check the clock is still running
  if ( .not. clock_running ) then
    write ( log_scratch_space, '(A)' ) &
      "State::set_clock::The LFRic IO clock has stopped."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

end subroutine set_clock

!> @brief    The jedi_state_type finalizer
!>
subroutine jedi_state_destructor( self )

  implicit none

  type( jedi_state_type ), intent(inout) :: self

  self%geometry => null()
  if ( allocated(self%fields ) ) deallocate(self%fields )
  if ( allocated(self%model_clock ) ) then
    call final_time( self%model_clock, self%calendar )
  endif

end subroutine jedi_state_destructor

!! For testing

!> Print the 1st point in each field
subroutine print_field( self )

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  real(real64), pointer :: atlas_data_ptr(:,:)
  integer(i_def)        :: ivar
  character(str_def)    :: iso_datetime

  ! Printing data
  call log_event( "State print ----", LOG_LEVEL_INFO )
  call self%state_time%to_string( iso_datetime )
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

end module jedi_state_mod
