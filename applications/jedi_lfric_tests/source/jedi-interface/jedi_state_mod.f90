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
!>          (Atlas fields here), and iii) storage of modeldb instance to be
!>          used for model time stepping.
module jedi_state_mod

  use, intrinsic :: iso_fortran_env,  only : real64
  use atlas_field_emulator_mod,       only : atlas_field_emulator_type
  use atlas_field_interface_mod,      only : atlas_field_interface_type
  use constants_mod,                  only : i_def, l_def, str_def
  use field_collection_mod,           only : field_collection_type
  use driver_modeldb_mod,             only : modeldb_type
  use io_context_mod,                 only : io_context_type
  use jedi_geometry_mod,              only : jedi_geometry_type
  use jedi_lfric_datetime_mod,        only : jedi_datetime_type
  use jedi_lfric_duration_mod,        only : jedi_duration_type
  use jedi_lfric_field_meta_mod,      only : jedi_lfric_field_meta_type
  use jedi_setup_field_meta_data_mod, only : setup_field_meta_data
  use log_mod,                        only : log_event,          &
                                             log_scratch_space,  &
                                             LOG_LEVEL_INFO,     &
                                             LOG_LEVEL_ERROR
  use model_clock_mod,                only : model_clock_type
  use namelist_collection_mod,        only : namelist_collection_type
  use namelist_mod,                   only : namelist_type

  implicit none

  private

type, public :: jedi_state_type
  private

  !> These fields emulate a set of Atlas fields are and used purely for testing
  type ( atlas_field_emulator_type ), allocatable :: fields(:)

  !> An object that stores the field meta-data associated with the fields
  type( jedi_lfric_field_meta_type )              :: field_meta_data

  !> Interface field linking the Atlas emulator fields and LFRic fields in the
  !> modeldb (to do field copies)
  type( atlas_field_interface_type ), allocatable :: fields_to_modeldb(:)

  !> Logical that indicates the state has an instance of a modeldb that has been created
  logical                                         :: has_a_modeldb = .false.

  !> Modeldb that stores the fields to propagate
  type( modeldb_type ),  public                   :: modeldb

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

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the LFRic fields in the modeldb
  procedure, private :: setup_interface_to_modeldb

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

  !> Update the curent time
  procedure, public :: update_time

  !> Print the field statistics to log out
  procedure, public :: print

  !> Copy the data in the LFRic fields stored in the modeldb to the internal
  !> Atlas field emulators
  procedure, public :: from_modeldb

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
!> @param [in] geometry          The geometry object required to construct the
!>                               state
!> @param [in] config            The configuration object including the
!>                               required information to construct a state and
!>                               read a file to initialise the fields
!> @param [in] modeldb_filename  The location of the modeldb configuration file
subroutine state_initialiser_read( self, &
                                   geometry, &
                                   configuration, &
                                   modeldb_filename )

  use jedi_lfric_nl_modeldb_driver_mod, only : initialise_modeldb

  implicit none

  class( jedi_state_type ),           intent(inout) :: self
  type( jedi_geometry_type ), target, intent(in)    :: geometry
  type( namelist_collection_type ),   intent(in)    :: configuration
  character(len=*),         optional, intent(in)    :: modeldb_filename

  ! Local
  type( namelist_type ), pointer :: jedi_state_config
  logical( l_def )               :: use_pseudo_model

  jedi_state_config => configuration%get_namelist('jedi_state')

  call self%state_initialiser( geometry, jedi_state_config )

  ! Initialise the Atlas field emulators via the modeldb or the
  ! io_collection
  call jedi_state_config%get_value( 'use_pseudo_model', use_pseudo_model )

  if ( use_pseudo_model ) then
    ! We are not running the non-linear model so read the model fields from
    ! a file obtained by running the nl model offline and do a copy from
    ! io_collection to the fields stored internally (self%fields).
    call self%read_file( self%state_time )
  else
    ! We are running the non-linear model so initialise the modeldb and
    ! do a copy from model fields stored inside the modeldb
    if ( .not. present(modeldb_filename) ) then
      ! A configuration file is required if the pseudo model is not being used.
      log_scratch_space = &
        'jedi_state_mod: modeldb configuration file name not supplied.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    endif
    call initialise_modeldb( "non-linear modeldb", modeldb_filename, &
                             geometry%get_mpi_comm(), self%modeldb )
    self%has_a_modeldb = .true.
    call self%setup_interface_to_modeldb()
    call self%from_modeldb()
  end if

end subroutine state_initialiser_read

!> @brief    Initialiser for jedi_state_type
!>
!> @param [in] geometry The geometry object required to construct the state
!> @param [in] config   A configuration object including the required
!>                      information to construct a state
subroutine state_initialiser( self, geometry, config )

  use fs_continuity_mod,    only : W3, Wtheta

  implicit none

  class( jedi_state_type ),        intent(inout) :: self
  type( jedi_geometry_type ), target, intent(in) :: geometry
  type( namelist_type ),              intent(in) :: config

  ! Local
  integer(i_def)                  :: n_horizontal
  integer(i_def)                  :: n_levels
  integer(i_def)                  :: n_layers
  integer(i_def)                  :: ivar
  integer(i_def)                  :: n_variables
  integer(i_def)                  :: fs_id
  logical(l_def)                  :: twod_field
  character(str_def)              :: state_time
  character(str_def), allocatable :: variables(:)

  ! Setup
  call config%get_value( 'state_time', state_time )
  call self%state_time%init( state_time )

  call config%get_value( 'variables', variables )
  call setup_field_meta_data( self%field_meta_data, variables )

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

    call self%fields(ivar)%initialise( &
                                n_levels, &
                                n_horizontal, &
                                self%geometry%get_mpi_comm(), &
                                self%field_meta_data%get_variable_name(ivar) )
  end do

  ! Setup the io_collection (empty ready for read/write operations)
  call self%io_collection%initialise(name = 'io_collection', table_len=100)

end subroutine state_initialiser

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
subroutine read_file( self, read_time )

  use jedi_lfric_io_update_mod, only: update_io_field_collection
  use lfric_xios_read_mod,      only: read_state

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  type( jedi_datetime_type ),  intent(in) :: read_time

  ! Local
  character( len=str_def ), allocatable :: variable_names(:)
  character( len=str_def )              :: current_datetime
  class( io_context_type ),     pointer :: context_ptr
  character( len=* ),         parameter :: file_prefix="read_"

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

  ! Copy modeldb fields to the Atlas field emulators
  call self%field_meta_data%get_variable_names(variable_names)

  call self%set_from_field_collection(variable_names, self%io_collection)

end subroutine read_file

!> Write model fields to file
!>
!> @param [in] write_time  The datetime to be write to the file
subroutine write_file( self, write_time )

  use jedi_lfric_io_update_mod, only : update_io_field_collection
  use lfric_xios_write_mod,     only : write_state

  implicit none

  class( jedi_state_type ),   intent(inout) :: self
  type( jedi_datetime_type ),    intent(in) :: write_time

  ! Local
  character( len=str_def ), allocatable :: variable_names(:)
  character( len=str_def )              :: current_datetime
  class( io_context_type ),     pointer :: context_ptr
  character( len=* ), parameter         :: file_prefix="write_"

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

!> @brief    Setup fields_to_modeldb variable that enables copying between
!>           Atlas field emulators and the LFRic fields in the modeldb
!>
subroutine setup_interface_to_modeldb( self )

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
  logical(l_def)            :: is_2d

  type( field_collection_type ), pointer :: depository

  nullify(depository)
  depository => self%modeldb%fields%get_field_collection("depository")

  n_variables = self%field_meta_data%get_n_variables()

  ! Allocate space for the interface fields
  if ( allocated( self%fields_to_modeldb ) ) then
    deallocate( self%fields_to_modeldb )
  endif
  allocate( self%fields_to_modeldb( n_variables ) )

  ! Link the Atlas emulator fields with lfric fields
  call self%geometry%get_horizontal_map(horizontal_map_ptr)
  do ivar=1, n_variables

    ! Get the required data
    call get_model_field( self%field_meta_data%get_variable_name(ivar), &
                          depository, lfric_field_ptr )

    atlas_data_ptr => self%fields(ivar)%get_data()

    is_2d = self%field_meta_data%get_variable_is_2d(ivar)

    call self%fields_to_modeldb(ivar)%initialise( atlas_data_ptr,     &
                                                  horizontal_map_ptr, &
                                                  lfric_field_ptr, is_2d=is_2d )

  end do

end subroutine setup_interface_to_modeldb

!> @brief    Copy from modeldb to the Atlas field emulators
!>
subroutine from_modeldb( self )

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def) :: ivar

  !> @todo Will need some sort of transform for winds and
  !>       possibly other higher order elements.
  !>
  !>       call transform_winds(...)

  ! copy to the Atlas emulator fields
  do ivar = 1, size(self%fields_to_modeldb)
    call self%fields_to_modeldb(ivar)%copy_from_lfric()
  end do

end subroutine from_modeldb

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
  logical(l_def)            :: is_2d

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
    is_2d = self%field_meta_data%get_variable_is_2d(ivar)
    call interface_fields(ivar)%initialise( atlas_data_ptr,     &
                                            horizontal_map_ptr, &
                                            lfric_field_ptr, is_2d=is_2d )
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
  do while ( new_time > self%state_time )
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

  use jedi_lfric_nl_modeldb_driver_mod, only : finalise_modeldb

  implicit none

  type( jedi_state_type ), intent(inout) :: self

  self%geometry => null()
  if ( allocated(self%fields ) ) deallocate(self%fields )
  if ( allocated( self%fields_to_modeldb ) ) deallocate( self%fields_to_modeldb )
  if (self%has_a_modeldb) then
    call finalise_modeldb( self%modeldb )
    self%has_a_modeldb = .false.
  endif
  call self%io_collection%clear()

end subroutine jedi_state_destructor

!> Print the field statistics to log out
!>
subroutine print( self )

  implicit none

  class( jedi_state_type ), target, intent(inout) :: self

  ! Local
  integer(i_def)                            :: ivar
  type (atlas_field_emulator_type), pointer :: atlas_field_ptr

  ! Printing data
  call log_event( "State print ----", LOG_LEVEL_INFO )
  call self%state_time%print()
  do ivar = 1, self%field_meta_data%get_n_variables()
    ! Get Atlas field pointer
    atlas_field_ptr => self%fields(ivar)
    ! Print current field rms, max and min
    write ( log_scratch_space, '(A,3(A,E22.15))' ) &
      trim(self%field_meta_data%get_variable_name(ivar)), &
      ", RMS: ", atlas_field_ptr%root_mean_square(), &
      ", Max: ", atlas_field_ptr%maximum(), &
      ", Min: ", atlas_field_ptr%minimum()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do

end subroutine print

end module jedi_state_mod
