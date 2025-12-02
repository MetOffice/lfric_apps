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

  use, intrinsic :: iso_fortran_env,  only : real64
  use atlas_field_emulator_mod,       only : atlas_field_emulator_type
  use atlas_field_interface_mod,      only : atlas_field_interface_type
  use constants_mod,                  only : i_def, l_def, str_def
  use field_collection_mod,           only : field_collection_type
  use io_context_mod,                 only : io_context_type
  use jedi_geometry_mod,              only : jedi_geometry_type
  use jedi_lfric_datetime_mod,        only : jedi_datetime_type
  use jedi_lfric_duration_mod,        only : jedi_duration_type
  use jedi_lfric_field_meta_mod,      only : jedi_lfric_field_meta_type
  use jedi_setup_field_meta_data_mod, only : setup_field_meta_data
  use log_mod,                        only : log_event,          &
                                             log_scratch_space,  &
                                             LOG_LEVEL_INFO,     &
                                             LOG_LEVEL_DEBUG,    &
                                             LOG_LEVEL_ERROR
  use model_clock_mod,                only : model_clock_type
  use namelist_collection_mod,        only : namelist_collection_type
  use namelist_mod,                   only : namelist_type

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
  procedure, private :: construct_increment
  procedure, private :: initialise_via_configuration
  procedure, private :: initialise_via_copy
  generic, public    :: initialise => initialise_via_configuration, &
                                      initialise_via_copy

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the fields in a LFRic field collection
  procedure, private :: setup_interface_to_field_collection

  !> Get Atlas field from the bundle using the fieldname
  procedure, private :: get_field_data

  !> Set the internal Atlas field emulators from a input field_collection
  procedure, public :: set_from_field_collection
  procedure, public :: set_from_field_collection_ad

  !> Get via the internal Atlas field emulators via a copy to a field_collection
  procedure, public :: get_to_field_collection
  procedure, public :: get_to_field_collection_ad

  !> Return the curent time
  procedure, public :: valid_time

  !> Read the Atlas fields from file into fields
  procedure, public :: read_file

  !> @todo Write the Atlas fields to file from the fields
  !> procedure, public :: write_file

  !> Zero the Atlas fields
  procedure, public :: zero

  !> Randomise the Atlas fields
  procedure, public :: random

  !> Compute the root mean square norm
  procedure, public :: norm

  !> Compute dot_product with a supplied input increment
  procedure, public :: dot_product_with

  !> Compute dot_product with self, with each field component scaled by itself
  !> and the relevant fields scaled by the same number.
  procedure, public :: scaled_dot_product_with_itself

  !> Update the curent time
  procedure, public :: update_time

  !> Print the field statistics to log out
  procedure, public :: print

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

!> @brief Construct and initialise fields to zero or read from file
!>
!> @param [in] geometry      The geometry object required to construct the
!>                           increment
!> @param [in] configuration The configuration object including the required
!>                           information to construct an increment and set the
!>                           data via either a file-read or setting the values
!>                           to zero
subroutine initialise_via_configuration( self, geometry, configuration )

  implicit none

  class( jedi_increment_type ),       intent(inout) :: self
  type( jedi_geometry_type ), target, intent(in)    :: geometry
  type( namelist_collection_type ),   intent(in)    :: configuration

  ! Local
  type( namelist_type ),  pointer :: increment_config
  logical( l_def )                :: initialise_via_read
  character(str_def)              :: inc_time_str
  character(str_def), allocatable :: variables(:)

  ! 1. Setup from configuration
  increment_config => configuration%get_namelist('jedi_increment')

  call increment_config%get_value( 'inc_time', inc_time_str )
  call self%inc_time%init( inc_time_str )

  call increment_config%get_value( 'variables', variables )
  call setup_field_meta_data( self%field_meta_data, variables )

  self%geometry => geometry

  ! 2. Create fields
  call self%construct_increment()

  ! 3. Initialise fields by either reading or zeroing the fields
  call increment_config%get_value( 'initialise_via_read', initialise_via_read )
  if ( initialise_via_read ) then
    if ( geometry%get_io_setup_increment() ) then
      call self%read_file( self%inc_time )
    else
      ! Fail because the file was not setup
      log_scratch_space = 'The Increment file has not been setup, check the configuration.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    endif
  else
    call self%zero()
  endif

end subroutine initialise_via_configuration

!> @brief Construct and initialise fields by copy for jedi_increment_type
!>
!> @param [in] rhs  An increment to copy construct from
!>
subroutine initialise_via_copy( self, rhs )

  implicit none

  class( jedi_increment_type ), intent(inout) :: self
  class( jedi_increment_type ),    intent(in) :: rhs

  ! Local
  integer( kind=i_def ) :: n_variables
  integer( kind=i_def ) :: ivar

  ! 1. Setup from rhs increment
  self%field_meta_data = rhs%field_meta_data
  self%inc_time = rhs%inc_time
  self%geometry => rhs%geometry

  ! 2. Create fields
  call self%construct_increment()

  ! 3. Initialise fields by copy from rhs
  n_variables = self%field_meta_data%get_n_variables()
  do ivar=1,n_variables
    self%fields(ivar) = rhs%fields(ivar)
  enddo

end subroutine initialise_via_copy

!> @brief    Construct increment fields
!>
subroutine construct_increment( self )

  use fs_continuity_mod,     only : W3, Wtheta

  implicit none

  class( jedi_increment_type ),       intent(inout) :: self

  ! Local
  integer( kind=i_def ) :: n_horizontal
  integer( kind=i_def ) :: n_levels
  integer( kind=i_def ) :: n_layers
  integer( kind=i_def ) :: ivar
  integer( kind=i_def ) :: n_variables
  integer( kind=i_def ) :: fs_id
  logical( kind=l_def ) :: twod_field

  n_variables = self%field_meta_data%get_n_variables()

  allocate( self%fields(n_variables) )

  n_horizontal=self%geometry%get_n_horizontal()
  n_layers=self%geometry%get_n_layers()

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

end subroutine construct_increment

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
subroutine read_file( self, read_time )

  use jedi_lfric_io_update_mod,      only : update_io_field_collection
  use lfric_xios_read_mod,           only : read_state

  implicit none

  class( jedi_increment_type ), intent(inout) :: self
  type( jedi_datetime_type ),      intent(in) :: read_time

  ! Local
  character( len=* ),         parameter :: file_prefix="read_inc_"
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
  integer( kind=i_def ) :: n_variables
  integer( kind=i_def ) :: ivar

  ! Zero the Atlas fields
  n_variables = self%field_meta_data%get_n_variables()
  do ivar=1,n_variables
    call self%fields(ivar)%zero()
  enddo

end subroutine zero

!> @brief A method to set all internal Atlas field emulators to random values
!>
subroutine random( self )

  implicit none

  class( jedi_increment_type ), intent(inout) :: self

  ! Local
  integer( kind=i_def ) :: n_variables
  integer( kind=i_def ) :: ivar

  ! Randomise the Atlas fields
  n_variables = self%field_meta_data%get_n_variables()
  do ivar=1,n_variables
    call self%fields(ivar)%random()
  enddo

end subroutine random

!> @brief Compute the root mean square norm
!>
function norm( self ) result(rms)

  implicit none

  class( jedi_increment_type ), intent(in) :: self
  real(real64)                             :: rms

  ! Local
  integer( kind=i_def ) :: n_variables
  integer( kind=i_def ) :: ivar
  integer( kind=i_def ) :: n_values

  n_variables = self%field_meta_data%get_n_variables()
  rms = 0.0_real64
  n_values = 0_i_def
  do ivar=1,n_variables
    rms = rms + self%fields(ivar)%sum_of_squares()
    n_values = n_values + self%fields(ivar)%number_of_points()
  enddo
  rms = sqrt( rms / real(n_values, kind=real64) )

end function norm

!> @brief Compute dot_product with a supplied input increment
!>
function dot_product_with( self, rhs ) result(dot_product)

  implicit none

  class( jedi_increment_type ), intent(in) :: self
  class( jedi_increment_type ), intent(in) :: rhs
  real( real64 )                           :: dot_product

  ! Local
  integer( kind=i_def ) :: n_variables
  integer( kind=i_def ) :: ivar

  n_variables = self%field_meta_data%get_n_variables()

  dot_product = 0.0_real64
  do ivar=1,n_variables
    dot_product = dot_product + self%fields(ivar)%dot_product_with(rhs%fields(ivar))
  enddo

end function dot_product_with

!> @brief Scaled dot product with itself.
!> @details Compute dot product with itself for each field. Then scale each
!>          field with the corresponding dot product, and also the dot
!>          products themselves (adding a tiny number to avoid divide-by-zero).
!> @returns dot_product Total scaled dot product.
function scaled_dot_product_with_itself( self ) result(dot_product)

  implicit none

  class( jedi_increment_type ), intent(inout) :: self
  real( real64 )                              :: dot_product

  ! Local
  integer( kind=i_def )     :: n_variables
  integer( kind=i_def )     :: ivar
  real( real64 )            :: dot_product_ivar
  real( real64 )            :: scale_factor
  real( real64 ), parameter :: eps = 1.0e-30_real64

  call log_event( "jedi_increment_type%scaled_dot_product_with_itself", LOG_LEVEL_DEBUG )

  n_variables = self%field_meta_data%get_n_variables()

  dot_product = 0.0_real64
  do ivar=1,n_variables
    dot_product_ivar = self%fields(ivar)%dot_product_with(self%fields(ivar))

    write( log_scratch_space, * ) "Field = ", self%fields(ivar)%get_field_name(), &
                                  ", dot product = ", dot_product_ivar
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

    ! eps avoids divide-by-zero
    scale_factor = 1.0_real64 / (dot_product_ivar + eps)
    call self%fields(ivar)%multiply_by(scale_factor)
    dot_product_ivar = dot_product_ivar * scale_factor
    dot_product = dot_product + dot_product_ivar
  end do

end function scaled_dot_product_with_itself

!------------------------------------------------------------------------------
! Local methods to support LFRic-JEDI implementation
!------------------------------------------------------------------------------

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
  integer( kind=i_def )          :: ivar
  type(field_type),      pointer :: lfric_field_ptr
  real( kind=real64 ),   pointer :: atlas_data_ptr(:,:)
  integer( kind=i_def ), pointer :: horizontal_map_ptr(:)
  integer( kind=i_def )          :: n_variables
  logical( kind=l_def )          :: all_variables_exists

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
    call self%get_field_data(variable_names(ivar), atlas_data_ptr)
    call atlas_lfric_interface_fields(ivar)%initialise( atlas_data_ptr,     &
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

  class( jedi_increment_type ), intent(inout) :: self
  character( len=str_def ),        intent(in) :: variable_name
  real( real64 ), pointer,        intent(out) :: atlas_data_ptr(:,:)

  ! Local
  integer( kind=i_def ) :: ivar
  integer( kind=i_def ) :: n_variables
  logical( kind=l_def ) :: field_found

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
!> @param [in]    variable_names    The list of variables to copy
!> @param [inout] field_collection  The fields to copy from
!>
subroutine set_from_field_collection( self, variable_names, field_collection )

  implicit none

  class( jedi_increment_type ),  intent(inout) :: self
  character( len=str_def ),      intent(in)    :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer( kind=i_def )                           :: ivar
  integer( kind=i_def )                           :: n_variables
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

!> @brief    Adjoint of: set the internal Atlas field emulators by copying the
!>           fields stored in a field_collection
!>
!> @param [in]    variable_names    The list of variables to copy
!> @param [inout] field_collection  The fields to adjoint copy from
!>
subroutine set_from_field_collection_ad( self, variable_names, field_collection )

  implicit none

  class( jedi_increment_type ),  intent(inout) :: self
  character( len=str_def ),      intent(in)    :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer( kind=i_def )                           :: ivar
  integer( kind=i_def )                           :: n_variables
  type( atlas_field_interface_type ), allocatable :: atlas_lfric_interface_fields(:)

  ! Allocate space for the interface fields
  n_variables = size(variable_names)
  allocate( atlas_lfric_interface_fields(n_variables) )

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_field_collection( atlas_lfric_interface_fields, &
                                                 variable_names, &
                                                 field_collection )

  ! Copy the LFRic fields in the field_collection to the Atlas field emulators
  do ivar = n_variables, 1, -1
    call atlas_lfric_interface_fields(ivar)%copy_from_lfric_ad()
  end do

end subroutine set_from_field_collection_ad

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
  integer( kind=i_def )                           :: ivar
  integer( kind=i_def )                           :: n_variables
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

!> @brief    Adjoint of: Get a copy of the internal Atlas field emulators by
!>           copying into a field_collection
!>
!> @param [in]    variable_names    The list of variables to copy
!> @param [inout] field_collection  The fields to copy to
!>
subroutine get_to_field_collection_ad( self, variable_names, field_collection )

  implicit none

  class( jedi_increment_type ),  intent(inout) :: self
  character( len=str_def ),         intent(in) :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer( kind=i_def )                           :: ivar
  integer( kind=i_def )                           :: n_variables
  type( atlas_field_interface_type ), allocatable :: atlas_lfric_interface_fields(:)

  ! Allocate space for the interface fields
  n_variables = size(variable_names)
  allocate( atlas_lfric_interface_fields(n_variables) )

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_field_collection( atlas_lfric_interface_fields, &
                                                 variable_names, &
                                                 field_collection )

  ! Copy the LFRic fields in the field_collection to the Atlas field emulators
  do ivar = n_variables, 1, -1
    call atlas_lfric_interface_fields(ivar)%copy_to_lfric_ad()
  end do

end subroutine get_to_field_collection_ad

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
  logical( kind=l_def  )            :: clock_running

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
  do while ( new_time > self%inc_time )
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

!> Print the field statistics to log out
!>
subroutine print( self )

  implicit none

  class( jedi_increment_type ), target, intent(inout) :: self

  ! Local
  integer(i_def)                            :: ivar
  type (atlas_field_emulator_type), pointer :: atlas_field_ptr

  ! Printing data
  call log_event( "Increment print ----", LOG_LEVEL_INFO )
  call self%inc_time%print()
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

end module jedi_increment_mod
