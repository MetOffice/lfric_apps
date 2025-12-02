!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the JEDI Geometry emulator class.
!>
!> @details This module holds a JEDI Geometry emulator class that includes only
!>           the functionality that is required to support interface emulation.
!>
module jedi_geometry_mod

  use, intrinsic :: iso_fortran_env, only : real64

  use calendar_mod,                  only : calendar_type
  use constants_mod,                 only : i_def, l_def, str_def, &
                                            r_second, i_timestep
  use extrusion_mod,                 only : extrusion_type, TWOD
  use io_context_mod,                only : io_context_type
  use jedi_lfric_driver_time_mod,    only : jedi_lfric_init_time, &
                                            jedi_lfric_final_time
  use jedi_lfric_mesh_interface_mod, only : is_mesh_cubesphere,        &
                                            get_cubesphere_resolution, &
                                            get_layer_ncells,          &
                                            get_lonlat,                &
                                            get_sigma_w3_levels,       &
                                            get_sigma_wtheta_levels,   &
                                            get_stretching_height
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_lfric_file_meta_mod,      only : jedi_lfric_file_meta_type
  use jedi_lfric_io_setup_mod,       only : initialise_io
  use lfric_mpi_mod,                 only : lfric_mpi_type, &
                                            lfric_comm_type
  use log_mod,                       only : log_event, LOG_LEVEL_ERROR
  use mesh_mod,                      only : mesh_type
  use mesh_collection_mod,           only : mesh_collection
  use model_clock_mod,               only : model_clock_type
  use namelist_collection_mod,       only : namelist_collection_type
  use namelist_mod,                  only : namelist_type

  implicit none

  private

type, public :: jedi_geometry_type
  private

  !> The data map between external field data and LFRic fields
  integer( kind = i_def ), allocatable  :: horizontal_map(:)
  !> the LFRic field dimensions
  integer( kind = i_def )               :: n_horizontal
  !> Comm and mesh_name
  integer( kind = i_def )               :: mpi_comm
  character( len=str_def )              :: mesh_name
  !> IO clock, context and calendar
  type( model_clock_type ), allocatable :: io_clock
  class( io_context_type ), allocatable :: io_context
  class( calendar_type ),   allocatable :: calendar
  logical( kind=l_def )                 :: io_setup_increment=.false.

contains

  !> Field initialiser.
  procedure, public :: initialise

  !> Getters
  procedure, public  :: get_clock
  procedure, public  :: get_mpi_comm
  procedure, public  :: get_mesh_name
  procedure, public  :: get_mesh
  procedure, public  :: get_twod_mesh
  procedure, public  :: get_io_context
  procedure, public  :: get_io_setup_increment
  procedure, public  :: get_n_horizontal
  procedure, public  :: get_n_layers
  procedure, public  :: get_horizontal_map

  !> IO Setup
  procedure, private :: setup_io

  !> Finalizer
  final             :: jedi_geometry_destructor

end type jedi_geometry_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_geometry_type
!>
subroutine initialise( self, mpi_comm, configuration )
  ! Access config directly until modeldb ready
  use driver_mesh_mod,           only: init_mesh
  use driver_config_mod,         only: init_config
  use jedi_lfric_mesh_setup_mod, only: initialise_mesh
  use jedi_lfric_tests_mod,      only: jedi_lfric_tests_required_namelists

  implicit none

  class( jedi_geometry_type ), intent(inout) :: self
  integer( kind=i_def ),          intent(in) :: mpi_comm
  type(namelist_collection_type), intent(in) :: configuration

  ! Local
  type(mesh_type), pointer     :: mesh
  type(lfric_mpi_type)         :: mpi_obj
  type(namelist_type), pointer :: geometry_configuration
  integer                      :: i_horizontal
  real(real64)                 :: domain_height
  real(real64)                 :: stretching_height
  real(real64), allocatable    :: lonlat(:,:),          &
                                  sigma_W3_levels(:),   &
                                  sigma_Wtheta_levels(:)

  ! Save the mpi_comm
  self%mpi_comm = mpi_comm

  ! Setup mesh
  mpi_obj = self%get_mpi_comm()
  call initialise_mesh( self%mesh_name, configuration, mpi_obj )

  geometry_configuration => configuration%get_namelist('jedi_geometry')
  ! Setup the IO
  call self%setup_io( geometry_configuration )

  ! @todo: The geometry should read some fields: orog, height, ancils

  ! The following is testing the mesh interface

  ! Set target mesh for all functions in the interface
  mesh => self%get_mesh()

  if ( .not. is_mesh_cubesphere(mesh) ) then
    call log_event( "Working mesh is not a cubesphere", LOG_LEVEL_ERROR )
  end if

  ! Get grid size and layers
  self%n_horizontal = get_layer_ncells(mesh)

  ! Create horizontal_map
  lonlat = get_lonlat(mesh)
  allocate( self%horizontal_map( self%n_horizontal ) )

  ! For mock purposes return sequential map
  do i_horizontal=1,self%n_horizontal
    self%horizontal_map( i_horizontal ) = i_horizontal
  end do

  ! Here JEDI deals with physical coordinates
  domain_height = mesh%get_domain_top()
  sigma_W3_levels = get_sigma_w3_levels(mesh)
  sigma_Wtheta_levels = get_sigma_wtheta_levels(mesh)
  stretching_height = get_stretching_height()

end subroutine initialise

!> @brief    Get the number of horizontal points
!>
!> @return n_horizontal The number of horizontal points
function get_n_horizontal(self) result(n_horizontal)

  implicit none

  class( jedi_geometry_type ), intent(in) :: self
  integer( kind=i_def )                   :: n_horizontal

  n_horizontal = self%n_horizontal

end function get_n_horizontal

!> @brief    Get the number of model layers
!>
!> @return n_layers The number of model layers
function get_n_layers(self) result(n_layers)

  implicit none

  class( jedi_geometry_type ), intent(in) :: self
  integer( kind=i_def )                   :: n_layers

  type(mesh_type), pointer :: mesh
  mesh => self%get_mesh()

  n_layers = mesh%get_nlayers()

end function get_n_layers

!> @brief    Get a pointer to the horizontal map
!>
!> @return horizontal_map A pointer to the map providing the horizontal index
!>                        of the Atlas field
subroutine get_horizontal_map(self, horizontal_map)

  implicit none

  class( jedi_geometry_type ), target, intent(in) :: self
  integer( kind=i_def ), pointer,   intent(inout) :: horizontal_map(:)

  horizontal_map => self % horizontal_map

end subroutine get_horizontal_map

!> @brief    Get the IO clock
!>
!> @return io_clock A pointer to the IO clock
function get_clock(self) result(io_clock)

  implicit none

    class( jedi_geometry_type ), target, intent(in) :: self
    type( model_clock_type ), pointer               :: io_clock

    io_clock => self%io_clock

end function get_clock

!> @brief    Get the mpi communicator
!>
!> @return mpi_comm The mpi communicator
function get_mpi_comm(self) result(mpi_obj)

implicit none

  class( jedi_geometry_type ), intent(in) :: self
  type( lfric_mpi_type )                  :: mpi_obj
  type(lfric_comm_type)                   :: lfric_comm

  call lfric_comm%set_comm_mpi_val( self%mpi_comm )
  call mpi_obj%initialise( lfric_comm )

end function get_mpi_comm

!> @brief    Get the mesh mesh_name
!>
!> @return mesh_name The mesh mesh_name
function get_mesh_name(self) result(mesh_name)

  implicit none

    class( jedi_geometry_type ), intent(in) :: self
    character( len=str_def )                :: mesh_name

    mesh_name = self%mesh_name

end function get_mesh_name

!> @brief    Get the 3D mesh object
!>
!> @return mesh The 3D mesh object
function get_mesh(self) result(mesh)

  implicit none

    class( jedi_geometry_type ), intent(in) :: self
    type( mesh_type ), pointer              :: mesh

    mesh => mesh_collection%get_mesh(self%mesh_name)

end function get_mesh

!> @brief    Get the 2D mesh object
!>
!> @return mesh The 2D mesh object
function get_twod_mesh(self) result(twod_mesh)

  implicit none

    class( jedi_geometry_type ), intent(in) :: self
    type( mesh_type ), pointer              :: mesh
    type( mesh_type ), pointer              :: twod_mesh

    mesh => mesh_collection%get_mesh(self%mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

end function get_twod_mesh

!> @brief    Get the IO context
!>
!> @return context_ptr Pointer to the current context object
function get_io_context(self) result(context_ptr)

   implicit none

   class( jedi_geometry_type ), target, intent(in) :: self
   class( io_context_type ), pointer :: context_ptr

   context_ptr => self%io_context

end function get_io_context

!> @brief    Get the io_setup_increment logical
!>
!> @return io_setup_increment Logical that defines if the increment io has been
!>                            setup
function get_io_setup_increment(self) result(io_setup_increment)

  implicit none

  class( jedi_geometry_type ), target, intent(in) :: self
  logical( kind=l_def )                           :: io_setup_increment

  io_setup_increment = self%io_setup_increment

end function get_io_setup_increment

!> @brief    Private method to setup the IO for the application
!>
!> @param [in] configuration A configuration object containing the IO options
subroutine setup_io(self, configuration)

  implicit none

  class( jedi_geometry_type ), intent(inout) :: self
  type( namelist_type ),       intent(in)    :: configuration

  ! Local
  real( kind=r_second )      :: time_step
  integer( kind=i_timestep ) :: duration
  character( len= str_def )  :: calender_start

  ! Local
  character( len=str_def ) :: xios_id
  character( len=str_def ) :: io_mode_str
  character( len=str_def ) :: field_group_id
  character( len=str_def ) :: file_name
  character( len=str_def ) :: context_name
  character( len=str_def ) :: io_path_state_write
  character( len=str_def ) :: io_path_state_read
  character( len=str_def ) :: io_path_inc_read
  character( len=str_def ) :: io_time_step_str
  character( len=str_def ) :: io_calender_start_str
  integer( kind=i_def )    :: freq

  type( jedi_lfric_file_meta_type ), allocatable :: file_meta_data(:)
  type( jedi_duration_type )                     :: io_time_step
  type( jedi_datetime_type )                     :: io_calender_start

  ! Create IO clock and setup IO
  call configuration%get_value( 'io_time_step', io_time_step_str )
  call io_time_step%init( io_time_step_str )
  call io_time_step%get_duration( duration )
  time_step = real( duration, r_second )

  call configuration%get_value( 'io_calender_start', io_calender_start_str )
  call io_calender_start%init( io_calender_start_str )
  call io_calender_start%to_string(calender_start)

  call jedi_lfric_init_time( time_step, calender_start, &
                             self%io_clock, self%calendar )

  ! Allocate file_meta depending on how many files are required
  call configuration%get_value( 'io_setup_increment', self%io_setup_increment )
  if ( self%io_setup_increment ) then
    allocate( file_meta_data(3) )
  else
    allocate( file_meta_data(2) )
  endif

  ! Setup IO files: i) state read, ii) state write and iii) increment read
  !                 (if requested)
  ! Note, xios_id is arbitrary but has to be unique within a context

  ! Read state
  call configuration%get_value( 'io_path_state_read', io_path_state_read )
  file_name = io_path_state_read
  xios_id = "read_model_data"
  io_mode_str = "read"
  field_group_id = "read_fields"
  freq=1_i_def
  call file_meta_data(1)%initialise( file_name,   &
                                     xios_id,     &
                                     io_mode_str, &
                                     freq,        &
                                     field_group_id )
  ! Write state
  call configuration%get_value( 'io_path_state_write', io_path_state_write )
  file_name = io_path_state_write
  xios_id = "write_model_data"
  io_mode_str = "write"
  field_group_id = "write_fields"
  freq=1_i_def
  call file_meta_data(2)%initialise( file_name,   &
                                     xios_id,     &
                                     io_mode_str, &
                                     freq,        &
                                     field_group_id )

  ! Read increment if required
  if ( self%io_setup_increment ) then
    call configuration%get_value( 'io_path_inc_read', io_path_inc_read )
    file_name = io_path_inc_read
    xios_id = "read_inc_model_data"
    io_mode_str = "read"
    field_group_id = "read_inc_fields"
    freq=1_i_def
    call file_meta_data(3)%initialise( file_name,   &
                                       xios_id,     &
                                       io_mode_str, &
                                       freq,        &
                                       field_group_id )
  endif

  ! Setup XIOS with the files defined by file_meta_data
  context_name = "jedi_context"
  call initialise_io( context_name,         &
                      self%get_mpi_comm(),  &
                      file_meta_data,       &
                      self%get_mesh_name(), &
                      self%calendar,        &
                      self%io_context,      &
                      self%io_clock )

  ! Tick out of initialisation state
  if ( .not. self%io_clock%tick() ) then
    call log_event( 'The LFRic IO has stopped.', LOG_LEVEL_ERROR )
  end if

end subroutine setup_io

!> @brief    Finalizer for jedi_geometry_type
!>
subroutine jedi_geometry_destructor(self)

  implicit none

  type(jedi_geometry_type), intent(inout)    :: self

  if ( allocated(self % horizontal_map) ) deallocate(self % horizontal_map)
  if ( allocated(self % io_context) ) deallocate(self % io_context)
  call jedi_lfric_final_time( self%io_clock, self%calendar )
  self%io_setup_increment = .false.

end subroutine jedi_geometry_destructor

end module jedi_geometry_mod
