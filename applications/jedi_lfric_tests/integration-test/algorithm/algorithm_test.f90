!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   The top level program for the da dev algorithm
!!         integration tests.
!>@details Sets up and runs the integration tests specified in
!!         algorithms_test.py. Currently only
!!         jedi_lfric_increment_alg_mod.x90 is tested, this algorithm
!!         takes a real field and adds one to its value.
program algorithm_test

  use add_mesh_map_mod,        only: assign_mesh_maps
  use configuration_mod,       only: final_configuration, &
                                     read_configuration
  use constants_mod,           only: i_def, r_def, str_def, l_def
  use create_mesh_mod,         only: create_extrusion, create_mesh
  use test_algorithm_mod,      only: test_algorithm_finalise,   &
                                     test_algorithm_initialise, &
                                     test_jedi_lfric_increment_alg
  use driver_collections_mod,  only: init_collections, final_collections
  use driver_mesh_mod,         only: init_mesh
  use extrusion_mod,           only: extrusion_type,         &
                                     uniform_extrusion_type, &
                                     PRIME_EXTRUSION, TWOD
  use halo_comms_mod,          only: initialise_halo_comms, &
                                     finalise_halo_comms
  use lfric_mpi_mod,           only: global_mpi, &
                                     create_comm, destroy_comm, &
                                     lfric_comm_type
  use log_mod,                 only: log_event,          &
                                     initialise_logging, &
                                     finalise_logging,   &
                                     LOG_LEVEL_ERROR,    &
                                     LOG_LEVEL_INFO
  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type

  use base_mesh_config_mod, only: GEOMETRY_SPHERICAL, &
                                  GEOMETRY_PLANAR

  implicit none

  ! MPI communicator
  type(lfric_comm_type) :: comm

  ! Number of processes and local rank
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  type(namelist_collection_type), save :: configuration

  ! Variables used for parsing command line arguments
  integer :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag

  ! Flags which determine the tests that will be carried out
  logical :: do_test_jedi_lfric_increment_alg_mod = .false.

  class(extrusion_type),        allocatable :: extrusion
  type(uniform_extrusion_type), allocatable :: extrusion_2d

  character(str_def) ::  prime_mesh_name

  logical(l_def) :: apply_partition_check

  integer(i_def) :: geometry
  integer(i_def) :: stencil_depth
  integer(i_def) :: method
  integer(i_def) :: number_of_layers
  real(r_def)    :: domain_bottom
  real(r_def)    :: domain_height
  real(r_def)    :: scaled_radius

  type(namelist_type), pointer :: base_mesh_nml => null()
  type(namelist_type), pointer :: planet_nml    => null()
  type(namelist_type), pointer :: extrusion_nml => null()

  integer(i_def) :: i
  integer(i_def), parameter :: one_layer = 1_i_def

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Setup for all tests
  ! Usage message to print
  character(len=256)       :: usage_message

  character(str_def)       :: base_mesh_names(1)
  character(str_def), allocatable :: twod_names(:)
  real(r_def), parameter   :: tolerance = 1.0e-3_r_def

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Communicators and Logging Setup
  ! Initialise MPI communicatios and get a valid communicator
  call create_comm(comm)

  ! Save the communicator for later use
  call global_mpi%initialise(comm)

  ! Initialise halo functionality
  call initialise_halo_comms(comm)

  total_ranks = global_mpi%get_comm_size()
  local_rank  = global_mpi%get_comm_rank()

  call initialise_logging( comm%get_comm_mpi_val(), 'jedi_lfric_test_alg' )

  call log_event( 'da dev alg testing running ...', LOG_LEVEL_INFO )

  ! Parse command line parameters
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()

  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 2) then
    write(usage_message,*) "Usage: ",trim(program_name), &
      " <namelist filename> "       // &
      " test_XXX with XXX in { "    // &
      " jedi_lfric_increment_alg_mod, " // &
      " } "
    call log_event( trim(usage_message), LOG_LEVEL_ERROR )
  end if

  call get_command_argument( 1, dummy, length, status )
  allocate( character(length) :: filename )
  call get_command_argument( 1, filename, length, status )

  call get_command_argument( 2, dummy, length, status )
  allocate(character(length)::test_flag)
  call get_command_argument( 2, test_flag, length, status )

  ! Choose test case depending on flag provided in the first command
  ! line argument
  select case (trim(test_flag))
  case ("test_jedi_lfric_increment_alg_mod")
    do_test_jedi_lfric_increment_alg_mod = .true.
  case default
    call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  ! Setup configuration, mesh, and fem
  call configuration%initialise( program_name, table_len=10 )
  call read_configuration( filename, configuration )

  call init_collections()

  !--------------------------------------
  ! 0.0 Extract namelist variables
  !--------------------------------------
  base_mesh_nml => configuration%get_namelist('base_mesh')
  planet_nml    => configuration%get_namelist('planet')
  extrusion_nml => configuration%get_namelist('extrusion')

  call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
  call base_mesh_nml%get_value( 'geometry', geometry )
  call extrusion_nml%get_value( 'method', method )
  call extrusion_nml%get_value( 'domain_height', domain_height )
  call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
  call planet_nml%get_value( 'scaled_radius', scaled_radius )

  base_mesh_nml => null()
  planet_nml    => null()
  extrusion_nml => null()

  !--------------------------------------
  ! 1.0 Create the meshes
  !--------------------------------------
  base_mesh_names(1) = prime_mesh_name
  allocate( twod_names, source=base_mesh_names )

  !--------------------------------------
  ! 1.1 Create the required extrusions
  !--------------------------------------
  select case (geometry)
  case (geometry_planar)
    domain_bottom = 0.0_r_def
  case (geometry_spherical)
    domain_bottom = scaled_radius
  case default
    call log_event("Invalid geometry for mesh initialisation", LOG_LEVEL_ERROR)
  end select
  allocate( extrusion, source=create_extrusion( method,           &
                                                domain_height,    &
                                                domain_bottom,    &
                                                number_of_layers, &
                                                PRIME_EXTRUSION ) )

  extrusion_2d = uniform_extrusion_type( domain_height, &
                                         domain_bottom, &
                                         one_layer, TWOD )

  !-------------------------------------------------------------------------
  ! 1.2 Create the required meshes
  !-------------------------------------------------------------------------
  stencil_depth = 1
  apply_partition_check = .false.
  call init_mesh( configuration,              &
                  local_rank, total_ranks,    &
                  base_mesh_names, extrusion, &
                  stencil_depth,              &
                  apply_partition_check )

  do i=1, size(twod_names)
    twod_names(i) = trim(twod_names(i))//'_2d'
  end do
  call create_mesh( base_mesh_names, extrusion_2d, &
                    alt_name=twod_names )
  call assign_mesh_maps(twod_names)


  !-------------------------------------------------------------------------
  ! Tests
  !-------------------------------------------------------------------------
  call test_algorithm_initialise(prime_mesh_name) ! fem

  if ( do_test_jedi_lfric_increment_alg_mod ) then
    call test_jedi_lfric_increment_alg(tolerance)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Finalise and close down

  call test_algorithm_finalise()

  if (allocated(program_name)) deallocate(program_name)
  if (allocated(filename))     deallocate(filename)
  if (allocated(test_flag))    deallocate(test_flag)

  call log_event( 'da dev alg functional testing completed ...', LOG_LEVEL_INFO )

  ! (note that the order of the calls of the following finalisers matters)

  call final_collections()
  call final_configuration()

  call finalise_halo_comms()
  call global_mpi%finalise()
  call configuration%clear()
  call destroy_comm()

  call finalise_logging()

end program algorithm_test
