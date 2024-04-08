!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Main program for running the CMA operator tests in cma_test_mod.x90.
!>@details To run a particular test, provide the appropriate command line
!>         argument of the form --test_XXX where a list of possible values
!>         for XXX is given below. Running with --test_all will carry out all
!>         tests.
!>
!>         Any configuration data is read from the namelist file
!>         cma_test_configuration.nml.
!>
!>         Note that the test of the diag_DhMDhT term currently only works
!>         reliably sequentially.

program cma_test

  use, intrinsic :: iso_fortran_env,  only : real64

  use check_configuration_mod,        only : get_required_stencil_depth
  use cma_test_algorithm_mod,         only : cma_test_init,                  &
                                             test_cma_apply_mass_p,          &
                                             test_cma_apply_mass_v,          &
                                             test_cma_apply_div_v,           &
                                             test_cma_multiply_div_v_mass_v, &
                                             test_cma_multiply_grad_v_div_v, &
                                             test_cma_add,                   &
                                             test_cma_apply_inv,             &
                                             test_cma_diag_DhMDhT
  use constants_mod,                  only : i_def, r_def, i_def, l_def, &
                                             r_solver, pi, str_def
  use derived_config_mod,             only : set_derived_config
  use extrusion_mod,                  only : extrusion_type
  use mpi_mod,                        only : global_mpi, &
                                             create_comm, destroy_comm
  use field_mod,                      only : field_type
  use fs_continuity_mod,              only : W0,W1,W2,W3
  use function_space_mod,             only : function_space_type
  use halo_comms_mod,                 only : initialise_halo_comms, &
                                             finalise_halo_comms
  use configuration_mod,              only : read_configuration, &
                                             ensure_configuration
  use driver_collections_mod,         only : init_collections, final_collections
  use driver_mesh_mod,                only : init_mesh
  use gungho_extrusion_mod,           only : create_extrusion
  use log_mod,                        only : log_event,          &
                                             log_scratch_space,  &
                                             initialise_logging, &
                                             finalise_logging,   &
                                             LOG_LEVEL_ERROR,    &
                                             LOG_LEVEL_INFO
  use mesh_mod,                       only : mesh_type
  use mesh_collection_mod,            only : mesh_collection
  use namelist_collection_mod,        only : namelist_collection_type
  use namelist_mod,                   only : namelist_type

  use base_mesh_config_mod,           only : GEOMETRY_SPHERICAL

  implicit none

  ! MPI communicator
  integer(kind=i_def) :: comm

  ! Number of processes and local rank
  integer(kind=i_def) :: total_ranks, local_rank

  ! Filename to read namelist from
  character(:), allocatable :: filename

  ! Variables used for parsing command line arguments
  integer :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag
  ! Usage message to print
  character(len=256) :: usage_message
  ! The following variables are required for checking that the mesh
  ! has a suitable size (i.e. \f$dx \approx dz\f$)
  ! Pointer to mesh
  type(mesh_type), pointer :: mesh => null()
  ! Number of cells of 2d mesh (local and global),
  integer(kind=i_def) :: ncells_2d_local, ncells_2d
  ! Number of vertical layers
  integer(kind=i_def) :: nlayers
  ! vertical domain size
  real   (kind=r_def) :: domain_top
  ! Grid spacing in horizontal and vertical
  real   (kind=r_def) :: dx, dz
  character(str_def)  :: base_mesh_names(1)
  ! Variables for reading configuration from namelist file
  character(*), parameter ::       &
       required_configuration(6) = &
      (/'finite_element      ',    &
        'base_mesh           ',    &
        'formulation         ',    &
        'planet              ',    &
        'extrusion           ',    &
        'partitioning        '/)
  logical              :: okay
  logical, allocatable :: success_map(:)
  integer :: i

  class(extrusion_type), allocatable :: extrusion

  ! Flags which determine the tests that will be carried out
  logical :: do_test_apply_mass_p = .false.
  logical :: do_test_apply_mass_v = .false.
  logical :: do_test_apply_div_v = .false.
  logical :: do_test_multiply_div_v_mass_v = .false.
  logical :: do_test_multiply_grad_v_div_v = .false.
  logical :: do_test_add = .false.
  logical :: do_test_apply_inv = .false.
  logical :: do_test_diag_dhmdht = .false.

  ! Namelist and configuration variables
  type(namelist_collection_type), save :: configuration
  type(namelist_type), pointer         :: nml_obj => null()

  integer(i_def)     :: stencil_depth
  character(str_def) :: file_prefix
  character(str_def) :: prime_mesh_name
  real(r_def)        :: radius
  integer(i_def)     :: geometry
  logical(l_def)     :: prepartitioned
  logical            :: check_partitions = .false.

  ! Error tolerance for tests
  ! Note: tolerance is for r_solver = real64
  !       Tolerance at r_solver = real32: generated using spacing() in
  !           ./gungho/source/algorithm/cma_test_algorithm_mod.x90
  real(kind=r_solver) :: tolerance
  tolerance = 1.0E-12_r_solver

  ! Initialise MPI communicatios and get a valid communicator
  call create_comm(comm)

  ! Save the communicator for later use
  call global_mpi%initialise(comm)

  ! Initialise halo functionality
  call initialise_halo_comms( comm )

  total_ranks = global_mpi%get_comm_size()
  local_rank  = global_mpi%get_comm_rank()

  call initialise_logging( comm, 'cma_test' )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  call log_event( 'CMA functional testing running ...', LOG_LEVEL_INFO )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parse command line parameters
  !
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()
  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 2) then
     write(usage_message,*) "Usage: ",trim(program_name), &
          " <namelist filename> " // &
          " test_XXX with XXX in { " // &
          "apply_mass_p, "             // &
          "apply_mass_v, "             // &
          "apply_div_v, "              // &
          "multiply_div_v_mass_v, "    // &
          "multiply_grad_v_div_v, "    // &
          "add, "                      // &
          "apply_inv, "                // &
          "diag_dhmdht, "              // &
          "all"                        // &
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
  case ("test_apply_mass_p")
     do_test_apply_mass_p = .true.
  case ("test_apply_mass_v")
     do_test_apply_mass_v = .true.
  case ("test_apply_div_v")
     do_test_apply_div_v = .true.
  case ("test_multiply_div_v_mass_v")
     do_test_multiply_div_v_mass_v = .true.
  case ("test_multiply_grad_v_div_v")
     do_test_multiply_grad_v_div_v = .true.
  case ("test_add")
     do_test_add = .true.
  case ("test_apply_inv")
     do_test_apply_inv = .true.
  case ("test_diag_dhmdht")
     do_test_diag_dhmdht = .true.
  case ("test_all")
     do_test_apply_mass_p = .true.
     do_test_apply_mass_v = .true.
     do_test_apply_div_v = .true.
     do_test_multiply_div_v_mass_v = .true.
     do_test_multiply_grad_v_div_v = .true.
     do_test_add = .true.
     do_test_apply_inv = .true.
     do_test_diag_dhmdht = .true.
  case default
     call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  call configuration%initialise( program_name, table_len=10 )

  deallocate(program_name)
  deallocate(test_flag)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Load configuration from file
  !
  write( log_scratch_space, '("Reading namelist file: ", A)' ) trim(filename)
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  allocate( success_map(size(required_configuration)) )
  call read_configuration( filename, configuration )

  okay = ensure_configuration( required_configuration, success_map )
  if (.not. okay) then
     write( log_scratch_space, '(A)' ) &
            'The following required namelists were not loaded:'
     do i = 1,size(required_configuration)
        if (.not. success_map(i)) &
             log_scratch_space = trim(log_scratch_space) // ' ' &
                                 // required_configuration(i)
     end do
     call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if
  deallocate( success_map )

  call set_derived_config(.false.)

  call init_collections()

  if (configuration%namelist_exists('planet')) then
    nml_obj => configuration%get_namelist('planet')
    call nml_obj%get_value( 'radius', radius )
  end if

  if (configuration%namelist_exists('base_mesh')) then
    nml_obj => configuration%get_namelist('base_mesh')
    call nml_obj%get_value( 'file_prefix', file_prefix )
    call nml_obj%get_value( 'prepartitioned', prepartitioned )
    call nml_obj%get_value( 'geometry', geometry )
    call nml_obj%get_value( 'prime_mesh_name', prime_mesh_name )
  end if


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise

  call log_event( 'Initialising harness', LOG_LEVEL_INFO )

  base_mesh_names(1) = prime_mesh_name
  allocate( extrusion, source=create_extrusion() )

  stencil_depth = get_required_stencil_depth()
  check_partitions = .false.
  call init_mesh( configuration,              &
                  local_rank, total_ranks,    &
                  base_mesh_names, extrusion, &
                  stencil_depth, check_partitions )

  ! Work out grid spacing, which should be of order 1
  mesh => mesh_collection%get_mesh(prime_mesh_name)
  ncells_2d_local = mesh%get_ncells_2d()

  ! Ensure that a spherical geometry is used (otherwise tests are too simple)
  if (geometry /= GEOMETRY_SPHERICAL) then
     call log_event( "Geometry has to be spherical", &
                     LOG_LEVEL_ERROR )
  end if

  ! Work out total number of cells
  call global_mpi%global_sum(ncells_2d_local, ncells_2d)

  ! Check that the grid spacings are of order 1 (between 0.1 and 10).
  ! Otherwise the derivative and mass terms in the
  ! Helmholtz-operator are not well balanced and errors could go unnoticed
  nlayers = mesh%get_nlayers()
  domain_top = mesh%get_domain_top()
  dx = sqrt(4.*pi/ncells_2d)*radius
  dz = domain_top / nlayers

  if ( ( dx < 0.1_r_def) .or. ( dx > 10.0_r_def) ) then
     call log_event( "Average dx has to be in range 0.1 ... 10.0", &
                     LOG_LEVEL_ERROR )
  end if
  if ( ( dz < 0.1_r_def) .or. ( dz > 10.0_r_def) ) then
     call log_event( "Average dz has to be in range 0.1 ... 10.0", &
                     LOG_LEVEL_ERROR )
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Do tests
  !
  call log_event( 'Initialising test', LOG_LEVEL_INFO )

  ! Initialise CMA test module
  call cma_test_init( mesh )

  ! Run all requested tests
  if (do_test_apply_mass_p)             call test_cma_apply_mass_p(tolerance)
  if (do_test_apply_mass_v)             call test_cma_apply_mass_v(tolerance)
  if (do_test_apply_div_v)              call test_cma_apply_div_v(tolerance)
  if (do_test_multiply_div_v_mass_v)    call test_cma_multiply_div_v_mass_v(tolerance)
  if (do_test_multiply_grad_v_div_v)    call test_cma_multiply_grad_v_div_v(tolerance)
  if (do_test_add)                      call test_cma_add(tolerance)
  if (do_test_apply_inv)                call test_cma_apply_inv(tolerance)
  if (do_test_diag_dhmdht)              call test_cma_diag_DhMDhT(tolerance)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Finalise and close down
  !
  call log_event( ' CMA functional testing completed ...', LOG_LEVEL_INFO )

  call final_collections ()

  ! Finalise halo functionality
  call finalise_halo_comms()

  ! Finalise MPI communications
  call global_mpi%finalise()
  call destroy_comm()

  ! Finalise the logging system
  call finalise_logging()

end program cma_test
