!-----------------------------------------------------------------------------
! (c) Crown copyright 2017-2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @mainpage solver
!> The solver API is abstract so cannot be tested unless with a particular implementation.
!! The mini-app makes an example linear operator, preconditioner and krylov subspace solver.
!! This can then be executed to test the solver api against different compilers.

program solver

  use add_mesh_map_mod,        only: assign_mesh_maps
  use constants_mod,           only: i_def, r_def, PRECISION_REAL, str_def
  use convert_to_upper_mod,    only: convert_to_upper
  use cli_mod,                 only: parse_command_line
  use create_mesh_mod,         only: create_mesh, create_extrusion
  use driver_collections_mod,  only: init_collections, final_collections
  use driver_config_mod,       only: init_config, final_config
  use driver_mesh_mod,         only: init_mesh
  use driver_fem_mod,          only: init_fem
  use driver_log_mod,          only: init_logger, final_logger
  use extrusion_mod,           only: extrusion_type,         &
                                     uniform_extrusion_type, &
                                     PRIME_EXTRUSION, TWOD
  use halo_comms_mod,          only: initialise_halo_comms, &
                                     finalise_halo_comms
  use init_solver_miniapp_mod, only: init_solver_miniapp
  use inventory_by_mesh_mod,   only: inventory_by_mesh_type
  use lfric_mpi_mod,           only: global_mpi, &
                                     create_comm, destroy_comm, &
                                     lfric_comm_type
  use field_mod,               only: field_type
  use sci_field_vector_mod,    only: field_vector_type
  use solver_miniapp_alg_mod,  only: solver_miniapp_alg
  use configuration_mod,       only: final_configuration
  use solver_miniapp_mod,      only: solver_required_namelists
  use log_mod,                 only: log_event,            &
                                     log_scratch_space,    &
                                     LOG_LEVEL_ALWAYS,     &
                                     LOG_LEVEL_ERROR,      &
                                     LOG_LEVEL_INFO
  use mesh_mod,                only: mesh_type
  use mesh_collection_mod,     only: mesh_collection
  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type
  use sci_checksum_alg_mod,    only: checksum_alg


  !------------------------------------
  ! Configuration modules
  !------------------------------------
  use base_mesh_config_mod, only: GEOMETRY_SPHERICAL, &
                                  GEOMETRY_PLANAR

  implicit none

  character(*), parameter :: program_name = 'solver'

  character(:), allocatable :: filename
  type(namelist_collection_type), SAVE :: configuration

  integer(i_def) :: total_ranks, local_rank
  type(lfric_comm_type) :: comm

  ! prognostic fields
  type(field_type),    pointer :: chi(:) => null()
  type(field_type),    pointer :: panel_id => null()
  type(field_type)             :: field_1, field_2
  type(field_vector_type)      :: fv_1

  type(mesh_type),     pointer :: mesh => null()
  type(inventory_by_mesh_type) :: chi_inventory
  type(inventory_by_mesh_type) :: panel_id_inventory

  character(str_def)              :: base_mesh_names(1)
  character(str_def), allocatable :: twod_names(:)

  class(extrusion_type),        allocatable :: extrusion
  type(uniform_extrusion_type), allocatable :: extrusion_2d

  type(namelist_type), pointer :: base_mesh_nml => null()
  type(namelist_type), pointer :: planet_nml    => null()
  type(namelist_type), pointer :: extrusion_nml => null()

  character(str_def) :: prime_mesh_name

  integer(i_def) :: stencil_depth
  integer(i_def) :: geometry
  integer(i_def) :: method
  integer(i_def) :: number_of_layers
  real(r_def)    :: domain_bottom
  real(r_def)    :: domain_height
  real(r_def)    :: scaled_radius
  logical        :: check_partitions

  integer(i_def) :: i
  integer(i_def), parameter :: one_layer = 1_i_def

  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  call parse_command_line( filename )

  ! Initialise MPI communicatios and get a valid communicator
  call create_comm(comm)

  ! Initialise halo functionality
  call initialise_halo_comms( comm )

  ! Save the commmunicator for later use
  call global_mpi%initialise(comm)

  total_ranks = global_mpi%get_comm_size()
  local_rank  = global_mpi%get_comm_rank()

  call configuration%initialise( program_name, table_len=10 )
  call init_config( filename, solver_required_namelists, &
                    configuration )
  call init_logger( comm, program_name )
  call init_collections()

  deallocate( filename )

  write(log_scratch_space,'(A)')                        &
      'Application built with '//trim(PRECISION_REAL)// &
      '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

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

  call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )


  !=======================================================================
  ! 1.0 Mesh
  !=======================================================================

  !-----------------------------------------------------------------------
  ! 1.1 Determine the required meshes
  !-----------------------------------------------------------------------
  base_mesh_names(1) = prime_mesh_name

  !-----------------------------------------------------------------------
  ! 1.2 Create the required extrusions
  !-----------------------------------------------------------------------
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

  extrusion_2d = uniform_extrusion_type( domain_bottom, &
                                         domain_bottom, &
                                         one_layer, TWOD )

  !-----------------------------------------------------------------------
  ! 1.2 Create the required meshes
  !-----------------------------------------------------------------------
  stencil_depth = 1
  check_partitions = .false.
  call init_mesh( configuration,              &
                  local_rank, total_ranks,    &
                  base_mesh_names, extrusion, &
                  stencil_depth, check_partitions )

  allocate( twod_names, source=base_mesh_names )
  do i=1, size(twod_names)
    twod_names(i) = trim(twod_names(i))//'_2d'
  end do
  call create_mesh( base_mesh_names, extrusion_2d, &
                    alt_name=twod_names )
  call assign_mesh_maps(twod_names)


  !=======================================================================
  ! 2.0 Build the FEM function spaces and coordinate fields
  !=======================================================================
  ! Create FEM specifics (function spaces and chi field)
  call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

  ! Create and initialise prognostic fields
  mesh => mesh_collection%get_mesh(prime_mesh_name)
  call init_solver_miniapp( mesh, fv_1 )

  ! Call an algorithm
  call chi_inventory%get_field_array(mesh, chi)
  call panel_id_inventory%get_field(mesh, panel_id)
  call solver_miniapp_alg( fv_1, chi, panel_id )

  ! Write out output file
  call log_event(program_name//": writing diagnostic output", LOG_LEVEL_INFO)

  ! pull the fields from the vector
  call fv_1%export_field( field_1, 1 )
  call fv_1%export_field( field_2, 2 )

  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------
  call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  ! Write checksums to file
  call checksum_alg(program_name, field_1, 'solver_field_1',field_2, 'solver_field_2')

  call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  nullify(chi, panel_id, mesh)

  ! Finalise global collections
  call final_collections()

  ! Finalise namelist configurations
  call final_configuration()

  ! Finalise halo functionality
  call finalise_halo_comms()

  call final_logger( program_name )
  call final_config()

  ! Finalise MPI communications
  call global_mpi%finalise()
  call destroy_comm()

end program solver
