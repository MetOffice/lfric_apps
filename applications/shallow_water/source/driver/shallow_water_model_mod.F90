!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!!        and the shallow_water model simulation.
module shallow_water_model_mod

  use add_mesh_map_mod,               only: assign_mesh_maps
  use assign_orography_field_mod,     only: assign_orography_field
  use check_configuration_mod,        only: get_required_stencil_depth
  use sci_checksum_alg_mod,           only: checksum_alg
  use conservation_algorithm_mod,     only: conservation_algorithm
  use constants_mod,                  only: i_def, str_def, r_def, &
                                            PRECISION_REAL, l_def
  use convert_to_upper_mod,           only: convert_to_upper
  use create_mesh_mod,                only: create_mesh, create_extrusion
  use derived_config_mod,             only: set_derived_config
  use driver_fem_mod,                 only: init_fem, final_fem
  use driver_io_mod,                  only: init_io, &
                                            filelist_populator
  use driver_modeldb_mod,             only: modeldb_type
  use driver_mesh_mod,                only: init_mesh
  use extrusion_mod,                  only: extrusion_type,         &
                                            uniform_extrusion_type, &
                                            PRIME_EXTRUSION, TWOD
  use field_mod,                      only: field_type
  use field_parent_mod,               only: write_interface
  use field_collection_mod,           only: field_collection_type
  use sci_geometric_constants_mod,    only: get_chi_inventory, &
                                            get_panel_id_inventory
  use inventory_by_mesh_mod,          only: inventory_by_mesh_type
  use lfric_xios_file_mod,            only: lfric_xios_file_type
  use linked_list_mod,                only: linked_list_type
  use log_mod,                        only: log_event,         &
                                            log_set_level,     &
                                            log_scratch_space, &
                                            LOG_LEVEL_INFO,    &
                                            LOG_LEVEL_ERROR
  use mesh_collection_mod,            only: mesh_collection
  use mesh_mod,                       only: mesh_type
  use minmax_tseries_mod,             only: minmax_tseries,      &
                                            minmax_tseries_init, &
                                            minmax_tseries_final
  use namelist_collection_mod,        only: namelist_collection_type
  use namelist_mod,                   only: namelist_type
  use runtime_constants_mod,          only: create_runtime_constants
  use shallow_water_setup_io_mod,     only: init_shallow_water_files
  use xios,                           only: xios_update_calendar

  ! Configuration modules
  use base_mesh_config_mod, only: GEOMETRY_PLANAR, &
                                  GEOMETRY_SPHERICAL

  implicit none

  private
  public :: initialise_infrastructure, &
            initialise_model,          &
            finalise_infrastructure,   &
            finalise_model

  contains

  !=============================================================================
  !> @brief Initialises the infrastructure and sets up constants used
  !!        by the model.
  !> @param[in]     program_name   Identifier given to the model being run
  !> @param[inout]  modeldb        The ModelDB object
  subroutine initialise_infrastructure( program_name, modeldb)

    implicit none

    character(*),           intent(in)    :: program_name
    class(modeldb_type),    intent(inout) :: modeldb

    type(inventory_by_mesh_type),  pointer :: chi_inventory      => null()
    type(inventory_by_mesh_type),  pointer :: panel_id_inventory => null()
    procedure(filelist_populator), pointer :: files_init_ptr     => null()

    character(len=*),   parameter   :: io_context_name = "shallow_water"

    character(str_def), allocatable :: base_mesh_names(:)
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

    integer(i_def), parameter :: one_layer = 1_i_def
    integer(i_def) :: i

    !=======================================================================
    ! 0.0 Extract configuration variables
    !=======================================================================
    base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
    planet_nml    => modeldb%configuration%get_namelist('planet')
    extrusion_nml => modeldb%configuration%get_namelist('extrusion')

    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    call base_mesh_nml%get_value( 'geometry', geometry )
    call extrusion_nml%get_value( 'method', method )
    call extrusion_nml%get_value( 'domain_height', domain_height )
    call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
    call planet_nml%get_value( 'scaled_radius', scaled_radius )

    base_mesh_nml => null()
    planet_nml    => null()
    extrusion_nml => null()

    !-------------------------------------------------------------------------
    ! Initialise aspects of the infrastructure
    !-------------------------------------------------------------------------
    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call set_derived_config( .true. )


    !-------------------------------------------------------------------------
    ! 1.0 Mesh
    !-------------------------------------------------------------------------

    !=======================================================================
    ! 1.1 Determine the required meshes
    !=======================================================================

    ! Meshes that require a prime/2d extrusion
    ! ---------------------------------------------------------
    allocate(base_mesh_names(1))
    base_mesh_names(1) = prime_mesh_name

    !=======================================================================
    ! 1.2 Generate required extrusions
    !=======================================================================

    ! Extrusions for prime/2d meshes
    ! ---------------------------------------------------------
    select case (geometry)
    case (GEOMETRY_PLANAR)
      domain_bottom = 0.0_r_def
    case (GEOMETRY_SPHERICAL)
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


    !=======================================================================
    ! 1.3 Initialise mesh objects and assign InterGrid maps
    !=======================================================================

    ! Initialise prime/2d meshes
    ! ---------------------------------------------------------
    check_partitions = .false.
    stencil_depth = get_required_stencil_depth()
    call init_mesh( modeldb%configuration,       &
                    modeldb%mpi%get_comm_rank(), &
                    modeldb%mpi%get_comm_size(), &
                    base_mesh_names, extrusion,  &
                    stencil_depth, check_partitions )


    allocate( twod_names, source=base_mesh_names )
    do i=1, size(twod_names)
      twod_names(i) = trim(twod_names(i))//'_2d'
    end do
    call create_mesh( base_mesh_names, extrusion_2d, &
                      alt_name=twod_names )
    call assign_mesh_maps(twod_names)

    !-------------------------------------------------------------------------
    ! 2.0 Build the FEM function spaces and coordinate fields
    !-------------------------------------------------------------------------
    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()
    call init_fem(mesh_collection, chi_inventory, panel_id_inventory)

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
    ! domain and context

    files_init_ptr => init_shallow_water_files
    call init_io( io_context_name, prime_mesh_name, modeldb, &
                  chi_inventory, panel_id_inventory,         &
                  populate_filelist=files_init_ptr )

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    call create_runtime_constants()

    deallocate(base_mesh_names)
    nullify(chi_inventory, panel_id_inventory)

  end subroutine initialise_infrastructure

  !=============================================================================
  !> @brief Initialises the shallow_water application.
  !> @param[in]     mesh       The primary mesh
  !> @param[in,out] modeldb The working data set for the model run
  subroutine initialise_model( mesh,    &
                               modeldb )

    use swe_timestep_alg_mod, only: swe_timestep_alg_init

    implicit none

    type(mesh_type),    pointer, intent(in)    :: mesh
    type(modeldb_type),          intent(inout) :: modeldb

    type(field_collection_type), pointer :: prognostic_fields => null()

    ! Prognostic fields
    type(field_type), pointer :: wind => null()
    type(field_type), pointer :: geopot => null()
    type(field_type), pointer :: buoyancy => null()
    type(field_type), pointer :: q => null()

    call log_event( 'shallow_water: Initialising miniapp ...', LOG_LEVEL_INFO )

    prognostic_fields => modeldb%fields%get_field_collection("prognostics")

    call prognostic_fields%get_field("wind", wind)
    call prognostic_fields%get_field("buoyancy", buoyancy)
    call prognostic_fields%get_field("geopot", geopot)
    call prognostic_fields%get_field("q", q)

    ! Initialise transport and shallow water model
    call swe_timestep_alg_init( mesh, wind, geopot, buoyancy, q )

    call log_event( 'shallow_water: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine initialise_model


  !=============================================================================
  !> @brief Finalises infrastructure and constants used by the model.
  subroutine finalise_infrastructure()

    implicit none

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------
    call final_fem()

  end subroutine finalise_infrastructure

  !=============================================================================
  !> @brief Finalise the shallow_water application.
  !> @param[in,out] modeldb   The working data set for the model run
  !> @param[in]     program_name An identifier given to the model begin run
  subroutine finalise_model( modeldb, &
                             program_name )

    use swe_timestep_alg_mod, only: swe_timestep_alg_final

    implicit none

    type(modeldb_type), intent(inout) :: modeldb
    character(*),       intent(in)    :: program_name

    type( field_collection_type ), pointer :: prognostic_fields => null()

    type(field_type), pointer :: wind => null()
    type(field_type), pointer :: geopot => null()
    type(field_type), pointer :: q => null()

    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    prognostic_fields => modeldb%fields%get_field_collection("prognostics")

    call prognostic_fields%get_field('wind', wind)
    call prognostic_fields%get_field('geopot', geopot)
    call prognostic_fields%get_field('q', q)

    ! Checksums
    call checksum_alg( program_name, wind, 'wind', geopot, 'geopot', q, 'q' )

    ! Finalise transport
    call swe_timestep_alg_final()

  end subroutine finalise_model

end module shallow_water_model_mod
