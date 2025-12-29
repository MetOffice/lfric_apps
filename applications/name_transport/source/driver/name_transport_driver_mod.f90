!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>
!> @brief Drives the execution of the name_transport miniapp.
!>
module name_transport_driver_mod

  use add_mesh_map_mod,                   only: assign_mesh_maps
  use check_configuration_mod,            only: get_required_stencil_depth
  use constants_mod,                      only: i_def, l_def, &
                                                r_def, r_second, str_def
  use create_mesh_mod,                    only: create_mesh, create_extrusion
  use driver_fem_mod,                     only: init_fem
  use driver_io_mod,                      only: init_io, final_io
  use driver_mesh_mod,                    only: init_mesh
  use driver_modeldb_mod,                 only: modeldb_type
  use derived_config_mod,                 only: set_derived_config
  use diagnostics_io_mod,                 only: write_scalar_diagnostic, &
                                                write_vector_diagnostic
  use extrusion_mod,                      only: extrusion_type,              &
                                                uniform_extrusion_type,      &
                                                shifted_extrusion_type,      &
                                                double_level_extrusion_type, &
                                                PRIME_EXTRUSION, TWOD,       &
                                                SHIFTED, DOUBLE_LEVEL
  use field_mod,                          only: field_type
  use fs_continuity_mod,                  only: W3, Wtheta
  use inventory_by_mesh_mod,              only: inventory_by_mesh_type
  use lfric_xios_context_mod,             only: lfric_xios_context_type
  use lfric_xios_action_mod,              only: advance
  use log_mod,                            only: log_event,         &
                                                log_scratch_space, &
                                                LOG_LEVEL_ERROR,   &
                                                LOG_LEVEL_INFO,    &
                                                LOG_LEVEL_TRACE
  use mesh_mod,                           only: mesh_type
  use mesh_collection_mod,                only: mesh_collection
  use model_clock_mod,                    only: model_clock_type
  use namelist_mod,                       only: namelist_type
  use runtime_constants_mod,              only: create_runtime_constants
  use sci_checksum_alg_mod,               only: checksum_alg
  use sci_geometric_constants_mod,        only: get_chi_inventory,      &
                                                get_panel_id_inventory, &
                                                get_height_fe
  use timing_mod,                         only: start_timing, stop_timing, &
                                                tik, LPROF

  ! Transport algorithms
  use name_transport_init_fields_alg_mod, only: name_transport_init_fields_alg
  use name_transport_control_alg_mod,     only: name_transport_prerun_setup, &
                                                name_transport_init,         &
                                                name_transport_step,         &
                                                name_transport_final

  ! Configuration modules
  use base_mesh_config_mod,               only: GEOMETRY_PLANAR, &
                                                GEOMETRY_SPHERICAL
  use name_options_config_mod,            only: transport_density

  implicit none

  private

  public :: initialise_name_transport
  public :: step_name_transport
  public :: finalise_name_transport

  ! Prognostic fields
  type(field_type) :: density
  type(field_type) :: tracer_con
  ! Prescribed wind field
  type(field_type) :: wind

contains

  !============================================================================
  !> @brief Sets up required state in preparation for run.
  !> @param[in]      program_name  Identifier given to the model being run
  !> @param[inout]   modeldb       The modeldb object
  subroutine initialise_name_transport( program_name, modeldb )

    implicit none

    character(*),            intent(in)    :: program_name
    type(modeldb_type),      intent(inout) :: modeldb

    character(len=*),         parameter   :: xios_ctx = "transport"
    integer(kind=i_def)                   :: num_base_meshes
    type(mesh_type),              pointer :: mesh
    type(inventory_by_mesh_type), pointer :: chi_inventory
    type(inventory_by_mesh_type), pointer :: panel_id_inventory
    character(len=str_def),   allocatable :: base_mesh_names(:)
    character(len=str_def),   allocatable :: extra_io_mesh_names(:)
    type(field_type),             pointer :: height_w3
    type(field_type),             pointer :: height_wth

    type(lfric_xios_context_type),         pointer :: io_context
    class(extrusion_type),             allocatable :: extrusion
    type(uniform_extrusion_type),      allocatable :: extrusion_2d
    type(shifted_extrusion_type),      allocatable :: extrusion_shifted
    type(double_level_extrusion_type), allocatable :: extrusion_double

    character(len=str_def), allocatable :: meshes_to_shift(:)
    character(len=str_def), allocatable :: meshes_to_double(:)
    character(len=str_def), allocatable :: twod_names(:)
    character(len=str_def), allocatable :: shifted_names(:)
    character(len=str_def), allocatable :: double_names(:)
    character(len=str_def)              :: prime_mesh_name

    logical(kind=l_def) :: prepartitioned
    logical(kind=l_def) :: apply_partition_check

    integer(kind=i_def) :: geometry
    integer(kind=i_def) :: stencil_depth
    real(kind=r_def)    :: domain_bottom
    real(kind=r_def)    :: domain_height
    real(kind=r_def)    :: scaled_radius
    integer(kind=i_def) :: method
    integer(kind=i_def) :: number_of_layers
    logical(kind=l_def) :: nodal_output_on_w3
    logical(kind=l_def) :: write_diag
    logical(kind=l_def) :: use_xios_io

    type(namelist_type), pointer :: base_mesh_nml
    type(namelist_type), pointer :: extrusion_nml
    type(namelist_type), pointer :: planet_nml
    type(namelist_type), pointer :: io_nml

    integer(i_def) :: i
    integer(i_def), parameter :: one_layer = 1_i_def

    !=======================================================================
    ! 0.0 Extract configuration variables
    !=======================================================================

    base_mesh_nml   => modeldb%configuration%get_namelist('base_mesh')
    extrusion_nml   => modeldb%configuration%get_namelist('extrusion')
    planet_nml      => modeldb%configuration%get_namelist('planet')
    io_nml          => modeldb%configuration%get_namelist('io')

    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    call base_mesh_nml%get_value( 'geometry', geometry )
    call base_mesh_nml%get_value( 'prepartitioned', prepartitioned )
    call extrusion_nml%get_value( 'method', method )
    call extrusion_nml%get_value( 'domain_height', domain_height )
    call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
    call planet_nml%get_value( 'scaled_radius', scaled_radius )
    call io_nml%get_value( 'nodal_output_on_w3', nodal_output_on_w3 )
    call io_nml%get_value( 'write_diag', write_diag )
    call io_nml%get_value( 'use_xios_io', use_xios_io )

    base_mesh_nml   => null()
    extrusion_nml   => null()
    planet_nml      => null()
    io_nml          => null()

    !-----------------------------------------------------------------------
    ! Initialise infrastructure
    !-----------------------------------------------------------------------
    call set_derived_config( .true. )

    !=======================================================================
    ! 1.0 Mesh
    !=======================================================================

    !=======================================================================
    ! 1.1 Determine the required meshes (only primal mesh for NAME)
    !=======================================================================
    num_base_meshes = 1

    ! 1.1a Meshes that require a prime/2d extrusion
    ! ---------------------------------------------------------
    allocate(base_mesh_names(num_base_meshes))
    base_mesh_names(1) = prime_mesh_name

    ! 1.1b Meshes the require a shifted extrusion
    ! ---------------------------------------------------------
    allocate(meshes_to_shift,  source=base_mesh_names)

    ! 1.1c Meshes that require a double-level extrusion
    ! ---------------------------------------------------------
    allocate(meshes_to_double, source=base_mesh_names)


    !=======================================================================
    ! 1.2 Generate required extrusions
    !=======================================================================

    ! 1.2a Extrusions for prime/2d meshes
    ! ---------------------------------------------------------
    select case (geometry)
    case (geometry_planar)
      domain_bottom = 0.0_r_def
    case (geometry_spherical)
      domain_bottom = scaled_radius
    case default
      call log_event( "Invalid geometry for mesh initialisation", &
                      log_level_error )
    end select

    allocate( extrusion, source=create_extrusion( method,           &
                                                  domain_height,    &
                                                  domain_bottom,    &
                                                  number_of_layers, &
                                                  prime_extrusion ) )

    allocate( twod_names, source=base_mesh_names )
    do i=1, size(twod_names)
      twod_names(i) = trim(twod_names(i))//'_2d'
    end do
    extrusion_2d = uniform_extrusion_type( domain_bottom, &
                                           domain_bottom, &
                                           one_layer, TWOD )

    ! 1.2b Extrusions for shifted meshes
    ! ---------------------------------------------------------
    if ( allocated(meshes_to_shift) ) then
      if ( size(meshes_to_shift) > 0 ) then
        extrusion_shifted = shifted_extrusion_type(extrusion)
      end if
    end if

    ! 1.2c Extrusions for double-level meshes
    ! ---------------------------------------------------------
    if ( allocated(meshes_to_double) ) then
      if ( size(meshes_to_double) > 0 ) then
        extrusion_double = double_level_extrusion_type(extrusion)
      end if
    end if

    !=======================================================================
    ! 1.3 Initialise mesh objects and assign InterGrid maps
    !=======================================================================

    ! 1.3a Initialise prime/2d meshes
    ! ---------------------------------------------------------
    stencil_depth = get_required_stencil_depth()
    apply_partition_check = .false.

    call init_mesh( modeldb%configuration,        &
                    modeldb%mpi%get_comm_rank(),  &
                    modeldb%mpi%get_comm_size(),  &
                    base_mesh_names,              &
                    extrusion, stencil_depth,     &
                    apply_partition_check )

    call create_mesh( base_mesh_names, extrusion_2d, &
                      alt_name=twod_names )
    call assign_mesh_maps(twod_names)

    ! 1.3b Initialise shifted meshes
    ! ---------------------------------------------------------
    if (allocated(meshes_to_shift)) then
      if (size(meshes_to_shift) > 0) then

        allocate( shifted_names, source=meshes_to_shift )
        do i=1, size(shifted_names)
          shifted_names(i) = trim(shifted_names(i))//'_shifted'
        end do
        call create_mesh( meshes_to_shift,   &
                          extrusion_shifted, &
                          alt_name=shifted_names )
        call assign_mesh_maps(shifted_names)

      end if
    end if

    ! 1.3c Initialise double-level meshes
    ! ---------------------------------------------------------
    if (allocated(meshes_to_double)) then
      if (size(meshes_to_double) > 0) then

        allocate( double_names, source=meshes_to_double )
        do i=1, size(double_names)
          double_names(i) = trim(double_names(i))//'_double'
        end do
        call create_mesh( meshes_to_double, &
                          extrusion_double, &
                          alt_name=double_names )
        call assign_mesh_maps(double_names)

      end if
    end if

    !=======================================================================
    ! 2.0 Initialise FEM / Coordinates / Fields
    !=======================================================================

    ! FEM initialisation
    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()

    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    call create_runtime_constants()

    ! Set up transport runtime collection type
    ! Transport on only one horizontal local mesh
    mesh => mesh_collection%get_mesh(prime_mesh_name)

    ! Set transport metadata for primal mesh
    call name_transport_prerun_setup( num_base_meshes )

    ! Initialise prognostic variables
    call name_transport_init_fields_alg( mesh, wind, density, tracer_con )

    ! Initialise all transport-only control algorithm
    call name_transport_init( density, tracer_con )

    !=======================================================================
    ! 3.0 Initialise IO / Clock
    !=======================================================================

    ! I/O initialisation
    call init_io( xios_ctx,           &
                  prime_mesh_name,    &
                  modeldb,            &
                  chi_inventory,      &
                  panel_id_inventory )

    ! Call clock initial step before initial conditions output
    ! This ensures that lfric_initial.nc will be written out
    if (modeldb%clock%is_initialisation() .and. use_xios_io) then
      call modeldb%io_contexts%get_io_context(xios_ctx, io_context)
      call advance(io_context, modeldb%clock)
    end if

    !=======================================================================
    ! 4.0 Initial output
    !=======================================================================

    ! Output initial conditions
    if (modeldb%clock%is_initialisation() .and. write_diag) then

      call write_vector_diagnostic( 'u', wind, modeldb%clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'rho', density, modeldb%clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer_con', tracer_con, modeldb%clock, &
                                    mesh, nodal_output_on_w3 )

      height_w3 => get_height_fe(W3, mesh%get_id())
      height_wth => get_height_fe(Wtheta, mesh%get_id())
      call write_scalar_diagnostic( 'height_w3', height_w3, modeldb%clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'height_wth', height_wth, modeldb%clock, &
                                    mesh, nodal_output_on_w3 )

    end if

    if (allocated(base_mesh_names))  deallocate(base_mesh_names)
    if (allocated(meshes_to_shift))  deallocate(meshes_to_shift)
    if (allocated(meshes_to_double)) deallocate(meshes_to_double)

    if (allocated(extra_io_mesh_names)) deallocate(extra_io_mesh_names)
    nullify(chi_inventory, panel_id_inventory, mesh)

  end subroutine initialise_name_transport

  !============================================================================
  !> @brief Performs a time step of the name_transport app.
  !>
  subroutine step_name_transport( model_clock )

    use base_mesh_config_mod,     only: prime_mesh_name
    use io_config_mod,            only: diagnostic_frequency, &
                                        nodal_output_on_w3,   &
                                        write_diag
    use sci_field_minmax_alg_mod, only: log_field_minmax

    implicit none

    class(model_clock_type), intent(in) :: model_clock

    type(mesh_type), pointer :: mesh
    integer(tik)             :: id

    ! Get mesh
    mesh => mesh_collection%get_mesh(prime_mesh_name)

    ! Print min/max of fields before transport step
    if (transport_density) then
      call log_field_minmax( LOG_LEVEL_INFO, 'rho', density )
    end if
    call log_field_minmax( LOG_LEVEL_INFO, 'tracer_con', tracer_con )

    write(log_scratch_space, '("/", A, "\ ")') repeat('*', 76)
    call log_event( log_scratch_space, LOG_LEVEL_TRACE )
    write( log_scratch_space, '(A,I0)' ) &
      'Start of timestep ', model_clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    if ( LPROF ) call start_timing( id, 'name_transport_step' )

    ! Transport field
    call name_transport_step( model_clock, wind, tracer_con, &
                              density, transport_density )

    if ( LPROF ) call stop_timing( id, 'name_transport_step' )

    ! Print min/max of fields after transport step
    if (transport_density) then
      call log_field_minmax( LOG_LEVEL_INFO, 'rho', density )
    end if
    call log_field_minmax( LOG_LEVEL_INFO, 'tracer_con', tracer_con )

    write( log_scratch_space, &
           '(A,I0)' ) 'End of timestep ', model_clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write(log_scratch_space, '("\", A, "/ ")') repeat('*', 76)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Output wind, density and tracer values
    if ( (mod( model_clock%get_step(), diagnostic_frequency ) == 0) &
         .and. write_diag ) then
      call write_vector_diagnostic( 'u', wind, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer_con', tracer_con, model_clock, &
                                    mesh, nodal_output_on_w3 )
      if (transport_density) then
        call write_scalar_diagnostic( 'rho', density, model_clock, &
                                      mesh, nodal_output_on_w3 )
      end if
    end if

    nullify(mesh)

  end subroutine step_name_transport

  !============================================================================
  !> @brief Tidies up after a run.
  !>
  subroutine finalise_name_transport( program_name, modeldb )

    implicit none

    character(*),        intent(in)    :: program_name
    class(modeldb_type), intent(inout) :: modeldb

    call name_transport_final( density, tracer_con, transport_density )

    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------

    ! Write checksums to file
    call checksum_alg( program_name, density, 'rho',  wind, 'u',  &
                       tracer_con, 'tracer' )

    call final_io(modeldb)

  end subroutine finalise_name_transport

end module name_transport_driver_mod
