!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls infrastructure related information used by the model

module gravity_wave_infrastructure_mod

  use add_mesh_map_mod,           only : assign_mesh_maps
  use driver_modeldb_mod,         only : modeldb_type
  use check_configuration_mod,    only : get_required_stencil_depth
  use constants_mod,              only : i_def,           &
                                         PRECISION_REAL,  &
                                         r_def, r_second, &
                                         l_def, str_def
  use convert_to_upper_mod,       only : convert_to_upper
  use create_mesh_mod,            only : create_mesh, create_extrusion
  use derived_config_mod,         only : set_derived_config
  use extrusion_mod,              only : extrusion_type,         &
                                         uniform_extrusion_type, &
                                         TWOD, PRIME_EXTRUSION
  use sci_geometric_constants_mod,        &
                                  only : get_chi_inventory,  &
                                         get_panel_id_inventory
  use inventory_by_mesh_mod,      only : inventory_by_mesh_type
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_ERROR
  use mesh_collection_mod,        only : mesh_collection
  use namelist_mod,               only : namelist_type
  use field_mod,                  only : field_type
  use driver_fem_mod,             only : init_fem, init_function_space_chains
  use driver_io_mod,              only : init_io, final_io
  use driver_mesh_mod,            only : init_mesh
  use runtime_constants_mod,      only : create_runtime_constants
  use remove_duplicates_mod,      only : remove_duplicates

  ! Configuration modules
  use base_mesh_config_mod,   only: GEOMETRY_PLANAR, &
                                    GEOMETRY_SPHERICAL

  implicit none

  private
  public initialise_infrastructure, finalise_infrastructure

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Initialises the infrastructure used by the model
  !> @param [in]      program_name  The name of the program being run
  !> @param [in,out]  modeldb       The modeldb object
  subroutine initialise_infrastructure( program_name,  &
                                        modeldb)

    implicit none

    character(*),       intent(in)    :: program_name
    type(modeldb_type), intent(inout) :: modeldb

    type(inventory_by_mesh_type), pointer :: chi_inventory => null()
    type(inventory_by_mesh_type), pointer :: panel_id_inventory => null()

    character(str_def), allocatable :: base_mesh_names(:)
    character(str_def), allocatable :: chain_mesh_tags(:)
    character(str_def), allocatable :: twod_names(:)

    character(str_def), allocatable :: tmp_mesh_names(:)

    class(extrusion_type),         allocatable :: extrusion
    type(uniform_extrusion_type),  allocatable :: extrusion_2d

    integer(i_def) :: start_index, end_index

    character(str_def) :: prime_mesh_name

    logical(l_def) :: l_multigrid
    logical(l_def) :: prepartitioned
    logical        :: apply_partition_check

    integer(i_def) :: stencil_depth
    integer(i_def) :: geometry
    integer(i_def) :: method
    integer(i_def) :: number_of_layers
    real(r_def)    :: domain_bottom
    real(r_def)    :: domain_height
    real(r_def)    :: scaled_radius

    type(namelist_type), pointer :: base_mesh_nml   => null()
    type(namelist_type), pointer :: formulation_nml => null()
    type(namelist_type), pointer :: extrusion_nml   => null()
    type(namelist_type), pointer :: planet_nml      => null()
    type(namelist_type), pointer :: multigrid_nml   => null()

    integer(i_def) :: i
    integer(i_def), parameter :: one_layer = 1_i_def

    !=======================================================================
    ! 0.0 Extract configuration variables
    !=======================================================================
    base_mesh_nml   => modeldb%configuration%get_namelist('base_mesh')
    formulation_nml => modeldb%configuration%get_namelist('formulation')
    extrusion_nml   => modeldb%configuration%get_namelist('extrusion')
    planet_nml      => modeldb%configuration%get_namelist('planet')

    call formulation_nml%get_value( 'l_multigrid', l_multigrid )

    if (l_multigrid) then
      multigrid_nml => modeldb%configuration%get_namelist('multigrid')
      call multigrid_nml%get_value( 'chain_mesh_tags', chain_mesh_tags )
      multigrid_nml => null()
    end if

    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    call base_mesh_nml%get_value( 'geometry', geometry )
    call base_mesh_nml%get_value( 'prepartitioned', prepartitioned )
    call extrusion_nml%get_value( 'method', method )
    call extrusion_nml%get_value( 'domain_height', domain_height )
    call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
    call planet_nml%get_value( 'scaled_radius', scaled_radius )

    base_mesh_nml   => null()
    extrusion_nml   => null()
    formulation_nml => null()
    planet_nml      => null()

    !-------------------------------------------------------------------------
    ! Initialise infrastructure
    !-------------------------------------------------------------------------
    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    call set_derived_config( .false. )


    !=======================================================================
    ! 1.0 Mesh
    !=======================================================================

    !=======================================================================
    ! 1.1 Determine the required meshes
    !=======================================================================
    if ( allocated(tmp_mesh_names) ) deallocate(tmp_mesh_names)

    allocate( tmp_mesh_names(100) )
    tmp_mesh_names(:) = ''
    tmp_mesh_names(1) = prime_mesh_name
    start_index = 2

    if ( l_multigrid ) then
      end_index = start_index + SIZE(chain_mesh_tags) - 1
      tmp_mesh_names(start_index:end_index) = chain_mesh_tags
      start_index = end_index + 1
    end if

    base_mesh_names = remove_duplicates( tmp_mesh_names(:)  )

    deallocate(tmp_mesh_names)


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
      call log_event("Invalid geometry for mesh initialisation", LOG_LEVEL_ERROR)
    end select
    allocate( extrusion, source=create_extrusion( method,           &
                                                  domain_height,       &
                                                  domain_bottom,    &
                                                  number_of_layers, &
                                                  PRIME_EXTRUSION ) )

    extrusion_2d = uniform_extrusion_type( domain_height,    &
                                           domain_bottom, &
                                           one_layer, TWOD )


    !=======================================================================
    ! 1.3 Initialise mesh objects and assign InterGrid maps
    !=======================================================================
    stencil_depth = get_required_stencil_depth()
    apply_partition_check = .false.
    if ( .not. prepartitioned .and. l_multigrid ) then
      apply_partition_check = .true.
    end if

    call init_mesh( modeldb%configuration,        &
                    modeldb%mpi%get_comm_rank(),  &
                    modeldb%mpi%get_comm_size(),  &
                    base_mesh_names,              &
                    extrusion, stencil_depth,     &
                    apply_partition_check )

    allocate( twod_names, source=base_mesh_names )
    do i=1, size(twod_names)
      twod_names(i) = trim(twod_names(i))//'_2d'
    end do
    call create_mesh( base_mesh_names, extrusion_2d, &
                      alt_name=twod_names )
    call assign_mesh_maps(twod_names)

    !=======================================================================
    ! 2.0 Create FEM specifics (function spaces and chi field)
    !=======================================================================
    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )
    if ( l_multigrid ) then
      call init_function_space_chains( mesh_collection, chain_mesh_tags )
    end if

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------
    call init_io( program_name, prime_mesh_name, &
                  modeldb, chi_inventory, panel_id_inventory)

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    call create_runtime_constants()


    nullify(chi_inventory, panel_id_inventory)
    deallocate(base_mesh_names)

  end subroutine initialise_infrastructure


  !> @brief Finalises infrastructure used by the model
  !> @param[inout] modeldb The model database object
  subroutine finalise_infrastructure(modeldb)

    implicit none

    class(modeldb_type), intent(inout) :: modeldb

    ! Finalise the IO and destroy the IO contexts
    call final_io(modeldb)

  end subroutine finalise_infrastructure

end module gravity_wave_infrastructure_mod
