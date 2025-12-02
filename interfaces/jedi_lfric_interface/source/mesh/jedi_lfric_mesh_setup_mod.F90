!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Provides a method to setup a  mesh for a JEDI-LFRIC
!>
module jedi_lfric_mesh_setup_mod

  use add_mesh_map_mod,        only: assign_mesh_maps
  use base_mesh_config_mod,    only: GEOMETRY_SPHERICAL, &
                                     GEOMETRY_PLANAR
  use constants_mod,           only: str_def, i_def, l_def, r_def
  use create_mesh_mod,         only: create_mesh
  use driver_mesh_mod,         only: init_mesh
  use extrusion_mod,           only: extrusion_type,         &
                                     uniform_extrusion_type, &
                                     TWOD
  use gungho_extrusion_mod,    only: create_extrusion
  use lfric_mpi_mod,           only: lfric_mpi_type
  use log_mod,                 only: log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_ERROR
  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type

  implicit none

  private

  public initialise_mesh

contains

  !> @brief Initialise the mesh and store it in the global mesh collection
  !>
  !> @param [out]   mesh_name     The name of the mesh being setup
  !> @param [in]    configuration The geometry configuration
  !> @param [inout] mpi_obj       The mpi communicator
  !> @param [in]    alt_mesh_name The name of an alternative mesh_name to setup
  subroutine initialise_mesh( mesh_name, configuration, mpi_obj, alt_mesh_name )

    implicit none

    character(len=*),              intent(out) :: mesh_name
    type(namelist_collection_type), intent(in) :: configuration
    !> @todo: This should be intent in but when calling the method I get
    !> a compiler failure
    class(lfric_mpi_type),       intent(inout) :: mpi_obj
    character(len=*),     optional, intent(in) :: alt_mesh_name

    ! Local
    type(namelist_type),              pointer :: base_mesh_nml
    type(namelist_type),              pointer :: extrusion_nml
    type(namelist_type),              pointer :: planet_nml
    class(extrusion_type),        allocatable :: extrusion
    type(uniform_extrusion_type), allocatable :: extrusion_2d

    character(str_def), allocatable :: twod_names(:)
    character(str_def)              :: base_mesh_names(1)
    character(str_def)              :: prime_mesh_name
    integer(i_def),       parameter :: one_layer = 1_i_def
    integer(i_def)                  :: geometry
    integer(i_def)                  :: extrusion_method
    integer(i_def)                  :: stencil_depth
    integer(i_def)                  :: number_of_layers
    integer(i_def)                  :: i
    real(r_def)                     :: domain_bottom
    real(r_def)                     :: domain_height
    real(r_def)                     :: scaled_radius
    logical(l_def)                  :: apply_partition_check

    !--------------------------------------
    ! 0.0 Extract namelist variables
    !--------------------------------------
    nullify(base_mesh_nml, planet_nml, extrusion_nml)
    base_mesh_nml => configuration%get_namelist('base_mesh')
    planet_nml    => configuration%get_namelist('planet')
    extrusion_nml => configuration%get_namelist('extrusion')

    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    call base_mesh_nml%get_value( 'geometry', geometry )
    call extrusion_nml%get_value( 'domain_height', domain_height )
    call extrusion_nml%get_value( 'method', extrusion_method )
    call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
    call planet_nml%get_value( 'scaled_radius', scaled_radius )

    !--------------------------------------
    ! 1.0 Create the meshes
    !--------------------------------------
    ! Choose mesh to use
    if ( present(alt_mesh_name) ) then
      !@todo: this is required until we can create multiple configurations
      base_mesh_names(1) = alt_mesh_name
    else
      base_mesh_names(1) = prime_mesh_name
    endif
    mesh_name = base_mesh_names(1)

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

    allocate( extrusion, source=create_extrusion( extrusion_method, &
                                                  geometry,         &
                                                  number_of_layers, &
                                                  domain_height,    &
                                                  scaled_radius ) )

    extrusion_2d = uniform_extrusion_type( domain_height, &
                                           domain_bottom, &
                                           one_layer, TWOD )

    !-------------------------------------------------------------------------
    ! 1.2 Create the required meshes
    !-------------------------------------------------------------------------
    stencil_depth = 2
    apply_partition_check = .false.
    call init_mesh( configuration,           &
                    mpi_obj%get_comm_rank(), &
                    mpi_obj%get_comm_size(), &
                    base_mesh_names,         &
                    extrusion,               &
                    stencil_depth,           &
                    apply_partition_check )

    allocate( twod_names, source=base_mesh_names )
    do i=1, size(twod_names)
      twod_names(i) = trim(twod_names(i))//'_2d'
    end do
    call create_mesh( base_mesh_names, extrusion_2d, &
                      alt_name=twod_names )
    call assign_mesh_maps(twod_names)

  end subroutine initialise_mesh

end module jedi_lfric_mesh_setup_mod
