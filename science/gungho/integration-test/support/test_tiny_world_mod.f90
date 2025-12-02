!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module test_tiny_world_mod

  use constants_mod,       only : i_def, i_timestep, r_def, r_second
  use field_mod,           only : field_type
  use global_mesh_mod,     only : global_mesh_type
  use lfric_mpi_mod,       only : global_mpi
  use local_mesh_mod,      only : local_mesh_type
  use mesh_mod,            only : mesh_type
  use mesh_collection_mod, only : mesh_collection
  use model_clock_mod,     only : model_clock_type

  implicit none

  private
  public :: initialise_tiny_world, finalise_tiny_world

  type(field_type)                 :: chi(3)
  type(global_mesh_type), pointer  :: global_mesh
  type(local_mesh_type),  pointer  :: local_mesh
  type(mesh_type), public, pointer :: mesh
  type(model_clock_type)           :: model_clock
  type(field_type)                 :: panel_id

  real(r_def), parameter :: planet_radius = 600.0_r_def

contains

  subroutine initialise_tiny_world()

    use function_space_collection_mod, only : function_space_collection_type, &
                                              function_space_collection
    use halo_comms_mod,                only : initialise_halo_comms
    use halo_routing_collection_mod,   only : halo_routing_collection_type, &
                                              halo_routing_collection
    use lfric_mpi_mod,                 only : lfric_comm_type, &
                                              create_comm
    use mesh_collection_mod,           only : mesh_collection_type, &
                                              mesh_collection
    use runtime_tools_mod,             only : init_hierarchical_mesh_id_list

    implicit none

    type(lfric_comm_type) :: communicator

    integer :: mesh_id

    call create_comm( communicator )
    call global_mpi%initialise( communicator )
    call initialise_halo_comms( communicator )

    mesh_collection = mesh_collection_type()
    function_space_collection = function_space_collection_type()
    halo_routing_collection = halo_routing_collection_type()

    call build_topology( planet_radius, mesh_id )
    call build_geometry( mesh, chi, panel_id )

    model_clock = model_clock_type( first=1_i_timestep,            &
                                    last=3_i_timestep,             &
                                    seconds_per_step=1.0_r_second, &
                                    spinup_period=0.0_r_second )

    ! Gung Ho finite element bits
    !
    call init_hierarchical_mesh_id_list( (/mesh_id/) )

  end subroutine initialise_tiny_world


  subroutine build_topology(planet_radius, mesh_id )

    use extrusion_mod,           only : extrusion_type,         &
                                        uniform_extrusion_type, &
                                        prime_extrusion
    use panel_decomposition_mod, only : custom_decomposition_type
    use partition_mod,           only : partition_type,        &
                                        partitioner_interface, &
                                        partitioner_cubedsphere
    use ugrid_mesh_data_mod,     only : ugrid_mesh_data_type

    implicit none

    procedure(partitioner_interface), pointer :: partitioner_proc
    real(r_def), intent(in)  :: planet_radius
    integer,     intent(out) :: mesh_id

    type(ugrid_mesh_data_type)                   :: ugrid_data
    class(extrusion_type), allocatable           :: extrusion
    type(custom_decomposition_type), allocatable :: decomposition
    type(partition_type)                         :: partitioner

    call ugrid_data%read_from_file('shared-resources/tiny_world.nc', 'tiny')
    allocate( global_mesh, source=global_mesh_type( ugrid_data ) )
    decomposition = custom_decomposition_type( num_xprocs=1_i_def, num_yprocs=1_i_def )
    partitioner_proc => partitioner_cubedsphere
    partitioner = partition_type( global_mesh,                           &
                                  partitioner_proc,                      &
                                  decomposition,                         &
                                  max_stencil_depth=1_i_def,             &
                                  generate_inner_halos=.true.,           &
                                  local_rank=global_mpi%get_comm_rank(), &
                                  total_ranks=global_mpi%get_comm_size() )
    allocate( local_mesh )
    call local_mesh%initialise( global_mesh, partitioner )
    call local_mesh%init_cell_owner()
    allocate( extrusion,                                   &
              source=uniform_extrusion_type(               &
                atmosphere_bottom=planet_radius,           &
                atmosphere_top=planet_radius + 30.0_r_def, &
                number_of_layers=3,                        &
                extrusion_id=prime_extrusion               &
              ) )
    allocate( mesh, source=mesh_type( local_mesh, extrusion ) )
    mesh_id = mesh_collection%add_new_mesh( mesh )

    if (allocated(extrusion)) deallocate( extrusion )

  end subroutine build_topology


  subroutine build_geometry( mesh, chi, panel_id )

    use coord_transform_mod,           only : identify_panel, xyz2llr
    use field_mod,                     only : field_proxy_type
    use function_space_mod,            only : function_space_type
    use function_space_collection_mod, only : function_space_collection
    use fs_continuity_mod,             only : W3, W0
    use reference_element_mod,         only : reference_element_type

    implicit none

    type(mesh_type),   intent(in),    target  :: mesh
    class(field_type), intent(inout)          :: chi(3)
    class(field_type), intent(inout)          :: panel_id

    type(function_space_type), pointer :: chi_fs
    type(function_space_type), pointer :: w3_fs

    class(reference_element_type), pointer :: reference_element

    type(field_proxy_type) :: chi_proxy(3)
    type(field_proxy_type) :: panel_id_proxy

    real(r_def), allocatable :: column_coords(:,:,:)
    real(r_def), allocatable :: vertex_coords(:, :)
    real(r_def), allocatable :: vertex_local_coords(:, :)

    real(r_def), pointer :: chi_dof_coords(:, :)

    real(r_def) :: x, y, z
    real(r_def) :: interp_weight
    real(r_def) :: v_x, v_y, v_z
    real(r_def) :: v_lon, v_lat, v_r
    real(r_def) :: longitude, latitude, radius

    integer(i_def), pointer :: chi_map(:, :)
    integer(i_def), pointer :: panel_id_map(:, :)

    integer(i_def) :: chi_undf, chi_ndf, chi_nlayers
    integer(i_def) :: panel_id_undf, panel_id_ndf, panel_id_nlayers
    integer(i_def) :: nverts

    integer :: cell
    integer :: vert
    integer :: i, k
    integer :: df
    integer :: dfk
    integer :: panel

    chi_fs => function_space_collection%get_fs( mesh,              &
                                                element_order_h=0, &
                                                element_order_v=0, &
                                                lfric_fs=W0 )
    do i = 1, 3
      call chi(i)%initialise( chi_fs )
    end do

    w3_fs => function_space_collection%get_fs( mesh,              &
                                               element_order_h=0, &
                                               element_order_v=0, &
                                               lfric_fs=W3 )
    call panel_id%initialise( w3_fs )

    ! Break encapsulation and get the proxy.
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    chi_undf = chi_proxy(1)%vspace%get_undf()
    chi_ndf = chi_proxy(1)%vspace%get_ndf()
    chi_nlayers = chi_proxy(1)%vspace%get_nlayers()
    chi_map => chi_proxy(1)%vspace%get_whole_dofmap()
    chi_dof_coords => chi_proxy(1)%vspace%get_nodes()

    panel_id_proxy = panel_id%get_proxy()
    panel_id_undf = panel_id_proxy%vspace%get_undf()
    panel_id_ndf = panel_id_proxy%vspace%get_ndf()
    panel_id_nlayers = panel_id_proxy%vspace%get_nlayers()
    panel_id_map => panel_id_proxy%vspace%get_whole_dofmap()

    reference_element => mesh%get_reference_element()
    call reference_element%get_vertex_coordinates( vertex_coords )
    nverts = reference_element%get_number_vertices()

    allocate( column_coords(3, nverts, chi_nlayers) )
    allocate( vertex_local_coords(3, nverts) )

    do cell = 1, chi_proxy(1)%vspace%get_ncell()

      call mesh%get_column_coords( cell, column_coords )

      ! Calculate penel IDs
      !
      x = 0.0_r_def
      y = 0.0_r_def
      z = 0.0_r_def
      interp_weight = 1.0_r_def / nverts

      do vert = 1,nverts
        ! evaluate at cell centre
        x = x + interp_weight * column_coords(1,vert,1)
        y = y + interp_weight * column_coords(2,vert,1)
        z = z + interp_weight * column_coords(3,vert,1)
      end do

      panel = identify_panel(x, y, z)
      panel_id_proxy%data(panel_id_map(1, cell)                          &
                          :panel_id_map(1, cell) + panel_id_nlayers - 1) &
        = real(panel, r_def)

      ! Assign the co-ordinates
      !
      do k = 0, chi_nlayers - 1
        vertex_local_coords(:,:) = column_coords(:,:,k+1)

        do df = 1, chi_ndf
          dfk = chi_map(df, cell) + k
          ! Compute interpolation weights
          longitude  = 0.0_r_def
          latitude = 0.0_r_def
          radius = 0.0_r_def
          do vert = 1, nverts
            interp_weight =                                                    &
              (1.0_r_def - abs(vertex_coords(vert,1) - chi_dof_coords(1,df)))  &
              *(1.0_r_def - abs(vertex_coords(vert,2) - chi_dof_coords(2,df))) &
              *(1.0_r_def - abs(vertex_coords(vert,3) - chi_dof_coords(3,df)))
            v_x = vertex_local_coords(1,vert)
            v_y = vertex_local_coords(2,vert)
            v_z = vertex_local_coords(3,vert)

            call xyz2llr(v_x,v_y,v_z,v_lon,v_lat,v_r)
            longitude = longitude + interp_weight*v_lon
            latitude = latitude + interp_weight*v_lat
            radius = radius + interp_weight*v_r
          end do

          chi_proxy(1)%data(dfk) = longitude
          chi_proxy(2)%data(dfk) = latitude
          chi_proxy(3)%data(dfk) = radius - planet_radius

        end do
      end do

    end do

  end subroutine build_geometry


  subroutine finalise_tiny_world

    use halo_comms_mod,        only : finalise_halo_comms
    use lfric_mpi_mod,         only : destroy_comm
    use runtime_tools_mod,     only : final_hierarchical_mesh_id_list
    use sci_fem_constants_mod, only : final_fem_constants

    implicit none

    call final_fem_constants()
    call final_hierarchical_mesh_id_list
    call finalise_halo_comms()
    if (associated(mesh)) deallocate( mesh )
    if (associated(local_mesh)) deallocate( local_mesh )
    if (associated(global_mesh)) deallocate( global_mesh )
    call destroy_comm()

  end subroutine finalise_tiny_world

end module test_tiny_world_mod
