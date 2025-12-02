!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief PSy-lite code required for invokation of the lookup table generation kernels.
!> @detail Halo exchange depths have been manually changed for `lookup_X_field` and
!!         `set_counts_X_field` to halo exchange to the required depth. This is
!!         required because these fields are accessed with a STENCIL, but are
!!         not marked as such. The reason being that PSyclone does not allow
!!         writing to STENCIL fields. We are working around this
!!         by passing a second set of "dummy" fields in READ mode with a STENCIL.
!!         See PSyclone issue #3003.
module psykal_lite_gen_lookup_tables_psy_mod
  use constants_mod, only: r_tran, r_solver, i_def, r_def
  use r_tran_field_mod, only: r_tran_field_type, r_tran_field_proxy_type
  use r_solver_field_mod, only: r_solver_field_type, r_solver_field_proxy_type
  use integer_field_mod, only: integer_field_type, integer_field_proxy_type
  use operator_mod, only: operator_type, operator_proxy_type
  implicit none
  contains

  subroutine invoke_gen_poly1d_lookup_kernel(tracer, coeff, reconstruction, lookup_poly1d_field, set_counts_poly1d_field, &
&lookup_poly1d_field_dummy, set_counts_poly1d_field_dummy, order, nindices, stencil_extent, loop_halo_depth)
    use gen_poly1d_lookup_kernel_mod, only: gen_poly1d_lookup_code
    use mesh_mod, only: mesh_type
    use stencil_dofmap_mod, only: STENCIL_CROSS
    use stencil_dofmap_mod, only: stencil_dofmap_type
    integer(kind=i_def), intent(in) :: order, nindices
    type(r_tran_field_type), intent(in) :: tracer, coeff, reconstruction
    type(integer_field_type), intent(in) :: lookup_poly1d_field, set_counts_poly1d_field, lookup_poly1d_field_dummy, &
&set_counts_poly1d_field_dummy
    integer(kind=i_def), intent(in) :: stencil_extent
    integer, intent(in) :: loop_halo_depth
    integer(kind=i_def) cell
    integer(kind=i_def) loop0_start, loop0_stop
    integer(kind=i_def) nlayers_tracer
    integer(kind=i_def), pointer, dimension(:) :: set_counts_poly1d_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_poly1d_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: set_counts_poly1d_field_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_poly1d_field_data
    type(integer_field_proxy_type) lookup_poly1d_field_proxy, set_counts_poly1d_field_proxy, lookup_poly1d_field_dummy_proxy, &
&set_counts_poly1d_field_dummy_proxy
    real(kind=r_tran), pointer, dimension(:) :: reconstruction_data
    real(kind=r_tran), pointer, dimension(:) :: coeff_data
    real(kind=r_tran), pointer, dimension(:) :: tracer_data
    type(r_tran_field_proxy_type) tracer_proxy, coeff_proxy, reconstruction_proxy
    integer(kind=i_def), pointer :: map_adspc1_tracer(:,:), map_adspc2_coeff(:,:), &
&map_adspc3_reconstruction(:,:), map_adspc4_lookup_poly1d_field(:,:), &
&map_adspc5_set_counts_poly1d_field(:,:)
    integer(kind=i_def) ndf_adspc1_tracer, undf_adspc1_tracer, ndf_adspc2_coeff, undf_adspc2_coeff, ndf_adspc3_reconstruction, &
&undf_adspc3_reconstruction, ndf_adspc4_lookup_poly1d_field, undf_adspc4_lookup_poly1d_field, ndf_adspc5_set_counts_poly1d_field, &
&undf_adspc5_set_counts_poly1d_field
    integer(kind=i_def) max_halo_depth_mesh
    type(mesh_type), pointer :: mesh
    integer(kind=i_def), pointer :: set_counts_poly1d_field_dummy_stencil_size(:)
    integer(kind=i_def), pointer :: set_counts_poly1d_field_dummy_stencil_dofmap(:,:,:)
    type(stencil_dofmap_type), pointer :: set_counts_poly1d_field_dummy_stencil_map
    integer(kind=i_def), pointer :: lookup_poly1d_field_dummy_stencil_size(:)
    integer(kind=i_def), pointer :: lookup_poly1d_field_dummy_stencil_dofmap(:,:,:)
    type(stencil_dofmap_type), pointer :: lookup_poly1d_field_dummy_stencil_map

    nullify( set_counts_poly1d_field_dummy_data, lookup_poly1d_field_dummy_data, set_counts_poly1d_field_data, &
             lookup_poly1d_field_data, reconstruction_data, coeff_data, tracer_data, &
             map_adspc1_tracer, map_adspc2_coeff, map_adspc4_lookup_poly1d_field, &
             mesh, set_counts_poly1d_field_dummy_stencil_size, set_counts_poly1d_field_dummy_stencil_dofmap, &
             set_counts_poly1d_field_dummy_stencil_map, lookup_poly1d_field_dummy_stencil_size, &
             lookup_poly1d_field_dummy_stencil_dofmap, lookup_poly1d_field_dummy_stencil_map)

    !
    ! Initialise field and/or operator proxies
    !
    tracer_proxy = tracer%get_proxy()
    tracer_data => tracer_proxy%data
    coeff_proxy = coeff%get_proxy()
    coeff_data => coeff_proxy%data
    reconstruction_proxy = reconstruction%get_proxy()
    reconstruction_data => reconstruction_proxy%data
    lookup_poly1d_field_proxy = lookup_poly1d_field%get_proxy()
    lookup_poly1d_field_data => lookup_poly1d_field_proxy%data
    set_counts_poly1d_field_proxy = set_counts_poly1d_field%get_proxy()
    set_counts_poly1d_field_data => set_counts_poly1d_field_proxy%data
    lookup_poly1d_field_dummy_proxy = lookup_poly1d_field_dummy%get_proxy()
    lookup_poly1d_field_dummy_data => lookup_poly1d_field_dummy_proxy%data
    set_counts_poly1d_field_dummy_proxy = set_counts_poly1d_field_dummy%get_proxy()
    set_counts_poly1d_field_dummy_data => set_counts_poly1d_field_dummy_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_tracer = tracer_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => tracer_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Initialise stencil dofmaps
    !
    lookup_poly1d_field_dummy_stencil_map => &
&lookup_poly1d_field_dummy_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS,stencil_extent)
    lookup_poly1d_field_dummy_stencil_dofmap => lookup_poly1d_field_dummy_stencil_map%get_whole_dofmap()
    lookup_poly1d_field_dummy_stencil_size => lookup_poly1d_field_dummy_stencil_map%get_stencil_sizes()
    set_counts_poly1d_field_dummy_stencil_map => &
&set_counts_poly1d_field_dummy_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS,stencil_extent)
    set_counts_poly1d_field_dummy_stencil_dofmap => set_counts_poly1d_field_dummy_stencil_map%get_whole_dofmap()
    set_counts_poly1d_field_dummy_stencil_size => set_counts_poly1d_field_dummy_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_tracer => tracer_proxy%vspace%get_whole_dofmap()
    map_adspc2_coeff => coeff_proxy%vspace%get_whole_dofmap()
    map_adspc3_reconstruction => reconstruction_proxy%vspace%get_whole_dofmap()
    map_adspc4_lookup_poly1d_field => lookup_poly1d_field_proxy%vspace%get_whole_dofmap()
    map_adspc5_set_counts_poly1d_field => set_counts_poly1d_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for adspc1_tracer
    !
    ndf_adspc1_tracer = tracer_proxy%vspace%get_ndf()
    undf_adspc1_tracer = tracer_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_coeff
    !
    ndf_adspc2_coeff = coeff_proxy%vspace%get_ndf()
    undf_adspc2_coeff = coeff_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc3_reconstruction
    !
    ndf_adspc3_reconstruction = reconstruction_proxy%vspace%get_ndf()
    undf_adspc3_reconstruction = reconstruction_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc4_lookup_poly1d_field
    !
    ndf_adspc4_lookup_poly1d_field = lookup_poly1d_field_proxy%vspace%get_ndf()
    undf_adspc4_lookup_poly1d_field = lookup_poly1d_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc5_set_counts_poly1d_field
    !
    ndf_adspc5_set_counts_poly1d_field = set_counts_poly1d_field_proxy%vspace%get_ndf()
    undf_adspc5_set_counts_poly1d_field = set_counts_poly1d_field_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_halo_cell(loop_halo_depth)
    !
    ! Call kernels and communication routines
    !
    if (tracer_proxy%is_dirty(depth=loop_halo_depth)) then
      call tracer_proxy%halo_exchange(depth=loop_halo_depth)
    end if
    if (coeff_proxy%is_dirty(depth=loop_halo_depth)) then
      call coeff_proxy%halo_exchange(depth=loop_halo_depth)
    end if
    if (reconstruction_proxy%is_dirty(depth=loop_halo_depth)) then
      call reconstruction_proxy%halo_exchange(depth=loop_halo_depth)
    end if
    if (lookup_poly1d_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call lookup_poly1d_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    if (set_counts_poly1d_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call set_counts_poly1d_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    do cell = loop0_start, loop0_stop, 1
      call gen_poly1d_lookup_code(nlayers_tracer, tracer_data, coeff_data, reconstruction_data, lookup_poly1d_field_data, &
&set_counts_poly1d_field_data, lookup_poly1d_field_dummy_data, lookup_poly1d_field_dummy_stencil_size(cell), &
&lookup_poly1d_field_dummy_stencil_dofmap(:,:,cell), set_counts_poly1d_field_dummy_data, &
&set_counts_poly1d_field_dummy_stencil_size(cell), set_counts_poly1d_field_dummy_stencil_dofmap(:,:,cell), order, nindices, &
&ndf_adspc1_tracer, undf_adspc1_tracer, map_adspc1_tracer(:,cell), ndf_adspc2_coeff, undf_adspc2_coeff, map_adspc2_coeff(:,cell), &
&ndf_adspc3_reconstruction, undf_adspc3_reconstruction, map_adspc3_reconstruction(:,cell), ndf_adspc4_lookup_poly1d_field, &
&undf_adspc4_lookup_poly1d_field, map_adspc4_lookup_poly1d_field(:,cell), ndf_adspc5_set_counts_poly1d_field, &
&undf_adspc5_set_counts_poly1d_field, map_adspc5_set_counts_poly1d_field(:,cell))
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call lookup_poly1d_field_proxy%set_dirty()
    call lookup_poly1d_field_proxy%set_clean(loop_halo_depth)
    call set_counts_poly1d_field_proxy%set_dirty()
    call set_counts_poly1d_field_proxy%set_clean(loop_halo_depth)
    !
    !
  end subroutine invoke_gen_poly1d_lookup_kernel

  subroutine invoke_gen_poly2d_lookup_kernel(tracer, coeff, reconstruction, lookup_poly2d_field, set_counts_poly2d_field, &
&lookup_poly2d_field_dummy, set_counts_poly2d_field_dummy, stencil_size, nindices, stencil_extent, loop_halo_depth)
    use gen_poly2d_lookup_kernel_mod, only: gen_poly2d_lookup_code
    use mesh_mod, only: mesh_type
    use stencil_dofmap_mod, only: STENCIL_REGION
    use stencil_dofmap_mod, only: stencil_dofmap_type
    integer(kind=i_def), intent(in) :: stencil_size, nindices
    type(r_tran_field_type), intent(in) :: tracer, coeff, reconstruction
    type(integer_field_type), intent(in) :: lookup_poly2d_field, set_counts_poly2d_field, lookup_poly2d_field_dummy, &
&set_counts_poly2d_field_dummy
    integer(kind=i_def), intent(in) :: stencil_extent
    integer, intent(in) :: loop_halo_depth
    integer(kind=i_def) cell
    integer(kind=i_def) loop0_start, loop0_stop
    integer(kind=i_def) nlayers_tracer
    integer(kind=i_def), pointer, dimension(:) :: set_counts_poly2d_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_poly2d_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: set_counts_poly2d_field_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_poly2d_field_data
    type(integer_field_proxy_type) lookup_poly2d_field_proxy, set_counts_poly2d_field_proxy, lookup_poly2d_field_dummy_proxy, &
&set_counts_poly2d_field_dummy_proxy
    real(kind=r_tran), pointer, dimension(:) :: reconstruction_data
    real(kind=r_tran), pointer, dimension(:) :: coeff_data
    real(kind=r_tran), pointer, dimension(:) :: tracer_data
    type(r_tran_field_proxy_type) tracer_proxy, coeff_proxy, reconstruction_proxy
    integer(kind=i_def), pointer :: map_adspc1_tracer(:,:), map_adspc2_coeff(:,:), &
&map_adspc3_reconstruction(:,:), map_adspc4_lookup_poly2d_field(:,:), &
&map_adspc5_set_counts_poly2d_field(:,:)
    integer(kind=i_def) ndf_adspc1_tracer, undf_adspc1_tracer, ndf_adspc2_coeff, undf_adspc2_coeff, ndf_adspc3_reconstruction, &
&undf_adspc3_reconstruction, ndf_adspc4_lookup_poly2d_field, undf_adspc4_lookup_poly2d_field, ndf_adspc5_set_counts_poly2d_field, &
&undf_adspc5_set_counts_poly2d_field
    integer(kind=i_def) max_halo_depth_mesh
    type(mesh_type), pointer :: mesh
    integer(kind=i_def), pointer :: set_counts_poly2d_field_dummy_stencil_size(:)
    integer(kind=i_def), pointer :: set_counts_poly2d_field_dummy_stencil_dofmap(:,:,:)
    type(stencil_dofmap_type), pointer :: set_counts_poly2d_field_dummy_stencil_map
    integer(kind=i_def), pointer :: lookup_poly2d_field_dummy_stencil_size(:)
    integer(kind=i_def), pointer :: lookup_poly2d_field_dummy_stencil_dofmap(:,:,:)
    type(stencil_dofmap_type), pointer :: lookup_poly2d_field_dummy_stencil_map

    nullify( set_counts_poly2d_field_dummy_data, lookup_poly2d_field_dummy_data, set_counts_poly2d_field_data, &
             lookup_poly2d_field_data, reconstruction_data, coeff_data, tracer_data, &
             map_adspc1_tracer, map_adspc2_coeff, map_adspc4_lookup_poly2d_field, &
             mesh, set_counts_poly2d_field_dummy_stencil_size, set_counts_poly2d_field_dummy_stencil_dofmap, &
             set_counts_poly2d_field_dummy_stencil_map, lookup_poly2d_field_dummy_stencil_size, &
             lookup_poly2d_field_dummy_stencil_dofmap, lookup_poly2d_field_dummy_stencil_map)

    !
    ! Initialise field and/or operator proxies
    !
    tracer_proxy = tracer%get_proxy()
    tracer_data => tracer_proxy%data
    coeff_proxy = coeff%get_proxy()
    coeff_data => coeff_proxy%data
    reconstruction_proxy = reconstruction%get_proxy()
    reconstruction_data => reconstruction_proxy%data
    lookup_poly2d_field_proxy = lookup_poly2d_field%get_proxy()
    lookup_poly2d_field_data => lookup_poly2d_field_proxy%data
    set_counts_poly2d_field_proxy = set_counts_poly2d_field%get_proxy()
    set_counts_poly2d_field_data => set_counts_poly2d_field_proxy%data
    lookup_poly2d_field_dummy_proxy = lookup_poly2d_field_dummy%get_proxy()
    lookup_poly2d_field_dummy_data => lookup_poly2d_field_dummy_proxy%data
    set_counts_poly2d_field_dummy_proxy = set_counts_poly2d_field_dummy%get_proxy()
    set_counts_poly2d_field_dummy_data => set_counts_poly2d_field_dummy_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_tracer = tracer_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => tracer_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Initialise stencil dofmaps
    !
    lookup_poly2d_field_dummy_stencil_map => &
&lookup_poly2d_field_dummy_proxy%vspace%get_stencil_dofmap(STENCIL_REGION,stencil_extent)
    lookup_poly2d_field_dummy_stencil_dofmap => lookup_poly2d_field_dummy_stencil_map%get_whole_dofmap()
    lookup_poly2d_field_dummy_stencil_size => lookup_poly2d_field_dummy_stencil_map%get_stencil_sizes()
    set_counts_poly2d_field_dummy_stencil_map => &
&set_counts_poly2d_field_dummy_proxy%vspace%get_stencil_dofmap(STENCIL_REGION,stencil_extent)
    set_counts_poly2d_field_dummy_stencil_dofmap => set_counts_poly2d_field_dummy_stencil_map%get_whole_dofmap()
    set_counts_poly2d_field_dummy_stencil_size => set_counts_poly2d_field_dummy_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_tracer => tracer_proxy%vspace%get_whole_dofmap()
    map_adspc2_coeff => coeff_proxy%vspace%get_whole_dofmap()
    map_adspc3_reconstruction => reconstruction_proxy%vspace%get_whole_dofmap()
    map_adspc4_lookup_poly2d_field => lookup_poly2d_field_proxy%vspace%get_whole_dofmap()
    map_adspc5_set_counts_poly2d_field => set_counts_poly2d_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for adspc1_tracer
    !
    ndf_adspc1_tracer = tracer_proxy%vspace%get_ndf()
    undf_adspc1_tracer = tracer_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_coeff
    !
    ndf_adspc2_coeff = coeff_proxy%vspace%get_ndf()
    undf_adspc2_coeff = coeff_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc3_reconstruction
    !
    ndf_adspc3_reconstruction = reconstruction_proxy%vspace%get_ndf()
    undf_adspc3_reconstruction = reconstruction_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc4_lookup_poly2d_field
    !
    ndf_adspc4_lookup_poly2d_field = lookup_poly2d_field_proxy%vspace%get_ndf()
    undf_adspc4_lookup_poly2d_field = lookup_poly2d_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc5_set_counts_poly2d_field
    !
    ndf_adspc5_set_counts_poly2d_field = set_counts_poly2d_field_proxy%vspace%get_ndf()
    undf_adspc5_set_counts_poly2d_field = set_counts_poly2d_field_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_halo_cell(loop_halo_depth)
    !
    ! Call kernels and communication routines
    !
    if (tracer_proxy%is_dirty(depth=loop_halo_depth)) then
      call tracer_proxy%halo_exchange(depth=loop_halo_depth)
    end if
    if (coeff_proxy%is_dirty(depth=loop_halo_depth)) then
      call coeff_proxy%halo_exchange(depth=loop_halo_depth)
    end if
    if (reconstruction_proxy%is_dirty(depth=loop_halo_depth)) then
      call reconstruction_proxy%halo_exchange(depth=loop_halo_depth)
    end if
    if (lookup_poly2d_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call lookup_poly2d_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    if (set_counts_poly2d_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call set_counts_poly2d_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    do cell = loop0_start, loop0_stop, 1
      call gen_poly2d_lookup_code(nlayers_tracer, tracer_data, coeff_data, reconstruction_data, lookup_poly2d_field_data, &
&set_counts_poly2d_field_data, lookup_poly2d_field_dummy_data, lookup_poly2d_field_dummy_stencil_size(cell), &
&lookup_poly2d_field_dummy_stencil_dofmap(:,:,cell), set_counts_poly2d_field_dummy_data, &
&set_counts_poly2d_field_dummy_stencil_size(cell), set_counts_poly2d_field_dummy_stencil_dofmap(:,:,cell), stencil_size, nindices, &
&ndf_adspc1_tracer, undf_adspc1_tracer, map_adspc1_tracer(:,cell), ndf_adspc2_coeff, undf_adspc2_coeff, map_adspc2_coeff(:,cell), &
&ndf_adspc3_reconstruction, undf_adspc3_reconstruction, map_adspc3_reconstruction(:,cell), ndf_adspc4_lookup_poly2d_field, &
&undf_adspc4_lookup_poly2d_field, map_adspc4_lookup_poly2d_field(:,cell), ndf_adspc5_set_counts_poly2d_field, &
&undf_adspc5_set_counts_poly2d_field, map_adspc5_set_counts_poly2d_field(:,cell))
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call lookup_poly2d_field_proxy%set_dirty()
    call lookup_poly2d_field_proxy%set_clean(loop_halo_depth)
    call set_counts_poly2d_field_proxy%set_dirty()
    call set_counts_poly2d_field_proxy%set_clean(loop_halo_depth)
    !
    !
  end subroutine invoke_gen_poly2d_lookup_kernel

  subroutine invoke_gen_poly_adv_upd_lookup_kernel(advective, reconstruction_big_halo, wind_big_halo, lookup_field, &
&set_counts_field, lookup_field_dummy, set_counts_field_dummy, nsets_max, nindices, stencil_extent, loop_halo_depth)
    use gen_poly_adv_upd_lookup_kernel_mod, only: gen_poly_adv_upd_lookup_code
    use mesh_mod, only: mesh_type
    use stencil_2D_dofmap_mod, only: stencil_2D_dofmap_type, STENCIL_2D_CROSS
    integer(kind=i_def), intent(in) :: nsets_max, nindices
    type(r_tran_field_type), intent(in) :: advective, reconstruction_big_halo, wind_big_halo
    type(integer_field_type), intent(in) :: lookup_field, set_counts_field, lookup_field_dummy, set_counts_field_dummy
    integer(kind=i_def), intent(in) :: stencil_extent
    integer, intent(in) :: loop_halo_depth
    integer(kind=i_def) cell
    integer(kind=i_def) loop0_start, loop0_stop
    integer(kind=i_def) nlayers_advective
    integer(kind=i_def), pointer, dimension(:) :: set_counts_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: set_counts_field_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_field_data
    type(integer_field_proxy_type) lookup_field_proxy, set_counts_field_proxy, lookup_field_dummy_proxy, &
&set_counts_field_dummy_proxy
    real(kind=r_tran), pointer, dimension(:) :: wind_big_halo_data
    real(kind=r_tran), pointer, dimension(:) :: reconstruction_big_halo_data
    real(kind=r_tran), pointer, dimension(:) :: advective_data
    type(r_tran_field_proxy_type) advective_proxy, reconstruction_big_halo_proxy, wind_big_halo_proxy
    integer(kind=i_def), pointer :: map_adspc1_reconstruction_big_halo(:,:), map_adspc2_lookup_field(:,:), &
&map_adspc3_set_counts_field(:,:), map_w2(:,:), map_wtheta(:,:)
    integer(kind=i_def) ndf_wtheta, undf_wtheta, ndf_adspc1_reconstruction_big_halo, undf_adspc1_reconstruction_big_halo, &
&ndf_w2, undf_w2, ndf_adspc2_lookup_field, undf_adspc2_lookup_field, ndf_adspc3_set_counts_field, undf_adspc3_set_counts_field
    integer(kind=i_def) max_halo_depth_mesh
    type(mesh_type), pointer :: mesh
    integer(kind=i_def) set_counts_field_dummy_max_branch_length
    integer(kind=i_def), pointer :: set_counts_field_dummy_stencil_size(:,:)
    integer(kind=i_def), pointer :: set_counts_field_dummy_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: set_counts_field_dummy_stencil_map
    integer(kind=i_def) lookup_field_dummy_max_branch_length
    integer(kind=i_def), pointer :: lookup_field_dummy_stencil_size(:,:)
    integer(kind=i_def), pointer :: lookup_field_dummy_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: lookup_field_dummy_stencil_map
    integer(kind=i_def) wind_big_halo_max_branch_length
    integer(kind=i_def), pointer :: wind_big_halo_stencil_size(:,:)
    integer(kind=i_def), pointer :: wind_big_halo_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: wind_big_halo_stencil_map
    integer(kind=i_def) reconstruction_big_halo_max_branch_length
    integer(kind=i_def), pointer :: reconstruction_big_halo_stencil_size(:,:)
    integer(kind=i_def), pointer :: reconstruction_big_halo_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: reconstruction_big_halo_stencil_map

    nullify( set_counts_field_dummy_data, lookup_field_dummy_data, set_counts_field_data, &
             lookup_field_data, wind_big_halo_data, reconstruction_big_halo_data, advective_data, &
             map_adspc1_reconstruction_big_halo, map_adspc2_lookup_field, map_w2, &
             map_wtheta, mesh, set_counts_field_dummy_stencil_size, set_counts_field_dummy_stencil_dofmap, &
             set_counts_field_dummy_stencil_map, lookup_field_dummy_stencil_size, &
             lookup_field_dummy_stencil_dofmap, lookup_field_dummy_stencil_map, wind_big_halo_stencil_size, &
             wind_big_halo_stencil_dofmap, wind_big_halo_stencil_map, reconstruction_big_halo_stencil_size, &
             reconstruction_big_halo_stencil_dofmap, reconstruction_big_halo_stencil_map )

    !
    ! Initialise field and/or operator proxies
    !
    advective_proxy = advective%get_proxy()
    advective_data => advective_proxy%data
    reconstruction_big_halo_proxy = reconstruction_big_halo%get_proxy()
    reconstruction_big_halo_data => reconstruction_big_halo_proxy%data
    wind_big_halo_proxy = wind_big_halo%get_proxy()
    wind_big_halo_data => wind_big_halo_proxy%data
    lookup_field_proxy = lookup_field%get_proxy()
    lookup_field_data => lookup_field_proxy%data
    set_counts_field_proxy = set_counts_field%get_proxy()
    set_counts_field_data => set_counts_field_proxy%data
    lookup_field_dummy_proxy = lookup_field_dummy%get_proxy()
    lookup_field_dummy_data => lookup_field_dummy_proxy%data
    set_counts_field_dummy_proxy = set_counts_field_dummy%get_proxy()
    set_counts_field_dummy_data => set_counts_field_dummy_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_advective = advective_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => advective_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Initialise stencil dofmaps
    !
    reconstruction_big_halo_stencil_map => &
&reconstruction_big_halo_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    reconstruction_big_halo_max_branch_length = stencil_extent + 1_i_def
    reconstruction_big_halo_stencil_dofmap => reconstruction_big_halo_stencil_map%get_whole_dofmap()
    reconstruction_big_halo_stencil_size => reconstruction_big_halo_stencil_map%get_stencil_sizes()
    wind_big_halo_stencil_map => wind_big_halo_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    wind_big_halo_max_branch_length = stencil_extent + 1_i_def
    wind_big_halo_stencil_dofmap => wind_big_halo_stencil_map%get_whole_dofmap()
    wind_big_halo_stencil_size => wind_big_halo_stencil_map%get_stencil_sizes()
    lookup_field_dummy_stencil_map => lookup_field_dummy_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    lookup_field_dummy_max_branch_length = stencil_extent + 1_i_def
    lookup_field_dummy_stencil_dofmap => lookup_field_dummy_stencil_map%get_whole_dofmap()
    lookup_field_dummy_stencil_size => lookup_field_dummy_stencil_map%get_stencil_sizes()
    set_counts_field_dummy_stencil_map => &
&set_counts_field_dummy_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    set_counts_field_dummy_max_branch_length = stencil_extent + 1_i_def
    set_counts_field_dummy_stencil_dofmap => set_counts_field_dummy_stencil_map%get_whole_dofmap()
    set_counts_field_dummy_stencil_size => set_counts_field_dummy_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_wtheta => advective_proxy%vspace%get_whole_dofmap()
    map_adspc1_reconstruction_big_halo => reconstruction_big_halo_proxy%vspace%get_whole_dofmap()
    map_w2 => wind_big_halo_proxy%vspace%get_whole_dofmap()
    map_adspc2_lookup_field => lookup_field_proxy%vspace%get_whole_dofmap()
    map_adspc3_set_counts_field => set_counts_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for wtheta
    !
    ndf_wtheta = advective_proxy%vspace%get_ndf()
    undf_wtheta = advective_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc1_reconstruction_big_halo
    !
    ndf_adspc1_reconstruction_big_halo = reconstruction_big_halo_proxy%vspace%get_ndf()
    undf_adspc1_reconstruction_big_halo = reconstruction_big_halo_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for w2
    !
    ndf_w2 = wind_big_halo_proxy%vspace%get_ndf()
    undf_w2 = wind_big_halo_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_lookup_field
    !
    ndf_adspc2_lookup_field = lookup_field_proxy%vspace%get_ndf()
    undf_adspc2_lookup_field = lookup_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc3_set_counts_field
    !
    ndf_adspc3_set_counts_field = set_counts_field_proxy%vspace%get_ndf()
    undf_adspc3_set_counts_field = set_counts_field_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_halo_cell(loop_halo_depth)
    !
    ! Call kernels and communication routines
    !
    if (advective_proxy%is_dirty(depth=loop_halo_depth)) then
      call advective_proxy%halo_exchange(depth=loop_halo_depth)
    end if
    if (reconstruction_big_halo_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call reconstruction_big_halo_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    if (wind_big_halo_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call wind_big_halo_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    if (lookup_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call lookup_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    if (set_counts_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call set_counts_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    do cell = loop0_start, loop0_stop, 1
      call gen_poly_adv_upd_lookup_code(nlayers_advective, advective_data, reconstruction_big_halo_data, &
&reconstruction_big_halo_stencil_size(:,cell), reconstruction_big_halo_max_branch_length, &
&reconstruction_big_halo_stencil_dofmap(:,:,:,cell), wind_big_halo_data, wind_big_halo_stencil_size(:,cell), &
&wind_big_halo_max_branch_length, wind_big_halo_stencil_dofmap(:,:,:,cell), lookup_field_data, set_counts_field_data, &
&lookup_field_dummy_data, lookup_field_dummy_stencil_size(:,cell), lookup_field_dummy_max_branch_length, &
&lookup_field_dummy_stencil_dofmap(:,:,:,cell), set_counts_field_dummy_data, set_counts_field_dummy_stencil_size(:,cell), &
&set_counts_field_dummy_max_branch_length, set_counts_field_dummy_stencil_dofmap(:,:,:,cell), nsets_max, nindices, ndf_wtheta, &
&undf_wtheta, map_wtheta(:,cell), ndf_adspc1_reconstruction_big_halo, undf_adspc1_reconstruction_big_halo, &
&map_adspc1_reconstruction_big_halo(:,cell), ndf_w2, undf_w2, map_w2(:,cell), ndf_adspc2_lookup_field, undf_adspc2_lookup_field, &
&map_adspc2_lookup_field(:,cell), ndf_adspc3_set_counts_field, undf_adspc3_set_counts_field, map_adspc3_set_counts_field(:,cell))
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call lookup_field_proxy%set_dirty()
    call lookup_field_proxy%set_clean(loop_halo_depth)
    call set_counts_field_proxy%set_dirty()
    call set_counts_field_proxy%set_clean(loop_halo_depth)
    !
    !
  end subroutine invoke_gen_poly_adv_upd_lookup_kernel

  !> @brief Additional edit required due to PSyclone issue #2934. Use of integer fields with LMA operators.
  subroutine invoke_gen_w3h_adv_upd_lookup_kernel(advective_increment, wind, m3_inv, lookup_field, set_counts_field, &
lookup_field_dummy, set_counts_field_dummy, nsets_max, nindices, stencil_extent, loop_halo_depth)
    use gen_w3h_adv_upd_lookup_kernel_mod, only: gen_w3h_adv_upd_lookup_code
    use mesh_mod, only: mesh_type
    use stencil_2D_dofmap_mod, only: stencil_2D_dofmap_type, STENCIL_2D_CROSS
    integer(kind=i_def), intent(in) :: nsets_max, nindices
    type(r_tran_field_type), intent(in) :: advective_increment, wind
    type(operator_type), intent(in) :: m3_inv
    type(integer_field_type), intent(in) :: lookup_field, set_counts_field, lookup_field_dummy, set_counts_field_dummy
    integer(kind=i_def), intent(in) :: stencil_extent
    integer, intent(in) :: loop_halo_depth
    integer(kind=i_def) cell
    integer(kind=i_def) loop0_start, loop0_stop
    integer(kind=i_def) nlayers_advective_increment
    real(kind=r_def), pointer, dimension(:,:,:) :: m3_inv_local_stencil
    type(operator_proxy_type) m3_inv_proxy
    integer(kind=i_def), pointer, dimension(:) :: set_counts_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: set_counts_field_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_field_data
    type(integer_field_proxy_type) lookup_field_proxy, set_counts_field_proxy, lookup_field_dummy_proxy, &
&set_counts_field_dummy_proxy
    real(kind=r_tran), pointer, dimension(:) :: wind_data
    real(kind=r_tran), pointer, dimension(:) :: advective_increment_data
    type(r_tran_field_proxy_type) advective_increment_proxy, wind_proxy
    integer(kind=i_def), pointer :: map_adspc1_lookup_field(:,:), map_adspc2_set_counts_field(:,:), &
&map_any_w2(:,:), map_w3(:,:)
    integer(kind=i_def) ndf_w3, undf_w3, ndf_any_w2, undf_any_w2, ndf_adspc1_lookup_field, undf_adspc1_lookup_field, &
&ndf_adspc2_set_counts_field, undf_adspc2_set_counts_field
    integer(kind=i_def) max_halo_depth_mesh
    type(mesh_type), pointer :: mesh
    integer(kind=i_def) set_counts_field_dummy_max_branch_length
    integer(kind=i_def), pointer :: set_counts_field_dummy_stencil_size(:,:)
    integer(kind=i_def), pointer :: set_counts_field_dummy_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: set_counts_field_dummy_stencil_map
    integer(kind=i_def) lookup_field_dummy_max_branch_length
    integer(kind=i_def), pointer :: lookup_field_dummy_stencil_size(:,:)
    integer(kind=i_def), pointer :: lookup_field_dummy_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: lookup_field_dummy_stencil_map
    integer(kind=i_def) wind_max_branch_length
    integer(kind=i_def), pointer :: wind_stencil_size(:,:)
    integer(kind=i_def), pointer :: wind_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: wind_stencil_map

    nullify( m3_inv_local_stencil, set_counts_field_dummy_data, lookup_field_dummy_data, &
             set_counts_field_data, lookup_field_data, wind_data, advective_increment_data, &
             map_adspc1_lookup_field, map_adspc2_set_counts_field, map_any_w2, map_w3, &
             mesh, &
             set_counts_field_dummy_stencil_size, set_counts_field_dummy_stencil_dofmap, set_counts_field_dummy_stencil_map, &
             lookup_field_dummy_stencil_size, lookup_field_dummy_stencil_dofmap, lookup_field_dummy_stencil_map, &
             wind_stencil_size, wind_stencil_dofmap, wind_stencil_map )

    !
    ! Initialise field and/or operator proxies
    !
    advective_increment_proxy = advective_increment%get_proxy()
    advective_increment_data => advective_increment_proxy%data
    wind_proxy = wind%get_proxy()
    wind_data => wind_proxy%data
    m3_inv_proxy = m3_inv%get_proxy()
    m3_inv_local_stencil => m3_inv_proxy%local_stencil
    lookup_field_proxy = lookup_field%get_proxy()
    lookup_field_data => lookup_field_proxy%data
    set_counts_field_proxy = set_counts_field%get_proxy()
    set_counts_field_data => set_counts_field_proxy%data
    lookup_field_dummy_proxy = lookup_field_dummy%get_proxy()
    lookup_field_dummy_data => lookup_field_dummy_proxy%data
    set_counts_field_dummy_proxy = set_counts_field_dummy%get_proxy()
    set_counts_field_dummy_data => set_counts_field_dummy_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_advective_increment = advective_increment_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => advective_increment_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Initialise stencil dofmaps
    !
    wind_stencil_map => wind_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    wind_max_branch_length = stencil_extent + 1_i_def
    wind_stencil_dofmap => wind_stencil_map%get_whole_dofmap()
    wind_stencil_size => wind_stencil_map%get_stencil_sizes()
    lookup_field_dummy_stencil_map => lookup_field_dummy_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    lookup_field_dummy_max_branch_length = stencil_extent + 1_i_def
    lookup_field_dummy_stencil_dofmap => lookup_field_dummy_stencil_map%get_whole_dofmap()
    lookup_field_dummy_stencil_size => lookup_field_dummy_stencil_map%get_stencil_sizes()
    set_counts_field_dummy_stencil_map => &
&set_counts_field_dummy_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    set_counts_field_dummy_max_branch_length = stencil_extent + 1_i_def
    set_counts_field_dummy_stencil_dofmap => set_counts_field_dummy_stencil_map%get_whole_dofmap()
    set_counts_field_dummy_stencil_size => set_counts_field_dummy_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_w3 => advective_increment_proxy%vspace%get_whole_dofmap()
    map_any_w2 => wind_proxy%vspace%get_whole_dofmap()
    map_adspc1_lookup_field => lookup_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_set_counts_field => set_counts_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for w3
    !
    ndf_w3 = advective_increment_proxy%vspace%get_ndf()
    undf_w3 = advective_increment_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for any_w2
    !
    ndf_any_w2 = wind_proxy%vspace%get_ndf()
    undf_any_w2 = wind_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc1_lookup_field
    !
    ndf_adspc1_lookup_field = lookup_field_proxy%vspace%get_ndf()
    undf_adspc1_lookup_field = lookup_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_set_counts_field
    !
    ndf_adspc2_set_counts_field = set_counts_field_proxy%vspace%get_ndf()
    undf_adspc2_set_counts_field = set_counts_field_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_halo_cell(loop_halo_depth)
    !
    ! Call kernels and communication routines
    !

    if (advective_increment_proxy%is_dirty(depth=loop_halo_depth)) then
      call advective_increment_proxy%halo_exchange(depth=loop_halo_depth)
    end if
    if (wind_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call wind_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    if (lookup_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call lookup_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    if (set_counts_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call set_counts_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    do cell = loop0_start, loop0_stop, 1
      call gen_w3h_adv_upd_lookup_code(cell, nlayers_advective_increment, advective_increment_data, wind_data, &
&wind_stencil_size(:,cell), wind_max_branch_length, wind_stencil_dofmap(:,:,:,cell), m3_inv_proxy%ncell_3d, &
&m3_inv_local_stencil, lookup_field_data, set_counts_field_data, &
&lookup_field_dummy_data, lookup_field_dummy_stencil_size(:,cell), lookup_field_dummy_max_branch_length, &
&lookup_field_dummy_stencil_dofmap(:,:,:,cell), set_counts_field_dummy_data, set_counts_field_dummy_stencil_size(:,cell), &
&set_counts_field_dummy_max_branch_length, set_counts_field_dummy_stencil_dofmap(:,:,:,cell), nsets_max, nindices, ndf_w3, &
&undf_w3, map_w3(:,cell), ndf_any_w2, undf_any_w2, map_any_w2(:,cell), ndf_adspc1_lookup_field, undf_adspc1_lookup_field, &
&map_adspc1_lookup_field(:,cell), ndf_adspc2_set_counts_field, undf_adspc2_set_counts_field, map_adspc2_set_counts_field(:,cell))
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call lookup_field_proxy%set_dirty()
    call lookup_field_proxy%set_clean(loop_halo_depth)
    call set_counts_field_proxy%set_dirty()
    call set_counts_field_proxy%set_clean(loop_halo_depth)
    !
    !
  end subroutine invoke_gen_w3h_adv_upd_lookup_kernel

  ! Subroutine name shortened from invoke_gen_apply_helmholtz_op_lookup_kernel because of ifort compiler used in JEDI
  subroutine invoke_gen_a_h_o_lookup_kernel(dummy_w3_big_halo, lookup_field, set_counts_field, lookup_field_dummy, &
&set_counts_field_dummy, nindices, stencil_extent, loop_halo_depth)
    use gen_apply_helmholtz_op_lookup_kernel_mod, only: gen_apply_helmholtz_op_lookup_code
    use mesh_mod, only: mesh_type
    use stencil_2D_dofmap_mod, only: stencil_2D_dofmap_type, STENCIL_2D_CROSS
    integer(kind=i_def), intent(in) :: nindices
    type(r_solver_field_type), intent(in) :: dummy_w3_big_halo
    type(integer_field_type), intent(in) :: lookup_field, set_counts_field, lookup_field_dummy, set_counts_field_dummy
    integer(kind=i_def), intent(in) :: stencil_extent
    integer, intent(in) :: loop_halo_depth
    integer(kind=i_def) cell
    integer(kind=i_def) loop0_start, loop0_stop
    integer(kind=i_def) nlayers_dummy_w3_big_halo
    integer(kind=i_def), pointer, dimension(:) :: set_counts_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_field_dummy_data
    integer(kind=i_def), pointer, dimension(:) :: set_counts_field_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_field_data
    type(integer_field_proxy_type) lookup_field_proxy, set_counts_field_proxy, lookup_field_dummy_proxy, &
&set_counts_field_dummy_proxy
    real(kind=r_solver), pointer, dimension(:) :: dummy_w3_big_halo_data
    type(r_solver_field_proxy_type) dummy_w3_big_halo_proxy
    integer(kind=i_def), pointer :: map_adspc1_lookup_field(:,:), map_adspc2_set_counts_field(:,:), &
&map_w3(:,:)
    integer(kind=i_def) ndf_w3, undf_w3, ndf_adspc1_lookup_field, undf_adspc1_lookup_field, ndf_adspc2_set_counts_field, &
&undf_adspc2_set_counts_field
    integer(kind=i_def) max_halo_depth_mesh
    type(mesh_type), pointer :: mesh
    integer(kind=i_def) set_counts_field_dummy_max_branch_length
    integer(kind=i_def), pointer :: set_counts_field_dummy_stencil_size(:,:)
    integer(kind=i_def), pointer :: set_counts_field_dummy_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: set_counts_field_dummy_stencil_map
    integer(kind=i_def) lookup_field_dummy_max_branch_length
    integer(kind=i_def), pointer :: lookup_field_dummy_stencil_size(:,:)
    integer(kind=i_def), pointer :: lookup_field_dummy_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: lookup_field_dummy_stencil_map
    integer(kind=i_def) dummy_w3_big_halo_max_branch_length
    integer(kind=i_def), pointer :: dummy_w3_big_halo_stencil_size(:,:)
    integer(kind=i_def), pointer :: dummy_w3_big_halo_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: dummy_w3_big_halo_stencil_map
    nullify( set_counts_field_dummy_data, lookup_field_dummy_data, set_counts_field_data, &
             lookup_field_data, dummy_w3_big_halo_data, map_adspc1_lookup_field, &
             map_adspc2_set_counts_field, mesh, set_counts_field_dummy_stencil_size, &
             set_counts_field_dummy_stencil_dofmap, set_counts_field_dummy_stencil_map, &
             lookup_field_dummy_stencil_size, lookup_field_dummy_stencil_dofmap, &
             lookup_field_dummy_stencil_map, dummy_w3_big_halo_stencil_size, &
             dummy_w3_big_halo_stencil_dofmap, dummy_w3_big_halo_stencil_map )
    !
    ! Initialise field and/or operator proxies
    !
    dummy_w3_big_halo_proxy = dummy_w3_big_halo%get_proxy()
    dummy_w3_big_halo_data => dummy_w3_big_halo_proxy%data
    lookup_field_proxy = lookup_field%get_proxy()
    lookup_field_data => lookup_field_proxy%data
    set_counts_field_proxy = set_counts_field%get_proxy()
    set_counts_field_data => set_counts_field_proxy%data
    lookup_field_dummy_proxy = lookup_field_dummy%get_proxy()
    lookup_field_dummy_data => lookup_field_dummy_proxy%data
    set_counts_field_dummy_proxy = set_counts_field_dummy%get_proxy()
    set_counts_field_dummy_data => set_counts_field_dummy_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_dummy_w3_big_halo = dummy_w3_big_halo_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => dummy_w3_big_halo_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Initialise stencil dofmaps
    !
    dummy_w3_big_halo_stencil_map => dummy_w3_big_halo_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    dummy_w3_big_halo_max_branch_length = stencil_extent + 1_i_def
    dummy_w3_big_halo_stencil_dofmap => dummy_w3_big_halo_stencil_map%get_whole_dofmap()
    dummy_w3_big_halo_stencil_size => dummy_w3_big_halo_stencil_map%get_stencil_sizes()
    lookup_field_dummy_stencil_map => lookup_field_dummy_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    lookup_field_dummy_max_branch_length = stencil_extent + 1_i_def
    lookup_field_dummy_stencil_dofmap => lookup_field_dummy_stencil_map%get_whole_dofmap()
    lookup_field_dummy_stencil_size => lookup_field_dummy_stencil_map%get_stencil_sizes()
    set_counts_field_dummy_stencil_map => &
&set_counts_field_dummy_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    set_counts_field_dummy_max_branch_length = stencil_extent + 1_i_def
    set_counts_field_dummy_stencil_dofmap => set_counts_field_dummy_stencil_map%get_whole_dofmap()
    set_counts_field_dummy_stencil_size => set_counts_field_dummy_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_w3 => dummy_w3_big_halo_proxy%vspace%get_whole_dofmap()
    map_adspc1_lookup_field => lookup_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_set_counts_field => set_counts_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for w3
    !
    ndf_w3 = dummy_w3_big_halo_proxy%vspace%get_ndf()
    undf_w3 = dummy_w3_big_halo_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc1_lookup_field
    !
    ndf_adspc1_lookup_field = lookup_field_proxy%vspace%get_ndf()
    undf_adspc1_lookup_field = lookup_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_set_counts_field
    !
    ndf_adspc2_set_counts_field = set_counts_field_proxy%vspace%get_ndf()
    undf_adspc2_set_counts_field = set_counts_field_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_halo_cell(loop_halo_depth)
    !
    ! Call kernels and communication routines
    !
    if (lookup_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call lookup_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    if (set_counts_field_proxy%is_dirty(depth=loop_halo_depth + stencil_extent)) then
      call set_counts_field_proxy%halo_exchange(depth=loop_halo_depth + stencil_extent)
    end if
    do cell = loop0_start, loop0_stop, 1
      call gen_apply_helmholtz_op_lookup_code(nlayers_dummy_w3_big_halo, dummy_w3_big_halo_data, &
&dummy_w3_big_halo_stencil_size(:,cell), dummy_w3_big_halo_max_branch_length, dummy_w3_big_halo_stencil_dofmap(:,:,:,cell), &
&lookup_field_data, set_counts_field_data, lookup_field_dummy_data, lookup_field_dummy_stencil_size(:,cell), &
&lookup_field_dummy_max_branch_length, lookup_field_dummy_stencil_dofmap(:,:,:,cell), set_counts_field_dummy_data, &
&set_counts_field_dummy_stencil_size(:,cell), set_counts_field_dummy_max_branch_length, &
&set_counts_field_dummy_stencil_dofmap(:,:,:,cell), nindices, ndf_w3, undf_w3, map_w3(:,cell), ndf_adspc1_lookup_field, &
&undf_adspc1_lookup_field, map_adspc1_lookup_field(:,cell), ndf_adspc2_set_counts_field, undf_adspc2_set_counts_field, &
&map_adspc2_set_counts_field(:,cell))
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call lookup_field_proxy%set_dirty()
    call lookup_field_proxy%set_clean(loop_halo_depth)
    call set_counts_field_proxy%set_dirty()
    call set_counts_field_proxy%set_clean(loop_halo_depth)
    !
    !
  end subroutine invoke_gen_a_h_o_lookup_kernel

end module psykal_lite_gen_lookup_tables_psy_mod
