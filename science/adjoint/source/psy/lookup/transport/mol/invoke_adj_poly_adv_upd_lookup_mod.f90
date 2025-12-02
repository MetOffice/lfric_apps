!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   PSyKAl lite code to compute the stencil adj_poly_adv_upd_lookup_kernel.
!!          PSyclone issue: #2932.
!> @details PSy-lite code required here as halo exchanges must be enforced
!!          for fields being used in calculations. The lookup table acts as
!!          a stencil operation.

module invoke_adj_poly_adv_upd_lookup_mod
  use constants_mod, only: r_tran, r_def, i_def
  use r_tran_field_mod, only: r_tran_field_type, r_tran_field_proxy_type
  use integer_field_mod, only: integer_field_type, integer_field_proxy_type
  implicit none
  contains

  !> @brief PSy-lite code required here as halo exchanges must be enforced
  !>        for fields being used in calculations. The lookup table acts as
  !>        a stencil operation.
  subroutine invoke_adj_poly_adv_upd_lookup(advective, reconstruction, &
                                            lookup_poly_adv_upd_field, num_sets_poly_adv_upd_field, &
                                            wind, wind_dir, &
                                            nsets, nindices, stencil_extent)
    use adj_poly_adv_upd_lookup_kernel_mod, only: adj_poly_adv_upd_lookup_code
    use mesh_mod, only: mesh_type
    use stencil_2D_dofmap_mod, only: stencil_2D_dofmap_type, STENCIL_2D_CROSS
    integer(kind=i_def), intent(in) :: nsets, nindices
    type(r_tran_field_type), intent(in) :: advective, reconstruction, wind, wind_dir
    type(integer_field_type), intent(in) :: lookup_poly_adv_upd_field, num_sets_poly_adv_upd_field
    integer(kind=i_def), intent(in) :: stencil_extent
    integer(kind=i_def) df
    integer(kind=i_def) cell
    integer(kind=i_def) loop1_start, loop1_stop
    integer(kind=i_def) loop0_start, loop0_stop
    integer(kind=i_def) nlayers_advective
    integer(kind=i_def), pointer, dimension(:) :: num_sets_poly_adv_upd_field_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_poly_adv_upd_field_data
    type(integer_field_proxy_type) lookup_poly_adv_upd_field_proxy, num_sets_poly_adv_upd_field_proxy
    real(kind=r_tran), pointer, dimension(:) :: wind_dir_data
    real(kind=r_tran), pointer, dimension(:) :: wind_data
    real(kind=r_tran), pointer, dimension(:) :: reconstruction_data
    real(kind=r_tran), pointer, dimension(:) :: advective_data
    type(r_tran_field_proxy_type) advective_proxy, reconstruction_proxy, wind_proxy, wind_dir_proxy
    integer(kind=i_def), pointer :: map_adspc1_reconstruction(:,:), &
&map_adspc2_lookup_poly_adv_upd_field(:,:), map_adspc3_num_sets_poly_adv_upd_field(:,:), &
&map_w2(:,:), map_wtheta(:,:)
    integer(kind=i_def) ndf_wtheta, undf_wtheta, ndf_adspc1_reconstruction, undf_adspc1_reconstruction, &
&ndf_adspc2_lookup_poly_adv_upd_field, undf_adspc2_lookup_poly_adv_upd_field, ndf_adspc3_num_sets_poly_adv_upd_field, &
&undf_adspc3_num_sets_poly_adv_upd_field, ndf_w2, undf_w2, ndf_aspc1_advective, undf_aspc1_advective
    integer(kind=i_def) max_halo_depth_mesh
    type(mesh_type), pointer :: mesh
    integer(kind=i_def) wind_max_branch_length
    integer(kind=i_def), pointer :: wind_stencil_size(:,:)
    integer(kind=i_def), pointer :: wind_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type), pointer :: wind_stencil_map

    nullify( num_sets_poly_adv_upd_field_data, lookup_poly_adv_upd_field_data, &
             wind_dir_data, wind_data, reconstruction_data, advective_data, &
             map_adspc1_reconstruction, map_adspc2_lookup_poly_adv_upd_field, &
             map_adspc3_num_sets_poly_adv_upd_field, map_w2, map_wtheta, &
             mesh, &
             wind_stencil_size, wind_stencil_dofmap, wind_stencil_map )

    !
    ! Initialise field and/or operator proxies
    !
    advective_proxy = advective%get_proxy()
    advective_data => advective_proxy%data
    reconstruction_proxy = reconstruction%get_proxy()
    reconstruction_data => reconstruction_proxy%data
    lookup_poly_adv_upd_field_proxy = lookup_poly_adv_upd_field%get_proxy()
    lookup_poly_adv_upd_field_data => lookup_poly_adv_upd_field_proxy%data
    num_sets_poly_adv_upd_field_proxy = num_sets_poly_adv_upd_field%get_proxy()
    num_sets_poly_adv_upd_field_data => num_sets_poly_adv_upd_field_proxy%data
    wind_proxy = wind%get_proxy()
    wind_data => wind_proxy%data
    wind_dir_proxy = wind_dir%get_proxy()
    wind_dir_data => wind_dir_proxy%data
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
    wind_stencil_map => wind_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    wind_max_branch_length = stencil_extent + 1_i_def
    wind_stencil_dofmap => wind_stencil_map%get_whole_dofmap()
    wind_stencil_size => wind_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_wtheta => advective_proxy%vspace%get_whole_dofmap()
    map_adspc1_reconstruction => reconstruction_proxy%vspace%get_whole_dofmap()
    map_adspc2_lookup_poly_adv_upd_field => lookup_poly_adv_upd_field_proxy%vspace%get_whole_dofmap()
    map_adspc3_num_sets_poly_adv_upd_field => num_sets_poly_adv_upd_field_proxy%vspace%get_whole_dofmap()
    map_w2 => wind_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for wtheta
    !
    ndf_wtheta = advective_proxy%vspace%get_ndf()
    undf_wtheta = advective_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc1_reconstruction
    !
    ndf_adspc1_reconstruction = reconstruction_proxy%vspace%get_ndf()
    undf_adspc1_reconstruction = reconstruction_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_lookup_poly_adv_upd_field
    !
    ndf_adspc2_lookup_poly_adv_upd_field = lookup_poly_adv_upd_field_proxy%vspace%get_ndf()
    undf_adspc2_lookup_poly_adv_upd_field = lookup_poly_adv_upd_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc3_num_sets_poly_adv_upd_field
    !
    ndf_adspc3_num_sets_poly_adv_upd_field = num_sets_poly_adv_upd_field_proxy%vspace%get_ndf()
    undf_adspc3_num_sets_poly_adv_upd_field = num_sets_poly_adv_upd_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for w2
    !
    ndf_w2 = wind_proxy%vspace%get_ndf()
    undf_w2 = wind_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for aspc1_advective
    !
    ndf_aspc1_advective = advective_proxy%vspace%get_ndf()
    undf_aspc1_advective = advective_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_edge_cell()
    loop1_start = 1
    loop1_stop = advective_proxy%vspace%get_last_dof_halo(1)
    !
    ! Call kernels and communication routines
    !
    if (wind_proxy%is_dirty(depth=stencil_extent)) then
      call wind_proxy%halo_exchange(depth=stencil_extent)
    end if

    ! EDIT(JC): Call halo exchange for stencil-like accesses.
    if (advective_proxy%is_dirty(depth=stencil_extent)) then
      call advective_proxy%halo_exchange(depth=stencil_extent)
    end if
    if (wind_dir_proxy%is_dirty(depth=stencil_extent)) then
      call wind_dir_proxy%halo_exchange(depth=stencil_extent)
    end if

    !$omp parallel default(shared), private(cell)
    !$omp do schedule(static)
    do cell = loop0_start, loop0_stop, 1
      call adj_poly_adv_upd_lookup_code(nlayers_advective, advective_data, reconstruction_data, lookup_poly_adv_upd_field_data, &
&num_sets_poly_adv_upd_field_data, wind_data, wind_stencil_size(:,cell), wind_max_branch_length, wind_stencil_dofmap(:,:,:,cell), &
&wind_dir_data, nsets, nindices, ndf_wtheta, undf_wtheta, map_wtheta(:,cell), ndf_adspc1_reconstruction, &
&undf_adspc1_reconstruction, map_adspc1_reconstruction(:,cell), ndf_adspc2_lookup_poly_adv_upd_field, &
&undf_adspc2_lookup_poly_adv_upd_field, map_adspc2_lookup_poly_adv_upd_field(:,cell), ndf_adspc3_num_sets_poly_adv_upd_field, &
&undf_adspc3_num_sets_poly_adv_upd_field, map_adspc3_num_sets_poly_adv_upd_field(:,cell), ndf_w2, undf_w2, map_w2(:,cell))
    end do
    !$omp end do
    !$omp end parallel

    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call advective_proxy%set_dirty()
    call reconstruction_proxy%set_dirty()
    !
    do df = loop1_start, loop1_stop, 1
      ! Built-in: setval_c (set a real-valued field to a real scalar value)
      advective_data(df) = 0.0_r_def
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call advective_proxy%set_dirty()
    call advective_proxy%set_clean(1)
    !
    !
  end subroutine invoke_adj_poly_adv_upd_lookup

end module invoke_adj_poly_adv_upd_lookup_mod
