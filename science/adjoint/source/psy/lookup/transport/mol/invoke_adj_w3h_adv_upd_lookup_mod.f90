!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   PSyKAl lite code to compute the stencil adj_w3h_advective_update_kernel.
!!          PSyclone issues: #2932, #2934.
!> @details #2934 Needed because attempting to parse the kernel code gives an error:
!!
!!          "Parse Error: In the LFRic API a kernel that has an LMA operator argument must
!!           only have field arguments with 'gh_real' data type but kernel
!!           'adj_w3h_adv_upd_lookup_kernel_type' has a field argument with 'gh_integer' data type."
!!
!!          Additionally halo exchanges must be enforced for fields being used in
!!          calculations. The lookup table acts as a stencil operation. Issue #2932.

module invoke_adj_w3h_adv_upd_lookup_mod

  implicit none

  public  :: invoke_adj_w3h_adv_upd_lookup_kernel_type

  contains

  !=============================================================================
  !> @brief PSy layer invoke code for adj_w3h_advective_update_kernel.
  !> @param[in]     advective_increment Advective update field
  !> @param[in,out] tracer                            Pointwise tracer field to advect stored on cell faces
  !> @param[in]     lookup_w3h_adv_upd                Lookup table between read and write indices
  !> @param[in]     set_count_w3h_adv_upd_field_data  Index set counters for the lookup table
  !> @param[in]     wind                              Wind field
  !> @param[in]     stencil_extent                    Extent of stencil
  !> @param[in]     m3_inv                            Inverse mass matrix for W3 space
  !> @param[in]     nsets                             Number of index sets per cell in lookup table
  !> @param[in]     nindices                          Number of indices per cell in lookup table
  subroutine invoke_adj_w3h_adv_upd_lookup_kernel_type( advective_increment,      &
                                                        tracer,                   &
                                                        lookup_w3h_adv_upd_field, &
                                                        set_count_w3h_adv_upd_field, &
                                                        wind,                     &
                                                        stencil_extent,           &
                                                        m3_inv,                   &
                                                        nsets,                    &
                                                        nindices )

    use adj_w3h_adv_upd_lookup_kernel_mod, only: adj_w3h_adv_upd_lookup_code
    use constants_mod,                     only: r_tran, r_def, i_def, l_def
    use mesh_mod,                          only: mesh_type
    use operator_mod,                      only: operator_type, &
                                                 operator_proxy_type
    use integer_field_mod,                 only: integer_field_type, &
                                                 integer_field_proxy_type
    use r_tran_field_mod,                  only: r_tran_field_type, &
                                                 r_tran_field_proxy_type
    use stencil_2D_dofmap_mod,             only: stencil_2D_dofmap_type, &
                                                 STENCIL_2D_CROSS

    implicit none

    ! Arguments
    type(r_tran_field_type),     intent(in) :: advective_increment
    type(r_tran_field_type),  intent(inout) :: tracer
    type(integer_field_type),    intent(in) :: lookup_w3h_adv_upd_field
    type(integer_field_type),    intent(in) :: set_count_w3h_adv_upd_field
    type(r_tran_field_type),     intent(in) :: wind
    integer(kind=i_def),         intent(in) :: stencil_extent
    type(operator_type),         intent(in) :: m3_inv
    integer(kind=i_def),         intent(in) :: nsets
    integer(kind=i_def),         intent(in) :: nindices

    ! Internal variables
    integer(kind=i_def)                         :: cell
    integer(kind=i_def)                         :: loop0_start
    integer(kind=i_def)                         :: loop0_stop
    integer(kind=i_def)                         :: nlayers_advective_increment
    real(kind=r_def), pointer, dimension(:,:,:) :: m3_inv_local_stencil
    type(operator_proxy_type)                   :: m3_inv_proxy
    real(kind=r_tran),    pointer, dimension(:) :: wind_data
    real(kind=r_tran),    pointer, dimension(:) :: tracer_data
    real(kind=r_tran),    pointer, dimension(:) :: advective_increment_data
    integer(kind=i_def),  pointer, dimension(:) :: lookup_w3h_adv_upd_field_data
    integer(kind=i_def),  pointer, dimension(:) :: set_count_w3h_adv_upd_field_data
    type(r_tran_field_proxy_type)               :: advective_increment_proxy
    type(r_tran_field_proxy_type)               :: tracer_proxy
    type(r_tran_field_proxy_type)               :: wind_proxy
    type(integer_field_proxy_type)              :: lookup_w3h_adv_upd_field_proxy
    type(integer_field_proxy_type)              :: set_count_w3h_adv_upd_field_proxy
    integer(kind=i_def),                pointer :: map_adspc1_tracer(:,:)
    integer(kind=i_def),                pointer :: map_adspc2_lookup_w3h_adv_upd_field(:,:)
    integer(kind=i_def),                pointer :: map_adspc2_set_count_w3h_adv_upd_field(:,:)
    integer(kind=i_def),                pointer :: map_any_w2(:,:)
    integer(kind=i_def),                pointer :: map_w3(:,:)
    integer(kind=i_def)                         :: ndf_w3
    integer(kind=i_def)                         :: undf_w3
    integer(kind=i_def)                         :: ndf_adspc1_tracer
    integer(kind=i_def)                         :: undf_adspc1_tracer
    integer(kind=i_def)                         :: ndf_any_w2
    integer(kind=i_def)                         :: undf_any_w2
    integer(kind=i_def)                         :: ndf_adspc2_lookup_w3h_adv_upd_field
    integer(kind=i_def)                         :: undf_adspc2_lookup_w3h_adv_upd_field
    integer(kind=i_def)                         :: ndf_adspc3_set_count_w3h_adv_upd_field
    integer(kind=i_def)                         :: undf_adspc3_set_count_w3h_adv_upd_field
    integer(kind=i_def)                         :: max_halo_depth_mesh
    type(mesh_type),                    pointer :: mesh
    integer(kind=i_def)                         :: wind_max_branch_length
    integer(kind=i_def),                pointer :: wind_stencil_size(:,:)
    integer(kind=i_def),                pointer :: wind_stencil_dofmap(:,:,:,:)
    type(stencil_2D_dofmap_type),       pointer :: wind_stencil_map

    nullify( m3_inv_local_stencil, &
             wind_data, tracer_data, advective_increment_data, &
             lookup_w3h_adv_upd_field_data, set_count_w3h_adv_upd_field_data, &
             map_adspc1_tracer, map_adspc2_lookup_w3h_adv_upd_field, &
             map_adspc2_set_count_w3h_adv_upd_field, map_any_w2, map_w3, &
             mesh, &
             wind_stencil_size, wind_stencil_dofmap, wind_stencil_map )

    !
    ! Initialise field and/or operator proxies
    !
    advective_increment_proxy = advective_increment%get_proxy()
    advective_increment_data => advective_increment_proxy%data
    tracer_proxy = tracer%get_proxy()
    tracer_data => tracer_proxy%data
    lookup_w3h_adv_upd_field_proxy = lookup_w3h_adv_upd_field%get_proxy()
    lookup_w3h_adv_upd_field_data => lookup_w3h_adv_upd_field_proxy%data
    set_count_w3h_adv_upd_field_proxy = set_count_w3h_adv_upd_field%get_proxy()
    set_count_w3h_adv_upd_field_data => set_count_w3h_adv_upd_field_proxy%data
    wind_proxy = wind%get_proxy()
    wind_data => wind_proxy%data
    m3_inv_proxy = m3_inv%get_proxy()
    m3_inv_local_stencil => m3_inv_proxy%local_stencil
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
    wind_stencil_map => wind_proxy%vspace%get_stencil_2D_dofmap( STENCIL_2D_CROSS, stencil_extent )
    wind_max_branch_length = stencil_extent + 1_i_def
    wind_stencil_dofmap => wind_stencil_map%get_whole_dofmap()
    wind_stencil_size => wind_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_w3 => advective_increment_proxy%vspace%get_whole_dofmap()
    map_adspc1_tracer => tracer_proxy%vspace%get_whole_dofmap()
    map_adspc2_lookup_w3h_adv_upd_field => lookup_w3h_adv_upd_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_set_count_w3h_adv_upd_field => set_count_w3h_adv_upd_field_proxy%vspace%get_whole_dofmap()
    map_any_w2 => wind_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for w3
    !
    ndf_w3 = advective_increment_proxy%vspace%get_ndf()
    undf_w3 = advective_increment_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc1_tracer
    !
    ndf_adspc1_tracer = tracer_proxy%vspace%get_ndf()
    undf_adspc1_tracer = tracer_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for lookup_w3h_adv_upd_field
    !
    ndf_adspc2_lookup_w3h_adv_upd_field = lookup_w3h_adv_upd_field_proxy%vspace%get_ndf()
    undf_adspc2_lookup_w3h_adv_upd_field = lookup_w3h_adv_upd_field_proxy%vspace%get_undf()

    ! Initialise number of DoFs for set_count_w3h_adv_upd_field
    !
    ndf_adspc3_set_count_w3h_adv_upd_field = set_count_w3h_adv_upd_field_proxy%vspace%get_ndf()
    undf_adspc3_set_count_w3h_adv_upd_field = set_count_w3h_adv_upd_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for any_w2
    !
    ndf_any_w2 = wind_proxy%vspace%get_ndf()
    undf_any_w2 = wind_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_edge_cell()
    !
    ! Call kernels and communication routines
    !
    if ( wind_proxy%is_dirty(depth=stencil_extent) ) then
      call wind_proxy%halo_exchange(depth=stencil_extent)
    end if
    if ( advective_increment_proxy%is_dirty(depth=stencil_extent) ) then
      call advective_increment_proxy%halo_exchange(depth=stencil_extent)
    end if

    !$omp parallel default(shared), private(cell)
    !$omp do schedule(static)
    do cell = loop0_start, loop0_stop
      call adj_w3h_adv_upd_lookup_code( cell,                                           &
                                        nlayers_advective_increment,                    &
                                        advective_increment_data,                       &
                                        tracer_data,                                    &
                                        lookup_w3h_adv_upd_field_data,                  &
                                        set_count_w3h_adv_upd_field_data,               &
                                        wind_data,                                      &
                                        wind_stencil_size(:,cell),                      &
                                        wind_max_branch_length,                         &
                                        wind_stencil_dofmap(:,:,:,cell),                &
                                        m3_inv_proxy%ncell_3d,                          &
                                        m3_inv_local_stencil,                           &
                                        nsets,                                          &
                                        nindices,                                       &
                                        ndf_w3,                                         &
                                        undf_w3,                                        &
                                        map_w3(:,cell),                                 &
                                        ndf_adspc1_tracer,                              &
                                        undf_adspc1_tracer,                             &
                                        map_adspc1_tracer(:,cell),                      &
                                        ndf_adspc2_lookup_w3h_adv_upd_field,            &
                                        undf_adspc2_lookup_w3h_adv_upd_field,           &
                                        map_adspc2_lookup_w3h_adv_upd_field(:,cell),    &
                                        ndf_adspc3_set_count_w3h_adv_upd_field,         &
                                        undf_adspc3_set_count_w3h_adv_upd_field,        &
                                        map_adspc2_set_count_w3h_adv_upd_field(:,cell), &
                                        ndf_any_w2,                                     &
                                        undf_any_w2,                                    &
                                        map_any_w2(:,cell) )

    end do
    !$omp end do
    !$omp end parallel

    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call tracer_proxy%set_dirty()
    !
    !
  end subroutine invoke_adj_w3h_adv_upd_lookup_kernel_type

end module invoke_adj_w3h_adv_upd_lookup_mod
