!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   PSyKAl lite code to compute the stencil adj_poly1d_recon_lookup_kernel.
!!          PSyclone issue: #2932.
!> @details PSy-lite code required here as halo exchanges must be enforced
!!          for fields being used in calculations. The lookup table acts as
!!          a stencil operation.

module invoke_adj_poly1d_recon_lookup_mod

  use constants_mod, only: r_tran, r_def, i_def
  use r_tran_field_mod, only: r_tran_field_type, r_tran_field_proxy_type
  use integer_field_mod, only: integer_field_type, integer_field_proxy_type
  implicit none
  contains

  subroutine invoke_adj_poly1d_recon_lookup(reconstruction, tracer, lookup_poly1d_field, &
                                            num_sets_poly1d_field, coeff, nsets, nindices, &
                                            stencil_extent)
    use adj_polynd_recon_zeroing_kernel_mod, only: adj_polynd_recon_zeroing_code
    use adj_poly1d_recon_lookup_kernel_mod, only: adj_poly1d_recon_lookup_code
    use mesh_mod, only: mesh_type
    integer(kind=i_def), intent(in) :: nsets, nindices, stencil_extent
    type(r_tran_field_type), intent(in) :: reconstruction, tracer, coeff
    type(integer_field_type), intent(in) :: lookup_poly1d_field, num_sets_poly1d_field
    integer(kind=i_def) cell
    integer(kind=i_def) loop1_start, loop1_stop
    integer(kind=i_def) loop0_start, loop0_stop
    integer(kind=i_def) nlayers_reconstruction
    integer(kind=i_def), pointer, dimension(:) :: num_sets_poly1d_field_data
    integer(kind=i_def), pointer, dimension(:) :: lookup_poly1d_field_data
    type(integer_field_proxy_type) lookup_poly1d_field_proxy, num_sets_poly1d_field_proxy
    real(kind=r_tran), pointer, dimension(:) :: coeff_data
    real(kind=r_tran), pointer, dimension(:) :: tracer_data
    real(kind=r_tran), pointer, dimension(:) :: reconstruction_data
    type(r_tran_field_proxy_type) reconstruction_proxy, tracer_proxy, coeff_proxy
    integer(kind=i_def), pointer :: map_adspc1_reconstruction(:,:), map_adspc2_tracer(:,:), &
&map_adspc3_lookup_poly1d_field(:,:), map_adspc4_num_sets_poly1d_field(:,:), map_adspc5_coeff(:,:)
    integer(kind=i_def) ndf_adspc1_reconstruction, undf_adspc1_reconstruction, ndf_adspc2_tracer, undf_adspc2_tracer, &
&ndf_adspc3_lookup_poly1d_field, undf_adspc3_lookup_poly1d_field, ndf_adspc4_num_sets_poly1d_field, &
&undf_adspc4_num_sets_poly1d_field, ndf_adspc5_coeff, undf_adspc5_coeff
    integer(kind=i_def) max_halo_depth_mesh
    type(mesh_type), pointer :: mesh

    nullify( num_sets_poly1d_field_data, lookup_poly1d_field_data, &
             coeff_data, tracer_data, reconstruction_data, &
             map_adspc1_reconstruction, map_adspc2_tracer, &
             map_adspc3_lookup_poly1d_field, map_adspc4_num_sets_poly1d_field, &
             map_adspc5_coeff, &
             mesh )

    !
    ! Initialise field and/or operator proxies
    !
    reconstruction_proxy = reconstruction%get_proxy()
    reconstruction_data => reconstruction_proxy%data
    tracer_proxy = tracer%get_proxy()
    tracer_data => tracer_proxy%data
    lookup_poly1d_field_proxy = lookup_poly1d_field%get_proxy()
    lookup_poly1d_field_data => lookup_poly1d_field_proxy%data
    num_sets_poly1d_field_proxy = num_sets_poly1d_field%get_proxy()
    num_sets_poly1d_field_data => num_sets_poly1d_field_proxy%data
    coeff_proxy = coeff%get_proxy()
    coeff_data => coeff_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_reconstruction = reconstruction_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => reconstruction_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_reconstruction => reconstruction_proxy%vspace%get_whole_dofmap()
    map_adspc2_tracer => tracer_proxy%vspace%get_whole_dofmap()
    map_adspc3_lookup_poly1d_field => lookup_poly1d_field_proxy%vspace%get_whole_dofmap()
    map_adspc4_num_sets_poly1d_field => num_sets_poly1d_field_proxy%vspace%get_whole_dofmap()
    map_adspc5_coeff => coeff_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for adspc1_reconstruction
    !
    ndf_adspc1_reconstruction = reconstruction_proxy%vspace%get_ndf()
    undf_adspc1_reconstruction = reconstruction_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_tracer
    !
    ndf_adspc2_tracer = tracer_proxy%vspace%get_ndf()
    undf_adspc2_tracer = tracer_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc3_lookup_poly1d_field
    !
    ndf_adspc3_lookup_poly1d_field = lookup_poly1d_field_proxy%vspace%get_ndf()
    undf_adspc3_lookup_poly1d_field = lookup_poly1d_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc4_num_sets_poly1d_field
    !
    ndf_adspc4_num_sets_poly1d_field = num_sets_poly1d_field_proxy%vspace%get_ndf()
    undf_adspc4_num_sets_poly1d_field = num_sets_poly1d_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc5_coeff
    !
    ndf_adspc5_coeff = coeff_proxy%vspace%get_ndf()
    undf_adspc5_coeff = coeff_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_edge_cell()
    loop1_start = 1
    loop1_stop = mesh%get_last_edge_cell()
    !
    ! Call kernels and communication routines
    !

    ! EDIT(JC): Call halo exchange on reconstruction field, as this is used like a stencil.
    if (reconstruction_proxy%is_dirty(depth=stencil_extent)) then
      call reconstruction_proxy%halo_exchange(depth=stencil_extent)
    end if
    if (coeff_proxy%is_dirty(depth=stencil_extent)) then
      call coeff_proxy%halo_exchange(depth=stencil_extent)
    end if

    !$omp parallel default(shared), private(cell)
    !$omp do schedule(static)
    do cell = loop0_start, loop0_stop, 1
      call adj_poly1d_recon_lookup_code(nlayers_reconstruction, reconstruction_data, tracer_data, lookup_poly1d_field_data, &
&num_sets_poly1d_field_data, coeff_data, nsets, nindices, ndf_adspc1_reconstruction, undf_adspc1_reconstruction, &
&map_adspc1_reconstruction(:,cell), ndf_adspc2_tracer, undf_adspc2_tracer, map_adspc2_tracer(:,cell), &
&ndf_adspc3_lookup_poly1d_field, undf_adspc3_lookup_poly1d_field, map_adspc3_lookup_poly1d_field(:,cell), &
&ndf_adspc4_num_sets_poly1d_field, undf_adspc4_num_sets_poly1d_field, map_adspc4_num_sets_poly1d_field(:,cell), ndf_adspc5_coeff, &
&undf_adspc5_coeff, map_adspc5_coeff(:,cell))
    end do
    !$omp end do
    !$omp end parallel

    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call tracer_proxy%set_dirty()
    !
    do cell = loop1_start, loop1_stop, 1
      call adj_polynd_recon_zeroing_code(nlayers_reconstruction, reconstruction_data, tracer_data, ndf_adspc1_reconstruction, &
&undf_adspc1_reconstruction, map_adspc1_reconstruction(:,cell), ndf_adspc2_tracer, undf_adspc2_tracer, map_adspc2_tracer(:,cell))
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call reconstruction_proxy%set_dirty()
    !
    !
  end subroutine invoke_adj_poly1d_recon_lookup
end module invoke_adj_poly1d_recon_lookup_mod
