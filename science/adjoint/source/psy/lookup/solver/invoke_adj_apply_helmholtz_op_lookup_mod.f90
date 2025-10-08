!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   PSyKAl lite code to compute the stencil invoke_adj_apply_helmholtz_operator_kernel.
!!          PSyclone issue: #2932.
!> @details PSy-lite code required here as halo exchanges must be enforced
!!          for fields being used in calculations. The lookup table acts as
!!          a stencil operation.

! Module name shortened from invoke_adj_apply_helmholtz_op_lookup_kernel_mod because of ifort compiler used in JEDI
MODULE invoke_adj_a_h_o_lookup_kernel_mod
  USE constants_mod, ONLY: r_solver, r_def, l_def, i_def
  USE field_mod, ONLY: field_type, field_proxy_type
  USE r_solver_field_mod, ONLY: r_solver_field_type, r_solver_field_proxy_type
  USE integer_field_mod, ONLY: integer_field_type, integer_field_proxy_type
  IMPLICIT NONE
  CONTAINS
  ! Subroutine name shortened from invoke_adj_apply_helmholtz_op_lookup_kernel because of ifort compiler used in JEDI
  SUBROUTINE invoke_adj_a_h_o_lookup_kernel(vector_mx, vector_amx, helmholtz_operator, lookup_field, &
&set_counts_field, stencil_extent, nindices, lam_mesh)
    USE adj_apply_helmholtz_op_lookup_kernel_mod, ONLY: adj_apply_helmholtz_op_lookup_code
    USE mesh_mod, ONLY: mesh_type
    INTEGER(KIND=i_def), intent(in) :: stencil_extent
    INTEGER(KIND=i_def), intent(in) :: nindices
    LOGICAL(KIND=l_def), intent(in) :: lam_mesh
    TYPE(r_solver_field_type), intent(in) :: vector_mx, vector_amx, helmholtz_operator(9)
    TYPE(integer_field_type), intent(in) :: lookup_field, set_counts_field
    INTEGER(KIND=i_def) cell
    INTEGER(KIND=i_def) i
    INTEGER(KIND=i_def) loop0_start, loop0_stop
    INTEGER(KIND=i_def) nlayers_vector_mx
    INTEGER(KIND=i_def), pointer, dimension(:) :: set_counts_field_data
    INTEGER(KIND=i_def), pointer, dimension(:) :: lookup_field_data
    TYPE(integer_field_proxy_type) lookup_field_proxy, set_counts_field_proxy
    REAL(KIND=r_solver), pointer, dimension(:) :: helmholtz_operator_1_data, helmholtz_operator_2_data, &
&helmholtz_operator_3_data, helmholtz_operator_4_data, helmholtz_operator_5_data, &
&helmholtz_operator_6_data, helmholtz_operator_7_data, helmholtz_operator_8_data, &
&helmholtz_operator_9_data
    REAL(KIND=r_solver), pointer, dimension(:) :: vector_amx_data
    REAL(KIND=r_solver), pointer, dimension(:) :: vector_mx_data
    TYPE(r_solver_field_proxy_type) vector_mx_proxy, vector_amx_proxy, helmholtz_operator_proxy(9)
    INTEGER(KIND=i_def), pointer :: map_adspc1_lookup_field(:,:), map_adspc2_set_counts_field(:,:), &
&map_w3(:,:)
    INTEGER(KIND=i_def) ndf_w3, undf_w3, ndf_adspc1_lookup_field, undf_adspc1_lookup_field, ndf_adspc2_set_counts_field, &
&undf_adspc2_set_counts_field
    INTEGER(KIND=i_def) max_halo_depth_mesh
    TYPE(mesh_type), pointer :: mesh

    nullify( set_counts_field_data, lookup_field_data, helmholtz_operator_1_data, &
             helmholtz_operator_2_data, helmholtz_operator_4_data, helmholtz_operator_5_data, &
             helmholtz_operator_7_data, helmholtz_operator_8_data, vector_amx_data, &
             vector_mx_data, map_adspc1_lookup_field, map_adspc2_set_counts_field, mesh )

    !
    ! Initialise field and/or operator proxies
    !
    vector_mx_proxy = vector_mx%get_proxy()
    vector_mx_data => vector_mx_proxy%data
    vector_amx_proxy = vector_amx%get_proxy()
    vector_amx_data => vector_amx_proxy%data
    helmholtz_operator_proxy(1) = helmholtz_operator(1)%get_proxy()
    helmholtz_operator_1_data => helmholtz_operator_proxy(1)%data
    helmholtz_operator_proxy(2) = helmholtz_operator(2)%get_proxy()
    helmholtz_operator_2_data => helmholtz_operator_proxy(2)%data
    helmholtz_operator_proxy(3) = helmholtz_operator(3)%get_proxy()
    helmholtz_operator_3_data => helmholtz_operator_proxy(3)%data
    helmholtz_operator_proxy(4) = helmholtz_operator(4)%get_proxy()
    helmholtz_operator_4_data => helmholtz_operator_proxy(4)%data
    helmholtz_operator_proxy(5) = helmholtz_operator(5)%get_proxy()
    helmholtz_operator_5_data => helmholtz_operator_proxy(5)%data
    helmholtz_operator_proxy(6) = helmholtz_operator(6)%get_proxy()
    helmholtz_operator_6_data => helmholtz_operator_proxy(6)%data
    helmholtz_operator_proxy(7) = helmholtz_operator(7)%get_proxy()
    helmholtz_operator_7_data => helmholtz_operator_proxy(7)%data
    helmholtz_operator_proxy(8) = helmholtz_operator(8)%get_proxy()
    helmholtz_operator_8_data => helmholtz_operator_proxy(8)%data
    helmholtz_operator_proxy(9) = helmholtz_operator(9)%get_proxy()
    helmholtz_operator_9_data => helmholtz_operator_proxy(9)%data
    lookup_field_proxy = lookup_field%get_proxy()
    lookup_field_data => lookup_field_proxy%data
    set_counts_field_proxy = set_counts_field%get_proxy()
    set_counts_field_data => set_counts_field_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_vector_mx = vector_mx_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => vector_mx_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Look-up dofmaps for each function space
    !
    map_w3 => vector_mx_proxy%vspace%get_whole_dofmap()
    map_adspc1_lookup_field => lookup_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_set_counts_field => set_counts_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for w3
    !
    ndf_w3 = vector_mx_proxy%vspace%get_ndf()
    undf_w3 = vector_mx_proxy%vspace%get_undf()
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
    loop0_stop = mesh%get_last_edge_cell()
    !
    ! Call kernels and communication routines
    !
    do i = LBOUND(helmholtz_operator_proxy, dim=1), UBOUND(helmholtz_operator_proxy, dim=1)
      if (helmholtz_operator_proxy(i)%is_dirty(depth=stencil_extent)) then
        call helmholtz_operator_proxy(i)%halo_exchange(depth=stencil_extent)
      end if
    end do
    call vector_mx_proxy%halo_exchange(depth=stencil_extent)

    !$omp parallel default(shared), private(cell)
    !$omp do schedule(static)
    DO cell = loop0_start, loop0_stop, 1
      CALL adj_apply_helmholtz_op_lookup_code(nlayers_vector_mx, vector_mx_data, vector_amx_data, helmholtz_operator_1_data, &
&helmholtz_operator_2_data, helmholtz_operator_3_data, helmholtz_operator_4_data, helmholtz_operator_5_data, &
&helmholtz_operator_6_data, helmholtz_operator_7_data, helmholtz_operator_8_data, helmholtz_operator_9_data, lookup_field_data, &
&set_counts_field_data, nindices, lam_mesh, ndf_w3, undf_w3, map_w3(:,cell), ndf_adspc1_lookup_field, undf_adspc1_lookup_field, &
&map_adspc1_lookup_field(:,cell), ndf_adspc2_set_counts_field, undf_adspc2_set_counts_field, map_adspc2_set_counts_field(:,cell))
    END DO
    !$omp end do
    !$omp end parallel
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    CALL vector_mx_proxy%set_dirty()
    !
    !
  END SUBROUTINE invoke_adj_a_h_o_lookup_kernel
END MODULE invoke_adj_a_h_o_lookup_kernel_mod
