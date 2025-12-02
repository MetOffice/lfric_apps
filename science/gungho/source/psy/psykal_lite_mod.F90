!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of the PSy layer

!> @details Contains hand-rolled versions of the PSy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_mod

  use field_mod,                    only : field_type, field_proxy_type
  use r_tran_field_mod,             only : r_tran_field_type, r_tran_field_proxy_type
  use integer_field_mod,            only : integer_field_type, integer_field_proxy_type
  use scalar_mod,                   only : scalar_type
  use operator_mod,                 only : operator_type, operator_proxy_type,                   &
                                           r_solver_operator_type, r_solver_operator_proxy_type, &
                                           r_tran_operator_type, r_tran_operator_proxy_type
  use constants_mod,                only : r_def, i_def, r_double, r_solver, r_tran, l_def, cache_block
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use mesh_mod,                     only : mesh_type
  use function_space_mod,           only : BASIS, DIFF_BASIS

  use quadrature_xyoz_mod,          only : quadrature_xyoz_type, &
                                           quadrature_xyoz_proxy_type
  use quadrature_face_mod,          only : quadrature_face_type, &
                                           quadrature_face_proxy_type
  use r_solver_field_mod,           only : r_solver_field_type, &
                                           r_solver_field_proxy_type

  implicit none
  public

contains

!> Non pointwise Kernels

!-------------------------------------------------------------------------------
  subroutine invoke_compute_dof_level_kernel(level)

  use sci_compute_dof_level_kernel_mod, only: compute_dof_level_code
  use mesh_mod,                     only: mesh_type
  implicit none

  type(field_type), intent(inout) :: level
  type(field_proxy_type) :: l_p
  integer :: cell, ndf, undf
  real(kind=r_def), pointer :: nodes(:,:) => null()
  integer, pointer :: map(:) => null()
  type(mesh_type), pointer :: mesh => null()
  l_p = level%get_proxy()
  undf = l_p%vspace%get_undf()
  ndf  = l_p%vspace%get_ndf()
  nodes => l_p%vspace%get_nodes( )

  mesh => l_p%vspace%get_mesh()
  do cell = 1,mesh%get_last_halo_cell(1)
    map => l_p%vspace%get_cell_dofmap(cell)
    call compute_dof_level_code(l_p%vspace%get_nlayers(),                 &
                                l_p%data,                                 &
                                ndf,                                      &
                                undf,                                     &
                                map,                                      &
                                nodes                                     &
                               )
  end do
  call l_p%set_dirty()

  end subroutine invoke_compute_dof_level_kernel

  !----------------------------------------------------------------------------
  ! This requires a stencil of horizontal cells for the operators
  ! see PSyclone #1103: https://github.com/stfc/PSyclone/issues/1103
  ! The LFRic infrastructure for this will be introduced in #2532
  subroutine invoke_helmholtz_operator_kernel_type(helmholtz_operator, hb_lumped_inv, stencil_depth, u_normalisation, div_star, &
                                                   t_normalisation, ptheta2v, compound_div, m3_exner_star, p3theta, w2_mask)
    use helmholtz_operator_kernel_mod, only: helmholtz_operator_code
    use mesh_mod, only: mesh_type
    use stencil_2d_dofmap_mod, only: stencil_2d_cross
    use stencil_2d_dofmap_mod, only: stencil_2d_dofmap_type
    use reference_element_mod, only: reference_element_type
    implicit none

    type(r_solver_field_type), intent(in) :: helmholtz_operator(9)
    type(r_solver_field_type), intent(in) ::  hb_lumped_inv, u_normalisation, t_normalisation, w2_mask
    type(r_solver_operator_type), intent(in) :: div_star, ptheta2v, compound_div, m3_exner_star, p3theta
    integer(kind=i_def), intent(in) :: stencil_depth
    integer(kind=i_def) :: stencil_size
    integer(kind=i_def) cell
    integer(kind=i_def) nlayers
    type(r_solver_operator_proxy_type) div_star_proxy, ptheta2v_proxy, compound_div_proxy, m3_exner_star_proxy, p3theta_proxy
    type(r_solver_field_proxy_type) helmholtz_operator_proxy(9)
    type(r_solver_field_proxy_type) hb_lumped_inv_proxy, u_normalisation_proxy, t_normalisation_proxy, &
                           w2_mask_proxy
    integer(kind=i_def), pointer :: map_w2(:,:) => null(), map_w3(:,:) => null(), map_wtheta(:,:) => null()
    integer(kind=i_def) ndf_w3, undf_w3, ndf_w2, undf_w2, ndf_wtheta, undf_wtheta
    type(mesh_type), pointer :: mesh => null()
    INTEGER(KIND=i_def) :: hb_lumped_inv_max_branch_length
    integer(kind=i_def), pointer :: hb_lumped_inv_stencil_sizes(:,:) => null()
    integer(kind=i_def), pointer :: hb_lumped_inv_stencil_dofmap(:,:,:,:) => null()
    type(stencil_2d_dofmap_type), pointer :: hb_lumped_inv_stencil_map => null()
    integer(kind=i_def) :: i,j
    integer(kind=i_def), allocatable :: cell_stencil(:)
    integer(kind=i_def) nfaces_re_h
    class(reference_element_type), pointer :: reference_element => null()
    !
    ! Initialise field and/or operator proxies
    !
    helmholtz_operator_proxy(1) = helmholtz_operator(1)%get_proxy()
    helmholtz_operator_proxy(2) = helmholtz_operator(2)%get_proxy()
    helmholtz_operator_proxy(3) = helmholtz_operator(3)%get_proxy()
    helmholtz_operator_proxy(4) = helmholtz_operator(4)%get_proxy()
    helmholtz_operator_proxy(5) = helmholtz_operator(5)%get_proxy()
    helmholtz_operator_proxy(6) = helmholtz_operator(6)%get_proxy()
    helmholtz_operator_proxy(7) = helmholtz_operator(7)%get_proxy()
    helmholtz_operator_proxy(8) = helmholtz_operator(8)%get_proxy()
    helmholtz_operator_proxy(9) = helmholtz_operator(9)%get_proxy()
    hb_lumped_inv_proxy = hb_lumped_inv%get_proxy()
    u_normalisation_proxy = u_normalisation%get_proxy()
    div_star_proxy = div_star%get_proxy()
    t_normalisation_proxy = t_normalisation%get_proxy()
    ptheta2v_proxy = ptheta2v%get_proxy()
    compound_div_proxy = compound_div%get_proxy()
    m3_exner_star_proxy = m3_exner_star%get_proxy()
    p3theta_proxy = p3theta%get_proxy()
    w2_mask_proxy = w2_mask%get_proxy()
    !
    ! Initialise number of layers
    !
    nlayers = helmholtz_operator_proxy(1)%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => helmholtz_operator_proxy(1)%vspace%get_mesh()
    !
    ! Initialise stencil dofmaps
    !
    hb_lumped_inv_stencil_map => hb_lumped_inv_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_depth)
    hb_lumped_inv_max_branch_length = stencil_depth + 1
    hb_lumped_inv_stencil_dofmap => hb_lumped_inv_stencil_map%get_whole_dofmap()
    hb_lumped_inv_stencil_sizes => hb_lumped_inv_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_w3 => helmholtz_operator_proxy(1)%vspace%get_whole_dofmap()
    map_w2 => hb_lumped_inv_proxy%vspace%get_whole_dofmap()
    map_wtheta => t_normalisation_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for w3
    !
    ndf_w3 = helmholtz_operator_proxy(1)%vspace%get_ndf()
    undf_w3 = helmholtz_operator_proxy(1)%vspace%get_undf()
    !
    ! Initialise number of DoFs for w2
    !
    ndf_w2 = hb_lumped_inv_proxy%vspace%get_ndf()
    undf_w2 = hb_lumped_inv_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for wtheta
    !
    ndf_wtheta = t_normalisation_proxy%vspace%get_ndf()
    undf_wtheta = t_normalisation_proxy%vspace%get_undf()
    !
    ! Call kernels and communication routines
    !
    if (hb_lumped_inv_proxy%is_dirty(depth=1)) then
      call hb_lumped_inv_proxy%halo_exchange(depth=1)
    end if
    if (u_normalisation_proxy%is_dirty(depth=1)) then
      call u_normalisation_proxy%halo_exchange(depth=1)
    end if
    if (t_normalisation_proxy%is_dirty(depth=1)) then
      call t_normalisation_proxy%halo_exchange(depth=1)
    end if
    if (w2_mask_proxy%is_dirty(depth=1)) then
      call w2_mask_proxy%halo_exchange(depth=1)
    end if
    !
    ! Get the reference element and query its properties
    !
    reference_element => mesh%get_reference_element()
    nfaces_re_h = reference_element%get_number_horizontal_faces()
    !
    ! Create cell stencil of the correct size
    stencil_size =  1 + nfaces_re_h*stencil_depth
    allocate( cell_stencil( stencil_size ) )

    !$omp parallel default(shared), private(cell, cell_stencil, i, j)
    !$omp do schedule(static)
    do cell=1,mesh%get_last_edge_cell()
      !
      ! Populate cell_stencil array used for operators
      ! (this is the id of each cell in the stencil)
      cell_stencil(:) = 0
      cell_stencil(1) = cell
      j=0
      do i = 1,nfaces_re_h
        if (mesh%get_cell_next(i, cell) /= 0)then
          j=j+1
          cell_stencil(j+1) = mesh%get_cell_next(i, cell)
        end if
      end do
      call helmholtz_operator_code(stencil_size,                     &
                                   cell_stencil, nlayers,            &
                                   helmholtz_operator_proxy(1)%data, &
                                   helmholtz_operator_proxy(2)%data, &
                                   helmholtz_operator_proxy(3)%data, &
                                   helmholtz_operator_proxy(4)%data, &
                                   helmholtz_operator_proxy(5)%data, &
                                   helmholtz_operator_proxy(6)%data, &
                                   helmholtz_operator_proxy(7)%data, &
                                   helmholtz_operator_proxy(8)%data, &
                                   helmholtz_operator_proxy(9)%data, &
                                   hb_lumped_inv_proxy%data, &
                                   hb_lumped_inv_stencil_sizes(:,cell), &
                                   hb_lumped_inv_max_branch_length,  &
                                   hb_lumped_inv_stencil_dofmap(:,:,:,cell), &
                                   u_normalisation_proxy%data, &
                                   div_star_proxy%ncell_3d, &
                                   div_star_proxy%local_stencil, &
                                   t_normalisation_proxy%data, &
                                   ptheta2v_proxy%ncell_3d, &
                                   ptheta2v_proxy%local_stencil, &
                                   compound_div_proxy%ncell_3d, &
                                   compound_div_proxy%local_stencil, &
                                   m3_exner_star_proxy%ncell_3d, &
                                   m3_exner_star_proxy%local_stencil, &
                                   p3theta_proxy%ncell_3d, &
                                   p3theta_proxy%local_stencil, &
                                   w2_mask_proxy%data, &
                                   ndf_w3, undf_w3, map_w3(:,cell), &
                                   ndf_w2, undf_w2, map_w2(:,cell), &
                                   ndf_wtheta, undf_wtheta, map_wtheta(:,cell))
    end do
    !$omp end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    !$omp master
    call helmholtz_operator_proxy(1)%set_dirty()
    call helmholtz_operator_proxy(2)%set_dirty()
    call helmholtz_operator_proxy(3)%set_dirty()
    call helmholtz_operator_proxy(4)%set_dirty()
    call helmholtz_operator_proxy(5)%set_dirty()
    call helmholtz_operator_proxy(6)%set_dirty()
    call helmholtz_operator_proxy(7)%set_dirty()
    call helmholtz_operator_proxy(8)%set_dirty()
    call helmholtz_operator_proxy(9)%set_dirty()
    !$omp end master
    !
    !$omp end parallel
    !
  end subroutine invoke_helmholtz_operator_kernel_type

  !----------------------------------------------------------------------------
  ! This requires a stencil of horizontal cells for the operators
  ! see PSyclone #1103: https://github.com/stfc/PSyclone/issues/1103
  ! The LFRic infrastructure for this will be introduced in #2532
  subroutine invoke_elim_helmholtz_operator_kernel_type(helmholtz_operator, hb_lumped_inv, stencil_depth, &
                                                   u_normalisation, div_star, &
                                                   m3_exner_star, Q32, &
                                                   w2_mask)
    use elim_helmholtz_operator_kernel_mod, only: elim_helmholtz_operator_code
    use mesh_mod, only: mesh_type
    use stencil_dofmap_mod, only: stencil_cross
    use stencil_dofmap_mod, only: stencil_dofmap_type
    use reference_element_mod, only: reference_element_type

    implicit none

    type(r_solver_field_type), intent(in) :: helmholtz_operator(9)
    type(r_solver_field_type), intent(in) :: hb_lumped_inv, u_normalisation, w2_mask
    type(r_solver_operator_type), intent(in) :: div_star, m3_exner_star, Q32
    integer(kind=i_def), intent(in) :: stencil_depth
    integer(kind=i_def) :: stencil_size
    integer(kind=i_def) cell
    integer(kind=i_def) nlayers
    type(r_solver_operator_proxy_type) div_star_proxy, m3_exner_star_proxy, Q32_proxy
    type(r_solver_field_proxy_type) helmholtz_operator_proxy(9)
    type(r_solver_field_proxy_type) hb_lumped_inv_proxy, u_normalisation_proxy, &
                           w2_mask_proxy
    integer(kind=i_def), pointer :: map_w2(:,:) => null(), map_w3(:,:) => null()
    integer(kind=i_def) ndf_w3, undf_w3, ndf_w2, undf_w2
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) hb_lumped_inv_stencil_size
    integer(kind=i_def), pointer :: hb_lumped_inv_stencil_dofmap(:,:,:) => null()
    type(stencil_dofmap_type), pointer :: hb_lumped_inv_stencil_map => null()
    integer(kind=i_def) :: i
    integer(kind=i_def), allocatable :: cell_stencil(:)
    integer(kind=i_def) nfaces_re_h
    class(reference_element_type), pointer :: reference_element => null()
    !
    ! Initialise field and/or operator proxies
    !
    helmholtz_operator_proxy(1) = helmholtz_operator(1)%get_proxy()
    helmholtz_operator_proxy(2) = helmholtz_operator(2)%get_proxy()
    helmholtz_operator_proxy(3) = helmholtz_operator(3)%get_proxy()
    helmholtz_operator_proxy(4) = helmholtz_operator(4)%get_proxy()
    helmholtz_operator_proxy(5) = helmholtz_operator(5)%get_proxy()
    helmholtz_operator_proxy(6) = helmholtz_operator(6)%get_proxy()
    helmholtz_operator_proxy(7) = helmholtz_operator(7)%get_proxy()
    helmholtz_operator_proxy(8) = helmholtz_operator(8)%get_proxy()
    helmholtz_operator_proxy(9) = helmholtz_operator(9)%get_proxy()
    hb_lumped_inv_proxy = hb_lumped_inv%get_proxy()
    u_normalisation_proxy = u_normalisation%get_proxy()
    div_star_proxy = div_star%get_proxy()
    m3_exner_star_proxy = m3_exner_star%get_proxy()
    Q32_proxy = Q32%get_proxy()
    w2_mask_proxy = w2_mask%get_proxy()
    !
    ! Initialise number of layers
    !
    nlayers = helmholtz_operator_proxy(1)%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => helmholtz_operator_proxy(1)%vspace%get_mesh()
    !
    ! Initialise stencil dofmaps
    !
    hb_lumped_inv_stencil_map => hb_lumped_inv_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS,stencil_depth)
    hb_lumped_inv_stencil_dofmap => hb_lumped_inv_stencil_map%get_whole_dofmap()
    hb_lumped_inv_stencil_size = hb_lumped_inv_stencil_map%get_size()
    !
    ! Look-up dofmaps for each function space
    !
    map_w3 => helmholtz_operator_proxy(1)%vspace%get_whole_dofmap()
    map_w2 => hb_lumped_inv_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for w3
    !
    ndf_w3 = helmholtz_operator_proxy(1)%vspace%get_ndf()
    undf_w3 = helmholtz_operator_proxy(1)%vspace%get_undf()
    !
    ! Initialise number of DoFs for w2
    !
    ndf_w2 = hb_lumped_inv_proxy%vspace%get_ndf()
    undf_w2 = hb_lumped_inv_proxy%vspace%get_undf()
    !
    ! Call kernels and communication routines
    !
    if (hb_lumped_inv_proxy%is_dirty(depth=1)) then
      call hb_lumped_inv_proxy%halo_exchange(depth=1)
    end if
    if (u_normalisation_proxy%is_dirty(depth=1)) then
      call u_normalisation_proxy%halo_exchange(depth=1)
    end if
    if (w2_mask_proxy%is_dirty(depth=1)) then
      call w2_mask_proxy%halo_exchange(depth=1)
    end if
    !
    ! Get the reference element and query its properties
    !
    reference_element => mesh%get_reference_element()
    nfaces_re_h = reference_element%get_number_horizontal_faces()
    !
    ! Create cell stencil of the correct size
    stencil_size =  1 + nfaces_re_h*stencil_depth
    allocate( cell_stencil( stencil_size ) )

    !$omp parallel default(shared), private(cell, cell_stencil, i)
    !$omp do schedule(static)
    do cell=1,mesh%get_last_edge_cell()
      !
      ! Populate cell_stencil array used for operators
      ! (this is the id of each cell in the stencil)
      cell_stencil(1) = cell
      do i = 1,nfaces_re_h
        cell_stencil(i+1) = mesh%get_cell_next(i, cell)
      end do
      call elim_helmholtz_operator_code(stencil_size,                &
                                   cell_stencil, nlayers,            &
                                   helmholtz_operator_proxy(1)%data, &
                                   helmholtz_operator_proxy(2)%data, &
                                   helmholtz_operator_proxy(3)%data, &
                                   helmholtz_operator_proxy(4)%data, &
                                   helmholtz_operator_proxy(5)%data, &
                                   helmholtz_operator_proxy(6)%data, &
                                   helmholtz_operator_proxy(7)%data, &
                                   helmholtz_operator_proxy(8)%data, &
                                   helmholtz_operator_proxy(9)%data, &
                                   hb_lumped_inv_proxy%data, &
                                   hb_lumped_inv_stencil_size, &
                                   hb_lumped_inv_stencil_dofmap(:,:,cell), &
                                   u_normalisation_proxy%data, &
                                   div_star_proxy%ncell_3d, &
                                   div_star_proxy%local_stencil, &
                                   m3_exner_star_proxy%ncell_3d, &
                                   m3_exner_star_proxy%local_stencil, &
                                   Q32_proxy%ncell_3d, &
                                   Q32_proxy%local_stencil, &
                                   w2_mask_proxy%data, &
                                   ndf_w3, undf_w3, map_w3(:,cell), &
                                   ndf_w2, undf_w2, map_w2(:,cell))
    end do
    !$omp end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    !$omp master
    call helmholtz_operator_proxy(1)%set_dirty()
    call helmholtz_operator_proxy(2)%set_dirty()
    call helmholtz_operator_proxy(3)%set_dirty()
    call helmholtz_operator_proxy(4)%set_dirty()
    call helmholtz_operator_proxy(5)%set_dirty()
    call helmholtz_operator_proxy(6)%set_dirty()
    call helmholtz_operator_proxy(7)%set_dirty()
    call helmholtz_operator_proxy(8)%set_dirty()
    call helmholtz_operator_proxy(9)%set_dirty()
    !$omp end master
    !
    !$omp end parallel
    !
  end subroutine invoke_elim_helmholtz_operator_kernel_type

!-------------------------------------------------------------------------------
!> Routine to perform higher-order reconstruction of a field on a fine mesh to
!! a coarse mesh. There is a field in this kernel that requires both a mesh
!! argument and a stencil argument, and PSyclone currently does not support
!! this. The issue to address this is #1983
  subroutine invoke_prolong_scalar_linear_kernel_type(target_field, source_field, stencil_extent)

    use sci_prolong_scalar_linear_kernel_mod, only: prolong_scalar_linear_kernel_code
    use mesh_map_mod,                     only: mesh_map_type
    use mesh_mod,                         only: mesh_type
    use stencil_dofmap_mod,               only: STENCIL_CROSS
    use stencil_dofmap_mod,               only: stencil_dofmap_type

    implicit none

    type(field_type),    intent(in) :: target_field, source_field
    integer(kind=i_def), intent(in) :: stencil_extent

    integer(kind=i_def) :: cell
    integer(kind=i_def) :: nlayers
    type(field_proxy_type) :: target_field_proxy, source_field_proxy
    integer(kind=i_def), pointer :: map_adspc1_target_field(:,:) => null(), map_adspc2_source_field(:,:) => null()
    integer(kind=i_def) :: ndf_adspc1_target_field, undf_adspc1_target_field, ndf_adspc2_source_field, undf_adspc2_source_field
    integer(kind=i_def) :: ncell_target_field, ncpc_target_field_source_field_x, ncpc_target_field_source_field_y
    integer(kind=i_def), pointer :: cell_map_source_field(:,:,:) => null()
    type(mesh_map_type), pointer :: mmap_target_field_source_field => null()
    type(mesh_type),     pointer :: mesh_target_field => null()
    type(mesh_type),     pointer :: mesh_source_field => null()
    integer(kind=i_def), pointer :: stencil_size(:) => null()
    integer(kind=i_def), pointer :: stencil_dofmap(:,:,:) => null()
    type(stencil_dofmap_type), pointer :: stencil_map => null()

    !
    ! Initialise field and/or operator proxies
    !
    target_field_proxy = target_field%get_proxy()
    source_field_proxy = source_field%get_proxy()
    !
    ! Initialise number of layers
    !
    nlayers = target_field_proxy%vspace%get_nlayers()
    !
    ! Look-up mesh objects and loop limits for inter-grid kernels
    !
    mesh_target_field => target_field_proxy%vspace%get_mesh()
    mesh_source_field => source_field_proxy%vspace%get_mesh()
    mmap_target_field_source_field => mesh_source_field%get_mesh_map(mesh_target_field)
    cell_map_source_field => mmap_target_field_source_field%get_whole_cell_map()
    ncell_target_field = mesh_target_field%get_last_halo_cell(depth=2)
    ncpc_target_field_source_field_x = mmap_target_field_source_field%get_ntarget_cells_per_source_x()
    ncpc_target_field_source_field_y = mmap_target_field_source_field%get_ntarget_cells_per_source_y()
    !
    ! Initialise stencil dofmaps
    !
    stencil_map => source_field_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS,stencil_extent)
    stencil_dofmap => stencil_map%get_whole_dofmap()
    stencil_size => stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_target_field => target_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_source_field => source_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for adspc1_target_field
    !
    ndf_adspc1_target_field = target_field_proxy%vspace%get_ndf()
    undf_adspc1_target_field = target_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_source_field
    !
    ndf_adspc2_source_field = source_field_proxy%vspace%get_ndf()
    undf_adspc2_source_field = source_field_proxy%vspace%get_undf()
    !
    ! Call kernels and communication routines
    !
    if (source_field_proxy%is_dirty(depth=stencil_extent)) then
      call source_field_proxy%halo_exchange(depth=stencil_extent)
    end if
    !
    do cell=1,mesh_source_field%get_last_edge_cell()
      !
      call prolong_scalar_linear_kernel_code(nlayers, cell_map_source_field(:,:,cell), ncpc_target_field_source_field_x, &
&ncpc_target_field_source_field_y, ncell_target_field, target_field_proxy%data, source_field_proxy%data, stencil_size(cell), &
stencil_dofmap(:,:,cell), ndf_adspc1_target_field, &
&undf_adspc1_target_field, map_adspc1_target_field, undf_adspc2_source_field, map_adspc2_source_field(:,cell))
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call target_field_proxy%set_dirty()
    !
    !
  end subroutine invoke_prolong_scalar_linear_kernel_type

!-------------------------------------------------------------------------------
!> Routine to perform injection of a multidata field on a fine mesh to
!> a coarse mesh. Intermesh kernels cannot currently take in integer arguments,
!> this has been raised as an issue https://github.com/stfc/PSyclone/issues/2504
  subroutine invoke_prolong_multidata_scalar_kernel_type(target_field, source_field, ndata)

    use sci_prolong_multidata_scalar_kernel_mod, only: prolong_multidata_scalar_kernel_code
    use mesh_map_mod,                        only: mesh_map_type
    use mesh_mod,                            only: mesh_type

    implicit none

    type(field_type),    intent(in) :: target_field, source_field
    integer(kind=i_def), intent(in) :: ndata

    integer(kind=i_def) :: cell
    integer(kind=i_def) :: nlayers
    type(field_proxy_type) :: target_field_proxy, source_field_proxy
    integer(kind=i_def), pointer :: map_adspc1_target_field(:,:) => null(), map_adspc2_source_field(:,:) => null()
    integer(kind=i_def) :: ndf_adspc1_target_field, undf_adspc1_target_field, ndf_adspc2_source_field, undf_adspc2_source_field
    integer(kind=i_def) :: ncell_target_field, ncpc_target_field_source_field_x, ncpc_target_field_source_field_y
    integer(kind=i_def), pointer :: cell_map_source_field(:,:,:) => null()
    type(mesh_map_type), pointer :: mmap_target_field_source_field => null()
    type(mesh_type),     pointer :: mesh_target_field => null()
    type(mesh_type),     pointer :: mesh_source_field => null()

    !
    ! Initialise field and/or operator proxies
    !
    target_field_proxy = target_field%get_proxy()
    source_field_proxy = source_field%get_proxy()
    !
    ! Initialise number of layers
    !
    nlayers = target_field_proxy%vspace%get_nlayers()
    !
    ! Look-up mesh objects and loop limits for inter-grid kernels
    !
    mesh_target_field => target_field_proxy%vspace%get_mesh()
    mesh_source_field => source_field_proxy%vspace%get_mesh()
    mmap_target_field_source_field => mesh_source_field%get_mesh_map(mesh_target_field)
    cell_map_source_field => mmap_target_field_source_field%get_whole_cell_map()
    ncell_target_field = mesh_target_field%get_last_halo_cell(depth=2)
    ncpc_target_field_source_field_x = mmap_target_field_source_field%get_ntarget_cells_per_source_x()
    ncpc_target_field_source_field_y = mmap_target_field_source_field%get_ntarget_cells_per_source_y()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_target_field => target_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_source_field => source_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of doFs for adspc1_target_field
    !
    ndf_adspc1_target_field = target_field_proxy%vspace%get_ndf()
    undf_adspc1_target_field = target_field_proxy%vspace%get_undf()
    !
    ! Initialise number of doFs for adspc2_source_field
    !
    ndf_adspc2_source_field = source_field_proxy%vspace%get_ndf()
    undf_adspc2_source_field = source_field_proxy%vspace%get_undf()
    !
    ! call kernels and communication routines
    !
    ! if (source_field_proxy%is_dirty(depth=1)) then
    !   call source_field_proxy%halo_exchange(depth=1)
    ! end if
    !
    do cell=1,mesh_source_field%get_last_edge_cell()
      !
      call prolong_multidata_scalar_kernel_code(nlayers, cell_map_source_field(:,:,cell), ncpc_target_field_source_field_x, &
&ncpc_target_field_source_field_y, ncell_target_field, &
&target_field_proxy%data, source_field_proxy%data, ndata, ndf_adspc1_target_field, &
&undf_adspc1_target_field, map_adspc1_target_field, undf_adspc2_source_field, map_adspc2_source_field(:,cell))
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call target_field_proxy%set_dirty()
    !
    !
  end subroutine invoke_prolong_multidata_scalar_kernel_type

!==============================================================================

!-------------------------------------------------------------------------------
!> Routine to perform averaging of a multidata field on a coarse mesh to
!> a fine mesh. Intermesh kernels cannot currently take in integer arguments,
!> this has been raised as an issue https://github.com/stfc/PSyclone/issues/2504
  subroutine invoke_restrict_multidata_scalar_kernel_type(target_field, source_field, ndata)
    use sci_restrict_multidata_scalar_kernel_mod, only: restrict_multidata_scalar_kernel_code
    use mesh_map_mod, only: mesh_map_type
    use mesh_mod, only: mesh_type
    type(field_type),    intent(in) :: target_field, source_field
    integer(kind=i_def), intent(in) :: ndata
    integer(kind=i_def) cell
    integer(kind=i_def) loop0_start, loop0_stop
    integer(kind=i_def) nlayers
    type(field_proxy_type) target_field_proxy, source_field_proxy
    integer(kind=i_def), pointer :: map_adspc1_target_field(:,:) => null(), map_adspc2_source_field(:,:) => null()
    integer(kind=i_def) ndf_adspc1_target_field, undf_adspc1_target_field, ndf_adspc2_source_field, undf_adspc2_source_field
    integer(kind=i_def) ncell_source_field, ncpc_source_field_target_field_x, ncpc_source_field_target_field_y
    integer(kind=i_def), pointer :: cell_map_target_field(:,:,:) => null()
    type(mesh_map_type), pointer :: mmap_source_field_target_field => null()
    integer(kind=i_def) max_halo_depth_mesh_target_field
    type(mesh_type), pointer :: mesh_target_field => null()
    integer(kind=i_def) max_halo_depth_mesh_source_field
    type(mesh_type), pointer :: mesh_source_field => null()
    !
    ! Initialise field and/or operator proxies
    !
    target_field_proxy = target_field%get_proxy()
    source_field_proxy = source_field%get_proxy()
    !
    ! Initialise number of layers
    !
    nlayers = target_field_proxy%vspace%get_nlayers()
    !
    ! Look-up mesh objects and loop limits for inter-grid kernels
    !
    mesh_source_field => source_field_proxy%vspace%get_mesh()
    max_halo_depth_mesh_source_field = mesh_source_field%get_halo_depth()
    mesh_target_field => target_field_proxy%vspace%get_mesh()
    max_halo_depth_mesh_target_field = mesh_target_field%get_halo_depth()
    mmap_source_field_target_field => mesh_target_field%get_mesh_map(mesh_source_field)
    cell_map_target_field => mmap_source_field_target_field%get_whole_cell_map()
    ncell_source_field = mesh_source_field%get_last_halo_cell(depth=2)
    ncpc_source_field_target_field_x = mmap_source_field_target_field%get_ntarget_cells_per_source_x()
    ncpc_source_field_target_field_y = mmap_source_field_target_field%get_ntarget_cells_per_source_y()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_target_field => target_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_source_field => source_field_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of doFs for adspc1_target_field
    !
    ndf_adspc1_target_field = target_field_proxy%vspace%get_ndf()
    undf_adspc1_target_field = target_field_proxy%vspace%get_undf()
    !
    ! Initialise number of doFs for adspc2_source_field
    !
    ndf_adspc2_source_field = source_field_proxy%vspace%get_ndf()
    undf_adspc2_source_field = source_field_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh_target_field%get_last_edge_cell()
    !
    ! call kernels and communication routines
    !
    do cell=loop0_start,loop0_stop
      !
      call restrict_multidata_scalar_kernel_code(nlayers, cell_map_target_field(:,:,cell), ncpc_source_field_target_field_x, &
&ncpc_source_field_target_field_y, ncell_source_field, target_field_proxy%data, &
&source_field_proxy%data, ndata, undf_adspc1_target_field, &
&map_adspc1_target_field(:,cell), ndf_adspc2_source_field, undf_adspc2_source_field, map_adspc2_source_field)
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call target_field_proxy%set_dirty()
    !
    !
  end subroutine invoke_restrict_multidata_scalar_kernel_type

  !-------------------------------------------------------------------------------
!> Routine to perform higher-order reconstruction of a multidata field on a fine mesh to
!> a coarse mesh. Intermesh kernels cannot currently take in integer arguments,
!> this has been raised as an issue https://github.com/stfc/PSyclone/issues/2504
  subroutine invoke_prolong_multidata_linear_kernel_type(target_field, source_field, source_mask, stencil_extent, ndata)

    use sci_prolong_multidata_linear_kernel_mod, only: prolong_multidata_linear_kernel_code
    use mesh_map_mod,                     only: mesh_map_type
    use mesh_mod,                         only: mesh_type
    use stencil_dofmap_mod,               only: STENCIL_CROSS
    use stencil_dofmap_mod,               only: stencil_dofmap_type

    implicit none

    type(field_type),    intent(in) :: target_field, source_field, source_mask
    integer(kind=i_def), intent(in) :: stencil_extent
    integer(kind=i_def), intent(in) :: ndata

    integer(kind=i_def) :: cell
    integer(kind=i_def) :: nlayers
    type(field_proxy_type) :: target_field_proxy, source_field_proxy, source_mask_proxy
    integer(kind=i_def), pointer :: map_adspc1_target_field(:,:) => null(), map_adspc2_source_field(:,:) => null(), &
                                    map_adspc3_source_mask(:,:) => null()
    integer(kind=i_def) :: ndf_adspc1_target_field, undf_adspc1_target_field, ndf_adspc2_source_field, undf_adspc2_source_field, &
                           ndf_adspc3_source_mask, undf_adspc3_source_mask
    integer(kind=i_def) :: ncell_target_field, ncpc_target_field_source_field_x, ncpc_target_field_source_field_y
    integer(kind=i_def), pointer :: cell_map_source_field(:,:,:) => null()
    type(mesh_map_type), pointer :: mmap_target_field_source_field => null()
    type(mesh_type),     pointer :: mesh_target_field => null()
    type(mesh_type),     pointer :: mesh_source_field => null()
    integer(kind=i_def), pointer :: stencil_size(:) => null()
    integer(kind=i_def), pointer :: stencil_dofmap(:,:,:) => null()
    type(stencil_dofmap_type), pointer :: stencil_map => null()
    integer(kind=i_def), pointer :: stencil_size_mask(:) => null()
    integer(kind=i_def), pointer :: stencil_dofmap_mask(:,:,:) => null()
    type(stencil_dofmap_type), pointer :: stencil_map_mask => null()

    !
    ! Initialise field and/or operator proxies
    !
    target_field_proxy = target_field%get_proxy()
    source_field_proxy = source_field%get_proxy()
    source_mask_proxy = source_mask%get_proxy()
    !
    ! Initialise number of layers
    !
    nlayers = target_field_proxy%vspace%get_nlayers()
    !
    ! Look-up mesh objects and loop limits for inter-grid kernels
    !
    mesh_target_field => target_field_proxy%vspace%get_mesh()
    mesh_source_field => source_field_proxy%vspace%get_mesh()
    mmap_target_field_source_field => mesh_source_field%get_mesh_map(mesh_target_field)
    cell_map_source_field => mmap_target_field_source_field%get_whole_cell_map()
    ncell_target_field = mesh_target_field%get_last_halo_cell(depth=2)
    ncpc_target_field_source_field_x = mmap_target_field_source_field%get_ntarget_cells_per_source_x()
    ncpc_target_field_source_field_y = mmap_target_field_source_field%get_ntarget_cells_per_source_y()
    !
    ! Initialise stencil dofmaps
    !
    stencil_map => source_field_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS,stencil_extent)
    stencil_dofmap => stencil_map%get_whole_dofmap()
    stencil_size => stencil_map%get_stencil_sizes()
    stencil_map_mask => source_mask_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS,stencil_extent)
    stencil_size_mask => stencil_map_mask%get_stencil_sizes()
    stencil_dofmap_mask => stencil_map_mask%get_whole_dofmap()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_target_field => target_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_source_field => source_field_proxy%vspace%get_whole_dofmap()
    map_adspc3_source_mask => source_mask_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of doFs for adspc1_target_field
    !
    ndf_adspc1_target_field = target_field_proxy%vspace%get_ndf()
    undf_adspc1_target_field = target_field_proxy%vspace%get_undf()
    !
    ! Initialise number of doFs for adspc2_source_field
    !
    ndf_adspc2_source_field = source_field_proxy%vspace%get_ndf()
    undf_adspc2_source_field = source_field_proxy%vspace%get_undf()
    !
    ! Initialise number of doFs for adspc2_source_field
    !
    ndf_adspc3_source_mask = source_mask_proxy%vspace%get_ndf()
    undf_adspc3_source_mask = source_mask_proxy%vspace%get_undf()
    ! call kernels and communication routines
    !
    if (source_field_proxy%is_dirty(depth=stencil_extent)) then
      call source_field_proxy%halo_exchange(depth=stencil_extent)
    end if
    !
    if (source_mask_proxy%is_dirty(depth=stencil_extent)) then
      call source_mask_proxy%halo_exchange(depth=stencil_extent)
    end if
    !
    do cell=1,mesh_source_field%get_last_edge_cell()
      !
      call prolong_multidata_linear_kernel_code(nlayers, cell_map_source_field(:,:,cell), &
ncpc_target_field_source_field_x, &
&ncpc_target_field_source_field_y, ncell_target_field, target_field_proxy%data, &
source_field_proxy%data, stencil_size(cell), &
stencil_dofmap(:,:,cell), source_mask_proxy%data, &
stencil_size_mask(cell), stencil_dofmap_mask(:,:,cell), &
ndata, ndf_adspc1_target_field, &
&undf_adspc1_target_field, map_adspc1_target_field, undf_adspc2_source_field, &
map_adspc2_source_field(:,cell), ndf_adspc3_source_mask, &
&undf_adspc3_source_mask, map_adspc3_source_mask)
    end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    call target_field_proxy%set_dirty()
    !
    !
  end subroutine invoke_prolong_multidata_linear_kernel_type

end module psykal_lite_mod
