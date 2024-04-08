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
  !> io_mod uses this routine. However, because io_mod is not a algorithm its currently
  !> not clear if it should call into PSy. #1253 will address this point and remove
  !> once decided.
  subroutine invoke_nodal_coordinates_kernel(nodal_coords, chi )
    use nodal_coordinates_kernel_mod, only: nodal_coordinates_code
    use mesh_mod,                     only: mesh_type
    implicit none

    type(field_type), intent(inout)      :: nodal_coords(3)
    type(field_type), intent(in)         :: chi(3)

    type(field_proxy_type) :: x_p(3), chi_p(3)

    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_x
    integer                 :: undf_chi, undf_x
    integer                 :: dim_chi
    integer                 :: df_x, df_chi

    integer, pointer        :: map_chi(:) => null()
    integer, pointer        :: map_x(:) => null()
    real(kind=r_def), pointer :: nodes_x(:,:) => null()

    real(kind=r_def), allocatable  :: basis_chi(:,:,:)
    integer :: i
    type(mesh_type), pointer :: mesh => null()

    do i = 1,3
      x_p(i)   = nodal_coords(i)%get_proxy()
      chi_p(i) = chi(i)%get_proxy()
    end do

    nlayers = x_p(1)%vspace%get_nlayers()

    ndf_x  = x_p(1)%vspace%get_ndf( )
    undf_x = x_p(1)%vspace%get_undf()
    nodes_x => x_p(1)%vspace%get_nodes()

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    dim_chi = chi_p(1)%vspace%get_dim_space( )

    ! Evaluate the basis function
    allocate(basis_chi(dim_chi, ndf_chi, ndf_x))
    do df_x = 1, ndf_x
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
      end do
    end do

    if (chi_p(1)%is_dirty(depth=1)) then
       call chi_p(1)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(2)%is_dirty(depth=1)) then
       call chi_p(2)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(3)%is_dirty(depth=1)) then
       call chi_p(3)%halo_exchange(depth=1)
    end if
    mesh => x_p(1)%vspace%get_mesh()

    do cell = 1, mesh%get_last_halo_cell(1)
       map_x   => x_p(1)%vspace%get_cell_dofmap( cell )
       map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
       call nodal_coordinates_code(nlayers, &
                                   x_p(1)%data, &
                                   x_p(2)%data, &
                                   x_p(3)%data, &
                                   chi_p(1)%data, &
                                   chi_p(2)%data, &
                                   chi_p(3)%data, &
                                   ndf_x, undf_x, map_x, &
                                   ndf_chi, undf_chi, map_chi, &
                                   basis_chi &
                                  )
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

    deallocate(basis_chi)
  end subroutine invoke_nodal_coordinates_kernel

  !-------------------------------------------------------------------------------
  !> io_mod uses this routine. However, because io_mod is not a algorithm its currently
  !> not clear if it should call into PSy. #1253 will address this point and remove
  !> once decided.
  subroutine invoke_nodal_xyz_coordinates_kernel(nodal_coords, chi, panel_id)
    use nodal_xyz_coordinates_kernel_mod, only: nodal_xyz_coordinates_code
    use mesh_mod,                     only: mesh_type
    implicit none

    type(field_type), intent(inout)      :: nodal_coords(3)
    type(field_type), intent(in)         :: chi(3)
    type(field_type), intent(in)         :: panel_id

    type(field_proxy_type) :: x_p(3), chi_p(3), panel_id_proxy

    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_x, ndf_pid
    integer                 :: undf_chi, undf_x, undf_pid
    integer                 :: dim_chi
    integer                 :: df_x, df_chi

    integer, pointer        :: map_pid(:) => null()
    integer, pointer        :: map_chi(:) => null()
    integer, pointer        :: map_x(:) => null()
    real(kind=r_def), pointer :: nodes_x(:,:) => null()

    real(kind=r_def), allocatable  :: basis_chi(:,:,:)
    integer :: i
    type(mesh_type), pointer :: mesh => null()

    do i = 1,3
      x_p(i)   = nodal_coords(i)%get_proxy()
      chi_p(i) = chi(i)%get_proxy()
    end do

    panel_id_proxy = panel_id%get_proxy()

    nlayers = x_p(1)%vspace%get_nlayers()

    ndf_x  = x_p(1)%vspace%get_ndf()
    undf_x = x_p(1)%vspace%get_undf()
    nodes_x => x_p(1)%vspace%get_nodes()

    ndf_chi  = chi_p(1)%vspace%get_ndf()
    undf_chi = chi_p(1)%vspace%get_undf()
    dim_chi = chi_p(1)%vspace%get_dim_space()

    ndf_pid  = panel_id_proxy%vspace%get_ndf()
    undf_pid = panel_id_proxy%vspace%get_undf()

    ! Evaluate the basis function
    allocate(basis_chi(dim_chi, ndf_chi, ndf_x))
    do df_x = 1, ndf_x
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
      end do
    end do

    if (chi_p(1)%is_dirty(depth=1)) then
       call chi_p(1)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(2)%is_dirty(depth=1)) then
       call chi_p(2)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(3)%is_dirty(depth=1)) then
       call chi_p(3)%halo_exchange(depth=1)
    end if
    if (panel_id_proxy%is_dirty(depth=1)) then
      call panel_id_proxy%halo_exchange(depth=1)
    end if

    mesh => x_p(1)%vspace%get_mesh()

    do cell = 1, mesh%get_last_halo_cell(1)
       map_x   => x_p(1)%vspace%get_cell_dofmap( cell )
       map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
       map_pid => panel_id_proxy%vspace%get_cell_dofmap( cell )
       call nodal_xyz_coordinates_code(nlayers, &
                                       x_p(1)%data, &
                                       x_p(2)%data, &
                                       x_p(3)%data, &
                                       chi_p(1)%data, &
                                       chi_p(2)%data, &
                                       chi_p(3)%data, &
                                       panel_id_proxy%data, &
                                       ndf_x, undf_x, map_x, &
                                       ndf_chi, undf_chi, map_chi, &
                                       basis_chi, &
                                       ndf_pid, undf_pid, map_pid &
                                      )
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

    deallocate(basis_chi)
  end subroutine invoke_nodal_xyz_coordinates_kernel

!-------------------------------------------------------------------------------
  subroutine invoke_compute_dof_level_kernel(level)

  use compute_dof_level_kernel_mod, only: compute_dof_level_code
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

    !>@brief Remap a scalar field from the standard cubed sphere mesh onto an extended
    !!       mesh
    !!       This routine loops only over halo cells (and always to the full
    !!       depth of the halo). It is not clear if this functionality will ever
    !!       be supported by psyclone or if it should be treated as a special
    !!       case (as with computation of coordinate fields). Issue 2300 has been
    !!       opened to investigate this.
    subroutine invoke_init_remap_on_extended_mesh_kernel_type(remap_weights, remap_indices, &
                                                              chi_ext, chi, chi_stencil_depth, &
                                                              panel_id, pid_stencil_depth, &
                                                              linear_remap, ndata)

      use init_remap_on_extended_mesh_kernel_mod, only: init_remap_on_extended_mesh_code
      use function_space_mod,                     only: BASIS, DIFF_BASIS
      use mesh_mod,                               only: mesh_type
      use stencil_2D_dofmap_mod,                  only: stencil_2D_dofmap_type, STENCIL_2D_CROSS

      implicit none

      type(r_tran_field_type), intent(in) :: remap_weights
      type(field_type), intent(in) :: chi_ext(3), chi(3), panel_id
      type(integer_field_type), intent(in) :: remap_indices
      logical(kind=l_def), intent(in) :: linear_remap
      integer(kind=i_def), intent(in) :: ndata
      integer(kind=i_def), intent(in) :: chi_stencil_depth, pid_stencil_depth
      integer(kind=i_def) :: cell
      integer(kind=i_def) :: df_nodal, df_wchi
      real(kind=r_def), allocatable :: basis_wchi(:,:,:)
      integer(kind=i_def) :: dim_wchi
      real(kind=r_def), pointer :: nodes_remap(:,:) => null()
      integer(kind=i_def) :: nlayers
      type(r_tran_field_proxy_type) :: remap_weights_proxy
      type(integer_field_proxy_type) :: remap_indices_proxy
      type(field_proxy_type) :: chi_ext_proxy(3), chi_proxy(3), panel_id_proxy
      integer(kind=i_def), pointer :: map_remap(:,:) => null(), map_panel_id(:,:) => null(), map_wchi(:,:) => null()
      integer(kind=i_def) :: ndf_remap, undf_remap, ndf_wchi, undf_wchi, ndf_panel_id, undf_panel_id
      type(mesh_type), pointer :: mesh => null()
      type(stencil_2d_dofmap_type), pointer :: stencil_map => null()
      integer(kind=i_def), pointer :: wchi_stencil_size(:,:) => null()
      integer(kind=i_def), pointer :: wchi_stencil_dofmap(:,:,:,:) => null()
      integer(kind=i_def)          :: wchi_stencil_max_branch_length
      integer(kind=i_def), pointer :: pid_stencil_size(:,:) => null()
      integer(kind=i_def), pointer :: pid_stencil_dofmap(:,:,:,:) => null()
      integer(kind=i_def)          :: pid_stencil_max_branch_length
      integer(kind=i_def)          :: cell_start, cell_end

      ! Initialise field and/or operator proxies
      remap_weights_proxy = remap_weights%get_proxy()
      remap_indices_proxy = remap_indices%get_proxy()
      chi_ext_proxy(1) = chi_ext(1)%get_proxy()
      chi_ext_proxy(2) = chi_ext(2)%get_proxy()
      chi_ext_proxy(3) = chi_ext(3)%get_proxy()
      chi_proxy(1) = chi(1)%get_proxy()
      chi_proxy(2) = chi(2)%get_proxy()
      chi_proxy(3) = chi(3)%get_proxy()
      panel_id_proxy = panel_id%get_proxy()

      ! Initialise number of layers
      nlayers = remap_weights_proxy%vspace%get_nlayers()

      ! Create a mesh object
      mesh => remap_weights_proxy%vspace%get_mesh()

      ! Initialise stencil dofmaps
      stencil_map => chi_ext_proxy(1)%vspace%get_stencil_2d_dofmap(STENCIL_2D_CROSS, chi_stencil_depth)
      wchi_stencil_max_branch_length = chi_stencil_depth + 1_i_def
      wchi_stencil_dofmap => stencil_map%get_whole_dofmap()
      wchi_stencil_size => stencil_map%get_stencil_sizes()

      stencil_map => panel_id_proxy%vspace%get_stencil_2d_dofmap(STENCIL_2D_CROSS, pid_stencil_depth)
      pid_stencil_max_branch_length = pid_stencil_depth + 1_i_def
      pid_stencil_dofmap => stencil_map%get_whole_dofmap()
      pid_stencil_size => stencil_map%get_stencil_sizes()

      ! Look-up dofmaps for each function space
      map_remap => remap_weights_proxy%vspace%get_whole_dofmap()
      map_wchi => chi_ext_proxy(1)%vspace%get_whole_dofmap()
      map_panel_id => panel_id_proxy%vspace%get_whole_dofmap()

      ! Initialise number of DoFs for remap
      ndf_remap = remap_weights_proxy%vspace%get_ndf()
      undf_remap = remap_weights_proxy%vspace%get_undf()

      ! Initialise number of DoFs for wchi
      ndf_wchi = chi_ext_proxy(1)%vspace%get_ndf()
      undf_wchi = chi_ext_proxy(1)%vspace%get_undf()

      ! Initialise number of DoFs for panel_id
      ndf_panel_id = panel_id_proxy%vspace%get_ndf()
      undf_panel_id = panel_id_proxy%vspace%get_undf()

      ! Initialise evaluator-related quantities for the target function spaces
      nodes_remap => remap_weights_proxy%vspace%get_nodes()

      ! Allocate basis/diff-basis arrays
      dim_wchi = chi_ext_proxy(1)%vspace%get_dim_space()
      allocate (basis_wchi(dim_wchi, ndf_wchi, ndf_remap))

      ! Compute basis/diff-basis arrays
      do df_nodal = 1,ndf_remap
        do df_wchi = 1,ndf_wchi
          basis_wchi(:,df_wchi,df_nodal) = chi_ext_proxy(1)%vspace%call_function(BASIS,df_wchi,nodes_remap(:,df_nodal))
        end do
      end do

      ! Call kernels and communication routines
      if (panel_id_proxy%is_dirty(depth=mesh%get_halo_depth())) THEN
        call panel_id_proxy%halo_exchange(depth=mesh%get_halo_depth())
      end if

      cell_start = mesh%get_last_edge_cell() + 1
      cell_end   = mesh%get_last_halo_cell(mesh%get_halo_depth())

      !$omp parallel default(shared), private(cell)
      !$omp do schedule(static)
      do cell = cell_start, cell_end
        call init_remap_on_extended_mesh_code(nlayers, &
                                         remap_weights_proxy%data, &
                                         remap_indices_proxy%data, &
                                         chi_ext_proxy(1)%data, &
                                         chi_ext_proxy(2)%data, &
                                         chi_ext_proxy(3)%data, &
                                         chi_proxy(1)%data, &
                                         chi_proxy(2)%data, &
                                         chi_proxy(3)%data, &
                                         wchi_stencil_size(:,cell), &
                                         wchi_stencil_dofmap(:,:,:,cell), &
                                         wchi_stencil_max_branch_length, &
                                         panel_id_proxy%data, &
                                         pid_stencil_size(:,cell), &
                                         pid_stencil_dofmap(:,:,:,cell), &
                                         pid_stencil_max_branch_length, &
                                         linear_remap, &
                                         ndata, &
                                         ndf_remap, &
                                         undf_remap, &
                                         map_remap(:,cell), &
                                         ndf_wchi, &
                                         undf_wchi,&
                                         map_wchi(:,cell), &
                                         basis_wchi, &
                                         ndf_panel_id, &
                                         undf_panel_id, map_panel_id(:,cell))
      end do
      !$omp end do

      ! Set halos dirty/clean for fields modified in the above loop
      !$omp master
      call remap_weights_proxy%set_clean(mesh%get_halo_depth())
      call remap_indices_proxy%set_clean(mesh%get_halo_depth())
      !$omp end master
      !
      !$omp end parallel

      ! Deallocate basis arrays
      deallocate (basis_wchi)

    end subroutine invoke_init_remap_on_extended_mesh_kernel_type


   !>@brief Remap a scalar field from the standard cubed sphere mesh onto an extended
    !!       mesh
    !!       This routine loops only over halo cells (and always to the full
    !!       depth of the halo). It is not clear if this functionality will ever
    !!       be supported by psyclone or if it should be treated as a special
    !!       case (as with computation of coordinate fields). Issue 2300 has been
    !!       opened to investigate this.
    subroutine invoke_remap_on_extended_mesh_kernel_type(remap_field, field, stencil_depth, &
                                                         remap_weights, remap_indices, &
                                                         panel_id, &
                                                         ndata, &
                                                         monotone, enforce_minvalue, minvalue, &
                                                         halo_compute_depth )

      use remap_on_extended_mesh_kernel_mod, only: remap_on_extended_mesh_code
      use mesh_mod,                          only: mesh_type
      use stencil_2D_dofmap_mod,             only: stencil_2D_dofmap_type, STENCIL_2D_CROSS
      implicit none

      type(r_tran_field_type), intent(in) :: remap_field, field, remap_weights
      type(integer_field_type), intent(in) :: remap_indices
      type(field_type), intent(in) :: panel_id
      integer(kind=i_def), intent(in) :: ndata
      logical(kind=l_def), intent(in) :: monotone
      logical(kind=l_def), intent(in) :: enforce_minvalue
      real(kind=r_tran),   intent(in) :: minvalue
      integer(kind=i_def), intent(in) :: halo_compute_depth
      integer(kind=i_def) :: cell, stencil_depth
      integer(kind=i_def) :: nlayers
      type(r_tran_field_proxy_type) :: remap_field_proxy, field_proxy, remap_weights_proxy
      type(integer_field_proxy_type) :: remap_indices_proxy
      type(field_proxy_type) :: panel_id_proxy
      integer(kind=i_def), pointer :: map_remap_field(:,:) => null(), map_panel_id(:,:) => null(), map_remap(:,:) => null()
      integer(kind=i_def) :: ndf_remap_field, undf_remap_field, ndf_remap, undf_remap, ndf_panel_id, undf_panel_id
      type(mesh_type), pointer :: mesh => null()
      type(stencil_2d_dofmap_type), pointer :: stencil_map => null()
      integer(kind=i_def), pointer :: stencil_size(:,:) => null()
      integer(kind=i_def), pointer :: stencil_dofmap(:,:,:,:) => null()
      integer(kind=i_def)          :: stencil_max_branch_length
      integer(kind=i_def)          :: cell_start, cell_end

      ! Initialise field and/or operator proxies
      remap_field_proxy = remap_field%get_proxy()
      field_proxy = field%get_proxy()
      remap_weights_proxy = remap_weights%get_proxy()
      remap_indices_proxy = remap_indices%get_proxy()
      panel_id_proxy = panel_id%get_proxy()

      ! Initialise number of layers
      nlayers = remap_field_proxy%vspace%get_nlayers()

      ! Create a mesh object
      mesh => remap_field_proxy%vspace%get_mesh()

      ! Initialise stencil dofmaps
      stencil_map => field_proxy%vspace%get_stencil_2d_dofmap(STENCIL_2D_CROSS, stencil_depth)
      stencil_max_branch_length = stencil_depth + 1_i_def
      stencil_dofmap => stencil_map%get_whole_dofmap()
      stencil_size => stencil_map%get_stencil_sizes()

      ! Look-up dofmaps for each function space
      map_remap_field => remap_field_proxy%vspace%get_whole_dofmap()
      map_remap => remap_weights_proxy%vspace%get_whole_dofmap()
      map_panel_id => panel_id_proxy%vspace%get_whole_dofmap()

      ! Initialise number of DoFs for remap_field
      ndf_remap_field = remap_field_proxy%vspace%get_ndf()
      undf_remap_field = remap_field_proxy%vspace%get_undf()

      ! Initialise number of DoFs for interpolation fields
      ndf_remap = remap_weights_proxy%vspace%get_ndf()
      undf_remap = remap_weights_proxy%vspace%get_undf()

      ! Initialise number of DoFs for panel_id
      ndf_panel_id = panel_id_proxy%vspace%get_ndf()
      undf_panel_id = panel_id_proxy%vspace%get_undf()

      ! Call kernels and communication routines
      if (field_proxy%is_dirty(depth=mesh%get_halo_depth())) THEN
        call field_proxy%halo_exchange(depth=mesh%get_halo_depth())
      end if
      if (panel_id_proxy%is_dirty(depth=halo_compute_depth)) THEN
        call panel_id_proxy%halo_exchange(depth=halo_compute_depth)
      end if
      cell_start = mesh%get_last_edge_cell() + 1
      cell_end   = mesh%get_last_halo_cell(halo_compute_depth)

      !$omp parallel default(shared), private(cell)
      !$omp do schedule(static)
      do cell = cell_start, cell_end
        call remap_on_extended_mesh_code(nlayers, &
                                         remap_field_proxy%data, &
                                         field_proxy%data, &
                                         stencil_size(:,cell), &
                                         stencil_dofmap(:,:,:,cell), &
                                         stencil_max_branch_length, &
                                         remap_weights_proxy%data, &
                                         remap_indices_proxy%data, &
                                         panel_id_proxy%data, &
                                         ndata, &
                                         monotone, &
                                         enforce_minvalue, &
                                         minvalue, &
                                         ndf_remap_field, &
                                         undf_remap_field, &
                                         map_remap_field(:,cell), &
                                         ndf_remap, &
                                         undf_remap, &
                                         map_remap(:,cell), &
                                         ndf_panel_id, &
                                         undf_panel_id, map_panel_id(:,cell))
      end do
      !$omp end do

      ! Set halos dirty/clean for fields modified in the above loop
      !$omp master
      call remap_field_proxy%set_clean(halo_compute_depth)
      !$omp end master
      !
      !$omp end parallel
    end subroutine invoke_remap_on_extended_mesh_kernel_type


    ! Psyclone does not currently have native support for builtins with mixed
    ! precision, this will be addressed in https://github.com/stfc/PSyclone/issues/1786
    ! Copy a r_tran_field_type to a field_type
    subroutine invoke_copy_rtran_to_rdef(rdef_field, field)

      use omp_lib,            only: omp_get_thread_num
      use omp_lib,            only: omp_get_max_threads
      use mesh_mod,           only: mesh_type
      use r_tran_field_mod,   only: r_tran_field_type, r_tran_field_proxy_type
      use field_mod,          only: field_type, field_proxy_type

      implicit none

      type(field_type),        intent(inout) :: rdef_field
      type(r_tran_field_type), intent(in)    :: field

      integer(kind=i_def)             :: df
      integer(kind=i_def)             :: loop0_start, loop0_stop
      type(r_tran_field_proxy_type)   :: field_proxy
      type(field_proxy_type)          :: rdef_field_proxy
      integer(kind=i_def)             :: max_halo_depth_mesh
      type(mesh_type), pointer        :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      rdef_field_proxy = rdef_field%get_proxy()
      field_proxy = field%get_proxy()
      !
      ! Create a mesh object
      !
      mesh => rdef_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      IF (field_proxy%is_dirty(depth=1)) THEN
        ! only copy the owned dofs
        loop0_stop = rdef_field_proxy%vspace%get_last_dof_annexed()
      ELSE
        ! copy the 1st halo row as well
        loop0_stop = rdef_field_proxy%vspace%get_last_dof_halo(1)
      END IF
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(df)
      !$omp do schedule(static)
      DO df=loop0_start,loop0_stop
        rdef_field_proxy%data(df) = real(field_proxy%data(df), r_def)
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL rdef_field_proxy%set_dirty()
      IF (.not. field_proxy%is_dirty(depth=1)) THEN
        CALL rdef_field_proxy%set_clean(1)
      END IF
      !
    end subroutine invoke_copy_rtran_to_rdef


    ! Psyclone does not currently have native support for builtins with mixed
    ! precision, this will be addressed in https://github.com/stfc/PSyclone/issues/1786
    ! Copy a field_type to a r_tran_field_type
    subroutine invoke_copy_to_rtran(r_tran_field, field)

      use omp_lib,            only: omp_get_thread_num
      use omp_lib,            only: omp_get_max_threads
      use mesh_mod,           only: mesh_type
      use r_tran_field_mod,   only: r_tran_field_type, r_tran_field_proxy_type
      use field_mod,          only: field_type, field_proxy_type

      implicit none

      type(r_tran_field_type), intent(inout) :: r_tran_field
      type(field_type),        intent(in)    :: field

      integer(kind=i_def)             :: df
      integer(kind=i_def)             :: loop0_start, loop0_stop
      type(r_tran_field_proxy_type)   :: r_tran_field_proxy
      type(field_proxy_type)          :: field_proxy
      integer(kind=i_def)             :: max_halo_depth_mesh
      type(mesh_type), pointer        :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      r_tran_field_proxy = r_tran_field%get_proxy()
      field_proxy = field%get_proxy()
      !
      ! Create a mesh object
      !
      mesh => r_tran_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      IF (field_proxy%is_dirty(depth=1)) THEN
        ! only copy the owned dofs
        loop0_stop = r_tran_field_proxy%vspace%get_last_dof_annexed()
      ELSE
        ! copy the 1st halo row as well
        loop0_stop = r_tran_field_proxy%vspace%get_last_dof_halo(1)
      END IF
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(df)
      !$omp do schedule(static)
      DO df=loop0_start,loop0_stop
        r_tran_field_proxy%data(df) = real(field_proxy%data(df), r_tran)
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL r_tran_field_proxy%set_dirty()
      IF (.not. field_proxy%is_dirty(depth=1)) THEN
        CALL r_tran_field_proxy%set_clean(1)
      END IF
      !
    end subroutine invoke_copy_to_rtran

    ! Psyclone does not currently recognize r_tran_operator_type,
    ! this will be addressed in https://github.com/stfc/PSyclone/issues/1998
    SUBROUTINE invoke_operator_setval_x_rtran(op_r_tran, op)
      USE operator_algebra_kernel_mod, ONLY: operator_setval_x_kernel_code
      USE mesh_mod, ONLY: mesh_type
      implicit none
      TYPE(r_tran_operator_type), intent(in) :: op_r_tran
      TYPE(operator_type), intent(in) :: op
      INTEGER colour
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop1_start
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      TYPE(operator_proxy_type) op_proxy
      TYPE(r_tran_operator_proxy_type) op_r_tran_proxy
      INTEGER(KIND=i_def) ndf_aspc1_op_r_tran, ndf_aspc2_op_r_tran
      INTEGER(KIND=i_def), allocatable :: last_halo_cell_all_colours(:,:)
      INTEGER(KIND=i_def) ncolour
      INTEGER(KIND=i_def), pointer :: cmap(:,:)
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      op_r_tran_proxy = op_r_tran%get_proxy()
      op_proxy = op%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = op_r_tran_proxy%fs_from%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => op_r_tran_proxy%fs_from%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Get the colourmap
      !
      ncolour = mesh%get_ncolours()
      cmap => mesh%get_colour_map()
      !
      ! Initialise number of DoFs for aspc1_op_r_tran
      !
      ndf_aspc1_op_r_tran = op_r_tran_proxy%fs_to%get_ndf()
      !
      ! Initialise number of DoFs for aspc2_op_r_tran
      !
      ndf_aspc2_op_r_tran = op_r_tran_proxy%fs_from%get_ndf()
      !
      ! Initialise mesh propertie
      !
      last_halo_cell_all_colours = mesh%get_last_halo_cell_all_colours()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = ncolour
      loop1_start = 1
      !
      ! Call kernels and communication routines
      !
      DO colour=loop0_start,loop0_stop
        !$omp parallel default(shared), private(cell)
        !$omp do schedule(static)
        DO cell=loop1_start,last_halo_cell_all_colours(colour,1)
          !
          CALL operator_setval_x_kernel_code(cmap(colour, cell), nlayers, op_r_tran_proxy%ncell_3d, &
&op_r_tran_proxy%local_stencil, op_proxy%ncell_3d, op_proxy%local_stencil, ndf_aspc1_op_r_tran, ndf_aspc2_op_r_tran)
        END DO
        !$omp end do
        !$omp end parallel
      END DO
      !
    end SUBROUTINE invoke_operator_setval_x_rtran

    ! Psyclone does not currently recognize r_tran_operator_type,
    ! this will be addressed in https://github.com/stfc/PSyclone/issues/1998
    SUBROUTINE invoke_dg_matrix_vector_rtran(lhs, x, matrix)
      USE dg_matrix_vector_kernel_mod, ONLY: dg_matrix_vector_code
      USE mesh_mod, ONLY: mesh_type
      implicit none
      TYPE(r_tran_field_type), intent(in) :: lhs, x
      TYPE(r_tran_operator_type), intent(in) :: matrix
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      TYPE(r_tran_operator_proxy_type) matrix_proxy
      TYPE(r_tran_field_proxy_type) lhs_proxy, x_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_lhs(:,:) => null(), map_aspc1_x(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_lhs, undf_adspc1_lhs, ndf_aspc1_x, undf_aspc1_x, &
&ndf_aspc1_lhs, undf_aspc1_lhs
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      lhs_proxy = lhs%get_proxy()
      x_proxy = x%get_proxy()
      matrix_proxy = matrix%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = lhs_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => lhs_proxy%vspace%get_mesh()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_lhs => lhs_proxy%vspace%get_whole_dofmap()
      map_aspc1_x => x_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_lhs
      !
      ndf_adspc1_lhs = lhs_proxy%vspace%get_ndf()
      undf_adspc1_lhs = lhs_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for aspc1_x
      !
      ndf_aspc1_x = x_proxy%vspace%get_ndf()
      undf_aspc1_x = x_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for aspc1_lhs
      !
      ndf_aspc1_lhs = lhs_proxy%vspace%get_ndf()
      undf_aspc1_lhs = lhs_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(cell)
      !$omp do schedule(static)
      DO cell=loop0_start,loop0_stop
        !
        CALL dg_matrix_vector_code(cell, nlayers, lhs_proxy%data, x_proxy%data, matrix_proxy%ncell_3d, &
&matrix_proxy%local_stencil, ndf_adspc1_lhs, undf_adspc1_lhs, map_adspc1_lhs(:,cell), &
&ndf_aspc1_x, undf_aspc1_x, map_aspc1_x(:,cell))
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! Set halos dirty/clean for fields modified in the above loop(s)
      !
      CALL lhs_proxy%set_dirty()

    END SUBROUTINE invoke_dg_matrix_vector_rtran

!-------------------------------------------------------------------------------
!> Routine to perform higher-order reconstruction of a field on a fine mesh to
!! a coarse mesh. There is a field in this kernel that requires both a mesh
!! argument and a stencil argument, and PSyclone currently does not support
!! this. The issue to address this is #1983
  subroutine invoke_prolong_scalar_linear_kernel_type(target_field, source_field, stencil_extent)

    use prolong_scalar_linear_kernel_mod, only: prolong_scalar_linear_kernel_code
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

!==============================================================================

  subroutine invoke_r32_field_min_max(field_min_norm, &
                                  field_max_norm, &
                                  r32_field)

    use scalar_r32_mod,     only: scalar_r32_type
    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use mesh_mod,           only: mesh_type
    use field_r32_mod,      only: field_r32_type, field_r32_proxy_type

    implicit none

    real(kind=real32),              intent(out)  :: field_min_norm
    real(kind=real32),              intent(out)  :: field_max_norm
    type(field_r32_type),            intent(in)  :: r32_field
    type(scalar_r32_type)                        :: global_min, global_max
    integer(kind=i_def)                          :: df
    real(kind=real32), allocatable, dimension(:) :: l_field_min_norm
    real(kind=real32), allocatable, dimension(:) :: l_field_max_norm
    real(kind=real32)                            :: minv, maxv
    integer(kind=i_def)                          :: th_idx
    integer(kind=i_def)                          :: loop0_start, loop0_stop
    integer(kind=i_def)                          :: nthreads
    type(field_r32_proxy_type)                   :: field_proxy
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = r32_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    ALLOCATE (l_field_min_norm(nthreads))
    ALLOCATE (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    DO df=loop0_start,loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                 field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                 field_proxy%data(df))
    END DO
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    DO th_idx=2,nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    END DO
    DEALLOCATE (l_field_min_norm, l_field_max_norm)
    global_min%value = field_min_norm
    global_max%value = field_max_norm
    field_min_norm = global_min%get_min()
    field_max_norm = global_max%get_max()
    !
  end subroutine invoke_r32_field_min_max

  !-------------------------------------------------------------------------------
  subroutine invoke_r64_field_min_max(field_min_norm, &
                                  field_max_norm, &
                                  r64_field)

    use scalar_r64_mod,     only: scalar_r64_type
    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use mesh_mod,           only: mesh_type
    use field_r64_mod,      only: field_r64_type, field_r64_proxy_type

    implicit none

    real(kind=real64),               intent(out) :: field_min_norm
    real(kind=real64),               intent(out) :: field_max_norm
    type(field_r64_type),             intent(in) :: r64_field
    type(scalar_r64_type)                        :: global_min, global_max
    integer(kind=i_def)                          :: df
    real(kind=real64), allocatable, dimension(:) :: l_field_min_norm
    real(kind=real64), allocatable, dimension(:) :: l_field_max_norm
    real(kind=real64)                            :: minv, maxv
    integer(kind=i_def)                          :: th_idx
    integer(kind=i_def)                          :: loop0_start, loop0_stop
    integer(kind=i_def)                          :: nthreads
    type(field_r64_proxy_type)                   :: field_proxy
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = r64_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    ALLOCATE (l_field_min_norm(nthreads))
    ALLOCATE (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    DO df=loop0_start,loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                 field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                 field_proxy%data(df))
    END DO
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    DO th_idx=2,nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    END DO
    DEALLOCATE (l_field_min_norm, l_field_max_norm)
    global_min%value = field_min_norm
    global_max%value = field_max_norm
    field_min_norm = global_min%get_min()
    field_max_norm = global_max%get_max()
    !
  end subroutine invoke_r64_field_min_max
  !-------------------------------------------------------------------------------
!> Routine to perform injection of a multidata field on a fine mesh to
!> a coarse mesh. Intermesh kernels cannot currently take in integer arguments,
!> this has been raised as an issue https://github.com/stfc/PSyclone/issues/2504
  subroutine invoke_prolong_multidata_scalar_kernel_type(target_field, source_field, ndata)

    use prolong_multidata_scalar_kernel_mod, only: prolong_multidata_scalar_kernel_code
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
    use restrict_multidata_scalar_kernel_mod, only: restrict_multidata_scalar_kernel_code
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

    use prolong_multidata_linear_kernel_mod, only: prolong_multidata_linear_kernel_code
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
