!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculate diffusion flux of momentum as detailed in UMDP28.
!>
module diffusion_momentum_flux_kernel_mod

  use argument_mod,                  only : arg_type,                          &
                                            GH_FIELD, GH_REAL, GH_INTEGER,     &
                                            GH_READ, GH_WRITE,                 &
                                            STENCIL, CROSS, CELL_COLUMN,       &
                                            ANY_SPACE_1, ANY_SPACE_2,          &
                                            ANY_DISCONTINUOUS_SPACE_3,         &
                                            ANY_DISCONTINUOUS_SPACE_9
  use constants_mod,                 only : r_def, i_def
  use fs_continuity_mod,             only : W3, W2, W1, Wtheta, W2H
  use kernel_mod,                    only : kernel_type
  use sci_face_selector_support_mod, only : face_from_face_selector
  use mixing_config_mod,             only : fullstress

  implicit none
  private

  type, public, extends(kernel_type) :: diffusion_momentum_flux_kernel_type
    private
    type(arg_type) :: meta_args(19) = (/                                       &
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                       &! uflux_w3
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W3),                       &! vflux_w3
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W1),                       &! uflux_w1
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W1),                       &! vflux_w1
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_SPACE_1),              &! wflux_sh_w2h
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2, STENCIL(CROSS)),       &! u_physics
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2, STENCIL(CROSS)),       &! dx_at_w2
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2, STENCIL(CROSS)),       &! dA_at_w2
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_SPACE_2),              &! dx_at_sh_w2
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_SPACE_2),              &! dA_at_sh_w2
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H, STENCIL(CROSS)),      &! visc_m_w2h
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_SPACE_1),              &! visc_m_sh_w2h
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                       &! visc_m_w3
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),                       &! rho_in_w3
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_SPACE_1),              &! rho_in_sh_w2h
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2H, STENCIL(CROSS)),      &! rho_in_w2h
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9, &
                                                            STENCIL(CROSS)),   &! panel_id
         arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),&! face_selector_ew
         arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &! face_selector_ns
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: diffusion_momentum_flux_code
  end type

  public :: diffusion_momentum_flux_code

contains

!> @brief Calculates momentum diffusion flux.
!> @param[in]     nlayers          Number of layers in the mesh
!> @param[in,out] uflux_w3         One component of flux of U
!> @param[in,out] vflux_w3         One component of flux of V
!> @param[in,out] uflux_w1         One component of flux of U
!> @param[in,out] vflux_w1         One component of flux of V
!> @param[in,out] wflux_sh_w2h     One component of flux of W
!> @param[in]     u_physics        Input wind field
!> @param[in]     smap_w2_size     Size of the stencil map for u_physics
!> @param[in]     smap_w2          Stencil map for u_physics
!> @param[in]     dx_at_w2         Distance between cell centres at W2 points
!> @param[in]     smap_dx_size     Size of the stencil map for dx_at_w2
!> @param[in]     smap_dx          Stencil map for dx_at_w2
!> @param[in]     dA_at_w2         Area of faces at W2 points
!> @param[in]     smap_dA_size     Size of the stencil map for dA_at_w2
!> @param[in]     smap_dA          Stencil map for dA_at_w2
!> @param[in]     dx_at_sh_w2      Distance between cell centres at shifted W2 points
!> @param[in]     dA_at_sh_w2      Area of faces at shifted W2 points
!> @param[in]     visc_m_w2h       Viscosity field in W2H space
!> @param[in]     smap_w2h_size    Size of the stencil map for visc_m_w2h
!> @param[in]     smap_w2h         Stencil map for visc_m_w2h
!> @param[in]     visc_m_sh_w2h    Viscosity field in shifted W2H space
!> @param[in]     visc_m_w3        Viscosity field in W3 space
!> @param[in]     rho_in_w3        Density field in W3 space
!> @param[in]     rho_in_sh_w2h    Density field in shifted W2H space
!> @param[in]     rho_in_w2h       Density field in W2H space
!> @param[in]     smap_rho_size    Size of the stencil map for rho_in_w2h
!> @param[in]     smap_rho         Stencil map for rho_in_w2h
!> @param[in]     panel_id         The ID number of the current panel
!> @param[in]     smap_pid_size    Size of the stencil map for panel_id
!> @param[in]     smap_pid         Stencil map for panel_id
!> @param[in]     face_selector_ew 2D field indicating which W/E faces to loop
!> @param[in]     face_selector_ns 2D field indicating which N/S faces to loop
!> @param[in]     ndf_w3           Number of DOFs for W3 space
!> @param[in]     undf_w3          Number of unique DOFs for W3 space
!> @param[in]     map_w3           Dofmap for the cell at the base of the column
!> @param[in]     ndf_w1           Number of DOFs for W1 space
!> @param[in]     undf_w1          Number of unique DOFs for W1 space
!> @param[in]     map_w1           Dofmap for the cell at the base of the column
!> @param[in]     ndf_sh_w2h       Number of DOFs for shifted W2H space
!> @param[in]     undf_sh_w2h      Number of unique DOFs for shifted W2H space
!> @param[in]     map_sh_w2h       Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2           Number of DOFs for W2 space
!> @param[in]     undf_w2          Number of unique DOFs for W2 space
!> @param[in]     map_w2           Dofmap for the cell at the base of the column
!> @param[in]     ndf_sh_w2        Number of DOFs for shifted W2 space
!> @param[in]     undf_sh_w2       Number of unique DOFs for shifted W2 space
!> @param[in]     map_sh_w2        Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2h          Number of DOFs for W2H space
!> @param[in]     undf_w2h         Number of unique DOFs for W2H space
!> @param[in]     map_w2h          Dofmap for the cell at the base of the column
!> @param[in]     ndf_pid          Number of DOFs for pid space
!> @param[in]     undf_pid         Number of unique DOFs for pid space
!> @param[in]     map_pid          Dofmap for the cell at the base of the column
!> @param[in]     ndf_w3_2d        Number of DOFs for 2D W3 space
!> @param[in]     undf_w3_2d       umber of unique DOFs for 2D W3 space
!> @param[in]     map_w3_2d        Dofmap for the cell at the base of the column
subroutine diffusion_momentum_flux_code( nlayers,                              &
                                         uflux_w3,                             &
                                         vflux_w3,                             &
                                         uflux_w1,                             &
                                         vflux_w1,                             &
                                         wflux_sh_w2h,                         &
                                         u_physics,                            &
                                         smap_w2_size, smap_w2,                &
                                         dx_at_w2,                             &
                                         smap_dx_size, smap_dx,                &
                                         dA_at_w2,                             &
                                         smap_dA_size, smap_dA,                &
                                         dx_at_sh_w2,                          &
                                         dA_at_sh_w2,                          &
                                         visc_m_w2h,                           &
                                         smap_w2h_size, smap_w2h,              &
                                         visc_m_sh_w2h,                        &
                                         visc_m_w3,                            &
                                         rho_in_w3,                            &
                                         rho_in_sh_w2h,                        &
                                         rho_in_w2h,                           &
                                         smap_rho_size, smap_rho,              &
                                         panel_id,                             &
                                         smap_pid_size, smap_pid,              &
                                         face_selector_ew,                     &
                                         face_selector_ns,                     &
                                         ndf_w3, undf_w3, map_w3,              &
                                         ndf_w1, undf_w1, map_w1,              &
                                         ndf_sh_w2h, undf_sh_w2h, map_sh_w2h,  &
                                         ndf_w2, undf_w2, map_w2,              &
                                         ndf_sh_w2, undf_sh_w2, map_sh_w2,     &
                                         ndf_w2h, undf_w2h, map_w2h,           &
                                         ndf_pid, undf_pid, map_pid,           &
                                         ndf_w3_2d, undf_w3_2d, map_w3_2d      &
                                        )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: ndf_w1, undf_w1
  integer(kind=i_def), intent(in) :: map_w1(ndf_w1)
  integer(kind=i_def), intent(in) :: ndf_sh_w2h, undf_sh_w2h
  integer(kind=i_def), intent(in) :: map_sh_w2h(ndf_sh_w2h)
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: ndf_sh_w2, undf_sh_w2
  integer(kind=i_def), intent(in) :: map_sh_w2(ndf_sh_w2)
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
  integer(kind=i_def), intent(in) :: smap_w2_size
  integer(kind=i_def), intent(in) :: smap_w2(ndf_w2,smap_w2_size)
  integer(kind=i_def), intent(in) :: smap_dx_size
  integer(kind=i_def), intent(in) :: smap_dx(ndf_w2,smap_dx_size)
  integer(kind=i_def), intent(in) :: smap_dA_size
  integer(kind=i_def), intent(in) :: smap_dA(ndf_w2,smap_dA_size)
  integer(kind=i_def), intent(in) :: smap_w2h_size
  integer(kind=i_def), intent(in) :: smap_w2h(ndf_w2h,smap_w2h_size)
  integer(kind=i_def), intent(in) :: smap_rho_size
  integer(kind=i_def), intent(in) :: smap_rho(ndf_w2h,smap_rho_size)
  integer(kind=i_def), intent(in) :: smap_pid_size
  integer(kind=i_def), intent(in) :: smap_pid(ndf_pid,smap_pid_size)

  real(kind=r_def), dimension(undf_w3),     intent(inout) :: uflux_w3,         &
                                                             vflux_w3
  real(kind=r_def), dimension(undf_w1),     intent(inout) :: uflux_w1,         &
                                                             vflux_w1
  real(kind=r_def), dimension(undf_sh_w2h), intent(inout) :: wflux_sh_w2h
  real(kind=r_def), dimension(undf_w2),     intent(in)    :: u_physics,        &
                                                             dx_at_w2,         &
                                                             dA_at_w2
  real(kind=r_def), dimension(undf_sh_w2),  intent(in)    :: dx_at_sh_w2,      &
                                                             dA_at_sh_w2
  real(kind=r_def), dimension(undf_sh_w2h), intent(in)    :: visc_m_sh_w2h,    &
                                                             rho_in_sh_w2h
  real(kind=r_def), dimension(undf_w2h),    intent(in)    :: visc_m_w2h,       &
                                                             rho_in_w2h
  real(kind=r_def), dimension(undf_w3),     intent(in)    :: visc_m_w3,        &
                                                             rho_in_w3
  real(kind=r_def), dimension(undf_pid),    intent(in)    :: panel_id
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)  :: face_selector_ew, &
                                                             face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: k, df, j
  integer(kind=i_def) :: st_df_xu, st_df_yv
  integer(kind=i_def) :: df_xu, df_yv
  integer(kind=i_def) :: grad_sign_x, grad_sign_y
  real(kind=r_def)    :: dx_at_w1, dy_at_w1
  real(kind=r_def)    :: dudx_w3, dvdy_w3
  real(kind=r_def)    :: dudy_w1, dvdx_w1
  real(kind=r_def)    :: dudz_sh_w2, dwdx_sh_w2, dwdy_sh_w2
  real(kind=r_def)    :: dz_at_sh_w2
  real(kind=r_def)    :: dAx_at_w3, dAy_at_w3
  real(kind=r_def)    :: dAy_at_w1, dAx_at_w1
  real(kind=r_def)    :: dAz_at_sh_w2
  real(kind=r_def)    :: visc_m_w1
  real(kind=r_def)    :: rho_in_w1

  integer(kind=i_def) :: smap_w2_true(ndf_w2,smap_w2_size)
  integer(kind=i_def) :: smap_w2h_true(ndf_w2h,smap_w2h_size)
  integer(kind=i_def) :: stencil_cell
  integer(kind=i_def) :: cell_panel, stencil_panel, panel_edge
  integer(kind=i_def) :: vec_dir(ndf_w2,smap_w2_size)

  integer(kind=i_def), parameter :: grad_sign(4) = (/-1, -1,  1, 1 /)
  integer(kind=i_def), parameter :: wind_sign(4) = (/1, -1,  1, -1 /)

  ! Assumed direction for derivatives in this kernel is:
  !  y
  !  ^
  !  |_> x
  !

  ! Layout of dofs for the stencil map
  ! dimensions of map are (ndf, ncell)
  ! Horizontally:
  !
  !   -- 4 --
  !   |     |
  !   1     3
  !   |     |
  !   -- 2 --
  !
  ! df = 5 is in the centre on the bottom face
  ! df = 6 is in the centre on the top face

  ! The layout of the cells in the stencil is:
  !
  !          -----
  !          |   |
  !          | 5 |
  !     ---------------
  !     |    |   |    |
  !     |  2 | 1 |  4 |
  !     ---------------
  !          |   |
  !          | 3 |
  !          -----

  ! If the full stencil isn't available, we must be at the domain edge.
  ! The increment is already 0, so we just exit the routine.
  if (smap_w2_size < 5_i_def) then
    return
  end if

  ! The W2H DoF values change in orientation when we cross over a panel
  ! Vector directions parallel to the boundary (i.e. the winds on faces
  ! perpendicular to the boundary) also flip sign
  ! We need to take this into account by adjusting the stencil map used for
  ! the wind field. Do this by looking at whether the panel changes
  ! for other cells in the stencil
  cell_panel = int(panel_id(map_pid(1)), i_def)

  do stencil_cell = 1, smap_w2_size
    stencil_panel = int(panel_id(smap_pid(1, stencil_cell)), i_def)
    ! Create panel_edge to check whether a panel is changing
    panel_edge = 10*cell_panel + stencil_panel

    select case (panel_edge)
    case (41, 32, 16, 25, 64, 53)
      ! Clockwise rotation of panel
      smap_w2_true(1, stencil_cell) = smap_w2(2, stencil_cell)
      smap_w2_true(2, stencil_cell) = smap_w2(3, stencil_cell)
      smap_w2_true(3, stencil_cell) = smap_w2(4, stencil_cell)
      smap_w2_true(4, stencil_cell) = smap_w2(1, stencil_cell)
      ! Vertical dofs unchanged
      smap_w2_true(5, stencil_cell) = smap_w2(5, stencil_cell)
      smap_w2_true(6, stencil_cell) = smap_w2(6, stencil_cell)
      ! Clockwise rotation of panel
      smap_w2h_true(1, stencil_cell) = smap_w2h(2, stencil_cell)
      smap_w2h_true(2, stencil_cell) = smap_w2h(3, stencil_cell)
      smap_w2h_true(3, stencil_cell) = smap_w2h(4, stencil_cell)
      smap_w2h_true(4, stencil_cell) = smap_w2h(1, stencil_cell)
      ! Flip direction of vectors if necessary
      vec_dir(1,stencil_cell) = -1_i_def
      vec_dir(2,stencil_cell) = 1_i_def
      vec_dir(3,stencil_cell) = -1_i_def
      vec_dir(4,stencil_cell) = 1_i_def
    case (14, 23, 61, 52, 46, 35)
      ! Anti-clockwise rotation of panel
      smap_w2_true(1, stencil_cell) = smap_w2(4, stencil_cell)
      smap_w2_true(2, stencil_cell) = smap_w2(1, stencil_cell)
      smap_w2_true(3, stencil_cell) = smap_w2(2, stencil_cell)
      smap_w2_true(4, stencil_cell) = smap_w2(3, stencil_cell)
      ! Vertical dofs unchanged
      smap_w2_true(5, stencil_cell) = smap_w2(5, stencil_cell)
      smap_w2_true(6, stencil_cell) = smap_w2(6, stencil_cell)
      ! Anti-clockwise rotation of panel
      smap_w2h_true(1, stencil_cell) = smap_w2h(4, stencil_cell)
      smap_w2h_true(2, stencil_cell) = smap_w2h(1, stencil_cell)
      smap_w2h_true(3, stencil_cell) = smap_w2h(2, stencil_cell)
      smap_w2h_true(4, stencil_cell) = smap_w2h(3, stencil_cell)
      ! Flip direction of vectors if necessary
      vec_dir(1,stencil_cell) = 1_i_def
      vec_dir(2,stencil_cell) = -1_i_def
      vec_dir(3,stencil_cell) = 1_i_def
      vec_dir(4,stencil_cell) = -1_i_def
    case default
      ! Same panel or crossing panel with no rotation, so stencil map is unchanged
      smap_w2_true(:, stencil_cell) = smap_w2(:, stencil_cell)
      smap_w2h_true(:, stencil_cell) = smap_w2h(:, stencil_cell)
      vec_dir(:,stencil_cell) = 1_i_def
    end select
  end do

  ! --------------------------------------------------------------------
  ! Compute horizontal flux on W3 point.
  ! --------------------------------------------------------------------
  do k = 0, nlayers - 1
    dudx_w3 = ( u_physics(map_w2(3)+k) - u_physics(map_w2(1)+k) ) &
              * 2.0_r_def / (dx_at_w2(map_w2(3)+k) + dx_at_w2(map_w2(1)+k))
    dAx_at_w3 = 0.5_r_def * (dA_at_w2(map_w2(3)+k) + dA_at_w2(map_w2(1)+k))
    ! X direction flux of U
    uflux_w3(map_w3(1)+k) = -rho_in_w3(map_w3(1)+k) * &
                             visc_m_w3(map_w3(1)+k) * dudx_w3 * dAx_at_w3

    dvdy_w3 = ( u_physics(map_w2(4)+k) - u_physics(map_w2(2)+k) ) &
              * 2.0_r_def / (dx_at_w2(map_w2(4)+k) + dx_at_w2(map_w2(2)+k))
    ! Note u_physics at 2nd and 4th dof indicate from north to south wind.
    dvdy_w3 = -1 * dvdy_w3
    dAy_at_w3 = 0.5_r_def * (dA_at_w2(map_w2(4)+k) + dA_at_w2(map_w2(2)+k))
    ! Y direction flux of V
    vflux_w3(map_w3(1)+k) = -rho_in_w3(map_w3(1)+k) * &
                             visc_m_w3(map_w3(1)+k) * dvdy_w3 * dAy_at_w3
  end do

  if (fullstress) then
    do k = 0, nlayers - 1
      uflux_w3(map_w3(1)+k) = 2.0_r_def * uflux_w3(map_w3(1)+k)
      vflux_w3(map_w3(1)+k) = 2.0_r_def * vflux_w3(map_w3(1)+k)
    end do
  end if

  ! --------------------------------------------------------------------
  ! Compute horizontal flux on W1 point.
  ! --------------------------------------------------------------------
  !
  !              u(1,5)              u(3,5)
  !                 |                   |
  !  v(4,2) --- flux(8)--- v(4) --- flux(7)--- v(4,4)
  !                 |                   |
  !               u(1)                u(3)
  !                 |                   |
  !  v(2,2) --- flux(5)--- v(2) --- flux(6)--- v(2,4)
  !                 |                   |
  !              u(1,3)              u(3,3)
  !
  do df = 1, 4
    ! Only calculate this dof if it hasn't already been done
    if (uflux_w1(map_w1(df+4)) == 0.0_r_def) then
      select case (df)
      case (1)
        grad_sign_x = -1
        grad_sign_y = -1
        st_df_xu = 3
        st_df_yv = 2
        df_xu = 1
        df_yv = 2
      case (2)
        grad_sign_x = 1
        grad_sign_y = -1
        st_df_xu = 3
        st_df_yv = 4
        df_xu = 3
        df_yv = 2
      case (3)
        grad_sign_x = 1
        grad_sign_y = 1
        st_df_xu = 5
        st_df_yv = 4
        df_xu = 3
        df_yv = 4
      case (4)
        grad_sign_x = -1
        grad_sign_y = 1
        st_df_xu = 5
        st_df_yv = 2
        df_xu = 1
        df_yv = 4
      end select

      do k = 0, nlayers - 1
        visc_m_w1 = 0.25_r_def * &
                   ( visc_m_w2h(smap_w2h_true(df_xu, 1       )+k) + &
                     visc_m_w2h(smap_w2h_true(df_yv, 1       )+k) + &
                     visc_m_w2h(smap_w2h_true(df_xu, st_df_xu)+k) + &
                     visc_m_w2h(smap_w2h_true(df_yv, st_df_yv)+k) )
        rho_in_w1 = 0.25_r_def * &
                   ( rho_in_w2h(smap_w2h_true(df_xu, 1       )+k) + &
                     rho_in_w2h(smap_w2h_true(df_yv, 1       )+k) + &
                     rho_in_w2h(smap_w2h_true(df_xu, st_df_xu)+k) + &
                     rho_in_w2h(smap_w2h_true(df_yv, st_df_yv)+k) )

        dy_at_w1 = 0.5_r_def * &
          (dx_at_w2(smap_w2_true(df_yv,1)+k) + &
           dx_at_w2(smap_w2_true(df_yv,st_df_yv)+k))

        dudy_w1 = grad_sign_y * &
                ( vec_dir(df_xu,st_df_xu) * &
                  u_physics(smap_w2_true(df_xu,st_df_xu)+k) - &
                  u_physics(smap_w2_true(df_xu,1)+k) ) / dy_at_w1
        ! Y direction flux of U
        uflux_w1(map_w1(df+4)+k) = -rho_in_w1 * &
                                    visc_m_w1 * dudy_w1

        dx_at_w1 = 0.5_r_def * &
          (dx_at_w2(smap_w2_true(df_xu,1)+k) + &
           dx_at_w2(smap_w2_true(df_xu,st_df_xu)+k))
        dvdx_w1 = grad_sign_x * &
                ( vec_dir(df_yv,st_df_yv) * &
                  u_physics(smap_w2_true(df_yv,st_df_yv)+k) - &
                  u_physics(smap_w2_true(df_yv,1)+k) ) / dx_at_w1
        ! Note u_physics at 2nd and 4th dof indicate from north to south wind.
        dvdx_w1 = -1 * dvdx_w1
        ! X direction flux of V
        vflux_w1(map_w1(df+4)+k) = -rho_in_w1 * &
                                    visc_m_w1 * dvdx_w1
      end do

      if (fullstress) then
        do k = 0, nlayers - 1
          uflux_w1(map_w1(df+4)+k) = uflux_w1(map_w1(df+4)+k) + &
                                     vflux_w1(map_w1(df+4)+k)
          vflux_w1(map_w1(df+4)+k) = uflux_w1(map_w1(df+4)+k)
        end do
      end if

      do k = 0, nlayers - 1
        dAy_at_w1 = 0.5_r_def * &
          (dA_at_w2(smap_w2_true(df_yv,1)+k) + &
           dA_at_w2(smap_w2_true(df_yv,st_df_yv)+k))
        ! Y direction flux of U * Area
        uflux_w1(map_w1(df+4)+k) = uflux_w1(map_w1(df+4)+k) * dAy_at_w1

        dAx_at_w1 = 0.5_r_def * &
          (dA_at_w2(smap_w2_true(df_xu,1)+k) + &
           dA_at_w2(smap_w2_true(df_xu,st_df_xu)+k))
        ! X direction flux of V * Area
        vflux_w1(map_w1(df+4)+k) = vflux_w1(map_w1(df+4)+k) * dAx_at_w1
      end do

    end if
  end do

  ! --------------------------------------------------------------------
  ! Compute horizontal flux on sfifted W2 point.
  ! --------------------------------------------------------------------
  do j = 1, ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))
    do k = 1, nlayers-1
      dwdx_sh_w2 = grad_sign(df) * &
                   ( u_physics(smap_w2_true(5,df+1)+k) - &
                     u_physics(smap_w2_true(5,1)+k) ) / &
                     dx_at_sh_w2(map_sh_w2(df)+k)
      ! X/Y direction flux of W * Area
      wflux_sh_w2h(map_sh_w2h(df)+k) = -rho_in_sh_w2h(map_sh_w2h(df)+k) * &
                                        visc_m_sh_w2h(map_sh_w2h(df)+k) * &
                                        dwdx_sh_w2 * &
                                        dA_at_sh_w2(map_sh_w2(df)+k)
    end do

    if (fullstress) then
      do k = 1, nlayers-1
        dz_at_sh_w2 = 0.5_r_def * &
          ( dx_at_w2(smap_w2_true(5,df+1)+k) + dx_at_w2(smap_w2_true(5,1)+k) )
        dudz_sh_w2 = ( u_physics(map_w2(df)+k) - &
                       u_physics(map_w2(df)+k-1) ) / dz_at_sh_w2
        ! Note u_physics at 2nd and 4th dof indicate from north to south wind.
        dudz_sh_w2 = wind_sign(df) * dudz_sh_w2
        ! X/Y direction flux of W * Area
        wflux_sh_w2h(map_sh_w2h(df)+k) = wflux_sh_w2h(map_sh_w2h(df)+k) - &
                                         rho_in_sh_w2h(map_sh_w2h(df)+k) * &
                                         visc_m_sh_w2h(map_sh_w2h(df)+k) * &
                                         dudz_sh_w2 * &
                                         dA_at_sh_w2(map_sh_w2(df)+k)
      end do
    end if

  end do

  ! --------------------------------------------------------------------
  ! Compute vertical flux on shifted W2 point.
  ! --------------------------------------------------------------------
  if (fullstress) then
    do df = 1,3,2
      ! Only calculate this dof if it hasn't already been done
      if (uflux_w1(map_w1(df)+1) == 0.0_r_def) then
        do k = 1, nlayers - 1
          dAz_at_sh_w2 = 0.5_r_def * &
            ( dA_at_w2(smap_w2_true(5,1)+k) + &
              dA_at_w2(smap_w2_true(5,df+1)+k) )
          dwdx_sh_w2 = grad_sign(df) * ( u_physics(smap_w2_true(5,df+1)+k) - &
                                         u_physics(smap_w2_true(5,1)+k) ) / &
                                         dx_at_sh_w2(map_sh_w2(df)+k)
          ! Z direction flux of U * Area
          uflux_w1(map_w1(df)+k) = -rho_in_sh_w2h(map_sh_w2h(df)+k) * &
                                    visc_m_sh_w2h(map_sh_w2h(df)+k) * &
                                    dwdx_sh_w2 * dAz_at_sh_w2
        end do
      end if
    end do

    do df = 2,4,2
      ! Only calculate this dof if it hasn't already been done
      if (vflux_w1(map_w1(df)+1) == 0.0_r_def) then
        do k = 1, nlayers - 1
          dAz_at_sh_w2 = 0.5_r_def * &
            ( dA_at_w2(smap_w2_true(5,1)+k) + &
              dA_at_w2(smap_w2_true(5,df+1)+k) )
          dwdy_sh_w2 = grad_sign(df) * ( u_physics(smap_w2_true(5,df+1)+k) - &
                                         u_physics(smap_w2_true(5,1)+k) ) / &
                                         dx_at_sh_w2(map_sh_w2(df)+k)
          ! Z direction flux of V * Area
          vflux_w1(map_w1(df)+k) = -rho_in_sh_w2h(map_sh_w2h(df)+k) * &
                                    visc_m_sh_w2h(map_sh_w2h(df)+k) * &
                                    dwdy_sh_w2 * dAz_at_sh_w2
        end do
      end if
    end do
  end if

end subroutine diffusion_momentum_flux_code

end module diffusion_momentum_flux_kernel_mod
