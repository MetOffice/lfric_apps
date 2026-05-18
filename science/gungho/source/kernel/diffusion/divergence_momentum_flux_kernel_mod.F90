!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculate divergence of momentum flux and update wind field
!>        as detailed in UMDP28.
!>
module divergence_momentum_flux_kernel_mod

  use argument_mod,                  only : arg_type,                          &
                                            GH_FIELD, GH_REAL, GH_INTEGER,     &
                                            GH_READ, GH_WRITE,                 &
                                            STENCIL, CROSS, CELL_COLUMN,       &
                                            ANY_SPACE_1,                       &
                                            ANY_DISCONTINUOUS_SPACE_3,         &
                                            ANY_DISCONTINUOUS_SPACE_9
  use constants_mod,                 only : r_def, i_def, l_def
  use fs_continuity_mod,             only : Wtheta, W3, W2, W2H, W1
  use kernel_mod,                    only : kernel_type
  use sci_face_selector_support_mod, only : face_from_face_selector
  use mixing_config_mod,             only : fullstress

  implicit none
  private

  type, public, extends(kernel_type) :: divergence_momentum_flux_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                       &
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE,  W2),                      &! u_inc
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   W3, STENCIL(CROSS)),      &! uflux_w3
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   W3, STENCIL(CROSS)),      &! vflux_w3
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   W1),                      &! uflux_w1
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   W1),                      &! vflux_w1
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   ANY_SPACE_1),             &! wflux_sh_w2h
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   W2),                      &! detj_at_w2
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   Wtheta),                  &! rho_in_wth
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   W2H),                     &! rho_in_w2h
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,   ANY_DISCONTINUOUS_SPACE_9,&
                                                             STENCIL(CROSS)),  &! panel_id
         arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),&! face_selector_ew
         arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &! face_selector_ns
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: divergence_momentum_flux_code
  end type

  public :: divergence_momentum_flux_code

contains

!> @brief Calculates dievrgence of momentum flux.
!> @param[in]     nlayers          Number of layers in the mesh
!> @param[in,out] u_n              Increment of wind field
!> @param[in]     uflux_w3         One component of flux of U
!> @param[in]     smap_w3_size     Size of the stencil map for uflux_w3
!> @param[in]     smap_w3          Stencil map for uflux_w3
!> @param[in]     vflux_w3         One component of flux of V
!> @param[in]     smap_w3v_size    Size of the stencil map for uflux_w3v
!> @param[in]     smap_w3v         Stencil map for uflux_w3v
!> @param[in]     uflux_w1         One component of flux of U
!> @param[in]     vflux_w1         One component of flux of V
!> @param[in]     wflux_sh_w2h     One component of flux of W
!> @param[in]     detj_at_w2       Cell volume at W2 points
!> @param[in]     rho_in_wth       Density field in Wtheta space
!> @param[in]     rho_in_w2h       Density field in W2H space
!> @param[in]     panel_id         The ID number of the current panel
!> @param[in]     smap_pid_size    Size of the stencil map for panel_id
!> @param[in]     smap_pid         Stencil map for panel_id
!> @param[in]     face_selector_ew 2D field indicating which W/E faces to loop
!> @param[in]     face_selector_ns 2D field indicating which N/S faces to loop
!> @param[in]     ndf_w2           Number of DOFs for W2 space
!> @param[in]     undf_w2          Number of unique DOFs for W2 space
!> @param[in]     map_w2           Dofmap for the cell at the base of the column
!> @param[in]     ndf_w3           Number of DOFs for W3 space
!> @param[in]     undf_w3          Number of unique DOFs for W3 space
!> @param[in]     map_w3           Dofmap for the cell at the base of the column
!> @param[in]     ndf_w1           Number of DOFs for W1 space
!> @param[in]     undf_w1          Number of unique DOFs for W1 space
!> @param[in]     map_w1           Dofmap for the cell at the base of the column
!> @param[in]     ndf_sh_w2h       Number of DOFs for shifted W2H space
!> @param[in]     undf_sh_w2h      Number of unique DOFs for shifted W2H space
!> @param[in]     map_sh_w2h       Dofmap for the cell at the base of the column
!> @param[in]     ndf_wt           Number of DOFs for Wtheta space
!> @param[in]     undf_wt          Number of unique DOFs for Wtheta space
!> @param[in]     map_wt           Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2h          Number of DOFs for W2H space
!> @param[in]     undf_w2h         Number of unique DOFs for W2H space
!> @param[in]     map_w2h          Dofmap for the cell at the base of the column
!> @param[in]     ndf_pid          Number of DOFs for pid space
!> @param[in]     undf_pid         Number of unique DOFs for pid space
!> @param[in]     map_pid          Dofmap for the cell at the base of the column
!> @param[in]     ndf_w3_2d        Number of DOFs for 2D W3 space
!> @param[in]     undf_w3_2d       umber of unique DOFs for 2D W3 space
!> @param[in]     map_w3_2d        Dofmap for the cell at the base of the column
subroutine divergence_momentum_flux_code( nlayers,                             &
                                          u_inc,                               &
                                          uflux_w3,                            &
                                          smap_w3_size, smap_w3,               &
                                          vflux_w3,                            &
                                          smap_w3v_size, smap_w3v,             &
                                          uflux_w1,                            &
                                          vflux_w1,                            &
                                          wflux_sh_w2h,                        &
                                          detj_at_w2,                          &
                                          rho_in_wth,                          &
                                          rho_in_w2h,                          &
                                          panel_id,                            &
                                          smap_pid_size, smap_pid,             &
                                          face_selector_ew,                    &
                                          face_selector_ns,                    &
                                          ndf_w2, undf_w2, map_w2,             &
                                          ndf_w3, undf_w3, map_w3,             &
                                          ndf_w1, undf_w1, map_w1,             &
                                          ndf_sh_w2h, undf_sh_w2h, map_sh_w2h, &
                                          ndf_wt, undf_wt, map_wt,             &
                                          ndf_w2h, undf_w2h, map_w2h,          &
                                          ndf_pid, undf_pid, map_pid,          &
                                          ndf_w3_2d, undf_w3_2d, map_w3_2d     &
                                         )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: ndf_w1, undf_w1
  integer(kind=i_def), intent(in) :: map_w1(ndf_w1)
  integer(kind=i_def), intent(in) :: ndf_sh_w2h, undf_sh_w2h
  integer(kind=i_def), intent(in) :: map_sh_w2h(ndf_sh_w2h)
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: map_wt(ndf_wt)
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
  integer(kind=i_def), intent(in) :: smap_w3_size
  integer(kind=i_def), intent(in) :: smap_w3(ndf_w3,smap_w3_size)
  integer(kind=i_def), intent(in) :: smap_w3v_size
  integer(kind=i_def), intent(in) :: smap_w3v(ndf_w3,smap_w3v_size)
  integer(kind=i_def), intent(in) :: smap_pid_size
  integer(kind=i_def), intent(in) :: smap_pid(ndf_pid,smap_pid_size)

  real(kind=r_def), dimension(undf_w2),     intent(inout) :: u_inc
  real(kind=r_def), dimension(undf_w3),     intent(in)    :: uflux_w3,         &
                                                             vflux_w3
  real(kind=r_def), dimension(undf_w1),     intent(in)    :: uflux_w1,         &
                                                             vflux_w1
  real(kind=r_def), dimension(undf_sh_w2h), intent(in)    :: wflux_sh_w2h
  real(kind=r_def), dimension(undf_w2),     intent(in)    :: detj_at_w2
  real(kind=r_def), dimension(undf_wt),     intent(in)    :: rho_in_wth
  real(kind=r_def), dimension(undf_w2h),    intent(in)    :: rho_in_w2h
  real(kind=r_def), dimension(undf_pid),    intent(in)    :: panel_id
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)  :: face_selector_ew, &
                                                             face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: k, df, j
  integer(kind=i_def) :: dfp3
  real(kind=r_def)    :: dflux
  real(kind=r_def), dimension(0:nlayers-1) :: r_volume

  integer(kind=i_def) :: stencil_cell
  integer(kind=i_def) :: cell_panel, stencil_panel, panel_edge
  integer(kind=i_def), dimension(4) :: vec_dir_x, vec_dir_y
  logical(kind=l_def), dimension(4) :: rotate

  integer(kind=i_def), parameter :: grad_sign(4) = (/-1, -1,  1, 1 /)

  ! If the full stencil isn't available, we must be at the domain edge.
  ! The increment is already 0, so we just exit the routine.
  if (smap_w3_size < 5_i_def) then
    return
  end if

  cell_panel = int(panel_id(map_pid(1)), i_def)

  do df = 1, 4
    stencil_cell = df + 1
    stencil_panel = int(panel_id(smap_pid(1, stencil_cell)), i_def)
    ! Create panel_edge to check whether a panel is changing
    panel_edge = 10*cell_panel + stencil_panel

    select case (panel_edge)
    case (41, 32, 16, 25, 64, 53)
      rotate(df) = .true.
      vec_dir_x(df) = 1_i_def
      vec_dir_y(df) = -1_i_def
    case (14, 23, 61, 52, 46, 35)
      rotate(df) = .true.
      vec_dir_x(df) = -1_i_def
      vec_dir_y(df) = 1_i_def
    case default
      rotate(df) = .false.
      vec_dir_x(df) = 1_i_def
      vec_dir_y(df) = 1_i_def
    end select
  end do

  do j = 1, ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))
    ! Calculate increment of U
    if (mod(df,2) == 1) then
      do k = 0, nlayers-1
        r_volume(k) = 1.0_r_def / detj_at_w2(map_w2(df)+k) / &
                                  rho_in_w2h(map_w2h(df)+k)
      end do

      ! Calculate increment of U by X direction divergence of flux
      if (rotate(df)) then
        do k = 0, nlayers - 1
          dflux = grad_sign(df) * &
                  ( vec_dir_x(df) * vflux_w3(smap_w3(1,df+1)+k) - &
                    uflux_w3(smap_w3(1,1)+k) )
          u_inc(map_w2(df)+k) = -dflux * r_volume(k)
        end do
      else
        do k = 0, nlayers - 1
          dflux = grad_sign(df) * &
                  ( uflux_w3(smap_w3(1,df+1)+k) - uflux_w3(smap_w3(1,1)+k) )
          u_inc(map_w2(df)+k) = -dflux * r_volume(k)
        end do
      end if

      ! Calculate increment of U by Y direction divergence of flux
      if (df == 1) then
        dfp3 = 8
      else
        dfp3 = df + 3
      end if
      do k = 0, nlayers - 1
        dflux = grad_sign(df) * &
                ( uflux_w1(map_w1(df+4)+k) - uflux_w1(map_w1(dfp3)+k) )
        u_inc(map_w2(df)+k) = u_inc(map_w2(df)+k) - dflux * r_volume(k)
      end do

      if (fullstress) then
        ! Calculate increment of U by Z direction divergence of flux
        do k = 0, nlayers - 1
          dflux = uflux_w1(map_w1(df)+k+1) - uflux_w1(map_w1(df)+k)
          u_inc(map_w2(df)+k) = u_inc(map_w2(df)+k) - dflux * r_volume(k)
        end do
      end if

    end if

    ! Calculate increment of V
    if (mod(df,2) == 0) then
      do k = 0, nlayers-1
        r_volume(k) = 1.0_r_def / detj_at_w2(map_w2(df)+k) / &
                                  rho_in_w2h(map_w2h(df)+k)
      end do

      ! Calculate increment of V by Y direction divergence of flux
      if (rotate(df)) then
        do k = 0, nlayers - 1
          dflux = grad_sign(df) * &
                  ( vec_dir_y(df) * uflux_w3(smap_w3(1,df+1)+k) - &
                    vflux_w3(smap_w3(1,1)+k) )
          ! Note u_inc at 2nd and 4th dof indicate from north to south wind.
          u_inc(map_w2(df)+k) = dflux * r_volume(k)
        end do
      else
        do k = 0, nlayers - 1
          dflux = grad_sign(df) * &
                  ( vflux_w3(smap_w3(1,df+1)+k) - vflux_w3(smap_w3(1,1)+k) )
          ! Note u_inc at 2nd and 4th dof indicate from north to south wind.
          u_inc(map_w2(df)+k) = dflux * r_volume(k)
        end do
      end if

      ! Calculate increment of V by X direction divergence of flux
      do k = 0, nlayers - 1
        dflux = grad_sign(df) * &
                ( vflux_w1(map_w1(df+3)+k) - vflux_w1(map_w1(df+4)+k) )
        ! Note u_inc at 2nd and 4th dof indicate from north to south wind.
        u_inc(map_w2(df)+k) = u_inc(map_w2(df)+k) + dflux * r_volume(k)
      end do

      if (fullstress) then
        ! Calculate increment of V by Z direction divergence of flux
        do k = 0, nlayers - 1
          dflux = vflux_w1(map_w1(df)+k+1) - vflux_w1(map_w1(df)+k)
          ! Note u_inc at 2nd and 4th dof indicate from north to south wind.
          u_inc(map_w2(df)+k) = u_inc(map_w2(df)+k) + dflux * r_volume(k)
        end do
      end if

    end if
  end do

  ! Calculate increment of W by X/Y direction divergence of flux
  do k = 1, nlayers - 1
    dflux = ( wflux_sh_w2h(map_sh_w2h(3)+k) - wflux_sh_w2h(map_sh_w2h(1)+k) + &
              wflux_sh_w2h(map_sh_w2h(4)+k) - wflux_sh_w2h(map_sh_w2h(2)+k) )
    u_inc(map_w2(5)+k) = - dflux / detj_at_w2(map_w2(5)+k) / &
                           rho_in_wth(map_wt(1)+k)
  end do

end subroutine divergence_momentum_flux_code

end module divergence_momentum_flux_kernel_mod
