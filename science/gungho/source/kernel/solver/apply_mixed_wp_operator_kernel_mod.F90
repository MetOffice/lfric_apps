!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the LHS of the semi-implicit system for the
!!         vertical velocity and pressure equations:
!!        (lhs_w) = norm_u*(Mu*u - P2t*t - grad*p),
!!        lhs_p = M3p*p - P3t*t + Q32*u,
!!        with t = -Mt^(-1) * Pt2*u
module apply_mixed_wp_operator_kernel_mod

use argument_mod,      only : arg_type,              &
                              GH_FIELD, GH_OPERATOR, &
                              GH_READ,               &
                              GH_WRITE,              &
                              GH_REAL, CELL_COLUMN
use constants_mod,     only : r_solver, i_def
use kernel_mod,        only : kernel_type
use fs_continuity_mod, only : W2, W3, Wtheta, W2h, W2v

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: apply_mixed_wp_operator_kernel_type
  private
  type(arg_type) :: meta_args(14) = (/                       &
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W2v),        & ! lhs_w
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),         & ! lhs_p
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2h),        & ! uv'
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2v),        & ! w'
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),         & ! exner'
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  Wtheta, W2), & ! Ptheta2
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),     & ! Mtheta^-1
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W2, W2),     & ! Mu^{c,d}
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W2, Wtheta), & ! P2theta
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W2, W3),     & ! grad
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2),         & ! norm_u
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W3),     & ! m3p
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W2),     & ! q32
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, Wtheta)  & ! p3t
       /)
  integer :: operates_on = CELL_COLUMN
  contains
  procedure, nopass :: apply_mixed_wp_operator_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: apply_mixed_wp_operator_code

contains

!> @brief Compute the LHS of the semi-implicit system
!> @param[in]     cell          Horizontal cell index
!> @param[in]     nlayers       Number of layers
!> @param[in,out] lhs_w         Mixed operator applied to the vertical momentum equation
!> @param[in,out] lhs_p         Mixed operator applied to the equation of state
!> @param[in]     wind_uv       Horizontal wind field
!> @param[in]     wind_w        Vertical wind field
!> @param[in]     exner         Exner pressure field
!> @param[in]     ncell0        Total number of cells for the pt2 operator
!> @param[in]     pt2           Projection operator from W2 to Wtheta
!> @param[in]     mt_lumped_inv Lumped inverse mass matrix for the Wtheta space
!> @param[in]     ncell1        Total number of cells for the mu_cd operator
!> @param[in]     mu_cd         Generalised mass matrix for the momentum equation
!> @param[in]     ncell2        Total number of cells for the p2t operator
!> @param[in]     p2t           Generalised projection matrix from Wtheta to W2
!> @param[in]     ncell3        Total number of cells for the grad operator
!> @param[in]     grad          Generalised gradient operator for the momentum equation
!> @param[in]     norm_u        Normalisation field for the momentum equation
!> @param[in]     ncell4        Total number of cells for the m3p operator
!> @param[in]     m3p           Weighted mass matrix for the W3 space
!> @param[in]     ncell5        Total number of cells for the q32 operator
!> @param[in]     q32           Projection operator from W2 to W3
!> @param[in]     ncell6        Total number of cells for the p3t operator
!> @param[in]     p3t           Projection operator from W2 to Wtheta
!> @param[in]     ndf_w2v       Number of degrees of freedom per cell for the vertical wind space
!> @param[in]     undf_w2v      Unique number of degrees of freedom for the vertical wind space
!> @param[in]     map_w2v       Dofmap for the cell at the base of the column for the vertical wind space
!> @param[in]     ndf_w3        Norm_umber of degrees of freedom per cell for the pressure space
!> @param[in]     ndf_w3        Unique number of degrees of freedom for the pressure space
!> @param[in]     map_w3        Dofmap for the cell at the base of the column for the pressure space
!> @param[in]     ndf_w2h       Number of degrees of freedom per cell for the horizontal wind space
!> @param[in]     undf_w2h      Unique number of degrees of freedom for the horizontal wind space
!> @param[in]     map_w2h       Dofmap for the cell at the base of the column for the horizontal wind space
!> @param[in]     undf_wt       Number of degrees of freedom per cell for the potential
!!                              temperature space
!> @param[in]     undf_wt       Unique number of degrees of freedom for the potential
!!                              temperature space
!> @param[in]     map_wt        Dofmap for the cell at the base of the column for the
!!                              potential temperature space
!> @param[in]     ndf_w2        Number of degrees of freedom per cell for the wind space
!> @param[in]     undf_w2       Unique number of degrees of freedom for the wind space
!> @param[in]     map_w2        Dofmap for the cell at the base of the column for the wind space
subroutine apply_mixed_wp_operator_code(cell,                       &
                                        nlayers,                    &
                                        lhs_w,                      &
                                        lhs_p,                      &
                                        wind_uv, wind_w, exner,     &
                                        ncell0, pt2,                &
                                        mt_lumped_inv,              &
                                        ncell1, mu_cd,              &
                                        ncell2, P2t,                &
                                        ncell3, grad,               &
                                        norm_u,                     &
                                        ncell4, m3p,                &
                                        ncell5, q32,                &
                                        ncell6, p3t,                &
                                        ndf_w2v, undf_w2v, map_w2v, &
                                        ndf_w3, undf_w3, map_w3,    &
                                        ndf_w2h, undf_w2h, map_w2h, &
                                        ndf_wt, undf_wt, map_wt,    &
                                        ndf_w2, undf_w2, map_w2)

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: cell, nlayers
  integer(kind=i_def),                     intent(in) :: ncell0, ncell1, ncell2, ncell3
  integer(kind=i_def),                     intent(in) :: ncell4, ncell5, ncell6
  integer(kind=i_def),                     intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                     intent(in) :: undf_w2h, ndf_w2h
  integer(kind=i_def),                     intent(in) :: undf_w2v, ndf_w2v
  integer(kind=i_def),                     intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def),                     intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
  integer(kind=i_def), dimension(ndf_w2v), intent(in) :: map_w2v
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

  ! Fields
  real(kind=r_solver), dimension(undf_w2v), intent(inout) :: lhs_w
  real(kind=r_solver), dimension(undf_w3),  intent(inout) :: lhs_p
  real(kind=r_solver), dimension(undf_w2h), intent(in)    :: wind_uv
  real(kind=r_solver), dimension(undf_w2v), intent(in)    :: wind_w
  real(kind=r_solver), dimension(undf_w2),  intent(in)    :: norm_u
  real(kind=r_solver), dimension(undf_wt),  intent(in)    :: mt_lumped_inv
  real(kind=r_solver), dimension(undf_w3),  intent(in)    :: exner

  ! Operators
  real(kind=r_solver), dimension(ncell0, ndf_wt, ndf_w2), intent(in) :: pt2
  real(kind=r_solver), dimension(ncell1, ndf_w2, ndf_w2), intent(in) :: mu_cd
  real(kind=r_solver), dimension(ncell2, ndf_w2, ndf_wt), intent(in) :: p2t
  real(kind=r_solver), dimension(ncell3, ndf_w2, ndf_w3), intent(in) :: grad
  real(kind=r_solver), dimension(ncell4, ndf_w3, ndf_w3), intent(in) :: m3p
  real(kind=r_solver), dimension(ncell5, ndf_w3, ndf_w2), intent(in) :: q32
  real(kind=r_solver), dimension(ncell6, ndf_w3, ndf_wt), intent(in) :: p3t

  ! Internal variables
  integer(kind=i_def)                                :: df, df2, ij,   &
                                                        nm1, iw3, iwt, &
                                                        iw2, iw2h, iw2v
  real(kind=r_solver), dimension(0:nlayers-1,ndf_w2) :: u_e
  real(kind=r_solver), dimension(0:nlayers)          :: t_col

  ! Set up some useful shorthands for indices
  ij = (cell-1)*nlayers + 1
  nm1 = nlayers-1
  iw3 = map_w3(1)
  iwt = map_wt(1)

  ! Create the element velocity field
  do df = 1, ndf_w2h
    iw2h = map_w2h(df)
    u_e(:,df) = wind_uv(iw2h:iw2h+nm1)
  end do
  do df = 1, ndf_w2v
    iw2v = map_w2v(df)
    u_e(:,ndf_w2h+df) = wind_w(iw2v:iw2v+nm1)
  end do

  ! Compute t for the column
  t_col(:) = 0.0_r_solver
  do df = 1, ndf_w2
    t_col(0:nm1)   = t_col(0:nm1)   - pt2(ij:ij+nm1, 1, df)*u_e(:,df)
    t_col(1:nm1+1) = t_col(1:nm1+1) - pt2(ij:ij+nm1, 2, df)*u_e(:,df)
  end do
  t_col(:) = t_col(:) * mt_lumped_inv(iwt:iwt+1+nm1)

  ! LHS W
  iw2v = map_w2v(1)
  lhs_w(iw2v:iw2v+nlayers) = 0.0_r_solver

  do df = 1, ndf_w2v
    iw2v = map_w2v(df)
    iw2  = map_w2(ndf_w2h+df)
    lhs_w(iw2v:iw2v+nm1) = lhs_w(iw2v:iw2v+nm1) &
                         + norm_u(iw2:iw2+nm1)*( &
                         - p2t(ij:ij+nm1, ndf_w2h+df, 1)*t_col(0:nm1)   &
                         - p2t(ij:ij+nm1, ndf_w2h+df, 2)*t_col(1:nm1+1) &
                         - grad(ij:ij+nm1, ndf_w2h+df, 1)*exner(iw3:iw3+nm1))

  end do
  do df2 = 1, ndf_w2
    do df = 1, ndf_w2v
      iw2v = map_w2v(df)
      iw2  = map_w2(ndf_w2h+df)
      lhs_w(iw2v:iw2v+nm1) = lhs_w(iw2v:iw2v+nm1) &
                           + norm_u(iw2:iw2+nm1)* &
                             mu_cd(ij:ij+nm1, ndf_w2h+df, df2)*u_e(:,df2)

    end do
  end do
  ! Set BC for lhs_w
  lhs_w(map_w2v(1)) = 0.0_r_solver
  lhs_w(map_w2v(2)+nlayers-1) = 0.0_r_solver

  ! LHS P
  lhs_p(iw3:iw3+nm1) = m3p(ij:ij+nm1, 1, 1)*exner(iw3:iw3+nm1) &
                     - p3t(ij:ij+nm1, 1, 1)*t_col(0:nm1)       &
                     - p3t(ij:ij+nm1, 1, 2)*t_col(1:nm1+1)
  do df = 1, ndf_w2
    lhs_p(iw3:iw3+nm1) = lhs_p(iw3:iw3+nm1) + q32(ij:ij+nm1, 1, df)*u_e(:,df)
  end do

end subroutine apply_mixed_wp_operator_code

end module apply_mixed_wp_operator_kernel_mod
