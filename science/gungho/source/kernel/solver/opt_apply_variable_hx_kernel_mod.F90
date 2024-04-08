!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @todO Create unit test for this kernel, see #2935
module opt_apply_variable_hx_kernel_mod

  use argument_mod,      only : arg_type,             &
                                GH_FIELD, GH_SCALAR,  &
                                GH_OPERATOR, GH_REAL, &
                                GH_READ, GH_WRITE,    &
                                CELL_COLUMN
  use constants_mod,     only : r_double, r_single, i_def
  use fs_continuity_mod, only : W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: opt_apply_variable_hx_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                           &
         arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),             &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2),             &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),         &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),             &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3,     W2),     &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3,     Wtheta), &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  Wtheta, W2),     &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3,     W3),     &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                  &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3)              &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: opt_apply_variable_hx_code

  interface opt_apply_variable_hx_code
    module procedure &
      opt_apply_variable_hx_code_r_single, &
      opt_apply_variable_hx_code_r_double
  end interface


contains

!> @brief Applies the component of the helmholtz operator that maps from velocity space
!>        to the pressure space as well as the constant in space part, optimised for lowest
!>        order elements with horizontally discontinuous temperature space
!> @details The Helmholtz operator can be summarised as:
!>          \f[
!>             H(p) = Mp + \nabla.\left( \nabla p \right) + \bar{ \nabla p }
!>         \f]
!>        For a given p & \f[ \nabla p \f] this kernel applies the
!>        divergence \f[ \nabla. X \f] and averaging  \f[ \bar{X} \f]
!>        operators as well as the application of the mass matrix M
!> @param[in] cell Horizontal cell index
!> @param[in] nlayers Number of layers
!> @param[in,out] lhs Pressure field with helmholtz operator applied to it
!> @param[in] x Gradient of the pressure field in the velocity space
!> @param[in] mt_inv Lumped inverse mass matrix for the temperature space
!> @param[in] pressure Field that helmholtz operator is being applied to
!> @param[in] ncell_3d_1 Total number of cells for divergence matrix
!> @param[in] div Generalised divergence matrix
!> @param[in] ncell_3d_2 Total number of cells for p3t matrix
!> @param[in] p3t Mapping from temperature space to pressure space
!> @param[in] ncell_3d_3 Total number of cells for pt2 matrix
!> @param[in] pt2 Mapping from velocity space to temperature space
!> @param[in] ncell_3d_4 Total number of cells for m3 matrix
!> @param[in] m3 Mass matrix for the pressure space
!> @param[in] sgn +/- 1 Weight
!> @param[in] rhs_p Pressure field in w3
!> @param[in] ndf_w3 Number of degrees of freedom per cell for the pressure space
!> @param[in] undf_w3 Unique number of degrees of freedom  for the pressure space
!> @param[in] map_w3 Dofmap for the cell at the base of the column for the pressure space
!> @param[in] ndf_w2 Number of degrees of freedom per cell for the velocity space
!> @param[in] undf_w2 Unique number of degrees of freedom  for the velocity space
!> @param[in] map_w2 Dofmap for the cell at the base of the column for the velocity space
!> @param[in] ndf_wt Number of degrees of freedom per cell for the temperature space
!> @param[in] undf_wt Unique number of degrees of freedom  for the temperature space
!> @param[in] map_wt Dofmap for the cell at the base of the column for the temperature space

! R_DOUBLE PRECISION
! ==================
subroutine opt_apply_variable_hx_code_r_double(cell,                    &
                                               nlayers,                 &
                                               lhs, x,                  &
                                               mt_inv,                  &
                                               pressure,                &
                                               ncell_3d_1,              &
                                               div,                     &
                                               ncell_3d_2,              &
                                               p3t,                     &
                                               ncell_3d_3,              &
                                               pt2,                     &
                                               ncell_3d_4,              &
                                               m3,                      &
                                               sgn,                     &
                                               rhs_p,                   &
                                               ndf_w3, undf_w3, map_w3, &
                                               ndf_w2, undf_w2, map_w2, &
                                               ndf_wt, undf_wt, map_wt)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell_3d_1, ncell_3d_2
  integer(kind=i_def),                    intent(in) :: ncell_3d_3, ncell_3d_4
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_double), dimension(undf_w2), intent(in)    :: x
  real(kind=r_double), dimension(undf_wt), intent(in)    :: mt_inv
  real(kind=r_double), dimension(undf_w3), intent(inout) :: lhs
  real(kind=r_double), dimension(undf_w3), intent(in)    :: pressure
  real(kind=r_double),                     intent(in)    :: sgn
  real(kind=r_double), dimension(undf_w3), intent(in)    :: rhs_p

  real(kind=r_double), dimension(1,6,ncell_3d_1), intent(in) :: div
  real(kind=r_double), dimension(2,6,ncell_3d_2), intent(in) :: pt2
  real(kind=r_double), dimension(1,2,ncell_3d_3), intent(in) :: p3t
  real(kind=r_double), dimension(1,1,ncell_3d_4), intent(in) :: m3

  ! Internal variables
  integer(kind=i_def)            :: k, ik
  real(kind=r_double), dimension(2) :: t_e
  real(kind=r_double)               :: div_u

  ! Compute D * u + P3t * Mt^-1 * ( Pt2 * u ) + rhs_p
  ! Hard wired optimisation for desired configuration (p=0 elements with pt2
  ! only acting on vertical components of u )
  k = 0
  ik = (cell-1)*nlayers + k + 1

  t_e(1) = mt_inv(map_wt(1)+k)*p3t(1,1,ik)                              &
          *(pt2(1,5,ik)*x(map_w2(5)+k) + pt2(1,6,ik)  *x(map_w2(6)+k))
  t_e(2) = mt_inv(map_wt(2)+k)*p3t(1,2,ik)                              &
          *(pt2(2,5,ik)*x(map_w2(5)+k) + pt2(1,5,ik+1)*x(map_w2(5)+k+1) &
          + pt2(2,6,ik)*x(map_w2(6)+k) + pt2(1,6,ik+1)*x(map_w2(6)+k+1))

  div_u = div(1,1,ik)*x(map_w2(1)+k) + div(1,2,ik)*x(map_w2(2)+k) &
        + div(1,3,ik)*x(map_w2(3)+k) + div(1,4,ik)*x(map_w2(4)+k) &
        + div(1,5,ik)*x(map_w2(5)+k) + div(1,6,ik)*x(map_w2(6)+k)
  lhs(map_w3(1)+k) = m3(1,1,ik)*pressure(map_w3(1)+k) &
                   + sgn*(div_u + (t_e(1) + t_e(2))) + rhs_p(map_w3(1)+k)

  do k = 1,nlayers-2
    ik = (cell-1)*nlayers + k + 1

    t_e(1) = mt_inv(map_wt(1)+k)*p3t(1,1,ik)                              &
            *(pt2(1,5,ik)*x(map_w2(5)+k) + pt2(2,5,ik-1)*x(map_w2(5)+k-1) &
            + pt2(1,6,ik)*x(map_w2(6)+k) + pt2(2,6,ik-1)*x(map_w2(6)+k-1))
    t_e(2) = mt_inv(map_wt(2)+k)*p3t(1,2,ik)                              &
            *(pt2(2,5,ik)*x(map_w2(5)+k) + pt2(1,5,ik+1)*x(map_w2(5)+k+1) &
            + pt2(2,6,ik)*x(map_w2(6)+k) + pt2(1,6,ik+1)*x(map_w2(6)+k+1))

    div_u = div(1,1,ik)*x(map_w2(1)+k) + div(1,2,ik)*x(map_w2(2)+k) &
          + div(1,3,ik)*x(map_w2(3)+k) + div(1,4,ik)*x(map_w2(4)+k) &
          + div(1,5,ik)*x(map_w2(5)+k) + div(1,6,ik)*x(map_w2(6)+k)
    lhs(map_w3(1)+k) = m3(1,1,ik)*pressure(map_w3(1)+k) &
                     + sgn*(div_u + (t_e(1) + t_e(2))) + rhs_p(map_w3(1)+k)
  end do

  k = nlayers-1
  ik = (cell-1)*nlayers + k + 1

  t_e(1) = mt_inv(map_wt(1)+k)*p3t(1,1,ik)                              &
          *(pt2(1,5,ik)*x(map_w2(5)+k) + pt2(2,5,ik-1)*x(map_w2(5)+k-1) &
          + pt2(1,6,ik)*x(map_w2(6)+k) + pt2(2,6,ik-1)*x(map_w2(6)+k-1))
  t_e(2) = mt_inv(map_wt(2)+k)*p3t(1,2,ik)                              &
          *(pt2(2,5,ik)*x(map_w2(5)+k) + pt2(2,6,ik)  *x(map_w2(6)+k))

  div_u = div(1,1,ik)*x(map_w2(1)+k) + div(1,2,ik)*x(map_w2(2)+k) &
        + div(1,3,ik)*x(map_w2(3)+k) + div(1,4,ik)*x(map_w2(4)+k) &
        + div(1,5,ik)*x(map_w2(5)+k) + div(1,6,ik)*x(map_w2(6)+k)
  lhs(map_w3(1)+k) = m3(1,1,ik)*pressure(map_w3(1)+k) &
                   + sgn*(div_u + (t_e(1) + t_e(2))) + rhs_p(map_w3(1)+k)

end subroutine opt_apply_variable_hx_code_r_double

! R_SINGLE PRECISION
! ==================
subroutine opt_apply_variable_hx_code_r_single(cell,                    &
                                               nlayers,                 &
                                               lhs, x,                  &
                                               mt_inv,                  &
                                               pressure,                &
                                               ncell_3d_1,              &
                                               div,                     &
                                               ncell_3d_2,              &
                                               p3t,                     &
                                               ncell_3d_3,              &
                                               pt2,                     &
                                               ncell_3d_4,              &
                                               m3,                      &
                                               sgn,                     &
                                               rhs_p,                   &
                                               ndf_w3, undf_w3, map_w3, &
                                               ndf_w2, undf_w2, map_w2, &
                                               ndf_wt, undf_wt, map_wt)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell_3d_1, ncell_3d_2
  integer(kind=i_def),                    intent(in) :: ncell_3d_3, ncell_3d_4
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_single), dimension(undf_w2), intent(in)    :: x
  real(kind=r_single), dimension(undf_wt), intent(in)    :: mt_inv
  real(kind=r_single), dimension(undf_w3), intent(inout) :: lhs
  real(kind=r_single), dimension(undf_w3), intent(in)    :: pressure
  real(kind=r_single),                     intent(in)    :: sgn
  real(kind=r_single), dimension(undf_w3), intent(in)    :: rhs_p

  real(kind=r_single), dimension(1,6,ncell_3d_1), intent(in) :: div
  real(kind=r_single), dimension(2,6,ncell_3d_2), intent(in) :: pt2
  real(kind=r_single), dimension(1,2,ncell_3d_3), intent(in) :: p3t
  real(kind=r_single), dimension(1,1,ncell_3d_4), intent(in) :: m3

  ! Internal variables
  integer(kind=i_def)               :: k, ik
  real(kind=r_single), dimension(2) :: t_e
  real(kind=r_single)               :: div_u

  ! Compute D * u + P3t * Mt^-1 * ( Pt2 * u ) + rhs_p
  ! Hard wired optimisation for desired configuration (p=0 elements with pt2
  ! only acting on vertical components of u )
  k = 0
  ik = (cell-1)*nlayers + k + 1

  t_e(1) = mt_inv(map_wt(1)+k)*p3t(1,1,ik)                              &
          *(pt2(1,5,ik)*x(map_w2(5)+k) + pt2(1,6,ik)  *x(map_w2(6)+k))
  t_e(2) = mt_inv(map_wt(2)+k)*p3t(1,2,ik)                              &
          *(pt2(2,5,ik)*x(map_w2(5)+k) + pt2(1,5,ik+1)*x(map_w2(5)+k+1) &
          + pt2(2,6,ik)*x(map_w2(6)+k) + pt2(1,6,ik+1)*x(map_w2(6)+k+1))

  div_u = div(1,1,ik)*x(map_w2(1)+k) + div(1,2,ik)*x(map_w2(2)+k) &
        + div(1,3,ik)*x(map_w2(3)+k) + div(1,4,ik)*x(map_w2(4)+k) &
        + div(1,5,ik)*x(map_w2(5)+k) + div(1,6,ik)*x(map_w2(6)+k)
  lhs(map_w3(1)+k) = m3(1,1,ik)*pressure(map_w3(1)+k) &
                   + sgn*(div_u + (t_e(1) + t_e(2))) + rhs_p(map_w3(1)+k)

  do k = 1,nlayers-2
    ik = (cell-1)*nlayers + k + 1

    t_e(1) = mt_inv(map_wt(1)+k)*p3t(1,1,ik)                              &
            *(pt2(1,5,ik)*x(map_w2(5)+k) + pt2(2,5,ik-1)*x(map_w2(5)+k-1) &
            + pt2(1,6,ik)*x(map_w2(6)+k) + pt2(2,6,ik-1)*x(map_w2(6)+k-1))
    t_e(2) = mt_inv(map_wt(2)+k)*p3t(1,2,ik)                              &
            *(pt2(2,5,ik)*x(map_w2(5)+k) + pt2(1,5,ik+1)*x(map_w2(5)+k+1) &
            + pt2(2,6,ik)*x(map_w2(6)+k) + pt2(1,6,ik+1)*x(map_w2(6)+k+1))

    div_u = div(1,1,ik)*x(map_w2(1)+k) + div(1,2,ik)*x(map_w2(2)+k) &
          + div(1,3,ik)*x(map_w2(3)+k) + div(1,4,ik)*x(map_w2(4)+k) &
          + div(1,5,ik)*x(map_w2(5)+k) + div(1,6,ik)*x(map_w2(6)+k)
    lhs(map_w3(1)+k) = m3(1,1,ik)*pressure(map_w3(1)+k) &
                     + sgn*(div_u + (t_e(1) + t_e(2))) + rhs_p(map_w3(1)+k)
  end do

  k = nlayers-1
  ik = (cell-1)*nlayers + k + 1

  t_e(1) = mt_inv(map_wt(1)+k)*p3t(1,1,ik)                              &
          *(pt2(1,5,ik)*x(map_w2(5)+k) + pt2(2,5,ik-1)*x(map_w2(5)+k-1) &
          + pt2(1,6,ik)*x(map_w2(6)+k) + pt2(2,6,ik-1)*x(map_w2(6)+k-1))
  t_e(2) = mt_inv(map_wt(2)+k)*p3t(1,2,ik)                              &
          *(pt2(2,5,ik)*x(map_w2(5)+k) + pt2(2,6,ik)  *x(map_w2(6)+k))

  div_u = div(1,1,ik)*x(map_w2(1)+k) + div(1,2,ik)*x(map_w2(2)+k) &
        + div(1,3,ik)*x(map_w2(3)+k) + div(1,4,ik)*x(map_w2(4)+k) &
        + div(1,5,ik)*x(map_w2(5)+k) + div(1,6,ik)*x(map_w2(6)+k)
  lhs(map_w3(1)+k) = m3(1,1,ik)*pressure(map_w3(1)+k) &
                   + sgn*(div_u + (t_e(1) + t_e(2))) + rhs_p(map_w3(1)+k)

end subroutine opt_apply_variable_hx_code_r_single

end module opt_apply_variable_hx_kernel_mod
