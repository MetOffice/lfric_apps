!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office GitHub Copilot Enterprise.
!-----------------------------------------------------------------------------

!> @brief Compute analytic Q32 contribution to the pressure RHS.
!> @details Computes:
!!          rhs = const_r*Q32_rho*rhs_u + const_t*Q32_theta*rhs_u.
module apply_analytic_q32_rhs_kernel_mod

  use argument_mod,      only : arg_type,                                      &
                                GH_FIELD, GH_OPERATOR,                         &
                                GH_READ, GH_WRITE,                             &
                                GH_SCALAR,                                     &
                                GH_REAL, CELL_COLUMN
  use constants_mod,     only : r_solver, i_def
  use kernel_mod,        only : kernel_type
  use fs_continuity_mod, only : W2, W3

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: apply_analytic_q32_rhs_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                        &
        arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),                          & ! rhs
        arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2),                          & ! rhs_u
        arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W2),                      & ! q32_rho
        arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W2),                      & ! q32_theta
        arg_type(GH_SCALAR,   GH_REAL, GH_READ),                               & ! const_r
        arg_type(GH_SCALAR,   GH_REAL, GH_READ)                                & ! const_t
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: apply_analytic_q32_rhs_code
  end type apply_analytic_q32_rhs_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: apply_analytic_q32_rhs_code

contains

!> @brief Compute analytic-elimination Q32 contribution for one column.
!> @param[in]     cell       Horizontal cell index
!> @param[in]     nlayers    Number of layers
!> @param[in,out] rhs        Output Q32 contribution field in W3
!> @param[in]     rhs_u      Input velocity RHS field in W2
!> @param[in]     ncell1     Total number of cells for q32_rho
!> @param[in]     q32_rho    Weighted projection operator from W2 to W3
!> @param[in]     ncell2     Total number of cells for q32_theta
!> @param[in]     q32_theta  Projection operator from W2 to W3 from
!!                          eliminated theta
!> @param[in]     const_r    tau_r*dt scaling constant applied to q32_rho
!> @param[in]     const_t    tau_t*dt scaling constant applied to q32_theta
!> @param[in]     ndf_w3     Number of degrees of freedom per cell for W3
!> @param[in]     undf_w3    Unique number of degrees of freedom for W3
!> @param[in]     map_w3     Dofmap for W3
!> @param[in]     ndf_w2     Number of degrees of freedom per cell for W2
!> @param[in]     undf_w2    Unique number of degrees of freedom for W2
!> @param[in]     map_w2     Dofmap for W2
subroutine apply_analytic_q32_rhs_code(cell,                                   &
                                       nlayers,                                &
                                       rhs,                                    &
                                       rhs_u,                                  &
                                       ncell1, q32_rho,                        &
                                       ncell2, q32_theta,                      &
                                       const_r, const_t,                       &
                                       ndf_w3, undf_w3, map_w3,                &
                                       ndf_w2, undf_w2, map_w2)

  implicit none

  integer(kind=i_def), intent(in)    :: cell, nlayers
  integer(kind=i_def), intent(in)    :: ncell1, ncell2
  integer(kind=i_def), intent(in)    :: undf_w3, ndf_w3
  integer(kind=i_def), intent(in)    :: undf_w2, ndf_w2
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)

  real(kind=r_solver), intent(inout) :: rhs(undf_w3)
  real(kind=r_solver), intent(in)    :: rhs_u(undf_w2)
  real(kind=r_solver), intent(in)    :: const_r
  real(kind=r_solver), intent(in)    :: const_t

  real(kind=r_solver), intent(in)    :: q32_rho(ncell1, ndf_w3, ndf_w2)
  real(kind=r_solver), intent(in)    :: q32_theta(ncell2, ndf_w3, ndf_w2)

  integer(kind=i_def) :: df1, df2, ik, i1, i2, nl

  nl = nlayers-1

  do df1 = 1, ndf_w3
    i1 = map_w3(df1)
    rhs(i1:i1+nl) = 0.0_r_solver
  end do

  ik = (cell-1)*nlayers + 1
  do df2 = 1, ndf_w2
    i2 = map_w2(df2)
    do df1 = 1, ndf_w3
      i1 = map_w3(df1)
      rhs(i1:i1+nl) = rhs(i1:i1+nl)                                            &
                    + const_r*q32_rho(ik:ik+nl,df1,df2)*rhs_u(i2:i2+nl)        &
                    + const_t*q32_theta(ik:ik+nl,df1,df2)*rhs_u(i2:i2+nl)
    end do
  end do

end subroutine apply_analytic_q32_rhs_code

end module apply_analytic_q32_rhs_kernel_mod