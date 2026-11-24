!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Adjoint of the application of the semi-implicit mixed operator
!>        to the equation of state.

module adj_apply_elim_mixed_lp_operator_kernel_mod

  use argument_mod,      only : arg_type,              &
                                GH_FIELD, GH_OPERATOR, &
                                GH_READ, GH_INC,       &
                                GH_READWRITE,          &
                                GH_SCALAR,             &
                                GH_REAL, CELL_COLUMN
  use constants_mod,     only : r_solver, i_def
  use kernel_mod,        only : kernel_type
  use fs_continuity_mod, only : W2, W3, Wtheta

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------

  type, public, extends(kernel_type) :: adj_apply_elim_mixed_lp_operator_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                            &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3),         & ! lhs_exner
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, Wtheta),     & ! theta'
         arg_type(GH_FIELD,    GH_REAL, GH_INC,       W2),         & ! u'
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3),         & ! exner'
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3, W3),     & ! M3^exner
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3, W2),     & ! Q32_rho
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3, Wtheta), & ! P3theta
         arg_type(GH_SCALAR,   GH_REAL, GH_READ)                   & ! const_r
         /)
    integer :: operates_on = CELL_COLUMN
  contains
  procedure, nopass :: adj_apply_elim_mixed_lp_operator_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: adj_apply_elim_mixed_lp_operator_code

contains

!> @brief Compute the adjoint of the LHS of the equation of state.
!> @param[in]     cell      Horizontal cell index
!> @param[in]     nlayers   Number of layers
!> @param[in,out] lhs_exner Mixed operator applied to the equation of state
!> @param[in,out] theta     Potential temperature field
!> @param[in,out] u         Wind field
!> @param[in,out] exner     Exner pressure field
!> @param[in]     ncell1    Total number of cells for the m3exner operator
!> @param[in]     m3exner   W3 mass matrix weighted by the reference Exner pressure
!> @param[in]     ncell2    Total number of cells for the m3rho operator
!> @param[in]     q32_rho   Projection matrix from W2 to W3 (density
!!                          contribution to Q32; the theta contribution is
!!                          carried instead through the p3theta*theta' term,
!!                          since theta' already carries the tau_t*dt scaling
!!                          when discretely eliminated -- see call sites)
!> @param[in]     ncell3    Total number of cells for the p3theta operator
!> @param[in]     p3theta   Projection from Wtheta to W3 weighted by the reference
!!                          potential temperature
!> @param[in]     const_r   tau_r*dt scaling constant applied to the Q32
!!                          (density) contribution to the equation of state
!> @param[in]     ndf_w3    Number of degrees of freedom per cell for the density space
!> @param[in]     undf_w3   Unique number of degrees of freedom for the density space
!> @param[in]     map_w3    Dofmap for the cell at the base of the column for the
!!                          density space
!> @param[in]     ndf_wt    Number of degrees of freedom per cell for the potential temperature space
!> @param[in]     undf_wt   Unique number of degrees of freedom for the potential temperature space
!> @param[in]     map_wt    Dofmap for the cell at the base of the column for the potential temperature space
!> @param[in]     ndf_w2    Number of degrees of freedom per cell for the wind space
!> @param[in]     undf_w2   Unique number of degrees of freedom for the wind space
!> @param[in]     map_w2    Dofmap for the cell at the base of the column for the
!!                          wind space
subroutine adj_apply_elim_mixed_lp_operator_code(cell,                    &
                                                 nlayers,                 &
                                                 lhs_exner,               &
                                                 theta, u, exner,         &
                                                 ncell1, m3exner,         &
                                                 ncell2, q32_rho,         &
                                                 ncell3, p3theta,         &
                                                 const_r,                 &
                                                 ndf_w3, undf_w3, map_w3, &
                                                 ndf_wt, undf_wt, map_wt, &
                                                 ndf_w2, undf_w2, map_w2)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell1, ncell2, ncell3
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  ! Fields
  real(kind=r_solver), dimension(undf_w3), intent(inout) :: lhs_exner
  real(kind=r_solver), dimension(undf_wt), intent(inout) :: theta
  real(kind=r_solver), dimension(undf_w3), intent(inout) :: exner
  real(kind=r_solver), dimension(undf_w2), intent(inout) :: u
  real(kind=r_solver),                     intent(in)    :: const_r

  ! Operators
  real(kind=r_solver), dimension(ncell1, ndf_w3, ndf_w3), intent(in) :: m3exner
  real(kind=r_solver), dimension(ncell2, ndf_w3, ndf_w2), intent(in) :: q32_rho
  real(kind=r_solver), dimension(ncell3, ndf_w3, ndf_wt), intent(in) :: p3theta

  ! Internal variables
  integer(kind=i_def)                    :: df, k, ik
  real(kind=r_solver), dimension(ndf_wt) :: t_e
  real(kind=r_solver), dimension(ndf_w3) :: p_e, lhs_e
  real(kind=r_solver), dimension(ndf_w2) :: u_e

  do k = nlayers-1, 0, -1

    ik = (cell-1)*nlayers + k + 1

    ! LHS for this element
    do df = ndf_w3, 1, -1
      lhs_e(df) = lhs_exner(map_w3(df)+k)
      lhs_exner(map_w3(df)+k) = 0.0_r_solver
    end do

    u_e = const_r*matmul(transpose(q32_rho(ik,:,:)), lhs_e)
    t_e = -matmul(transpose(p3theta(ik,:,:)), lhs_e)
    p_e = matmul(transpose(m3exner(ik,:,:)), lhs_e)

    do df = ndf_w2, 1, -1
      u(map_w2(df)+k) = u(map_w2(df)+k) + u_e(df)
    end do

    do df = ndf_w3, 1, -1
      exner(map_w3(df)+k) = exner(map_w3(df)+k) + p_e(df)
    end do

    do df = ndf_wt, 1, -1
      theta(map_wt(df)+k) = theta(map_wt(df)+k) + t_e(df)
    end do

  end do

end subroutine adj_apply_elim_mixed_lp_operator_code

end module adj_apply_elim_mixed_lp_operator_kernel_mod
