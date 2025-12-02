!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Apply the semi-implicit mixed operator to the equation of state.

module apply_mixed_lp_operator_kernel_mod

use argument_mod,            only : arg_type,              &
                                    GH_FIELD, GH_OPERATOR, &
                                    GH_READ, GH_WRITE,     &
                                    GH_REAL, CELL_COLUMN
use constants_mod,           only : r_solver, i_def
use kernel_mod,              only : kernel_type
use fs_continuity_mod,       only : W3, Wtheta

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: apply_mixed_lp_operator_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                       &
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),        & ! lhs_exner
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),    & ! theta'
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),        & ! rho'
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),        & ! exner'
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W3),    & ! M3^exner
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W3),    & ! M3^rho
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, Wtheta) & ! P3t
       /)
  integer :: operates_on = CELL_COLUMN
  contains
  procedure, nopass :: apply_mixed_lp_operator_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: apply_mixed_lp_operator_code

contains
!> @brief Compute the LHS of the equation of state:
!!        lhs_exner = m3exner*exner - m3rho*rho - p3theta*theta
!! @param[in] cell Horizontal cell index
!! @param[in] nlayers Number of layers
!! @param[in,out] lhs_exner Mixed operator applied to the equation of state
!! @param[in] theta Potential temperature field
!! @param[in] rho Density field
!! @param[in] exner Exner pressure field
!! @param[in] ncell1 Total number of cells for the m3exner operator
!! @param[in] m3exner W3 mass matrix weighted by the reference Exner pressure
!! @param[in] ncell2 Total number of cells for the m3rho operator
!! @param[in] m3rho W3 mass matrix weighted by the reference density
!! @param[in] ncell3 Total number of cells for the p3theta operator
!! @param[in] p3theta Projection from Wtheta to W3 weighted by the reference
!!                    potential temperature
!! @param[in] ndf_w3 Number of degrees of freedom per cell for the density space
!! @param[in] undf_w3 Unique number of degrees of freedom for the density space
!! @param[in] map_w3 Dofmap for the cell at the base of the column for the
!!                   density space
!! @param[in] ndf_wt Number of degrees of freedom per cell for the potential temperature space
!! @param[in] undf_wt Unique number of degrees of freedom for the potential temperature space
!! @param[in] map_wt Dofmap for the cell at the base of the column for the potential temperature space
subroutine apply_mixed_lp_operator_code(cell,                    &
                                        nlayers,                 &
                                        lhs_exner,               &
                                        theta, rho, exner,       &
                                        ncell1, m3exner,         &
                                        ncell2, m3rho,           &
                                        ncell3, p3theta,         &
                                        ndf_w3, undf_w3, map_w3, &
                                        ndf_wt, undf_wt, map_wt)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell1, ncell2, ncell3
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  ! Fields
  real(kind=r_solver), dimension(undf_w3), intent(inout) :: lhs_exner
  real(kind=r_solver), dimension(undf_wt), intent(in)    :: theta
  real(kind=r_solver), dimension(undf_w3), intent(in)    :: rho, exner

  ! Operators
  real(kind=r_solver), dimension(ncell1, ndf_w3, ndf_w3), intent(in) :: m3exner
  real(kind=r_solver), dimension(ncell2, ndf_w3, ndf_w3), intent(in) :: m3rho
  real(kind=r_solver), dimension(ncell3, ndf_w3, ndf_wt), intent(in) :: p3theta

  ! Internal variables
  integer(kind=i_def)                    :: df, k, ik
  real(kind=r_solver), dimension(ndf_wt) :: t_e
  real(kind=r_solver), dimension(ndf_w3) :: r_e, p_e, lhs_e

  do k = 0, nlayers-1
    do df = 1, ndf_wt
      t_e(df) = theta(map_wt(df)+k)
    end do
    do df = 1, ndf_w3
      r_e(df) = rho(map_w3(df)+k)
      p_e(df) = exner(map_w3(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    ! lhs for this element
    lhs_e = matmul(m3exner(ik,:,:),p_e) - matmul(m3rho(ik,:,:),r_e)   &
            - matmul(p3theta(ik,:,:),t_e)
    do df = 1, ndf_w3
      lhs_exner(map_w3(df)+k) = lhs_e(df)
    end do

  end do

end subroutine apply_mixed_lp_operator_code

end module apply_mixed_lp_operator_kernel_mod
