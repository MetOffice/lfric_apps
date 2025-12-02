!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Apply the semi-implicit mixed operator to the momentum equation.
module apply_mixed_lu_operator_kernel_mod

use argument_mod,      only : arg_type,              &
                              GH_FIELD, GH_OPERATOR, &
                              GH_READ, GH_INC,       &
                              GH_REAL, CELL_COLUMN
use constants_mod,     only : r_solver, i_def
use kernel_mod,        only : kernel_type
use fs_continuity_mod, only : W2, W3, Wtheta

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: apply_mixed_lu_operator_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/                       &
       arg_type(GH_FIELD,    GH_REAL, GH_INC,  W2),         & ! lhs_u
       arg_type(GH_FIELD,    GH_REAL, GH_READ, W2),         & ! u'
       arg_type(GH_FIELD,    GH_REAL, GH_READ, Wtheta),     & ! theta'
       arg_type(GH_FIELD,    GH_REAL, GH_READ, W3),         & ! exner'
       arg_type(GH_OPERATOR, GH_REAL, GH_READ, W2, W2),     & ! Mu^{c,d}
       arg_type(GH_OPERATOR, GH_REAL, GH_READ, W2, Wtheta), & ! P2theta
       arg_type(GH_OPERATOR, GH_REAL, GH_READ, W2, W3),     & ! grad
       arg_type(GH_FIELD,    GH_REAL, GH_READ, W2)          & ! norm_u
       /)
  integer :: operates_on = CELL_COLUMN
  contains
  procedure, nopass :: apply_mixed_lu_operator_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: apply_mixed_lu_operator_code

contains

!> @brief Compute the LHS of the momentum equation lhs_u = norm_u*(Mu*u - P2t*t - grad*p).
!! @param[in] cell Horizontal cell index
!! @param[in] nlayers number of layers
!! @param[in,out] lhs_u Mixed operator applied to the momentum equation
!! @param[in] wind Wind field
!! @param[in] theta Potential temperature field
!! @param[in] exner Exner pressure field
!! @param[in] ncell1 Total number of cells for the mu_cd operator
!! @param[in] mu_cd Generalised mass matrix for the momentum equation
!! @param[in] ncell2 Total number of cells for the p2theta operator
!! @param[in] p2theta Generalised projection matrix from Wtheta to W2
!! @param[in] ncell3 Total number of cells for the grad operator
!! @param[in] grad Generalised gradient operator for the momentum equation
!! @param[in] norm_u Normalisation field for the momentum equation
!! @param[in] ndf_w2 number of degrees of freedom per cell for the wind space
!! @param[in] undf_w2 Unique number of degrees of freedom for the wind space
!! @param[in] map_w2 Dofmap for the cell at the base of the column for the wind space
!! @param[in] ndf_wt number of degrees of freedom per cell for the potential
!!                   temperature space
!! @param[in] undf_wt Unique number of degrees of freedom for the potential
!!                    temperature space
!! @param[in] map_wt Dofmap for the cell at the base of the column for the
!!                   potential temperature space
!! @param[in] ndf_w3 norm_umber of degrees of freedom per cell for the pressure space
!! @param[in] undf_w3 Unique number of degrees of freedom for the pressure space
!! @param[in] map_w3 Dofmap for the cell at the base of the column for the pressure space

subroutine apply_mixed_lu_operator_code(cell,                    &
                                        nlayers,                 &
                                        lhs_u,                   &
                                        wind, theta, exner,      &
                                        ncell1, mu_cd,           &
                                        ncell2, P2theta,         &
                                        ncell3, grad,            &
                                        norm_u,                  &
                                        ndf_w2, undf_w2, map_w2, &
                                        ndf_wt, undf_wt, map_wt, &
                                        ndf_w3, undf_w3, map_w3)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell1, ncell2, ncell3
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  ! Fields
  real(kind=r_solver), dimension(undf_w2), intent(inout) :: lhs_u
  real(kind=r_solver), dimension(undf_w2), intent(in)    :: wind, norm_u
  real(kind=r_solver), dimension(undf_wt), intent(in)    :: theta
  real(kind=r_solver), dimension(undf_w3), intent(in)    :: exner

  ! Operators
  real(kind=r_solver), dimension(ncell1, ndf_w2, ndf_w2), intent(in) :: mu_cd
  real(kind=r_solver), dimension(ncell2, ndf_w2, ndf_wt), intent(in) :: p2theta
  real(kind=r_solver), dimension(ncell3, ndf_w2, ndf_w3), intent(in) :: grad

  ! Internal variables
  integer(kind=i_def)                       :: k, ik, df1, df2
  real(kind=r_solver), dimension(0:nlayers) :: lhs_t

    do df1 = 1, ndf_w2
      lhs_t = 0.0_r_solver

      do df2 = 1, ndf_w2
        do k = 0, nlayers-1
          ik = (cell-1)*nlayers + k + 1
          lhs_t(k) = lhs_t(k) + mu_cd(ik,df1,df2) * wind(map_w2(df2)+k)
        end do
      end do

      do df2 = 1, ndf_wt
        do k = 0, nlayers-1
          ik = (cell-1)*nlayers + k + 1
          lhs_t(k) = lhs_t(k) - p2theta(ik,df1,df2) * theta(map_wt(df2)+k)
        end do
      end do

      do df2 = 1, ndf_w3
        do k = 0, nlayers-1
          ik = (cell-1)*nlayers + k + 1
          lhs_t(k) = lhs_t(k) - grad(ik,df1,df2) * exner(map_w3(df2)+k)
        end do
      end do

      do k = 0, nlayers-1
        lhs_u(map_w2(df1)+k) =  lhs_u(map_w2(df1)+k) + &
            norm_u(map_w2(df1)+k) * lhs_t(k)
      end do
    end do

end subroutine apply_mixed_lu_operator_code

end module apply_mixed_lu_operator_kernel_mod
