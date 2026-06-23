!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the adjoint of LHS of the semi-implicit system
!!        for the horizontal velocity.
module adj_apply_mixed_u_operator_kernel_mod

use argument_mod,      only : arg_type,              &
                              GH_FIELD, GH_OPERATOR, &
                              GH_READ,               &
                              GH_WRITE,              &
                              GH_REAL, CELL_COLUMN
use constants_mod,     only : r_solver, i_def
use kernel_mod,        only : kernel_type
use fs_continuity_mod, only : W2, W3, W2h, W2v, W2broken

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: adj_apply_mixed_u_operator_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                       &
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W2broken),  & ! lhs_uv
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W2h),       & ! uv'
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W2v),       & ! w'
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),        & ! exner'
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W2, W2),    & ! Mu^{c,d}
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W2, W3),    & ! grad
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2)         & ! norm_u
       /)
  integer :: operates_on = CELL_COLUMN
  contains
  procedure, nopass :: adj_apply_mixed_u_operator_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_apply_mixed_u_operator_code

contains

!> @brief Compute the adjoint of the LHS of the semi-implicit system
!> @param[in]     cell          Horizontal cell index
!> @param[in]     nlayers       Number of layers
!> @param[in,out] lhs_uv        Mixed operator applied to the horizontal momentum equation
!> @param[in,out] wind_uv       Horizontal wind field
!> @param[in,out] wind_w        Vertical wind field
!> @param[in,out] exner         Exner pressure field
!> @param[in]     ncell1        Total number of cells for the mu_cd operator
!> @param[in]     mu_cd         Generalised mass matrix for the momentum equation
!> @param[in]     ncell2        Total number of cells for the grad operator
!> @param[in]     grad          Generalised gradient operator for the momentum equation
!> @param[in]     norm_u        Normalisation field for the momentum equation
!> @param[in]     ndf_w2hb      number of degrees of freedom per cell for the broken horizontal wind space
!> @param[in]     undf_w2hb     unique number of degrees of freedom for the broken horizontal wind space
!> @param[in]     map_w2hb      dofmap for the cell at the base of the column for the broken horizontal wind space
!> @param[in]     ndf_w2h       number of degrees of freedom per cell for the horizontal wind space
!> @param[in]     undf_w2h      unique number of degrees of freedom for the horizontal wind space
!> @param[in]     map_w2h       dofmap for the cell at the base of the column for the horizontal wind space
!> @param[in]     ndf_w2v       Number of degrees of freedom per cell for the vertical wind space
!> @param[in]     undf_w2v      Unique number of degrees of freedom for the vertical wind space
!> @param[in]     map_w2v       Dofmap for the cell at the base of the column for the vertical wind space
!> @param[in]     ndf_w3        Norm_umber of degrees of freedom per cell for the pressure space
!> @param[in]     ndf_w3        Unique number of degrees of freedom for the pressure space
!> @param[in]     map_w3        Dofmap for the cell at the base of the column for the pressure space
!> @param[in]     ndf_w2        Number of degrees of freedom per cell for the 3d wind space
!> @param[in]     undf_w2       Unique number of degrees of freedom for the 3d wind space
!> @param[in]     map_w2        Dofmap for the cell at the base of the column for the 3d wind space

subroutine adj_apply_mixed_u_operator_code(cell,                          &
                                           nlayers,                       &
                                           lhs_uv,                        &
                                           wind_uv, wind_w, exner,        &
                                           ncell1, mu_cd,                 &
                                           ncell2, grad,                  &
                                           norm_u,                        &
                                           ndf_w2hb, undf_w2hb, map_w2hb, &
                                           ndf_w2h, undf_w2h, map_w2h,    &
                                           ndf_w2v, undf_w2v, map_w2v,    &
                                           ndf_w3, undf_w3, map_w3,       &
                                           ndf_w2, undf_w2, map_w2)

  implicit none

  ! Arguments
  integer(kind=i_def),                      intent(in) :: cell, nlayers
  integer(kind=i_def),                      intent(in) :: ncell1, ncell2
  integer(kind=i_def),                      intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                      intent(in) :: undf_w2h, ndf_w2h
  integer(kind=i_def),                      intent(in) :: undf_w2hb, ndf_w2hb
  integer(kind=i_def),                      intent(in) :: undf_w2v, ndf_w2v
  integer(kind=i_def),                      intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def), dimension(ndf_w2hb), intent(in) :: map_w2hb
  integer(kind=i_def), dimension(ndf_w2h),  intent(in) :: map_w2h
  integer(kind=i_def), dimension(ndf_w2v),  intent(in) :: map_w2v
  integer(kind=i_def), dimension(ndf_w2),   intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3),   intent(in) :: map_w3

  ! Fields
  real(kind=r_solver), dimension(undf_w2hb), intent(inout) :: lhs_uv
  real(kind=r_solver), dimension(undf_w2h),  intent(inout) :: wind_uv
  real(kind=r_solver), dimension(undf_w2v),  intent(inout) :: wind_w
  real(kind=r_solver), dimension(undf_w2),   intent(in)    :: norm_u
  real(kind=r_solver), dimension(undf_w3),   intent(inout) :: exner

  ! Operators
  real(kind=r_solver), dimension(ncell1, ndf_w2, ndf_w2), intent(in) :: mu_cd
  real(kind=r_solver), dimension(ncell2, ndf_w2, ndf_w3), intent(in) :: grad

  ! Internal variables
  integer(kind=i_def) :: df, df2, ij, &
                         nm1, iw3,    &
                         iw2, iw2hb,  &
                         iw2h2, iw2v2

  ! Set up some useful shorthands for indices
  ij = (cell-1)*nlayers + 1
  nm1 = nlayers-1
  iw3 = map_w3(1)

  do df2 = ndf_w2v, 1, -1
    do df = ndf_w2h, 1, -1
      iw2hb = map_w2hb(df)
      iw2   = map_w2(df)
      iw2v2 = map_w2v(df2)
      wind_w(iw2v2:iw2v2+nm1) = wind_w(iw2v2:iw2v2+nm1) &
                              + norm_u(iw2:iw2+nm1)*  &
                                mu_cd(ij:ij+nm1, df, ndf_w2h+df2)*lhs_uv(iw2hb:iw2hb+nm1)
    end do
  end do

  do df2 = ndf_w2h, 1, -1
    do df = ndf_w2h, 1, -1
      iw2hb = map_w2hb(df)
      iw2   = map_w2(df)
      iw2h2 = map_w2h(df2)
      wind_uv(iw2h2:iw2h2+nm1) = wind_uv(iw2h2:iw2h2+nm1) &
                               + norm_u(iw2:iw2+nm1)*  &
                                 mu_cd(ij:ij+nm1, df, df2)*lhs_uv(iw2hb:iw2hb+nm1)
    end do
  end do

  do df = ndf_w2h, 1, -1
    iw2hb = map_w2hb(df)
    iw2   = map_w2(df)
    exner(iw3:iw3+nm1) = exner(iw3:iw3+nm1) &
                       - norm_u(iw2:iw2+nm1)   &
                         *grad(ij:ij+nm1, df, 1)*lhs_uv(iw2hb:iw2hb+nm1)
    lhs_uv(iw2hb:iw2hb+nm1) = 0.0_r_solver
  end do

end subroutine adj_apply_mixed_u_operator_code

end module adj_apply_mixed_u_operator_kernel_mod
