!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the backsubstitution for the wind field from the Schur
!!        complement preconditioner when using a split W2 field (W2h, W2v).
module schur_backsub_kernel_mod

  use argument_mod,      only : arg_type,              &
                                GH_FIELD, GH_OPERATOR, &
                                GH_READ, GH_INC,       &
                                GH_WRITE,              &
                                GH_REAL, CELL_COLUMN,  &
                                GH_SCALAR, GH_LOGICAL
  use constants_mod,     only : r_solver, i_def, l_def
  use fs_continuity_mod, only : W2h, W2v, W2, W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  ! Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: schur_backsub_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                       &
         arg_type(GH_FIELD,    GH_REAL,    GH_INC,   W2h),    & ! lhs_h
         arg_type(GH_FIELD,    GH_REAL,    GH_WRITE, W2v),    & ! lhs_v
         arg_type(GH_FIELD,    GH_REAL,    GH_READ,  W2),     & ! rhs
         arg_type(GH_FIELD,    GH_REAL,    GH_READ,  W3),     & ! exner_inc
         arg_type(GH_OPERATOR, GH_REAL,    GH_READ,  W2, W3), & ! div
         arg_type(GH_FIELD,    GH_REAL,    GH_READ,  W2),     & ! norm
         arg_type(GH_SCALAR,   GH_LOGICAL, GH_READ),          & ! lam
         arg_type(GH_FIELD,    GH_REAL,    GH_READ,  W2)      & ! mask
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: schur_backsub_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: schur_backsub_code

contains

!> @details Computes (lhs_h,lhs_v) = rhs + norm*Hb_inv*div(exner_inc)
!!          where lhs_h and lhs_v are the W2h and W2v components of a W2 space.
!!          The weighted div matrix maps from the pressure space (W3) space to
!!          the full lhs space (W2) and norm and Hb_inv are
!!          fields in the same space as and rhs.
!> @param[in]     cell      Horizontal cell index
!> @param[in]     nlayers   Number of layers
!> @param[in,out] lhs_h     Horizontal output lhs (A*x)_h
!> @param[in,out] lhs_v     Vertical output lhs (A*x)_v
!> @param[in]     rhs       3D rhs field for computing lhs from
!> @param[in]     exner_inc Pressure increment
!> @param[in]     ncell_3d  Total number of cells
!> @param[in]     div       Local matrix assembly form of the weighted diverence operator
!> @param[in]     norm      Normalisation field to scale output by
!> @param[in]     Hb_inv    Second field to scale output by
!> @param[in]     lam       Flag to apply lateral boundary conditions
!> @param[in]     mask      Boundary condition mask for lam domains
!> @param[in]     ndfh      Number of degrees of freedom per cell for the horizontal output field
!> @param[in]     undfh     Unique number of degrees of freedom  for the horizontal output field
!> @param[in]     maph      Dofmap for the cell at the base of the column for the horizontal output field
!> @param[in]     ndfv      Number of degrees of freedom per cell for the vertical output field
!> @param[in]     undfv     Unique number of degrees of freedom for the vertical output field
!> @param[in]     mapv      Dofmap for the cell at the base of the column for the vertical output field
!> @param[in]     ndf1      Number of degrees of freedom per cell for the norm, Hb_inv and rhs
!!                          fields (W2)
!> @param[in]     undf1     Unique number of degrees of freedom for the norm, Hb_inv and rhs
!!                          fields (W2)
!> @param[in]     map1      Dofmap for the cell at the base of the column for the norm, Hb_inv and rhs
!!                          fields (W2)
!> @param[in]     ndf2      Number of degrees of freedom per cell for the exner_inc field (W3)
!> @param[in]     undf2     Unique number of degrees of freedom for the exner_inc field (W3)
!> @param[in]     map2      Dofmap for the cell at the base of the column for the exner_inc field (W3)
subroutine schur_backsub_code(cell,              &
                              nlayers,           &
                              lhs_h, lhs_v,      &
                              rhs, exner_inc,    &
                              ncell_3d,          &
                              div,               &
                              norm,              &
                              lam, mask,         &
                              ndfh, undfh, maph, &
                              ndfv, undfv, mapv, &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)

  implicit none

  ! Arguments
  integer(kind=i_def),                   intent(in) :: cell, nlayers, ncell_3d
  integer(kind=i_def),                   intent(in) :: undfh, ndfh
  integer(kind=i_def),                   intent(in) :: undfv, ndfv
  integer(kind=i_def),                   intent(in) :: undf1, ndf1
  integer(kind=i_def),                   intent(in) :: undf2, ndf2
  integer(kind=i_def), dimension(ndf1),  intent(in) :: map1
  integer(kind=i_def), dimension(ndf2),  intent(in) :: map2
  integer(kind=i_def), dimension(ndfh),  intent(in) :: maph
  integer(kind=i_def), dimension(ndfv),  intent(in) :: mapv

  logical(kind=l_def), intent(in) :: lam

  real(kind=r_solver), dimension(undf2),              intent(in)    :: exner_inc
  real(kind=r_solver), dimension(undfh),              intent(inout) :: lhs_h
  real(kind=r_solver), dimension(undfv),              intent(inout) :: lhs_v
  real(kind=r_solver), dimension(ncell_3d,ndf1,ndf2), intent(in)    :: div
  real(kind=r_solver), dimension(undf1),              intent(in)    :: norm
  real(kind=r_solver), dimension(undf1),              intent(in)    :: rhs
  real(kind=r_solver), dimension(undf1),              intent(in)    :: mask

  ! Internal variables
  integer(kind=i_def) :: df, df2, ij, nl, i1, i2, ih, iv

  ij = (cell-1)*nlayers + 1
  nl = nlayers - 1

  ! Horizontal terms
  do df = 1, ndfh
    ih = maph(df)
    i1 = map1(df)
    lhs_h(ih:ih+nl) = lhs_h(ih:ih+nl) + 0.5_r_solver*rhs(i1:i1+nl)
    do df2 = 1, ndf2
      i2 = map2(df2)
      lhs_h(ih:ih+nl) = lhs_h(ih:ih+nl) &
                      + div(ij:ij+nl, df, df2)*exner_inc(i2:i2+nl)*norm(i1:i1+nl)
    end do
  end do

  ! Vertical terms
  iv = mapv(1)
  i1 = map1(ndfh+1)
  lhs_v(iv:iv+nl+1) = rhs(i1:i1+nl+1)

  do df = 1, ndfv
    iv = mapv(df)
    i1 = map1(ndfh+df)
    do df2 = 1, ndf2
      i2 = map2(df2)
      lhs_v(iv:iv+nl) = lhs_v(iv:iv+nl) &
                      + div(ij:ij+nl, ndfh + df, df2)*exner_inc(i2:i2+nl)*norm(i1:i1+nl)
    end do
  end do

  ! Vertical boundary conditions
  lhs_v(mapv(1))    = 0.0_r_solver
  lhs_v(mapv(2)+nl) = 0.0_r_solver

  ! Lateral boundary conditions
  if ( lam ) then
    do df = 1, ndfh
      ih = maph(df)
      i1 = map1(df)
      lhs_h(ih:ih+nl) = lhs_h(ih:ih+nl)*mask(i1:i1+nl)
    end do
    do df = 1, ndfv
      iv = mapv(df)
      i1 = map1(ndfh+df)
      lhs_v(iv:iv+nl) = lhs_v(iv:iv+nl)*mask(i1:i1+nl)
    end do
  end if

end subroutine schur_backsub_code

end module schur_backsub_kernel_mod
