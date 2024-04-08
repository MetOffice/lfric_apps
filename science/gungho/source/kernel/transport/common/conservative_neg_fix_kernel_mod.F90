!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! TODO: this should use r_tran diagonal mass matrix for W3
!-------------------------------------------------------------------------------

!> @brief A negative fix that preserves columnar mass
!> @details Only written for the lowest-order function space.

module conservative_neg_fix_kernel_mod
use argument_mod,            only : arg_type, GH_FIELD, GH_READ,     &
                                    CELL_COLUMN, GH_REAL, GH_READWRITE
use fs_continuity_mod,       only : W3
use constants_mod,           only : r_tran, i_def, EPS_R_TRAN
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: conservative_neg_fix_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                    &
       arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3), & ! field
       arg_type(GH_FIELD,    GH_REAL, GH_READ,      W3)  & ! diag mass matrix
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: conservative_neg_fix_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public conservative_neg_fix_code

contains

!> @brief A negative fix that preserves columnar mass
!! @param[in]     nlayers     Number of layers
!! @param[in,out] field       Density field to be adjusted
!! @param[in,out] mm_w3_diag  Diagonal W3 mass matrix field
!! @param[in]     ndf_w3      Num of DoFs per cell for the density field
!! @param[in]     undf_w3     Num of DoFs per partition for the density field
!! @param[in]     map_w3      Base cell DoF-map for the density field
subroutine conservative_neg_fix_code(nlayers,                 &
                                     field, mm_w3_diag,       &
                                     ndf_w3, undf_w3, map_w3  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_w3, ndf_w3
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  real(kind=r_tran),   intent(inout) :: field(undf_w3)
  real(kind=r_tran),   intent(in)    :: mm_w3_diag(undf_w3)

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_tran)   :: old_mass, new_mass, min_value, scaling
  real(kind=r_tran)   :: field_local(nlayers)
  real(kind=r_tran)   :: mm_local(nlayers)

  ! Extract values to put into local array
  do k = 0, nlayers-1
    field_local(k+1) = field(map_w3(1)+k)
    mm_local(k+1) = mm_w3_diag(map_w3(1)+k)
  end do

  min_value = minval(field_local)

  ! Only do adjustment if min value is negative
  if (min_value < 0.0_r_tran) then
    old_mass = 0.0_r_tran
    new_mass = 0.0_r_tran

    ! Calculate old mass, adjust field and calculate new mass within same loop
    do k = 1, nlayers
      old_mass = old_mass + mm_local(k)*field_local(k)
      field_local(k) = max(field_local(k), 0.0_r_tran)
      new_mass = new_mass + mm_local(k)*field_local(k)
    end do

    ! Conservative scaling factor
    old_mass = max(0.0_r_tran, old_mass)
    new_mass = max(EPS_R_TRAN, new_mass)
    scaling = old_mass / new_mass

    ! Return local values to field object
    do k = 0, nlayers - 1
      field(map_w3(1)+k) = scaling * field_local(k+1)
    end do
  end if

end subroutine conservative_neg_fix_code

end module conservative_neg_fix_kernel_mod