!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel that computes the maximum horizontal and 3D CFL by summing
!!        the maximum CFLs in a cell.
!> @details The maximum CFL for each cell is computed for the horizontal and
!!          in 3D. This can be used to compute the maximum CFL used in the
!!          substepping for MOL transport. This kernel assumes lowest order
!!          elements.

module calc_max_cfl_kernel_mod

use argument_mod,      only : arg_type,          &
                              GH_FIELD, GH_REAL, &
                              GH_WRITE, GH_READ, &
                              CELL_COLUMN
use constants_mod,     only : r_tran, i_def
use fs_continuity_mod, only : W3, W2h, W2v
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: calc_max_cfl_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                    &
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),     & ! sum_cfl
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),     & ! sum_horizontal_cfl
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2h),    & ! horizontal_cfl
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2v)     & ! vertical_cfl
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: calc_max_cfl_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: calc_max_cfl_code
contains

!> @brief Computes the sum of the maximum CFLs in each direction in a cell.
!> @param[in]     nlayers             Number of layers
!> @param[in,out] sum_cfl             Sum of 3D CFL in a cell
!> @param[in,out] sum_horizontal_cfl  Sum of horizontal CFL in a cell
!> @param[in]     horizontal_cfl      Horizontal CFL values
!> @param[in]     vertical_cfl        Vertical CFL values
!> @param[in]     ndf_w3              Number of degrees of freedom per cell
!> @param[in]     undf_w3             Number of unique degrees of freedom for the
!!                                    advective_update field
!> @param[in]     map_w3              Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2h             Number of degrees of freedom per cell for the horizontal wind fields
!> @param[in]     undf_w2h            Number of unique degrees of freedom for the horizontal wind fields
!> @param[in]     map_w2h             Dofmap for the cell at the base of the column for the horizontal wind fields
!> @param[in]     ndf_w2v             Number of degrees of freedom per cell for the vertical wind fields
!> @param[in]     undf_w2v            Number of unique degrees of freedom for the vertical wind fields
!> @param[in]     map_w2v             Dofmap for the cell at the base of the column for the vertical wind fields
subroutine calc_max_cfl_code( nlayers,              &
                              sum_cfl,              &
                              sum_horizontal_cfl,   &
                              horizontal_cfl,       &
                              vertical_cfl,         &
                              ndf_w3,               &
                              undf_w3,              &
                              map_w3,               &
                              ndf_w2h,              &
                              undf_w2h,             &
                              map_w2h,              &
                              ndf_w2v,              &
                              undf_w2v,             &
                              map_w2v               &
                              )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                     :: nlayers
  integer(kind=i_def), intent(in)                     :: ndf_w3, ndf_w2h, ndf_w2v
  integer(kind=i_def), intent(in)                     :: undf_w3, undf_w2h, undf_w2v
  integer(kind=i_def), dimension(ndf_w3), intent(in)  :: map_w3
  integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
  integer(kind=i_def), dimension(ndf_w2v), intent(in) :: map_w2v

  real(kind=r_tran), dimension(undf_w3), intent(inout) :: sum_cfl
  real(kind=r_tran), dimension(undf_w3), intent(inout) :: sum_horizontal_cfl
  real(kind=r_tran), dimension(undf_w2h), intent(in)   :: horizontal_cfl
  real(kind=r_tran), dimension(undf_w2v), intent(in)   :: vertical_cfl

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_tran)    :: max_cfl_h, max_cfl_v

  ! W2h dof map
  !
  !    4
  ! 1     3
  !    2
  !
  ! W2v dof map
  !
  !   2
  !
  !   1

  do k = 0, nlayers - 1
    max_cfl_h = max( abs(horizontal_cfl(map_w2h(1)+k)), abs(horizontal_cfl(map_w2h(3)+k)) ) &
              + max( abs(horizontal_cfl(map_w2h(2)+k)), abs(horizontal_cfl(map_w2h(4)+k)) )
    max_cfl_v = max( abs(vertical_cfl(map_w2v(1)+k)), abs(vertical_cfl(map_w2v(2)+k)) )
    sum_cfl(map_w3(1)+k) = max_cfl_h + max_cfl_v
    sum_horizontal_cfl(map_w3(1)+k) = max_cfl_h
  end do

end subroutine calc_max_cfl_code

end module calc_max_cfl_kernel_mod
