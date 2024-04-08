!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to sort the vertical departure distances
!> @details This sorts the vertical departure distances, so that the vertical
!!          departure points are vertically increasing
!!          Note that this method loops upwards, so does not treat upward and
!!          downward winds symmetrically
module vert_deppt_sorting_kernel_mod

use argument_mod,                only : arg_type, GH_FIELD,    &
                                        GH_REAL, GH_READWRITE, &
                                        CELL_COLUMN
use fs_continuity_mod,           only : W2v
use constants_mod,               only : r_tran, i_def
use kernel_mod,                  only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vert_deppt_sorting_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                   &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, W2v)  & ! dep_dist
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vert_deppt_sorting_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vert_deppt_sorting_code

contains

!> @brief Kernel to cap vertical departure points
!> @param[in]     nlayers        The number of layers in the mesh
!> @param[in,out] dep_dist_z     Field with vertical departure distances
!> @param[in]     ndf_w2v        Number of DoFs per cell for W2V
!> @param[in]     undf_w2v       Number of W2V DoFs in memory for this partition
!> @param[in]     map_w2v        Map of lowest-cell W2V DoFs
subroutine vert_deppt_sorting_code( nlayers,             &
                                    dep_dist_z,          &
                                    ndf_w2v,             &
                                    undf_w2v,            &
                                    map_w2v )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  real(kind=r_tran),   intent(inout) :: dep_dist_z(undf_w2v)

  ! Internal variables
  integer(kind=i_def) :: k

  ! Loop through internal wind values
  do k = 1, nlayers - 2

    ! If the departure points have crossed, set the new departure point to be
    ! slightly above the old departure point
    if (dep_dist_z(map_w2v(1)+k+1) > dep_dist_z(map_w2v(1)+k) + 1.0_r_tran) then
      dep_dist_z(map_w2v(1)+k+1) =                                             &
        MAX(dep_dist_z(map_w2v(1)+k)                                           &
            + 1.0_r_tran + 1.0_r_tran / REAL(nlayers, r_tran),                 &
            REAL(k - nlayers, r_tran))
    end if

  end do

end subroutine vert_deppt_sorting_code

end module vert_deppt_sorting_kernel_mod
