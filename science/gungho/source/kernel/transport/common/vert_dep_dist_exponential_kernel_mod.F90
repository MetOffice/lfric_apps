!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to cap vertical departure points using an exponential method
!> @details This caps a vertical departure points field, to prevent departure
!!          points from leaving the model domain. Capped departure points are
!!          limited using an exponential method.
!!          This follows equation (12) of Wood, Staniforth and White (2009):
!!          z_dep = z_arr * exp(-dt*w/z_arr)
!!                = z_arr * exp(-dimensionless_dep_dist)

module vert_dep_dist_exponential_kernel_mod

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
type, public, extends(kernel_type) :: vert_dep_dist_exponential_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                   &
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, W2v)  & ! dep_dist
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vert_dep_dist_exponential_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vert_dep_dist_exponential_code

contains

!> @brief Kernel to cap vertical departure points using an exponential method
!> @param[in]     nlayers        The number of layers in the mesh
!> @param[in,out] dep_dist_z     Field with vertical departure distances
!> @param[in]     ndf_w2v        Number of DoFs per cell for W2V
!> @param[in]     undf_w2v       Number of W2V DoFs in memory for this partition
!> @param[in]     map_w2v        Map of lowest-cell W2V DoFs
subroutine vert_dep_dist_exponential_code( nlayers,             &
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

  ! Loop through all faces
  do k = 0, nlayers

    ! Check departure points haven't gone out of the domain bottom
    if (dep_dist_z(map_w2v(1)+k) > REAL(k, r_tran)) then
      ! z_dep = z_arr * exp(-dep_dist) where dep_dist is dimensionless
      ! The departure distance with dimensions of length is
      ! z_arr - z_dep
      ! = z_arr - z_arr * exp(-dep_dist)
      ! So the new dimensionless departure distance is found by multiplying by
      ! k / z_arr, giving the new dimensionless departure distance as
      ! dep_dist_new = k * (1 - exp(-dep_dist_old))
      dep_dist_z(map_w2v(1)+k) = REAL(k, r_tran) * (1.0_r_tran - EXP(-dep_dist_z(map_w2v(1)+k)))

    ! Check departure points haven't gone out of the domain top
    else if (dep_dist_z(map_w2v(1)+k) < REAL(k - nlayers, r_tran)) then
      ! This calculation values the above, with some signs changed
      dep_dist_z(map_w2v(1)+k) = (EXP(dep_dist_z(map_w2v(1)+k)) - 1.0_r_tran) * REAL(nlayers - k, r_tran)

    end if

  end do

end subroutine vert_dep_dist_exponential_code

end module vert_dep_dist_exponential_kernel_mod
