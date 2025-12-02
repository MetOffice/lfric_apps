!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to compute the vertical departure distances for FFSL
!> @details This code calculates the distance which is swept through a cell
!!          in the z direction during one time step of ffsl transport. The
!!          arrival point is the cell face and the departure point is calculated
!!          and stored as a field.
!!          The departure points are a dimensionless displacement, corresponding
!!          to the number of cells moved by a fluid parcel. The part of the wind
!!          corresponding to the fractional part of the departure points is
!!          also computed.

module vert_dep_dist_ffsl_kernel_mod

use argument_mod,                only : arg_type, GH_SCALAR,   &
                                        GH_FIELD, GH_REAL,     &
                                        GH_WRITE, GH_READ,     &
                                        CELL_COLUMN
use fs_continuity_mod,           only : W3, W2v
use constants_mod,               only : r_tran, i_def
use kernel_mod,                  only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vert_dep_dist_ffsl_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                  &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2v),    & ! dep_dist
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2v),    & ! frac_wind
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2v),    & ! wind
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3)      & ! detj
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vert_dep_dist_ffsl_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vert_dep_dist_ffsl_code

contains

!> @brief Computes the vertical departure distances using an Eulerian method for
!!        FFSL transport.
!> @param[in]     nlayers        The number of layers in the mesh
!> @param[in,out] dep_dist_z     Field with vertical departure distances
!> @param[in,out] frac_wind      Fractional part of the wind to be computed here
!> @param[in]     wind           The vertical wind
!> @param[in]     detj           Det(J) at W3: the volume of cells
!> @param[in]     ndf_w2v        Number of DoFs per cell for W2V
!> @param[in]     undf_w2v       Number of W2V DoFs in memory for this partition
!> @param[in]     map_w2v        Map of lowest-cell W2V DoFs
!> @param[in]     ndf_w3         Number of DoFs per cell for W3
!> @param[in]     undf_w3        Number of W3 DoFs in memory for this partition
!> @param[in]     map_w3         Map of lowest-cell W3 DoFs
subroutine vert_dep_dist_ffsl_code( nlayers,             &
                                    dep_dist_z,          &
                                    frac_wind,           &
                                    wind,                &
                                    detj,                &
                                    ndf_w2v,             &
                                    undf_w2v,            &
                                    map_w2v,             &
                                    ndf_w3,              &
                                    undf_w3,             &
                                    map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2v, ndf_w3
  integer(kind=i_def), intent(in)    :: undf_w2v, undf_w3
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  real(kind=r_tran),   intent(in)    :: detj(undf_w3)
  real(kind=r_tran),   intent(in)    :: wind(undf_w2v)
  real(kind=r_tran),   intent(inout) :: dep_dist_z(undf_w2v)
  real(kind=r_tran),   intent(inout) :: frac_wind(undf_w2v)

  integer(kind=i_def) :: j, k, cell_idx, max_num_cells, num_int_cells
  integer(kind=i_def) :: sign_flux, offset
  real(kind=r_tran)   :: flux_face_k, running_flux_k
  real(kind=r_tran)   :: frac_dep_dist, frac_flux_face_k

  ! Factor used to ensure consistency of fractional flux if
  ! there are any floating point rounding errors
  integer(kind=i_def) :: rounding_factor

  ! Set the bottom values
  frac_wind(map_w2v(1)) = 0.0_r_tran
  dep_dist_z(map_w2v(1)) = 0.0_r_tran

  do k = 1, nlayers - 1

    flux_face_k = wind(map_w2v(1)+k)

    ! Get an upper limit on the number of cells to step through
    if (flux_face_k >= 0.0_r_tran) then
      max_num_cells = k
      offset = 0
      sign_flux = 1
    else
      max_num_cells = nlayers - k
      offset = -1
      sign_flux = -1
    end if

    num_int_cells = 0
    running_flux_k = 0.0_r_tran

    ! Step backwards through flux to find departure point
    do j = 1, max_num_cells
      ! Index of cell to look at
      cell_idx = map_w3(1) + k - sign_flux*j + offset

      ! We have found integer number of cells, so exit do-loop
      if (running_flux_k + detj(cell_idx) > abs(flux_face_k)) EXIT

      ! Increment running values, if we have not exited loop
      running_flux_k = running_flux_k + detj(cell_idx)
      num_int_cells = num_int_cells + 1
    end do

    ! Set fractional flux. running_flux_k is now the integer flux
    frac_flux_face_k = sign_flux*(abs(flux_face_k) - running_flux_k)

    ! Determine fractional distance
    cell_idx = map_w3(1) + k - sign_flux*(num_int_cells + 1) + offset
    ! Cap cell_idx to prevent it exceeding the column
    cell_idx = MIN(MAX(cell_idx, map_w3(1)), map_w3(1) + nlayers - 1)
    frac_dep_dist = frac_flux_face_k / detj(cell_idx)

    ! Set the values of the output fields
    dep_dist_z(map_w2v(1)+k) = real(sign_flux*num_int_cells, r_tran) + frac_dep_dist

    ! If floating point rounding puts the dep_dist into the next integer value
    ! then set the fractional part back to zero
    rounding_factor = 1 - (abs(int(dep_dist_z(map_w2v(1)+k), i_def)) - num_int_cells)
    frac_wind(map_w2v(1)+k) = frac_flux_face_k*rounding_factor

  end do

  ! Set the top values
  frac_wind(map_w2v(1)+nlayers) = 0.0_r_tran
  dep_dist_z(map_w2v(1)+nlayers) = 0.0_r_tran

end subroutine vert_dep_dist_ffsl_code

end module vert_dep_dist_ffsl_kernel_mod
