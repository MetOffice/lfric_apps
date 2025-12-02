!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to Socrates for calculation of topographic parameters

module set_topography_kernel_mod

use argument_mod,  only : arg_type, CELL_COLUMN, &
                          GH_FIELD, GH_SCALAR, &
                          GH_REAL, GH_INTEGER, &
                          GH_READ, GH_WRITE, GH_READWRITE, &
                          ANY_DISCONTINUOUS_SPACE_1, &
                          ANY_DISCONTINUOUS_SPACE_2, &
                          ANY_DISCONTINUOUS_SPACE_3
use constants_mod, only : r_def, i_def, pi
use kernel_mod,    only : kernel_type

implicit none

private

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: set_topography_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/ &
    arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! grad_x_orog
    arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! grad_y_orog
    arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! slope_angle
    arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! slope_aspect
    arg_type(GH_FIELD, GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! skyview
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! horizon_angle
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_3), & ! horizon_aspect
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ), & ! n_horizon_angle
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ)  & ! n_horizon_layer
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: set_topography_code
end type

public :: set_topography_code

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers         Number of layers
!> @param[in]     grad_x_orog     Ancillary surface X gradient
!> @param[in]     grad_y_orog     Ancillary surface Y gradient
!> @param[in,out] slope_angle     Surface slope angle
!> @param[in,out] slope_aspect    Surface slope aspect
!> @param[in,out] skyview         Skyview factor
!> @param[in,out] horizon_angle   Horizon angles
!> @param[in,out] horizon_aspect  Horizon angle directions
!> @param[in]     n_horizon_angle Number of horizon angle directions
!> @param[in]     n_horizon_layer Number of layers of horizon angles
!> @param[in]     ndf_2d          No. of degrees of freedom per cell for 2D space
!> @param[in]     undf_2d         No. unique of degrees of freedom for 2D space
!> @param[in]     map_2d          Dofmap for cell at base of column for 2D space
!> @param[in]     ndf_h_ang       No. of degrees of freedom per cell for h_ang space
!> @param[in]     undf_h_ang      No. unique of degrees of freedom for h_ang space
!> @param[in]     map_h_ang       Dofmap for cell at base of column for h_ang space
!> @param[in]     ndf_h_asp       No. of degrees of freedom per cell for h_asp space
!> @param[in]     undf_h_asp      No. unique of degrees of freedom for h_asp space
!> @param[in]     map_h_asp       Dofmap for cell at base of column for h_asp space
subroutine set_topography_code(nlayers, &
                               grad_x_orog, grad_y_orog, &
                               slope_angle, slope_aspect, skyview, &
                               horizon_angle, horizon_aspect, &
                               n_horizon_angle, n_horizon_layer, &
                               ndf_2d, undf_2d, map_2d, &
                               ndf_h_ang, undf_h_ang, map_h_ang, &
                               ndf_h_asp, undf_h_asp, map_h_asp)

  use socrates_set_topography, only: set_topography

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, n_horizon_angle, n_horizon_layer
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: ndf_h_ang, undf_h_ang
  integer(i_def), intent(in) :: ndf_h_asp, undf_h_asp

  integer(i_def), dimension(ndf_2d), intent(in) :: map_2d
  integer(i_def), dimension(ndf_h_ang), intent(in) :: map_h_ang
  integer(i_def), dimension(ndf_h_asp), intent(in) :: map_h_asp

  real(r_def), dimension(undf_2d), intent(in) :: grad_x_orog, grad_y_orog
  real(r_def), dimension(undf_2d), intent(inout) :: slope_angle, slope_aspect
  real(r_def), dimension(undf_2d), intent(inout) :: skyview
  real(r_def), dimension(undf_h_ang), intent(inout) :: horizon_angle
  real(r_def), dimension(undf_h_asp), intent(inout) :: horizon_aspect

  ! Local variables
  integer :: n_profile = 1
  integer :: h_ang_1, h_ang_last
  integer :: h_asp_1, h_asp_last


  h_ang_1 = map_h_ang(1)
  h_ang_last = map_h_ang(1) + max(n_horizon_angle*n_horizon_layer - 1, 0)
  h_asp_1 = map_h_asp(1)
  h_asp_last = map_h_asp(1) + max(n_horizon_angle - 1, 0)

  call set_topography(n_profile, pi, &
                      n_horiz_ang    = n_horizon_angle, &
                      n_horiz_layer  = n_horizon_layer, &
                      grad_x         = grad_x_orog(map_2d(1):map_2d(1)), &
                      grad_y         = grad_y_orog(map_2d(1):map_2d(1)), &
                      horizon_angle  = horizon_angle(h_ang_1:h_ang_last), &
                      horizon_aspect = horizon_aspect(h_asp_1:h_asp_last), &
                      slope_angle    = slope_angle(map_2d(1):map_2d(1)), &
                      slope_aspect   = slope_aspect(map_2d(1):map_2d(1)), &
                      skyview        = skyview(map_2d(1):map_2d(1)))

end subroutine set_topography_code
end module set_topography_kernel_mod
