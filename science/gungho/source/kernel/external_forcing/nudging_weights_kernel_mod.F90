!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes weights for scaling the spectral nudging increment.
!> @details This kernel computes the height-dependent weights used in spectral
!!          nudging, to allow the nudging to be targeted at certain vertical
!!          heights in the atmosphere.
!!          The weights are currently ramped linearly from 0 at the surface to
!!          2km, and then held constant above this height
!!          Only implemented for the lowest-order elements
module nudging_weights_kernel_mod

  use argument_mod,              only: arg_type,                               &
                                       GH_FIELD, GH_SCALAR,                    &
                                       GH_REAL, GH_INTEGER,                    &
                                       GH_READ, GH_WRITE,                      &
                                       CELL_COLUMN
  use constants_mod,             only: r_def, i_def
  use fs_continuity_mod,         only: W3
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: nudging_weights_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3),                         &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                               &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: nudging_weights_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: nudging_weights_code

contains

!> @brief Computes weights for scaling the spectral nudging increment.
!> @param[in]     nlayers             Number of layers
!> @param[in,out] weights             The weights to be computed
!> @param[in]     bottom_nudge_level  Level below which nudging is turned off
!> @param[in]     top_nudge_level     Level above which nudging is turned off
!> @param[in]     bottom_nudge_width  Width of the ramp at the bottom
!> @param[in]     top_nudge_width     Width of the ramp at the top
!> @param[in]     ndf                 Num DoFs per cell
!> @param[in]     undf                Num DoFs in this partition
!> @param[in]     map                 DoFmap for this column
subroutine nudging_weights_code(nlayers,                                       &
                                weights,                                       &
                                bottom_nudge_level,                            &
                                top_nudge_level,                               &
                                bottom_nudge_width,                            &
                                top_nudge_width,                               &
                                ndf,                                           &
                                undf,                                          &
                                map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf
  integer(kind=i_def), intent(in)    :: undf
  integer(kind=i_def), intent(in)    :: map(ndf)
  integer(kind=i_def), intent(in)    :: bottom_nudge_level
  integer(kind=i_def), intent(in)    :: top_nudge_level
  integer(kind=i_def), intent(in)    :: bottom_nudge_width
  integer(kind=i_def), intent(in)    :: top_nudge_width
  real(kind=r_def),    intent(inout) :: weights(undf)


  ! Local variables
  integer(kind=i_def) :: k, idx, nl

  nl = nlayers - 2 + ndf  ! nlayers for Wtheta, nlayers - 1 for W3
  idx = map(1)

  do k = 0, bottom_nudge_level-bottom_nudge_width-1
    weights(idx + k) = 0.0_r_def
  end do

  do k = bottom_nudge_level-bottom_nudge_width, bottom_nudge_level+bottom_nudge_width
    weights(idx + k) = real(k - (bottom_nudge_level - bottom_nudge_width), r_def) / &
                       real(MAX(1, 2*bottom_nudge_width), r_def)
  end do

  do k = bottom_nudge_level+bottom_nudge_width, top_nudge_level-top_nudge_width-1
    weights(idx + k) = 1.0_r_def
  end do

  do k = top_nudge_level-top_nudge_width+1, top_nudge_level+top_nudge_width
    weights(idx + k) = 1.0_r_def - real(k - (top_nudge_level - top_nudge_width), r_def) / &
                       real(MAX(1, 2*top_nudge_width), r_def)
  end do

  do k = top_nudge_level+top_nudge_width+1, nl
    weights(idx + k) = 0.0_r_def
  end do

end subroutine nudging_weights_code

end module nudging_weights_kernel_mod
