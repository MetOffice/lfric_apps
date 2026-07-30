!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes weights for scaling the spectral nudging increment.
!> @details This kernel computes the height-dependent weights used in spectral
!!          nudging, to allow the nudging to be targeted at certain vertical
!!          heights in the atmosphere.
!!          Only implemented for the lowest-order elements
module nudging_weights_kernel_mod

  use argument_mod,              only: arg_type,                               &
                                       GH_FIELD, GH_SCALAR,                    &
                                       GH_REAL, GH_INTEGER,                    &
                                       GH_READ, GH_WRITE,                      &
                                       ANY_DISCONTINUOUS_SPACE_1,              &
                                       CELL_COLUMN
  use constants_mod,             only: r_def, i_def, EPS
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
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  &
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
!> @param[in]     bottom_nudge_level  Level at and above which nudging is full
!!                                    strength
!> @param[in]     top_nudge_level     Level at and below which nudging is full
!!                                    strength
!> @param[in]     taper_bottom_level  Level at and below which nudging is 0
!> @param[in]     taper_top_level     Level at and above which nudging is 0
!> @param[in]     ndf                 Num DoFs per cell
!> @param[in]     undf                Num DoFs in this partition
!> @param[in]     map                 DoFmap for this column
subroutine nudging_weights_code(nlayers,                                       &
                                weights,                                       &
                                bottom_nudge_level,                            &
                                top_nudge_level,                               &
                                taper_bottom_level,                            &
                                taper_top_level,                               &
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
  integer(kind=i_def), intent(in)    :: taper_bottom_level
  integer(kind=i_def), intent(in)    :: taper_top_level
  real(kind=r_def),    intent(inout) :: weights(undf)


  ! Local variables
  integer(kind=i_def) :: k, idx, nl
  real(kind=r_def)    :: k_real, w3_offset
  real(kind=r_def)    :: taper_bottom, nudging_bottom, nudging_top, taper_top

  nl = nlayers + ndf - 2  ! nlayers for W3, nlayers+1 for Wtheta
  w3_offset = 0.5_r_def * real(2 - ndf, r_def)  ! 0.5 for W3, 0.0 for Wtheta
  idx = map(1)

  ! Convert levels to reals for comparison with k_real
  taper_bottom = real(taper_bottom_level, r_def)
  nudging_bottom = real(bottom_nudge_level, r_def)
  nudging_top = real(top_nudge_level, r_def)
  taper_top = real(taper_top_level, r_def)

  do k = 0, nl
    k_real = real(k, r_def) + w3_offset

    ! Below start of taper, weights = 0
    if (k_real < taper_bottom + EPS) then
      weights(idx + k) = 0.0_r_def

    ! Ramp up from 0 to 1 between taper_bottom and nudging_bottom
    else if (k_real < nudging_bottom + EPS) then
      weights(idx + k) = (k_real - taper_bottom) /                             &
                         max(1.0_r_def, nudging_bottom - taper_bottom)

    ! Full nudging between nudging_bottom and nudging_top
    else if (k_real < nudging_top + EPS) then
      weights(idx + k) = 1.0_r_def

    ! Ramp down from 1 to 0 between nudging_top and taper_top
    else if (k_real < taper_top + EPS) then
      weights(idx + k) = 1.0_r_def - (k_real - nudging_top) /                  &
                                     max(1.0_r_def, taper_top - nudging_top)

    ! Above taper_top, weights = 0
    else
      weights(idx + k) = 0.0_r_def

    end if
  end do

end subroutine nudging_weights_code

end module nudging_weights_kernel_mod
