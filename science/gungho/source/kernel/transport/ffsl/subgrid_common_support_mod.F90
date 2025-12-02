!------------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Common routines for horizontal and vertical subgrid reconstructions
!!          for use in FFSL
!------------------------------------------------------------------------------
module subgrid_common_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN
use transport_enumerated_types_mod, only: monotone_strict,                     &
                                          monotone_relaxed,                    &
                                          monotone_positive,                   &
                                          monotone_qm_pos,                     &
                                          monotone_none
implicit none

private

public :: monotonic_edge
public :: subgrid_quadratic_recon

contains

  !----------------------------------------------------------------------------
  !> @brief  Applies monotonicity constraints to the edge values needed
  !!         for the quadratic subgrid reconstruction (i.e. PPM or Nirvana).
  !!         If a monotonic limiter is specified, this generally involves
  !!         bounding the edge values by the field values that the edge lies
  !!         between.
  !!
  !> @param[in]     field          An array of field values and their neighbours
  !> @param[in]     monotone       Monotone option to ensure no over/undershoots
  !> @param[in]     min_val        Minimum value to enforce edge value to be for
  !!                               quasi-monotone limiter
  !> @param[in,out] edge_left      Field value at left edge of cell
  !> @param[in,out] edge_right     Field value at right edge of cell
  !> @param[in]     order          Order of the reconstruction
  !> @param[in]     bottom_layer   Index of layer to apply monotonicity from
  !> @param[in]     nlayers        Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine monotonic_edge(field, monotone, min_val, edge_left, edge_right,   &
                            order, bottom_layer, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: bottom_layer
    integer(kind=i_def), intent(in)    :: order
    real(kind=r_tran),   intent(in)    :: field(nlayers,2*order+1)
    integer(kind=i_def), intent(in)    :: monotone
    real(kind=r_tran),   intent(in)    :: min_val
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers), edge_right(nlayers)

    real(kind=r_tran)   :: t1(nlayers)
    integer(kind=i_def) :: i, b, t

    i = order + 1
    b = bottom_layer
    t = nlayers

    select case (monotone)
    case (monotone_relaxed, monotone_strict)
      ! Ensure that edge values lie between those of neighbouring cells
      edge_left(b:t) = MIN(                                                    &
          MAX(field(b:t,i-1), field(b:t,i)),                                   &
          MAX(edge_left(b:t), MIN(field(b:t,i-1), field(b:t,i)))               &
      )
      edge_right(b:t) = MIN(                                                   &
          MAX(field(b:t,i), field(b:t,i+1)),                                   &
          MAX(edge_right(b:t), MIN(field(b:t,i), field(b:t,i+1)))              &
      )

    case (monotone_positive)
      ! Ensure that edges also lie above specified minimum value
      edge_left(b:t)  = MAX(edge_left(b:t), 0.0_r_tran)
      edge_right(b:t) = MAX(edge_right(b:t), 0.0_r_tran)

    case (monotone_qm_pos)
      ! Only bound edges when the central field is outside its neighbours
      ! The t1 variable is a test for this
      ! This is 1 if the cell's value is outside of its neighbours', 0 otherwise
      t1(b:t) = 0.5_r_tran - SIGN(0.5_r_tran,                                  &
          (field(b:t,i) - field(b:t,i-1)) * (field(b:t,i) - field(b:t,i+1))    &
      )

      ! Ensure that edge values lie between those of neighbouring cells
      edge_left(b:t) = t1(b:t)*edge_left(b:t) + (1.0_r_tran - t1(b:t))*MAX(    &
          MIN(field(b:t,i-1), field(b:t,i)),                                   &
          MIN(edge_left(b:t), MAX(field(b:t,i-1), field(b:t,i)))               &
      )
      edge_right(b:t) = t1(b:t)*edge_right(b:t) + (1.0_r_tran - t1(b:t))*MAX(  &
          MIN(field(b:t,i), field(b:t,i+1)),                                   &
          MIN(edge_right(b:t), MAX(field(b:t,i), field(b:t,i+1)))              &
      )

      ! Ensure that edges also lie above specified minimum value
      edge_left(b:t)  = MAX(edge_left(b:t), min_val)
      edge_right(b:t) = MAX(edge_right(b:t), min_val)
    end select

  end subroutine monotonic_edge

  !----------------------------------------------------------------------------
  !> @brief  Returns the horizontal subgrid quadratic reconstruction.
  !!         This can be used to compute the flux as:
  !!         flux = u * reconstruction
  !!         The reconstruction depends upon the cell edge values, the field
  !!         value, and the departure distance. The reconstruction is
  !!         third-order, and is based on the quadratic subgrid reconstruction
  !!         of PPM and Nirvana.
  !!         The reconstruction corresponds to the integral of the field in the
  !!         interval (0, 1) between 1-dep and 1. In the absence of a monotonic
  !!         limiter, the quadratic is determined by the constraints:
  !!         (a) the integral of the quadratic from 0 to 1 is equal to the
  !!             integral of the field value from 0 to 1 (which is just the
  !!             field value) so that mass is conserved
  !!         (b) the quadratic should equal the edge reconstruction at the left
  !!             (upwind) and right (downwind) edges of the cell. Monotonic
  !!             limiters may relax this constraint to ensure that the quadratic
  !!             remains bounded.
  !!         Note that the departure distance is always the absolute value of
  !!         the fractional part of the distance. This ensures that the left
  !!         edge is always upwind, and the reconstruction in this routine does
  !!         not need to consider the sign of the wind.
  !!
  !> @param[in,out] reconstruction The quadratic reconstruction
  !> @param[in]     dep            The absolute value of the fractional
  !!                               departure distance for the cell
  !> @param[in]     field          Field value in the cell
  !> @param[in]     edge_left      Field value at left edge of cell
  !> @param[in]     edge_right     Field value at right edge of cell
  !> @param[in]     monotone       Monotone option to ensure no over/undershoots
  !> @param[in]     order          Order of the reconstruction
  !> @param[in]     nlayers        Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine subgrid_quadratic_recon(reconstruction, dep, field, edge_left,    &
                                     edge_right, monotone, order, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: order
    real(kind=r_tran),   intent(inout) :: reconstruction(nlayers)
    real(kind=r_tran),   intent(in)    :: dep(nlayers)
    real(kind=r_tran),   intent(in)    :: field(nlayers,2*order+1)
    real(kind=r_tran),   intent(in)    :: edge_left(nlayers)
    real(kind=r_tran),   intent(in)    :: edge_right(nlayers)
    integer(kind=i_def), intent(in)    :: monotone

    real(kind=r_tran) :: cm(nlayers), cc(nlayers), cp(nlayers)
    real(kind=r_tran) :: a1(nlayers), a2(nlayers)
    real(kind=r_tran) :: tau(nlayers), sigma(nlayers), upsilon(nlayers)
    real(kind=r_tran) :: ts(nlayers), tl(nlayers), tr(nlayers)

    integer(kind=i_def) :: i

    i = order + 1

    select case (monotone)
    ! No limiter ---------------------------------------------------------------
    case (monotone_none)
      cp = 1.0_r_tran - 2.0_r_tran*dep + dep**2
      cc = 3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = -dep + dep**2

    ! Positivity ---------------------------------------------------------------
    case (monotone_positive)
      ! Get quadratic coefficients
      a1 = -4.0_r_tran*edge_left - 2.0_r_tran*edge_right + 6.0_r_tran*field(:,i)
      a2 = 3.0_r_tran*edge_left + 3.0_r_tran*edge_right - 6.0_r_tran*field(:,i)

      ! This is 1 if the turning point would be in the cell, 0 otherwise
      tau = 0.5_r_tran + SIGN(0.5_r_tran,                                      &
          (2.0_r_tran*edge_left + edge_right - 3.0_r_tran*field(:,i))          &
          * (edge_left + 2.0_r_tran*edge_right - 3.0_r_tran*field(:,i))        &
      )

      ! This is the turning point multiplied by a common denominator, to avoid
      ! potential division by zero. If the turning point is negative, this is 0
      ! while if the turning point is positive this is 1
      sigma = 0.5_r_tran + SIGN(0.5_r_tran,                                    &
          a2*(4.0_r_tran*edge_left*a2 - a1**2)                                 &
      )

      ! Combine the switches, so this is 1 if the turning point is in the cell
      ! and is negative
      tau(:) = tau(:)*(1 - sigma(:))

      ! If tau is 1, revert to constant reconstruction. Otherwise do normal
      ! quadratic reconstruction
      cp = (1.0_r_tran - tau)*(1.0_r_tran - 2.0_r_tran*dep + dep**2)
      cc = (1.0_r_tran - tau)*(3.0_r_tran*dep - 2.0_r_tran*dep**2) + tau
      cm = (1.0_r_tran - tau)*(-dep + dep**2)

    ! Strict limiter -----------------------------------------------------------
    case (monotone_strict)
      ! This is 1 if the turning point would be in the cell, 0 otherwise
      tau = 0.5_r_tran + SIGN(0.5_r_tran,                                      &
          (2.0_r_tran*edge_left + edge_right - 3.0_r_tran*field(:,i))          &
          * (edge_left + 2.0_r_tran*edge_right - 3.0_r_tran*field(:,i))        &
      )

      ! If the tau switch is 1, we revert to constant reconstruction. Otherwise
      ! do the normal quadratic reconstruction
      cm = (1.0_r_tran - tau) * (-dep + dep**2)
      cc = tau + (1.0_r_tran - tau) * (3.0_r_tran*dep - 2.0_r_tran*dep**2)
      cp = (1.0_r_tran - tau) * (1.0_r_tran - 2.0_r_tran*dep + dep**2)

    ! Relaxed limiters ---------------------------------------------------------
    case (monotone_relaxed, monotone_qm_pos)
      ! This is 1 if the turning point would be in the cell, 0 otherwise
      tau = 0.5_r_tran + SIGN(0.5_r_tran,                                      &
          (2.0_r_tran*edge_left + edge_right - 3.0_r_tran*field(:,i))          &
          * (edge_left + 2.0_r_tran*edge_right - 3.0_r_tran*field(:,i))        &
      )

      ! This is 1 if we would need to revert to a constant reconstruction
      sigma = 0.5_r_tran + SIGN(0.5_r_tran,                                    &
          (field(:,i) - edge_right) * (field(:,i) - edge_left)                 &
      )

      ! This is 1 if field_cell is closer to field_edge_right than
      ! field_edge_left, and 0 otherwise
      upsilon = 0.5_r_tran + SIGN(0.5_r_tran,                                  &
          ABS(field(:,i) - edge_left) - ABS(field(:,i) - edge_right)           &
      )

      ! Monotonicity switches, which turn on/off appropriate reconstructions:
      ! ts: apply "strict" monotonicity (reverting to a constant reconstruction)
      ts = tau*sigma
      ! tr: shift turning point to right
      tr = tau*(1.0_r_tran - sigma)*upsilon
      ! tl: shift turning point to left
      tl = tau*(1.0_r_tran - sigma)*(1.0_r_tran - upsilon)

      ! Calculate coefficients corresponding to facet reconstruction
      cp = (1.0_r_tran - ts)*(1.0_r_tran - tl)*(                               &
          1.0_r_tran - 2.0_r_tran*(1.0_r_tran - tr)*dep                        &
          + (1.0_r_tran - 2*tr)*dep**2                                         &
      )
      cc = ts + (1.0_r_tran - ts)*(                                            &
          3.0_r_tran*tl + 3.0_r_tran*(1.0_r_tran - 2*tl)*(1.0_r_tran - tr)*dep &
          + (3*tl + 3*tr - 2)*dep**2                                           &
      )
      cm = (1.0_r_tran - ts)*(1.0_r_tran - tr)*(                               &
          -2.0_r_tran*tl + (4*tl - 1.0_r_tran)*dep                             &
          + (1.0_r_tran - 2*tl)*dep**2                                         &
      )
    end select

    ! Apply weights to field and field edge values
    reconstruction(:) = cm(:)*edge_left(:) + cc(:)*field(:,i) + cp(:)*edge_right(:)

  end subroutine subgrid_quadratic_recon

end module subgrid_common_support_mod
