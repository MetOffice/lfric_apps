!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Routines for calculating semi-Lagrangian interpolation coefficients.
!!
!> @details This module contains functions and subroutines which allow the
!!          semi-Lagrangian transport scheme to pre-compute interpolation
!!          coefficients, and to apply monotonicity to the transported field.
!------------------------------------------------------------------------------

module sl_support_mod

  use constants_mod,                   only: i_def, r_tran
  use transport_enumerated_types_mod,  only: monotone_none,                    &
                                             monotone_strict,                  &
                                             monotone_relaxed,                 &
                                             vertical_monotone_order_constant, &
                                             vertical_monotone_order_linear,   &
                                             vertical_monotone_order_high

  implicit none

  private

  public :: monotone_cubic_sl
  public :: monotone_quintic_sl
  public :: compute_linear_coeffs
  public :: compute_cubic_coeffs
  public :: compute_quintic_coeffs
  public :: compute_cubic_hermite_coeffs

  contains

  !-----------------------------------------------------------------------------
  !> @brief   Limits the cubic semi-Lagrangian advective transport.
  !> @details Identifies the non-monotone cubic-interpolated values
  !!          and modifies them using the specified order/option.
  !> @param[in,out] field_dep       The calculated field at departure points
  !> @param[in]     field_local     Grid-data field around the departure points
  !> @param[in]     linear_coef_1   First linear coefficient
  !> @param[in]     linear_coef_2   Second linear coefficient
  !> @param[in]     vertical_monotone
  !!                                Option to identify non-monotone points
  !> @param[in]     vertical_monotone_order
  !!                                Specifies the monotonic modification
  !> @param[in]     nl              Number of data points in the vertical
  !-----------------------------------------------------------------------------
  subroutine monotone_cubic_sl( field_dep, field_local,       &
                                linear_coef_1, linear_coef_2, &
                                vertical_monotone,            &
                                vertical_monotone_order,      &
                                nl )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nl
    integer(kind=i_def), intent(in)    :: vertical_monotone
    integer(kind=i_def), intent(in)    :: vertical_monotone_order
    real(kind=r_tran),   intent(inout) :: field_dep(nl)
    real(kind=r_tran),   intent(in)    :: field_local(nl,4)
    real(kind=r_tran),   intent(in)    :: linear_coef_1(nl)
    real(kind=r_tran),   intent(in)    :: linear_coef_2(nl)

    ! Local variables
    real(kind=r_tran)   :: tau(nl), sigma(nl)

    ! Determine if monotonicity needs applying to each point -------------------
    select case (vertical_monotone)
    case (monotone_strict)
      ! 1 if there is a new extremum between the central points, 0 otherwise
      tau = 0.5_r_tran + SIGN(0.5_r_tran,                                      &
        (field_local(:,2) - field_dep(:)) * (field_local(:,3) - field_dep(:))  &
      )

    case (monotone_relaxed)
      ! This switch is equivalent to the following logical process:
      ! is_new_extremum AND NOT (expect_new_extremum AND correct_extremum)
      tau = (0.5_r_tran + SIGN(0.5_r_tran,                                     &
        ! is_new_extremum: 1 if new extremum between central points, 0 otherwise
        (field_local(:,2) - field_dep(:)) * (field_local(:,3) - field_dep(:))  &
      )) * (1.0_r_tran - (0.5_r_tran + SIGN(0.5_r_tran,                        &
        ! expect_new_extremum: 1 if cell should be extremum, 0 otherwise
        (field_local(:,2) - field_local(:,1))                                  &
        * (field_local(:,3) - field_local(:,4))                                &
      )) * (0.5_r_tran + SIGN(0.5_r_tran,                                      &
        ! correct_extremum: 1 if extremum is in correct direction, 0 otherwise
        ! NB: this switch is only relevant if expecting new extremum. Don't need
        ! to also test with field_local(:,3) and field_local(:,4) because in
        ! that situation we would get the same answer
        (field_dep(:) - field_local(:,2))                                      &
        * (field_local(:,2) - field_local(:,1))                                &
      )))
    end select

    ! Apply monotonicity -------------------------------------------------------
    if (vertical_monotone /= monotone_none) then
      select case (vertical_monotone_order)
      case (vertical_monotone_order_constant)
        ! if monotonicity needs applying, just bound the field by its neighbours
        field_dep(:) = (1.0_r_tran - tau) * field_dep + tau * MIN(             &
          MAX(field_local(:,2), field_local(:,3)),                             &
          MAX(field_dep(:), MIN(field_local(:,2), field_local(:,3)))           &
        )

      case (vertical_monotone_order_linear)
        ! if monotonicity needs applying, revert to linear reconstruction
        field_dep(:) = (1.0_r_tran - tau)*field_dep(:) + tau * (               &
          linear_coef_1(:)*field_local(:,2)                                    &
          + linear_coef_2(:)*field_local(:,3)                                  &
        )

      case (vertical_monotone_order_high)
        ! sigma is 1 if the dep field value is closer to right field value than
        ! the left, and 0 if closer to left field value than the right
        sigma = 0.5_r_tran + SIGN(0.5_r_tran,                                  &
          (field_dep(:) - field_local(:,2))                                    &
          * (field_local(:,3) - field_local(:,2))                              &
        )
        ! Interpolated field combines sigma and tau switches:
        ! tau switches on monotonicity, sigma chooses if shifting left or right
        field_dep(:) = (1.0_r_tran - tau)*field_dep(:) + tau * (               &
        ! Monotonicity uses a bounded cubic polynomial
          (1.0_r_tran - sigma) * (field_local(:,2)                             &
            + (field_local(:,3) - field_local(:,2))*linear_coef_2(:)**3        &
          )                                                                    &
          + sigma * (field_local(:,3)                                          &
            - (field_local(:,3) - field_local(:,2))*linear_coef_1(:)**3        &
          )                                                                    &
        )
      end select
    end if

  end subroutine monotone_cubic_sl

  !-----------------------------------------------------------------------------
  !> @brief   Limits the quintic semi-Lagrangian advective transport.
  !> @details Identifies the non-monotone quintic-interpolated values
  !!          and modifies them using the specified order/option.
  !> @param[in,out] field_dep       The calculated field at departure points
  !> @param[in]     field_local     Grid-data field around the departure points
  !> @param[in]     linear_coef_1   First linear coefficient
  !> @param[in]     linear_coef_2   Second linear coefficient
  !> @param[in]     vertical_monotone
  !!                                Option to identify non-monotone points
  !> @param[in]     vertical_monotone_order
  !!                                Specifies the monotonic modification
  !> @param[in]     nl              Number of data points in the vertical
  !-----------------------------------------------------------------------------
  subroutine monotone_quintic_sl( field_dep, field_local,       &
                                  linear_coef_1, linear_coef_2, &
                                  vertical_monotone,            &
                                  vertical_monotone_order,      &
                                  nl )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nl
    integer(kind=i_def), intent(in)    :: vertical_monotone
    integer(kind=i_def), intent(in)    :: vertical_monotone_order
    real(kind=r_tran),   intent(inout) :: field_dep(nl)
    real(kind=r_tran),   intent(in)    :: field_local(6,nl)
    real(kind=r_tran),   intent(in)    :: linear_coef_1(nl)
    real(kind=r_tran),   intent(in)    :: linear_coef_2(nl)

    ! Local variables
    real(kind=r_tran)   :: tau(nl), sigma(nl)

    ! Determine if monotonicity needs applying to each point -------------------
    select case (vertical_monotone)
    case (monotone_strict)
      ! 1 if there is a new extremum between the central points, 0 otherwise
      tau = 0.5_r_tran + SIGN(0.5_r_tran,                                      &
        (field_local(:,3) - field_dep(:)) * (field_local(:,4) - field_dep(:))  &
      )

    case (monotone_relaxed)
      ! This switch is equivalent to the following logical process:
      ! is_new_extremum AND NOT (expect_new_extremum AND correct_extremum)
      tau = (0.5_r_tran + SIGN(0.5_r_tran,                                     &
        ! is_new_extremum: 1 if new extremum between central points, 0 otherwise
        (field_local(:,3) - field_dep(:)) * (field_local(:,4) - field_dep(:))  &
      )) * (1.0_r_tran - (0.5_r_tran + SIGN(0.5_r_tran,                        &
        ! expect_new_extremum: 1 if cell should be extremum, 0 otherwise
        (field_local(:,3) - field_local(:,2))                                  &
        * (field_local(:,4) - field_local(:,5))                                &
      )) * (0.5_r_tran + SIGN(0.5_r_tran,                                      &
        ! correct_extremum: 1 if extremum is in correct direction, 0 otherwise
        ! NB: this switch is only relevant if expecting new extremum. Don't need
        ! to also test with field_local(:,3) and field_local(:,4) because in
        ! that situation we would get the same answer
        (field_dep(:) - field_local(:,3))                                      &
        * (field_local(:,3) - field_local(:,2))                                &
      )))
    end select

    ! Apply monotonicity -------------------------------------------------------
    if (vertical_monotone /= monotone_none) then
      select case (vertical_monotone_order)
      case (vertical_monotone_order_constant)
        ! if monotonicity needs applying, just bound the field by its neighbours
        field_dep(:) = (1.0_r_tran - tau) * field_dep + tau * MIN(             &
          MAX(field_local(:,3), field_local(:,4)),                             &
          MAX(field_dep(:), MIN(field_local(:,3), field_local(:,4)))           &
        )

      case (vertical_monotone_order_linear)
        ! if monotonicity needs applying, revert to linear reconstruction
        field_dep(:) = (1.0_r_tran - tau)*field_dep(:) + tau * (               &
          linear_coef_1(:)*field_local(:,3)                                    &
          + linear_coef_2(:)*field_local(:,4)                                  &
        )

      case (vertical_monotone_order_high)
        ! sigma is 1 if the dep field value is closer to right field value than
        ! the left, and 0 if closer to left field value than the right
        sigma = 0.5_r_tran + SIGN(0.5_r_tran,                                  &
          (field_dep(:) - field_local(:,3))                                    &
          * (field_local(:,4) - field_local(:,3))                              &
        )
        ! Interpolated field combines sigma and tau switches:
        ! tau switches on monotonicity, sigma chooses if shifting left or right
        field_dep(:) = (1.0_r_tran - tau)*field_dep(:) + tau * (               &
        ! Monotonicity uses a bounded quintic polynomial
          (1.0_r_tran - sigma) * (field_local(:,3)                             &
            + (field_local(:,4) - field_local(:,3))*linear_coef_2(:)**5        &
          )                                                                    &
          + sigma * (field_local(:,4)                                          &
            - (field_local(:,4) - field_local(:,3))*linear_coef_1(:)**5        &
          )                                                                    &
        )
      end select
    end if

  end subroutine monotone_quintic_sl

  !-----------------------------------------------------------------------------
  !> @brief   This subroutine computes the cubic-Lagrange weights.
  !> @details Compute the cubic-Lagrange interpolation weights for vertical
  !!          semi-Lagrangian advective transport.
  !> @param[in]  k_dep   Indices of the grid points below the departure points
  !> @param[in]  z_dep   The departure (interpolation) points
  !> @param[in]  z_arr   The grid points where the data is located
  !> @param[out] indices Column index for interpolation for each point
  !> @param[out] coeffs  The cubic interpolation weights
  !> @param[in]  nl      Number of values in a column
  !-----------------------------------------------------------------------------
  subroutine compute_cubic_coeffs(k_dep, z_dep, z_arr, dz, indices, coeffs, nl)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nl
    integer(kind=i_def), intent(in)  :: k_dep(nl)
    real(kind=r_tran),   intent(in)  :: z_dep(nl)
    real(kind=r_tran),   intent(in)  :: z_arr(nl)
    real(kind=r_tran),   intent(in)  :: dz(nl)
    real(kind=r_tran),   intent(out) :: coeffs(nl,4)
    integer(kind=i_def), intent(out) :: indices(nl,4)

    ! Local variables
    real(kind=r_tran)   :: z1, z2, z3, z4, xi
    real(kind=r_tran)   :: d1, d2, d3, d4
    real(kind=r_tran)   :: n1, n2, n3, n4
    integer(kind=i_def) :: k, km

    do k = 1, nl
      km = k_dep(k)
      xi = (z_dep(k) - z_arr(km))/dz(km)

      indices(k,1) = max(1, km - 1)
      indices(k,2) = km
      indices(k,3) = min(nl, km + 1)
      indices(k,4) = min(nl, km + 2)

      z1 = 0.0_r_tran
      z2 = z1 + dz(indices(k,1))
      z3 = z2 + dz(indices(k,2))
      z4 = z3 + dz(indices(k,3))
      xi = z2 + xi*(z3-z2)

      d1 = (z1-z2)*(z1-z3)*(z1-z4)
      d2 = (z2-z1)*(z2-z3)*(z2-z4)
      d3 = (z3-z1)*(z3-z2)*(z3-z4)
      d4 = (z4-z1)*(z4-z2)*(z4-z3)

      n1 = (xi-z2)*(xi-z3)*(xi-z4)
      n2 = (xi-z1)*(xi-z3)*(xi-z4)
      n3 = (xi-z1)*(xi-z2)*(xi-z4)
      n4 = (xi-z1)*(xi-z2)*(xi-z3)

      ! cubic weights
      coeffs(k,1) = n1/d1
      coeffs(k,2) = n2/d2
      coeffs(k,3) = n3/d3
      coeffs(k,4) = n4/d4

      ! Next to boundaries there are not enough points for cubic
      ! so revert to linear

      if (indices(k,1) == indices(k,2) .or. indices(k,3) == indices(k,4)) then
        coeffs(k,1) = 0.0_r_tran
        coeffs(k,2) = (z3-xi)/(z3-z2)
        coeffs(k,3) = 1.0_r_tran - coeffs(k,2)
        coeffs(k,4) = 0.0_r_tran
      end if
    end do

  end subroutine compute_cubic_coeffs

  !-----------------------------------------------------------------------------
  !> @brief   This subroutine computes cubic-Hermite weights.
  !> @details Compute the cubic-Hermite interpolation weights for vertical
  !!          semi-Lagrangian advective transport.
  !> @param[in]  k_dep   Indices of the grid points below the departure points
  !> @param[in]  z_dep   The departure (interpolation) points
  !> @param[in]  z_arr   The grid points where the data is located
  !> @param[out] indices Column index for interpolation for each point
  !> @param[out] coeffs  The cubic interpolation weights
  !> @param[in]  nl      Number of values in a column
  !-----------------------------------------------------------------------------
  subroutine compute_cubic_hermite_coeffs(k_dep, z_dep, z_arr, dz, indices, coeffs, nl)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nl
    integer(kind=i_def), intent(in)  :: k_dep(nl)
    real(kind=r_tran),   intent(in)  :: z_dep(nl)
    real(kind=r_tran),   intent(in)  :: z_arr(nl)
    real(kind=r_tran),   intent(in)  :: dz(nl)
    real(kind=r_tran),   intent(out) :: coeffs(nl,4)
    integer(kind=i_def), intent(out) :: indices(nl,4)

    ! Local variables
    real(kind=r_tran)   :: xi, alpha, beta, inv_1p_alpha, inv_1p_beta
    real(kind=r_tran)   :: xip2, xip3, c1, c2, c3, c4
    integer(kind=i_def) :: k, km

    do k = 1, nl
      km = k_dep(k)
      xi = (z_dep(k) - z_arr(km))/dz(km)

      indices(k,1) = max(1, km - 1)
      indices(k,2) = km
      indices(k,3) = min(nl, km + 1)
      indices(k,4) = min(nl, km + 2)

      alpha = dz(indices(k,1))/dz(indices(k,2))
      beta = dz(indices(k,3))/dz(indices(k,2))
      inv_1p_alpha = 1.0_r_tran/(1.0_r_tran + alpha)
      inv_1p_beta = 1.0_r_tran/(1.0_r_tran + beta)
      xip2 = xi**2
      xip3 = xi**3
      c1 = 2.0_r_tran*xip3 - 3.0_r_tran*xip2 + 1.0_r_tran
      c2 = 1.0_r_tran - c1
      c3 = xip3 - 2.0_r_tran*xip2 + xi
      c4 = xip3 - xip2

      ! Cubic-hermite weights
      if (indices(k,1) == indices(k,2)) then
        coeffs(k,1) = 0.0_r_tran
        coeffs(k,2) = c1 - c3 - c4*inv_1p_beta
        coeffs(k,3) = c2 + c3
        coeffs(k,4) = c4*inv_1p_beta
      else if (indices(k,3) == indices(k,4)) then
        coeffs(k,1) = -c3*inv_1p_alpha
        coeffs(k,2) = c1 - c4
        coeffs(k,3) = c2 + c4 + c3*inv_1p_alpha
        coeffs(k,4) = 0.0_r_tran
      else
        coeffs(k,1) = -c3*inv_1p_alpha
        coeffs(k,2) = c1 - c4*inv_1p_beta
        coeffs(k,3) = c2 + c3*inv_1p_alpha
        coeffs(k,4) = c4*inv_1p_beta
      end if

    end do

  end subroutine compute_cubic_hermite_coeffs

  !-----------------------------------------------------------------------------
  !> @brief   This subroutine computes cubic-Hermite weights.
  !> @details Compute the cubic-Hermite interpolation weights for vertical
  !!          semi-Lagrangian advective transport.
  !> @param[in]  k_dep   Indices of the grid points below the departure points
  !> @param[in]  z_dep   The departure (interpolation) points
  !> @param[in]  z_arr   The grid points where the data is located
  !> @param[out] coeffs  The linear interpolation weights
  !> @param[in]  nl      Number of values in a column
  !-----------------------------------------------------------------------------
  subroutine compute_linear_coeffs(k_dep, z_dep, z_arr, dz, coeffs, nl)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nl
    integer(kind=i_def), intent(in)  :: k_dep(nl)
    real(kind=r_tran),   intent(in)  :: z_dep(nl)
    real(kind=r_tran),   intent(in)  :: z_arr(nl)
    real(kind=r_tran),   intent(in)  :: dz(nl)
    real(kind=r_tran),   intent(out) :: coeffs(nl,2)

    ! Local variables
    real(kind=r_tran)   :: xi
    integer(kind=i_def) :: k, km

    do k = 1, nl
      km = k_dep(k)
      xi = (z_dep(k) - z_arr(km))/dz(km)

      ! linear weights
      coeffs(k,2) = xi
      coeffs(k,1) = 1.0_r_tran - coeffs(k,2)
    end do

  end subroutine compute_linear_coeffs

  !-------------------------------------------------------------------------------
  !> @details This subroutine computes the quintic-Lagrange weights.
  !> @details Compute the quintic-Lagrange interpolation weights for vertical
  !!          semi-Lagrangian advective transport.
  !> @param[in]  k_dep   Indices of the grid points below the departure points
  !> @param[in]  z_dep   The departure (interpolation) points
  !> @param[in]  z_arr   The grid points where the data is located
  !> @param[out] indices Column index for interpolation for each point
  !> @param[out] coeffs  The cubic interpolation weights
  !> @param[in]  nl      Number of values in a column
  !-------------------------------------------------------------------------------
  subroutine compute_quintic_coeffs(k_dep, z_dep, z_arr, dz, indices, coeffs, nl)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nl
    integer(kind=i_def), intent(in)  :: k_dep(nl)
    real(kind=r_tran),   intent(in)  :: z_dep(nl)
    real(kind=r_tran),   intent(in)  :: z_arr(nl)
    real(kind=r_tran),   intent(in)  :: dz(nl)
    real(kind=r_tran),   intent(out) :: coeffs(nl,6)
    integer(kind=i_def), intent(out) :: indices(nl,6)

    ! Local variables
    real(kind=r_tran)   :: z1, z2, z3, z4, z5, z6, xi
    real(kind=r_tran)   :: d1, d2, d3, d4, d5, d6
    real(kind=r_tran)   :: n1, n2, n3, n4, n5, n6
    integer(kind=i_def) :: k, km

    do k = 1, nl
      km = k_dep(k)
      xi = (z_dep(k) - z_arr(km))/dz(km)

      indices(k,1) = max(1, km - 2)
      indices(k,2) = max(1, km - 1)
      indices(k,3) = km
      indices(k,4) = min(nl, km + 1)
      indices(k,5) = min(nl, km + 2)
      indices(k,6) = min(nl, km + 3)

      z1 = 0.0_r_tran
      z2 = z1 + dz(indices(k,1))
      z3 = z2 + dz(indices(k,2))
      z4 = z3 + dz(indices(k,3))
      z5 = z4 + dz(indices(k,4))
      z6 = z5 + dz(indices(k,5))
      xi = z3 + xi*(z4-z3)

      d1 = (z1-z2)*(z1-z3)*(z1-z4)*(z1-z5)*(z1-z6)
      d2 = (z2-z1)*(z2-z3)*(z2-z4)*(z2-z5)*(z2-z6)
      d3 = (z3-z1)*(z3-z2)*(z3-z4)*(z3-z5)*(z3-z6)
      d4 = (z4-z1)*(z4-z2)*(z4-z3)*(z4-z5)*(z4-z6)
      d5 = (z5-z1)*(z5-z2)*(z5-z3)*(z5-z4)*(z5-z6)
      d6 = (z6-z1)*(z6-z2)*(z6-z3)*(z6-z4)*(z6-z5)

      n1 = (xi-z2)*(xi-z3)*(xi-z4)*(xi-z5)*(xi-z6)
      n2 = (xi-z1)*(xi-z3)*(xi-z4)*(xi-z5)*(xi-z6)
      n3 = (xi-z1)*(xi-z2)*(xi-z4)*(xi-z5)*(xi-z6)
      n4 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z5)*(xi-z6)
      n5 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z4)*(xi-z6)
      n6 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z4)*(xi-z5)

      ! quintic weights
      coeffs(k,1) = n1/d1
      coeffs(k,2) = n2/d2
      coeffs(k,3) = n3/d3
      coeffs(k,4) = n4/d4
      coeffs(k,5) = n5/d5
      coeffs(k,6) = n6/d6

      if (indices(k,1) == indices(k,2) .or. indices(k,5) == indices(k,6)) then
        ! Revert to cubic weights
        d1 = (z2-z3)*(z2-z4)*(z2-z5)
        d2 = (z3-z2)*(z3-z4)*(z3-z5)
        d3 = (z4-z2)*(z4-z3)*(z4-z5)
        d4 = (z5-z2)*(z5-z3)*(z5-z4)

        n1 = (xi-z3)*(xi-z4)*(xi-z5)
        n2 = (xi-z2)*(xi-z4)*(xi-z5)
        n3 = (xi-z2)*(xi-z3)*(xi-z5)
        n4 = (xi-z2)*(xi-z3)*(xi-z4)

        coeffs(k,1) = 0.0_r_tran
        coeffs(k,2) = n1/d1
        coeffs(k,3) = n2/d2
        coeffs(k,4) = n3/d3
        coeffs(k,5) = n4/d4
        coeffs(k,6) = 0.0_r_tran
      end if

      if (indices(k,2) == indices(k,3) .or. indices(k,4) == indices(k,5)) then
        ! Revert to linear weights
        coeffs(k,1) = 0.0_r_tran
        coeffs(k,2) = 0.0_r_tran
        coeffs(k,3) = (z4-xi)/(z4-z3)
        coeffs(k,4) = 1.0_r_tran - coeffs(k,3)
        coeffs(k,5) = 0.0_r_tran
        coeffs(k,6) = 0.0_r_tran
      end if
    end do

  end subroutine compute_quintic_coeffs

end module sl_support_mod