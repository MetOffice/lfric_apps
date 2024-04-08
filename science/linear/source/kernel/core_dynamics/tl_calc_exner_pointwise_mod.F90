!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module tl_calc_exner_pointwise_mod

use constants_mod,     only : r_def
use planet_config_mod, only : kappa, Rd, p_zero

implicit none

private

public :: tl_calc_exner_pointwise

contains
!-------------------------------------------------------------------------------
! Contained subroutines
!
! These subroutines are called by the tangent linear kernels, which are then
! passed through PSyAD to create the adjoint. At the present time, PSyAD can
! only deal with subroutines and not functions. This is because a function
! with 1 output and 2 inputs would require the adjoint function to have 2
! outputs and 1 input, which is not possible. However it is possible for a
! subroutine to have 2 outputs.
!
! In summary, the following subroutines should NOT be changed to functions,
! even though the corresponding nonlinear gungho code uses functions.
!
!-------------------------------------------------------------------------------

!> @brief Compute the change in exner pressure from the tangent
!>        linear of the equation of state.
!> @details The nonlinear equation of state is:
!>          exner = ( Rd/p0 * rho * theta ) ^ (  k / ( 1 - k ) )
!>          The tangent linear is:
!>          exner = (k/(1-k)) * ls_exner * ( rho / ls_rho + theta / ls_theta )
!! @param[in,out] exner    Change in pressure
!! @param[in]     rho      Change in density
!! @param[in]     theta    Change in potential temperature
!! @param[in]     ls_rho   Linearisation state for density
!! @param[in]     ls_theta Linearisation state for potential temperature
subroutine tl_calc_exner_pointwise(exner, rho, theta, ls_rho, ls_theta)

  implicit none

  real(kind=r_def), intent(inout) :: exner
  real(kind=r_def), intent(in)    :: rho, theta
  real(kind=r_def), intent(in)    :: ls_rho, ls_theta
  real(kind=r_def)                :: ls_exner

  ls_exner = ( ( Rd / p_zero ) * ls_rho * ls_theta ) ** &
             ( kappa / ( 1.0_r_def - kappa ) )

  exner = ( kappa / ( 1.0_r_def - kappa ) ) * ls_exner * &
          ( ( rho / ls_rho ) + ( theta / ls_theta )  )

end subroutine tl_calc_exner_pointwise

end module tl_calc_exner_pointwise_mod
