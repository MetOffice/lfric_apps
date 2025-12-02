!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------

!> @brief Contains forcing terms for use in the held-suarez kernels

!> @details Kernel adds a Held-Suarez forcing based on Wedi and Smolarkiewicz 2009:
!> Wedi, N. P. and Smolarkiewicz, P. K. (2009), A framework for testing global
!> non-hydrostatic models. Q.J.R. Meteorol. Soc., 135: 469-484. doi: 10.1002/qj.377

module held_suarez_forcings_mod

use constants_mod,            only: r_def
use planet_config_mod,        only: scaling_factor, kappa
implicit none

private

!-------------------------------------------------------------------------------
! local parameters
!-------------------------------------------------------------------------------
! Held-Suarez parameters
real(kind=r_def), parameter :: SIGMA_B = 0.7_r_def  ! non-dimensional pressure threshold
! relaxation and damping coefficients
real(kind=r_def), parameter :: KF = 1._r_def/86400._r_def ! 1 day-1
real(kind=r_def), parameter :: KA = KF/40.0_r_def         ! 1/40 day-1
real(kind=r_def), parameter :: KS = KF/4.0_r_def          ! 1/4 day-1

real(kind=r_def), parameter :: T_MIN            = 200.0_r_def ! Minimum/Stratospheric temperature
real(kind=r_def), parameter :: T_SURF           = 315.0_r_def ! surface temperature
real(kind=r_def), parameter :: DT_EQ_POLE       = 60.0_r_def  ! Equator-Pole Temp diff (deltaT)_y
real(kind=r_def), parameter :: STATIC_STABILITY = 10.0_r_def  ! Static Stability temperature (delta \theta)_z

public held_suarez_newton_frequency
public held_suarez_damping
public held_suarez_equilibrium_theta

contains

!> @brief Function to calculate equilibrium theta profile for held-suarez
!! @param[in] exner exner pressure
!! @param[in] lat latitude
function held_suarez_equilibrium_theta(exner, lat) result(theta_eq)

  implicit none

  real(kind=r_def), intent(in)  :: exner, lat
  real(kind=r_def) :: theta_eq         ! Equilibrium theta

  theta_eq = max(T_MIN/exner, (T_SURF - DT_EQ_POLE*sin(lat)*sin(lat) &
     - STATIC_STABILITY*log(exner)*cos(lat)*cos(lat)/kappa))

end function held_suarez_equilibrium_theta

!> @brief Function to calculate the newton relaxation frequency for held-suarez
!! @param[in] sigma nondimensional pressure p/p_surf
!! @param[in] lat latitude
function held_suarez_newton_frequency(sigma, lat) result(held_suarez_frequency)

  implicit none

  real(kind=r_def), intent(in) :: sigma
  real(kind=r_def), intent(in) :: lat
  real(kind=r_def)             :: held_suarez_frequency
  real(kind=r_def)             :: sigma_func

  sigma_func = max((sigma - SIGMA_B)/(1.0_r_def - SIGMA_B), 0.0_r_def)
  held_suarez_frequency = KA + (KS - KA)*sigma_func*(cos(lat)**4)

  ! If running on a scaled planet, then reduce the timescale...
  held_suarez_frequency = held_suarez_frequency*scaling_factor

end function held_suarez_newton_frequency

!> @brief Function to calculate the damping coefficent for held-suarez
!! @param[in] sigma nondimensional pressure p/p_surf
function held_suarez_damping(sigma) result(held_suarez_damping_rate)

  implicit none

  real(kind=r_def), intent(in) :: sigma
  real(kind=r_def)             :: held_suarez_damping_rate
  real(kind=r_def) :: sigma_func

  sigma_func = max((sigma - SIGMA_B)/(1.0_r_def - SIGMA_B), 0.0_r_def)
  held_suarez_damping_rate = -KF*sigma_func

  ! If running on a scaled planet, then reduce the timescale...
  held_suarez_damping_rate = held_suarez_damping_rate*scaling_factor

end function held_suarez_damping

end module held_suarez_forcings_mod
