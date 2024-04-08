!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains forcing terms for the Tidally Locked Earth (TLE) case.
!>
!> @details Contains functions to calculate a reference profile and timescale
!!          that are then used by the corresponding TLE kernel, which adds
!!          increments to the temperature field (similar to the Held-Suarez
!!          test).
!!          References:
!!          * Mayne, N. J., Baraffe, I., Acreman, D. M., Smith, C., Wood, N.,
!!          Amundsen, D. S., Thuburn, J., and Jackson, D. R.: Using the UM
!!          dynamical cores to reproduce idealised 3-D flows, Geosci. Model
!!          Dev., 7, 3059-3087, https://doi.org/10.5194/gmd-7-3059-2014, 2014.
!!          * Merlis, T. M., and Schneider, T. (2010), Atmospheric Dynamics of
!!          Earth-Like Tidally Locked Aquaplanets, J. Adv. Model. Earth Syst.,
!!          2, 13, doi:10.3894/JAMES.2010.2.13.

module tidally_locked_earth_forcings_mod

  use constants_mod,     only: r_def, i_def, pi
  use planet_config_mod, only: kappa

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Local parameters
  !---------------------------------------------------------------------------
  ! Stratospheric temperature
  real(kind=r_def), parameter :: T_STRA = 200.0_r_def
  ! Surface temperature at equator
  real(kind=r_def), parameter :: T_SURF = 315.0_r_def
  ! Equator-to-pole temperature difference
  real(kind=r_def), parameter :: DT_EQ_POLE = 60.0_r_def
  ! Static stability
  real(kind=r_def), parameter :: DT_Z = 10.0_r_def

  !---------------------------------------------------------------------------
  ! Local variables
  !---------------------------------------------------------------------------
  real(kind=r_def) :: theta1, theta2 ! Potential temperatures

  public :: tidally_locked_earth_equilibrium_theta

contains

!> @brief Calculates equilibrium theta profile for the TLE case.
!>
!> @details Calculates the equilibrium (reference) profile of air potential
!!          temperature using eqns. 17 and 18 in Mayne et al. (2014).
!> @param[in] exner         Exner function
!> @param[in] lon           Longitude
!> @param[in] lat           Latitude
!> @return    theta_eq      Equilibrium theta profile
function tidally_locked_earth_equilibrium_theta(exner, lon, lat) &
  result(theta_eq)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in) :: exner    ! Exner function
  real(kind=r_def), intent(in) :: lon, lat ! Longitude and latitude
  real(kind=r_def)             :: theta_eq ! Equilibrium theta

  ! We are operating in terms of potential temperature, which is related to
  ! the real temperature via Exner function:
  ! theta * exner = temp
  theta1 = T_STRA / exner
  ! and log(exner) = kappa * log(p/p0)
  theta2 = T_SURF + DT_EQ_POLE * cos(lon) * cos(lat) &
           - DT_Z * (log(exner) / kappa) * cos(lat) * cos(lat)
  ! Note that the substellar point (maximum forcing) is at lon=0,
  ! and not at lon=180 (as it was in the UM).
  theta_eq = max(theta1, theta2)

end function tidally_locked_earth_equilibrium_theta

end module tidally_locked_earth_forcings_mod
