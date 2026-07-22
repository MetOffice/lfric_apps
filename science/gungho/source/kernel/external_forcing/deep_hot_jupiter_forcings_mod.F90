!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains forcing terms for use in the deep hot Jupiter kernels.
!>
!> @details Support functions for a kernel that adds the HD209458b test
!!          based on Menou & Rauscher (2009),
!!          Atmospheric Circulation of Hot Jupiters: A Shallow Three-Dimensional Model,
!!          ApJ, 700, 887-897, 2009, DOI: 10.1088/0004-637X/700/1/887.
!!          Also performed in Mayne et al., (2014),
!!          Using the UM dynamical cores to reproduce idealised 3-D flows,
!!          Geoscientific Model Development, Volume 7, Issue 6, 2014, pp. 3059-3087,
!!          DOI: 10.5194/gmd-7-3059-2014.

module deep_hot_jupiter_forcings_mod

  use constants_mod,     only: i_def, r_def, pi
  use planet_config_mod, only: omega, p_zero, kappa

  implicit none

  private

  public :: deep_hot_jupiter_newton_frequency
  public :: deep_hot_jupiter_equilibrium_theta

contains

!> @brief Function to calculate equilibrium theta profile for deep hot Jupiter temperature forcing.
!> @param[in] exner         Exner pressure
!> @param[in] lat           Latitude
!> @param[in] lon           Longitude
!> @return    theta_eq      Equilibrium theta
function deep_hot_jupiter_equilibrium_theta(exner, lat, lon, height) result(theta_eq)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in)    :: exner, lat, lon, height

  ! Local variables
  real(kind=r_def)                :: t_night, t_day
  real(kind=r_def)                :: t_eq, theta_eq ! Equilibrium temperature and theta theta

  ! Calculate night and day side temperature model variables
  t_night = night_side_temp(exner, height)
  t_day = day_side_temp(exner, height)

  ! Calculate equilibrium temperature according to equation 29 in Mayne et. al. 2014
  if (lon >= -1.0_r_def*pi/2.0_r_def .and. lon <= pi/2.0_r_def) then
    t_eq = (t_night**4_i_def + (t_day**2_i_def + t_night**2_i_def) * (t_day + t_night) * (t_day - t_night) &
                      * cos(lon) * cos(lat))**0.25_r_def
  else
    t_eq = t_night
  end if

  ! Recall using potential temperature
  ! Therefore, must convert the temperature to potential
  ! temperature using the Exner function (exner*theta=Temp)
  theta_eq = t_eq / exner

end function deep_hot_jupiter_equilibrium_theta

!> @brief Function to calculate the Newton relaxation frequency for deep hot Jupiter idealised test case.
!> @return    jupiter_like_frequency      Newton cooling relaxation frequency
function deep_hot_jupiter_newton_frequency(height) result(jupiter_like_frequency)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in) :: height

  ! Parameters
  real(kind=r_def), parameter :: h_low = 2833333.3333333246_r_def
  real(kind=r_def), parameter :: h_high = 10166666.666666657_r_def

  ! Local variables
  real(kind=r_def) :: log_height
  real(kind=r_def) :: jupiter_like_frequency

  ! Calculate the relaxation frequency (inverse of Newton radiative cooling timescale)
  if (height <= h_low) then
    jupiter_like_frequency = 0.0_r_def
  else
    if (height >= h_high) then
      log_height = log10(h_high / h_high)
    else
      log_height = log10(height / h_high)
    end if
    jupiter_like_frequency = 1.0_r_def /                                       &
                               (10.0_r_def **                                  &
                                 ( 3.64396238_r_def                            &
                                  -6.62330063_r_def * log_height               &
                                  -15.19454913_r_def * log_height ** 2_i_def   &
                                  -58.13043383_r_def * log_height ** 3_i_def   &
                                  -50.34836486_r_def * log_height ** 4_i_def   &
                                 )                                             &
                               )
  end if

end function deep_hot_jupiter_newton_frequency

!> @brief Function to calculate the day side temperature model variable
!> @param[in] exner         Exner pressure
!> @return    t_day         Day side temperature model variable
function day_side_temp(exner, height) result(t_day)

   implicit none

   ! Arguments
   real(kind=r_def), intent(in) :: exner, height

   ! Parameters
   real(kind=r_def), parameter :: p_low = 100.0_r_def, p_high = 1.0e6_r_def
   real(kind=r_def), parameter :: t_low = 1000.0_r_def, alpha = 0.015_r_def
   real(kind=r_def), parameter :: beta = -120.0_r_def
   real(kind=r_def), parameter :: h_low = 2833333.3333333246_r_def
   real(kind=r_def), parameter :: h_high = 10166666.666666657_r_def

   ! Local variables
   real(kind=r_def) :: pressure
   real(kind=r_def) :: t_day_active, t_day, log_height

   pressure = pressure_from_exner(exner)

   if (height <= h_low) then
     log_height = log10(h_low / h_high)
   else if (height >= h_high) then
     log_height = log10(h_high / h_high)
   else
     log_height = log10(height / h_high)
   end if

   t_day_active = 1.47255194e+03_r_def                            &
                  - 2.10792942e+03_r_def * log_height             &
                  - 2.14940583e+04_r_def * log_height ** 2_i_def  &
                  - 5.70631448e+05_r_def * log_height ** 3_i_def  &
                  - 4.27527770e+06_r_def * log_height ** 4_i_def  &
                  - 1.37981342e+07_r_def * log_height ** 5_i_def  &
                  - 1.55828322e+07_r_def * log_height ** 6_i_def  &
                  + 2.55667194e+07_r_def * log_height ** 7_i_def  &
                  + 1.02354499e+08_r_def * log_height ** 8_i_def  &
                  + 1.22122184e+08_r_def * log_height ** 9_i_def  &
                  + 5.28052053e+07_r_def * log_height ** 10_i_def

   if (pressure >= p_high) then
     t_day = t_day_active + beta * (1.0_r_def - exp(log10(p_high / pressure)))
   else if (pressure < p_low) then
     t_day = max(t_day_active * exp(alpha * log10(pressure / p_low)), t_low)
   else
     t_day = t_day_active
   end if

end function day_side_temp

!> @brief Function to calculate the night side temperature model variable
!> @param[in] exner         Exner pressure
!> @return    t_night       Night side temperature model variable
function night_side_temp(exner, height) result(t_night)

   implicit none

   ! Arguments
   real(kind=r_def), intent(in) :: exner, height

   ! Parameters
   real(kind=r_def), parameter :: p_low = 100.0_r_def, p_high = 1.0e6_r_def
   real(kind=r_def), parameter :: t_low = 250.0_r_def, alpha = 0.10_r_def
   real(kind=r_def), parameter :: beta = 100.0_r_def
   real(kind=r_def), parameter :: h_low = 2833333.3333333246_r_def
   real(kind=r_def), parameter :: h_high = 10166666.666666657_r_def

   ! Local variables
   real(kind=r_def) :: pressure
   real(kind=r_def) :: t_night_active, t_night, log_height

   pressure = pressure_from_exner(exner)

   if (height <= h_low) then
     log_height = log10(h_low / h_high)
   else if (height >= h_high) then
     log_height = log10(h_high / h_high)
   else
     log_height = log10(height / h_high)
   end if

   t_night_active = 4.75007738e+02_r_def                             &
                    - 2.33169450e+03_r_def * log_height              &
                    - 2.15104001e+04_r_def * log_height ** 2_i_def   &
                    - 3.12941407e+05_r_def * log_height ** 3_i_def   &
                    - 2.19075899e+04_r_def * log_height ** 4_i_def   &
                    + 1.54282163e+07_r_def * log_height ** 5_i_def   &
                    + 9.65080096e+07_r_def * log_height ** 6_i_def   &
                    + 2.85808438e+08_r_def * log_height ** 7_i_def   &
                    + 4.67082601e+08_r_def * log_height ** 8_i_def   &
                    + 4.06668435e+08_r_def * log_height ** 9_i_def   &
                    + 1.47801814e+08_r_def * log_height ** 10_i_def

   if (pressure >= p_high) then
     t_night = t_night_active + beta * (1.0_r_def - exp(log10(p_high / pressure)))
   else if (pressure < p_low) then
     t_night = max(t_night_active * exp(alpha * log10(pressure / p_low)), t_low)
   else
     t_night = t_night_active
   end if

end function night_side_temp

!> @brief Function to calculate pressure from exner function
!> @param[in] exner         Exner pressure
!> @return    pressure      Pressure
function pressure_from_exner(exner) result(pressure)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in) :: exner

  ! Local variables
  real(kind=r_def) :: pressure

  pressure = p_zero * exner ** (1.0_r_def / kappa)

end function pressure_from_exner

end module deep_hot_jupiter_forcings_mod
