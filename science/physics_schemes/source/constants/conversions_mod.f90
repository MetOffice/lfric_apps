!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic interface module for UM code (conversions_mod)
!----------------------------------------------------------------------------

module conversions_mod

  use, intrinsic :: iso_fortran_env, only: real32
  use constants_mod, only: r_um, i_um, r_bl,   &
                           degrees_to_radians, &
                           radians_to_degrees, &
                           lfric_pi => pi

  use science_conversions_mod, only: seconds_per_day,           &
                                     seconds_per_hour,          &
                                     seconds_per_minute,        &
                                     hours_per_day,             &
                                     seconds_to_hours,          &
                                     hours_to_days,             &
                                     zero_degrees_celsius,      &
                                     knots_to_metre_per_second, &
                                     feet_to_metres

  implicit none

  private
  public :: isec_per_day, isec_per_hour, pi, pi_over_180, recip_pi_over_180,   &
            rhour_per_day, rhour_per_sec, rsec_per_day, rsec_per_hour,         &
            zerodegc, zerodegc_32, pi_bl

  ! The following variables have been hidden as they are not currently
  ! required to build the extracted UM code. They have been left in
  ! in case they are required as more UM code is drawn into the lfric_atm
  ! build. Should they be required at a later date, they should simply be
  ! added to the public statement above.

  ! Disabled variables:
  !   isec_per_min, ihour_per_day, rday_per_hour, kt2ms, ft2m


  integer(i_um), parameter :: isec_per_day  = int(seconds_per_day,    i_um)
  integer(i_um), parameter :: isec_per_hour = int(seconds_per_hour,   i_um)
  integer(i_um), parameter :: isec_per_min  = int(seconds_per_minute, i_um)
  integer(i_um), parameter :: ihour_per_day = int(hours_per_day,      i_um)

  real(r_um),    parameter :: rsec_per_day  = real(seconds_per_day,   r_um)
  real(r_um),    parameter :: rsec_per_hour = real(seconds_per_hour,  r_um)
  real(r_um),    parameter :: rhour_per_day = real(hours_per_day,     r_um)

  real(r_um),    parameter :: rhour_per_sec = real(seconds_to_hours,  r_um)
  real(r_um),    parameter :: rday_per_hour = real(hours_to_days,     r_um)

  real(r_um),    parameter :: pi = real( lfric_pi, r_um )
  real(r_bl),    parameter :: pi_bl = real( lfric_pi, r_bl )

  ! Conversion factor degrees to radians
  real(r_um),    parameter :: pi_over_180 = real(degrees_to_radians, r_um)

  ! Conversion factor radians to degrees
  real(r_um),    parameter :: recip_pi_over_180 = real(radians_to_degrees, r_um)

  ! zerodegc is a conversion between degrees centigrade and kelvin
  real(r_um),    parameter :: zerodegc     = real(zero_degrees_celsius, r_um)
  real(real32),  parameter :: zerodegc_32  = real(zero_degrees_celsius, real32)

  ! Knots to m/s conversion
  real(r_um),    parameter :: kt2ms = real(knots_to_metre_per_second, r_um)

  ! Feet to metres conversion
  real(r_um),    parameter :: ft2m  = real(feet_to_metres, r_um)

end module conversions_mod
