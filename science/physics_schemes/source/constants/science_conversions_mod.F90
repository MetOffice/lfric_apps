!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic conversions module
!----------------------------------------------------------------------------

module science_conversions_mod

  use constants_mod,     only : r_def
  use sci_constants_mod, only : zero_C_in_K

  implicit none

  private
  public :: seconds_per_day, seconds_per_hour, seconds_per_minute, &
            hours_per_day, seconds_to_hours, hours_to_days,        &
            zero_degrees_celsius, knots_to_metre_per_second,       &
            feet_to_metres

  !> @name Time conversions (Conversion factors)
  !> @{
  real(r_def), parameter :: seconds_per_day    = 86400.0_r_def
  real(r_def), parameter :: seconds_per_hour   = 3600.0_r_def
  real(r_def), parameter :: seconds_per_minute = 60.0_r_def
  real(r_def), parameter :: hours_per_day      = 24.0_r_def
  real(r_def), parameter :: seconds_to_hours   = 1.0_r_def/seconds_per_hour
  real(r_def), parameter :: hours_to_days      = 1.0_r_def/hours_per_day
  !> @}

  ! Celsius to Kelvin (Conversion offset)
  real(r_def), parameter :: zero_degrees_celsius = zero_C_in_K

  ! Knots to m/s (Conversion factor)
  real(r_def), parameter :: knots_to_metre_per_second = 1852.0_r_def / &
                                                        seconds_per_hour
  ! Feet to metres (Conversion factor)
  real(r_def), parameter :: feet_to_metres = 0.3048_r_def

end module science_conversions_mod
