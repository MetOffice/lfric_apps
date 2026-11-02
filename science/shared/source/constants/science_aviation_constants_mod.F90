! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
!> @brief LFRic aviation constants module
!>
module science_aviation_constants_mod

  use constants_mod, only: r_def
  use science_conversions_mod, only: feet_to_metres

  implicit none

  private
  public :: mtokft, isa_lapse_ratel, isa_lapse_rateu, &
      isa_press_bot, isa_press_mid, isa_press_top, &
      isa_temp_bot, isa_temp_top, gpm1, gpm2

  ! Meters to thousands of feet.
  real(r_def), parameter :: mtokft = 1.0_r_def / (feet_to_metres * 1000.0_r_def)

  ! International Standard Atmosphere (ISA) within troposphere for:
  ! Values taken from:
  ! https://en.wikipedia.org/wiki/international_standard_atmosphere
  ! Lapse rate (degreeC/km) for levels below 11,000 gpm.
  real(r_def), parameter :: isa_lapse_ratel = 6.5e-03_r_def

  ! Lapse rate (degreeC/km) for levels above 20,000 gpm.
  real(r_def), parameter :: isa_lapse_rateu = -1.0e-03_r_def

  ! Surface pressure (Pa).
  real(r_def), parameter :: isa_press_bot = 101325.0_r_def

  ! Pressure (Pa) at 11,000 gpm.
  real(r_def), parameter :: isa_press_mid = 22632.0_r_def

  ! Pressure (Pa) at 20,000 gpm.
  real(r_def), parameter :: isa_press_top = 5474.87_r_def

  ! Surface temperature (K).
  real(r_def), parameter :: isa_temp_bot = 288.15_r_def

  ! Temperature (K) of isothermal layer - at tropopause.
  real(r_def), parameter :: isa_temp_top = 216.65_r_def

  ! Height limit (gpm) for standard lower lapse rate.
  real(r_def), parameter :: gpm1 = 11000.0_r_def

  ! Height (gpm) of top of isothermal layer.
  real(r_def), parameter :: gpm2 = 20000.0_r_def

end module science_aviation_constants_mod

