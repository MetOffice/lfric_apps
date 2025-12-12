!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic aviation constants module
!----------------------------------------------------------------------------
!
MODULE science_aviation_constants_mod

  USE constants_mod, ONLY: r_def
  USE science_conversions_mod, ONLY: feet_to_metres

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mtokft, isa_Lapse_RateL, isa_Lapse_RateU, &
      isa_Press_Bot, isa_Press_Mid, isa_Press_Top, &
      isa_Temp_Bot, isa_Temp_Top, Gpm1, Gpm2

  ! meters to thousands of feet
  REAL(r_def), PARAMETER :: mtokft = 1.0_r_def / (feet_to_metres * 1000.0_r_def)

  ! International Standard Atmosphere (isa) within troposhpere for:
  ! Values taken from:
  ! https://en.wikipedia.org/wiki/International_Standard_Atmosphere
  ! Lapse rate (degreeC/km) for levels below 11,000 gpm
  REAL(r_def), PARAMETER :: isa_Lapse_RateL = 6.5e-03_r_def

  ! Lapse rate (degreeC/km) for levels above 11,000 gpm
  REAL(r_def), PARAMETER :: isa_Lapse_RateU = -1.0e-03_r_def

  ! surface pressure (Pa)
  REAL(r_def), PARAMETER :: isa_Press_Bot = 101325.0_r_def

  ! pressure (Pa) at 11,000 gpm
  REAL(r_def), PARAMETER :: isa_Press_Mid = 22632.0_r_def

  ! pressure (Pa) at 20,000 gpm
  REAL(r_def), PARAMETER :: isa_Press_Top = 5474.87_r_def

  ! Surface temperature (K)
  REAL(r_def), PARAMETER :: isa_Temp_Bot = 288.15_r_def

  ! Temperature (K) of isothermal layer - at tropopause
  REAL(r_def), PARAMETER :: isa_Temp_Top = 216.65_r_def

  ! Height limit (gpm) for std lower lapse rate
  REAL(r_def), PARAMETER :: Gpm1 = 11000.0_r_def

  ! Height (gpm) of top of isothermal layer
  REAL(r_def), PARAMETER :: Gpm2 = 20000.0_r_def

END MODULE science_aviation_constants_mod

