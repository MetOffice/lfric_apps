! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code.
! -----------------------------------------------------------------------------
!> @brief LFRic aviation constants module
!>
module science_aviation_constants_mod

  use constants_mod, only: r_def
  use science_conversions_mod, only: feet_to_metres

  implicit none

  private
  public :: mtokft, isa_lapse_ratel, isa_lapse_rateu,              &
      isa_press_bot, isa_press_mid, isa_press_top,                 &
      isa_temp_bot, isa_temp_top, gpm1, gpm2,                      &
      parcel_mixed_layer_depth, parcel_surface_layer_depth,        &
      parcel_min_cape, parcel_min_cin, parcel_qsat_delta_t,        &
      parcel_screen_height, parcel_top_start_pressure,             &
      inflow_min_cape, inflow_min_cin,                             &
      bolton_a, bolton_b, bolton_c, bolton_d

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

  ! Undilute parcel ascents used for CAPE and CIN diagnostics.
  ! Depth (m) of the near-surface mixed layer averaged to form the parcel.
  real(r_def), parameter :: parcel_mixed_layer_depth = 500.0_r_def

  ! Depth (m) of the surface layer represented by a 1.5 m parcel, taken as
  ! twice the screen height.
  real(r_def), parameter :: parcel_surface_layer_depth = 3.0_r_def

  ! CAPE (J/kg) above which a parcel is treated as buoyant, so that the
  ! next negatively buoyant layer marks its level of neutral buoyancy.
  real(r_def), parameter :: parcel_min_cape = 1.0_r_def

  ! Most negative CIN (J/kg); an ascent stops once this is exceeded.
  real(r_def), parameter :: parcel_min_cin = -1000.0_r_def

  ! Temperature step (K) used to estimate the rate of change of saturation
  ! mixing ratio with temperature along a saturated ascent.
  real(r_def), parameter :: parcel_qsat_delta_t = 0.2_r_def

  ! Height (m) of the screen-level parcel above the surface.
  real(r_def), parameter :: parcel_screen_height = 1.5_r_def

  ! Model levels at or below this pressure (Pa), that is at or above it in
  ! height, are not used as starting levels for the most unstable parcel.
  real(r_def), parameter :: parcel_top_start_pressure = 50000.0_r_def

  ! Minimum CAPE (J/kg) and CIN (J/kg) of a parcel for its starting level to
  ! lie in the effective inflow layer of a convective storm, from Thompson,
  ! Mead and Edwards (2007), Wea. Forecasting 22, 102-115.
  real(r_def), parameter :: inflow_min_cape = 100.0_r_def
  real(r_def), parameter :: inflow_min_cin = -250.0_r_def

  ! Empirical coefficients for the temperature at the lifting condensation
  ! level, from equation 21 of Bolton (1980), Mon. Wea. Rev. 108, 1046-1053:
  ! T_lcl = bolton_b / (bolton_c ln(T) - ln(e) - bolton_d) + bolton_a,
  ! with temperatures in K and vapour pressure e in hPa.
  real(r_def), parameter :: bolton_a = 55.0_r_def
  real(r_def), parameter :: bolton_b = 2840.0_r_def
  real(r_def), parameter :: bolton_c = 3.5_r_def
  real(r_def), parameter :: bolton_d = 4.805_r_def

end module science_aviation_constants_mod

