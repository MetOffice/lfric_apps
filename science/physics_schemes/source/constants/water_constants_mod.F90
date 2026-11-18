!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic interface module for UM code (water_constants_mod)
!----------------------------------------------------------------------------

module water_constants_mod

  use, intrinsic :: iso_fortran_env, only: real32
  use constants_mod, only : r_um, r_bl
  use driver_water_constants_mod, only: t_freeze_h2o_sea,             &
                                        t_freeze_h2o,                 &
                                        density_h2o,                  &
                                        density_h2o_sea,              &
                                        density_ice,                  &
                                        latent_heat_h2o_condensation, &
                                        latent_heat_h2o_fusion,       &
                                        heat_capacity_h2o_vapour,     &
                                        heat_capacity_h2o,            &
                                        heat_capacity_ice,            &
                                        jules_dpsidt

  implicit none

  private
  public :: hcapi, hcapw, hcapv, lc, lf, rho_ice, rho_water, rhosea, tfs, tm, &
            dpsidt, tm_bl, lc_bl, tm_32b, lc_32b

  !-----------------------------------------------------------------------
  ! Parameters names contained here are fixed, as these names are
  ! what the UM code base with look for.
  !-----------------------------------------------------------------------

  ! Temperature at which sea water freezes, [K]
  real(r_um), parameter :: tfs = real(t_freeze_h2o_sea, r_um)

  ! Temperature at which fresh water freezes and ice melts, [K]
  real(r_um), parameter :: tm  = real(t_freeze_h2o, r_um)
  real(real32), parameter :: tm_32b = real(tm, real32)
  real(r_bl), parameter :: tm_bl  = real(t_freeze_h2o, r_bl)

  ! Density of pure water [kg/m3]
  real(r_um), parameter :: rho_water = real(density_h2o, r_um)

  ! Density of sea water [kg/m3]
  real(r_um), parameter :: rhosea = real(density_h2o_sea, r_um)

  ! Density of ice [kg/m3]
  real(r_um), parameter :: rho_ice = real(density_ice, r_um)

  ! Latent heat of condensation of water at 0 degC [J/kg]
  real(r_um), parameter :: lc = real(latent_heat_h2o_condensation, r_um)
  real(real32), parameter :: lc_32b = real(lc, real32)
  real(r_bl), parameter :: lc_bl = real(latent_heat_h2o_condensation, r_bl)

  ! Latent heat of fusion of water at 0 degC [J/kg]
  real(r_um), parameter :: lf = real(latent_heat_h2o_fusion, r_um)

  ! Specific heat capacity of water vapour [J/kg/K]
  real(r_um), parameter :: hcapv = real(heat_capacity_h2o_vapour, r_um)

  ! Specific heat capacity of water [J/kg/K]
  real(r_um), parameter :: hcapw = real(heat_capacity_h2o, r_um)

  ! Specific heat capacity of ice [J/kg/K]
  real(r_um), parameter :: hcapi = real(heat_capacity_ice, r_um)

  ! Rate of change of soil matrix potential with temperature at
  ! equilibrium between water and ice in partially frozen soil [m/K].
  real(r_um), parameter :: dpsidt = real(jules_dpsidt, r_um)

end module water_constants_mod
