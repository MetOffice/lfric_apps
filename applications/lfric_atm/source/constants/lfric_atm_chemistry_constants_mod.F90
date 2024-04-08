!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic chemistry constants module
!----------------------------------------------------------------------------
!
module lfric_atm_chemistry_constants_mod

  use constants_mod, only: r_def

  implicit none

  private
  public :: avogadro, boltzmann, stefan_boltzmann, planck_constant, &
            speed_of_light, density_so4, mean_free_path_reference,  &
            temperature_mean_free_path_reference,                   &
            pressure_mean_free_path_reference

  !> Number of molecules in 1 mole
  real(r_def), parameter :: avogadro = 6.022e23_r_def

  !> Boltzmanns constant (J/K)
  real(r_def), parameter :: boltzmann = 1.3804e-23_r_def

  !> Stefan-Boltzmann constant (W/m^2/K^4).
  real(r_def), parameter :: stefan_boltzmann = 5.670374419e-08_r_def

  !> Planck constant (J s)
  real(r_def), parameter :: planck_constant = 6.626070150e-34_r_def

  !> Speed of light in a vacuum (m/s)
  real(r_def), parameter :: speed_of_light = 2.99792458e+08_r_def

  !> Density of SO4 particle (kg/m^3)
  real(r_def), parameter :: density_so4 = 1769.0_r_def

  !> @name Mean Free Path constants
  !> @{
  real(r_def), parameter :: mean_free_path_reference = 6.6e-8_r_def
  !< [m]
  real(r_def), parameter :: temperature_mean_free_path_reference = 293.15_r_def
  !< [K]
  real(r_def), parameter :: pressure_mean_free_path_reference = 1.01325e5_r_def
  !< [Pa]
  !> @}

end module lfric_atm_chemistry_constants_mod

