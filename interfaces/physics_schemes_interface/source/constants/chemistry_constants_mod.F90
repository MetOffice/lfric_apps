!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic inteface module for UM code (chemistry_constants_mod)
!----------------------------------------------------------------------------

module chemistry_constants_mod

  use constants_mod, only: r_um
  use science_chemistry_constants_mod, only:                                   &
                                         lfric_avogadro => avogadro,           &
                                         lfric_boltzmann => boltzmann,         &
                                         density_so4,                          &
                                         mean_free_path_reference,             &
                                         temperature_mean_free_path_reference, &
                                         pressure_mean_free_path_reference

  implicit none

  private
  public :: avogadro, boltzmann, rho_so4, mfp_ref, pref_mfp, tref_mfp

  ! Number of molecules in 1 mole
  real(r_um), parameter :: avogadro = real(lfric_avogadro, r_um)

  ! Boltzmanns constant [J/K]
  real(r_um), parameter :: boltzmann = real(lfric_boltzmann, r_um)

  ! Density of so4 particle [kg/m^3]
  real(r_um), parameter :: rho_so4 = real(density_so4, r_um)

  ! Mean free path
  ! Ref value [m]
  real(r_um), parameter :: mfp_ref =  real(mean_free_path_reference, r_um)
  ! Ref temperature [K]
  real(r_um), parameter :: tref_mfp= real(temperature_mean_free_path_reference,&
                                          r_um)
  ! Ref pressure [Pa]
  real(r_um), parameter :: pref_mfp= real(pressure_mean_free_path_reference,   &
                                          r_um)

end module chemistry_constants_mod

