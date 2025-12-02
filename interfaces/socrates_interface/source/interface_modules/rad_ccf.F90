!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Radiation Constant Configuration File

module rad_ccf

  use constants_mod,        only: lfric_pi => pi
  use science_chemistry_constants_mod, &
                            only: avogadro, &
                                  boltzmann, &
                                  lfric_stefan_boltzmann => stefan_boltzmann, &
                                  planck_constant, &
                                  speed_of_light
  use science_conversions_mod, &
                            only: lfric_seconds_per_day => seconds_per_day
  use science_rel_mol_mass_mod, &
                            only: relative_molecular_mass_dry_air
  use driver_water_constants_mod, &
                            only: density_h2o
  use realtype_rd,          only: RealK

  implicit none
  private
  public :: k_boltzmann, n_avogadro, stefan_boltzmann, pi, seconds_per_day, &
            mol_weight_air, repsilon, rho_water, &
            h_planck, c_light
  public :: set_socrates_constants

  real(RealK), parameter :: k_boltzmann  = real(boltzmann, RealK)
  real(RealK), parameter :: n_avogadro   = real(avogadro, RealK)
  real(RealK), parameter :: stefan_boltzmann &
                                         = real(lfric_stefan_boltzmann, RealK)
  real(RealK), parameter :: h_planck     = real(planck_constant, RealK)
  real(RealK), parameter :: c_light      = real(speed_of_light, RealK)
  real(RealK), parameter :: pi           = real(lfric_pi, RealK)
  real(RealK), parameter :: rho_water    = real(density_h2o, RealK)
  real(RealK), parameter :: seconds_per_day &
                              = real(lfric_seconds_per_day, RealK)
  real(RealK), parameter :: mol_weight_air &
                              = real(relative_molecular_mass_dry_air, RealK)

  ! Ratio of molecular weights of water and dry air
  real(RealK), protected :: repsilon

contains

subroutine set_socrates_constants()

  use planet_constants_mod, only: lfric_repsilon => repsilon

  implicit none

  repsilon = real(lfric_repsilon, RealK)

end subroutine set_socrates_constants

end module rad_ccf
