!----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief  Sets scalar gas concentrations for GHGs and other lower boundary
!!         conditions required by UKCA chemistry schemes.
!!         Values can be specified in namelist or obtained from file.

module ukca_scalar_gas_init_mod

  use constants_mod, only : r_def

  implicit none

  ! Currently set to Y2000 mean values (MMR)
  real(r_def), parameter, public :: atmospheric_ccl4     = 4.442e-10_r_def
  real(r_def), parameter, public :: atmospheric_cfc113   = 4.71e-10_r_def
  real(r_def), parameter, public :: atmospheric_cfc114   = 9.62e-11_r_def
  real(r_def), parameter, public :: atmospheric_cfc115   = 4.49e-11_r_def
  real(r_def), parameter, public :: atmospheric_cfc11    = 1.11e-09_r_def
  real(r_def), parameter, public :: atmospheric_cfc12    = 2.18e-09_r_def
  real(r_def), parameter, public :: atmospheric_ch2br2   = 18.0186e-12_r_def
  real(r_def), parameter, public :: atmospheric_chbr3    = 6.97904e-12_r_def
  real(r_def), parameter, public :: atmospheric_ch4      = 9.697e-07_r_def
  real(r_def), parameter, public :: atmospheric_co2      = 5.58e-04_r_def
  real(r_def), parameter, public :: atmospheric_csul     = 1.000e-09_r_def
  real(r_def), parameter, public :: atmospheric_h1202    = 3.788e-13_r_def
  real(r_def), parameter, public :: atmospheric_h1211    = 2.14e-11_r_def
  real(r_def), parameter, public :: atmospheric_h1301    = 1.70e-11_r_def
  real(r_def), parameter, public :: atmospheric_h2       = 3.453e-08_r_def
  real(r_def), parameter, public :: atmospheric_h2402    = 3.87e-12_r_def
  real(r_def), parameter, public :: atmospheric_hfc125   = 6.937e-12_r_def
  real(r_def), parameter, public :: atmospheric_hfc134a  = 5.158e-12_r_def
  real(r_def), parameter, public :: atmospheric_hcfc141b = 4.092e-11_r_def
  real(r_def), parameter, public :: atmospheric_hcfc142b = 3.637e-11_r_def
  real(r_def), parameter, public :: atmospheric_hcfc22   = 2.851e-10_r_def
  real(r_def), parameter, public :: atmospheric_mebr     = 3.009e-11_r_def
  real(r_def), parameter, public :: atmospheric_meccl3   = 2.775e-10_r_def
  real(r_def), parameter, public :: atmospheric_mecl     = 9.378e-10_r_def
  real(r_def), parameter, public :: atmospheric_n2       = 7.547e-01_r_def
  real(r_def), parameter, public :: atmospheric_n2o      = 4.800e-07_r_def
  real(r_def), parameter, public :: atmospheric_o2       = 2.314e-01_r_def

end module ukca_scalar_gas_init_mod
