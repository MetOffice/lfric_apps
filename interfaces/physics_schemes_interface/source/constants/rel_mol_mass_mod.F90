!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic interface module for UM code (rel_mol_mass_mod)
!----------------------------------------------------------------------------

module rel_mol_mass_mod

  use constants_mod, only: r_um
  use science_rel_mol_mass_mod, only: relative_molecular_mass_s,    &
                                      relative_molecular_mass_h2o2, &
                                      relative_molecular_mass_o3,   &
                                      relative_molecular_mass_h2o,  &
                                      relative_molecular_mass_dry_air

  implicit none

  private

  ! The following variables have been hidden as they are not currently
  ! required to build the extracted UM code. They have been left in
  ! in case they are required as more UM code is drawn into the lfric_atm
  ! build. Should they be required at a later date, they should simply be
  ! made public.

  ! Disabled variables:
  !   rmm_s, rmm_h2o2, rmm_o3, rmm_air, rmm_w

  ! Relative Molecular Mass [kg/mole]
  ! Sulphur
  real(r_um), parameter :: rmm_s    = real(relative_molecular_mass_s,    r_um)
  ! Hydrogen Peroxide
  real(r_um), parameter :: rmm_h2o2 = real(relative_molecular_mass_h2o2, r_um)
  ! Ozone
  real(r_um), parameter :: rmm_o3   = real(relative_molecular_mass_o3,   r_um)
  ! Dry Air
  real(r_um), parameter :: rmm_air  = real(relative_molecular_mass_h2o,  r_um)
  ! Water
  real(r_um), parameter :: rmm_w    = real(relative_molecular_mass_dry_air, r_um)

end module rel_mol_mass_mod
