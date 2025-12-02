!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Burrows & Sharp (1999) gas abundances for hot Jupiters
!> @details The Burrows & Sharp (1999) chemistry scheme is applicable
!>          to hot Jupiter atmospheres and calculates equilibrium
!>          gas abundances depending on pressure and temperature.

module flexchem_bs1999_kernel_mod

use argument_mod,      only : arg_type, &
                              GH_FIELD, &
                              GH_REAL, GH_READ, GH_READWRITE, &
                              CELL_COLUMN
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type
use log_mod,           only : log_event, log_scratch_space, LOG_LEVEL_ERROR
! SOCRATES module with gas molar masses
use gas_list_pcf,      only : molar_weight, ip_h2o, ip_co, ip_ch4, ip_nh3,     &
                              ip_tio, ip_vo, ip_h2, ip_he, ip_na, ip_k, ip_li, &
                              ip_rb, ip_cs

implicit none

private
!-------------------------------------------------------------------------------
! Private types used in this module
!-------------------------------------------------------------------------------
! Parameters
! Number fraction of hydrogen
real(r_def), parameter :: hydrogen_num_frac  = 0.91183_r_def
! Number fraction of helium
real(r_def), parameter :: helium_num_frac = 1.0_r_def - hydrogen_num_frac
! Universal gas constant
real(r_def), parameter :: r_gas = 8.3143_r_def
! Ideal gas constant in cal/(mol*K)
real(r_def), parameter :: r_gas_cal = 1.9858775_r_def
! Number of atmospheres per Pa
real(r_def), parameter :: atm_per_pa = 1.0_r_def / 1.01325e+05_r_def
! Average number of oxygen atoms removed per silicon atom
real(r_def), parameter :: x_si = 3.28_r_def
! Temperature below which oxygen is removed by silicon atoms
real(r_def), parameter :: t_remove_o = 1500_r_def
! Characteristic transition scale in T over which O is removed
real(r_def), parameter :: t_trans_scale_remove_o = 10.0_r_def
! Characteristic transition scale in T over which the abundance changes
real(r_def), parameter ::                                                      &
    t_trans_scale_tio = 10.0_r_def,                                            &
    t_trans_scale_vo  = 10.0_r_def,                                            &
    t_trans_scale_na  = 20.0_r_def,                                            &
    t_trans_scale_k   = 20.0_r_def,                                            &
    t_trans_scale_li  = 20.0_r_def,                                            &
    t_trans_scale_rb  = 20.0_r_def,                                            &
    t_trans_scale_cs  = 20.0_r_def

! Elemental abundances as log(epsilon_z) = log(N_z/N_H) + 12, where N_z is the
! number of atoms of species z
real(r_def), parameter ::                                                      &
    log_eps_c  = 9.50_r_def,                                                   &
    log_eps_n  = 8.86_r_def,                                                   &
    log_eps_o  = 9.76_r_def,                                                   &
    log_eps_na = 7.24_r_def,                                                   &
    log_eps_si = 8.51_r_def,                                                   &
    log_eps_k  = 6.03_r_def,                                                   &
    log_eps_ti = 5.95_r_def,                                                   &
    log_eps_v  = 4.93_r_def,                                                   &
    log_eps_li = 4.26_r_def,                                                   &
    log_eps_rb = 3.52_r_def,                                                   &
    log_eps_cs = 2.08_r_def

! Convert elemental abundances to number of atoms relative to H:
! num_atoms_z = N_z/N_H
real(r_def), parameter ::                                                      &
    num_atoms_c  = (10.0_r_def)**(log_eps_c  - 12.0_r_def),                    &
    num_atoms_n  = (10.0_r_def)**(log_eps_n  - 12.0_r_def),                    &
    num_atoms_o  = (10.0_r_def)**(log_eps_o  - 12.0_r_def),                    &
    num_atoms_na = (10.0_r_def)**(log_eps_na - 12.0_r_def),                    &
    num_atoms_si = (10.0_r_def)**(log_eps_si - 12.0_r_def),                    &
    num_atoms_k  = (10.0_r_def)**(log_eps_k  - 12.0_r_def),                    &
    num_atoms_ti = (10.0_r_def)**(log_eps_ti - 12.0_r_def),                    &
    num_atoms_v  = (10.0_r_def)**(log_eps_v  - 12.0_r_def),                    &
    num_atoms_li = (10.0_r_def)**(log_eps_li - 12.0_r_def),                    &
    num_atoms_rb = (10.0_r_def)**(log_eps_rb - 12.0_r_def),                    &
    num_atoms_cs = (10.0_r_def)**(log_eps_cs - 12.0_r_def)

! Polynomial fits to condensation curves
real(r_def), parameter ::                                                      &
  ! Titanium and Vanadium Oxides
  poly_highp_tio_vo(2) =  [                                                    &
    -3.96274342e-05_r_def,                                                     &
      5.20741797e-04_r_def ],                                                  &
  ! Sodium
  poly_highp_na(2) = [                                                         &
    -5.58256919e-05_r_def,                                                     &
      8.81851644e-04_r_def ],                                                  &
  ! Potassium
  poly_highp_k(2) = [                                                          &
    -5.46977180e-05_r_def,                                                     &
      8.19104478e-04_r_def ],                                                  &
  ! Lithium
  poly_highp_li(2) = [                                                         &
    -3.50995394e-05_r_def,                                                     &
      6.51993843e-04_r_def ],                                                  &
  ! Rubidium
  poly_highp_rb(2) = [                                                         &
    -6.06654087e-05_r_def,                                                     &
      8.09569974e-04_r_def ],                                                  &
  ! Cesium
  poly_highp_cs(2) = [                                                         &
    -5.29210264e-05_r_def,                                                     &
      7.71577097e-04_r_def ]

real(r_def), parameter ::                                                      &
  ! Titanium and Vanadium Oxides
  poly_lowp_tio_vo(2)  = [                                                     &
    -2.91788901e-05_r_def,                                                     &
      5.11801742e-04_r_def ],                                                  &
  ! Sodium
  poly_lowp_na(2) = [                                                          &
    -6.69921629e-05_r_def,                                                     &
      8.90116829e-04_r_def ],                                                  &
  ! Potassium
  poly_lowp_k(2) = [                                                           &
    -6.46633572e-05_r_def,                                                     &
      8.29549449e-04_r_def ],                                                  &
  ! Lithium
  poly_lowp_li(2) = [                                                          &
    -3.55469185e-05_r_def,                                                     &
      6.52116945e-04_r_def ],                                                  &
  ! Rubidium
  poly_lowp_rb(2) = [                                                          &
    -3.19328287e-05_r_def,                                                     &
    8.69542964e-04_r_def ],                                                    &
  ! Cesium
  poly_lowp_cs(2) = [                                                          &
    -3.85306167e-05_r_def,                                                     &
      7.63040762e-04_r_def ]

! Transition point between polynomial fits
real(r_def), parameter ::                                                      &
  ! Titanium and Vanadium Oxides
  log_p_trans_tio_vo =  1.0_r_def,                                             &
  ! Sodium
  log_p_trans_na =      1.0_r_def,                                             &
  ! Potassium
  log_p_trans_k =       1.0_r_def,                                             &
  ! Lithium
  log_p_trans_li =      1.0_r_def,                                             &
  ! Rubidium
  log_p_trans_rb =     -2.0_r_def,                                             &
  ! Cesium
  log_p_trans_cs =      1.0_r_def

! Maximum value to take the exponential of - returns about 5e21
real(r_def), parameter :: max_exp = 50.0_r_def

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: flexchem_bs1999_kernel_type
  private
  type(arg_type) :: meta_args(15) = (/ &
    arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), & ! pressure_in_wth
    arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), & ! temperature_in_wth
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! h2o
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! co
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! ch4
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! nh3
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! h2
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! he
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! tio
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! vo
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! na
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! k
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! li
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), & ! rb
    arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta)  & ! cs
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: flexchem_bs1999_code
end type

public :: flexchem_bs1999_code

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

  !> @brief Kernel interface to the Burrows & Sharp (1999) scheme
  !> @param[in]     nlayers            Number of layers
  !> @param[in]     pressure_in_wth    Pressure in wth space
  !> @param[in]     temperature_in_wth Temperature in wth space
  !> @param[in,out] h2o                Water
  !> @param[in,out] co                 Carbon monoxide
  !> @param[in,out] ch4                Methane
  !> @param[in,out] nh3                Ammonia
  !> @param[in,out] h2                 Hydrogen
  !> @param[in,out] he                 Helium
  !> @param[in,out] tio                Titanium oxide
  !> @param[in,out] vo                 Vanadium oxide
  !> @param[in,out] na                 Sodium
  !> @param[in,out] k                  Potassium
  !> @param[in,out] li                 Lithium
  !> @param[in,out] rb                 Rubidium
  !> @param[in,out] cs                 Cesium
  !> @param[in]     ndf_wth            No. DOFs per cell for wth space
  !> @param[in]     undf_wth           No. unique of DOFs for wth space
  !> @param[in]     map_wth            Dofmap for wth space column base cell
  subroutine flexchem_bs1999_code(nlayers,                                     &
                                  pressure_in_wth, temperature_in_wth,         &
                                  h2o, co, ch4, nh3, h2, he,                   &
                                  tio, vo, na, k, li, rb, cs,                  &
                                  ndf_wth, undf_wth, map_wth)

    implicit none

    ! Arguments
    integer(i_def), intent(in) :: nlayers
    integer(i_def), intent(in) :: ndf_wth, undf_wth

    integer(i_def), dimension(ndf_wth), intent(in) :: map_wth

    real(r_def), dimension(undf_wth), intent(in) :: pressure_in_wth
    real(r_def), dimension(undf_wth), intent(in) :: temperature_in_wth

    real(r_def), dimension(undf_wth), intent(inout) ::                         &
      h2o, co, ch4, nh3, h2, he, tio, vo, na, k, li, rb, cs

    ! Loop indices
    integer(i_def) :: k_layer, j
    real(r_def) :: temperature_capped

    do k_layer = 0, nlayers
      j = map_wth(1) + k_layer

      ! Cap temperature at a minimum value of 300K
      temperature_capped = max(temperature_in_wth(j), 300.0_r_def)
      ! Hydrogen
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_h2,   &
        h2(j))
      ! Helium
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_he,   &
        he(j))
      ! Carbon monoxide
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_co,   &
        co(j))
      ! Methane
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_ch4,  &
        ch4(j))
      ! Water
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_h2o,  &
        h2o(j))
      ! Ammonia
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_nh3,  &
        nh3(j))
      ! Cesium
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_cs,   &
        cs(j))
      ! Potassium
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_k,    &
        k(j))
      ! Lithium
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_li,   &
        li(j))
      ! Sodium
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_na,   &
        na(j))
      ! Rubidium
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_rb,   &
        rb(j))
      ! Titanium oxide
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_tio,  &
        tio(j))
      ! Vanadium oxide
      call calc_mmr_bs1999(pressure_in_wth(j), temperature_capped, ip_vo,   &
        vo(j))
    end do

  end subroutine flexchem_bs1999_code

  !> @brief Burrows & Sharp (1999) scheme
  !> @param[in]  p    Pressure (Pa)
  !> @param[in]  t    Temperature (K)
  !> @param[in]  iump Identifier of molecule
  !> @param[out] mmr  Mass mixing ratio
  subroutine calc_mmr_bs1999(p, t, iump, mmr)

    implicit none

    real(r_def), intent(in) :: p, t
    integer(i_def), intent(in) :: iump

    real(r_def), intent(out) :: mmr

    ! Local variables
    ! Equilibrium coefficients
    real(r_def) :: eq_coef
    ! Condensation temperature
    real(r_def) :: t_cond
    ! Number of oxygen atoms once depleted
    real(r_def) :: oxygen_atoms_once_depleted
    ! Working variable in the calculation of CO equilibrium abundance
    real(r_def) :: co_work
    ! Working variable in the calculation of NH3 equilibrium abundance
    real(r_def) :: nh3_work
    ! Log (base 10) of pressure in bars
    real(r_def) :: log_p_bar
    ! Fraction of oxygen available to deplete that is depleted
    real(r_def) :: oxygen_frac_remove
    ! Total mass density
    real(r_def) :: tot_mass_dens
    ! H2 and He partial pressures
    real(r_def) :: pp_h2, pp_he

    log_p_bar = log10(p / 1.0e+05_r_def)

    ! Calculate fraction of available oxygen to deplete
    ! that have depleted given temperature in layer
    oxygen_frac_remove = 1.0_r_def /                                           &
      ( exp(min((t - t_remove_o) / t_trans_scale_remove_o,max_exp)) + 1.0_r_def )

    ! Calculate H2 and He partial pressures
    pp_h2 = p * (hydrogen_num_frac / 2.0_r_def) / &
      (hydrogen_num_frac / 2.0_r_def + helium_num_frac)
    pp_he = p * helium_num_frac / &
      (hydrogen_num_frac / 2.0_r_def + helium_num_frac)

    ! Total atmospheric mass density
    ! NB: All mass is assumed to be held in H2 and He. Should possibly use
    ! a mean molecular mass for all molecules instead to find total mass
    ! of atmosphere. pp_h2 converted from Pa to bar.
    tot_mass_dens = 1.0e-3_r_def * (molar_weight(ip_h2) * pp_h2 +              &
      molar_weight(ip_he) * pp_he) / (r_gas * t)

    ! Calculate partial pressures of gases relative to pp_h2. The unit of
    ! pp_h2 in Burrows & Sharp, ApJ, 1999 are in atmospheres, must convert
    ! pp_h2 to atmospheres before use in these formulas.
    select case(iump)

    case (ip_h2)
      ! H2 partial pressure relative to H2 partial pressure
      mmr = 1.0_r_def ! pp_h2/pp_h2

    case (ip_he)
      ! He partial pressure relative to H2 partial pressure
      mmr = pp_he / pp_h2

    case (ip_co, ip_ch4, ip_h2o)
      ! CO partial pressure relative to H2 partial pressure
      eq_coef = calc_k1(t)
      oxygen_atoms_once_depleted = num_atoms_o                                 &
        - x_si * num_atoms_si * oxygen_frac_remove

      co_work = num_atoms_c + oxygen_atoms_once_depleted +                     &
        (pp_h2 * atm_per_pa)**2 / (2.0_r_def * eq_coef)

      mmr = co_work -                                                          &
        sqrt(co_work**2 - 4.0_r_def * num_atoms_c * oxygen_atoms_once_depleted)

      ! CH4 partial pressure relative to H2 partial pressure
      if (iump == ip_ch4) then
        mmr = 2.0_r_def * num_atoms_c - mmr ! i.e. minus CO mmr
      end if

      ! H2O partial pressure relative to H2 partial pressure
      if (iump == ip_h2o) then
        mmr = 2.0_r_def * oxygen_atoms_once_depleted - mmr ! i.e. minus CO mmr
      end if

    case (ip_nh3)
      ! NH3 partial pressure relative to H2 partial pressure
      eq_coef = calc_k2(t)
      nh3_work = (pp_h2*atm_per_pa)**2 / (8.0_r_def * eq_coef)

      if ( num_atoms_n / nh3_work < 1.0e-10_r_def ) then
        mmr = 2.0_r_def * num_atoms_n
      else
        mmr = 2.0_r_def * (                                                    &
          sqrt( nh3_work * (2.0_r_def * num_atoms_n + nh3_work) ) - nh3_work )
      end if

    case (ip_tio)
      ! TiO partial pressure relative to H2 partial pressure
      t_cond = calc_cond_temp(log_p_bar, log_p_trans_tio_vo,                   &
        poly_highp_tio_vo, poly_lowp_tio_vo)
      mmr = 2.0_r_def * num_atoms_ti /                                         &
            (exp(min(-(t - t_cond)/t_trans_scale_tio,max_exp)) + 1.0_r_def)

    case (ip_vo)
      ! VO partial pressure relative to H2 partial pressure
      t_cond = calc_cond_temp(log_p_bar, log_p_trans_tio_vo,                   &
        poly_highp_tio_vo, poly_lowp_tio_vo)
      mmr = 2.0_r_def * num_atoms_V /                                          &
            (exp(min(-(t - t_cond) / t_trans_scale_vo,max_exp)) + 1.0_r_def)

    case (ip_na)
      ! Sodium partial pressure relative to H2 partial pressure
      t_cond = calc_cond_temp(log_p_bar, log_p_trans_na,                       &
        poly_highp_na, poly_lowp_na)
      mmr = 2.0_r_def * num_atoms_na /                                         &
        (exp(min(-(t - t_cond) / t_trans_scale_na,max_exp)) + 1.0_r_def)

    case (ip_k)
      ! Potassium partial pressure relative to H2 partial pressure
      t_cond = calc_cond_temp(log_p_bar, log_p_trans_k,                        &
        poly_highp_k, poly_lowp_k)
      mmr = 2.0_r_def * num_atoms_k /                                          &
        (exp(min(-(t - t_cond) / t_trans_scale_k,max_exp)) + 1.0_r_def)

    case (ip_li)
      ! Lithium partial pressure relative to H2 partial pressure
      t_cond = calc_cond_temp(log_p_bar, log_p_trans_li,                       &
        poly_highp_li, poly_lowp_li)
      mmr = 2.0_r_def * num_atoms_li /                                         &
        (exp(min(-(t - t_cond) / t_trans_scale_li,max_exp)) + 1.0_r_def)

    case (ip_rb)
      ! Rubidium partial pressure relative to H2 partial pressure
      t_cond = calc_cond_temp(log_p_bar, log_p_trans_rb,                       &
        poly_highp_rb, poly_lowp_rb)
      mmr = 2.0_r_def * num_atoms_rb /                                         &
        (exp(min(-(t - t_cond) / t_trans_scale_rb,max_exp)) + 1.0_r_def)

    case (ip_cs)
      ! Cesium partial pressure relative to H2 partial pressure
      t_cond = calc_cond_temp(log_p_bar, log_p_trans_cs,                       &
        poly_highp_cs, poly_lowp_cs)
      mmr = 2.0_r_def * num_atoms_cs /                                         &
        (exp(min(-(t - t_cond) / t_trans_scale_cs,max_exp)) + 1.0_r_def)

    case default
      write( log_scratch_space, '(A,I3,A)' )                                  &
        'Gas ', iump, ' is not part of the BS1999 scheme'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

    ! Convert to actual pressure (not relative to H2 partial pressure), and
    ! convert to mass density
    mmr = mmr * pp_h2 * molar_weight(iump) * 1.0e-3_r_def / (r_gas * t)

    ! Convert to molecular mass fraction
    mmr = max(0.0_r_def, mmr / tot_mass_dens)
  end subroutine calc_mmr_bs1999

  !> @brief Calculate condensation temperature from polynomial fits
  !> @param[in] log_p_bar   Log (base 10) of the total gas pressure in bar
  !> @param[in] log_p_trans Log (base 10) of transition between polynomial fits
  !> @param[in] poly_highp  Polynomial fit for high pressure
  !> @param[in] poly_lowp   Polynomial fit for low pressure
  !> @return    t_cond      Condensation temperature (K)
  pure function calc_cond_temp(log_p_bar, log_p_trans, poly_highp, poly_lowp)  &
    result(t_cond)

    implicit none

    real(r_def), intent(in) :: log_p_bar
    real(r_def), intent(in) :: poly_highp(2)
    real(r_def), intent(in) :: poly_lowp(2)
    real(r_def), intent(in) :: log_p_trans
    real(r_def) :: t_cond
    ! Local variables
    real(r_def) :: work

    work = 1.0_r_def /                                                         &
      (exp(min(-(log_p_bar - log_p_trans) / 0.1_r_def,max_exp)) + 1.0_r_def)
    t_cond = 1.0_r_def / ( (poly_lowp(1) * log_p_bar + poly_lowp(2))) *        &
      (1.0_r_def - work) +                                                     &
      1.0_r_def / ( (poly_highp(1)*log_p_bar + poly_highp(2))) * work
  end function calc_cond_temp


  !> @brief Calculate the equilibrium coefficient k_1
  !> @param[in] t Temperature (K)
  !> @return    k Equilibrium coefficient k_1
  pure function calc_k1(t) result (k)

    implicit none

    real(r_def), intent(in) :: t
    real(r_def)             :: k
    ! Equilibrium coefficients from Burrows & Sharp, ApJ, 1999
    real(r_def), parameter :: eqcoeff_1(5) = [                                 &
      1.106131e+06_r_def,                                                      &
      -5.6895e+04_r_def,                                                       &
      62.565_r_def,                                                            &
      -5.81396e-04_r_def,                                                      &
      2.346515e-08_r_def ]

    k = exp( (eqcoeff_1(1) / t + eqcoeff_1(2)                                  &
            + eqcoeff_1(3) * t                                                 &
            + eqcoeff_1(4) * t**2                                              &
            + eqcoeff_1(5) * t**3)                                             &
            / (r_gas_cal * t) )

  end function calc_k1

  !> @brief Calculate the equilibrium coefficient k_2
  !> @param[in] t Temperature (K)
  !> @return    k Equilibrium coefficient k_2
  pure function calc_k2(t) result(k)

    implicit none

    real(r_def), intent(in) :: t
    real(r_def)             :: k
    ! Equilibrium coefficients from Burrows & Sharp, ApJ, 1999
    real(r_def), parameter :: eqcoeff_2(5) = [                                 &
      8.16413e+05_r_def,                                                       &
      -2.9109e+04_r_def,                                                       &
      58.5878_r_def,                                                           &
      -7.8284e-04_r_def,                                                       &
      4.729048e-08_r_def ]

    k = exp( (eqcoeff_2(1) / t + eqcoeff_2(2)                                  &
            + eqcoeff_2(3) * t                                                 &
            + eqcoeff_2(4) * t**2                                              &
            + eqcoeff_2(5) * t**3)                                             &
            / (r_gas_cal * t) )
  end function calc_k2
end module flexchem_bs1999_kernel_mod
