! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsprec_mod

use um_types, only: real_lsprec

!Use statements for parameters
!General atmosphere modules
use conversions_mod,      only:                                                &
  zerodegc_orig => zerodegc,                                                   &
  pi_orig       => pi

use water_constants_mod,  only:                                                &
  lc_orig => lc,                                                               &
  lf_orig => lf

use planet_constants_mod, only:                                                &
  rv_orig => rv

use pc2_constants_mod,    only:                                                &
  wind_shear_factor_orig => wind_shear_factor

!Microphysics modules
use mphys_ice_mod,        only:                                                &
  t_scaling_orig => t_scaling,                                                 &
  qcf0_orig      => qcf0,                                                      &
  t_agg_min_orig => t_agg_min,                                                 &
  qcfmin_orig    => qcfmin,                                                    &
  thomo_orig     => thomo

use lsp_dif_mod,          only:                                                &
  cpwr_orig   => cpwr,                                                         &
  tw1_orig    => tw1,                                                          &
  tw2_orig    => tw2,                                                          &
  tw3_orig    => tw3,                                                          &
  tw4_orig    => tw4,                                                          &
  tw5_orig    => tw5,                                                          &
  tcor1_orig  => tcor1,                                                        &
  tcor2_orig  => tcor2

use lsp_autoc_consts_mod, only:                                                &
  auto_graup_qcf_thresh_orig  => auto_graup_qcf_thresh,                        &
  auto_graup_t_thresh_orig    => auto_graup_t_thresh,                          &
  auto_graup_coeff_orig       => auto_graup_coeff,                             &
  inhomog_rate_orig           => inhomog_rate,                                 &
  inhomog_lim_orig            => inhomog_lim,                                  &
  r_thresh_orig               => r_thresh,                                     &
  r_auto_orig                 => r_auto,                                       &
  n_auto_orig                 => n_auto,                                       &
  power_droplet_auto_orig     => power_droplet_auto,                           &
  power_qcl_auto_orig         => power_qcl_auto,                               &
  power_rho_auto_orig         => power_rho_auto,                               &
  consts_auto_orig            => consts_auto,                                  &
  aut_pref_orig               => aut_pref,                                     &
  aut_nc_orig                 => aut_nc

use mphys_constants_mod,  only:                                                &
  acc_pref_orig   => acc_pref,                                                 &
  acc_qc_orig     => acc_qc,                                                   &
  acc_qr_orig     => acc_qr,                                                   &
  mprog_min_orig  => mprog_min,                                                &
  mprog_abs_orig  => mprog_abs,                                                &
  ntot_land_orig  => ntot_land,                                                &
  ntot_sea_orig   => ntot_sea,                                                 &
  max_as_enh_orig => max_as_enh

use mphys_ice_mod,        only:                                                &
  m0_orig => m0

use mphys_radar_mod,      only:                                                &
  kliq_orig         => kliq,                                                   &
  kice_orig         => kice,                                                   &
  mm6m3_orig        => mm6m3,                                                  &
  ref_lim_orig      => ref_lim,                                                &
  ref_lim_lin_orig  => ref_lim_lin,                                            &
  mr_lim_orig       => mr_lim,                                                 &
  cf_lim_orig       => cf_lim,                                                 &
  nd_lim_orig       => nd_lim,                                                 &
  ref_mom_orig      => ref_mom

implicit none

character(len=*), parameter, private :: ModuleName='LSPREC_MOD'

!Variables to be copied at lsprec precision, allowing full precision to be
!kept in their native modules as some are used all over the UM.

!Common magic numbers
real (kind=real_lsprec), parameter ::                                          &
  zero = 0.0_real_lsprec,                                                      &
  half = 0.5_real_lsprec,                                                      &
  one  = 1.0_real_lsprec,                                                      &
  two  = 2.0_real_lsprec,                                                      &
  ten  = 10.0_real_lsprec

! Very small number for safety-checks to avoid div-by-zero
real (kind=real_lsprec), parameter :: small_number = sqrt(tiny(zero))

!Other PARAMETERS
!lsp_dif_mod
real (kind=real_lsprec), parameter ::                                          &
  cpwr  = real(cpwr_orig, kind=real_lsprec),                                   &
  tw1   = real(tw1_orig, kind=real_lsprec),                                    &
  tw2   = real(tw2_orig, kind=real_lsprec),                                    &
  tw3   = real(tw3_orig, kind=real_lsprec),                                    &
  tw4   = real(tw4_orig, kind=real_lsprec),                                    &
  tw5   = real(tw5_orig, kind=real_lsprec),                                    &
  tcor1 = real(tcor1_orig, kind=real_lsprec),                                  &
  tcor2 = real(tcor2_orig, kind=real_lsprec)

!lsp_autoc_consts_mod
real (kind=real_lsprec), parameter ::                                          &
  auto_graup_qcf_thresh = real(auto_graup_qcf_thresh_orig, kind=real_lsprec),  &
  auto_graup_t_thresh   = real(auto_graup_t_thresh_orig, kind=real_lsprec),    &
  auto_graup_coeff      = real(auto_graup_coeff_orig, kind=real_lsprec),       &
  inhomog_rate          = real(inhomog_rate_orig, kind=real_lsprec),           &
  inhomog_lim           = real(inhomog_lim_orig, kind=real_lsprec),            &
  r_thresh              = real(r_thresh_orig, kind=real_lsprec),               &
  r_auto                = real(r_auto_orig, kind=real_lsprec),                 &
  n_auto                = real(n_auto_orig, kind=real_lsprec),                 &
  power_droplet_auto    = real(power_droplet_auto_orig, kind=real_lsprec),     &
  power_qcl_auto        = real(power_qcl_auto_orig, kind=real_lsprec),         &
  power_rho_auto        = real(power_rho_auto_orig, kind=real_lsprec),         &
  consts_auto           = real(consts_auto_orig, kind=real_lsprec),            &
  aut_pref              = real(aut_pref_orig, kind=real_lsprec),               &
  aut_nc                = real(aut_nc_orig, kind=real_lsprec)

!mphys_constants_mod
real (kind=real_lsprec), parameter ::                                          &
  acc_pref    = real(acc_pref_orig, kind=real_lsprec),                         &
  acc_qc      = real(acc_qc_orig, kind=real_lsprec),                           &
  acc_qr      = real(acc_qr_orig, kind=real_lsprec),                           &
  mprog_min   = real(mprog_min_orig, kind=real_lsprec),                        &
  mprog_abs   = real(mprog_abs_orig, kind=real_lsprec),                        &
  ntot_land   = real(ntot_land_orig, kind=real_lsprec),                        &
  ntot_sea    = real(ntot_sea_orig, kind=real_lsprec),                         &
  max_as_enh  = real(max_as_enh_orig, kind=real_lsprec)

!conversions_mod
real (kind=real_lsprec), parameter ::                                          &
  zerodegc = real(zerodegc_orig, kind=real_lsprec),                            &
  pi       = real(pi_orig, kind=real_lsprec)

!mphys_ice_mod
real (kind=real_lsprec), parameter ::                                          &
  m0        = real(m0_orig, kind=real_lsprec),                                 &
  t_scaling = real(t_scaling_orig, kind=real_lsprec),                          &
  qcf0      = real(qcf0_orig, kind=real_lsprec),                               &
  t_agg_min = real(t_agg_min_orig, kind=real_lsprec),                          &
  qcfmin    = real(qcfmin_orig, kind=real_lsprec),                             &
  thomo     = real(thomo_orig, kind=real_lsprec)

!water_constants_mod
real (kind=real_lsprec), parameter ::                                          &
  lc = real(lc_orig, kind=real_lsprec),                                        &
  lf = real(lf_orig, kind=real_lsprec)

!planet_constants_mod
real (kind=real_lsprec), parameter ::                                          &
  rv = real(rv_orig, kind=real_lsprec)

!pc2_constants_mod
real (kind=real_lsprec), parameter ::                                          &
  wind_shear_factor = real(wind_shear_factor_orig, kind=real_lsprec)

!mphys_radar_mod
real (kind=real_lsprec), parameter ::                                          &
  kliq        = real(kliq_orig, kind=real_lsprec),                             &
  kice        = real(kice_orig, kind=real_lsprec),                             &
  mm6m3       = real(mm6m3_orig, kind=real_lsprec),                            &
  ref_lim     = real(ref_lim_orig, kind=real_lsprec),                          &
  ref_lim_lin = real(ref_lim_lin_orig, kind=real_lsprec),                      &
  mr_lim      = real(mr_lim_orig, kind=real_lsprec),                           &
  cf_lim      = real(cf_lim_orig, kind=real_lsprec),                           &
  nd_lim      = real(nd_lim_orig, kind=real_lsprec),                           &
  ref_mom     = real(ref_mom_orig, kind=real_lsprec)



!Non-parameter variables
!General atmosphere modules
real (kind=real_lsprec) ::                                                     &
  !planet_constants_mod
  lcrcp, lfrcp, cp, r, repsilon, recip_epsilon, one_minus_epsilon, g,          &
  !timestep_mod
  timestep,                                                                    &
  !cloud_inputs_mod
  ice_width, cff_spread_rate,                                                  &
  !cderived_mod
  delta_lambda, delta_phi,                                                     &
  !fsd_parameters_mod
  f_cons(3), fsd_eff_lam, fsd_eff_phi,                                         &
  !rad_input_mod
  rad_mcica_sigma, two_d_fsd_factor

!Microphysics modules
real (kind=real_lsprec) ::                                                     &
  !mphys_constants_mod
  timestep_mp, cx(200), constp(200), qclmin_rime, area_ratio_prefac,           &
  area_ratio_expn, rho_q_veloc, lam_evap_enh, ec_auto, tnuc,                   &
  !lsp_dif_mod
  apb1, apb2, apb3, apb4, apb5, apb6,                                          &
  !mphys_inputs_mod
  ai, bi, nscalesf, aut_qc

contains

subroutine lsprec_set_reals()

! Purpose:
!Module to manage conversions of all reals USEd below ls_ppnc to lsprec
!precision.

! Method:
! PARAMETERs are copied to PARAMETERs
! Normal REALS are copied to REALs

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: To follow

!All PARAMETERs are dealt with at the module level.
!All others are dealt with at the subroutine level.


!General atmosphere modules- not parameter
use planet_constants_mod, only:                                                &
  lcrcp_orig              => lcrcp,                                            &
  lfrcp_orig              => lfrcp,                                            &
  cp_orig                 => cp,                                               &
  r_orig                  => r,                                                &
  repsilon_orig           => repsilon,                                         &
  recip_epsilon_orig      => recip_epsilon,                                    &
  one_minus_epsilon_orig  => one_minus_epsilon,                                &
  g_orig                  => g

use timestep_mod,         only:                                                &
  timestep_orig => timestep

use cloud_inputs_mod,     only:                                                &
  ice_width_orig        => ice_width,                                          &
  cff_spread_rate_orig  => cff_spread_rate

use cderived_mod,         only:                                                &
  delta_lambda_orig   => delta_lambda,                                         &
  delta_phi_orig      => delta_phi

use fsd_parameters_mod,   only:                                                &
  f_cons_orig         => f_cons,                                               &
  fsd_eff_lam_orig    => fsd_eff_lam,                                          &
  fsd_eff_phi_orig    => fsd_eff_phi

use rad_input_mod,        only:                                                &
  rad_mcica_sigma_orig  => rad_mcica_sigma,                                    &
  two_d_fsd_factor_orig => two_d_fsd_factor

!Microphysics modules- not parameter
use mphys_constants_mod,  only:                                                &
  timestep_mp_orig        => timestep_mp,                                      &
  cx_orig                 => cx,                                               &
  constp_orig             => constp,                                           &
  qclmin_rime_orig        => qclmin_rime,                                      &
  area_ratio_prefac_orig  => area_ratio_prefac,                                &
  area_ratio_expn_orig    => area_ratio_expn,                                  &
  rho_q_veloc_orig        => rho_q_veloc,                                      &
  lam_evap_enh_orig       => lam_evap_enh,                                     &
  ec_auto_orig            => ec_auto,                                          &
  tnuc_orig               => tnuc

use lsp_dif_mod,          only:                                                &
  apb1_orig => apb1,                                                           &
  apb2_orig => apb2,                                                           &
  apb3_orig => apb3,                                                           &
  apb4_orig => apb4,                                                           &
  apb5_orig => apb5,                                                           &
  apb6_orig => apb6

use mphys_inputs_mod,     only:                                                &
  ai_orig             => ai,                                                   &
  bi_orig             => bi,                                                   &
  nscalesf_orig       => nscalesf,                                             &
  aut_qc_orig         => aut_qc

implicit none

!===============================================================================
!Use in reals in lsprec precision, both microphysics related and general atmos
!use lsprec_mod, only:

!- logicals and integers
!===============================================================================

!In general, variables should be declared in the module so they can be used
!elsewhere

!End of header

!General atmosphere modules
  !planet_constants_mod
lcrcp             = real(lcrcp_orig, kind=real_lsprec)
lfrcp             = real(lfrcp_orig, kind=real_lsprec)
cp                = real(cp_orig, kind=real_lsprec)
r                 = real(r_orig, kind=real_lsprec)
repsilon          = real(repsilon_orig, kind=real_lsprec)
recip_epsilon     = real(recip_epsilon_orig, kind=real_lsprec)
one_minus_epsilon = real(one_minus_epsilon_orig, kind=real_lsprec)
g                 = real(g_orig, kind=real_lsprec)

!timestep_mod
timestep = real(timestep_orig, kind=real_lsprec)

!cloud_inputs_mod
ice_width       = real(ice_width_orig, kind=real_lsprec)
cff_spread_rate = real(cff_spread_rate_orig, kind=real_lsprec)

!cderived_mod
delta_lambda  = real(delta_lambda_orig, kind=real_lsprec)
delta_phi     = real(delta_phi_orig, kind=real_lsprec)

!fsd_parameters_mod
f_cons(:)   = real(f_cons_orig(:), kind=real_lsprec)
fsd_eff_lam = real(fsd_eff_lam_orig, kind=real_lsprec)
fsd_eff_phi = real(fsd_eff_phi_orig, kind=real_lsprec)

!rad_input_mod
rad_mcica_sigma   = real(rad_mcica_sigma_orig, kind=real_lsprec)
two_d_fsd_factor  = real(two_d_fsd_factor_orig, kind=real_lsprec)

!Microphysics modules
  !mphys_constants_mod
timestep_mp       = real(timestep_mp_orig, kind=real_lsprec)
cx(:)             = real(cx_orig(:), kind=real_lsprec)
constp(:)         = real(constp_orig(:), kind=real_lsprec)
qclmin_rime       = real(qclmin_rime_orig, kind=real_lsprec)
area_ratio_prefac = real(area_ratio_prefac_orig, kind=real_lsprec)
area_ratio_expn   = real(area_ratio_expn_orig, kind=real_lsprec)
rho_q_veloc       = real(rho_q_veloc_orig, kind=real_lsprec)
lam_evap_enh      = real(lam_evap_enh_orig, kind=real_lsprec)
ec_auto           = real(ec_auto_orig, kind=real_lsprec)
tnuc              = real(tnuc_orig, kind=real_lsprec)

!lsp_dif_mod
apb1  = real(apb1_orig, kind=real_lsprec)
apb2  = real(apb2_orig, kind=real_lsprec)
apb3  = real(apb3_orig, kind=real_lsprec)
apb4  = real(apb4_orig, kind=real_lsprec)
apb5  = real(apb5_orig, kind=real_lsprec)
apb6  = real(apb6_orig, kind=real_lsprec)

!mphys_inputs_mod
ai        = real(ai_orig, kind=real_lsprec)
bi        = real(bi_orig, kind=real_lsprec)
nscalesf  = real(nscalesf_orig, kind=real_lsprec)
aut_qc    = real(aut_qc_orig, kind=real_lsprec)

end subroutine lsprec_set_reals
end module lsprec_mod
