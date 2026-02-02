! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module comorph_diags_um_mod

use um_types, only: real_umphys

implicit none

contains



!----------------------------------------------------------------
! Subroutine to set the diagnostic requests in CoMorph's
! diagnostics system, based on the switches set for the
! convection diagnostics in the UM.
!----------------------------------------------------------------
subroutine comorph_diags_um_reqs( comorph_diags,                               &
                                  pressure_incr_env )

use comorph_diags_type_mod, only: comorph_diags_type

use cv_diagnostic_array_mod, only:                                             &
      up_flux_half, dwn_flux_half,                                             &
      entrain_up, detrain_up, entrain_dwn, detrain_dwn,                        &
      cape_out, par_radius_up, par_radius_dwn,                                 &
      par_meanup_dtv, par_meandn_dtv ,par_meanup_rhl, par_meandn_rhl,          &
      par_meanup_q, par_meandn_q, par_meanup_qcl, par_meandn_qcl,              &
      par_meanup_qcf, par_meandn_qcf, par_meanup_qrain, par_meandn_qrain,      &
      par_meanup_qgr, par_meandn_qgr, par_meanup_qsnow, par_meandn_qsnow,      &
      par_meanup_cfl, par_meandn_cfl, par_meanup_cff, par_meandn_cff,          &
      par_meanup_cfb, par_meandn_cfb,                                          &
      par_meanup_u, par_meandn_u, par_meanup_v, par_meandn_v,                  &
      par_meanup_w, par_meandn_w, par_meanup_t, par_meandn_t,                  &
      par_coreup_dtv, par_coredn_dtv ,par_coreup_rhl, par_coredn_rhl,          &
      dry_fraction, liq_fraction, mix_fraction, icr_fraction,                  &
      dry_temp, liq_temp, mix_temp, icr_temp ,dry_q_vap, liq_q_vap,            &
      mix_q_vap, icr_q_vap ,dry_rhl, liq_rhl, mix_rhl, icr_rhl,                &
      turb_pert_u, turb_pert_v, turb_pert_w, turb_pert_t, turb_pert_q,         &
      turb_radius

use cv_stash_flg_mod, only:                                                    &
      flg_up_flx_half, flg_dwn_flx_half,                                       &
      flg_entr_up, flg_detr_up, flg_entr_dwn, flg_detr_dwn,                    &
      flg_par_radius_up, flg_par_radius_dwn,                                   &
      flg_par_meanup_dtv, flg_par_meandn_dtv,flg_par_meanup_rhl,               &
      flg_par_meandn_rhl, flg_par_coreup_dtv, flg_par_coredn_dtv,              &
      flg_par_meanup_q, flg_par_meandn_q, flg_par_meanup_qcl,                  &
      flg_par_meandn_qcl, flg_par_meanup_qcf, flg_par_meandn_qcf,              &
      flg_par_meanup_qrain, flg_par_meandn_qrain, flg_par_meanup_qgr,          &
      flg_par_meandn_qgr, flg_par_meanup_qsnow, flg_par_meandn_qsnow,          &
      flg_par_meanup_cfl, flg_par_meandn_cfl, flg_par_meanup_cff,              &
      flg_par_meandn_cff, flg_par_meanup_cfb, flg_par_meandn_cfb,              &
      flg_par_meanup_u, flg_par_meandn_u, flg_par_meanup_v, flg_par_meandn_v,  &
      flg_par_meanup_w, flg_par_meandn_w, flg_par_meanup_t, flg_par_meandn_t,  &
      flg_par_coreup_rhl, flg_par_coredn_rhl,                                  &
      flg_dry_fraction, flg_liq_fraction, flg_mix_fraction, flg_icr_fraction,  &
      flg_dry_temp, flg_liq_temp, flg_mix_temp, flg_icr_temp,                  &
      flg_dry_q_vap, flg_liq_q_vap, flg_mix_q_vap, flg_icr_q_vap,              &
      flg_dry_rhl, flg_liq_rhl, flg_mix_rhl, flg_icr_rhl,                      &
      flg_turb_pert_u, flg_turb_pert_v, flg_turb_pert_w, flg_turb_pert_t,      &
      flg_turb_pert_q, flg_turb_radius


use atm_fields_bounds_mod, only: tdims

use cloud_inputs_mod, only: i_cld_vn, l_pc2_homog_conv_pressure
use pc2_constants_mod, only: i_cld_pc2

implicit none

! Structure containing flags and meta-data for all the CoMorph
! diagnostics, and pointers to the output arrays
type(comorph_diags_type), intent(in out) :: comorph_diags

! Pressure change following parcels in the environment
! (needed to calculate the PC2 cloud response to convective
!  subsidence)
real(kind=real_umphys), target, intent(in out) :: pressure_incr_env            &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)
! This is handled as a diagnostic inside CoMorph,
! but is needed elsewhere in the model.


! Note: the following code asigns pointers from the comorph_diags
! structure for each requested diagnostic.  The pointers must
! be nullified again later, in the routine comorph_diags_um_null
! below.
! Therefore, if you add extra diag requests to this routine,
! you MUST also add the corresponding nullifications in
! the comorph_diags_um_null routine.

! For each requested diagnostic, we need to:
! a) Set its request flag for the appropriate domain to true.
! b) Point the pointer for the appropriate domain at the
!    output diagnostic array.

! Updraft and downdraft mass-fluxes:
! CoMorph calculates the mass-fluxes on rho-levels.
! If they have been requested on theta-levels, we'll need to
! interpolate them onto theta-levels afterwards.
! The logic in cv_stash_flg_mod sets the flags for mass-fluxes
! on rho-levels to true even if they are only requested on
! theta-levels, to ensure they are available to be interpolated.
if ( flg_up_flx_half ) then
  comorph_diags % updraft % par % massflux_d                                   &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % massflux_d                                   &
                                % field_3d => up_flux_half
end if
if ( flg_dwn_flx_half ) then
  comorph_diags % dndraft % par % massflux_d                                   &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % massflux_d                                   &
                                % field_3d => dwn_flux_half
end if

! Set requests and assign pointers for entrainment and
! detrainment diagnostics
if ( flg_entr_up ) then
  comorph_diags % updraft % plume_model % ent_mass_d                           &
                                % request % x_y_z = .true.
  comorph_diags % updraft % plume_model % ent_mass_d                           &
                                % field_3d => entrain_up
end if
if ( flg_detr_up ) then
  comorph_diags % updraft % plume_model % det_mass_d                           &
                                % request % x_y_z = .true.
  comorph_diags % updraft % plume_model % det_mass_d                           &
                                % field_3d => detrain_up
end if

if ( flg_entr_dwn ) then
  comorph_diags % dndraft % plume_model % ent_mass_d                           &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % plume_model % ent_mass_d                           &
                                % field_3d => entrain_dwn
end if
if ( flg_detr_dwn ) then
  comorph_diags % dndraft % plume_model % det_mass_d                           &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % plume_model % det_mass_d                           &
                                % field_3d => detrain_dwn
end if

! Average parcel radius
if ( flg_par_radius_up ) then
  comorph_diags % updraft % par % radius                                       &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % radius                                       &
                                % field_3d => par_radius_up
end if
if ( flg_par_radius_dwn ) then
  comorph_diags % dndraft % par % radius                                       &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % radius                                       &
                                % field_3d => par_radius_dwn
end if

! Parcel mean properties
if ( flg_par_meanup_dtv ) then
  comorph_diags % updraft % par % mean % virt_temp_excess                      &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % virt_temp_excess                      &
                                % field_3d => par_meanup_dtv
end if
if ( flg_par_meandn_dtv ) then
  comorph_diags % dndraft % par % mean % virt_temp_excess                      &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % virt_temp_excess                      &
                                % field_3d => par_meandn_dtv
end if
if ( flg_par_meanup_rhl ) then
  comorph_diags % updraft % par % mean % rel_hum_liq                           &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % rel_hum_liq                           &
                                % field_3d => par_meanup_rhl
end if
if ( flg_par_meandn_rhl ) then
  comorph_diags % dndraft % par % mean % rel_hum_liq                           &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % rel_hum_liq                           &
                                % field_3d => par_meandn_rhl
end if

if ( flg_par_meanup_q ) then
  comorph_diags % updraft % par % mean % q_vap                                 &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % q_vap                                 &
                                % field_3d => par_meanup_q
end if
if ( flg_par_meandn_q ) then
  comorph_diags % dndraft % par % mean % q_vap                                 &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % q_vap                                 &
                                % field_3d => par_meandn_q
end if

if ( flg_par_meanup_qcl ) then
  comorph_diags % updraft % par % mean % q_cl                                  &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % q_cl                                  &
                                % field_3d => par_meanup_qcl
end if
if ( flg_par_meandn_qcl ) then
  comorph_diags % dndraft % par % mean % q_cl                                  &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % q_cl                                  &
                                % field_3d => par_meandn_qcl
end if

if ( flg_par_meanup_qcf ) then
  comorph_diags % updraft % par % mean % q_cf                                  &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % q_cf                                  &
                                % field_3d => par_meanup_qcf
end if
if ( flg_par_meandn_qcf ) then
  comorph_diags % dndraft % par % mean % q_cf                                  &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % q_cf                                  &
                                % field_3d => par_meandn_qcf
end if

if ( flg_par_meanup_qrain ) then
  comorph_diags % updraft % par % mean % q_rain                                &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % q_rain                                &
                                % field_3d => par_meanup_qrain
end if
if ( flg_par_meandn_qrain ) then
  comorph_diags % dndraft % par % mean % q_rain                                &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % q_rain                                &
                                % field_3d => par_meandn_qrain
end if

if ( flg_par_meanup_qgr ) then
  comorph_diags % updraft % par % mean % q_graup                               &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % q_graup                               &
                                % field_3d => par_meanup_qgr
end if
if ( flg_par_meandn_qgr ) then
  comorph_diags % dndraft % par % mean % q_graup                               &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % q_graup                               &
                                % field_3d => par_meandn_qgr
end if

if ( flg_par_meanup_qsnow ) then
  comorph_diags % updraft % par % mean % q_snow                                &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % q_snow                                &
                                % field_3d => par_meanup_qsnow
end if
if ( flg_par_meandn_qsnow ) then
  comorph_diags % dndraft % par % mean % q_snow                                &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % q_snow                                &
                                % field_3d => par_meandn_qsnow
end if

if ( flg_par_meanup_cfl ) then
  comorph_diags % updraft % par % mean % cf_liq                                &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % cf_liq                                &
                                % field_3d => par_meanup_cfl
end if
if ( flg_par_meandn_cfl ) then
  comorph_diags % dndraft % par % mean % cf_liq                                &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % cf_liq                                &
                                % field_3d => par_meandn_cfl
end if

if ( flg_par_meanup_cff ) then
  comorph_diags % updraft % par % mean % cf_ice                                &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % cf_ice                                &
                                % field_3d => par_meanup_cff
end if
if ( flg_par_meandn_cff ) then
  comorph_diags % dndraft % par % mean % cf_ice                                &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % cf_ice                                &
                                % field_3d => par_meandn_cff
end if

if ( flg_par_meanup_cfb ) then
  comorph_diags % updraft % par % mean % cf_bulk                               &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % cf_bulk                               &
                                % field_3d => par_meanup_cfb
end if
if ( flg_par_meandn_cfb ) then
  comorph_diags % dndraft % par % mean % cf_bulk                               &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % cf_bulk                               &
                                % field_3d => par_meandn_cfb
end if

if ( flg_par_meanup_u ) then
  comorph_diags % updraft % par % mean % wind_u                                &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % wind_u                                &
                                % field_3d => par_meanup_u
end if
if ( flg_par_meandn_u ) then
  comorph_diags % dndraft % par % mean % wind_u                                &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % wind_u                                &
                                % field_3d => par_meandn_u
end if

if ( flg_par_meanup_v ) then
  comorph_diags % updraft % par % mean % wind_v                                &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % wind_v                                &
                                % field_3d => par_meanup_v
end if
if ( flg_par_meandn_v ) then
  comorph_diags % dndraft % par % mean % wind_v                                &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % wind_v                                &
                                % field_3d => par_meandn_v
end if

if ( flg_par_meanup_w ) then
  comorph_diags % updraft % par % mean % wind_w                                &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % wind_w                                &
                                % field_3d => par_meanup_w
end if
if ( flg_par_meandn_w ) then
  comorph_diags % dndraft % par % mean % wind_w                                &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % wind_w                                &
                                % field_3d => par_meandn_w
end if

if ( flg_par_meanup_t ) then
  comorph_diags % updraft % par % mean % temperature                           &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % mean % temperature                           &
                                % field_3d => par_meanup_t
end if
if ( flg_par_meandn_t ) then
  comorph_diags % dndraft % par % mean % temperature                           &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % mean % temperature                           &
                                % field_3d => par_meandn_t
end if


! Parcel core properties
if ( flg_par_coreup_dtv ) then
  comorph_diags % updraft % par % core % virt_temp_excess                      &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % core % virt_temp_excess                      &
                                % field_3d => par_coreup_dtv
end if
if ( flg_par_coredn_dtv ) then
  comorph_diags % dndraft % par % core % virt_temp_excess                      &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % core % virt_temp_excess                      &
                                % field_3d => par_coredn_dtv
end if
if ( flg_par_coreup_rhl ) then
  comorph_diags % updraft % par % core % rel_hum_liq                           &
                                % request % x_y_z = .true.
  comorph_diags % updraft % par % core % rel_hum_liq                           &
                                % field_3d => par_coreup_rhl
end if
if ( flg_par_coredn_rhl ) then
  comorph_diags % dndraft % par % core % rel_hum_liq                           &
                                % request % x_y_z = .true.
  comorph_diags % dndraft % par % core % rel_hum_liq                           &
                                % field_3d => par_coredn_rhl
end if

! Turbulent perturbations
if ( flg_turb_pert_u ) then
  comorph_diags % turb_fields_pert % wind_u                                    &
                                % request % x_y_z = .true.
  comorph_diags % turb_fields_pert % wind_u                                    &
                                % field_3d => turb_pert_u
end if
if ( flg_turb_pert_v ) then
  comorph_diags % turb_fields_pert % wind_v                                    &
                                % request % x_y_z = .true.
  comorph_diags % turb_fields_pert % wind_v                                    &
                                % field_3d => turb_pert_v
end if
if ( flg_turb_pert_w ) then
  comorph_diags % turb_fields_pert % wind_w                                    &
                                % request % x_y_z = .true.
  comorph_diags % turb_fields_pert % wind_w                                    &
                                % field_3d => turb_pert_w
end if
if ( flg_turb_pert_t ) then
  comorph_diags % turb_fields_pert % temperature                               &
                                % request % x_y_z = .true.
  comorph_diags % turb_fields_pert % temperature                               &
                                % field_3d => turb_pert_t
end if
if ( flg_turb_pert_q ) then
  comorph_diags % turb_fields_pert % q_vap                                     &
                                % request % x_y_z = .true.
  comorph_diags % turb_fields_pert % q_vap                                     &
                                % field_3d => turb_pert_q
end if
if ( flg_turb_radius ) then
  comorph_diags % turb_radius % request % x_y_z = .true.
  comorph_diags % turb_radius % field_3d => turb_radius
end if


! Comorph division of gridbox into sub-regions
if ( flg_dry_fraction ) then
  comorph_diags % genesis_diags % subregion_diags(1) % frac                    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(1) % frac                    &
                                % field_3d => dry_fraction
end if
if ( flg_liq_fraction ) then
  comorph_diags % genesis_diags % subregion_diags(2) % frac                    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(2) % frac                    &
                                % field_3d => liq_fraction
end if
if ( flg_mix_fraction ) then
  comorph_diags % genesis_diags % subregion_diags(3) % frac                    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(3) % frac                    &
                                % field_3d => mix_fraction
end if
if ( flg_icr_fraction ) then
  comorph_diags % genesis_diags % subregion_diags(4) % frac                    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(4) % frac                    &
                                % field_3d => icr_fraction
end if
if ( flg_dry_temp ) then
  comorph_diags % genesis_diags % subregion_diags(1) % fields % temperature    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(1) % fields % temperature    &
                                % field_3d => dry_temp
end if
if ( flg_liq_temp ) then
  comorph_diags % genesis_diags % subregion_diags(2) % fields % temperature    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(2) % fields % temperature    &
                                % field_3d => liq_temp
end if
if ( flg_mix_temp ) then
  comorph_diags % genesis_diags % subregion_diags(3) % fields % temperature    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(3) % fields % temperature    &
                                % field_3d => mix_temp
end if
if ( flg_icr_temp ) then
  comorph_diags % genesis_diags % subregion_diags(4) % fields % temperature    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(4) % fields % temperature    &
                                % field_3d => icr_temp
end if
if ( flg_dry_q_vap ) then
  comorph_diags % genesis_diags % subregion_diags(1) % fields % q_vap          &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(1) % fields % q_vap          &
                                % field_3d => dry_q_vap
end if
if ( flg_liq_q_vap ) then
  comorph_diags % genesis_diags % subregion_diags(2) % fields % q_vap          &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(2) % fields % q_vap          &
                                % field_3d => liq_q_vap
end if
if ( flg_mix_q_vap ) then
  comorph_diags % genesis_diags % subregion_diags(3) % fields % q_vap          &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(3) % fields % q_vap          &
                                % field_3d => mix_q_vap
end if
if ( flg_icr_q_vap ) then
  comorph_diags % genesis_diags % subregion_diags(4) % fields % q_vap          &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(4) % fields % q_vap          &
                                % field_3d => icr_q_vap
end if
if ( flg_dry_rhl ) then
  comorph_diags % genesis_diags % subregion_diags(1) % fields % rel_hum_liq    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(1) % fields % rel_hum_liq    &
                                % field_3d => dry_rhl
end if
if ( flg_liq_rhl ) then
  comorph_diags % genesis_diags % subregion_diags(2) % fields % rel_hum_liq    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(2) % fields % rel_hum_liq    &
                                % field_3d => liq_rhl
end if
if ( flg_mix_rhl ) then
  comorph_diags % genesis_diags % subregion_diags(3) % fields % rel_hum_liq    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(3) % fields % rel_hum_liq    &
                                % field_3d => mix_rhl
end if
if ( flg_icr_rhl ) then
  comorph_diags % genesis_diags % subregion_diags(4) % fields % rel_hum_liq    &
                                % request % x_y_z = .true.
  comorph_diags % genesis_diags % subregion_diags(4) % fields % rel_hum_liq    &
                                % field_3d => icr_rhl
end if


! Always request the CAPE diagnostic, as it can be used by
! other parts of the UM, regardless of whether it is
! requested as an output diagostic.
! Request the mass-flux-weighted CAPE, not straight CAPE,
! since CAPE itself can be noisy / not representative
! (e.g. when most of the mass-flux is shallow, but there
!  is a tiny mass-flux going deep with high CAPE).
comorph_diags % updraft_diags_2d % mfw_cape                                    &
                                % request % x_y = .true.
comorph_diags % updraft_diags_2d % mfw_cape                                    &
                                % field_2d => cape_out

! Set request and assign pointer for environment pressure
! change due to convective subsidence, if needed.
if ( i_cld_vn == i_cld_pc2 .and. l_pc2_homog_conv_pressure ) then
  comorph_diags % pressure_incr_env                                            &
                                % request % x_y_z = .true.
  comorph_diags % pressure_incr_env                                            &
                                % field_3d => pressure_incr_env
end if


return
end subroutine comorph_diags_um_reqs



!----------------------------------------------------------------
! Subroutine to nullify the pointers used by the CoMorph
! diagnostics system, after the call to CoMorph is complete.
!----------------------------------------------------------------
! Some of these pointers will be pointing at allocatable arrays;
! if any of these get deallocated while the pointers are still
! pointing at them, it may cause problems, so need to nullify
! the pointers first.
subroutine comorph_diags_um_null( comorph_diags )

use comorph_diags_type_mod, only: comorph_diags_type

use cv_stash_flg_mod, only:                                                    &
      flg_up_flx_half, flg_dwn_flx_half,                                       &
      flg_entr_up, flg_detr_up, flg_entr_dwn, flg_detr_dwn,                    &
      flg_par_radius_up, flg_par_radius_dwn,                                   &
      flg_par_meanup_dtv, flg_par_meandn_dtv,flg_par_meanup_rhl,               &
      flg_par_meanup_q, flg_par_meandn_q, flg_par_meanup_qcl,                  &
      flg_par_meandn_qcl, flg_par_meanup_qcf, flg_par_meandn_qcf,              &
      flg_par_meanup_qrain, flg_par_meandn_qrain, flg_par_meanup_qgr,          &
      flg_par_meandn_qgr, flg_par_meanup_qsnow, flg_par_meandn_qsnow,          &
      flg_par_meanup_cfl, flg_par_meandn_cfl, flg_par_meanup_cff,              &
      flg_par_meandn_cff, flg_par_meanup_cfb, flg_par_meandn_cfb,              &
      flg_par_meanup_u, flg_par_meandn_u, flg_par_meanup_v, flg_par_meandn_v,  &
      flg_par_meanup_w, flg_par_meandn_w, flg_par_meanup_t, flg_par_meandn_t,  &
      flg_par_meandn_rhl, flg_par_coreup_dtv, flg_par_coredn_dtv,              &
      flg_par_coreup_rhl, flg_par_coredn_rhl,                                  &
      flg_dry_fraction, flg_liq_fraction, flg_mix_fraction, flg_icr_fraction,  &
      flg_dry_temp, flg_liq_temp, flg_mix_temp, flg_icr_temp,                  &
      flg_dry_q_vap, flg_liq_q_vap, flg_mix_q_vap, flg_icr_q_vap,              &
      flg_dry_rhl, flg_liq_rhl, flg_mix_rhl, flg_icr_rhl,                      &
      flg_turb_pert_u, flg_turb_pert_v, flg_turb_pert_w, flg_turb_pert_t,      &
      flg_turb_pert_q, flg_turb_radius

use cloud_inputs_mod, only: i_cld_vn
use pc2_constants_mod, only: i_cld_pc2

implicit none

! Structure containing flags and meta-data for all the CoMorph
! diagnostics, and pointers to the output arrays
type(comorph_diags_type), intent(in out) :: comorph_diags


! Nullify pointers for all used diagnostics.
! This mirrors the structure of comorph_diags_um_reqs above.

if ( flg_up_flx_half )   comorph_diags % updraft % par                         &
                         % massflux_d % field_3d => null()
if ( flg_dwn_flx_half )  comorph_diags % dndraft % par                         &
                         % massflux_d % field_3d => null()

if ( flg_entr_up )   comorph_diags % updraft % plume_model                     &
                     % ent_mass_d % field_3d => null()
if ( flg_detr_up )   comorph_diags % updraft % plume_model                     &
                     % det_mass_d % field_3d => null()
if ( flg_entr_dwn )  comorph_diags % dndraft % plume_model                     &
                     % ent_mass_d % field_3d => null()
if ( flg_detr_dwn )  comorph_diags % dndraft % plume_model                     &
                     % det_mass_d % field_3d => null()

if ( flg_par_radius_up )  comorph_diags % updraft % par                        &
                     % radius % field_3d => null()
if ( flg_par_radius_dwn )  comorph_diags % dndraft % par                       &
                     % radius % field_3d => null()

if ( flg_par_meanup_dtv )  comorph_diags % updraft % par % mean                &
                     % virt_temp_excess % field_3d => null()
if ( flg_par_meandn_dtv )  comorph_diags % dndraft % par % mean                &
                     % virt_temp_excess % field_3d => null()
if ( flg_par_meanup_rhl )  comorph_diags % updraft % par % mean                &
                     % rel_hum_liq % field_3d => null()
if ( flg_par_meandn_rhl )  comorph_diags % dndraft % par % mean                &
                     % rel_hum_liq % field_3d => null()
if ( flg_par_meanup_q )  comorph_diags % updraft % par % mean                  &
                     % q_vap % field_3d => null()
if ( flg_par_meandn_q )  comorph_diags % dndraft % par % mean                  &
                     % q_vap % field_3d => null()
if ( flg_par_meanup_qcl )  comorph_diags % updraft % par % mean                &
                     % q_cl % field_3d => null()
if ( flg_par_meandn_qcl )  comorph_diags % dndraft % par % mean                &
                     % q_cl % field_3d => null()
if ( flg_par_meanup_qcf )  comorph_diags % updraft % par % mean                &
                     % q_cf % field_3d => null()
if ( flg_par_meandn_qcf )  comorph_diags % dndraft % par % mean                &
                     % q_cf % field_3d => null()
if ( flg_par_meanup_qrain )  comorph_diags % updraft % par % mean              &
                     % q_rain % field_3d => null()
if ( flg_par_meandn_qrain )  comorph_diags % dndraft % par % mean              &
                     % q_rain % field_3d => null()
if ( flg_par_meanup_qgr )  comorph_diags % updraft % par % mean                &
                     % q_graup % field_3d => null()
if ( flg_par_meandn_qgr )  comorph_diags % dndraft % par % mean                &
                     % q_graup % field_3d => null()
if ( flg_par_meanup_qsnow )  comorph_diags % updraft % par % mean              &
                     % q_snow % field_3d => null()
if ( flg_par_meandn_qsnow )  comorph_diags % dndraft % par % mean              &
                     % q_snow % field_3d => null()
if ( flg_par_meanup_cfl )  comorph_diags % updraft % par % mean                &
                     % cf_liq % field_3d => null()
if ( flg_par_meandn_cfl )  comorph_diags % dndraft % par % mean                &
                     % cf_liq % field_3d => null()
if ( flg_par_meanup_cff )  comorph_diags % updraft % par % mean                &
                     % cf_ice % field_3d => null()
if ( flg_par_meandn_cff )  comorph_diags % dndraft % par % mean                &
                     % cf_ice % field_3d => null()
if ( flg_par_meanup_cfb )  comorph_diags % updraft % par % mean                &
                     % cf_bulk % field_3d => null()
if ( flg_par_meandn_cfb )  comorph_diags % dndraft % par % mean                &
                     % cf_bulk % field_3d => null()

if ( flg_par_meanup_u )  comorph_diags % updraft % par % mean                  &
                     % wind_u % field_3d => null()
if ( flg_par_meandn_u )  comorph_diags % dndraft % par % mean                  &
                     % wind_u % field_3d => null()
if ( flg_par_meanup_v )  comorph_diags % updraft % par % mean                  &
                     % wind_v % field_3d => null()
if ( flg_par_meandn_v )  comorph_diags % dndraft % par % mean                  &
                     % wind_v % field_3d => null()
if ( flg_par_meanup_w )  comorph_diags % updraft % par % mean                  &
                     % wind_w % field_3d => null()
if ( flg_par_meandn_w )  comorph_diags % dndraft % par % mean                  &
                     % wind_w % field_3d => null()
if ( flg_par_meanup_t )  comorph_diags % updraft % par % mean                  &
                     % temperature % field_3d => null()
if ( flg_par_meandn_t )  comorph_diags % dndraft % par % mean                  &
                     % temperature % field_3d => null()

if ( flg_par_coreup_dtv )  comorph_diags % updraft % par % core                &
                     % virt_temp_excess % field_3d => null()
if ( flg_par_coredn_dtv )  comorph_diags % dndraft % par % core                &
                     % virt_temp_excess % field_3d => null()
if ( flg_par_coreup_rhl )  comorph_diags % updraft % par % core                &
                     % rel_hum_liq % field_3d => null()
if ( flg_par_coredn_rhl )  comorph_diags % dndraft % par % core                &
                     % rel_hum_liq % field_3d => null()

if ( flg_turb_pert_u )  comorph_diags % turb_fields_pert                       &
                     % wind_u % field_3d => null()
if ( flg_turb_pert_v )  comorph_diags % turb_fields_pert                       &
                     % wind_v % field_3d => null()
if ( flg_turb_pert_w )  comorph_diags % turb_fields_pert                       &
                     % wind_w % field_3d => null()
if ( flg_turb_pert_t )  comorph_diags % turb_fields_pert                       &
                     % temperature % field_3d => null()
if ( flg_turb_pert_q )  comorph_diags % turb_fields_pert                       &
                     % q_vap % field_3d => null()
if ( flg_turb_radius )  comorph_diags % turb_radius % field_3d => null()

if ( flg_dry_fraction )  comorph_diags % genesis_diags % subregion_diags(1)    &
                     % frac % field_3d => null()
if ( flg_liq_fraction )  comorph_diags % genesis_diags % subregion_diags(2)    &
                     % frac % field_3d => null()
if ( flg_mix_fraction )  comorph_diags % genesis_diags % subregion_diags(3)    &
                     % frac % field_3d => null()
if ( flg_icr_fraction )  comorph_diags % genesis_diags % subregion_diags(4)    &
                     % frac % field_3d => null()
if ( flg_dry_temp )  comorph_diags % genesis_diags % subregion_diags(1)        &
                     % fields % temperature % field_3d => null()
if ( flg_liq_temp )  comorph_diags % genesis_diags % subregion_diags(2)        &
                     % fields % temperature % field_3d => null()
if ( flg_mix_temp )  comorph_diags % genesis_diags % subregion_diags(3)        &
                     % fields % temperature % field_3d => null()
if ( flg_icr_temp )  comorph_diags % genesis_diags % subregion_diags(4)        &
                     % fields % temperature % field_3d => null()
if ( flg_dry_q_vap )  comorph_diags % genesis_diags % subregion_diags(1)       &
                     % fields % q_vap % field_3d => null()
if ( flg_liq_q_vap )  comorph_diags % genesis_diags % subregion_diags(2)       &
                     % fields % q_vap % field_3d => null()
if ( flg_mix_q_vap )  comorph_diags % genesis_diags % subregion_diags(3)       &
                     % fields % q_vap % field_3d => null()
if ( flg_icr_q_vap )  comorph_diags % genesis_diags % subregion_diags(4)       &
                     % fields % q_vap % field_3d => null()
if ( flg_dry_rhl )  comorph_diags % genesis_diags % subregion_diags(1)         &
                     % fields % rel_hum_liq % field_3d => null()
if ( flg_liq_rhl )  comorph_diags % genesis_diags % subregion_diags(2)         &
                     % fields % rel_hum_liq % field_3d => null()
if ( flg_mix_rhl )  comorph_diags % genesis_diags % subregion_diags(3)         &
                     % fields % rel_hum_liq % field_3d => null()
if ( flg_icr_rhl )  comorph_diags % genesis_diags % subregion_diags(4)         &
                     % fields % rel_hum_liq % field_3d => null()


comorph_diags % updraft_diags_2d % mfw_cape % field_2d => null()

if ( i_cld_vn == i_cld_pc2 ) comorph_diags % pressure_incr_env                 &
                                  % field_3d => null()

return
end subroutine comorph_diags_um_null



!----------------------------------------------------------------
! Subroutine to process the diagnostic fields output by the
! CoMorph convection scheme
!----------------------------------------------------------------
subroutine comorph_diags_um_proc( n_conv_levels, z_rho, z_theta, lcbase )

use cv_diagnostic_array_mod, only:                                             &
      up_flux, up_flux_half, dwn_flux, dwn_flux_half,                          &
      entrain_up, detrain_up, entrain_dwn, detrain_dwn, freq_up, freq_dwn
use cv_stash_flg_mod, only:                                                    &
      flg_up_flx, flg_up_flx_half, flg_dwn_flx, flg_dwn_flx_half,              &
      flg_entr_up, flg_detr_up, flg_entr_dwn, flg_detr_dwn,                    &
      flg_freq_up, flg_freq_dwn

use atm_fields_bounds_mod, only: pdims, tdims
use nlsizes_namelist_mod, only: row_length, rows, model_levels
use planet_constants_mod, only: g

implicit none

! Highest model-level where convection is allowed
integer, intent(in) :: n_conv_levels

! Model-level heights above surface
real(kind=real_umphys), intent(in) :: z_rho   ( row_length, rows, model_levels )
real(kind=real_umphys), intent(in) :: z_theta ( row_length, rows, model_levels )

! Model-level of lowest convective cloud-base
integer, intent(in) :: lcbase ( row_length, rows )


! Weight used for vertical interpolation
real(kind=real_umphys) :: interp

! Loop counters
integer :: i, j, k

character(len=*), parameter :: RoutineName                                     &
                               = "COMORPH_DIAGS_UM_PROC"

! Mass-fluxes and entrainment / detrainment rates need to be
! converted from kg m-2 s-1 (as calculated in CoMorph) to
! Pa s-1 (the units stated on the output diagsnostics in both
! STASH and the SCM diags system).  This just entails scaling
! by the acceleration due to gravity, g.
! Note that there is a slight inconsistency here, as the
! CoMorph mass-fluxes are in kg m-2 s-1 of dry-mass only,
! whereas the UM diagnostics are presumably meant to be fluxes
! of total-mass.
! However, the entrainment and detrainment do not form a closed
! budget of total-mass, since they neglect the mass transfered
! by precipitation.  This issue is just swept under the
! carpet in the 6A convection scheme.  For the sake of outputting
! fluxes which do form a closed mass budget, we keep these
! diagnostics in units of dry-mass.

! Another inconsistency is that comorph outputs mass-flux diagnostics
! which include the dry convection in the sub-cloud layer, whereas
! the rest of the UM system expects the updraft mass-flux fields from
! convection to only include the in-cloud mass-flux.
! We therefore need to reset the mass-fluxes to zero below cloud-base.

!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k, interp )                        &
!$OMP SHARED( pdims, tdims, n_conv_levels, g, lcbase, z_rho, z_theta,          &
!$OMP         flg_freq_up, freq_up, flg_freq_dwn, freq_dwn,                    &
!$OMP         flg_up_flx_half, up_flux_half, flg_up_flx, up_flux,              &
!$OMP         flg_dwn_flx_half, dwn_flux_half, flg_dwn_flx, dwn_flux,          &
!$OMP         flg_entr_up, entrain_up, flg_detr_up, detrain_up,                &
!$OMP         flg_entr_dwn, entrain_dwn, flg_detr_dwn, detrain_dwn )

! If total updraft mass-flux requested
! (on theta-levels or rho-levels)
if ( flg_up_flx_half ) then

  if (flg_freq_up) then
    ! Using CoMorph mass flux before modified so know complete set of levels
    ! with non-zero values
!$OMP DO SCHEDULE(STATIC)
    do k = 1, n_conv_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          if ( up_flux_half(i,j,k) > 0.0) then
            freq_up(i,j,k) = 1.0
          else
            freq_up(i,j,k) = 0.0
          end if
        end do
      end do
    end do
!$OMP END DO
  end if

  ! Scale output flux on rho-levels by g to convert to Pa s-1
!$OMP DO SCHEDULE(STATIC)
  do k = 1, n_conv_levels + 1
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        up_flux_half(i,j,k) = up_flux_half(i,j,k) * g
        if ( k<=lcbase(i,j) .or. lcbase(i,j)==0 )  up_flux_half(i,j,k) = 0.0
      end do
    end do
  end do
!$OMP END DO

  ! If flux requested on theta-levels
  if ( flg_up_flx ) then
    ! Interpolate flux onto theta-levels
!$OMP DO SCHEDULE(STATIC)
    do k = 1, n_conv_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          interp = ( z_theta(i,j,k) - z_rho(i,j,k) )                           &
                 / ( z_rho(i,j,k+1) - z_rho(i,j,k) )
          up_flux(i,j,k) = (1.0-interp) * up_flux_half(i,j,k)                  &
                         +      interp  * up_flux_half(i,j,k+1)
          if ( k<=lcbase(i,j) .or. lcbase(i,j)==0 )  up_flux(i,j,k) = 0.0
        end do
      end do
    end do
!$OMP END DO NOWAIT
  end if

end if  ! ( flg_up_flx_half )


! If total downdraft mass-flux requested
! (on theta-levels or rho-levels)
if ( flg_dwn_flx_half ) then

  if (flg_freq_dwn) then
    ! Using CoMorph mass flux before modified so know complete set of levels
    ! with non-zero values
!$OMP DO SCHEDULE(STATIC)
    do k = 1, n_conv_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          if ( dwn_flux_half(i,j,k) > 0.0) then
            freq_dwn(i,j,k) = 1.0
          else
            freq_dwn(i,j,k) = 0.0
          end if
        end do
      end do
    end do
!$OMP END DO
  end if

  ! Scale output flux on rho-levels by g to convert to Pa s-1
!$OMP DO SCHEDULE(STATIC)
  do k = 1, n_conv_levels + 1
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        dwn_flux_half(i,j,k) = dwn_flux_half(i,j,k) * g
      end do
    end do
  end do
!$OMP END DO

  ! If flux requested on theta-levels
  if ( flg_dwn_flx ) then
    ! Interpolate flux onto theta-levels
!$OMP DO SCHEDULE(STATIC)
    do k = 1, n_conv_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          interp = ( z_theta(i,j,k) - z_rho(i,j,k) )                           &
                 / ( z_rho(i,j,k+1) - z_rho(i,j,k) )
          dwn_flux(i,j,k) = (1.0-interp) * dwn_flux_half(i,j,k)                &
                          +      interp  * dwn_flux_half(i,j,k+1)
        end do
      end do
    end do
!$OMP END DO NOWAIT
  end if

end if  ! ( flg_dwn_flx_half )


! Convert entrainment and detrainment diagnostics to Pa s-1
! if they are requested
if ( flg_entr_up ) then
!$OMP DO SCHEDULE(STATIC)
  do k = 1, n_conv_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        entrain_up(i,j,k) = entrain_up(i,j,k) * g
        if ( k<=lcbase(i,j) .or. lcbase(i,j)==0 )  entrain_up(i,j,k) = 0.0
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if
if ( flg_detr_up ) then
!$OMP DO SCHEDULE(STATIC)
  do k = 1, n_conv_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        detrain_up(i,j,k) = detrain_up(i,j,k) * g
        if ( k<=lcbase(i,j) .or. lcbase(i,j)==0 )  detrain_up(i,j,k) = 0.0
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if

if ( flg_entr_dwn ) then
!$OMP DO SCHEDULE(STATIC)
  do k = 1, n_conv_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        entrain_dwn(i,j,k) = entrain_dwn(i,j,k) * g
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if
if ( flg_detr_dwn ) then
!$OMP DO SCHEDULE(STATIC)
  do k = 1, n_conv_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        detrain_dwn(i,j,k) = detrain_dwn(i,j,k) * g
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if

!$OMP END PARALLEL


return
end subroutine comorph_diags_um_proc



end module comorph_diags_um_mod
#endif
