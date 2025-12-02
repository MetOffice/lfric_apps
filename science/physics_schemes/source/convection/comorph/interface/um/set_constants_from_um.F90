! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_constants_from_um_mod

implicit none

contains

! Subroutine to set switches and constants in comorph based on the
! settings in the UM.  This overrides comorph's default settings.
subroutine set_constants_from_um( n_conv_levels, ntra_fld, i_tr_vars )

! UM settings and constants used from the various UM modules
use timestep_mod, only: um_timestep => timestep
use nlsizes_namelist_mod, only: row_length, rows, bl_levels
use cloud_inputs_mod, only: i_cld_vn
use pc2_constants_mod, only: i_cld_pc2
use ukca_option_mod, only: l_ukca, l_ukca_plume_scav
use model_domain_mod, only: model_type, mt_single_column
use planet_constants_mod, only: g, r, rv, cp
use water_constants_mod, only: tm, lc, lf,                                     &
                               um_rho_liq => rho_water,                        &
                               um_rho_ice => rho_ice
use gen_phys_inputs_mod, only: l_mr_physics
use mphys_inputs_mod, only: l_mcr_precfrac, l_subgrid_graupel_frac

use comorph_um_namelist_mod, only:                                             &
                      par_radius_knob,                                         &
                      par_radius_evol_method_um => par_radius_evol_method,     &
                      n_dndraft_types_um        => n_dndraft_types,            &
                      l_core_ent_cmr_um         => l_core_ent_cmr,             &
                      core_ent_fac_um           => core_ent_fac,               &
                      par_gen_pert_fac_um       => par_gen_pert_fac,           &
                      par_gen_rhpert_um         => par_gen_rhpert,             &
                      cf_conv_fac_um            => cf_conv_fac,                &
                      autoc_opt_um              => autoc_opt,                  &
                      coef_auto_um              => coef_auto,                  &
                      q_cl_auto_um              => q_cl_auto,                  &
                      drag_coef_par_um          => drag_coef_par,              &
!
                      par_gen_mass_fac_um       => par_gen_mass_fac,           &
                      wind_w_fac_um             => wind_w_fac,                 &
                      wind_w_buoy_fac_um        => wind_w_buoy_fac,            &
                      ass_min_radius_um         => ass_min_radius,             &
                      par_gen_core_fac_um       => par_gen_core_fac,           &
                      ent_coef_um               => ent_coef,                   &
                      overlap_power_um          => overlap_power,              &
                      rho_rim_um                => rho_rim,                    &
                      hetnuc_temp_um            => hetnuc_temp,                &
                      drag_coef_cond_um         => drag_coef_cond,             &
                      vent_factor_um            => vent_factor,                &
                      col_eff_coef_um           => col_eff_coef,               &
                      r_fac_tdep_n

! comorph settings and constants set by this routine
use comorph_constants_mod, only: real_cvprec, nx_full, ny_full,                &
                                 k_bot_conv, k_top_conv, k_top_init,           &
                                 n_tracers, n_dndraft_types,                   &
                                 l_cv_cloudfrac, l_tracer_scav,                &
                                 l_calc_cape, l_calc_mfw_cape,                 &
                                 l_spherical_coord, l_approx_dry_adiabat,      &
                                 comorph_timestep,                             &
                                 gravity, melt_temp, R_dry, R_vap,             &
                                 cp_dry, cp_vap, cp_liq, cp_ice,               &
                                 L_con_ref, L_fus_ref, rho_liq, rho_ice,       &
                                 indi_thresh, zero, one,                       &
                                 params_cl, params_rain,                       &
                                 params_cf, params_snow, params_graup,         &
                                 i_sg_homog, i_sg_frac_liq, i_sg_frac_ice,     &
                                 i_sg_frac_prec, tracer_positive,              &
                                 par_gen_pert_fac, par_gen_rhpert,             &
                                 par_gen_radius_fac, par_radius_evol_method,   &
                                 l_core_ent_cmr, core_ent_fac, cf_conv_fac,    &
                                 autoc_opt, coef_auto, q_cl_auto,              &
                                 drag_coef_par, rho_rim,                       &
                                 par_gen_mass_fac, wind_w_fac, wind_w_buoy_fac,&
                                 ass_min_radius, par_gen_core_fac, ent_coef,   &
                                 overlap_power, fac_tdep_n, hetnuc_temp,       &
                                 drag_coef_cond, vent_factor, col_eff_coef

implicit none

! Number of convection levels, set other_conv_ctl and passed in
integer, intent(in) :: n_conv_levels

! Number of tracer fields
integer, intent(in) :: ntra_fld

! Address of the first user-defined free tracer in the tracer super-array
integer, intent(in) :: i_tr_vars

! Loop counter
integer :: i_field


! Set CoMorph convection timestep equal to model timestep
! (no sub-stepping implemented so-far)
comorph_timestep = real(um_timestep,real_cvprec)

! Set CoMorph array sizes
nx_full = row_length
ny_full = rows
k_bot_conv = 1
k_top_conv = n_conv_levels
k_top_init = bl_levels - 1
n_tracers = ntra_fld

! For now, force CoMorph's cloud-fraction switch to be
! consistent with whether or not PC2 is on in the UM.
! If used, the cloud fractions are treated as primary
! prognostic fields within CoMorph.
l_cv_cloudfrac = ( i_cld_vn == i_cld_pc2 )
! Note: l_cv_cloudrac does not affect the convection
! simulated by CoMorph at all; it only switches on/off the
! calculation of cloud-fraction increments.
! Under some settings, CoMorph is also sensitive to the input
! cloud-fractions.  For this purpose, the cloud-fractions
! should be added to the separate input structure cloudfracs,
! in which the cloud fractions are treated as diagnostic
! inputs instead of primary fields.

! The UM does not account for the effect of moisture on R/cp,
! so set switch to use approximation R/cp = R_dry/cp_dry
l_approx_dry_adiabat = .true.

! Set switch for scavenging of tracers by convective precip production;
! currently only implemented for UKCA aerosol and chemistry fields
l_tracer_scav = l_ukca .and. l_ukca_plume_scav

! Set flag to calculate mass-flux-weighted CAPE (for SKEB and diagnostics)
l_calc_mfw_cape = .true.
! SCM also outputs straight CAPE as a diagnostic
l_calc_cape = (model_type==mt_single_column)

! Set gravitational acceleration
gravity = real( g, real_cvprec )

! Set thermodynamics constants
melt_temp = real( tm, real_cvprec )
R_dry = real( r, real_cvprec )
R_vap = real( rv, real_cvprec )
cp_dry = real( cp, real_cvprec )
L_con_ref = real( lc, real_cvprec )
L_fus_ref = real( lf, real_cvprec )
rho_liq = real( um_rho_liq, real_cvprec )
rho_ice = real( um_rho_ice, real_cvprec )

! Set specific heat capacities of the water species;
! this determines how the actual latent heats vary with
! temperature.
! In its latent heating calculations, the UM neglects
! both the temperature dependence of the latent heats, and the
! dependence of the heat capacity on the amount of water in the
! parcel (due to each water phase having a different heat
! capacity to dry air).  These appproximations can only be
! reproduced in CoMorph by setting the heat capacities of all
! 3 water phases equal to eachother.
if ( l_mr_physics ) then
  ! If using mixing-ratio physics, UM-consistent setting is
  ! for all 3 heat capacities to be zero
  cp_vap = zero
  cp_liq = zero
  cp_ice = zero
  ! Yields dT = L_ref/cp_dry dq
  ! where dq is a mixing-ratio increment.
else
  ! Otherwise, UM-consistent setting is
  ! for all 3 heat capacities to be equal to that of dry air
  cp_vap = real( cp, real_cvprec )
  cp_liq = real( cp, real_cvprec )
  cp_ice = real( cp, real_cvprec )
  ! Yields dT = L_ref/cp_dry dq
  ! where dq is a specific humidity increment.
end if

! Set which condensed water fields exist within which sub-grid fractions...

! Liquid and ice cloud exist in the liquid and ice cloud fractions
params_cl % i_sg = i_sg_frac_liq
params_cf % i_sg = i_sg_frac_ice
! UM's 2nd ice category also resides in the ice-cloud fraction if used
params_snow % i_sg = i_sg_frac_ice

if ( l_mcr_precfrac ) then
  ! If using prognostic precip fraction...
  ! Rain exists within the precip fraction
  params_rain % i_sg = i_sg_frac_prec
  ! Graupel optionally in the precip fraction or homogeneous within grid-box
  if ( l_subgrid_graupel_frac ) then
    params_graup % i_sg = i_sg_frac_prec
  else
    params_graup % i_sg = i_sg_homog
  end if
else
  ! No precip fraction field; assume rain and graupel are homogeneous
  params_rain  % i_sg = i_sg_homog
  params_graup % i_sg = i_sg_homog
end if

! Set list of flags for whether tracers are allowed to go negative
allocate( tracer_positive(n_tracers) )
! Nearly all UM tracers are not allowed to have negative values:
do i_field = 1, n_tracers
  tracer_positive(i_field) = .true.
end do

! Number of downdraughts types
n_dndraft_types = n_dndraft_types_um

! Scale default parcel radius factor by tuning knob from the UM
par_gen_radius_fac = par_gen_radius_fac * real( par_radius_knob, real_cvprec )

! Switch controlling how parcel radius evolves with height in the plume
par_radius_evol_method = par_radius_evol_method_um

! Overwrite turbulent parcel perturbation factor with UM value
par_gen_pert_fac = real( par_gen_pert_fac_um, real_cvprec )

! Overwrite neutrally-buoyant moisture perturbation with UM value
par_gen_rhpert = real( par_gen_rhpert_um, real_cvprec )

! Overwrite parameter controlling parcel core dilution
l_core_ent_cmr = l_core_ent_cmr_um
core_ent_fac   = real( core_ent_fac_um, real_cvprec )

! Parameter scaling convective cloud fraction
cf_conv_fac    = real( cf_conv_fac_um, real_cvprec )

! In-parcel cloud-to-rain autoconversion option
autoc_opt      = autoc_opt_um

! In-parcel cloud-to-rain autoconversion rate coefficient
coef_auto      = real( coef_auto_um, real_cvprec )

! In-parcel cloud-to-rain autoconversion threshold liquid water content
q_cl_auto      = real( q_cl_auto_um, real_cvprec )

! Drag coefficient applied to the parcel u,v winds
drag_coef_par  = real( drag_coef_par_um, real_cvprec )

! Dimensionless constant scaling the initiation mass-sources
par_gen_mass_fac = real(par_gen_mass_fac_um, real_cvprec )

! Prescribed draft vertical velocity excess
wind_w_fac       = real(wind_w_fac_um, real_cvprec )

! Tuning constant for buoyancy-dependent convective fraction
wind_w_buoy_fac  = real(wind_w_buoy_fac_um, real_cvprec )

! Minimum parcel initial radius
ass_min_radius   = real(ass_min_radius_um, real_cvprec )

! Scaling factor for par_gen core perturbations relative to
! the parcel mean properties (used if l_par_core = .true.)
par_gen_core_fac = real(par_gen_core_fac_um, real_cvprec )

! Entrainment
! Mixing entrainment rate (m-1) = ent_coef / parcel radius
ent_coef         = real(ent_coef_um, real_cvprec )

! Power for overlap between liquid and ice cloud fractions
! inside the parcel
overlap_power    = real(overlap_power_um, real_cvprec )

! Density of rimed ice (used for graupel)
rho_rim  = real(rho_rim_um, real_cvprec )

! Temperature-dependent ice number concentration slope
! The number concentration n(T) will be given by:
! n(T) = n0 exp( fac_tdep_n ( T - Tmelt ) )
! ( but limited above Tmelt and below T_homnuc)
fac_tdep_n  = -one/real( r_fac_tdep_n, real_cvprec )

! Heterogeneous nucleation temeprature / K
hetnuc_temp = real(hetnuc_temp_um, real_cvprec )

! Asymptotic drag coefficient for a sphere at high Reynolds
! number limit
drag_coef_cond = real(drag_coef_cond_um, real_cvprec )

! Coefficent scaling a term added onto the vapour and heat diffusion,
! for additional exchange due to fall-speed ventilation
vent_factor = real(vent_factor_um, real_cvprec )

! Coefficient for reduction of collection efficiency by
! deflection flow around hydrometeors
col_eff_coef = real(col_eff_coef_um, real_cvprec )

! Set threshold for using indirect indexing versus straight do-loops
! over all points in various calculations inside comorph.
! A value somewhere between 0 and 1 should optimise performance,
! however it breaks bit-reproducibility on different processor
! decompositions, as the indirect indexing loops aren't vectorised,
! whereas the straight do-loops are, and vectorisation can cause
! bit-level changes in answers.
! Tests with CCE have indicated that faster run-times are actually obtained
! by never using indirect indexing (i.e. always performing straight do-loops
! over all points, even when only needed at a few points), rather than
! always using indrect indexing, as the speed-up from vectorisation
! exceeds the penalty from doing redundant calculations.
! Therefore, set this to zero to disable indirect indexing:
indi_thresh = zero


return
end subroutine set_constants_from_um

end module set_constants_from_um_mod
