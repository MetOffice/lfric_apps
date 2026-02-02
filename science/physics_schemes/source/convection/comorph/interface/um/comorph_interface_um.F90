! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module comorph_interface_um_mod

use um_types, only: real_umphys

implicit none

contains

! Subroutine to interface the CoMorph convection scheme into the
! UM physics.  Does some interpolation of fields, and conversion
! of specific humidities to and from mixing ratios if needed.
! Also sets up the various derived-type structures to pass
! in and out of comorph_ctl, and calls that routine.
subroutine comorph_interface_um(                                               &
             n_conv_levels, cycleno,                                           &
             l_tracer_on, l_tracer, ntra_fld, i_tr_vars,                       &
             z_rho, z_theta, p_layer_boundaries, p_layer_centres,              &
             rho_dry_th, rho_wet_th, rho_dry, rho_wet,                         &
             exner_theta_levels,                                               &
             taux_p, tauy_p, ftl, fqw, rhokm, bl_w_var,                        &
             zh, dzh, zhnl, zhsc, bl_type_7, fb_surf, u_s,                     &
             ls_rain, ls_snow,                                                 &
             u_p, v_p, w,                                                      &
             theta_n, m_v, m_cl, m_cf, m_cf2, m_r, m_gr,                       &
             precfrac_n, cf_liquid_n, cf_frozen_n, bulk_cf_n,                  &
             ustar_p, vstar_p, r_w,                                            &
             theta_star, q_star, qcl_star, qcf_star, qcf2_star,                &
             qrain_star, qgraup_star,                                          &
             precfrac_star, cf_liquid_star, cf_frozen_star,                    &
             bulk_cf_star,                                                     &
             ccb0, cct0, lcbase0, ccb, cct, lcbase, lctop,                     &
             cca0, ccw0, cclwp0, cca, ccw, cclwp, cca_2d, lcca,                &
             aerosol, dust_div1, dust_div2,                                    &
             dust_div3, dust_div4, dust_div5, dust_div6,                       &
             so2, so4_aitken, so4_accu, so4_diss,                              &
             dms, nh3, soot_new, soot_aged, soot_cld,                          &
             bmass_new, bmass_aged, bmass_cld,                                 &
             ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss,               &
             co2, free_tracers, tracer_ukca, ozone_tracer,                     &
             theta_inc, q_inc, qcl_inc, qcf_inc,                               &
             qcf2_inc, qrain_inc, qgraup_inc,                                  &
             cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,                        &
             pressure_incr_env, dubydt_pout, dvbydt_pout )

use nlsizes_namelist_mod, only: row_length, rows, model_levels,                &
                                bl_levels, tr_vars, tr_ukca
use atm_fields_bounds_mod, only: pdims, pdims_s, wdims, wdims_s,               &
                                 tdims, tdims_s, tdims_l,                      &
                                 array_dims
use model_domain_mod, only: model_type, mt_single_column
use umPrintMgr, only: newline
!$ USE omp_lib, ONLY: omp_get_max_threads             ! Note OpenMP sentinel

use tuning_segments_mod, only: a_convect_seg_size, a_convect_trac_seg_size,    &
                               l_conv_1_seg_per_thread, l_autotune_segments
use gen_phys_inputs_mod, only: l_mr_physics
use mphys_inputs_mod, only: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,             &
                            l_mcr_precfrac
use ukca_option_mod, only: l_ukca_mode
use ukca_d1_defs, only: l_ukca_plume_diags
use s_main_force, only: gridbox_area
use horiz_grid_mod, only: xi1_u, xi2_v
use metric_terms_mod, only: h1_p_eta, h2_p_eta

use autotune_conv_seg_mod, only: autotune_conv_seg_start,                      &
                                 autotune_conv_seg_end
use ukca_scavenging_diags_mod, only: ukca_plume_scav_initial,                  &
                                     ukca_plume_scav_final
use q_to_mix_nocopy_mod, only: q_to_mix_nocopy
use level_heights_mod, only: r_theta_levels, r_rho_levels

! Comorph interface routines
use comorph_conv_cloud_extras_mod, only: comorph_conv_cloud_extras
use set_constants_from_um_mod, only: set_constants_from_um
use calc_conv_incs_mod, only: calc_conv_incs,                                  &
                              i_call_save_before_conv, i_call_diff_to_get_incs
use calc_qcf2_incs_mod, only: calc_qcf2_incs, i_call_combine_in_qcf2,          &
                              i_call_subtract_qcf, i_call_repartition
use fracs_consistency_mod, only: fracs_consistency
use interp_turb_mod, only: interp_turb
use calc_turb_len_mod, only: calc_turb_len
use limit_turb_perts_mod, only: limit_turb_perts
use conv_update_precfrac_mod, only: conv_update_precfrac
use assign_fields_mod, only: assign_fields
use assign_tracers_mod, only: assign_tracers

! Comorph itself
use comorph_ctl_mod, only: comorph_ctl

! Comorph modules
use comorph_diags_um_mod, only: comorph_diags_um_reqs,                         &
                                comorph_diags_um_null,                         &
                                comorph_diags_um_proc
use comorph_diags_scm_mod, only: comorph_diags_scm_reqs,                       &
                                 comorph_diags_scm_proc
use fields_type_mod, only: fields_type, fields_nullify
use grid_type_mod, only: grid_type, grid_nullify
use turb_type_mod, only: turb_type, turb_nullify
use cloudfracs_type_mod, only: cloudfracs_type, cloudfracs_nullify
use comorph_diags_type_mod, only: comorph_diags_type
use raise_error_mod, only: raise_fatal
use comorph_constants_mod, only: l_init_constants, real_hmprec,                &
                                 l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,    &
                                 l_turb_par_gen,                               &
                                 i_convcloud, i_convcloud_liqonly

implicit none

! Highest model-level where convection is allowed
integer, intent(in) :: n_conv_levels
! Model outer loop iteration count
integer, intent(in) :: cycleno

! Flag for whether any tracers are in use
logical, intent(in) :: l_tracer_on
! Flag for whether tracer transport is needed this call
logical, intent(in) :: l_tracer
! Number of tracer species
integer, intent(in) :: ntra_fld
! Address of the first user-defined free tracer in the tracer super-array
integer, intent(in) :: i_tr_vars

! Model-level heights above surface
real(kind=real_umphys), target, intent(in out) :: z_rho                        &
                            ( row_length, rows, model_levels )
real(kind=real_umphys), target, intent(in) :: z_theta                          &
                            ( row_length, rows, model_levels )
! Note: z_rho is INOUT because we temporarily reset the height
! of z_rho(:,:,1) to the surface (zero), as it is used as the
! lower interface of the bottom level.
! This is consistent with the setup of p_layer_boundaries in
! atmos_physics2.

! Pressure at top/bottom interfaces and centres of theta-levels
real(kind=real_umphys), target, intent(in) :: p_layer_boundaries               &
                            ( row_length, rows, 1:model_levels+1)
real(kind=real_umphys), target, intent(in) :: p_layer_centres                  &
                            ( row_length, rows, 0:model_levels )
! Note: for some reason in atmos_physics2 p_layer_boundaries is
! set so that p_layer_boundaries(k) is the pressure at the
! layer-interface just ABOVE p_theta_levels(k), whereas in the UM
! z_rho(k) is defined as the height of the layer interface just
! BELOW z_theta(k).  To avoid this discrepency, the array
! has been declared here with vertical dimensions
! 1:model_levels+1 instead of 0:model_levels, so that
! p_layer_boundaries(k) is defined at the same height as z_rho(k).

! Dry-density on theta-levels
real(kind=real_umphys), target, intent(in) :: rho_dry_th                       &
                            ( row_length, rows, model_levels-1 )
! Wet-density on theta-levels
real(kind=real_umphys), intent(in) :: rho_wet_th                               &
                            ( row_length, rows, model_levels-1 )
! Note: these fields do not exist on the top model-level!
! This means convection can only operate up to model_levels-1

! Dry-density on rho-levels
real(kind=real_umphys), intent(in) ::                                          &
                    rho_dry(pdims%i_end,pdims%j_end,pdims%k_end)

! Wet-density on rho-levels
real(kind=real_umphys), intent(in) ::                                          &
                    rho_wet(pdims%i_end,pdims%j_end,pdims%k_end)

! Exner pressure (p/p0)^R/cp on theta-levels
real(kind=real_umphys), intent(in) :: exner_theta_levels                       &
                    (tdims_s%i_start:tdims_s%i_end,                            &
                     tdims_s%j_start:tdims_s%j_end,                            &
                     tdims_s%k_start:tdims_s%k_end)


! Turbulence fields passed in from the boundary-layer scheme...

! Turbulent wind stress on p grid (on theta-levels) N m-2
real(kind=real_umphys), intent(in) ::                                          &
                    taux_p ( pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,                        &
                             0:bl_levels-1 )
real(kind=real_umphys), intent(in) ::                                          &
                    tauy_p ( pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,                        &
                             0:bl_levels-1 )

! Sensible heat flux / K kg m-2 s-1  (NOT including factor of cp)
real(kind=real_umphys), target, intent(in out) ::                              &
                                ftl ( pdims%i_start:pdims%i_end,               &
                                      pdims%j_start:pdims%j_end,               &
                                      bl_levels )
! Total water flux / kg m-2 s-1
real(kind=real_umphys), target, intent(in out) ::                              &
                                fqw ( pdims%i_start:pdims%i_end,               &
                                      pdims%j_start:pdims%j_end,               &
                                      bl_levels )
! (intent inout so we can change units to remove factor of rho
!  and then change back again).

! rho * turbulent diffusivity for momentum, on theta-grid
real(kind=real_umphys), intent(in) ::                                          &
                    rhokm    ( tdims_s%i_start:tdims_s%i_end,                  &
                               tdims_s%j_start:tdims_s%j_end,                  &
                               0:bl_levels-1 )

! Sub-grid turbulent variance in vertical velocity
real(kind=real_umphys), intent(in) ::                                          &
                    bl_w_var ( tdims%i_start:tdims%i_end,                      &
                               tdims%j_start:tdims%j_end,                      &
                               1:tdims%k_end )

! Boundary layer height / m
real(kind=real_umphys), intent(in) ::                                          &
                    zh       ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
! Resolved inversion thickness / m
real(kind=real_umphys), intent(in) ::                                          &
                    dzh      ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
! Surface-driven non-local BL-top height / m
real(kind=real_umphys), intent(in) ::                                          &
                    zhnl     ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
! Height of decoupled stratocumulus top / m
real(kind=real_umphys), intent(in) ::                                          &
                    zhsc     ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
! Indicator for shear-dominated boundary-layers
real(kind=real_umphys), intent(in) ::                                          &
                    bl_type_7( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
! Surface buoyancy flux / m2 s-3
real(kind=real_umphys), intent(in) ::                                          &
                    fb_surf  ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
! Surface friction velocity / m s-1
real(kind=real_umphys), intent(in) ::                                          &
                    u_s      ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )

! Surface precipitation rates from the "large-scale" microphysics scheme
real(kind=real_umphys), intent(in) ::                                          &
                    ls_rain  ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )
real(kind=real_umphys), intent(in) ::                                          &
                    ls_snow  ( pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end )

! Start-of-timestep primary fields...

! Wind velocity components / ms-1
real(kind=real_umphys), intent(in) :: u_p(pdims%i_end,pdims%j_end,pdims%k_end)
real(kind=real_umphys), intent(in) :: v_p(pdims%i_end,pdims%j_end,pdims%k_end)
real(kind=real_umphys), target, intent(in) ::                                  &
                            w(wdims_s%i_start:wdims_s%i_end,                   &
                              wdims_s%j_start:wdims_s%j_end,                   &
                              wdims_s%k_start:wdims_s%k_end)
! Note: u_p, v_p are interpolated onto the p-grid but are
! still on the layer-interfaces.  CoMorph needs them co-located
! with the other fields, so we will interpolate them onto
! theta-levels.

! Potential temperature / K
real(kind=real_umphys), intent(in) ::                                          &
                    theta_n(tdims_s%i_start:tdims_s%i_end,                     &
                            tdims_s%j_start:tdims_s%j_end,                     &
                            tdims_s%k_start:tdims_s%k_end)

! Start-of-timestep mixing ratios of water species
! Note: in SCM runs, these are not set unless the flag l_mr_conv
!       in cv_run_mod is set to true.
                    ! Water vapour
real(kind=real_umphys), target, intent(in) ::                                  &
                            m_v  (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
                    ! Cloud liquid water
real(kind=real_umphys), target, intent(in) ::                                  &
                            m_cl (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
                    ! Cloud ice and snow
real(kind=real_umphys), target, intent(in) ::                                  &
                            m_cf (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
                    ! 2nd cloud ice category (optional)
real(kind=real_umphys), target, intent(in out) ::                              &
                            m_cf2(tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
                    ! Prognostic rain water (optional)
real(kind=real_umphys), target, intent(in) ::                                  &
                            m_r  (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
                    ! Prognostic graupel (optional)
real(kind=real_umphys), target, intent(in) ::                                  &
                            m_gr (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)

! Rain fraction
real(kind=real_umphys), target, intent(in) :: precfrac_n                       &
                                 (tdims_l%i_start:tdims_l%i_end,               &
                                  tdims_l%j_start:tdims_l%j_end,               &
                                  tdims_l%k_start:tdims_l%k_end)

! Cloud fraction fields
real(kind=real_umphys), target, intent(in) :: cf_liquid_n                      &
                                 (tdims_l%i_start:tdims_l%i_end,               &
                                  tdims_l%j_start:tdims_l%j_end,               &
                                  tdims_l%k_start:tdims_l%k_end)
real(kind=real_umphys), target, intent(in) :: cf_frozen_n                      &
                                 (tdims_l%i_start:tdims_l%i_end,               &
                                  tdims_l%j_start:tdims_l%j_end,               &
                                  tdims_l%k_start:tdims_l%k_end)
real(kind=real_umphys), target, intent(in) :: bulk_cf_n                        &
                                 (tdims_l%i_start:tdims_l%i_end,               &
                                  tdims_l%j_start:tdims_l%j_end,               &
                                  tdims_l%k_start:tdims_l%k_end)

! Latest primary fields (to be updated by convection)

! u and v updated with all increments so far and
! interpolated to p-grid (on rho-levels)
real(kind=real_umphys), intent(in) ::                                          &
                    ustar_p ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              pdims%k_start:pdims%k_end )
real(kind=real_umphys), intent(in) ::                                          &
                    vstar_p ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              pdims%k_start:pdims%k_end )

! Vertical velocity increment:
! IN:  Increment so far this timestep
! OUT: Increment with contribution from convection added as well
real(kind=real_umphys), target, intent(in out) :: r_w                          &
                                 (wdims%i_start:wdims%i_end,                   &
                                  wdims%j_start:wdims%j_end,                   &
                                  1:wdims%k_end)

! Potential temperature
real(kind=real_umphys), target, intent(in out) :: theta_star                   &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)

! Water species
real(kind=real_umphys), target, intent(in out) :: q_star                       &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)
real(kind=real_umphys), target, intent(in out) :: qcl_star                     &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)
real(kind=real_umphys), target, intent(in out) :: qcf_star                     &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)
real(kind=real_umphys), target, intent(in out) :: qcf2_star                    &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)
real(kind=real_umphys), target, intent(in out) :: qrain_star                   &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)
real(kind=real_umphys), target, intent(in out) :: qgraup_star                  &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)

! Rain fraction
real(kind=real_umphys), target, intent(in out) :: precfrac_star                &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)

! Cloud fractions
real(kind=real_umphys), target, intent(in out) :: cf_liquid_star               &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)
real(kind=real_umphys), target, intent(in out) :: cf_frozen_star               &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)
real(kind=real_umphys), target, intent(in out) :: bulk_cf_star                 &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)

! Convective cloud arrays...

! Variables with 0 are prognostics passed around to radiation,
! use depends on whether using PC2, convective cores etc.
! Others are diagnostics.
! NOTE: in CoMorph, none of these are used separately by
! radiation; instead, they are added onto the large-scale cloud
! variables in this routine.  They need to be prognostics so
! that they can be subtracted off again during the next timestep.

! Prognostics for integer level corresponding to convective
! cloud base and top for the highest convecting layer
integer, intent(in out) :: ccb0    (row_length,rows)
integer, intent(in out) :: cct0    (row_length,rows)
! Prognostic for model-level of base of lowest convecting layer
integer, intent(in out) :: lcbase0 (row_length,rows)

! Diagnostics for integer level corresponding to convective
! cloud base and top for the highest convecting layer
integer, intent(in out) :: ccb    (row_length,rows)
integer, intent(in out) :: cct    (row_length,rows)
! Diagnostics for model-level of base and top of lowest layer
integer, intent(in out) :: lcbase (row_length,rows)
integer, intent(in out) :: lctop  (row_length,rows)

! Prognostics for 3-D convective cloud amount and cloud water content
real(kind=real_umphys), intent(in out) :: cca0                                 &
                                    ( row_length, rows, model_levels )
real(kind=real_umphys), intent(in out) :: ccw0                                 &
                                    ( row_length, rows, model_levels )
! Prognostic for vertically-integrated convective liquid water path
real(kind=real_umphys), intent(in out) :: cclwp0  (row_length,rows)

! Diagnostics for 3-D convective cloud amount and cloud water content
real(kind=real_umphys), target, intent(in out) :: cca                          &
                                    ( row_length, rows, model_levels )
real(kind=real_umphys), target, intent(in out) :: ccw                          &
                                    ( row_length, rows, model_levels )
! Diagnostics for 2-D convective cloud amount and
! vertically-integrated water path
real(kind=real_umphys), intent(in out) :: cclwp  (row_length,rows)
real(kind=real_umphys), intent(in out) :: cca_2d (row_length,rows)
! Diagnostic of convective cloud fraction on the lowest
! model-level containing convective cloud
real(kind=real_umphys), intent(in out) :: lcca   (row_length,rows)

! Note: CoMorph only calculates the 3-D fields CCA and CCW.
! The rest of the convective cloud fields must be calculated
! from those after the call to CoMorph.


! Tracer Prognostics - these are time level n plus all
! increments up to this point.  Transport of these fields
! should only be done on the final solver outer loop cycle.

! Aerosol tracer (called "murk" in atm_step).
real(kind=real_umphys), target, intent(in out) :: aerosol                      &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)

! Dust fields
real(kind=real_umphys), target, intent(in out) :: dust_div1                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: dust_div2                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: dust_div3                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: dust_div4                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: dust_div5                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: dust_div6                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)

! Chemistry species
real(kind=real_umphys), target, intent(in out) :: so2                          &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: so4_aitken                   &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: so4_accu                     &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: so4_diss                     &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: dms                          &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: nh3                          &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: soot_new                     &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: soot_aged                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: soot_cld                     &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: bmass_new                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: bmass_aged                   &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: bmass_cld                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: ocff_new                     &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: ocff_aged                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: ocff_cld                     &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: nitr_acc                     &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: nitr_diss                    &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in out) :: co2                          &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)

! Free tracers
real(kind=real_umphys), target, intent(in out) :: free_tracers                 &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end,               &
                                  tr_vars)

! UKCA tracers
real(kind=real_umphys), target, intent(in out) :: tracer_ukca                  &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end,               &
                                  tr_ukca)

! Ozone tracer
real(kind=real_umphys), target, intent(in out) :: ozone_tracer                 &
                                 (tdims_s%i_start:tdims_s%i_end,               &
                                  tdims_s%j_start:tdims_s%j_end,               &
                                  tdims_s%k_start:tdims_s%k_end)

! Convective increments to primary fields, needed by PC2 and for diagnostics
real(kind=real_umphys), intent(out) :: theta_inc                               &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(out) :: q_inc                                   &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(out) :: qcl_inc                                 &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(out) :: qcf_inc                                 &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(in out) :: qcf2_inc                             &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(in out) :: qrain_inc                            &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(in out) :: qgraup_inc                           &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(in out) :: cf_liquid_inc                        &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(in out) :: cf_frozen_inc                        &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
real(kind=real_umphys), intent(in out) :: bulk_cf_inc                          &
                     ( tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                                   1:tdims%k_end )
! Optional hydrometeor species and cloud fraction increment arrays are outputs,
! but need intent(in out) because they are only minimally allocated in
! other_conv_ctl if not used, and in this case intent(out) would imply
! writing them out-of-bounds.

! Pressure change following parcels in the environment
! (needed to calculate the PC2 cloud response to convective
!  subsidence)
real(kind=real_umphys), target, intent(in out) :: pressure_incr_env            &
                                 (tdims%i_start:tdims%i_end,                   &
                                  tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end)

! Convective u,v tendencies on rho-levels, output to atmos_physics2
! to update u,v on their native staggered grids
real(kind=real_umphys), intent(out) :: dubydt_pout                             &
                     ( pdims_s%i_start:pdims_s%i_end,                          &
                       pdims_s%j_start:pdims_s%j_end,                          &
                       pdims_s%k_start:pdims_s%k_end )
real(kind=real_umphys), intent(out) :: dvbydt_pout                             &
                     ( pdims_s%i_start:pdims_s%i_end,                          &
                       pdims_s%j_start:pdims_s%j_end,                          &
                       pdims_s%k_start:pdims_s%k_end )


! Inputs to comorph_ctl:

! Derived type structures storing pointers to the primary fields,
! at start-of-timestep and latest fields
type(fields_type) :: fields_n
type(fields_type) :: fields_np1

! Structure storing pointers to the model-grid fields and dry-rho
type(grid_type) :: grid

! Structure containing pointers to turbulence fields passed
! in from the boundary-layer scheme
type(turb_type) :: turb

! Structure containing diagnostic area fractions for
! cloud and rain
type(cloudfracs_type) :: cloudfracs

! Structure storing meta-data and pointers to fields for all
! diagnostics calculated by CoMorph
type(comorph_diags_type) :: comorph_diags

! Number of columns per segment
integer :: segment_size
! Number of segments
integer :: n_segments


! Store for contents of bottom layer of z_rho, which are
! temporarily set to zero for CoMorph
real(kind=real_umphys) :: z_rho_1( row_length, rows )


! Winds interpolated onto theta-levels
real(kind=real_umphys), target ::                                              &
                u_th_n( tdims%i_start:tdims%i_end,                             &
                        tdims%j_start:tdims%j_end,                             &
                        1:tdims%k_end-1)
real(kind=real_umphys), target ::                                              &
                v_th_n( tdims%i_start:tdims%i_end,                             &
                        tdims%j_start:tdims%j_end,                             &
                        1:tdims%k_end-1)
real(kind=real_umphys), target ::                                              &
                u_th_np1( tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end,                           &
                          1:tdims%k_end-1)
real(kind=real_umphys), target ::                                              &
                v_th_np1( tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end,                           &
                          1:tdims%k_end-1)
! Note: these do not exist at the top theta-level.

! Work array for handling grid-mean vertical velocity inside comorph.
real(kind=real_umphys), target :: w_work                                       &
                                  ( wdims%i_start:wdims%i_end,                 &
                                    wdims%j_start:wdims%j_end,                 &
                                    1:wdims%k_end )

! Array for start-of-timestep temperature
! (needs to be converted from theta to T)
real(kind=real_umphys), target ::                                              &
                temperature_n( tdims%i_start:tdims%i_end,                      &
                               tdims%j_start:tdims%j_end,                      &
                               1:tdims%k_end)

! Work arrays for water species if they are on in CoMorph
! but off in the UM
real(kind=real_umphys), target, allocatable :: q_rain_work(:,:,:)
real(kind=real_umphys), target, allocatable :: q_snow_work(:,:,:)
real(kind=real_umphys), target, allocatable :: q_graup_work(:,:,:)

! Arrays for turbulence fields interpolated onto rho-levels
real(kind=real_umphys), target, allocatable :: w_var_rh(:,:,:)
real(kind=real_umphys), target, allocatable :: fu_rh(:,:,:)
real(kind=real_umphys), target, allocatable :: fv_rh(:,:,:)

! Turbulence lengthscale
real(kind=real_umphys), target :: turb_len                                     &
                                  ( pdims%i_start:pdims%i_end,                 &
                                    pdims%j_start:pdims%j_end,                 &
                                    1:bl_levels )

! Horizontal grid-lengths
real(kind=real_umphys) :: delta_x ( pdims%i_start:pdims%i_end,                 &
                                    pdims%j_start:pdims%j_end )
real(kind=real_umphys) :: delta_y ( pdims%i_start:pdims%i_end,                 &
                                    pdims%j_start:pdims%j_end )

! "Effective" boundary-layer-top height
real(kind=real_umphys) :: zh_eff                                               &
                                  ( pdims%i_start:pdims%i_end,                 &
                                    pdims%j_start:pdims%j_end )
! BL height up-to-which to homogenize convective tendencies inside comorph
real(kind=real_umphys), target :: zh_homog                                     &
                                  ( pdims%i_start:pdims%i_end,                 &
                                    pdims%j_start:pdims%j_end )

! Model-level closest to zh_eff
integer :: k_zh_eff

! Store to save qrain+qgraup before convection, for calculating precip
! fraction increment from convection
real(kind=real_umphys), allocatable :: q_prec_b4(:,:,:)

! Separate convective bulk cloud fraction
! (only used if cca is liquid-cloud only)
real(kind=real_umphys), allocatable, target :: frac_bulk_conv(:,:,:)

! CCA used for calculating other diagnostics
! (points to frac_bulk_conv if cca is liquid-only, or cca otherwise).
real(kind=real_umphys), pointer :: cca_bulk(:,:,:)

! Flags for treating condensed water species diagnostically when
! they are off in the UM but on in CoMorph
logical :: l_temporary_rain
logical :: l_temporary_snow
logical :: l_temporary_graup

! Flag passed into q_to_mix_nocopy to tell it whether to convert
! from specifics to mixing-ratios or vice-versa
logical :: l_reverse

! Structure containing array dims, to pass into q_to_mix_nocopy
type(array_dims) :: dims

! Parcel radius amplification factor passed into comorph
real(kind=real_umphys), target :: par_radius_amp_um(row_length,rows)

! Loop counters
integer :: i, j, k

character(len=*), parameter :: routinename = "COMORPH_INTERFACE_UM"


! Raise an error if real_hmprec has not been set to the same kind
! as the input fields
if ( .not. real_hmprec == kind(z_rho) ) then
  call raise_fatal( routinename,                                               &
         "Host-model precision specified in the CoMorph "     //               &
         "convection scheme does not match the "     //newline//               &
         "kind of the input fields." )
end if


!----------------------------------------------------------------
! 1) Set timestepping, segmenting, array-dimension and moist
!    thermodynamics constants for CoMorph consistent with the UM
!----------------------------------------------------------------

! Set stuff in the CoMorph constants module, if not already set:
if ( .not. l_init_constants ) then
  ! Comorph sets l_init_constants to true the first time it is called,
  ! so this block of code will only be entered before the first
  ! call to comorph, i.e. on the first timestep.

  ! Call routine to set comorph's switches and constants consistent with the UM
  call set_constants_from_um( n_conv_levels, ntra_fld, i_tr_vars )

end if


!----------------------------------------------------------------
! 2) Save values of fields before convection, for use in
!    calculating convective increments needed elsewhere in the UM
!----------------------------------------------------------------

call calc_conv_incs      ( i_call_save_before_conv, z_theta, z_rho,            &
                           u_p, v_p, ustar_p, vstar_p, w, w_work,              &
                           u_th_n, v_th_n, u_th_np1, v_th_np1,                 &
                           theta_star, q_star, qcl_star, qcf_star,             &
                           qcf2_star, qrain_star, qgraup_star,                 &
                           cf_liquid_star, cf_frozen_star, bulk_cf_star,       &
                           dubydt_pout, dvbydt_pout, r_w,                      &
                           theta_inc, q_inc, qcl_inc, qcf_inc,                 &
                           qcf2_inc, qrain_inc, qgraup_inc,                    &
                           cf_liquid_inc, cf_frozen_inc, bulk_cf_inc )


!----------------------------------------------------------------
! 3) Conversions of input fields
!----------------------------------------------------------------

! Convert potential temperature to actual temperature
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE( i, j, k )            &
!$OMP SHARED( tdims, temperature_n, theta_n, theta_star, exner_theta_levels )
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      temperature_n(i,j,k) = theta_n(i,j,k)                                    &
                           * exner_theta_levels(i,j,k)
      theta_star(i,j,k) = theta_star(i,j,k)                                    &
                        * exner_theta_levels(i,j,k)
    end do
  end do
end do
!$OMP END PARALLEL DO

! If the rest of the physics is not using mixing-ratios, the
! _star fields for the water species need to be converted,
! as CoMorph is written using mixing-ratios.
if ( .not. l_mr_physics ) then

  ! Set dimensions of q_star fields; same as tdims, but without
  ! the zero-levels
  dims = tdims
  dims % k_start = 1
  dims % k_len = dims % k_end

  ! Call routine to convert to mixing-ratios
  l_reverse = .false.
  call q_to_mix_nocopy( l_reverse, dims, dims,                                 &
                        q_star, qcl_star, qcf_star,                            &
                        qcf2_star, qrain_star, qgraup_star )

end if

! The "latest" (_star) fields used by convection are values
! interpolated to departure points by SL advection.
! The interpolation is not guaranteed to preserve consistency
! between the cloud fraction and cloud water fields.
! However, various things can go wrong within CoMorph if they
! are inconsistent; especially the routine calc_env_regions,
! which attempts to calculate the internal properties of
! the in-cloud and clear sub-regions of the grid-box.
! E.g. a common problem is if qcf is large but CFF is small,
! the in-cloud qcf can get implausibly large, which can
! cause the phase-change calculations to yield nonsense.
! Therefore, apply safety checks to the cloud-fractions
! before passing them into convection...
call fracs_consistency      ( qcl_star, qcf_star, qcf2_star,                   &
                              qrain_star, qgraup_star,                         &
                              cf_liquid_star, cf_frozen_star, bulk_cf_star,    &
                              precfrac_star )

! Note: if not using PC2, the cloud-fraction _star fields do exist
! but are just set to zero in atm_step_4a.
! Use the PC2 option l_cloud_call_b4_conv to ensure the latest
! cloud fractions are actually calculated before this point.

if ( i_convcloud == i_convcloud_liqonly ) then
  ! CCA / CCW contain liquid-only convective cloud.
  ! Allocate separate array for bulk convective cloud fraction
  allocate( frac_bulk_conv(tdims%i_start:tdims%i_end,                          &
                           tdims%j_start:tdims%j_end,                          &
                           1:tdims%k_end) )
  ! Point CCA used for computing other things at bulk convective cloud amount
  cca_bulk => frac_bulk_conv
else
  ! Minimal allocation when not used
  allocate( frac_bulk_conv(1,1,1) )
  ! Main CCA array already contains bulk convective cloud amount
  cca_bulk => cca
end if

! If prognostic precip fraction is in use:
if ( l_mcr_qrain .and. l_mcr_precfrac ) then
  ! Allocate array to store precip mixing ratio before convection
  allocate( q_prec_b4(tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end,                               &
                      1:tdims%k_end) )
  ! Save precip mass before convection, for use in updating precfrac later
  call conv_update_precfrac    ( i_call_save_before_conv,                      &
                                 qrain_star, qgraup_star,                      &
                                 cca_bulk, q_prec_b4, precfrac_star )
end if

!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k, k_zh_eff )                      &
!$OMP SHARED( row_length, rows, model_levels, z_rho_1, z_rho,                  &
!$OMP         zh_eff, zh_homog, zhnl, dzh, zh, bl_type_7, zhsc,                &
!$OMP         model_type, gridbox_area, h1_p_eta, h2_p_eta, xi1_u, xi2_v,      &
!$OMP         delta_x, delta_y )

! For conservation purposes on theta-levels, the bottom rho-level
! is actually the surface (so that the bottom theta-level extends
! from the surface up to the next rho-level).
! Save the actual height of the bottom rho-level and replace it
! with zero for the remaining calculations:
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    z_rho_1(i,j) = z_rho(i,j,1)
    z_rho(i,j,1) = 0.0
  end do
end do
!$OMP END DO NOWAIT
! Note: this is correct for the interpolation of turbulence
! fields below, since ftl(:,:,1) is actually at the surface,
! while the rest of the array is on rho-levels.

! Set an "effective" boundary-layer height, below-which comorph's
! increments are modified to avoid double-counting the boundary-layer
! scheme's non-local fluxes
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    ! Top height of surface-driven mixing is the surface-driven non-local
    ! BL-top height zhnl plus the inversion thickness dzh
    ! (note dzh defaults to a large negative number when not used
    !  due to the inversion being sub-grid; using MAX to remove these).
    ! Also sometimes the local mixing based BL-top height exceeds zhnl and
    ! replaces it, and is stored in zh in that case, so use zh if it is higher.
    zh_eff(i,j) = max( zhnl(i,j) + max( dzh(i,j), 0.0 ), zh(i,j) )
    if ( bl_type_7(i,j) == 1.0 ) then
      ! For shear-dominated layers, the Sc-top is essentially coupled to the
      ! surface mixed-layer by shear-driven turbulence.  In this case, use
      ! the Sc-top height.  Comorph suppresses the increments from convection
      ! that only occurs within the layer below zh_eff, so this effectively
      ! disables convection within the shear-dominated layer.
      zh_eff(i,j) = max( zh_eff(i,j), zhsc(i,j) )
    end if
    ! Find model-level straddling zh_eff
    k_zh_eff = 0
    do k = 1, model_levels-1
      if ( zh_eff(i,j)>=z_rho(i,j,k) .and. zh_eff(i,j)<z_rho(i,j,k+1) ) then
        k_zh_eff = k
      end if
    end do
    ! Fully homogenize the model-level straddling the BL-top,
    ! so round zh_eff up to the next rho-level
    zh_homog(i,j) = z_rho(i,j,k_zh_eff+1)
  end do
end do
!$OMP END DO NOWAIT

! Compute horizontal grid-sizes
if ( model_type == mt_single_column ) then
  ! Single-column model runs
  do j = 1, rows
    do i = 1, row_length
      ! Convert SCM grid-box area (in km2) to m2
      ! Assume a square so just take sqrt
      delta_x(i,j) = sqrt( gridbox_area(i,j) * 1.0e6 )
      delta_y(i,j) = delta_x(i,j)
    end do
  end do
else
  ! 3D model runs
!$OMP DO SCHEDULE(STATIC)
  do j = 1, rows
    do i = 1, row_length
      ! Compute grid dx and dy at surface from UM metric and coordinate arrays
      delta_x(i,j) = h1_p_eta(i,j,0) * (xi1_u(i) - xi1_u(i-1))
      delta_y(i,j) = h2_p_eta(i,j,0) * (xi2_v(j) - xi2_v(j-1))
    end do
  end do
!$OMP END DO NOWAIT
end if

!$OMP END PARALLEL

! If using turbulence-based parcel perturbations, need to convert
! the turbulence fields into the required units, and interpolate
! some of them onto rho-levels
if ( l_turb_par_gen ) then

  ! Allocate arrays for momentum diffusivities and fluxes
  ! interpolated onto rho-levels, and turbulence lengthscale
  allocate( w_var_rh ( pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end,                              &
                       1:bl_levels ) )
  allocate( fu_rh    ( pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end,                              &
                       1:bl_levels ) )
  allocate( fv_rh    ( pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end,                              &
                       1:bl_levels ) )

  ! Interpolate momentum diffusivity and fluxes onto rho-levels
  call interp_turb( z_rho, z_theta,                                            &
                    bl_w_var, fb_surf, zh, u_s, taux_p, tauy_p,                &
                    w_var_rh, fu_rh, fv_rh )

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE( i, j, k )            &
!$OMP SHARED( row_length, rows, bl_levels, fu_rh, fv_rh, rho_wet,              &
!$OMP         ftl, fqw, rho_dry )
  do k = 1, bl_levels
    do j = 1, rows
      do i = 1, row_length
        ! Convert wind stresses -rho <w'u'> / N m-2
        ! into just <w'u'> / m2 s-2
        fu_rh(i,j,k) = -fu_rh(i,j,k) / rho_wet(i,j,k)
        fv_rh(i,j,k) = -fv_rh(i,j,k) / rho_wet(i,j,k)
        ! Convert sensible heat flux rho <w'Tl'> / K kg m-2
        ! into just <w'Tl'>
        ftl(i,j,k) = ftl(i,j,k) / rho_wet(i,j,k)
        ! Convert total water flux rho_dry <w'q'> / kg m-2
        ! into just <w'q'>
        ! (where q denotes water mixing-ratio)
        fqw(i,j,k) = fqw(i,j,k) / rho_dry(i,j,k)
        ! Note: dividing by dry-density here instead of wet, so
        ! that the flux is expressed in terms of mixing-ratio
        ! rather than specific humidity
      end do
    end do
  end do
!$OMP END PARALLEL DO

  ! Calculate turbulence lengthscale used to set parcel initial radius
  call calc_turb_len    ( zh_eff, z_theta, z_rho, rho_wet_th, m_v,             &
                          rhokm, bl_w_var, ls_rain, ls_snow, w,                &
                          delta_x, delta_y,                                    &
                          turb_len, par_radius_amp_um )

  ! Check for instances of fluxes too big relative to the turbulent
  ! w-variance (causes excessive parcel perturbations);
  ! increase the w-variance where needed to avoid the problem
  call limit_turb_perts    ( z_theta, z_rho, p_layer_centres, theta_star,      &
                             ftl, fqw, fu_rh, fv_rh, w_var_rh )

end if  ! ( l_turb_par_gen )


!----------------------------------------------------------------
! 4) Initialise temporary arrays needed by comorph
!----------------------------------------------------------------

! Throw an error if any water species are switched on in the UM but
! switched off in comorph, since in this case comorph will not
! transport them consistently (e.g. qcf2 needs to be transported
! along with cf_frozen if it is used).
if ( ( .not. l_cv_cf ) .or.                                                    &
     ! Ice-cloud is always on in the UM, so must be on in comorph
     ( l_mcr_qrain .and. ( .not. l_cv_rain ) ) .or.                            &
     ! UM prognostic rain requires rain to be on in comorph
     ( l_mcr_qgraup .and. ( .not. l_cv_graup ) ) ) then
     ! UM prognostic graupel requires graupel to be on in comorph
  call raise_fatal( routinename,                                               &
         "At least one condensed water species is switched "  //               &
         "on in the UM but switched off in CoMorph." //newline//               &
         "Code has not yet been implemented to handle this "  //               &
         "combination consistently." )
end if

! Temporary fields needed for water species if they are switched
! on in comorph, but switched off in the UM...
l_temporary_rain = l_cv_rain .and. (.not. l_mcr_qrain)
l_temporary_snow = l_cv_snow .and. (.not. l_mcr_qcf2)
l_temporary_graup = l_cv_graup .and. (.not. l_mcr_qgraup)

! In this case, CoMorph will operate using an initial rain / snow / graupel
! mixing ratio which is zero everywhere.  Any rain / snow / graupel
! produced by CoMorph will be either rained-out this timestep or added to
! the cloud-water, to let the large-scale microphysics rain it out.

! Allocate and initialise temporary water species to zero if needed...
if ( l_temporary_rain ) then
  allocate( q_rain_work(tdims%i_start:tdims%i_end,                             &
                        tdims%j_start:tdims%j_end,                             &
                        1:tdims%k_end) )
end if
if ( l_temporary_snow ) then
  allocate( q_snow_work(tdims%i_start:tdims%i_end,                             &
                        tdims%j_start:tdims%j_end,                             &
                        1:tdims%k_end) )
end if
if ( l_temporary_graup ) then
  allocate( q_graup_work(tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                         1:tdims%k_end) )
end if

!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k )                                &
!$OMP SHARED( tdims, row_length, rows, model_levels,                           &
!$OMP         l_temporary_rain, q_rain_work, l_temporary_snow, q_snow_work,    &
!$OMP         l_temporary_graup, q_graup_work, i_convcloud, frac_bulk_conv )

if ( l_temporary_rain ) then
!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        q_rain_work(i,j,k) = 0.0
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if
if ( l_temporary_snow ) then
!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        q_snow_work(i,j,k) = 0.0
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if
if ( l_temporary_graup ) then
!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        q_graup_work(i,j,k) = 0.0
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if

if ( i_convcloud == i_convcloud_liqonly ) then
  ! Initialise local array for bulk conv cloud fraction to zero
!$OMP DO SCHEDULE(STATIC)
  do k = 1, model_levels
    do j = 1, rows
      do i = 1, row_length
        frac_bulk_conv(i,j,k) = 0.0
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if

!$OMP END PARALLEL

if ( l_mcr_qcf2 .and. ( .not. l_cv_snow ) ) then
  ! 2nd ice category switched on in the UM but not in comorph
  ! Put all the ice-cloud mass in the qcf2 field to pass into comorph
  call calc_qcf2_incs( i_call_combine_in_qcf2,                                 &
                       m_cf, m_cf2, qcf_star, qcf2_star, qcf_inc, qcf2_inc )
end if

!----------------------------------------------------------------
! 5) Assign pointers to fields to pass into comorph
!----------------------------------------------------------------

! All the array arguments to this routine are passed into comorph_ctl,
! but we pass them in via pointers contained in the derived-type
! structures grid, turb, cloudfracs, fields_n, fields_np1.
! This routine assigns the pointers to the arrays.
call assign_fields      ( z_theta, z_rho, p_layer_centres, p_layer_boundaries, &
                          r_theta_levels,                                      &
                          rho_dry_th, w_var_rh, ftl, fqw, fu_rh, fv_rh,        &
                          turb_len, par_radius_amp_um, zh_homog,               &
                          cca, ccw, frac_bulk_conv,                            &
                          u_th_n, v_th_n, w, temperature_n,                    &
                          m_v, m_cl, m_cf,                                     &
                          m_cf2, m_r, m_gr,                                    &
                          cf_liquid_n, cf_frozen_n, bulk_cf_n,                 &
                          u_th_np1, v_th_np1, w_work, theta_star,              &
                          q_star, qcl_star, qcf_star,                          &
                          qcf2_star, qrain_star, qgraup_star,                  &
                          cf_liquid_star, cf_frozen_star, bulk_cf_star,        &
                          precfrac_star, l_temporary_snow,                     &
                          l_temporary_rain, l_temporary_graup,                 &
                          q_snow_work, q_rain_work, q_graup_work,              &
                          grid, turb, cloudfracs, fields_n, fields_np1 )


!----------------------------------------------------------------
! 6) Setup tracers if required
!----------------------------------------------------------------
if ( l_tracer ) then

  if ( tr_ukca > 0 .and. l_ukca_mode .and. l_ukca_plume_diags ) then
    ! Save values of UKCA tracers before convection, for calculating
    ! diagnostics of increments from plume scavenging
    call ukca_plume_scav_initial(tr_ukca, tracer_ukca)
  end if

  ! Allocate list of pointers to tracer fields
  allocate( fields_np1 % tracers(ntra_fld) )
  ! Call routine to point pointers at the tracer arrays
  call assign_tracers(                                                         &
             aerosol, dust_div1, dust_div2,                                    &
             dust_div3, dust_div4, dust_div5, dust_div6,                       &
             so2, so4_aitken, so4_accu, so4_diss,                              &
             dms, nh3, soot_new, soot_aged, soot_cld,                          &
             bmass_new, bmass_aged, bmass_cld,                                 &
             ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss,               &
             co2, free_tracers, tracer_ukca, ozone_tracer,                     &
             fields_np1 )

end if ! ( l_tracer )


!----------------------------------------------------------------
! 7) Diagnostic requests
!----------------------------------------------------------------
! CoMorph has an internal diagnostic system; meta-data for
! each diagnostic, and pointers to the output arrays, are
! contained in the derived-type structure comorph_diags.

! This routine USEs the UM's convection diagnostic arrays and
! flags from modules cv_diagnostic_array_mod, cv_stash_flg_mod
! and uses them to setup the pointers and diag request flags
! contained in the comorph_diags structure.
call comorph_diags_um_reqs( comorph_diags,                                     &
                            pressure_incr_env )

! Additional diagnostics for the SCM.  This is stored in a
! module which contains its own super-arrays for SCM diagnostics
if ( model_type == mt_single_column ) then
  call comorph_diags_scm_reqs( comorph_diags, l_tracer )
end if


!----------------------------------------------------------------
! 8) Call CoMorph
!----------------------------------------------------------------

if ( l_autotune_segments ) then
  ! Autotune the segment size used in comorph
  call autotune_conv_seg_start( l_tracer_on, l_tracer, cycleno,                &
                                row_length*rows, segment_size )

else
  ! Set segment size from UM namelist
  if ( l_tracer .and. a_convect_trac_seg_size > 0 ) then
    segment_size = a_convect_trac_seg_size
  else if ( a_convect_seg_size > 0 ) then
    segment_size = a_convect_seg_size
  else
    segment_size = row_length * rows
  end if
end if

! Set number of segments needed to get the desired segment size
n_segments = ceiling( real(row_length*rows) / real(segment_size) )

if ( l_conv_1_seg_per_thread ) then
  ! Override to just run one segment per OMP thread...

  ! If not using OMP, just run with a single segment for all columns.
  n_segments = 1
  ! Under OMP sentinel, set number of segments equal to number of threads
!$  n_segments = omp_get_max_threads()
end if

! The moment you've been waiting for!
call comorph_ctl( l_tracer, n_segments,                                        &
                  grid, turb, cloudfracs, fields_n, fields_np1,                &
                  comorph_diags )

if ( l_autotune_segments ) then
  ! Finish autotuning
  call autotune_conv_seg_end( l_tracer_on, l_tracer, cycleno )
end if

! Nullify pointers used to pass array arguments through CoMorph
call grid_nullify( grid )
call turb_nullify( turb )
call cloudfracs_nullify( cloudfracs )
call fields_nullify( fields_n )
call fields_nullify( fields_np1 )

! Deallocate list of pointers to tracer fields
if ( l_tracer )  deallocate( fields_np1 % tracers )

! Nullify diagnostic pointers
call comorph_diags_um_null( comorph_diags )

if ( l_tracer ) then
  if ( tr_ukca > 0 .and. l_ukca_mode .and. l_ukca_plume_diags ) then
    ! Use saved values of UKCA tracers before convection to calculate
    ! diagnostics of increments from plume scavenging
    call ukca_plume_scav_final(tr_ukca, tracer_ukca)
  end if
end if


!----------------------------------------------------------------
! 9) Conversions of output fields back into UM "format"
!----------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k )                                &
!$OMP SHARED( row_length, rows, tdims, z_rho, z_rho_1,                         &
!$OMP         theta_star, exner_theta_levels, qcl_star, qcf_star,              &
!$OMP         l_temporary_rain, q_rain_work, l_temporary_snow, q_snow_work,    &
!$OMP         l_temporary_graup, q_graup_work )

! Restore bottom level of z_rho, ready for whatever comes next...
!$OMP DO SCHEDULE(STATIC)
do j = 1, rows
  do i = 1, row_length
    z_rho(i,j,1) = z_rho_1(i,j)
  end do
end do
!$OMP END DO NOWAIT

! Convert final temperature back into potential temperature
!$OMP DO SCHEDULE(STATIC)
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      theta_star(i,j,k) = theta_star(i,j,k)                                    &
                        / exner_theta_levels(i,j,k)
    end do
  end do
end do
!$OMP END DO NOWAIT

! If any condensed water species are on in CoMorph but off
! in the UM, scatter any mixing-ratio of these species produced
! into the UM's equivalent fields

if ( l_temporary_rain ) then
  ! Add rain to the liquid cloud if prognostic rain is off
!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qcl_star(i,j,k) = qcl_star(i,j,k) + q_rain_work(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if

if ( l_temporary_snow ) then
  ! Add snow to single ice field qcf if 2nd ice category is off
!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qcf_star(i,j,k) = qcf_star(i,j,k) + q_snow_work(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if

if ( l_temporary_graup ) then
  ! Add graupel to qcf field if prognostic graupel is off
!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qcf_star(i,j,k) = qcf_star(i,j,k) + q_graup_work(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if

!$OMP END PARALLEL

! Note: it might be preferable to diagnostically rain-out these
! temporary fields instead of adding them to qcl or qcf, but
! that functionality hasn't been implemented yet.

! Deallocate temporary arrays for hydrometeor species
if ( l_temporary_graup )  deallocate( q_graup_work )
if ( l_temporary_snow )   deallocate( q_snow_work )
if ( l_temporary_rain )   deallocate( q_rain_work )

if ( l_mcr_qcf2 .and. ( .not. l_cv_snow ) ) then
  ! 2nd ice category switched on in the UM but not in comorph
  ! Subtract qcf off from qcf2 again
  call calc_qcf2_incs( i_call_subtract_qcf,                                    &
                       m_cf, m_cf2, qcf_star, qcf2_star, qcf_inc, qcf2_inc )
end if

! If prognostic precip fraction is in use:
if ( l_mcr_qrain .and. l_mcr_precfrac ) then
  ! Update the precip fraction using the convective rain and graupel increments
  call conv_update_precfrac    ( i_call_diff_to_get_incs,                      &
                                 qrain_star, qgraup_star,                      &
                                 cca_bulk, q_prec_b4, precfrac_star )
  ! Deallocate saved precip mass before convection, now we're done
  deallocate( q_prec_b4 )
end if

! Note: currently CoMorph does not compute any cloud fraction
! increment associated with transfer of condensate by
! precipitation.  If the parcel contains ice, some of which
! falls out into the environment, but detrainment is zero,
! this will yield zero ice-cloud fraction but non-zero ice
! mixing-ratio, which shouldn't be allowed.
! Pending implementation of code in CoMorph to compute the
! appropriate cloud-fraction associated with precip fall,
! here there is a temporary fix to avoid the problem.
! A minimum limit is applied to the cloud-fractions after
! convection, proportional to the mixing-ratio.
call fracs_consistency      ( qcl_star, qcf_star, qcf2_star,                   &
                              qrain_star, qgraup_star,                         &
                              cf_liquid_star, cf_frozen_star, bulk_cf_star,    &
                              precfrac_star )

! Convert final mixing ratios back to specific moisture contents
if ( .not. l_mr_physics ) then

  ! Call routine to convert back to specific moisture contents
  l_reverse = .true.
  call q_to_mix_nocopy( l_reverse, dims, dims,                                 &
                        q_star, qcl_star, qcf_star,                            &
                        qcf2_star, qrain_star, qgraup_star,                    &
                        q_cl_conv = ccw )
  ! Note also converting the convective cloud water mixing-ratio output
  ! by comorph to a specific humidity.

end if

!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k )                                &
!$OMP SHARED( row_length, rows, bl_levels,                                     &
!$OMP         fqw, ftl, rho_dry, rho_wet )

! If using turbulence-based parcel perturbations
if ( l_turb_par_gen ) then
  ! Convert scalar fluxes back from <w'q'> / kg kg-1 ms-1
  ! into rho <w'q'> / kg m-2 s-1
!$OMP DO SCHEDULE(STATIC)
  do k = 1, bl_levels
    do j = 1, rows
      do i = 1, row_length
        fqw(i,j,k) = fqw(i,j,k) * rho_dry(i,j,k)
        ftl(i,j,k) = ftl(i,j,k) * rho_wet(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
end if  ! ( l_turb_par_gen )

!$OMP END PARALLEL

if ( l_turb_par_gen ) then
  ! Deallocate the interpolated momentum diffusivity and fluxes
  ! on rho-levels
  deallocate( fv_rh )
  deallocate( fu_rh )
  deallocate( w_var_rh )
end if


! Calculate extra convective cloud fields from cca0, ccw0
call comorph_conv_cloud_extras(                                                &
             n_conv_levels, rho_dry_th, rho_wet_th,                            &
             r_theta_levels, r_rho_levels,                                     &
             cca, ccw, cca_bulk,                                               &
             cclwp, cca_2d, lcca, ccb, cct, lcbase, lctop,                     &
             cca0, ccw0, cclwp0, ccb0, cct0, lcbase0 )


!----------------------------------------------------------------
! 10) Calculate convective increments required elsewhere in the UM
!----------------------------------------------------------------

call calc_conv_incs      ( i_call_diff_to_get_incs, z_theta, z_rho,            &
                           u_p, v_p, ustar_p, vstar_p, w, w_work,              &
                           u_th_n, v_th_n, u_th_np1, v_th_np1,                 &
                           theta_star, q_star, qcl_star, qcf_star,             &
                           qcf2_star, qrain_star, qgraup_star,                 &
                           cf_liquid_star, cf_frozen_star, bulk_cf_star,       &
                           dubydt_pout, dvbydt_pout, r_w,                      &
                           theta_inc, q_inc, qcl_inc, qcf_inc,                 &
                           qcf2_inc, qrain_inc, qgraup_inc,                    &
                           cf_liquid_inc, cf_frozen_inc, bulk_cf_inc )

if ( l_mcr_qcf2 .and. ( .not. l_cv_snow ) ) then
  ! 2nd ice category switched on in the UM but not in comorph
  ! Repartition the ice-cloud increment between crystals and aggregates.
  call calc_qcf2_incs( i_call_repartition,                                     &
                       m_cf, m_cf2, qcf_star, qcf2_star, qcf_inc, qcf2_inc )
end if


!----------------------------------------------------------------
! 11) Processing of output diagnostics
!----------------------------------------------------------------

call comorph_diags_um_proc( n_conv_levels, z_rho, z_theta, lcbase )

! Finished with array for bulk convective cloud amount
cca_bulk => null()
deallocate( frac_bulk_conv )

! Additional diagnostics for the SCM
if ( model_type == mt_single_column ) then
  call comorph_diags_scm_proc( )
end if


return
end subroutine comorph_interface_um


end module comorph_interface_um_mod
#endif
