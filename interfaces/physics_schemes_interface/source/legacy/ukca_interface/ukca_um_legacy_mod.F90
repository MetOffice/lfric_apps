! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Module to satisfy legacy UM dependencies in UKCA.
!
! Method:
!
!  This module provides access from UKCA to UM variables, parameters and
!  procedures that are required by UKCA when running in the UM and are
!  not currently provided via the UKCA API.
!  It acts as a collation point for these UM legacy items that are defined in
!  various UM modules rather than including any definitions itself.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

module ukca_um_legacy_mod

! In alphabetical order of module name
use atm_fields_bounds_mod, only: pdims, tdims
use atm_fields_mod, only:                                                      &
  gc_nd_ait_sol,                                                               &
  gc_ait_sol_su,                                                               &
  gc_ait_sol_bc,                                                               &
  gc_ait_sol_oc,                                                               &
  gc_nd_acc_sol,                                                               &
  gc_acc_sol_su,                                                               &
  gc_acc_sol_bc,                                                               &
  gc_acc_sol_oc,                                                               &
  gc_acc_sol_ss,                                                               &
  gc_nd_cor_sol,                                                               &
  gc_cor_sol_su,                                                               &
  gc_cor_sol_bc,                                                               &
  gc_cor_sol_oc,                                                               &
  gc_cor_sol_ss,                                                               &
  gc_nd_ait_ins,                                                               &
  gc_ait_ins_bc,                                                               &
  gc_ait_ins_oc
use c_sulchm_mod, only: chi, rad_acc, rad_ait, sigma
use calc_surf_area_mod, only: calc_surf_area
use cderived_mod, only: delta_lambda, delta_phi
use copydiag_3d_mod, only: copydiag_3d
use copydiag_mod, only: copydiag
use cstash_mod, only: idom_b, iopl_d, isec_b, item_b, modl_b, ndiag
! Module for JULES-based atmospheric deposition
use deposition_from_ukca_chemistry_mod, only: deposition_from_ukca_chemistry
use dust_parameters_mod, only: l_dust, l_twobin_dust, drep, rhop
use get_emdiag_stash_mod, only: get_emdiag_stash
use glomap_clim_option_mod, only: i_gc_activation_arg,                         &
                                  i_glomap_clim_activation_scheme,             &
                                  i_glomap_clim_setup, l_glomap_clim_radaer,   &
                                  i_glomap_clim_tune_bc
use planet_constants_mod, only: c_virtual, cp, g, kappa, planet_radius, pref,  &
                                r, vkman
use set_pseudo_list_mod, only: set_pseudo_list
use stash_array_mod, only: len_stlist, num_stash_levels, num_stash_pseudo,     &
                           sf, si, si_last, stash_levels, stash_pseudo_levels, &
                           stindex, stlist
use atm_step_local, only: stashwork34, stashwork38, stashwork50
use stparam_mod, only: st_levels_model_theta, st_levels_single
use submodel_mod, only: submodel_for_sm, atmos_im
use trsrce_mod, only: trsrce
use um_stashcode_mod, only:                                                    &
  stashcode_bc_acc_sol,                                                        &
  stashcode_bc_acc_sol_load,                                                   &
  stashcode_bc_ait_insol,                                                      &
  stashcode_bc_ait_insol_load,                                                 &
  stashcode_bc_ait_sol,                                                        &
  stashcode_bc_ait_sol_load,                                                   &
  stashcode_bc_cor_sol,                                                        &
  stashcode_bc_cor_sol_load,                                                   &
  stashcode_bc_total_load,                                                     &
  stashcode_du_acc_insol,                                                      &
  stashcode_du_acc_insol_load,                                                 &
  stashcode_du_acc_sol,                                                        &
  stashcode_du_acc_sol_load,                                                   &
  stashcode_du_cor_insol,                                                      &
  stashcode_du_cor_insol_load,                                                 &
  stashcode_du_sup_insol,                                                      &
  stashcode_du_sup_insol_load,                                                 &
  stashcode_du_cor_sol,                                                        &
  stashcode_du_cor_sol_load,                                                   &
  stashcode_du_total_load,                                                     &
  stashcode_glomap_sec,                                                        &
  stashcode_h2o_acc_sol,                                                       &
  stashcode_h2o_acc_sol_load,                                                  &
  stashcode_h2o_ait_sol,                                                       &
  stashcode_h2o_ait_sol_load,                                                  &
  stashcode_h2o_cor_sol,                                                       &
  stashcode_h2o_cor_sol_load,                                                  &
  stashcode_h2o_mmr,                                                           &
  stashcode_h2o_nuc_sol,                                                       &
  stashcode_h2o_nuc_sol_load,                                                  &
  stashcode_h2o_total,                                                         &
  stashcode_h2o_total_load,                                                    &
  stashcode_mp_acc_sol,                                                        &
  stashcode_mp_acc_sol_load,                                                   &
  stashcode_mp_acc_insol,                                                      &
  stashcode_mp_acc_insol_load,                                                 &
  stashcode_mp_ait_sol,                                                        &
  stashcode_mp_ait_sol_load,                                                   &
  stashcode_mp_ait_insol,                                                      &
  stashcode_mp_ait_insol_load,                                                 &
  stashcode_mp_cor_sol,                                                        &
  stashcode_mp_cor_sol_load,                                                   &
  stashcode_mp_cor_insol,                                                      &
  stashcode_mp_cor_insol_load,                                                 &
  stashcode_mp_sup_insol,                                                      &
  stashcode_mp_sup_insol_load,                                                 &
  stashcode_mp_total_load,                                                     &
  stashcode_n_acc_insol,                                                       &
  stashcode_n_acc_sol,                                                         &
  stashcode_n_ait_insol,                                                       &
  stashcode_n_ait_sol,                                                         &
  stashcode_n_cor_insol,                                                       &
  stashcode_n_sup_insol,                                                       &
  stashcode_n_cor_sol,                                                         &
  stashcode_n_nuc_sol,                                                         &
  stashcode_nh4_acc_sol,                                                       &
  stashcode_nh4_acc_sol_load,                                                  &
  stashcode_nh4_ait_sol,                                                       &
  stashcode_nh4_ait_sol_load,                                                  &
  stashcode_nh4_cor_sol,                                                       &
  stashcode_nh4_cor_sol_load,                                                  &
  stashcode_nh4_total_load,                                                    &
  stashcode_nn_acc_sol,                                                        &
  stashcode_nn_acc_sol_load,                                                   &
  stashcode_nn_cor_sol,                                                        &
  stashcode_nn_cor_sol_load,                                                   &
  stashcode_nn_total_load,                                                     &
  stashcode_no3_acc_sol,                                                       &
  stashcode_no3_acc_sol_load,                                                  &
  stashcode_no3_ait_sol,                                                       &
  stashcode_no3_ait_sol_load,                                                  &
  stashcode_no3_cor_sol,                                                       &
  stashcode_no3_cor_sol_load,                                                  &
  stashcode_no3_total_load,                                                    &
  stashcode_oc_acc_sol,                                                        &
  stashcode_oc_acc_sol_load,                                                   &
  stashcode_oc_ait_insol,                                                      &
  stashcode_oc_ait_insol_load,                                                 &
  stashcode_oc_ait_sol,                                                        &
  stashcode_oc_ait_sol_load,                                                   &
  stashcode_oc_cor_sol,                                                        &
  stashcode_oc_cor_sol_load,                                                   &
  stashcode_oc_nuc_sol,                                                        &
  stashcode_oc_nuc_sol_load,                                                   &
  stashcode_oc_total_load,                                                     &
  stashcode_pm10_bc,                                                           &
  stashcode_pm10_dry,                                                          &
  stashcode_pm10_du,                                                           &
  stashcode_pm10_mp,                                                           &
  stashcode_pm10_nh4,                                                          &
  stashcode_pm10_nn,                                                           &
  stashcode_pm10_no3,                                                          &
  stashcode_pm10_oc,                                                           &
  stashcode_pm10_so4,                                                          &
  stashcode_pm10_ss,                                                           &
  stashcode_pm10_wet,                                                          &
  stashcode_pm2p5_bc,                                                          &
  stashcode_pm2p5_dry,                                                         &
  stashcode_pm2p5_du,                                                          &
  stashcode_pm2p5_mp,                                                          &
  stashcode_pm2p5_nh4,                                                         &
  stashcode_pm2p5_nn,                                                          &
  stashcode_pm2p5_no3,                                                         &
  stashcode_pm2p5_oc,                                                          &
  stashcode_pm2p5_so4,                                                         &
  stashcode_pm2p5_ss,                                                          &
  stashcode_pm2p5_wet,                                                         &
  stashcode_so4_acc_sol,                                                       &
  stashcode_so4_acc_sol_load,                                                  &
  stashcode_so4_ait_sol,                                                       &
  stashcode_so4_ait_sol_load,                                                  &
  stashcode_so4_cor_sol,                                                       &
  stashcode_so4_cor_sol_load,                                                  &
  stashcode_so4_nuc_sol,                                                       &
  stashcode_so4_nuc_sol_load,                                                  &
  stashcode_so4_total_load,                                                    &
  stashcode_ss_acc_sol,                                                        &
  stashcode_ss_acc_sol_load,                                                   &
  stashcode_ss_cor_sol,                                                        &
  stashcode_ss_cor_sol_load,                                                   &
  stashcode_ss_total_load,                                                     &
  stashcode_ukca_atmos_cfc11,                                                  &
  stashcode_ukca_atmos_cfc12,                                                  &
  stashcode_ukca_atmos_ch3br,                                                  &
  stashcode_ukca_atmos_ch4,                                                    &
  stashcode_ukca_atmos_co,                                                     &
  stashcode_ukca_atmos_h2,                                                     &
  stashcode_ukca_atmos_n2o,                                                    &
  stashcode_ukca_chem_diag,                                                    &
  stashcode_ukca_h_plus,                                                       &
  stashcode_ukca_nat,                                                          &
  stashcode_ukca_so4_sad,                                                      &
  stashcode_ukca_strat_ch4,                                                    &
  stashcode_ukca_strt_ch4_lss,                                                 &
  stashcode_ukca_trop_ch4,                                                     &
  stashcode_ukca_trop_o3,                                                      &
  stashcode_ukca_trop_oh,                                                      &
  stashcode_ukca_plumeria_height
use ukca_d1_defs, only: code, imode_first, istrat_first, item1_mode_diags,     &
                        item1_nitrate_diags, item1_nitrate_noems,              &
                        itemn_nitrate_diags,                                   &
                        item1_dust3mode_diags, itemN_dust3mode_diags,          &
                        item1_microplastic_diags, itemn_microplastic_diags,    &
                        l_ukca_mode_diags, l_ukca_plume_diags,                 &
                        l_ukca_stratflux, n_mode_diags, n_strat_fluxdiags,     &
                        nukca_d1items, ukca_diag_sect, ukca_item_sulpc,        &
                        ukcad1codes

use ukca_volcanic_so2_mod, only: ukca_volcanic_so2
use um_parcore, only: mype
use um_parvars, only: datastart
use vectlib_mod, only: cubrt_v, exp_v, log_v, oneover_v, powr_v
use umerf_mod, only: umerf

implicit none

! Flags to indicate UM code availability
logical, parameter :: l_um_infrastructure = .true.
  ! UM infrastructure code including STASH support, grid parameters
logical, parameter :: l_um_emissions_updates = .true.
  ! trsrc subroutine for emissions updating of tracers
logical, parameter :: l_um_calc_surf_area = .true.
  ! calc_surf_area subroutine

end module ukca_um_legacy_mod


