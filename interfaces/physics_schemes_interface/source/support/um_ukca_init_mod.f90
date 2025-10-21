!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief UKCA initialisation subroutine for UM science configuration

module um_ukca_init_mod

  ! LFRic namelists which have been read
  use aerosol_config_mod,        only: glomap_mode,                            &
                                       glomap_mode_ukca,                       &
                                       glomap_mode_dust_and_clim,              &
                                       emissions, emissions_GC3, emissions_GC5,&
                                       easyaerosol_cdnc, ukca_mode_seg_size
  use section_choice_config_mod, only: aerosol, aerosol_um
  use chemistry_config_mod,      only: chem_scheme, chem_scheme_offline_ox,    &
                                       chem_scheme_strattrop, chem_scheme_none,&
                                       chem_scheme_strat_test,                 &
                                       chem_scheme_flexchem,                   &
                                       l_ukca_ro2_ntp,                         &
                                       l_ukca_asad_full,                       &
                                       i_chem_timestep_halvings,               &
                                       l_ukca_quasinewton,                     &
                                       l_ukca_linox_scaling,                   &
                                       lightnox_scale_fac,                     &
                                       i_ukca_chem_version,                    &
                                       top_bdy_opt, top_bdy_opt_no_overwrt,    &
                                       top_bdy_opt_overwrt_top_two_lev,        &
                                       top_bdy_opt_overwrt_only_top_lev,       &
                                       top_bdy_opt_overwrt_co_no_o3_top,       &
                                       top_bdy_opt_overwrt_co_no_o3_h2o_top,   &
                                       fjx_solcyc_type

  ! Other LFRic modules used

  use log_mod, only : log_event,                                               &
                      log_scratch_space,                                       &
                      LOG_LEVEL_ERROR, LOG_LEVEL_INFO

  use constants_mod, only : i_def, r_um, i_um, imdi
  use water_constants_mod, only: tfs, rho_water, rhosea, latent_heat_cw => lc
  use chemistry_constants_mod, only: avogadro, boltzmann, rho_so4

  ! Callback module to interface with parent-specific routines needed
  ! by UKCA
  use lfric_ukca_callback_mod, only: bl_tracer_mix

  ! UM modules containing things that needs setting and setup routines
  use mphys_diags_mod, only: l_praut_diag, l_pracw_diag, l_piacw_diag,         &
                             l_psacw_diag
  use bl_diags_mod, only: bl_diag
  use ukca_nmspec_mod, only: nm_spec_active, nmspec_len
  use ukca_option_mod, only: i_mode_nucscav, l_ukca_plume_scav
  use ukca_scavenging_mod, only: ukca_set_conv_indices, tracer_info

  ! Other UM modules used
  use dms_flux_mod_4a,      only: i_liss_merlivat

  ! JULES modules used
  use jules_soil_mod,          only: dzsoil_io
  use c_rmol,                  only: rmol

  ! UKCA API module
  use ukca_api_mod, only: ukca_setup,                                          &
                          ukca_get_ntp_varlist,                                &
                          ukca_get_tracer_varlist,                             &
                          ukca_get_envgroup_varlists,                          &
                          ukca_get_emission_varlist,                           &
                          ukca_register_emission,                              &
                          ukca_chem_off,                                       &
                          ukca_chem_offline,                                   &
                          ukca_chem_strattrop,                                 &
                          ukca_age_reset_by_level,                             &
                          ukca_activation_arg,                                 &
                          ukca_activation_off,                                 &
                          ukca_maxlen_fieldname,                               &
                          ukca_maxlen_message,                                 &
                          ukca_maxlen_procname,                                &
                          ukca_maxlen_emiss_long_name,                         &
                          ukca_maxlen_emiss_tracer_name,                       &
                          ukca_maxlen_emiss_var_name,                          &
                          ukca_maxlen_emiss_vert_fact,                         &
                          ukca_strat_lbc_env,                                  &
                          ukca_get_photol_reaction_data,                       &
                          ukca_photol_varname_len

  ! Photolysis module
  use photol_api_mod,  only: photol_jlabel_len,                                &
    !  Max sizes for spectral file data
    photol_max_miesets, photol_n_solcyc_av, photol_sw_band_aer,                &
    photol_sw_phases, photol_max_wvl, photol_max_crossec,                      &
    photol_wvl_intervals, photol_num_tvals

  implicit none

  private
  public :: um_ukca_init

  ! Number of entrainment levels considered by UM routine tr_mix
  integer(i_um), parameter, public :: nlev_ent_tr_mix = 3

  ! Number of emission field entries required for the UKCA configuration
  integer(i_um) :: n_emiss_slots

  ! Number of emission entries for each emitted species
  integer, pointer :: n_slots(:)

  ! --------------------------------------------------------------------------
  ! List of dust_only required chemistry emissions
  ! This can be retired when GLOMAP can be run without chemistry

  character(len=ukca_maxlen_emiss_tracer_name), parameter ::   &
      required_2d_emissions(4) = [ 'DMS       ', 'Monoterp  ', &
                                   'SO2_low   ', 'SO2_high  ' ]

  character(len=ukca_maxlen_emiss_tracer_name), parameter :: &
      required_emissions(5) = [ 'DMS       ', 'Monoterp  ',  &
                                'SO2_low   ', 'SO2_high  ', 'SO2_nat   ' ]

  ! UKCA option choices

  ! UKCA tracers top boundary (topmost levels) overwrite method
  ! 0: Do not overwrite any levels
  ! 1: Overwrite top 2 levels with 3rd (except H2O)
  ! 2: Overwrite top level with level below
  ! 3: Overwrite top level of CO, NO, O3 with ACE-FTS climatology
  ! 4: Overwrite top level of CO, NO, O3, H2O with ACE-FTS climatology
  integer, parameter, public :: i_no_overwrt = 0
  integer, parameter, public :: i_overwrt_top_two_lev = 1
  integer, parameter, public :: i_overwrt_only_top_lev = 2
  integer, parameter, public :: i_overwrt_co_no_o3_top = 3
  integer, parameter, public :: i_overwrt_co_no_o3_h2o_top = 4
  ! --------------------------------------------------------------------------
  ! UKCA field name definitions follow here.
  ! All UKCA field names used in LFRic should be defined below and should be
  ! used from this module as required. This is intended to be the master copy
  ! that can be updated in the event of any changes to UKCA names in future
  ! UKCA versions.
  ! --------------------------------------------------------------------------

  ! -- UKCA field names - Tracers --
  character(len=*), parameter, public :: fldname_o3 = 'O3'
  character(len=*), parameter, public :: fldname_no3 = 'NO3'
  character(len=*), parameter, public :: fldname_oh = 'OH'
  character(len=*), parameter, public :: fldname_ho2 = 'HO2'

  character(len=*), parameter, public :: fldname_o3p = 'O(3P)'
  character(len=*), parameter, public :: fldname_o1d = 'O(1D)'
  character(len=*), parameter, public :: fldname_n = 'N'
  character(len=*), parameter, public :: fldname_no = 'NO'
  character(len=*), parameter, public :: fldname_no2 = 'NO2'
  character(len=*), parameter, public :: fldname_n2o5 = 'N2O5'
  character(len=*), parameter, public :: fldname_ho2no2 = 'HO2NO2'
  character(len=*), parameter, public :: fldname_hono2 = 'HONO2'
  character(len=*), parameter, public :: fldname_ch4 = 'CH4'
  character(len=*), parameter, public :: fldname_co = 'CO'
  character(len=*), parameter, public :: fldname_hcho = 'HCHO'
  character(len=*), parameter, public :: fldname_meoo = 'MeOO'
  character(len=*), parameter, public :: fldname_meooh = 'MeOOH'
  character(len=*), parameter, public :: fldname_h = 'H'
  character(len=*), parameter, public :: fldname_ch2o = 'H2O'
  character(len=*), parameter, public :: fldname_cl = 'Cl'
  character(len=*), parameter, public :: fldname_cl2o2 = 'Cl2O2'
  character(len=*), parameter, public :: fldname_clo = 'ClO'
  character(len=*), parameter, public :: fldname_oclo = 'OClO'
  character(len=*), parameter, public :: fldname_br = 'Br'
  character(len=*), parameter, public :: fldname_bro = 'BrO'
  character(len=*), parameter, public :: fldname_brcl = 'BrCl'
  character(len=*), parameter, public :: fldname_brono2 = 'BrONO2'
  character(len=*), parameter, public :: fldname_n2o = 'N2O'
  character(len=*), parameter, public :: fldname_hcl = 'HCl'
  character(len=*), parameter, public :: fldname_hocl = 'HOCl'
  character(len=*), parameter, public :: fldname_hbr = 'HBr'
  character(len=*), parameter, public :: fldname_hobr = 'HOBr'
  character(len=*), parameter, public :: fldname_clono2 = 'ClONO2'
  character(len=*), parameter, public :: fldname_cfcl3 = 'CFCl3'
  character(len=*), parameter, public :: fldname_cf2cl2 = 'CF2Cl2'
  character(len=*), parameter, public :: fldname_mebr = 'MeBr'
  character(len=*), parameter, public :: fldname_hono = 'HONO'
  character(len=*), parameter, public :: fldname_c2h6 = 'C2H6'
  character(len=*), parameter, public :: fldname_etoo = 'EtOO'
  character(len=*), parameter, public :: fldname_etooh = 'EtOOH'
  character(len=*), parameter, public :: fldname_mecho = 'MeCHO'
  character(len=*), parameter, public :: fldname_meco3 = 'MeCO3'
  character(len=*), parameter, public :: fldname_pan = 'PAN'
  character(len=*), parameter, public :: fldname_c3h8 = 'C3H8'
  character(len=*), parameter, public :: fldname_n_proo = 'n-PrOO'
  character(len=*), parameter, public :: fldname_i_proo = 'i-PrOO'
  character(len=*), parameter, public :: fldname_n_prooh = 'n-PrOOH'
  character(len=*), parameter, public :: fldname_i_prooh = 'i-PrOOH'
  character(len=*), parameter, public :: fldname_etcho = 'EtCHO'
  character(len=*), parameter, public :: fldname_etco3 = 'EtCO3'
  character(len=*), parameter, public :: fldname_me2co = 'Me2CO'
  character(len=*), parameter, public :: fldname_mecoch2oo = 'MeCOCH2OO'
  character(len=*), parameter, public :: fldname_mecoch2ooh = 'MeCOCH2OOH'
  character(len=*), parameter, public :: fldname_ppan = 'PPAN'
  character(len=*), parameter, public :: fldname_meono2 = 'MeONO2'
  character(len=*), parameter, public :: fldname_c5h8 = 'C5H8'
  character(len=*), parameter, public :: fldname_iso2 = 'ISO2'
  character(len=*), parameter, public :: fldname_isooh = 'ISOOH'
  character(len=*), parameter, public :: fldname_ison = 'ISON'
  character(len=*), parameter, public :: fldname_macr = 'MACR'
  character(len=*), parameter, public :: fldname_macro2 = 'MACRO2'
  character(len=*), parameter, public :: fldname_macrooh = 'MACROOH'
  character(len=*), parameter, public :: fldname_mpan = 'MPAN'
  character(len=*), parameter, public :: fldname_hacet = 'HACET'
  character(len=*), parameter, public :: fldname_mgly = 'MGLY'
  character(len=*), parameter, public :: fldname_nald = 'NALD'
  character(len=*), parameter, public :: fldname_hcooh = 'HCOOH'
  character(len=*), parameter, public :: fldname_meco3h = 'MeCO3H'
  character(len=*), parameter, public :: fldname_meco2h = 'MeCO2H'
  character(len=*), parameter, public :: fldname_h2 = 'H2'
  character(len=*), parameter, public :: fldname_meoh = 'MeOH'
  character(len=*), parameter, public :: fldname_msa = 'MSA'
  character(len=*), parameter, public :: fldname_nh3 = 'NH3'
  character(len=*), parameter, public :: fldname_cs2 = 'CS2'
  character(len=*), parameter, public :: fldname_csul = 'COS'
  character(len=*), parameter, public :: fldname_h2s = 'H2S'
  character(len=*), parameter, public :: fldname_so3 = 'SO3'
  character(len=*), parameter, public :: fldname_passive_o3 = 'PASSIVE O3'
  character(len=*), parameter, public :: fldname_age_of_air = 'AGE OF AIR'
  character(len=*), parameter, public :: fldname_co2 = 'CO2'
  character(len=*), parameter, public :: fldname_o2 = 'O2'
  character(len=*), parameter, public :: fldname_n2 = 'N2'

  ! Lumped tracers - names contain the species these are represented as
  character(len=*), parameter, public :: fldname_lumped_n = 'NO2'
  character(len=*), parameter, public :: fldname_lumped_cl = 'HCl'
  character(len=*), parameter, public :: fldname_lumped_br = 'BrO'

  character(len=*), parameter, public :: fldname_h2o2 = 'H2O2'
  character(len=*), parameter, public :: fldname_dms = 'DMS'
  character(len=*), parameter, public :: fldname_so2 = 'SO2'
  character(len=*), parameter, public :: fldname_h2so4 = 'H2SO4'
  character(len=*), parameter, public :: fldname_dmso = 'DMSO'
  character(len=*), parameter, public :: fldname_monoterpene = 'Monoterp'
  character(len=*), parameter, public :: fldname_secondary_organic = 'Sec_Org'
  character(len=*), parameter, public :: fldname_n_nuc_sol = 'Nuc_SOL_N'
  character(len=*), parameter, public :: fldname_nuc_sol_su = 'Nuc_SOL_SU'
  character(len=*), parameter, public :: fldname_nuc_sol_om = 'Nuc_SOL_OM'
  character(len=*), parameter, public :: fldname_n_ait_sol = 'Ait_SOL_N'
  character(len=*), parameter, public :: fldname_ait_sol_su = 'Ait_SOL_SU'
  character(len=*), parameter, public :: fldname_ait_sol_bc = 'Ait_SOL_BC'
  character(len=*), parameter, public :: fldname_ait_sol_om = 'Ait_SOL_OM'
  character(len=*), parameter, public :: fldname_n_acc_sol = 'Acc_SOL_N'
  character(len=*), parameter, public :: fldname_acc_sol_su = 'Acc_SOL_SU'
  character(len=*), parameter, public :: fldname_acc_sol_bc = 'Acc_SOL_BC'
  character(len=*), parameter, public :: fldname_acc_sol_om = 'Acc_SOL_OM'
  character(len=*), parameter, public :: fldname_acc_sol_ss = 'Acc_SOL_SS'
  character(len=*), parameter, public :: fldname_acc_sol_du = 'Acc_SOL_DU'
  character(len=*), parameter, public :: fldname_n_cor_sol = 'Cor_SOL_N'
  character(len=*), parameter, public :: fldname_cor_sol_su = 'Cor_SOL_SU'
  character(len=*), parameter, public :: fldname_cor_sol_bc = 'Cor_SOL_BC'
  character(len=*), parameter, public :: fldname_cor_sol_om = 'Cor_SOL_OM'
  character(len=*), parameter, public :: fldname_cor_sol_ss = 'Cor_SOL_SS'
  character(len=*), parameter, public :: fldname_cor_sol_du = 'Cor_SOL_DU'
  character(len=*), parameter, public :: fldname_n_ait_ins = 'Ait_INS_N'
  character(len=*), parameter, public :: fldname_ait_ins_bc = 'Ait_INS_BC'
  character(len=*), parameter, public :: fldname_ait_ins_om = 'Ait_INS_OM'
  character(len=*), parameter, public :: fldname_n_acc_ins = 'Acc_INS_N'
  character(len=*), parameter, public :: fldname_acc_ins_du = 'Acc_INS_DU'
  character(len=*), parameter, public :: fldname_n_cor_ins = 'Cor_INS_N'
  character(len=*), parameter, public :: fldname_cor_ins_du = 'Cor_INS_DU'

  ! -- UKCA field names - Non-transported prognostics --
  character(len=*), parameter, public :: fldname_cloud_drop_no_conc =          &
                                         'cdnc'
  character(len=*), parameter, public :: fldname_surfarea =                    &
                                         'surfarea'
  character(len=*), parameter, public :: fldname_drydp_ait_sol =               &
                                         'drydiam_ait_sol'
  character(len=*), parameter, public :: fldname_drydp_acc_sol =               &
                                         'drydiam_acc_sol'
  character(len=*), parameter, public :: fldname_drydp_cor_sol =               &
                                         'drydiam_cor_sol'
  character(len=*), parameter, public :: fldname_drydp_ait_ins =               &
                                         'drydiam_ait_insol'
  character(len=*), parameter, public :: fldname_drydp_acc_ins =               &
                                         'drydiam_acc_insol'
  character(len=*), parameter, public :: fldname_drydp_cor_ins =               &
                                         'drydiam_cor_insol'
  character(len=*), parameter, public :: fldname_wetdp_ait_sol =               &
                                         'wetdiam_ait_sol'
  character(len=*), parameter, public :: fldname_wetdp_acc_sol =               &
                                         'wetdiam_acc_sol'
  character(len=*), parameter, public :: fldname_wetdp_cor_sol =               &
                                         'wetdiam_cor_sol'
  character(len=*), parameter, public :: fldname_rhopar_ait_sol =              &
                                         'aerdens_ait_sol'
  character(len=*), parameter, public :: fldname_rhopar_acc_sol =              &
                                         'aerdens_acc_sol'
  character(len=*), parameter, public :: fldname_rhopar_cor_sol =              &
                                         'aerdens_cor_sol'
  character(len=*), parameter, public :: fldname_rhopar_ait_ins =              &
                                         'aerdens_ait_insol'
  character(len=*), parameter, public :: fldname_rhopar_acc_ins =              &
                                         'aerdens_acc_insol'
  character(len=*), parameter, public :: fldname_rhopar_cor_ins =              &
                                         'aerdens_cor_insol'
  character(len=*), parameter, public :: fldname_pvol_su_ait_sol =             &
                                         'pvol_su_ait_sol'
  character(len=*), parameter, public :: fldname_pvol_bc_ait_sol =             &
                                         'pvol_bc_ait_sol'
  character(len=*), parameter, public :: fldname_pvol_om_ait_sol =             &
                                         'pvol_oc_ait_sol'
  character(len=*), parameter, public :: fldname_pvol_wat_ait_sol =            &
                                         'pvol_h2o_ait_sol'
  character(len=*), parameter, public :: fldname_pvol_su_acc_sol =             &
                                         'pvol_su_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_bc_acc_sol =             &
                                         'pvol_bc_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_om_acc_sol =             &
                                         'pvol_oc_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_ss_acc_sol =             &
                                         'pvol_ss_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_du_acc_sol =             &
                                         'pvol_du_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_wat_acc_sol =            &
                                         'pvol_h2o_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_su_cor_sol =             &
                                         'pvol_su_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_bc_cor_sol =             &
                                         'pvol_bc_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_om_cor_sol =             &
                                         'pvol_oc_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_ss_cor_sol =             &
                                         'pvol_ss_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_du_cor_sol =             &
                                         'pvol_du_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_wat_cor_sol =            &
                                         'pvol_h2o_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_bc_ait_ins =             &
                                         'pvol_bc_ait_insol'
  character(len=*), parameter, public :: fldname_pvol_om_ait_ins =             &
                                         'pvol_oc_ait_insol'
  character(len=*), parameter, public :: fldname_pvol_du_acc_ins =             &
                                         'pvol_du_acc_insol'
  character(len=*), parameter, public :: fldname_pvol_du_cor_ins =             &
                                         'pvol_du_cor_insol'

  ! -- UKCA field names - Environmental drivers --

  ! - Drivers in scalar group -
  character(len=*), parameter, public :: fldname_sin_declination =             &
                                         'sin_declination'
  character(len=*), parameter, public :: fldname_equation_of_time =            &
                                         'equation_of_time'
  character(len=*), parameter, public :: fldname_atmos_ccl4 = 'atmospheric_ccl4'
  character(len=*), parameter, public :: fldname_atmos_cfc113 = 'atmospheric_cfc113'
  character(len=*), parameter, public :: fldname_atmos_cfc114 = 'atmospheric_cfc114'
  character(len=*), parameter, public :: fldname_atmos_cfc115 = 'atmospheric_cfc115'
  character(len=*), parameter, public :: fldname_atmos_cfc11 = 'atmospheric_cfc11'
  character(len=*), parameter, public :: fldname_atmos_cfc12 = 'atmospheric_cfc12'
  character(len=*), parameter, public :: fldname_atmos_ch2br2 = 'atmospheric_ch2br2'
  character(len=*), parameter, public :: fldname_atmos_chbr3 = 'atmospheric_chbr3'
  character(len=*), parameter, public :: fldname_atmos_ch4 = 'atmospheric_ch4'
  character(len=*), parameter, public :: fldname_atmos_co2 = 'atmospheric_co2'
  character(len=*), parameter, public :: fldname_atmos_csul = 'atmospheric_cos'
  character(len=*), parameter, public :: fldname_atmos_h1202 = 'atmospheric_h1202'
  character(len=*), parameter, public :: fldname_atmos_h1211 = 'atmospheric_h1211'
  character(len=*), parameter, public :: fldname_atmos_h1301 = 'atmospheric_h1301'
  character(len=*), parameter, public :: fldname_atmos_h2 = 'atmospheric_h2'
  character(len=*), parameter, public :: fldname_atmos_h2402 = 'atmospheric_h2402'
  character(len=*), parameter, public :: fldname_atmos_hfc125 = 'atmospheric_hfc125'
  character(len=*), parameter, public :: fldname_atmos_hfc134a =               &
                                         'atmospheric_hfc134a'
  character(len=*), parameter, public :: fldname_atmos_hcfc141b =              &
                                         'atmospheric_hcfc141b'
  character(len=*), parameter, public :: fldname_atmos_hcfc142b =              &
                                         'atmospheric_hcfc142b'
  character(len=*), parameter, public :: fldname_atmos_hcfc22 = 'atmospheric_hcfc22'
  character(len=*), parameter, public :: fldname_atmos_mebr = 'atmospheric_mebr'
  character(len=*), parameter, public :: fldname_atmos_meccl3 = 'atmospheric_meccl3'
  character(len=*), parameter, public :: fldname_atmos_mecl = 'atmospheric_mecl'
  character(len=*), parameter, public :: fldname_atmos_n2 = 'atmospheric_n2'
  character(len=*), parameter, public :: fldname_atmos_n2o = 'atmospheric_n2o'
  character(len=*), parameter, public :: fldname_atmos_o2 = 'atmospheric_o2'

  ! - Drivers in flat grid groups (integer, real & logical) -

  ! Photolysis & emissions-related drivers (integer)
  character(len=*), parameter, public :: fldname_kent = 'kent'
  character(len=*), parameter, public :: fldname_kent_dsc = 'kent_dsc'
  character(len=*), parameter, public :: fldname_cv_base = 'conv_cloud_base'
  character(len=*), parameter, public :: fldname_cv_top = 'conv_cloud_top'
  character(len=*), parameter, public :: fldname_cv_cloud_lwp = 'conv_cloud_lwp'

  ! General purpose drivers (real)
  character(len=*), parameter, public :: fldname_latitude = 'latitude'
  character(len=*), parameter, public :: fldname_longitude = 'longitude'
  character(len=*), parameter, public :: fldname_sin_latitude = 'sin_latitude'
  character(len=*), parameter, public :: fldname_cos_latitude = 'cos_latitude'
  character(len=*), parameter, public :: fldname_tan_latitude = 'tan_latitude'
  character(len=*), parameter, public :: fldname_frac_seaice = 'seaice_frac'
  character(len=*), parameter, public :: fldname_tstar = 'tstar'
  character(len=*), parameter, public :: fldname_pstar = 'pstar'
  character(len=*), parameter, public :: fldname_rough_length = 'rough_length'
  character(len=*), parameter, public :: fldname_ustar = 'u_s'
  character(len=*), parameter, public :: fldname_surf_hf = 'surf_hf'
  character(len=*), parameter, public :: fldname_zbl = 'zbl'
  character(len=*), parameter, public :: fldname_surf_albedo  = 'surf_albedo'
  character(len=*), parameter, public :: fldname_grid_surf_area =              &
                                         'grid_surf_area'
  ! Emissions-related drivers (real)
  character(len=*), parameter, public :: fldname_u_scalar_10m =                &
                                         'u_scalar_10m'
  character(len=*), parameter, public :: fldname_chloro_sea =                  &
                                         'chloro_sea'
  character(len=*), parameter, public :: fldname_dms_sea_conc =                &
                                         'dms_sea_conc'
  character(len=*), parameter, public :: fldname_dust_flux_div1 =              &
                                         'dust_flux_div1'
  character(len=*), parameter, public :: fldname_dust_flux_div2 =              &
                                         'dust_flux_div2'
  character(len=*), parameter, public :: fldname_dust_flux_div3 =              &
                                         'dust_flux_div3'
  character(len=*), parameter, public :: fldname_dust_flux_div4 =              &
                                         'dust_flux_div4'
  character(len=*), parameter, public :: fldname_dust_flux_div5 =              &
                                         'dust_flux_div5'
  character(len=*), parameter, public :: fldname_dust_flux_div6 =              &
                                         'dust_flux_div6'
  character(len=*), parameter, public :: fldname_zhsc =                        &
                                         'zhsc'
  character(len=*), parameter, public :: fldname_surf_wetness =                &
                                         'surf_wetness'
  ! General purpose drivers (logical)
  character(len=*), parameter, public :: fldname_l_land = 'land_sea_mask'

  ! - Drivers in flat grid plant functional type tile group -

  character(len=*), parameter, public :: fldname_stcon = 'stcon'

  ! - Drivers in full-height grid group -

  ! General purpose drivers
  character(len=*), parameter, public :: fldname_theta = 'theta'
  character(len=*), parameter, public :: fldname_exner_theta_lev =             &
                                         'exner_theta_levels'
  character(len=*), parameter, public :: fldname_p_theta_lev = 'p_theta_levels'
  character(len=*), parameter, public :: fldname_p_rho_lev = 'p_rho_levels'
  character(len=*), parameter, public :: fldname_q = 'q'
  character(len=*), parameter, public :: fldname_qcl = 'qcl'
  character(len=*), parameter, public :: fldname_qcf = 'qcf'
  character(len=*), parameter, public :: fldname_bulk_cloud_frac = 'cloud_frac'
  character(len=*), parameter, public :: fldname_ls_rain3d = 'ls_rain3d'
  character(len=*), parameter, public :: fldname_ls_snow3d = 'ls_snow3d'
  character(len=*), parameter, public :: fldname_conv_rain3d = 'conv_rain3d'
  character(len=*), parameter, public :: fldname_conv_snow3d = 'conv_snow3d'
  character(len=*), parameter, public :: fldname_rho_r2 = 'rho_r2'
  character(len=*), parameter, public :: fldname_area_cloud_frac =             &
                                         'area_cloud_fraction'
  character(len=*), parameter, public :: fldname_conv_cloud_amount =           &
                                         'conv_cloud_amount'
  character(len=*), parameter, public :: fldname_grid_volume = 'grid_volume'
  character(len=*), parameter, public :: fldname_grid_airmass = 'grid_airmass'
  character(len=*), parameter, public :: fldname_rel_humid = 'rel_humid_frac'
  ! Offline chemical fields
  character(len=*), parameter, public :: fldname_o3_offline = 'O3'
  character(len=*), parameter, public :: fldname_no3_offline = 'NO3'
  character(len=*), parameter, public :: fldname_oh_offline = 'OH'
  character(len=*), parameter, public :: fldname_ho2_offline = 'HO2'
  character(len=*), parameter, public :: fldname_h2o2_limit = 'H2O2'
  ! GLOMAP-specific drivers
  character(len=*), parameter, public :: fldname_liq_cloud_frac =              &
                                         'cloud_liq_frac'
  character(len=*), parameter, public :: fldname_autoconv = 'autoconv'
  character(len=*), parameter, public :: fldname_accretion = 'accretion'
  character(len=*), parameter, public :: fldname_rim_cry = 'rim_cry'
  character(len=*), parameter, public :: fldname_rim_agg = 'rim_agg'
  character(len=*), parameter, public :: fldname_rel_humid_clear_sky =         &
                                         'rel_humid_frac_clr'
  ! Activate-specific drivers
  character(len=*), parameter, public :: fldname_vertvel = 'vertvel'
  character(len=*), parameter, public :: fldname_svp = 'qsvp'

  ! - Drivers in full-height plus level 0 grid group -

  character(len=*), parameter, public :: fldname_interf_z = 'interf_z'

  ! - Drivers in full-height plus one grid group -

  character(len=*), parameter, public :: fldname_exner_rho_lev =               &
                                         'exner_rho_levels'

  ! - Drivers in boundary levels group -

  ! General purpose drivers
  character(len=*), parameter, public :: fldname_rhokh_rdz = 'rhokh_rdz'
  character(len=*), parameter, public :: fldname_dtrdz = 'dtrdz'
  ! Activate-specific drivers
  character(len=*), parameter, public :: fldname_bl_tke = 'bl_tke'

  ! - Drivers in entrainment levels group -

  character(len=*), parameter, public :: fldname_we_lim = 'we_lim'
  character(len=*), parameter, public :: fldname_t_frac = 't_frac'
  character(len=*), parameter, public :: fldname_zrzi = 'zrzi'
  character(len=*), parameter, public :: fldname_we_lim_dsc = 'we_lim_dsc'
  character(len=*), parameter, public :: fldname_t_frac_dsc = 't_frac_dsc'
  character(len=*), parameter, public :: fldname_zrzi_dsc = 'zrzi_dsc'

  ! - Drivers in land point group -

  character(len=*), parameter, public :: fldname_frac_land = 'fland'
  character(len=*), parameter, public :: fldname_soil_moisture_layer1 =        &
                                         'soil_moisture_layer1'

  ! - Drivers in land-point tile groups -

  character(len=*), parameter, public :: fldname_frac_surft = 'frac_types'
  character(len=*), parameter, public :: fldname_tstar_surft = 'tstar_tile'
  character(len=*), parameter, public :: fldname_z0_surft = 'z0tile_lp'
  character(len=*), parameter, public :: fldname_l_active_surft =              &
                                         'l_tile_active'

  ! - Drivers in land-point plant functional type tile group -

  character(len=*), parameter, public :: fldname_lai_pft = 'laift_lp'
  character(len=*), parameter, public :: fldname_canht_pft = 'canhtft_lp'

  ! - Photolysis rates - 4-D real
  character(len=*), parameter, public :: fldname_photol_rates = 'photol_rates'

  ! ------------------------------------------
  ! Pointers for accessing UKCA variable lists
  ! ------------------------------------------

  ! List of tracers required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public :: tracer_names(:)

  ! List of non-transported prognostics required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public :: ntp_names(:)

  ! List of species involved in photolytic reactions
  character(len=ukca_photol_varname_len), pointer, public :: ratj_varnames(:)
  ! Number of photolytic species ('jppj' elsewhere)
  integer(kind=i_um), save, public :: n_phot_spc

  ! Photolysis reaction data from UKCA
  character(len=10), save, pointer,public :: ratj_data(:,:)
  ! Variables for holding Photolysis / FastJX spectral data
  !- Fields containing spectral cross-section data (jvspec_file/ fjx_rd_xxx)
  integer(kind=i_um) :: njval               ! No of species to read x-sections for
  integer(kind=i_um) :: nw1, nw2            ! Max and min wavelength bins
  real(kind=r_um), allocatable :: fl(:)     ! TOA solar flux
  real(kind=r_um), allocatable :: q1d(:,:)  ! Photol rates for O(1D) at 3 temperatures
  real(kind=r_um), allocatable :: qo2(:,:)  ! Photol rates for O2 at 3 temperatures
  real(kind=r_um), allocatable :: qo3(:,:)  ! Photol rates for O3 at 3 temperatures
  real(kind=r_um), allocatable :: qqq(:,:,:) ! Photol rates for all other species
  real(kind=r_um), allocatable :: qrayl(:)  ! Rayleigh parameters
  real(kind=r_um), allocatable :: tqq(:,:)  ! Temperature values corresponding to rates
  real(kind=r_um), allocatable :: wl(:)     ! Effective wavelengths

  !- Fields containing (Mie) scattering parameters (jscat_file/ fjx_rd_mie)
  integer(kind=i_um) :: jtaumx              ! Max number of cloud sub layers
  integer(kind=i_um) :: naa                 ! Number of aerosol/ cloud data types
  real(kind=r_um) :: atau                   ! Cloud sub-layer factor
  real(kind=r_um) :: atau0                  ! min dtau
  real(kind=r_um), allocatable :: daa(:)    ! Density of scattering type
  real(kind=r_um), allocatable :: paa(:,:,:)  ! Phases of scattering types
  real(kind=r_um), allocatable :: qaa(:,:)  ! Q of scattering types
  real(kind=r_um), allocatable :: raa(:)    ! Effective radius of scattering type
  real(kind=r_um), allocatable :: saa(:,:)  ! Single Scattering Albedos
  real(kind=r_um), allocatable :: waa(:,:)  ! Wavelengths for scattering coefficients

  !- Fields containing Solar cycle information (solcyc_file)
  integer(kind=i_um), parameter :: n_solcyc_ts = 1476 ! Total months in solar cycle data
  real(kind=r_um), allocatable :: solcyc_av(:)  ! Average solar cycle
  real(kind=r_um), allocatable :: solcyc_quanta(:) ! Quanta component of solar cycle
  real(kind=r_um), allocatable :: solcyc_ts(:)  ! Obs. time series of solar cycle
  real(kind=r_um), allocatable :: solcyc_spec(:)  ! Spectral component of solar cyle

  integer(kind=i_um), allocatable :: jind(:)     ! Index of species from files
  character(len=photol_jlabel_len), allocatable :: jlabel(:)  ! Copy of species
                                    ! names to match those from files
  real(kind=r_um), pointer :: ratj_jfacta(:)   ! Quantum yield                                    
  real(kind=r_um), allocatable :: jfacta(:)    ! copy of quantum yield in correct units

  character(len=photol_jlabel_len), allocatable :: titlej(:)
                                    ! Names of species as read from the spec file

  ! Lists of environmental driver fields required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_scalar_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_integer(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_logical(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flatpft_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullht_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullht0_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullhtp1_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_bllev_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_entlev_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_land_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landtile_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landtile_logical(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landpft_real(:)
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullhtphot_real(:)

  ! Lists of emissions for the UKCA configuration
  character(len=ukca_maxlen_emiss_tracer_name), pointer ::                     &
    emiss_names(:)
  character(len=ukca_maxlen_emiss_var_name), allocatable, public ::            &
    emiss_names_flat(:)
  character(len=ukca_maxlen_emiss_var_name), allocatable, public ::            &
    emiss_names_fullht(:)

contains

  subroutine um_ukca_init(ncells_ukca)
  !> @brief Set up the UKCA model

  use nlsizes_namelist_mod, only: bl_levels, rows, model_levels

  use jules_surface_types_mod, only: ntype, npft,                              &
                                     brd_leaf, ndl_leaf,                       &
                                     c3_grass, c4_grass,                       &
                                     shrub, urban, lake, soil, ice

  ! UM modules used

  use timestep_mod,         only: timestep
  use cv_run_mod,           only: l_param_conv

  implicit none

    integer(i_def), intent(in) :: ncells_ukca

    integer(i_um) :: row_length_ukca

    ! Nullify pointers declared globally in this module
    nullify(n_slots)
    nullify(tracer_names)
    nullify(ntp_names)
    nullify(env_names_scalar_real)
    nullify(env_names_flat_integer)
    nullify(env_names_flat_real)
    nullify(env_names_flat_logical)
    nullify(env_names_flatpft_real)
    nullify(env_names_fullht_real)
    nullify(env_names_fullht0_real)
    nullify(env_names_fullhtp1_real)
    nullify(env_names_bllev_real)
    nullify(env_names_entlev_real)
    nullify(env_names_land_real)
    nullify(env_names_landtile_real)
    nullify(env_names_landtile_logical)
    nullify(env_names_landpft_real)
    nullify(env_names_fullhtphot_real)
    nullify(emiss_names)
    nullify(ratj_varnames)
    nullify(ratj_data)

    row_length_ukca = int( ncells_ukca, i_um )

    if ( ( aerosol == aerosol_um ) .and.                                       &
         ( glomap_mode == glomap_mode_dust_and_clim ) ) then

        call aerosol_ukca_dust_only_init( row_length_ukca, rows, model_levels, &
                                          bl_levels, timestep, l_param_conv )

    else if ( ( aerosol == aerosol_um .and.                                    &
                glomap_mode == glomap_mode_ukca ) .or.                         &
              ( chem_scheme == chem_scheme_offline_ox .or.                     &
                chem_scheme == chem_scheme_strat_test .or.                     &
                chem_scheme == chem_scheme_strattrop ) ) then

        call ukca_init( row_length_ukca, rows, model_levels, bl_levels,        &
                        ntype, npft, brd_leaf, ndl_leaf, c3_grass, c4_grass,   &
                        shrub, urban, lake, soil, ice,                         &
                        dzsoil_io(1), timestep, l_param_conv )

    end if

  end subroutine um_ukca_init


  !> @brief Set up UKCA for combination of Chemistry and Aerosol schemes
  !> @details Configure UKCA using options generally consistent with a GA9 or
  !>          UKESM run. Where possible, all values are taken from GA9 or UKESM
  !>          science settings.
  !>          Also register each emission field to be supplied to UKCA and
  !>          obtain lists of UKCA tracers, non-transported prognostics and
  !>          environmental drivers required for the selected configuration.
  !> @param[in] row_length      X dimension of model grid
  !> @param[in] rows            Y dimension of model grid
  !> @param[in] model_levels    Z dimension of model grid
  !> @param[in] bl_levels       No. of Z levels in boundary layer
  !> @param[in] ntype           No. of surface types considered in interactive dry deposition
  !> @param[in] npft            No. of plant functional types
  !> @param[in] i_brd_leaf      Index of surface type 'broad-leaf tree'
  !> @param[in] i_ndl_leaf      Index of surface type 'needle-leaf tree'
  !> @param[in] i_c3_grass      Index of surface type 'c3 grass'
  !> @param[in] i_c4_grass      Index of surface type 'c4 grass'
  !> @param[in] i_shrub         Index of surface type 'shrub'
  !> @param[in] i_urban         Index of surface type 'urban'
  !> @param[in] i_lake          Index of surface type 'lake'
  !> @param[in] i_soil          Index of surface type 'soil'
  !> @param[in] i_ice           Index of surface type 'ice'
  !> @param[in] dzsoil_layer1   Thickness of surface soil layer (m)
  !> @param[in] timestep        Model time step (s)
  !> @param[in] l_param_conv    True if convection is parameterized

  subroutine ukca_init( row_length, rows, model_levels, bl_levels,             &
                        ntype, npft, i_brd_leaf, i_ndl_leaf, i_c3_grass,       &
                        i_c4_grass, i_shrub, i_urban, i_lake, i_soil, i_ice,   &
                        dzsoil_layer1, timestep, l_param_conv )

    use ukca_mode_setup, only: i_ukca_bc_tuned
    use ukca_photol_param_mod, only: jppj
    use fastjx_inphot_mod, only: fastjx_inphot

    implicit none

    integer, intent(in) :: row_length      ! X dimension of model grid
    integer, intent(in) :: rows            ! Y dimension of model grid
    integer, intent(in) :: model_levels    ! Z dimension of model grid
    integer, intent(in) :: bl_levels       ! No. of Z levels in boundary layer
    integer, intent(in) :: ntype           ! No. of surface types considered in
                                           ! interactive dry deposition
    integer, intent(in) :: npft            ! No. of plant functional types
    integer, intent(in) :: i_brd_leaf      ! Index of type 'broad-leaf tree'
    integer, intent(in) :: i_ndl_leaf      ! Index of type 'needle-leaf tree'
    integer, intent(in) :: i_c3_grass      ! Index of type 'c3 grass'
    integer, intent(in) :: i_c4_grass      ! Index of type 'c4 grass'
    integer, intent(in) :: i_shrub         ! Index of type 'shrub'
    integer, intent(in) :: i_urban         ! Index of type 'urban'
    integer, intent(in) :: i_lake          ! Index of type 'lake'
    integer, intent(in) :: i_soil          ! Index of type 'soil'
    integer, intent(in) :: i_ice           ! Index of type 'ice'
    real(r_um), intent(in) :: dzsoil_layer1! Thickness of surface soil layer (m)
    real(r_um), intent(in) :: timestep     ! Model time step (s)
    logical, intent(in) :: l_param_conv    ! True if convection is parameterized

    ! Local variables

    integer :: n
    integer :: i

    ! List of photolysis file/variable names required for photolysis
    ! calculations
    character(len=10), save, pointer :: ratj_data(:,:)
    
    character(len=*), parameter  :: emiss_units = 'kg m-2 s-1'

    ! Local copies of ukca configuration variables to allow setting up chemistry
    ! and aerosol schemes through single call to ukca_setup.
    ! These should be obtained from namelists (chemistry_config) eventually.
    ! Some of these will not be considered if a top-level switch is turned off.

    integer :: i_tmp_ukca_chem=ukca_chem_off
    integer :: i_tmp_ukca_activation_scheme
    integer :: i_ukca_top_boundary_opt = i_no_overwrt

    logical :: l_ukca_ageair = .false.
    logical :: l_ukca_set_trace_gases = .false.
    logical :: l_chem_environ_gas_scalars = .false.
    logical :: l_ukca_photolysis = .false.
    logical :: l_ukca_prescribech4 = .false.
    logical :: l_ukca_chem_aero = .false.
    logical :: l_ukca_ro2_perm = .false.
    logical :: l_ukca_mode = .false.
    logical :: l_use_gridbox_mass = .false.

    ! Some default values - either standard values in UM, or the only ones
    ! currently supported in the LFRic-side implementation.

    integer :: i_ukca_light_param=1            ! Internal Price-Rind scheme
    integer :: i_ukca_quasinewton_start=2, i_ukca_quasinewton_end=3
    integer :: i_ukca_scenario=ukca_strat_lbc_env
    integer::  i_ukca_mode_seg_size=4          ! GLOMAP-mode segment size
    real(r_um) :: linox_scale_in    
    
    character(len=ukca_photol_varname_len) :: adjusted_fname  ! intermediate spc/ filename copy
    
    ! Variables for UKCA error handling
    integer :: ukca_errcode
    character(len=ukca_maxlen_message) :: ukca_errmsg
    character(len=ukca_maxlen_procname) :: ukca_errproc

    ! Set up GAL configuration based on GA9 or UKESM1.
    ! The ASAD Newton-Raphson Offline Oxidants scheme (ukca_chem_offline)
    ! is substituted for the Explicit backward-Euler scheme used in GA9
    ! which cannot be called by columns. Hence, configuration values for
    ! nrsteps and l_ukca_asad_columns are included.
    ! Other configuration values, with the exception of temporary logicals
    ! specifiying fixes, are set to match GA9/ UKESM science settings or taken
    ! from the LFRic context. Unlike GA9, all fixes that are controlled by
    ! temporary logicals will be on. (i.e. the defaults for these
    ! logicals, .true. by convention in UKCA, are not overridden.)

    ! Initialise values depending on configuration/ scheme
    ! The strat_test option runs Offline_ox scheme but with selective parts
    ! of strattrop scheme to be turned On for development.
    ! Return if 'no chemistry' option is chosen, as UKCA requires at least
    ! one option to be active.
    if ( chem_scheme == chem_scheme_strattrop ) then
       l_ukca_photolysis = .true.
       i_tmp_ukca_chem = ukca_chem_strattrop
       l_ukca_set_trace_gases = .true.
       l_chem_environ_gas_scalars = .true.
       l_ukca_prescribech4 = .true.
       l_ukca_ro2_perm = l_ukca_ro2_ntp     ! Reduces no. of transported fields
                         ! ro2_perm usually 'On' when ro2_ntp is 'On'.

       ! Top boundary overwrite method - convert to integer option
       select case( top_bdy_opt )
         case ( top_bdy_opt_no_overwrt )
           i_ukca_top_boundary_opt = i_no_overwrt
         case ( top_bdy_opt_overwrt_top_two_lev )
           i_ukca_top_boundary_opt =  i_overwrt_top_two_lev
         case ( top_bdy_opt_overwrt_only_top_lev )
           i_ukca_top_boundary_opt = i_overwrt_only_top_lev
         case ( top_bdy_opt_overwrt_co_no_o3_top )
           i_ukca_top_boundary_opt = i_overwrt_co_no_o3_top
         case ( top_bdy_opt_overwrt_co_no_o3_h2o_top )
           i_ukca_top_boundary_opt = i_overwrt_co_no_o3_h2o_top
         case default
           call log_event('Unknown option - UKCA tracer top boundary handling', &
                           LOG_LEVEL_ERROR)
       end select

    else if ( chem_scheme == chem_scheme_offline_ox .or.  &
              chem_scheme == chem_scheme_strat_test ) then
       i_tmp_ukca_chem = ukca_chem_offline
    else     ! chem_scheme_none or chem_scheme_flexchem
      call log_event('No Chemical scheme chosen for UKCA', LOG_LEVEL_INFO)
      return
    end if
    if (aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca) then
      l_ukca_chem_aero = .true.
      l_ukca_mode = .true.
    end if

    ! If the Easy Aerosol climatology is being used to set CDNC values
    ! then there is no need to calculate CDNCs via an activation scheme in UKCA
    if (easyaerosol_cdnc) then
      i_tmp_ukca_activation_scheme = ukca_activation_off
    else
      i_tmp_ukca_activation_scheme = ukca_activation_arg
    end if

    ! Logical to determine whether UKCA expects the mass of air in gridbox from
    ! parent model. For aerosols, UKCA contains a simplified calculation of a
    ! 'relative mass' sufficient for deposition and sedimentation processes
    ! and so does not require the parent to provide this field (=false)
    if ( chem_scheme == chem_scheme_strattrop .or.  &
         chem_scheme == chem_scheme_strat_test ) then
      l_use_gridbox_mass = .true.
    end if

    ! Set default lightning NOx scale factor to 1.0 if no factor supplied
    if ( l_ukca_linox_scaling ) then
      linox_scale_in = lightnox_scale_fac
    else
      linox_scale_in = 1.0_r_um
    end if

    ! Set default GLOMAP segment size if no factor supplied
    if ( ukca_mode_seg_size /= imdi ) THEN
      i_ukca_mode_seg_size = ukca_mode_seg_size 
    end if

    call ukca_setup( ukca_errcode,                                             &
           ! Context information
           row_length=row_length,                                              &
           rows=rows,                                                          &
           model_levels=model_levels,                                          &
           bl_levels=bl_levels,                                                &
           nlev_ent_tr_mix=nlev_ent_tr_mix,                                    &
           ntype=ntype,                                                        &
           npft=npft,                                                          &
           i_brd_leaf=i_brd_leaf,                                              &
           i_ndl_leaf=i_ndl_leaf,                                              &
           i_c3_grass=i_c3_grass,                                              &
           i_c4_grass=i_c4_grass,                                              &
           i_shrub=i_shrub,                                                    &
           i_urban=i_urban,                                                    &
           i_lake=i_lake,                                                      &
           i_soil=i_soil,                                                      &
           i_ice=i_ice,                                                        &
           dzsoil_layer1=dzsoil_layer1,                                        &
           timestep=timestep,                                                  &
           ! General UKCA configuration options
           i_ukca_chem=i_tmp_ukca_chem,                                        &
           l_ukca_chem_aero=l_ukca_chem_aero,                                  &
           l_ukca_mode=l_ukca_mode,                                            &
           l_fix_tropopause_level=.true.,                                      &
           l_ukca_persist_off=.true.,                                          &
           l_ukca_ageair=l_ukca_ageair,                                        &
           ! Chemistry configuration options
           i_ukca_chem_version=i_ukca_chem_version,                            &
           chem_timestep=3600,                                                 &
           i_chem_timestep_halvings=i_chem_timestep_halvings,                  &
           nrsteps=45,                                                         &
           l_ukca_asad_columns=.true.,                                         &
           l_ukca_asad_full=l_ukca_asad_full,                                  &
           l_ukca_intdd=.true.,                                                &
           l_ukca_ddep_lev1=.false.,                                           &
           l_ukca_ddepo3_ocean=.false.,                                        &
           l_ukca_dry_dep_so2wet=.true.,                                       &
           l_ukca_quasinewton=l_ukca_quasinewton,                              &
           i_ukca_quasinewton_start=i_ukca_quasinewton_start,                  &
           i_ukca_quasinewton_end=i_ukca_quasinewton_end,                      &
           l_use_photolysis=l_ukca_photolysis,                                 &
           i_ukca_topboundary=i_ukca_top_boundary_opt,                         &
           i_ukca_light_param=i_ukca_light_param,                              &
           l_ukca_linox_scaling=l_ukca_linox_scaling,                          &
           lightnox_scale_fac=linox_scale_in,                                  &
           l_chem_environ_gas_scalars=l_chem_environ_gas_scalars,              &
           l_ukca_prescribech4=l_ukca_prescribech4,                            &
           i_strat_lbc_source=i_ukca_scenario,                                 &
           l_ukca_ro2_ntp = l_ukca_ro2_ntp,                                    &
           l_ukca_ro2_perm = l_ukca_ro2_perm,                                  &
           l_use_gridbox_mass= l_use_gridbox_mass,                             &
           ! UKCA emissions configuration options
           mode_parfrac=2.5_r_um,                                              &
           l_ukca_enable_seadms_ems=.true.,                                    &
           i_ukca_dms_flux=i_liss_merlivat,                                    &
           l_ukca_scale_seadms_ems=.false.,                                    &
           l_ukca_scale_soa_yield_mt=.true.,                                   &
           soa_yield_scaling_mt=2.0_r_um,                                      &
           l_ukca_scale_soa_yield_isop=.false.,                                &
           soa_yield_scaling_isop=1.0_r_um,                                    &
           l_support_ems_vertprof=.true.,                                      &
           i_primss_method=2,                                                  &
           ! UKCA environmental driver configuration options
           l_param_conv=l_param_conv,                                          &
           l_ctile=.true.,                                                     &
           l_environ_rel_humid=.true.,                                         &
           ! General GLOMAP configuration options
           i_mode_nzts=15,                                                     &
           i_mode_setup=8,                                                     &
           l_mode_bhn_on=.true.,                                               &
           l_mode_bln_on=.false.,                                              &
           i_mode_nucscav=i_mode_nucscav,                                      &
           mode_activation_dryr=37.5_r_um,                                     &
           mode_incld_so2_rfrac=0.25_r_um,                                     &
           ukca_mode_seg_size=i_ukca_mode_seg_size,                            &
           l_cv_rainout=.not.(l_ukca_plume_scav),                              &
           l_dust_mp_slinn_impc_scav=.true.,                                   &
           l_dust_mp_ageing=.false.,                                           &
           ! GLOMAP emissions configuration options
           l_ukca_primsu=.true.,                                               &
           l_ukca_primss=.true.,                                               &
           l_ukca_primdu=.true.,                                               &
           l_ukca_primbcoc=.true.,                                             &
           l_ukca_prim_moc=.true.,                                             &
           l_bcoc_bf=.true.,                                                   &
           l_bcoc_bm=.true.,                                                   &
           l_bcoc_ff=.true.,                                                   &
           l_ukca_scale_biom_aer_ems=.true.,                                   &
           biom_aer_ems_scaling=2.0_r_um,                                      &
           ! GLOMAP feedback configuration options
           l_ukca_radaer=.true.,                                               &
           i_ukca_tune_bc=i_ukca_bc_tuned,                                     &
           i_ukca_activation_scheme=i_tmp_ukca_activation_scheme,              &
           i_ukca_nwbins=20,                                                   &
           ! Constants
           const_rmol=rmol,                                                    &
           const_tfs=tfs,                                                      &
           const_rho_water=rho_water,                                          &
           const_rhosea=rhosea,                                                &
           const_lc=latent_heat_cw,                                            &
           const_avogadro=avogadro,                                            &
           const_boltzmann=boltzmann,                                          &
           const_rho_so4=rho_so4,                                              &
           ! Callback procedures
           proc_bl_tracer_mix = bl_tracer_mix,                                 &
           ! UKCA temporary logicals
           l_fix_ukca_hygroscopicities=.false.,                                &
           l_fix_ukca_water_content=.true.,                                    &
           ! Return status information
           error_message=ukca_errmsg,                                          &
           error_routine=ukca_errproc)

    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Retrieve the lists of required fields for the configuration
    call set_ukca_field_lists()

    ! Set up indexing data needed for plume scavenging of UKCA tracers in the
    ! GR convection scheme
    if (l_ukca_plume_scav) then
      n = size(tracer_names)
      allocate(nm_spec_active(n))
      do i = 1, n
        nm_spec_active(i) = tracer_names(i)(1:nmspec_len)
      end do
      call ukca_set_conv_indices()
      tracer_info%i_ukca_first = 1
      tracer_info%i_ukca_last = n
    end if

    ! Register emissions required for this run
    call ukca_emiss_init()

    ! Photolysis initialisation
    ! Obtain and store photolysis reactions information, ensuring that
    ! number of species match those expected by photol_param_mod
    if ( l_ukca_photolysis .and. (chem_scheme==chem_scheme_strattrop .or.      &
      chem_scheme==chem_scheme_strat_test)) then
      call ukca_get_photol_reaction_data(ratj_data, ratj_varnames,             &
                                        jfacta_ptr=ratj_jfacta)
      n_phot_spc = SIZE(ratj_varnames)

      if ( n_phot_spc /= jppj ) then
        write( log_scratch_space, '(A,2I6,A)' )                                &
          'Mismatch in expected and registered photolysis reactions: ',        &
           n_phot_spc, jppj,'. Check definitions in ukca_photol_param_mod'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if    

      ! Read spectral data files, allocate arrays and set up
      call allocate_fastjx_filevars()

      adjusted_fname = ' '
      if (.not. allocated(jlabel)) allocate(jlabel(n_phot_spc))
      jlabel(:) = ''
      if (.not. allocated(jfacta)) allocate(jfacta(n_phot_spc))
      jfacta(:) = 0.0
      if (.not. allocated(jind)) allocate(jind(n_phot_spc))

      do i = 1, n_phot_spc
        jfacta(i)=ratj_jfacta(i)/100.0e0_r_um
        adjusted_fname=trim(adjustl(ratj_varnames(i)))
        jlabel(i)=adjusted_fname(1:photol_jlabel_len)
        write(log_scratch_space,'(A,I6,E12.3,A12)')'FJX_JFACTA ', i,jfacta(i),jlabel(i)
        call log_event(log_scratch_space, LOG_LEVEL_INFO)        
      end do

     ! call wrapper routine that reads FastJX spectral and solar cycle data
     call fastjx_inphot(                                                       &
                 ! (Max) data dimensions in files
                 photol_max_miesets, photol_n_solcyc_av, photol_sw_band_aer,   &
                 photol_sw_phases, photol_max_wvl, photol_wvl_intervals,       &
                 photol_max_crossec, photol_num_tvals, n_phot_spc,             &
                 ! Variables used for/ set from cross-section data
                 njval, nw1, nw2, fl, q1d, qo2, qo3, qqq, qrayl, tqq, wl,      &
                 titlej, jlabel, jfacta, jind,                                 &
                 ! Variables used for/ set from scatterrer data
                 jtaumx, naa, atau, atau0, daa, paa, qaa, raa, saa, waa,       &
                 ! Variables used for/ set from solar cycle data
                 n_solcyc_ts, solcyc_av, solcyc_quanta, solcyc_ts,         &
                 solcyc_spec )

      ! Deallocate fastjx spectral data as not currently used - placeholder for
      ! future API implementation
      call deallocate_fastjx_filevars()
    end if   ! Photolysis

    ! Switch on optional UM microphysics diagnostics required by UKCA
    if (any(env_names_fullht_real(:) == fldname_autoconv))                     &
      l_praut_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_accretion))                    &
      l_pracw_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_rim_cry))                      &
      l_piacw_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_rim_agg))                      &
      l_psacw_diag = .true.

    ! Switch on optional UM boundary layer diagnostics required by UKCA
    if (any(env_names_bllev_real(:) == fldname_bl_tke))                        &
      bl_diag%l_request_tke = .true.

  end subroutine ukca_init
  subroutine ukca_emiss_init()
  !> @brief Set up the emissions for UKCA model
    implicit none
    ! Emissions registration:
    ! Multiple emission fields can be provided for each active emission
    ! species (as for GA9).
    ! Register 2D emissions first, followed by 3D emissions; emissions
    ! must be registered in order of dimensionality for compatibility with
    ! the UKCA time step call which expects an array of 2D emission fields
    ! and an array of 3D emission fields. These arrays are expected to
    ! contain fields corresponding to consecutive emission id numbers,
    ! starting at 1, with the 3D fields having the highest numbers.
    ! The names of the emission species corresponding to each field array
    ! are given in separate reference arrays created below.

    ! 'Dictionary' of emission species and long_names
    integer, parameter :: num_tot = 23 !  Max for Strattrop+GLOMAP
    integer :: num_2d, num_3d

    character(len=ukca_maxlen_emiss_long_name) :: long_name

    integer(i_um) :: emiss_id
    integer :: n_emissions
    integer :: n2d, n3d
    integer :: i, j, k
    integer :: anc_per_emiss(num_tot)

    logical :: l_three_dim, l_found
    integer :: low_lev, hi_lev
    character(len=ukca_maxlen_emiss_vert_fact) :: vert_fact
    character(len=ukca_maxlen_emiss_var_name), allocatable :: tmp_names(:)

    character(len=ukca_maxlen_emiss_var_name) :: field_varname
    character(len=*), parameter  :: emiss_units = 'kg m-2 s-1'

    ! Populate dictionary of names. Using DATA as length of strings varies.
    ! Element 1 is the emission name
    ! Element 2 is the long name
    character(len=ukca_maxlen_emiss_long_name) :: emname_map(num_tot,2)
    data emname_map(1,:) /'C2H6', 'C2H6 surf emissions'/
    data emname_map(2,:) /'C3H8', 'C3H8 surf emissions'/
    data emname_map(3,:) /'C5H8', 'Biogenic C5H8 surf emissions'/
    data emname_map(4,:) /'CH4', 'CH4 surf emissions'/
    data emname_map(5,:) /'CO', 'CO surf emissions'/
    data emname_map(6,:) /'HCHO', 'HCHO surf emissions'/
    data emname_map(7,:) /'Me2CO', 'Me2CO surf emissions'/
    data emname_map(8,:) /'MeCHO', 'MeCHO surf emissions'/
    data emname_map(9,:) /'MeOH',        &
                           'Biogenic surface methanol (CH3OH) emissions'/
    data emname_map(10,:) /'NH3', 'NH3 surf emissions'/
    data emname_map(11,:) /'NO', 'NOx surf emissions'/
    data emname_map(12,:) /'BC_biofuel', 'BC biofuel surf emissions'/
    data emname_map(13,:) /'BC_fossil', 'BC fossil fuel surf emissions'/
    data emname_map(14,:) /'DMS', 'DMS emissions expressed as sulfur'/
    data emname_map(15,:) /'Monoterp',    &
                           'Monoterpene surf emissions expressed as carbon'/
    data emname_map(16,:) /'OM_biofuel',  &
                           'OC biofuel surf emissions expressed as carbon'/
    data emname_map(17,:) /'OM_fossil',   &
                          'OC fossil fuel surf emissions expressed as carbon'/
    data emname_map(18,:) /'SO2_high',    &
                           'SO2 high level emissions expressed as sulfur'/
    data emname_map(19,:) /'SO2_low',     &
                           'SO2 low level emissions'/
    ! could be 2D or could be 3D
    data emname_map(20,:) /'BC_biomass', 'BC biomass emissions'/
    data emname_map(21,:) /'OM_biomass', &
                           'OC biomass emissions expressed as carbon'/
    ! 3-D emissions
    data emname_map(22,:) /'NO_aircrft', 'NOx aircraft emissions'/
    data emname_map(23,:) /'SO2_nat',    &
                           'SO2 natural emissions'/

    ! End of Header

    n_emissions = size(emiss_names)
    anc_per_emiss(:) = 1
    if (emissions == emissions_GC3) then
      num_2d=19
      num_3d=4
      ! SO2_low and SO2_nat expressed as S in GC3 emissions
      emname_map(19,2) = trim(emname_map(19,2)) // ' expressed as sulfur'
      emname_map(23,2) = trim(emname_map(23,2)) // ' expressed as sulfur'
    else if (emissions == emissions_GC5) then
      num_2d=21
      num_3d=2
      ! BC_biomass and OM_biomass have _low and _high variants in GC5 set
      anc_per_emiss(20) = 2
      anc_per_emiss(21) = 2
    end if

    if ( n_emissions > (num_2d+num_3d) ) then
      write(log_scratch_space,'(2(A,I0),A)')' Number of expected emissions: ', &
        n_emissions,' higher than that currently registered: ',(num_2d+num_3d),&
        ' Update registered emissions in um_ukca_init_mod/ukca_emiss_init().'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Check that all emission species returned by UKCA are recognised here
    do i = 1, n_emissions
      l_found = .false.
      do j = 1, SIZE(emname_map, DIM=1)
        if ( emiss_names(i) == emname_map(j,1) ) l_found = .true.
      end do
      if ( .not. l_found ) then
        write( log_scratch_space, '(A,A,A)' ) ' Unknown Emission species ',    &
           emiss_names(i)//' returned by get_emission_varlist. ',              &
          ' Check registered emissions in um_ukca_init_mod/ukca_emiss_init().'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
    end do

    ! Determine number of offline emission entries required.
    ! Allow for multiple emission fields for each active emission species
    ! (consistent with GA9 requirements).
    n_emiss_slots = 0
    do i = 1, n_emissions
      do j = 1, num_2d+num_3d
        if ( emiss_names(i) == emname_map(j,1) ) then
          n_emiss_slots = n_emiss_slots + n_slots(i) * anc_per_emiss(j)
        end if
      end do
    end do

    ! === Register 2-D emissions ===

    allocate(tmp_names(sum(anc_per_emiss(1:num_2d))))  ! Max possible
    ! counter
    n2d = 0
    do i = 1, n_emissions
      do j = 1, num_2d
        if ( emiss_names(i) == emname_map(j,1) ) then
          do k = 1, anc_per_emiss(j)
            ! species found in dictionary, set up parameters
            long_name = emname_map(j,2)
            l_three_dim = .false.  ! A few default values.
            vert_fact = 'surface'
            low_lev = 1            ! Can pass as arguments in all calls since
            hi_lev = 1             ! these are ignored for vert_fact= surface
            n2d = n2d + 1
            field_varname = 'emissions_' // emiss_names(i)

            ! Assign long name as in GA9 ancil file.
            ! *** WARNING ***
            ! This is currently hard-wired in the absence of functionality
            ! to handle the NetCDF variable attributes that provide the long names
            ! in the UM. There is therefore no guarantee of compatibility with
            ! arbitrary ancil files and care must be taken to ensure that the
            ! data in the ancil files provided are consistent with the names
            ! defined here.
            if (emiss_names(i) == 'SO2_high') then
              vert_fact = 'high_level'
              low_lev = 8     ! Corresponding to L70 and L85
              hi_lev = 8
            end if
            if ( anc_per_emiss(j) > 1 ) then
              if (k == 1) then
                vert_fact = 'high_level'
                low_lev = 1
                hi_lev = 20
                field_varname = trim(field_varname) // '_high'
                long_name = trim(long_name) // ' high level'
              else
                field_varname = trim(field_varname) // '_low'
                long_name = trim(long_name) // ' low level'
              end if
            end if
            tmp_names(n2d) = field_varname
            call ukca_register_emission( n_emiss_slots, field_varname,         &
                                         emiss_names(i), emiss_units,          &
                                         l_three_dim, emiss_id,                &
                                         long_name=long_name,                  &
                                         vert_fact=vert_fact,                  &
                                         lowest_lev=low_lev, highest_lev=hi_lev)
            if (emiss_id /= int( n2d, i_um )) then
              write( log_scratch_space, '(A,I0,A,I0)' )                        &
                'Unexpected id (', emiss_id,                                   &
                ') assigned on registering a 2D UKCA emission. Expected ', n2d
              call log_event( log_scratch_space, LOG_LEVEL_ERROR )
            end if
          end do ! Loop over ancillaries per emission
        end if  ! species found in emnames ?

      end do    ! Loop over emnames_map 2-D
    end do      ! Loop over n_emissions

    ! Create reference array of 2D emission names
    allocate(emiss_names_flat(n2d))
    do i = 1, n2d
      emiss_names_flat(i) = tmp_names(i)
    end do
    deallocate(tmp_names)

    ! Register 3-D emissions
    l_three_dim = .true.
    vert_fact = 'all_levels'
    allocate(tmp_names(num_3d))  ! Max possible
    n3d = 0
    do i = 1, n_emissions
      do j = num_2d+1, SIZE(emname_map, DIM=1)
        if ( emiss_names(i) == emname_map(j,1) ) then
          long_name = emname_map(j,2)

          ! species found in dictionary, set up parameters
          field_varname = 'emissions_' // emiss_names(i)
          n3d = n3d + 1
          tmp_names(n3d) = field_varname
          ! Register with parameters as  in GA9 ancil file.
          call ukca_register_emission( n_emiss_slots, field_varname,           &
                                       emiss_names(i), emiss_units,            &
                                       l_three_dim, emiss_id,                  &
                                       long_name=long_name,                    &
                                       vert_fact=vert_fact )

          if (emiss_id /= int( n3d+n2d, i_um )) then
            write( log_scratch_space, '(A,I0,A,I0)' )                          &
              'Unexpected id (', emiss_id,                                     &
              ') assigned on registering a 3D UKCA emission. Expected ', n3d+n2d
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
          end if
        end if     ! species found in 3-D list

      end do       ! loop over emnames_map/3D
    end do         ! Loop over n_emissions

    ! Create reference array of 3D emission names
    allocate(emiss_names_fullht(n3d))
    do i = 1, n3d
      emiss_names_fullht(i) = tmp_names(i)
    end do
    deallocate(tmp_names)

    ! Print out number of emissions species for this run
    write( log_scratch_space, '(A,2(I0,A))' ) 'UKCA emissions: ',n2d,          &
      ' 2-d and ', n3d,' 3-d species active in this run '
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine ukca_emiss_init


  !>@brief Set the lists of fields required for the UKCA configuration
  subroutine set_ukca_field_lists()

    implicit none

    ! Local variables

    integer :: i
    integer :: n
    integer :: n_tot

    ! Variables for UKCA error handling
    integer :: ukca_errcode
    character(len=ukca_maxlen_message) :: ukca_errmsg
    character(len=ukca_maxlen_procname) :: ukca_errproc

    ! Get list of tracers required by the current UKCA configuration
    call ukca_get_tracer_varlist( tracer_names, ukca_errcode,                  &
                                  error_message=ukca_errmsg,                   &
                                  error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'Tracers required (', size(tracer_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(tracer_names)
      write( log_scratch_space, '(A)' ) tracer_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    ! Get list of NTPs required by the current UKCA configuration
    call ukca_get_ntp_varlist( ntp_names, ukca_errcode,                        &
                               error_message=ukca_errmsg,                      &
                               error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'NTPs required (', size(ntp_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(ntp_names)
      write( log_scratch_space, '(A)' ) ntp_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    ! Get lists of environmental drivers required by the current UKCA
    ! configuration.
    ! Note that these group lists derived from UKCA's master list must
    ! be requested here to force UKCA to set up the group lists prior to
    ! their use in kernel calls. (Setup within kernels is not thread-safe.)

    call ukca_get_envgroup_varlists(                                           &
           ukca_errcode,                                                       &
           varnames_scalar_real_ptr=env_names_scalar_real,                     &
           varnames_flat_integer_ptr=env_names_flat_integer,                   &
           varnames_flat_real_ptr=env_names_flat_real,                         &
           varnames_flat_logical_ptr=env_names_flat_logical,                   &
           varnames_flatpft_real_ptr=env_names_flatpft_real,                   &
           varnames_fullht_real_ptr=env_names_fullht_real,                     &
           varnames_fullht0_real_ptr=env_names_fullht0_real,                   &
           varnames_fullhtp1_real_ptr=env_names_fullhtp1_real,                 &
           varnames_bllev_real_ptr=env_names_bllev_real,                       &
           varnames_entlev_real_ptr=env_names_entlev_real,                     &
           varnames_land_real_ptr=env_names_land_real,                         &
           varnames_landtile_real_ptr=env_names_landtile_real,                 &
           varnames_landtile_logical_ptr=env_names_landtile_logical,           &
           varnames_landpft_real_ptr=env_names_landpft_real,                   &
           varnames_fullhtphot_real_ptr=env_names_fullhtphot_real,             &
           error_message=ukca_errmsg, error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    n_tot = 0

    n = size(env_names_scalar_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in scalar real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_scalar_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_integer)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat integer group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_flat_integer(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_flat_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_logical)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat logical group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_flat_logical(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flatpft_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat PFT real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_flatpft_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullht_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_fullht_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullht0_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height + level 0 real group required (',  &
      n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_fullht0_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullhtp1_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height + 1 real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_fullhtp1_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_bllev_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in boundary layer real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_bllev_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_entlev_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in entrainment levels real group required (',     &
      n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_entlev_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_land_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_land_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landtile_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land tile real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_landtile_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landtile_logical)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land tile logical group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_landtile_logical(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landpft_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land PFT real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_landpft_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullhtphot_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers as photolysis rates required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_fullhtphot_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    write( log_scratch_space, '(A,I0)' )                                      &
      'Total number of environmental drivers required: ', n_tot
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Get list of active offline emissions species in the current UKCA
    ! configuration
    call ukca_get_emission_varlist( emiss_names, n_slots, ukca_errcode,        &
                                    error_message=ukca_errmsg,                 &
                                    error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'Offline emissions active (', size(emiss_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(emiss_names)
      write( log_scratch_space, '(A)' ) emiss_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

  end subroutine set_ukca_field_lists

  subroutine aerosol_ukca_dust_only_init( row_length, rows, model_levels,      &
                                          bl_levels, timestep, l_param_conv )

    implicit none

    integer, intent(in) :: row_length      ! X dimension of model grid
    integer, intent(in) :: rows            ! Y dimension of model grid
    integer, intent(in) :: model_levels    ! Z dimension of model grid
    integer, intent(in) :: bl_levels       ! No. of Z levels in boundary layer
    real(r_um), intent(in) :: timestep     ! Model time step (s)
    logical, intent(in) :: l_param_conv    ! True if convection is parameterized

    ! Local variables

    integer(i_um) :: emiss_id
    integer(i_um) :: n_emissions
    integer(i_um) :: n
    integer(i_um) :: n_previous
    integer(i_um) :: i, j

    logical :: l_three_dim

    character(len=ukca_maxlen_emiss_long_name) :: long_name
    character(len=ukca_maxlen_emiss_var_name), allocatable :: tmp_names(:)

    character(len=ukca_maxlen_emiss_var_name) :: field_varname
    character(len=*), parameter  :: emiss_units = 'kg m-2 s-1'

    ! Variables for UKCA error handling
    integer :: ukca_errcode
    character(len=ukca_maxlen_message) :: ukca_errmsg
    character(len=ukca_maxlen_procname) :: ukca_errproc

    ! Set up proto-GA configuration based on GA9.
    ! The ASAD Newton-Raphson Offline Oxidants scheme (ukca_chem_offline)
    ! is substituted for the Explicit backward-Euler scheme used in GA9
    ! which cannot be called by columns. Hence, configuration values for
    ! nrsteps and l_ukca_asad_columns are included.
    ! Other configuration values, with the exception of temporary logicals
    ! specifiying fixes, are set to match GA9 science settings or taken
    ! from the LFRic context. Unlike GA9, all fixes that are controlled by
    ! temporary logicals will be on. (i.e. the defaults for these
    ! logicals, .true. by convention in UKCA, are not overridden.)

    call ukca_setup( ukca_errcode,                                             &

           ! Switch to skip setting up constants: these will have already
           ! been set when GLOMAP-clim is set up
           !
           l_skip_const_setup=.true.,                                          &

           ! Context information
           !
           row_length=row_length,                                              &
           rows=rows,                                                          &
           model_levels=model_levels,                                          &
           bl_levels=bl_levels,                                                &
           nlev_ent_tr_mix=nlev_ent_tr_mix,                                    &
           timestep=timestep,                                                  &

           ! General UKCA configuration options
           !
           i_ukca_chem=ukca_chem_off,                                          &
           l_ukca_mode=.true.,                                                 &
           l_fix_tropopause_level=.true.,                                      &
           l_ukca_persist_off=.true.,                                          &

           ! Chemistry configuration options
           !
           chem_timestep=3600,                                                 &
           i_chem_timestep_halvings=0,                                         &

           ! UKCA environmental driver configuration options
           !
           l_param_conv=l_param_conv,                                          &
           l_environ_rel_humid=.true.,                                         &

           ! General GLOMAP configuration options
           !
           i_mode_nzts=15,                                                     &
           i_mode_setup=6,                                                     &
           i_mode_nucscav=i_mode_nucscav,                                      &
           l_cv_rainout=.not.(l_ukca_plume_scav),                              &
           l_dust_mp_slinn_impc_scav=.true.,                                   &
           l_dust_mp_ageing=.false.,                                           &

           ! GLOMAP emissions configuration options
           !
           l_ukca_primdu=.true.,                                               &

           ! GLOMAP feedback configuration options
           !
           l_ukca_radaer=.false.,                                              &

           ! Callback procedures
           !
           proc_bl_tracer_mix = bl_tracer_mix,                                 &

           ! Return status information
           !
           error_message=ukca_errmsg,                                          &
           error_routine=ukca_errproc )

    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A)' )                                   &
           'UKCA error ', ukca_errcode, ' in ' // trim(ukca_errproc) // ' : '  &
           // trim(ukca_errmsg)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Retrieve the lists of required fields for the configuration
    call set_ukca_field_lists()

    ! Set up indexing data needed for plume scavenging of UKCA tracers in the
    ! GR convection scheme
    if (l_ukca_plume_scav) then
      n = size(tracer_names)
      allocate(nm_spec_active(n))
      do i = 1, n
        nm_spec_active(i) = tracer_names(i)(1:nmspec_len)
      end do

      call ukca_set_conv_indices()
      tracer_info%i_ukca_first = 1
      tracer_info%i_ukca_last  = n

    end if

    ! Switch on optional UM microphysics diagnostics required by UKCA
    if (any(env_names_fullht_real(:) == fldname_autoconv))                     &
      l_praut_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_accretion))                    &
      l_pracw_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_rim_cry))                      &
      l_piacw_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_rim_agg))                      &
      l_psacw_diag = .true.

    ! Switch on optional UM boundary layer diagnostics required by UKCA
    if (any(env_names_bllev_real(:) == fldname_bl_tke))                        &
      bl_diag%l_request_tke = .true.

    ! Emissions registration:
    ! One emission field is provided for each active emission species (as for
    ! GA9).
    ! Register 2D emissions first, followed by 3D emissions; emissions
    ! must be registered in order of dimensionality for compatibility with
    ! the UKCA time step call which expects an array of 2D emission fields
    ! and an array of 3D emission fields. These arrays are expected to
    ! contain fields corresponding to consecutive emission id numbers,
    ! starting at 1, with the 3D fields having the highest numbers.
    ! The names of the emission species corresponding to each field array
    ! are given in separate reference arrays created below.

    n_emissions = size(emiss_names)
    allocate(tmp_names(n_emissions))

    ! Determine number of offline emission entries required.
    ! Allow for one emission field for each active emission species
    ! (consistent with GA9 requirements).
    n_emiss_slots = 0
    do i = 1, size(emiss_names)
      do j = 1, size(required_emissions)
        if ( required_emissions(j) == emiss_names(i) ) then
          n_emiss_slots = n_emiss_slots + n_slots(i)
        end if
      end do
    end do

    ! Register 2D emissions
    n = 0
    do i = 1, n_emissions
      if ( any( required_2d_emissions(:) == emiss_names(i) ) ) then

        field_varname = 'emissions_' // emiss_names(i)
        l_three_dim = .false.

        ! Assign long name as in GA9 ancil file.
        ! *** WARNING ***
        ! This is currently hard-wired in the absence of functionality
        ! to handle the NetCDF variable attributes that provide the long names
        ! in the UM. There is therefore no guarantee of compatibility with
        ! arbitrary ancil files and care must be taken to ensure that the
        ! data in the ancil files provided are consistent with the names
        ! defined here.

        select case( emiss_names(i) )
        case('DMS')
          long_name = 'DMS emissions expressed as sulfur'

        case('Monoterp')
          long_name = 'Monoterpene surf emissions expressed as carbon'

        case('SO2_low')
          long_name = 'SO2 low level emissions expressed as sulfur'

        case('SO2_high')
          long_name = 'SO2 high level emissions expressed as sulfur'

        end select

        if (emiss_names(i) == 'SO2_high') then
          call ukca_register_emission( n_emiss_slots, field_varname, &
                                       emiss_names(i), emiss_units,  &
                                       l_three_dim, emiss_id,        &
                                       long_name=long_name,          &
                                       vert_fact='high_level',       &
                                       lowest_lev=8, highest_lev=8 )
        else
          call ukca_register_emission( n_emiss_slots, field_varname, &
                                       emiss_names(i), emiss_units,  &
                                       l_three_dim, emiss_id,        &
                                       long_name=long_name,          &
                                       vert_fact='surface' )
        end if

        n = n + 1

        if (emiss_id /= int( n, i_um )) then
          write( log_scratch_space, '(2(A,I0))' ) &
            'Unexpected id (', emiss_id,          &
            ') assigned on registering a 2D UKCA emission. Expected ', n
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if

        tmp_names(n) = field_varname

      end if
    end do

    ! Create reference array of 2D emission names
    allocate(emiss_names_flat(n))
    do i = 1, n
      emiss_names_flat(i) = tmp_names(i)
    end do

    ! Register 3D emissions
    n_previous = n
    n = 0

    do i = 1, n_emissions
      if ( emiss_names(i) == 'SO2_nat' ) then

        field_varname = 'emissions_' // emiss_names(i)
        l_three_dim = .true.

        ! Assign long name as in GA9 ancil file.
        ! *** WARNING ***
        ! This is currently hard-wired in the absence of functionality
        ! to handle the NetCDF variable attributes. There is therefore
        ! no guarantee of compatibility with arbitrary ancil files.

        long_name = 'SO2 natural emissions expressed as sulfur'

        call ukca_register_emission( n_emiss_slots, field_varname,             &
                                     emiss_names(i), emiss_units, l_three_dim, &
                                     emiss_id, long_name=long_name,            &
                                     vert_fact='all_levels' )

        n = n + 1

        if (emiss_id /= int( n + n_previous, i_um )) then
          write( log_scratch_space, '(2(A,I0))' ) &
            'Unexpected id (', emiss_id,          &
            ') assigned on registering 3D UKCA emission. Expected ', n
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if

        tmp_names(n) = field_varname

      end if
    end do

    ! Create reference array of 3D emission names
    allocate(emiss_names_fullht(n))
    do i = 1, n
      emiss_names_fullht(i) = tmp_names(i)
    end do

  end subroutine aerosol_ukca_dust_only_init

! ----------------------------------------------------------------------
subroutine allocate_fastjx_filevars()
! ----------------------------------------------------------------------
! Description:
!
! allocate arrays that will hold data from FastJX spectral files. 
! ----------------------------------------------------------------------

implicit none

! Ensure arrays are not already allocated
call deallocate_fastjx_filevars()

allocate(fl(photol_max_wvl))
allocate(q1d(photol_max_wvl, photol_num_tvals))
allocate(qo2(photol_max_wvl, photol_num_tvals))
allocate(qo3(photol_max_wvl, photol_num_tvals))
allocate(qqq(photol_max_wvl, 2, photol_max_crossec))
allocate(qrayl(photol_max_wvl+1))
allocate(tqq(photol_num_tvals, photol_max_crossec))
allocate(wl(photol_max_wvl))

allocate(daa(photol_max_miesets))
allocate(paa(photol_sw_phases, photol_sw_band_aer, photol_max_miesets))
allocate(qaa(photol_sw_band_aer, photol_max_miesets))
allocate(raa(photol_max_miesets))
allocate(saa(photol_sw_band_aer, photol_max_miesets))
allocate(waa(photol_sw_band_aer, photol_max_miesets))

allocate(solcyc_av(photol_n_solcyc_av))
allocate(solcyc_quanta(photol_wvl_intervals))
allocate(solcyc_spec(photol_max_wvl))
allocate(solcyc_ts(n_solcyc_ts))

allocate(titlej(photol_max_crossec))

return
end subroutine allocate_fastjx_filevars

subroutine deallocate_fastjx_filevars()
  ! ----------------------------------------------------------------------
  ! Description:
  !
  ! Deallocate arrays that hold data from FJX spectral files.
  ! ----------------------------------------------------------------------

implicit none

if (allocated(titlej)) deallocate(titlej)
if (allocated(solcyc_ts)) deallocate(solcyc_ts)
if (allocated(solcyc_spec)) deallocate(solcyc_spec)
if (allocated(solcyc_quanta)) deallocate(solcyc_quanta)
if (allocated(solcyc_av)) deallocate(solcyc_av)

if (allocated(waa)) deallocate(waa)
if (allocated(saa)) deallocate(saa)
if (allocated(raa)) deallocate(raa)
if (allocated(qaa)) deallocate(qaa)
if (allocated(paa)) deallocate(paa)
if (allocated(daa)) deallocate(daa)

if (allocated(jlabel)) deallocate(jlabel)
if (allocated(jind)) deallocate(jind)
if (allocated(jfacta)) deallocate(jfacta)

if (allocated(wl)) deallocate(wl)
if (allocated(tqq)) deallocate(tqq)
if (allocated(qrayl)) deallocate(qrayl)
if (allocated(qqq)) deallocate(qqq)
if (allocated(qo3)) deallocate(qo3)
if (allocated(qo2)) deallocate(qo2)
if (allocated(q1d)) deallocate(q1d)
if (allocated(fl)) deallocate(fl)

return
end subroutine deallocate_fastjx_filevars


end module um_ukca_init_mod
