!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of variables for UM physics schemes, which
!>         are either fixed in LFRic or derived from LFRic inputs

module um_physics_init_mod

  ! LFRic namelists which have been read
  use aerosol_config_mod,        only : activation_scheme,                     &
                                        activation_scheme_jones,               &
                                        glomap_mode,                           &
                                        glomap_mode_climatology,               &
                                        glomap_mode_ukca,                      &
                                        glomap_mode_dust_and_clim,             &
                                        glomap_mode_radaer_test,               &
                                        glomap_mode_off,                       &
                                        aclw_file,                             &
                                        acsw_file,                             &
                                        anlw_file,                             &
                                        answ_file,                             &
                                        crlw_file,                             &
                                        crsw_file,                             &
                                        prec_file,                             &
                                        l_radaer,                              &
                                        easyaerosol_sw,                        &
                                        easyaerosol_lw,                        &
                                        murk, murk_prognostic

  use blayer_config_mod,         only : a_ent_shr, a_ent_2_in => a_ent_2,      &
                                        cbl_opt, cbl_opt_conventional,         &
                                        cbl_opt_standard, cbl_opt_adjustable,  &
                                        cbl_mix_fac,                           &
                                        dec_thres_cloud_in => dec_thres_cloud, &
                                        dyn_diag, dyn_diag_zi_l_sea,           &
                                        dyn_diag_ri_based, dyn_diag_zi_l_cu,   &
                                        near_neut_z_on_l_in=>near_neut_z_on_l, &
                                        flux_bc_opt_in => flux_bc_opt,         &
                                        flux_bc_opt_interactive,               &
                                        flux_bc_opt_specified_scalars,         &
                                        flux_bc_opt_specified_scalars_tstar,   &
                                        fric_heating_in => fric_heating,       &
                                        free_atm_mix, free_atm_mix_to_sharp,   &
                                        free_atm_mix_ntml_corrected,           &
                                        free_atm_mix_free_trop_layer,          &
                                        interp_local, interp_local_gradients,  &
                                        interp_local_cf_dbdz,                  &
                                        new_kcloudtop, p_unstable,             &
                                        reduce_fa_mix, noice_in_turb,          &
                                        reduce_fa_mix_inv_and_cu_lcl,          &
                                        reduce_fa_mix_inv_only,                &
                                        relax_sc_over_cu_in=>relax_sc_over_cu, &
                                        sbl_opt, sbl_opt_sharpest,             &
                                        sbl_opt_sharp_sea_mes_land,            &
                                        sg_orog_mixing_in => sg_orog_mixing,   &
                                        sg_orog_mixing_none,                   &
                                        sg_orog_mixing_shear_plus_lambda,      &
                                        zhloc_depth_fac_in => zhloc_depth_fac, &
                                        bl_levels_in => bl_levels

  use cloud_config_mod,          only : scheme, scheme_smith, scheme_pc2,     &
                                        scheme_bimodal,                       &
                                        rh_crit, rh_crit_opt,                 &
                                        rh_crit_opt_namelist, rh_crit_opt_tke,&
                                        pc2ini, pc2ini_smith,                 &
                                        pc2ini_bimodal,                       &
                                        cff_spread_rate_in => cff_spread_rate,&
                                        falliceshear_method_in =>             &
                                        falliceshear_method,                  &
                                        falliceshear_method_real,             &
                                        falliceshear_method_constant,         &
                                        subgrid_qv, ice_width_in => ice_width,&
                                        ez_subcrit, ez_max,  &
                                        two_d_fsd_factor_in => two_d_fsd_factor

  use convection_config_mod,     only : cv_scheme,                    &
                                        cv_scheme_gregory_rowntree,   &
                                        cv_scheme_lambert_lewis,      &
                                        number_of_convection_substeps,&
                                        cape_timescale_in => cape_timescale

  use extrusion_config_mod,      only : domain_top, number_of_layers

  use formulation_config_mod,    only : moisture_formulation,    &
                                        moisture_formulation_dry

  use microphysics_config_mod,   only : a_ratio_exp_in => a_ratio_exp,       &
                                        a_ratio_fac_in => a_ratio_fac,       &
                                        graupel_scheme, graupel_scheme_none, &
                                        graupel_scheme_modified,             &
                                        droplet_tpr, shape_rime,             &
                                        qcl_rime,                            &
                                        ndrop_surf_in => ndrop_surf,         &
                                        z_surf_in => z_surf,                 &
                                        turb_gen_mixph,                      &
                                        orog_rain, orog_rime,                &
                                        prog_tnuc, orog_block,               &
                                        fcrit_in => fcrit,                   &
                                        nsigmasf_in => nsigmasf,             &
                                        nscalesf_in => nscalesf,             &
                                        microphysics_casim

  use mixing_config_mod,         only : smagorinsky,                 &
                                        mixing_method => method,     &
                                        method_3d_smag,              &
                                        method_2d_smag,              &
                                        method_blend_smag_fa,        &
                                        method_blend_1dbl_fa,        &
                                        mix_factor_in => mix_factor, &
                                        leonard_term

  use radiation_config_mod,      only : topography, topography_horizon

  use section_choice_config_mod, only : aerosol,           &
                                        aerosol_um,        &
                                        boundary_layer,    &
                                        boundary_layer_um, &
                                        convection,        &
                                        convection_um,     &
                                        cloud,             &
                                        cloud_um,          &
                                        electric,          &
                                        electric_um,       &
                                        microphysics,      &
                                        microphysics_um,   &
                                        orographic_drag,   &
                                        orographic_drag_um,&
                                        radiation,         &
                                        radiation_socrates,&
                                        spectral_gwd,      &
                                        spectral_gwd_um,   &
                                        stochastic_physics,&
                                        stochastic_physics_um, &
                                        surface,           &
                                        surface_jules

  use spectral_gwd_config_mod,   only :                                       &
                                 ussp_launch_factor_in => ussp_launch_factor, &
                                 wavelstar_in => wavelstar,                   &
                                 add_cgw_in => add_cgw,                       &
                                 cgw_scale_factor_in => cgw_scale_factor

  use socrates_init_mod, only: n_sw_band,                                      &
                               n_lw_band

  use stochastic_physics_config_mod, only: use_random_parameters,              &
                                           rp_bl_a_ent_1,                      &
                                           rp_bl_a_ent_shr_max,                &
                                           rp_bl_a_ent_shr,                    &
                                           rp_bl_cbl_mix_fac,                  &
                                           rp_bl_cld_top_diffusion,            &
                                           rp_bl_stable_ri_coef,               &
                                           rp_bl_min_mix_length,               &
                                           rp_bl_neutral_mix_length,           &
                                           rp_bl_ricrit,                       &
                                           rp_bl_smag_coef,                    &
                                           rp_callfreq,                        &
                                           rp_cycle_in,                        &
                                           rp_cycle_out,                       &
                                           rp_cycle_tm,                        &
                                           rp_decorr_ts,                       &
                                           rp_lsfc_lai_mult_max,               &
                                           rp_lsfc_lai_mult_min,               &
                                           rp_lsfc_lai_mult,                   &
                                           rp_lsfc_z0hm_pft_max,               &
                                           rp_lsfc_z0hm_pft_min,               &
                                           rp_lsfc_z0hm_pft,                   &
                                           rp_mp_ndrop_surf,                   &
                                           rp_ran_max

  use orographic_drag_config_mod, only:  include_moisture,          &
                                         include_moisture_lowmoist, &
                                         include_moisture_moist,    &
                                         include_moisture_dry


  ! Other LFRic modules used
  use constants_mod,        only : i_def, l_def, r_um, i_um, rmdi, r_def, r_bl
  use log_mod,              only : log_event,         &
                                   log_scratch_space, &
                                   LOG_LEVEL_ERROR,   &
                                   LOG_LEVEL_INFO

  use mr_indices_mod,       only : nummr_to_transport

  ! UM modules used
  use cderived_mod,         only : delta_lambda, delta_phi
  use nlsizes_namelist_mod, only : bl_levels, model_levels

  implicit none

  integer(i_def), protected :: n_radaer_mode
  integer(i_def), protected :: n_aer_mode_sw
  integer(i_def), protected :: n_aer_mode_lw
  integer(i_def), protected :: mode_dimen
  integer(i_def), protected :: sw_band_mode
  integer(i_def), protected :: lw_band_mode

  private
  public :: um_physics_init, &
            n_radaer_mode, mode_dimen, sw_band_mode, lw_band_mode, &
            n_aer_mode_sw, n_aer_mode_lw

contains

  !>@brief Initialise UM physics variables which are either fixed in LFRic
  !>        or derived from LFRic inputs or JULES variables
  !>@details This file sets many parameters and switches which are currently
  !>          in the UM namelists. Many of these will never be promoted to the
  !>          LFRic namelist as they are legacy options not fit for future use.
  !>          Hence we set them here until such time as we can retire them
  !>          from the UM code.
  !>        Other parameters and switches which are genuinely input variables,
  !>         via the LFRic namelists, are also set here for the UM code.

  subroutine um_physics_init()

    ! UM modules containing things that need setting and setup routines
    use bl_option_mod, only: i_bl_vn, sbl_op, ritrans,                     &
         cbl_op, lambda_min_nml, local_fa, keep_ri_fa,                     &
         sg_orog_mixing, fric_heating, idyndiag,                           &
         zhloc_depth_fac, flux_grad, entr_smooth_dec,                      &
         relax_sc_over_cu, bl_res_inv, blending_option,                    &
         a_ent_shr_nml, alpha_cd, puns, pstb, kprof_cu,                    &
         non_local_bl, flux_bc_opt, i_bl_vn_9c, sharp_sea_mes_land,        &
         lem_conven, to_sharp_across_1km, off, on, DynDiag_Ribased,        &
         DynDiag_ZL_corrn, blend_allpoints, ng_stress, lem_std, lem_adjust,&
         interactive_fluxes, specified_fluxes_only, except_disc_inv,       &
         ntml_level_corrn, free_trop_layers, sharpest, sg_shear_enh_lambda,&
         l_new_kcloudtop, buoy_integ, l_reset_dec_thres, DynDiag_ZL_CuOnly,&
         var_diags_opt, i_interp_local, i_interp_local_gradients,          &
         split_tke_and_inv, l_noice_in_turb, l_use_var_fixes,              &
         i_interp_local_cf_dbdz, tke_diag_fac, a_ent_2, dec_thres_cloud,   &
         dec_thres_cu, near_neut_z_on_l, blend_gridindep_fa,               &
         specified_fluxes_tstar
    use cloud_inputs_mod, only: i_cld_vn, forced_cu, i_rhcpt, i_cld_area,  &
         rhcrit, ice_fraction_method,falliceshear_method, cff_spread_rate, &
         l_subgrid_qv, ice_width, min_liq_overlap, i_eacf, not_mixph,      &
         i_pc2_checks_cld_frac_method, l_ensure_min_in_cloud_qcf,          &
         i_pc2_init_logic, dbsdtbs_turb_0,                                 &
         i_pc2_erosion_method, i_pc2_homog_g_method, i_pc2_init_method,    &
         check_run_cloud, l_pc2_implicit_erosion,                          &
         forced_cu_fac, i_pc2_conv_coupling, allicetdegc, starticetkelvin, &
         l_bm_ez_subcrit_only, ez_max_bm
    use cloud_config_mod, only: cld_fsd_hill
    use cv_run_mod, only: icvdiag, cvdiag_inv, cvdiag_sh_wtest,            &
         limit_pert_opt, tv1_sd_opt, iconv_congestus, iconv_deep,          &
         ent_fac_dp, cldbase_opt_dp, cldbase_opt_sh, w_cape_limit,         &
         l_param_conv, i_convection_vn, l_ccrad, l_mom, adapt, amdet_fac,  &
         bl_cnv_mix, anvil_factor, ccw_for_precip_opt, cldbase_opt_md,     &
         cnv_wat_load_opt,dd_opt, deep_cmt_opt, eff_dcff, eff_dcfl,        &
         ent_dp_power, ent_fac_md, ent_opt_dp, ent_opt_md, fac_qsat,       &
         ent_fac_md, iconv_deep, iconv_mid, iconv_shallow,l_3d_cca,        &
         l_anvil, l_cmt_heating, l_cv_conserve_check, l_safe_conv,         &
         mdet_opt_dp, mdet_opt_md, mid_cnv_pmin, mparwtr, qlmin,           &
         n_conv_calls, efrac,                                              &
         qstice, r_det, sh_pert_opt,t_melt_snow, termconv, tice,           &
         thpixs_mid,                                                       &
         tower_factor,ud_factor, fdet_opt, anv_opt, cape_timescale,        &
         cca2d_dp_opt,cca2d_md_opt,cca2d_sh_opt,                           &
         cca_dp_knob,cca_md_knob,cca_sh_knob,                              &
         ccw_dp_knob,ccw_md_knob,ccw_sh_knob,                              &
         cnv_cold_pools,dil_plume_water_load,l_cloud_deep, mid_cmt_opt,    &
         plume_water_load, rad_cloud_decay_opt, cape_bottom, cape_top,     &
         cape_min, i_convection_vn_6a, i_cv_llcs, midtrig_opt,             &
         llcs_cloud_precip, llcs_opt_all_rain, llcs_rhcrit, llcs_timescale,&
         check_run_convection, l_fcape, cape_ts_min, cape_ts_max,          &
         cpress_term, pr_melt_frz_opt, llcs_opt_crit_condens,              &
         llcs_detrain_coef, l_prog_pert, md_pert_opt, l_jules_flux,        &
         l_reset_neg_delthvu,                                              &
         l_conv_prog_precip, l_conv_prog_dtheta, l_conv_prog_dq,           &
         adv_conv_prog_dtheta, adv_conv_prog_dq,                           &
         tau_conv_prog_precip, tau_conv_prog_dtheta, tau_conv_prog_dq,     &
         prog_ent_grad, prog_ent_int, prog_ent_max, prog_ent_min,          &
         ent_fac_sh, c_mass_sh, orig_mdet_fac
    use cv_param_mod, only: mtrig_ntml, md_pert_efrac
    use cv_stash_flg_mod, only: set_convection_output_flags
    use cv_set_dependent_switches_mod, only: cv_set_dependent_switches
    use dust_parameters_mod, only: i_dust, i_dust_off, i_dust_flux,        &
         dust_veg_emiss, us_am, sm_corr, horiz_d, l_fix_size_dist,         &
         l_twobin_dust, dust_parameters_load,                              &
         dust_parameters_unload
    use electric_inputs_mod, only: electric_method, no_lightning, em_mccaul, &
         k1, k2, gwp_thresh, tiwp_thresh, storm_definition, graupel_and_ice
    use fsd_parameters_mod, only: f_cons
    use glomap_clim_option_mod, only: i_glomap_clim_setup,                     &
                                      l_glomap_clim_aie2,                      &
                                      i_glomap_clim_tune_bc
    use g_wave_input_mod, only: ussp_launch_factor, wavelstar, l_add_cgw,  &
         cgw_scale_factor, i_moist, scale_aware, middle, var
    use mphys_bypass_mod, only: mphys_mod_top,                               &
         qcf2_idims_start, qcf2_idims_end, qcf2_jdims_start, qcf2_jdims_end, &
         qcf2_kdims_end
    use mphys_constants_mod, only: cx, constp
    use mphys_inputs_mod, only: ai, ar, bi, c_r_correl, ci_input, cic_input, &
        di_input, dic_input, i_mcr_iter, l_diff_icevt,                       &
        l_mcr_qrain, l_psd, l_rain, l_warm_new, timestep_mp_in, x1r, x2r,    &
        sediment_loc, i_mcr_iter_tstep, all_sed_start,                       &
        check_run_precip, graupel_option, no_graupel, gr_srcols, a_ratio_exp,&
        a_ratio_fac, l_droplet_tpr, qclrime, l_shape_rime, ndrop_surf,       &
        z_surf, l_fsd_generator, mp_dz_scal, l_subgrid_qcl_mp, aut_qc,       &
        l_mphys_nonshallow, casim_iopt_act, casim_iopt_inuc,                 &
        casim_aerosol_couple_choice,l_casim,                                 &
        casim_aerosol_process_level,casim_moments_choice,                    &
        casim_aerosol_option,                                                &
        l_separate_process_rain, l_mcr_qcf2,                                 &
        l_mcr_qgraup, casim_max_sed_length, fixed_number, wvarfac,           &
        l_orograin, l_orogrime, l_orograin_block,                            &
        fcrit, nsigmasf, nscalesf, l_progn_tnuc, mp_czero, mp_tau_lim
    use mphys_psd_mod, only: x1g, x2g, x4g, x1gl, x2gl, x4gl
    use mphys_switches, only: set_mphys_switches,            &
        max_step_length, max_sed_length,                     &
        iopt_inuc, iopt_act, process_level, l_separate_rain, &
        l_ukca_casim, l_abelshipway, l_warm,                 &
        l_cfrac_casim_diag_scheme, l_prf_cfrac
    use murk_inputs_mod, only: l_murk_advect
    use casim_switches, only: its, ite, jts, jte, kts, kte,              &
                              ils, ile, jls, jle, kls, kle,              &
                              irs, ire, jrs, jre, krs, kre,              &
                              casim_moments_option, n_casim_tracers,     &
                              l_casim_warm_only,                         &
                              l_ukca_aerosol, no_aerosol_modes
    use casim_set_dependent_switches_mod, only:                                &
          casim_set_dependent_switches,                                        &
          casim_print_dependent_switches
    use casim_parent_mod, only: casim_parent, parent_um
    use initialize, only: mphys_init
    use generic_diagnostic_variables, only: casdiags
    use pc2_constants_mod, only: i_cld_off, i_cld_smith, i_cld_pc2,        &
         i_cld_bimodal, rhcpt_off, acf_off, real_shear, rhcpt_tke_based,   &
         pc2eros_exp_rh,pc2eros_hybrid_sidesonly,                          &
         original_but_wrong, acf_cusack, cbl_and_cu, pc2init_smith,        &
         pc2init_logic_original, pc2init_bimodal, i_pc2_homog_g_cf
    use rad_input_mod, only: two_d_fsd_factor
    use science_fixes_mod, only:  i_fix_mphys_drop_settle, second_fix,      &
         l_pc2_homog_turb_q_neg, l_fix_ccb_cct, l_fix_conv_precip_evap,     &
         l_fix_dyndiag, l_fix_pc2_cnv_mix_phase, l_fix_riming,              &
         l_fix_tidy_rainfracs, l_fix_zh, l_fix_incloud_qcf,                 &
         l_fix_mcr_frac_ice, l_fix_gr_autoc, l_improve_cv_cons,             &
         l_pc2_checks_sdfix
    use solinc_data, only: l_skyview
    use stochastic_physics_run_mod, only: a_ent_1_rp, a_ent_shr_rp,         &
         a_ent_shr_rp_max, cbl_mix_fac_rp, cs_rp, g0_rp, g1_rp,             &
         i_rp_scheme, l_rp2, l_rp2_cycle_in, l_rp2_cycle_out,               &
         lai_mult_rp_max, lai_mult_rp_min, lambda_min_rp, ndrop_surf_rp,    &
         lai_mult_rp, par_mezcla_rp, ran_max, ricrit_rp, rp2_callfreq,      &
         rp2_cycle_tm, rp2_decorr_ts, rp_max_idx, rp_min_idx,               &
         z0hm_pft_rp, z0hm_pft_rp_max, z0hm_pft_rp_min
    use turb_diff_mod, only: l_subfilter_horiz, l_subfilter_vert,        &
         mix_factor, turb_startlev_vert, turb_endlev_vert, l_leonard_term
    use ukca_mode_setup, only: ukca_mode_sussbcocdu_7mode, i_ukca_bc_tuned
    use glomap_clim_mode_setup_interface_mod, only:                            &
        glomap_clim_mode_setup_interface
    use ukca_config_specification_mod, only: i_sussbcocdu_7mode
    use ukca_option_mod, only: l_ukca, l_ukca_plume_scav, mode_aitsol_cvscav, &
         l_ukca_aie2, l_ukca_dust
    use ukca_scavenging_mod, only: ukca_mode_scavcoeff

    implicit none

    integer(i_def) :: k, n, n_pft
    logical(l_def) :: l_fix_nacl_density
    logical(l_def) :: l_fix_ukca_hygroscopicities
    logical(l_def) :: dust_loaded = .false.

    ! ----------------------------------------------------------------
    ! UM aerosol scheme settings
    ! For GLOMAP-mode climatology scheme (GLOMAP-clim), these are
    ! contained in glomap_clim_option_mod
    ! For prognostic GLOMAP-mode aerosols (using UKCA with aerosol
    ! precursor chemistry represented by an Offline Oxidants scheme),
    ! these are set via a UKCA API call in subroutine um_ukca_init.
    ! ----------------------------------------------------------------
    if ( aerosol == aerosol_um ) then

      if ( l_radaer .and. .not. ( radiation == radiation_socrates ) ) then
        !  RADAER only calculates required input fields for Socrates.
        !  There is no other output of RADAER, so this setting
        !  is not recommended.
        write(log_scratch_space,'(A)')                                         &
          "It is not recommended to run RADAER without Socrates."
        call log_event( log_scratch_space, LOG_LEVEL_ERROR)
      end if

      ! Options which are bespoke to the aerosol scheme chosen
      select case (glomap_mode)

        case(glomap_mode_climatology, glomap_mode_dust_and_clim)
          ! Set up the correct mode and components for GLOMAP-mode:
          ! 5 mode with SU SS OM BC components
          i_glomap_clim_setup = i_sussbcocdu_7mode
          l_fix_nacl_density = .true.
          l_fix_ukca_hygroscopicities = .false.
          i_glomap_clim_tune_bc = i_ukca_bc_tuned
          call glomap_clim_mode_setup_interface( i_glomap_clim_setup,          &
                                                 l_radaer,                     &
                                                 i_glomap_clim_tune_bc,        &
                                                 l_fix_nacl_density,           &
                                                 l_fix_ukca_hygroscopicities )
        case(glomap_mode_ukca)
          ! UKCA initialisation (via a call to um_ukca_init) is deferred
          ! until after that for JULES since JULES settings are required
          ! for configuring dry deposition.

        case(glomap_mode_radaer_test)
          ! Set up the correct mode and components for RADAER:
          ! 7 mode with SU SS OM BC DU components
          i_glomap_clim_setup = i_sussbcocdu_7mode
          l_fix_nacl_density = .true.
          l_fix_ukca_hygroscopicities = .false.
          i_glomap_clim_tune_bc = i_ukca_bc_tuned
          call glomap_clim_mode_setup_interface( i_glomap_clim_setup,          &
                                                 l_radaer,                     &
                                                 i_glomap_clim_tune_bc,        &
                                                 l_fix_nacl_density,           &
                                                 l_fix_ukca_hygroscopicities )
        case(glomap_mode_off)
          ! Do Nothing

        case default
          write( log_scratch_space, '(A,I0)' )                                 &
             'Invalid aerosol option, stopping', glomap_mode
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )

      end select

      ! Initialisation of RADAER fields
      if ( l_radaer ) then
        n_radaer_mode = 6_i_def
      else
        n_radaer_mode = 0_i_def
      end if

      if ( murk ) then
        l_murk_advect = murk_prognostic
      end if

    else ! if ( aerosol == aerosol_um ) then
      ! Initialisation of RADAER fields
      n_radaer_mode = 0_i_def
    end if ! if ( aerosol == aerosol_um ) then

    if (easyaerosol_sw) then
       n_aer_mode_sw = n_radaer_mode + 1_i_def
    else
       n_aer_mode_sw = n_radaer_mode
    end if
    if (easyaerosol_lw) then
       n_aer_mode_lw = n_radaer_mode + 1_i_def
    else
       n_aer_mode_lw = n_radaer_mode
    end if

    ! Initialisation of aerosol optical property fields
    mode_dimen   = max(  n_aer_mode_sw, n_aer_mode_lw, 1_i_def )
    sw_band_mode = max( (n_aer_mode_sw*n_sw_band) , 1_i_def )
    lw_band_mode = max( (n_aer_mode_lw*n_lw_band) , 1_i_def )

    ! ----------------------------------------------------------------
    ! UM boundary layer scheme settings - contained in UM module bl_option_mod
    ! ----------------------------------------------------------------
    if ( boundary_layer == boundary_layer_um ) then

      if ( surface /= surface_jules ) then
        write( log_scratch_space, '(A)' )                                   &
            'Jules surface is required for UM boundary layer - please switch on'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      a_ent_shr_nml = real(a_ent_shr, r_bl)
      a_ent_2       = real(a_ent_2_in, r_bl)
      bl_levels     = bl_levels_in
      if(allocated(alpha_cd))deallocate(alpha_cd)
      allocate(alpha_cd(bl_levels))
      alpha_cd      = 1.5_r_um
      alpha_cd(1)   = 2.0_r_um
      bl_res_inv    = on

      select case (cbl_opt)
        case(cbl_opt_conventional)
          cbl_op = lem_conven
        case(cbl_opt_standard)
          cbl_op = lem_std
        case(cbl_opt_adjustable)
          cbl_op = lem_adjust
      end select

      dec_thres_cloud = real(dec_thres_cloud_in, r_bl)
      dec_thres_cu = 0.5_r_bl * dec_thres_cloud
      entr_smooth_dec = on

      select case (flux_bc_opt_in)
        case(flux_bc_opt_interactive)
          flux_bc_opt = interactive_fluxes
        case(flux_bc_opt_specified_scalars)
          flux_bc_opt = specified_fluxes_only
        case(flux_bc_opt_specified_scalars_tstar)
          flux_bc_opt = specified_fluxes_tstar
      end select

      flux_grad = off

      if (fric_heating_in) then
        fric_heating = on
      else
        fric_heating = off
      end if

      i_bl_vn = i_bl_vn_9c

      select case (dyn_diag)
        case(dyn_diag_zi_l_sea)
          idyndiag = DynDiag_ZL_corrn
        case(dyn_diag_zi_l_cu)
          idyndiag = DynDiag_ZL_CuOnly
        case(dyn_diag_ri_based)
          idyndiag = DynDiag_Ribased
      end select
      near_neut_z_on_l = real(near_neut_z_on_l_in, r_bl)

      ! Interpolate the vertical gradients of sl,qw and calculate
      ! stability dbdz and Kh on theta-levels
      select case (interp_local)
        case(interp_local_gradients)
          i_interp_local = i_interp_local_gradients
        case(interp_local_cf_dbdz)
          i_interp_local = i_interp_local_cf_dbdz
        end select

      select case (reduce_fa_mix)
        case(reduce_fa_mix_inv_and_cu_lcl)
          keep_ri_fa = on
        case(reduce_fa_mix_inv_only)
          keep_ri_fa = except_disc_inv
      end select

      kprof_cu = buoy_integ
      l_noice_in_turb = noice_in_turb
      l_new_kcloudtop   = new_kcloudtop
      l_reset_dec_thres = .true.
      lambda_min_nml    = 40.0_r_um

      select case (free_atm_mix)
        case(free_atm_mix_to_sharp)
          local_fa = to_sharp_across_1km
        case(free_atm_mix_ntml_corrected)
          local_fa = ntml_level_corrn
        case(free_atm_mix_free_trop_layer)
          local_fa = free_trop_layers
      end select

      pstb = 2.0_r_um
      puns = real(p_unstable, r_um)

      if (relax_sc_over_cu_in) then
        relax_sc_over_cu = on
      else
        relax_sc_over_cu = off
      end if

      ritrans = 0.1_r_bl

      select case (sbl_opt)
        case(sbl_opt_sharpest)
          sbl_op = sharpest
        case(sbl_opt_sharp_sea_mes_land)
          sbl_op = sharp_sea_mes_land
      end select

      select case (sg_orog_mixing_in)
        case(sg_orog_mixing_none)
          sg_orog_mixing = off
        case(sg_orog_mixing_shear_plus_lambda)
          sg_orog_mixing = sg_shear_enh_lambda
      end select

      ! Switch for alternative TKE and variance diagnostics
      var_diags_opt = split_tke_and_inv
      tke_diag_fac  = 1.0_r_bl
      l_use_var_fixes = .true.
      zhloc_depth_fac = real(zhloc_depth_fac_in, r_bl)

      if (topography == topography_horizon) then
        ! Set control logical for use of skyview factor in JULES
        l_skyview = .true.
      end if

    end if

    ! ----------------------------------------------------------------
    ! UM convection scheme settings - contained in UM module cv_run_mod
    ! ----------------------------------------------------------------

    ! The following are needed by conv_diag regardless of whether
    ! convection is actually called or not
    ! Possibly these should vary with vertical level set??
    cape_bottom          = 5
    cape_top             = 50
    cldbase_opt_dp       = 8
    cldbase_opt_sh       = 0
    cvdiag_inv           = 0
    cvdiag_sh_wtest      = 0.02_r_um
    dil_plume_water_load = 0
    ent_fac_dp           = 1.0_r_um
    iconv_congestus      = 0
    iconv_deep           = 0
    icvdiag              = 1
    l_jules_flux         = .true.
    limit_pert_opt       = 2
    plume_water_load     = 0
    tv1_sd_opt           = 2
    w_cape_limit         = 0.4_r_um
    l_reset_neg_delthvu  = .true.

    if ( convection == convection_um ) then

      ! Options needed by all convection schemes
      l_param_conv = .true.
      fac_qsat     = 0.350_r_um
      mparwtr      = 1.0000e-3_r_um
      qlmin        = 4.0000e-4_r_um

      ! Options which are bespoke to the choice of scheme
      select case (cv_scheme)

      case(cv_scheme_gregory_rowntree)

      if ( boundary_layer /= boundary_layer_um ) then
        write( log_scratch_space, '(A)' )                                   &
            'UM boundary layer is required for GR convection - please switch on'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

        i_convection_vn     = i_convection_vn_6a
        adapt               = 8
        adv_conv_prog_dtheta  = .false.
        adv_conv_prog_dq      = .false.
        amdet_fac           = 3.0_r_um
        anv_opt             = 0
        anvil_factor        = 1.0000_r_um
        bl_cnv_mix          = 1
        cape_min            = 0.5_r_um
        cape_timescale      = cape_timescale_in
        cape_ts_max         = 14400.0_r_um
        cape_ts_min         = 1800.0_r_um
        cca2d_dp_opt        = 2
        cca2d_md_opt        = 2
        cca2d_sh_opt        = 2
        cca_dp_knob         = 0.80_r_um
        cca_md_knob         = 0.80_r_um
        cca_sh_knob         = 0.40_r_um
        ccw_dp_knob         = 1.00_r_um
        ccw_for_precip_opt  = 4
        ccw_md_knob         = 1.00_r_um
        ccw_sh_knob         = 1.00_r_um
        cldbase_opt_md      = 2
        cnv_cold_pools      = 0
        cnv_wat_load_opt    = 0
        cpress_term         = 0.3_r_um
        dd_opt              = 1
        deep_cmt_opt        = 6
        eff_dcff            = 3.0_r_um
        eff_dcfl            = 1.0_r_um
        efrac               = 1.0_r_um
        ent_dp_power        = 1.00_r_um
        ent_fac_md          = 1.00_r_um
        ent_opt_dp          = 7
        ent_opt_md          = 6
        fdet_opt            = 2
        iconv_mid           = 1
        iconv_shallow       = 1
        l_cloud_deep        = .true.
        l_3d_cca            = .true.
        l_anvil             = .true.
        l_ccrad             = .true.
        l_cmt_heating       = .true.
        l_conv_prog_precip  = .true.
        l_conv_prog_dtheta  = .true.
        l_conv_prog_dq      = .true.
        l_cv_conserve_check = .true.
        l_fcape             = .true.
        l_mom               = .true.
        l_prog_pert         = .false.
        l_safe_conv         = .true.
        md_pert_opt         = md_pert_efrac
        mdet_opt_dp         = 1
        mdet_opt_md         = 0
        mid_cmt_opt         = 2
        mid_cnv_pmin        = 10000.00_r_um
        midtrig_opt         = mtrig_ntml
        n_conv_calls        = number_of_convection_substeps
        pr_melt_frz_opt     = 2
        qstice              = 3.5000e-3_r_um
        r_det               = 0.5000_r_um
        rad_cloud_decay_opt = 0
        sh_pert_opt         = 1
        t_melt_snow         = 276.15_r_um
        termconv            = 2
        tice                = 263.1500_r_um
        thpixs_mid          = 0.5_r_um
        tower_factor        = 1.0000_r_um
        ud_factor           = 1.0000_r_um
        tau_conv_prog_precip = 10800.0_r_um
        tau_conv_prog_dtheta = 2700.0_r_um
        tau_conv_prog_dq    =  2700.0_r_um
        prog_ent_grad       = -1.1_r_um
        prog_ent_int        = -2.9_r_um
        prog_ent_min        = 0.5_r_um
        prog_ent_max        = 2.5_r_um
        ent_fac_sh          = 1.0_r_um
        c_mass_sh           = 0.03_r_um
        orig_mdet_fac       = 1.0_r_um

        ! Derived switches and parameters are set here based on the options
        ! above
        call cv_set_dependent_switches( )
        ! Flags for diagnostic output are set here
        call set_convection_output_flags( )

      case(cv_scheme_lambert_lewis)
        i_convection_vn   = i_cv_llcs
        non_local_bl      = off
        ng_stress         = off
        bl_res_inv        = off
        if ( scheme == scheme_pc2 ) then
          ! If pc2, we detrain some cloud
          llcs_cloud_precip = llcs_opt_crit_condens
          llcs_detrain_coef = 0.6_r_um
        else
          ! We just rain everything out
          llcs_cloud_precip = llcs_opt_all_rain
        end if
        llcs_rhcrit       = 0.8_r_um
        llcs_timescale    = 3600.0_r_um

      case default
        write( log_scratch_space, '(A,I5)' )  &
             'Invalid convection scheme option, stopping', cv_scheme
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )

      end select
    else ! convection /= convection_um
      ! Need to set the version of the convection diagnosis that we want to use
      i_convection_vn = i_convection_vn_6a
    end if

    ! Check the contents of the convection parameters module
    call check_run_convection()

    ! ----------------------------------------------------------------
    ! UM convection scheme settings for plume scavenging of UKCA
    ! aerosol tracers - contained in UM module ukca_option_mod
    ! ----------------------------------------------------------------

    if ( aerosol == aerosol_um .and.              &
         ( glomap_mode == glomap_mode_ukca ) .or. &
         ( glomap_mode == glomap_mode_dust_and_clim ) ) then
      l_ukca = .true.
      l_ukca_plume_scav = .true.
      if (l_ukca_plume_scav) then
        mode_aitsol_cvscav = 0.5    ! Plume scavenging fraction for soluble
                                    ! Aitken mode aerosol
        call ukca_mode_scavcoeff()
      end if
    end if

    ! ----------------------------------------------------------------
    ! UM cloud scheme settings - contained in UM module cloud_inputs_mod
    ! ----------------------------------------------------------------
    ! needed for visibility diags even without cloud scheme
    rhcrit(1) = 0.96_r_um
    ice_fraction_method = min_liq_overlap
    i_eacf = not_mixph
    if ( cloud == cloud_um ) then

      if ( moisture_formulation == moisture_formulation_dry ) then
        write( log_scratch_space, '(A)' )                                   &
            'moisture_formulation /= ''dry'' is required for UM cloud'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Options needed by all cloud schemes
      cff_spread_rate = real(cff_spread_rate_in, r_um)
      select case (falliceshear_method_in)
        case(falliceshear_method_constant)
          falliceshear_method = original_but_wrong
        case(falliceshear_method_real)
          falliceshear_method = real_shear
      end select
      select case (rh_crit_opt)
        case(rh_crit_opt_namelist)
          i_rhcpt = rhcpt_off
        case(rh_crit_opt_tke)
          i_rhcpt = rhcpt_tke_based
      end select
      ice_width           = real(ice_width_in, r_um)
      l_bm_ez_subcrit_only= ez_subcrit
      ez_max_bm           = ez_max
      ! l_add_cca_to_mcica is unused in LFRic, its functionality
      ! ... being replaced by the cloud_representation option in
      ! ... the radiation namelist (T=combined, F=liquid_and_ice).
      l_subgrid_qv               = subgrid_qv
      rhcrit(1:number_of_layers) = real(rh_crit, r_um)
      ! Options which are bespoke to the choice of scheme
      select case (scheme)

      case(scheme_smith)
        i_cld_vn   = i_cld_smith
        forced_cu  = off
        i_cld_area = acf_cusack

      case(scheme_pc2)
        i_cld_vn                     = i_cld_pc2
        allicetdegc                  = -20.0_r_um
        dbsdtbs_turb_0               = 1.50e-4_r_um
        forced_cu                    = cbl_and_cu
        forced_cu_fac                = 0.5_r_um
        i_cld_area                   = acf_off
        i_pc2_checks_cld_frac_method = 2
        i_pc2_conv_coupling          = 3
        i_pc2_erosion_method         = pc2eros_hybrid_sidesonly
        i_pc2_homog_g_method         = i_pc2_homog_g_cf
        l_ensure_min_in_cloud_qcf    = .false.
        i_pc2_init_logic             = pc2init_logic_original
        starticetkelvin              = 263.15_r_um
        if (pc2ini == pc2ini_smith)   i_pc2_init_method = pc2init_smith
        if (pc2ini == pc2ini_bimodal) i_pc2_init_method = pc2init_bimodal
        l_pc2_implicit_erosion       = .true.

      case(scheme_bimodal)
        i_cld_vn   = i_cld_bimodal
        forced_cu  = off
        i_cld_area = acf_off

      case default
        write( log_scratch_space, '(A,I3)' )  &
             'Invalid cloud scheme option, stopping', scheme
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )

      end select

      ! Check the contents of the cloud parameters module
      call check_run_cloud()

    else ! cloud /= cloud_um
      ! Set switch for no cloud scheme in UM
      i_cld_vn  = i_cld_off
      forced_cu = off
    end if

    ! ----------------------------------------------------------------
    ! Classic dust scheme - contained in dust_parameters_mod
    ! ----------------------------------------------------------------
    ! This is not used in LFRic but potentially called from UM code.
    ! Hence its inputs and options need setting according to the
    ! scheme being either off or running in diagnostic mode to calculate
    ! dust emissions only if these are potentially required by UKCA.

    if ( aerosol == aerosol_um .and.              &
         ( glomap_mode == glomap_mode_ukca ) .or. &
         ( glomap_mode == glomap_mode_dust_and_clim ) ) then
      i_dust = i_dust_flux
      dust_veg_emiss = 1
      us_am = 1.45              ! Increases friction velocity
      sm_corr = 0.5             ! Reduces soil moisture
      horiz_d = 2.25
      l_fix_size_dist = .false.
      l_twobin_dust = .false.
    else
      i_dust = i_dust_off
    end if
    if( dust_loaded ) then
      call dust_parameters_unload( )
      dust_loaded = .false.
    end if
    call dust_parameters_load()
    dust_loaded = .true.

    ! ----------------------------------------------------------------
    ! UM microphysics settings
    ! ----------------------------------------------------------------

    ! The following are needed by the bimodal cloud scheme, hence we initialise
    ! them even when microphysics isn't used
    ai             = 2.5700e-2_r_um
    bi             = 2.00_r_um
    cx(84)         = 1.0_r_um
    constp(35)     = 1.0_r_um
    mp_czero       = 10.0_r_um
    mp_tau_lim     = 1200.0_r_um
    ! The following are needed for the visibility diagnostic, hence we
    ! initialise them even when microphysics isn't used.  They are used in
    ! beta_precip which can still have convective rain/snow.
    x1r            = 2.2000e-1_r_um
    x2r            = 2.2000_r_um
    ! The following are also needed by COSP. They are initialised here but
    ! may be reset within microphysics.
    x1g = x1gl
    x2g = x2gl
    x4g = x4gl

    if ( microphysics == microphysics_um ) then

      if ( cloud /= cloud_um ) then
        write( log_scratch_space, '(A)' )                                   &
            'UM cloud is required for UM microphysics - please switch on'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      l_mcr_qrain    = .true.
      l_mphys_nonshallow = .true.
      l_rain         = .true.

      ! Domain top used in microphysics - contained in mphys_bypass_mod
      mphys_mod_top  = real(domain_top, r_um)

      ! Options only relevent to old microphysics scheme
      if (.not. microphysics_casim) then

        ! Set l_ukca_dust=true as only one dust option and
        ! prog_tnuc so lsp_prognostic_tnuc_kernel can interface
        !  with UM code. prog_tnuc is set from namelist
        l_ukca_dust  = .true.
        l_progn_tnuc = prog_tnuc

        ! These need to be set so that the cdnc value from ukca or glomap_clim
        ! is used by the microphysics scheme. No other option is available
        ! here
        ! The aie1 equivalents are not required for radiation, which is
        ! controlled via the droplet_effective_radius namelist option
        l_ukca_aie2 = .true.
        l_glomap_clim_aie2 = .true.

        select case (graupel_scheme)
        case (graupel_scheme_none)
          graupel_option = no_graupel
        case (graupel_scheme_modified)
          graupel_option = gr_srcols
          nummr_to_transport = 5_i_def
        end select

        a_ratio_exp    = real(a_ratio_exp_in, r_um)
        a_ratio_fac    = real(a_ratio_fac_in, r_um)
        ar             = 1.00_r_um
        c_r_correl     = 0.9_r_um
        ci_input       = 14.3_r_um
        cic_input      = 1024.0_r_um
        di_input       = 0.416_r_um
        dic_input      = 1.0_r_um
        i_mcr_iter     = i_mcr_iter_tstep
        l_diff_icevt   = .true.
        l_droplet_tpr  = droplet_tpr
        l_fsd_generator= cld_fsd_hill
        l_psd          = .true.
        l_shape_rime   = shape_rime
        l_subgrid_qcl_mp = turb_gen_mixph
        l_warm_new     = .true.
        mp_dz_scal     = 2.0_r_um
        ndrop_surf     = real(ndrop_surf_in, r_um)
        qclrime        = real(qcl_rime, r_um)
        sediment_loc   = all_sed_start
        timestep_mp_in = 120
        z_surf         = real(z_surf_in, r_um)
        aut_qc         = 2.47_r_um
        !     Needed by the Seeder Feeder scheme
        l_orograin     = orog_rain
        l_orogrime     = orog_rime
        l_orograin_block = orog_block
        nsigmasf       = real(nsigmasf_in, r_um)
        nscalesf       = real(nscalesf_in, r_um)
        fcrit          = real(fcrit_in, r_um)

      end if

      ! UM options needed if CASIM is being used
      if (microphysics_casim) then

        l_casim = .true.
        l_psd          = .false.
        graupel_option = 2_i_um
        l_mcr_qcf2 = .true.
        l_mcr_qgraup = .true.
        wvarfac = 1.0_r_um

        ! Transport all prognostic variables (include graupel and snow)
        nummr_to_transport = 6_i_def


        casim_moments_option = 22222   ! all double moment
        casim_iopt_act = 0_i_um     ! 'fixed number'
        casim_aerosol_option = 0_i_um   ! no soluble or insoluble aerosol modes
        casim_aerosol_process_level = 0_i_um
        casim_aerosol_couple_choice = 0_i_um
        l_ukca_aerosol = .false.

        casim_moments_choice = 1_i_um
        CALL casim_set_dependent_switches

        ! Tell CASIM that its parent model is the UM. This allows for any UM-specific
        ! operations to take place within CASIM.
        casim_parent = parent_um

        ! Tell CASIM that it needs to output rain and snowfall rates
        ! Required on every single UM timestep for JULES
        casdiags % l_surface_rain  = .true.
        casdiags % l_surface_snow  = .true.
        casdiags % l_surface_graup = .true.
        casdiags % l_rainfall_3d = .true.
        casdiags % l_snowonly_3d = .true.

        max_step_length = real(timestep_mp_in, r_um)

        !---------------------------------------------------------------------
        ! Set up microphysics sedimentation substep
        !---------------------------------------------------------------------
        max_sed_length = casim_max_sed_length  ! from parameter

        !---------------------------------------------------------------------
        ! Set up options for activation and deposition
        !---------------------------------------------------------------------
        iopt_act  = casim_iopt_act
        iopt_inuc = casim_iopt_inuc    ! from parameter

        !---------------------------------------------------------------------
        ! Set up microphysics switches
        !---------------------------------------------------------------------
        if (l_casim_warm_only) l_warm = .true.

        ! Set aerosol processing level directly on the CASIM side using the value
        ! obtained from the run_precip namelist.
        process_level = casim_aerosol_process_level

        ! Set up separate rain processing category for active aerosol on the
        ! CASIM side.
        l_separate_rain = l_separate_process_rain

        ! If using UKCA aerosol, set up switch on the CASIM side
        l_ukca_casim = l_ukca_aerosol

        ! Set up option for Abel and Shipway (2007) rain fall speeds
        ! (already default for Wilson and Ballard microphysics)
        l_abelshipway = .false.

        ! Set up options for CASIM cloud fraction scheme
        l_cfrac_casim_diag_scheme = .false.

        ! Set up cloud fraction scheme coupling (Smith or PC2)
        if ( i_cld_vn == i_cld_smith .or. i_cld_vn == i_cld_pc2                &
                                     .or. i_cld_vn == i_cld_bimodal) then
          l_prf_cfrac   = .true.
          l_abelshipway = .true.
        end if ! i_cld_vn == i_cld_smith/pc2/bimodal

        ! Call set mphys_switches with the options passed in from the run_precip
        ! namelist directly to CASIM.
        call set_mphys_switches(casim_moments_option, casim_aerosol_option)

        !---------------------------------------------------------------------
        ! Initialise and allocate the space required (CASIM repository)
        !---------------------------------------------------------------------
        ! Casim is written k-first
        its = 1_i_um
        ite = 1_i_um
        jts = 1_i_um
        jte = 1_i_um
        kts = 1_i_um
        kte = number_of_layers

        ils = its
        ile = ite
        jls = jts
        jle = jte
        kls = kts
        if ( i_cld_vn == i_cld_smith .or. i_cld_vn == i_cld_pc2                &
                                     .or. i_cld_vn == i_cld_bimodal) then
          kle = kte
        else
          kle = kte - 2
        end if ! i_cld_vn

        irs = its
        ire = ite
        jrs = jts
        jre = jte
        krs = kts
        if ( i_cld_vn == i_cld_smith .or. i_cld_vn == i_cld_pc2                &
                                     .or. i_cld_vn == i_cld_bimodal) then
          kre = kte
        else
          kre = kte - 2
        end if

        call mphys_init( its, ite, jts, jte, kts, kte,                         &
                         ils, ile, jls, jle, kls, kle,                         &
                         l_tendency=.false. )

      end if ! microphysics_casim
    end if ! microphysics == microphysics_um

    !---------------------------------------------------------
    ! UM electric (lightning) scheme settings
    !---------------------------------------------------------

    if ( electric == electric_um ) then
      ! Turn on mcccaul lightning scheme
      electric_method = em_mccaul
      storm_definition = graupel_and_ice
      gwp_thresh = 0.2_r_um
      tiwp_thresh = 1.0_r_um

      if (l_mcr_qcf2) then
        qcf2_idims_start = 1_i_um
        qcf2_jdims_start = 1_i_um
        qcf2_idims_end = 1_i_um
        qcf2_jdims_end = 1_i_um
        qcf2_kdims_end = model_levels
      end if

      if ( microphysics_casim ) then
        ! Use CASIM microphysics-specific settings
        k1 = 0.21_r_um
        k2 = 0.60_r_um
      else
        ! Use Wilson-Ballard microphysics settings
        k1 = 0.042_r_um
        k2 = 0.20_r_um
      end if

    else
      ! Switch off the lightning scheme
      electric_method = no_lightning
    end if


    if ( microphysics == microphysics_um                                    &
         .or. radiation == radiation_socrates ) then
      ! Options for the subgrid cloud variability parametrization used
      ! in microphysics and radiation but living elsewhere in the UM
      ! ... contained in rad_input_mod
      two_d_fsd_factor = two_d_fsd_factor_in
      ! ... contained in fsd_parameters_mod

      if ( cld_fsd_hill ) then
        ! Parameters for fractional standard deviation (fsd) of condensate taken
        ! from part of equation 3 in Hill et al (2015) DOI: 10.1002/qj.2506,
        ! i.e. phi(x,c) = R21 (xc) ^ 1/3 { (0.016 xc)^2.76 + 1 } ^ -0.09
        f_cons(1)      =  0.016
        f_cons(2)      =  2.76
        f_cons(3)      = -0.09
      end if

    end if

    ! Check the contents of the microphysics parameters module
    call check_run_precip()

    ! ----------------------------------------------------------------
    ! UM spectral gravity wave drag options - contained in g_wave_input_mod
    ! ----------------------------------------------------------------
    if ( spectral_gwd == spectral_gwd_um ) then

      cgw_scale_factor = real(cgw_scale_factor_in, r_um)
      l_add_cgw = add_cgw_in
      ussp_launch_factor = real(ussp_launch_factor_in, r_um)
      wavelstar = real(wavelstar_in, r_um)

    end if

    if ( orographic_drag == orographic_drag_um ) then
      scale_aware = .false.
      middle = 0.42_r_um
      var = 0.18_r_um
      select case (include_moisture)
        case(include_moisture_dry)
          i_moist = 0
        case(include_moisture_lowmoist)
          i_moist = 1
        case(include_moisture_moist)
          i_moist = 2
      end select
    end if

    ! ----------------------------------------------------------------
    ! Temporary logicals used to fix bugs in the UM - contained in science_fixes
    ! ----------------------------------------------------------------
    i_fix_mphys_drop_settle = second_fix ! This is a better fix than the
                                         ! original one.
    l_fix_ccb_cct           = .true.
    l_fix_conv_precip_evap  = .true.
    l_fix_dyndiag           = .true.
    l_fix_gr_autoc          = .true.
    l_fix_incloud_qcf       = .true.
    l_fix_mcr_frac_ice      = .true.
    l_fix_pc2_cnv_mix_phase = .true.
    l_fix_riming            = .true.
    l_fix_tidy_rainfracs    = .true.
    l_fix_zh                = .true.
    l_improve_cv_cons       = .true.
    l_pc2_checks_sdfix      = .true.
    l_pc2_homog_turb_q_neg  = .true.

    !-----------------------------------------------------------------------
    ! Smagorinsky mixing options - contained in turb_diff_mod and
    !                              turb_diff_ctl_mod
    !-----------------------------------------------------------------------
    if ( smagorinsky ) then

      ! The following are needed regardless of which mixing option is used
      mix_factor = real(mix_factor_in, r_um)
      turb_startlev_vert  = 2
      turb_endlev_vert    = bl_levels

      ! Options which are bespoke to the choice of scheme
      select case ( mixing_method )

      case( method_3d_smag )
        l_subfilter_horiz = .true.
        l_subfilter_vert  = .true.
        blending_option   = off
        non_local_bl      = off
        ng_stress         = off
      case( method_2d_smag )
        l_subfilter_horiz = .true.
        l_subfilter_vert  = .false.
        blending_option   = off
      case( method_blend_smag_fa )
        l_subfilter_horiz = .true.
        l_subfilter_vert  = .true.
        blending_option   = blend_allpoints
      case( method_blend_1dbl_fa )
        l_subfilter_horiz = .true.
        l_subfilter_vert  = .true.
        blending_option   = blend_gridindep_fa
      end select

    else ! not Smagorinsky

      ! Switches for Smagorinsky being off
      blending_option   = off
      l_subfilter_horiz = .false.
      l_subfilter_vert  = .false.

    end if

    !-----------------------------------------------------------------------
    ! Leonard terms on or off
    !-----------------------------------------------------------------------
    l_leonard_term = leonard_term

    !-----------------------------------------------------------------------
    ! UM Random Parameter scheme settings
    !-----------------------------------------------------------------------
    if ( stochastic_physics == stochastic_physics_um ) then

      ! Switches for the RP scheme
      l_rp2 = use_random_parameters
      i_rp_scheme = 1_i_def ! Use RP2b scheme
      l_rp2_cycle_in = rp_cycle_in
      l_rp2_cycle_out = rp_cycle_out
      rp2_cycle_tm = rp_cycle_tm

      ! Settings for RP algorithm
      rp2_callfreq = rp_callfreq
      rp2_decorr_ts = rp_decorr_ts
      ran_max = rp_ran_max

      ! Parameters
      a_ent_1_rp = rp_bl_a_ent_1
      a_ent_shr_rp = rp_bl_a_ent_shr
      a_ent_shr_rp_max = rp_bl_a_ent_shr_max
      cbl_mix_fac_rp = rp_bl_cbl_mix_fac
      cs_rp = rp_bl_smag_coef
      g0_rp = rp_bl_stable_ri_coef
      g1_rp = rp_bl_cld_top_diffusion
      lambda_min_rp = rp_bl_min_mix_length
      ndrop_surf_rp = rp_mp_ndrop_surf
      par_mezcla_rp = rp_bl_neutral_mix_length
      ricrit_rp = rp_bl_ricrit

      n_pft = size(rp_lsfc_z0hm_pft)
      do n = 1, n_pft
        lai_mult_rp(n) = rp_lsfc_lai_mult(n)
        z0hm_pft_rp(n) = rp_lsfc_z0hm_pft(n)
      end do

    end if ! if ( stochastic_physics == stochastic_physics_um ) then


  end subroutine um_physics_init

end module um_physics_init_mod
