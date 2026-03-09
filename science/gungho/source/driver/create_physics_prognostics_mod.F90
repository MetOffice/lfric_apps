!-------------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief create physics prognostics
!> @details Creates the physics prognostic fields
module create_physics_prognostics_mod

  use energy_correction_config_mod,   only : encorr_usage,                      &
                                             encorr_usage_none
  use clock_mod,                      only : clock_type
  use constants_mod,                  only : i_def, l_def
  use field_mod,                      only : field_type
  use integer_field_mod,              only : integer_field_type
  use field_spec_mod,                 only : main_coll_dict,                    &
                                             adv_coll_dict,                     &
                                             processor_type,                    &
                                             make_spec,                         &
                                             if_adv => if_advected
  use field_mapper_mod,               only : field_mapper_type
  use field_maker_mod,                only : field_maker_type
  use lfric_xios_diag_mod,            only : field_is_enabled
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W2, W3, Wtheta, W2H
  use function_space_mod,             only : function_space_type
  use log_mod,                        only : log_event,                         &
                                             LOG_LEVEL_INFO,                    &
                                             LOG_LEVEL_WARNING,                 &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use mixing_config_mod,              only : smagorinsky
  use physics_config_mod,             only : stochastic_physics_placement,      &
                                             stochastic_physics_placement_fast
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use radiation_config_mod,           only : n_radstep, cloud_representation,   &
                                             cloud_representation_combined,     &
                                             cloud_representation_conv_strat_liq_ice, &
                                             cloud_representation_split,        &
                                             topography, topography_flat,       &
                                             topography_horizon,                &
                                             n_horiz_layer, n_horiz_ang,        &
                                             l_inc_radstep
  use cosp_config_mod,                only : l_cosp, n_cosp_step
  use aerosol_config_mod,             only : glomap_mode,                       &
                                             glomap_mode_climatology,           &
                                             glomap_mode_dust_and_clim,         &
                                             glomap_mode_ukca, n_radaer_step,   &
                                             l_radaer, emissions,               &
                                             emissions_GC3, emissions_GC5,      &
                                             easyaerosol_sw, easyaerosol_lw,    &
                                             murk_prognostic, murk
  use section_choice_config_mod,      only : cloud, cloud_um,                   &
                                             aerosol, aerosol_um,               &
                                             radiation, radiation_socrates,     &
                                             boundary_layer,                    &
                                             boundary_layer_um,                 &
                                             electric, electric_um,             &
                                             iau_sst,                           &
                                             surface, surface_jules,            &
                                             orographic_drag,                   &
                                             orographic_drag_um,                &
                                             convection, convection_um,         &
                                             stochastic_physics,                &
                                             stochastic_physics_um
  use cloud_config_mod,               only : scheme,                            &
                                             scheme_pc2
  use convection_config_mod,          only : cv_scheme, cv_scheme_comorph
  use microphysics_config_mod,        only : microphysics_casim
  use multires_coupling_config_mod,   only : coarse_rad_aerosol

  use jules_surface_config_mod,       only : srf_ex_cnv_gust, l_vary_z0m_soil, &
                                             l_urban2t
  use jules_radiation_config_mod,     only : l_albedo_obs, l_sea_alb_var_chl
  use specified_surface_config_mod,   only : surf_temp_forcing, &
                                             surf_temp_forcing_int_flux
  use spectral_gwd_config_mod,        only : add_cgw
  use microphysics_config_mod,        only : turb_gen_mixph
  use derived_config_mod,             only : l_couple_ocean, l_couple_sea_ice
  use chemistry_config_mod,           only : chem_scheme, chem_scheme_none,    &
                                             chem_scheme_strattrop,            &
                                             chem_scheme_strat_test,           &
                                             chem_scheme_offline_ox,           &
                                             l_ukca_ro2_ntp
  use radiative_gases_config_mod,     only : &
    ch4_rad_opt, ch4_rad_opt_ancil, ch4_rad_opt_prognostic, &
    co_rad_opt, co_rad_opt_ancil, co_rad_opt_prognostic, &
    co2_rad_opt, co2_rad_opt_ancil, co2_rad_opt_prognostic, &
    cs_rad_opt, cs_rad_opt_ancil, cs_rad_opt_prognostic, &
    h2_rad_opt, h2_rad_opt_ancil, h2_rad_opt_prognostic, &
    h2o_rad_opt, h2o_rad_opt_ancil, h2o_rad_opt_prognostic, &
    hcn_rad_opt, hcn_rad_opt_ancil, hcn_rad_opt_prognostic, &
    he_rad_opt, he_rad_opt_ancil, he_rad_opt_prognostic, &
    k_rad_opt, k_rad_opt_ancil, k_rad_opt_prognostic, &
    li_rad_opt, li_rad_opt_ancil, li_rad_opt_prognostic, &
    n2_rad_opt, n2_rad_opt_ancil, n2_rad_opt_prognostic, &
    n2o_rad_opt, n2o_rad_opt_ancil, n2o_rad_opt_prognostic, &
    na_rad_opt, na_rad_opt_ancil, na_rad_opt_prognostic, &
    nh3_rad_opt, nh3_rad_opt_ancil, nh3_rad_opt_prognostic, &
    o2_rad_opt, o2_rad_opt_ancil, o2_rad_opt_prognostic, &
    o3_rad_opt, o3_rad_opt_prognostic, &
    rb_rad_opt, rb_rad_opt_ancil, rb_rad_opt_prognostic, &
    so2_rad_opt, so2_rad_opt_ancil, so2_rad_opt_prognostic, &
    tio_rad_opt, tio_rad_opt_ancil, tio_rad_opt_prognostic, &
    vo_rad_opt, vo_rad_opt_ancil, vo_rad_opt_prognostic
  use formulation_config_mod,         only : moisture_formulation,    &
                                             moisture_formulation_dry
  use stochastic_physics_config_mod,  only : blpert_type, blpert_type_off

#ifdef UM_PHYSICS
  use multidata_field_dimensions_mod, only :                                    &
    get_ndata_val => get_multidata_field_dimension
  use cv_run_mod,                     only:  l_conv_prog_precip,                &
                                             l_conv_prog_dtheta,                &
                                             l_conv_prog_dq,                    &
                                             adv_conv_prog_dtheta,              &
                                             adv_conv_prog_dq
  use mphys_inputs_mod, only: casim_iopt_act, l_mcr_precfrac
  use bl_option_mod, only: l_calc_tau_at_p
  use cloud_inputs_mod, only: l_pc2_homog_conv_pressure
  use io_config_mod,                  only : checkpoint_read, checkpoint_write
  use initialization_config_mod,      only : init_option,                       &
                                             init_option_checkpoint_dump
#endif


  implicit none

  private
  public :: create_physics_prognostics, process_physics_prognostics

contains

  !> @brief Iterate over active model fields and apply an arbitrary
  !> processor to the field specifiers.
  !> @details To be used by create_physics_prognostics to create fields and by
  !> gungho_model_mod / initialise_infrastrucure to enable checkpoint fields.
  !> These two operations have to be separate because they have to happen
  !> at different times. The enabling of checkpoint fields has to happen
  !> before the io context closes, and this is too early for field creation.
  !> @param  process Processor to be applied to selected field specifiers
  subroutine process_physics_prognostics(processor)
    use field_spec_mod,            only : main => main_coll_dict,               &
                                          adv => adv_coll_dict
    implicit none
    class(processor_type) :: processor
#ifdef UM_PHYSICS
    logical(l_def) :: checkpoint_flag, checkpoint_GC3, checkpoint_GC5
    logical(l_def) :: checkpoint_couple
    logical(l_def) :: advection_flag
    logical(l_def) :: advection_flag_dust
    logical(l_def) :: is_empty
    logical(l_def) :: is_rad ! Flag for chemistry fields
                             ! that are radiatively active
    logical(l_def) :: sst_pert_flag
#endif

    class(clock_type), pointer :: clock

    clock => processor%get_clock()

    !========================================================================
    ! Fields derived from the FE dynamical fields for use in physics
    !========================================================================

    ! Wtheta fields
    call processor%apply(make_spec('velocity_w2v', main%derived, Wtheta))
    call processor%apply(make_spec('w_in_wth', main%derived, Wtheta))
    call processor%apply(make_spec('rho_in_wth', main%derived, Wtheta))
    call processor%apply(make_spec('wetrho_in_wth', main%derived, Wtheta))
    call processor%apply(make_spec('exner_in_wth', main%derived, Wtheta))
    call processor%apply(make_spec('exner_wth_n', main%derived, Wtheta))

    if ( boundary_layer == boundary_layer_um .or.                              &
         convection     == convection_um ) then

      call processor%apply(make_spec('theta_star', main%derived, Wtheta))

    end if

    if ( boundary_layer == boundary_layer_um .or.                              &
         convection     == convection_um     .or.                              &
         smagorinsky ) then

      call processor%apply(make_spec('shear', main%derived, Wtheta, &
                                     empty = (.not. smagorinsky) ))
      call processor%apply(make_spec('visc_h', main%derived, Wtheta, &
                                     empty = (.not. smagorinsky) ))
      call processor%apply(make_spec('visc_m', main%derived, Wtheta, &
                                     empty = (.not. smagorinsky) ))

    end if

    ! W3 fields
    call processor%apply(make_spec('u_in_w3', main%derived, W3))
    call processor%apply(make_spec('v_in_w3', main%derived, W3))
    call processor%apply(make_spec('w_in_w3', main%derived, W3))
    call processor%apply(make_spec('theta_in_w3', main%derived, W3))
    call processor%apply(make_spec('wetrho_in_w3', main%derived, W3))

    if ( boundary_layer               == boundary_layer_um .or.                &
         convection                   == convection_um     .or.                &
         stochastic_physics_placement == stochastic_physics_placement_fast ) then

      call processor%apply(make_spec('u_in_w3_star', main%derived, W3))
      call processor%apply(make_spec('v_in_w3_star', main%derived, W3))

    end if

    ! W2 fields
    call processor%apply(make_spec('u_physics', main%derived, W2))
    call processor%apply(make_spec('u_star', main%derived, W2))
    call processor%apply(make_spec('wetrho_in_w2', main%derived, W2))

    ! W2H fields
    call processor%apply(make_spec('u_in_w2h', main%derived, W2H))
    call processor%apply(make_spec('v_in_w2h', main%derived, W2H))

    ! 2D fields
    if ( encorr_usage /= encorr_usage_none ) then
      call processor%apply(make_spec('accumulated_fluxes', main%derived, W3,   &
          twod=.true., empty = (encorr_usage == encorr_usage_none) ))
    end if

#ifdef UM_PHYSICS
    !========================================================================
    ! Fields owned by the radiation scheme
    !========================================================================

    ! 2D fields, might need checkpointing
    if (surface == surface_jules .and. l_albedo_obs) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('albedo_obs_vis', main%radiation,           &
        ckp=checkpoint_flag, empty = (.not. l_albedo_obs) ))
    call processor%apply(make_spec('albedo_obs_nir', main%radiation,           &
        ckp=checkpoint_flag, empty = (.not. l_albedo_obs) ))

    ! 3D fields, need checkpointing
    call processor%apply(make_spec('ozone', main%radiation,  ckp=.true.))

    ! 2D fields, don't need checkpointing
    call processor%apply(make_spec('lw_down_surf', main%radiation, W3, twod=.true.))
    call processor%apply(make_spec('sw_down_surf', main%radiation, W3, twod=.true.))
    call processor%apply(make_spec('sw_direct_surf', main%radiation, W3,        &
        twod=.true.))
    call processor%apply(make_spec('sw_down_blue_surf', main%radiation, W3,     &
        twod=.true.))
    call processor%apply(make_spec('sw_direct_blue_surf', main%radiation, W3,   &
        twod=.true.))

    call processor%apply(make_spec('cos_zenith_angle', main%radiation, W3,      &
        twod=.true.))
    call processor%apply(make_spec('lit_fraction', main%radiation, W3, twod=.true.))

    ! 3D fields, don't need checkpointing
    call processor%apply(make_spec('sw_heating_rate', main%radiation, Wtheta))
    call processor%apply(make_spec('lw_heating_rate', main%radiation, Wtheta))
    call processor%apply(make_spec('dtheta_rad', main%radiation, Wtheta))
    call processor%apply(make_spec('dmv_pc2_rad', main%radiation, Wtheta))

    ! Fields on surface tiles, don't need checkpointing
    call processor%apply(make_spec('lw_up_tile', main%radiation, W3,            &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('sw_up_tile', main%radiation, W3,            &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('sw_up_blue_tile', main%radiation, W3,       &
        mult='surface_tiles', twod=.true.))

    ! Fields that need checkpointing for the topographic correction scheme
    if (radiation == radiation_socrates .and.                                   &
        topography /= topography_flat) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('slope_angle', main%radiation,               &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('slope_aspect', main%radiation,              &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('skyview', main%radiation,                   &
        ckp=checkpoint_flag))
    if (radiation == radiation_socrates .and.                                   &
        topography == topography_horizon) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('horizon_angle', main%radiation,             &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('horizon_aspect', main%radiation,            &
          ckp=checkpoint_flag))

    ! Fields which need checkpointing for radiation timestepping
    !
    !>       involving passing the clock down here.
    !>
    if (radiation == radiation_socrates) then
      ! Checkpoint unless both the first timestep of this run and the
      ! first timestep of the next run are radiation timesteps
      checkpoint_flag =                                                        &
        mod(clock%get_first_step()-1, n_radstep) /= 0 .or.                     &
        mod(clock%get_last_step(),    n_radstep) /= 0

      if (checkpoint_read .or.                                                 &
          init_option == init_option_checkpoint_dump) then
        ! If the first timestep of this run IS a radiation timestep, but the
        ! first timestep of the next run IS NOT, then checkpoint_flag
        ! must be false to allow model to start running, as the _rts
        ! prognostics will not be in the initial dump
        if (mod(clock%get_first_step()-1, n_radstep) == 0 .and.                &
            mod(clock%get_last_step(),    n_radstep) /= 0) then
          checkpoint_flag = .false.
          if (checkpoint_write) then
            ! Any dump written will be incomplete, and the following run
            ! will need to start on a radiation timestep - print user a
            ! warning
            call log_event('Danger: start of this run is a radiation ' //      &
                           'timestep, but start of next run is not. ' //       &
                           'Written dump will be incomplete. Next run ' //     &
                           'must start with a radiation timestep',             &
                           LOG_LEVEL_WARNING)
          end if
        endif
      end if
    else
      checkpoint_flag = .false.
    end if

    ! 2D fields
    call processor%apply(make_spec('lw_down_surf_rts', main%radiation,          &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_down_surf_rts', main%radiation,          &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_direct_surf_rts', main%radiation,        &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_down_blue_surf_rts', main%radiation,     &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_direct_blue_surf_rts', main%radiation,   &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('lw_up_surf_rts', main%radiation,            &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_up_surf_rts', main%radiation,            &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('lw_up_toa_rts', main%radiation,             &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_up_toa_rts', main%radiation,             &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_direct_toa_rts', main%radiation,         &
         ckp=checkpoint_flag))

    call processor%apply(make_spec('cos_zenith_angle_rts', main%radiation,      &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('lit_fraction_rts', main%radiation,          &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('stellar_irradiance_rts', main%radiation,    &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('sin_stellar_declination_rts', main%radiation,  &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('stellar_eqn_of_time_rts', main%radiation,   &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('orographic_correction_rts', main%radiation, &
         ckp=checkpoint_flag))

    ! 3D fields
    call processor%apply(make_spec('sw_heating_rate_rts', main%radiation,       &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('lw_heating_rate_rts', main%radiation,       &
        ckp=checkpoint_flag))

    ! Fields on surface tiles
    call processor%apply(make_spec('lw_up_tile_rts', main%radiation,            &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_up_tile_rts', main%radiation,            &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('sw_up_blue_tile_rts', main%radiation,       &
          ckp=checkpoint_flag))


    ! Fields needed for radiation incremental timestepping
    if (l_inc_radstep) then

      ! 2D fields
      call processor%apply(make_spec('lw_down_surf_rtsi', main%radiation,      &
         ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_down_surf_rtsi', main%radiation,      &
         ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_direct_surf_rtsi', main%radiation,    &
         ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_down_blue_surf_rtsi', main%radiation, &
         ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_direct_blue_surf_rtsi', main%radiation,&
         ckp=checkpoint_flag))
      call processor%apply(make_spec('lw_up_surf_rtsi', main%radiation,        &
         ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_up_surf_rtsi', main%radiation,        &
         ckp=checkpoint_flag))
      call processor%apply(make_spec('lw_up_toa_rtsi', main%radiation,         &
        ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_up_toa_rtsi', main%radiation,         &
        ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_direct_toa_rtsi', main%radiation,     &
         ckp=checkpoint_flag))

      ! 3D fields
      call processor%apply(make_spec('sw_heating_rate_rtsi', main%radiation,   &
        ckp=checkpoint_flag))
      call processor%apply(make_spec('lw_heating_rate_rtsi', main%radiation,   &
        ckp=checkpoint_flag))

      ! Fields on surface tiles
      call processor%apply(make_spec('lw_up_tile_rtsi', main%radiation,        &
          ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_up_tile_rtsi', main%radiation,        &
          ckp=checkpoint_flag))
      call processor%apply(make_spec('sw_up_blue_tile_rtsi', main%radiation,   &
          ckp=checkpoint_flag))

    end if

    ! Fields which need checkpointing when superstepping radaer
    ! Checkpoint unless both the first timestep of this run and the
    ! first timestep of the next run are radaer timesteps
    if (l_radaer) then
      checkpoint_flag =                                                        &
        mod(clock%get_first_step()-1, n_radaer_step*n_radstep) /= 0 .or.       &
        mod(clock%get_last_step(),    n_radaer_step*n_radstep) /= 0

      if (checkpoint_read .or.                                                 &
          init_option == init_option_checkpoint_dump) then
        ! If the first timestep of this run IS a radaer timestep, but the
        ! first timestep of the next run IS NOT, then checkpoint_flag
        ! must be false to allow model to start running, as the radaer
        ! prognostics will not be in the initial dump
        if (mod(clock%get_first_step()-1, n_radaer_step*n_radstep) == 0 .and.  &
            mod(clock%get_last_step(),    n_radaer_step*n_radstep) /= 0) then
          checkpoint_flag = .false.
          if (checkpoint_write) then
            call log_event('Danger: start of this run is a radaer ' //         &
                           'timestep, but start of next run is not. ' //       &
                           'Written dump will be incomplete. Next run ' //     &
                           'must start with a radaer timestep',                &
                           LOG_LEVEL_WARNING)
          end if
        endif
      end if
    else
      checkpoint_flag = .false.
    end if

    ! vector_space=>function_space_collection%get_fs(mesh,0,0,Wtheta,
    !     get_ndata_val('aero_modes') )
    call processor%apply(make_spec('aer_mix_ratio', main%radiation,             &
         ckp=checkpoint_flag))
    ! vector_space=>function_space_collection%get_fs(mesh,0,0,Wtheta,
    !     get_ndata_val('sw_bands_aero_modes') )
    call processor%apply(make_spec('aer_sw_absorption', main%radiation,         &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('aer_sw_scattering', main%radiation,         &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('aer_sw_asymmetry', main%radiation,          &
         ckp=checkpoint_flag))
    ! vector_space=>function_space_collection%get_fs(mesh,0,0,Wtheta,
    !     get_ndata_val('lw_bands_aero_modes') )
    call processor%apply(make_spec('aer_lw_absorption', main%radiation,         &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('aer_lw_scattering', main%radiation,         &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('aer_lw_asymmetry', main%radiation,          &
         ckp=checkpoint_flag))

    ! Fields which need checkpointing for time sampling of COSP diagnostics.
    ! Checkpoint unless both the first timestep of this run and the
    ! first timestep of the next run are COSP timesteps
    if (l_cosp) then
      checkpoint_flag = &
        mod(clock%get_first_step()-1, n_cosp_step) /= 0 .or. &
        mod(clock%get_last_step(),    n_cosp_step) /= 0

      if (checkpoint_read .or. init_option == init_option_checkpoint_dump) then
        ! If the first timestep of this run IS a COSP timestep, but the
        ! first timestep of the next run IS NOT, then checkpoint_flag
        ! must be false to allow model to start running, as the COSP
        ! prognostics will not be in the initial dump
        if (mod(clock%get_first_step()-1, n_cosp_step) == 0 .and. &
            mod(clock%get_last_step(),    n_cosp_step) /= 0) then
          checkpoint_flag = .false.
          if (checkpoint_write) then
            call log_event('Danger: start of this run is a COSP ' // &
                           'timestep, but start of next run is not. ' // &
                           'Written dump will be incomplete. Next run ' // &
                           'must start with a COSP timestep', &
                           LOG_LEVEL_WARNING)
          end if
        endif
      end if
    else
      checkpoint_flag = .false.
    end if

    call processor%apply(make_spec('lit_fraction_cosp', main%radiation, &
      W3, twod=.true., ckp=checkpoint_flag))

    !========================================================================
    ! Fields owned by the microphysics scheme
    !========================================================================

    ! 3D fields, need checkpointing
    call processor%apply(make_spec('precfrac', main%microphysics,              &
        adv_coll=if_adv(l_mcr_precfrac, adv%all_adv), ckp=l_mcr_precfrac,      &
        empty = (.not. l_mcr_precfrac) ))

    ! Fields for CASIM (Cloud-AeroSol Interacting Microphysics)
    checkpoint_flag = microphysics_casim
    advection_flag = microphysics_casim

    call processor%apply(make_spec('nl_mphys', main%microphysics,              &
        adv_coll=if_adv((advection_flag .and. casim_iopt_act /= 0_i_def),      &
        adv%last_adv), ckp=checkpoint_flag, empty = (.not. microphysics_casim)))
    call processor%apply(make_spec('nr_mphys', main%microphysics,              &
        adv_coll=if_adv(advection_flag, adv%last_adv), ckp=checkpoint_flag,    &
        empty = (.not. microphysics_casim) ))
    call processor%apply(make_spec('ni_mphys', main%microphysics,              &
        adv_coll=if_adv(advection_flag, adv%last_adv), ckp=checkpoint_flag,    &
        empty = (.not. microphysics_casim) ))
    call processor%apply(make_spec('ns_mphys', main%microphysics,              &
        adv_coll=if_adv(advection_flag, adv%last_adv), ckp=checkpoint_flag,    &
        empty = (.not. microphysics_casim) ))
    call processor%apply(make_spec('ng_mphys', main%microphysics,              &
        adv_coll=if_adv(advection_flag, adv%last_adv), ckp=checkpoint_flag,    &
        empty = (.not. microphysics_casim) ))

    ! 2D fields, don't need checkpointing
    call processor%apply(make_spec('ls_rain', main%microphysics, W3, twod=.true.))
    call processor%apply(make_spec('ls_snow', main%microphysics, W3, twod=.true.))
    call processor%apply(make_spec('lsca_2d', main%microphysics, W3, twod=.true.))

    call processor%apply(make_spec('ls_graup', main%microphysics, W3, twod=.true.))
    call processor%apply(make_spec('tnuc_nlcl', main%microphysics, W3, twod=.true.))

    ! 3D fields, don't need checkpointing
    call processor%apply(make_spec('dtheta_mphys', main%microphysics, Wtheta))
    call processor%apply(make_spec('dmv_mphys', main%microphysics, Wtheta))
    call processor%apply(make_spec('ls_rain_3d', main%microphysics, Wtheta))
    call processor%apply(make_spec('ls_snow_3d', main%microphysics, Wtheta))
    call processor%apply(make_spec('autoconv', main%microphysics, Wtheta))
    call processor%apply(make_spec('accretion', main%microphysics, Wtheta))
    call processor%apply(make_spec('rim_cry', main%microphysics, Wtheta))
    call processor%apply(make_spec('rim_agg', main%microphysics, Wtheta))
    call processor%apply(make_spec('tnuc', main%microphysics, Wtheta))
    !========================================================================
    ! Fields owned by the electric scheme
    !========================================================================

    if ( electric == electric_um ) then

      ! 2D lightning potential field. Assuming needs checkpointing for now.
      ! This field doesn't need to be advected.
      call processor%apply(make_spec('flash_potential', main%electric,         &
           ckp=.true., empty = (electric /= electric_um) ))

    end if

    !========================================================================
    ! Fields owned by the orographic drag schemes
    !========================================================================

    ! 2D fields, might need checkpointing
    if (boundary_layer == boundary_layer_um .or.                                &
         surface == surface_jules           .or.                                &
         orographic_drag == orographic_drag_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('sd_orog', main%orography,                   &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('grad_xx_orog', main%orography,              &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('grad_xy_orog', main%orography,              &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('grad_yy_orog', main%orography,              &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('peak_to_trough_orog', main%orography,       &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('silhouette_area_orog', main%orography,      &
         ckp=checkpoint_flag))

    !========================================================================
    ! Fields owned by the turbulence scheme
    !========================================================================

    ! 2D fields, might need checkpointing

    if (boundary_layer == boundary_layer_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('zh', main%turbulence,                       &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('wvar', main%turbulence,                     &
        ckp=turb_gen_mixph))
    call processor%apply(make_spec('gradrinr', main%turbulence, Wtheta))

    ! 2D fields, don't need checkpointing
    call processor%apply(make_spec('ntml', main%turbulence, W3, twod=.true.,    &
        is_int=.true.))
    call processor%apply(make_spec('cumulus', main%turbulence, W3, twod=.true., &
        is_int=.true.))
    call processor%apply(make_spec('z_lcl', main%turbulence, W3, twod=.true.))
    call processor%apply(make_spec('inv_depth', main%turbulence, W3, twod=.true.))
    call processor%apply(make_spec('qcl_at_inv_top', main%turbulence, W3,       &
        twod=.true.))
    call processor%apply(make_spec('blend_height_tq', main%turbulence, W3,      &
        twod=.true., is_int=.true.))
    call processor%apply(make_spec('zh_nonloc', main%turbulence, W3, twod=.true.))
    call processor%apply(make_spec('zhsc', main%turbulence, W3, twod=.true.))
    call processor%apply(make_spec('bl_weight_1dbl', main%turbulence, W3,       &
        twod=.true.))
    call processor%apply(make_spec('level_ent', main%turbulence, W3, twod=.true.,     &
        is_int=.true.))
    call processor%apply(make_spec('level_ent_dsc', main%turbulence, W3, twod=.true., &
        is_int=.true.))

    ! Space for the 7 BL types
    ! vector_space => function_space_collection%get_fs(twod_mesh, 0, 0, W3,
    !     get_ndata_val('boundary_layer_types'))
    call processor%apply(make_spec('bl_type_ind', main%turbulence, W3,          &
        mult='boundary_layer_types', twod=.true., is_int=.true.))

    ! 3D fields, don't need checkpointing
    call processor%apply(make_spec('bq_bl', main%turbulence, Wtheta))
    call processor%apply(make_spec('bt_bl', main%turbulence, Wtheta))
    call processor%apply(make_spec('lmix_bl', main%turbulence, Wtheta))
    call processor%apply(make_spec('dsldzm', main%turbulence, Wtheta))
    call processor%apply(make_spec('mix_len_bm', main%turbulence, Wtheta))
    call processor%apply(make_spec('tke_bl', main%turbulence, Wtheta))
    call processor%apply(make_spec('rhokm_bl', main%turbulence, Wtheta))
    call processor%apply(make_spec('dtrdz_tq_bl', main%turbulence, Wtheta))
    call processor%apply(make_spec('dw_bl', main%turbulence, Wtheta))
    call processor%apply(make_spec('thetal_inc_leonard', main%turbulence, Wtheta))
    call processor%apply(make_spec('mt_inc_leonard', main%turbulence, Wtheta))

    ! 3D fields on W3 (rho) levels
    call processor%apply(make_spec('moist_flux_bl', main%turbulence, W3))
    call processor%apply(make_spec('heat_flux_bl', main%turbulence, W3))
    call processor%apply(make_spec('rhokh_bl', main%turbulence, W3))
    call processor%apply(make_spec('taux', main%turbulence, W3, &
                         empty=(.not. l_calc_tau_at_p) ))
    call processor%apply(make_spec('tauy', main%turbulence, W3, &
                         empty=(.not. l_calc_tau_at_p) ))

    ! W2 fields, don't need checkpointing
    call processor%apply(make_spec('rhokm_w2', main%turbulence, W2))
    call processor%apply(make_spec('tau_w2', main%turbulence, W2))

    ! Fields on entrainment levels, don't need checkpointing
    ! vector_space => function_space_collection%get_fs(twod_mesh, 0, 0, W3,
    !     get_ndata_val('entrainment_levels'))
    call processor%apply(make_spec('ent_we_lim', main%turbulence, W3,           &
        mult='entrainment_levels', twod=.true.))
    call processor%apply(make_spec('ent_t_frac', main%turbulence, W3,           &
        mult='entrainment_levels', twod=.true.))
    call processor%apply(make_spec('ent_zrzi', main%turbulence, W3,             &
        mult='entrainment_levels', twod=.true.))
    call processor%apply(make_spec('ent_we_lim_dsc', main%turbulence, W3,       &
        mult='entrainment_levels', twod=.true.))
    call processor%apply(make_spec('ent_t_frac_dsc', main%turbulence, W3,       &
        mult='entrainment_levels', twod=.true.))
    call processor%apply(make_spec('ent_zrzi_dsc', main%turbulence, W3,         &
        mult='entrainment_levels', twod=.true.))

    !========================================================================
    ! Fields owned by the convection scheme
    !========================================================================

    ! 2D fields, might need checkpointing
    if (convection == convection_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('conv_rain', main%convection,                &
        ckp=(checkpoint_flag .and. add_cgw)))
    call processor%apply(make_spec('conv_snow', main%convection,                &
        ckp=(checkpoint_flag .and. add_cgw)))
    call processor%apply(make_spec('dd_mf_cb', main%convection,                 &
        ckp=(checkpoint_flag .and. srf_ex_cnv_gust)))

    ! 3D fields, might need checkpointing
    if (convection == convection_um) then
      select case (cloud_representation)
      case (cloud_representation_combined,                                      &
            cloud_representation_conv_strat_liq_ice,                            &
            cloud_representation_split)
        checkpoint_flag = .true.
      case default
        checkpoint_flag = .false.
      end select
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('cca', main%convection,                      &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('ccw', main%convection,                      &
        ckp=checkpoint_flag))

    ! 2D fields, don't need checkpointing
    call processor%apply(make_spec('cca_2d', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('shallow_flag', main%convection, W3, twod=.true.,  &
        is_int=.true.))
    call processor%apply(make_spec('uw0_flux', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('vw0_flux', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('lcl_height', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('parcel_top', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('level_parcel_top', main%convection, W3,     &
        twod=.true., is_int=.true.))
    call processor%apply(make_spec('wstar', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('thv_flux', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('parcel_buoyancy', main%convection, W3,      &
        twod=.true.))
    call processor%apply(make_spec('qsat_at_lcl', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('cape_diluted', main%convection, W3,         &
        twod=.true.))

    ! Convective cloud base and top. Diagnostics in reality but defined as
    ! Prognostics as these are driving fields for UKCA
    call processor%apply(make_spec('cv_base', main%convection, W3, twod=.true.))
    call processor%apply(make_spec('cv_top', main%convection, W3, twod=.true.))

    ! 3D fields, don't need checkpointing
    call processor%apply(make_spec('dt_conv', main%convection, Wtheta))
    call processor%apply(make_spec('dmv_conv', main%convection, Wtheta))
    call processor%apply(make_spec('dmcl_conv', main%convection, Wtheta))
    call processor%apply(make_spec('dms_conv', main%convection, Wtheta))
    call processor%apply(make_spec('dcfl_conv', main%convection, Wtheta))
    call processor%apply(make_spec('dcff_conv', main%convection, Wtheta))
    call processor%apply(make_spec('dbcf_conv', main%convection, Wtheta))
    call processor%apply(make_spec('massflux_up', main%convection, Wtheta))
    call processor%apply(make_spec('massflux_down', main%convection, Wtheta))
    call processor%apply(make_spec('conv_rain_3d', main%convection, Wtheta))
    call processor%apply(make_spec('conv_snow_3d', main%convection, Wtheta))
    call processor%apply(make_spec('pressure_inc_env', main%convection, &
                         Wtheta, empty=(.not. l_pc2_homog_conv_pressure) ))

    ! 3D fields on W3 (rho) levels
    call processor%apply(make_spec('du_conv', main%convection, W3))
    call processor%apply(make_spec('dv_conv', main%convection, W3))

    call processor%apply(make_spec('conv_prog_dtheta', main%convection,        &
        adv_coll=if_adv((l_conv_prog_dtheta .and. adv_conv_prog_dtheta),       &
        adv%all_adv), ckp=l_conv_prog_dtheta, empty=(.not. l_conv_prog_dtheta)))
    call processor%apply(make_spec('conv_prog_dmv', main%convection,           &
        adv_coll=if_adv((l_conv_prog_dq .and. adv_conv_prog_dq), adv%all_adv), &
        ckp=l_conv_prog_dq, empty=(.not. l_conv_prog_dq)))

    call processor%apply(make_spec('conv_prog_precip', main%convection,        &
        adv_coll=if_adv(l_conv_prog_precip, adv%all_adv),                      &
        ckp=l_conv_prog_precip, empty=(.not. l_conv_prog_precip)))

    !========================================================================
    ! Fields owned by the cloud scheme
    !========================================================================

    ! 3D fields, might need checkpointing
    if (cloud == cloud_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('area_fraction', main%cloud,                 &
        ckp=(checkpoint_flag .and. radiation==radiation_socrates)))

    ! 3D fields, might need advecting
    if ( scheme == scheme_pc2 ) then
      advection_flag=.true.
    else
      advection_flag=.false.
    endif
    call processor%apply(make_spec('liquid_fraction', main%cloud,               &
        adv_coll=if_adv(advection_flag, adv%all_adv), ckp=checkpoint_flag))
    call processor%apply(make_spec('frozen_fraction', main%cloud,               &
        adv_coll=if_adv(advection_flag, adv%all_adv), ckp=checkpoint_flag))
    call processor%apply(make_spec('bulk_fraction', main%cloud,                 &
        adv_coll=if_adv(advection_flag, adv%all_adv), ckp=checkpoint_flag))

    call processor%apply(make_spec('rh_crit', main%cloud, Wtheta))
    call processor%apply(make_spec('departure_exner_wth', main%cloud, Wtheta,   &
        adv_coll=if_adv(advection_flag, adv%last_adv)))
    call processor%apply(make_spec('sigma_ml', main%cloud, Wtheta))
    call processor%apply(make_spec('sigma_mi', main%cloud, Wtheta))

    ! Fields for bimodal cloud scheme
    call processor%apply(make_spec('tau_dec_bm', main%cloud, Wtheta))
    call processor%apply(make_spec('tau_hom_bm', main%cloud, Wtheta))
    call processor%apply(make_spec('tau_mph_bm', main%cloud, Wtheta))

    is_empty = (cv_scheme /= cv_scheme_comorph)
    call processor%apply(make_spec('cf_liq_n', main%cloud, Wtheta, &
         empty = is_empty))
    call processor%apply(make_spec('cf_fro_n', main%cloud, Wtheta, &
         empty = is_empty))
    call processor%apply(make_spec('cf_bulk_n', main%cloud, Wtheta, &
         empty = is_empty))

    !========================================================================
    ! Fields owned by the surface exchange scheme
    !========================================================================

    ! 2D fields, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if

    ! Coupling fields might need checkpointing
    if (surface == surface_jules .and.                                          &
       (l_couple_sea_ice .or. l_couple_ocean)) then
      checkpoint_couple = .true.
    else
      checkpoint_couple = .false.
    end if

    call processor%apply(make_spec('z0msea', main%surface,                      &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('surface_conductance', main%surface,         &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('chloro_sea', main%surface,                  &
        ckp=(checkpoint_flag .and. l_sea_alb_var_chl)))

    ! Fields on surface tiles, might need checkpointing
    call processor%apply(make_spec('tile_fraction', main%surface,               &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('tile_temperature', main%surface,            &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('screen_temperature', main%surface,          &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('time_since_transition', main%surface,       &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('canopy_water', main%surface,                &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('tile_lw_grey_albedo', main%surface,         &
          ckp=checkpoint_flag))

    ! vector_space=>function_space_collection%get_fs(twod_mesh, 0, 0, W3,
    !     get_ndata_val('land_tile_rad_band'))
    call processor%apply(make_spec('albedo_obs_scaling', main%surface,         &
        ckp=(checkpoint_flag .and. l_albedo_obs), empty = (.not. l_albedo_obs)))

    ! Fields on plant functional types, might need checkpointing
    call processor%apply(make_spec('leaf_area_index', main%surface,             &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('canopy_height', main%surface,               &
          ckp=checkpoint_flag))


    ! Sea-ice category fields, might need checkpointing
    call processor%apply(make_spec('sea_ice_thickness', main%surface,           &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('sea_ice_temperature', main%surface,         &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('sea_ice_conductivity', main%surface, W3,    &
        mult='sea_ice_categories', twod=.true., ckp=.false.))
    call processor%apply(make_spec('melt_pond_fraction', main%surface, W3,      &
        mult='sea_ice_categories', twod=.true., ckp=.false.))
    call processor%apply(make_spec('melt_pond_depth', main%surface, W3,         &
        mult='sea_ice_categories', twod=.true., ckp=.false.))
    call processor%apply(make_spec('sea_ice_pensolar_frac_direct', main%surface, W3, &
        mult='sea_ice_categories', twod=.true., ckp=.false.))
    call processor%apply(make_spec('sea_ice_pensolar_frac_diffuse', main%surface, W3, &
        mult='sea_ice_categories', twod=.true., ckp=.false.))

    ! Sea surface velocity vector components provided via the coupler.
    call processor%apply(make_spec('sea_u_current', main%surface, W3, twod=.true.,    &
        ckp=.false.))
    call processor%apply(make_spec('sea_v_current', main%surface, W3, twod=.true.,    &
        ckp=.false.))
    call processor%apply(make_spec('sea_current_w2', main%surface, W2, ckp=.false.))

    ! Sea surface temperature perturbation (set in IAU, applied in coupling)
    ! Needs checkpointing if coupled and iau_sst is enabled
    if ((checkpoint_couple) .and. (iau_sst)) then
      sst_pert_flag = .true.
    else
      sst_pert_flag = .false.
    end if
    call processor%apply(make_spec('sea_surf_temp_pert', main%surface, W3,      &
        twod=.true., ckp=sst_pert_flag, empty = (.not. sst_pert_flag)))

    ! Fields on surface tiles, don't need checkpointing
    call processor%apply(make_spec('tile_heat_flux', main%surface, W3,          &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('alpha1_tile', main%surface, W3,             &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('ashtf_prime_tile', main%surface, W3,        &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('dtstar_tile', main%surface, W3,             &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('fracaero_t_tile', main%surface, W3,               &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('fracaero_s_tile', main%surface, W3,               &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('z0h_tile', main%surface, W3,                &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('z0m_tile', main%surface, W3,                &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('rhokh_tile', main%surface, W3,              &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('chr1p5m_tile', main%surface, W3,            &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('resfs_tile', main%surface, W3,              &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('canhc_tile', main%surface, W3,              &
        mult='surface_tiles', twod=.true.))
    call processor%apply(make_spec('gc_tile', main%surface, W3, mult='surface_tiles', &
        twod=.true.))

    ! Fields on surface tiles used by coupler, need checkpointing in coupled models
    call processor%apply(make_spec('tile_moisture_flux', main%surface, W3,      &
        mult='surface_tiles', twod=.true., ckp=checkpoint_couple))
    call processor%apply(make_spec('snowice_melt', main%surface, W3,            &
        mult='surface_tiles', twod=.true., ckp=checkpoint_couple))
    call processor%apply(make_spec('surf_ht_flux', main%surface, W3,            &
        mult='surface_tiles', twod=.true., ckp=checkpoint_couple))
    call processor%apply(make_spec('snowice_sublimation', main%surface, W3,     &
        mult='surface_tiles', twod=.true., ckp=checkpoint_couple))
    call processor%apply(make_spec('sea_ice_pensolar', main%surface, W3,        &
        mult='sea_ice_categories', twod=.true., ckp=checkpoint_couple))

    ! 2D fields, don't need checkpointing
    call processor%apply(make_spec('ocn_cpl_point', main%surface, W3,           &
        twod=.true., is_int=.true.))
    call processor%apply(make_spec('z0m_eff', main%surface, W3, twod=.true.))
    call processor%apply(make_spec('net_prim_prod', main%surface, W3, twod=.true.))
    call processor%apply(make_spec('taux_ssi', main%surface, W3, twod=.true.))
    call processor%apply(make_spec('tauy_ssi', main%surface, W3, twod=.true.))
    call processor%apply(make_spec('z0m', main%surface, W3, twod=.true.))
    call processor%apply(make_spec('ustar', main%surface, W3, twod=.true.))
    call processor%apply(make_spec('wspd10m', main%surface, W3, twod=.true.))
    call processor%apply(make_spec('urbdisp', main%surface, W3, twod=.true., &
                                   empty = (.not. l_urban2t) ))
    call processor%apply(make_spec('urbztm', main%surface, W3, twod=.true., &
                                   empty = (.not. l_urban2t) ))
    call processor%apply(make_spec('urbalbwl', main%surface, W3, twod=.true., &
                                   empty = (.not. l_urban2t) ))
    call processor%apply(make_spec('urbalbrd', main%surface, W3, twod=.true., &
                                   empty = (.not. l_urban2t) ))
    call processor%apply(make_spec('urbemisw', main%surface, W3, twod=.true., &
                                   empty = (.not. l_urban2t) ))
    call processor%apply(make_spec('urbemisr', main%surface, W3, twod=.true., &
                                   empty = (.not. l_urban2t) ))
    call processor%apply(make_spec('urbemisc', main%surface, W3, twod=.true., &
                                   empty = (.not. l_urban2t)))
    ! 2D fields, need checkpointing for urban-2-tile schemes
    call processor%apply(make_spec('urbwrr', main%surface, twod=.true., &
                                   ckp=l_urban2t, empty = (.not. l_urban2t) ))
    call processor%apply(make_spec('urbhwr', main%surface, twod=.true., &
                                   ckp=l_urban2t, empty = (.not. l_urban2t) ))
    call processor%apply(make_spec('urbhgt', main%surface, twod=.true., &
                                   ckp=l_urban2t, empty = (.not. l_urban2t) ))
    ! 2D field, need checkpointing if the internal flux scheme is used
    call processor%apply(make_spec('internal_flux', main%surface, &
        W3, twod=.true., &
        ckp=(surf_temp_forcing == surf_temp_forcing_int_flux), &
        empty = (surf_temp_forcing /= surf_temp_forcing_int_flux) ))

    ! Space for variables required for regridding to cell faces
    ! vector_space => function_space_collection%get_fs(twod_mesh, 0, 0, W2,
    !     get_ndata_val('surface_regrid_vars'))
    call processor%apply(make_spec('surf_interp_w2', main%surface, W2,          &
        mult='surface_regrid_vars', twod=.true.))
    call processor%apply(make_spec('surf_interp', main%surface, W3,            &
        mult='surface_regrid_vars', twod=.true.))

    ! 2D fields at W2 points
    ! vector_space => function_space_collection%get_fs(twod_mesh, 0, 0,
    !     W2)
    call processor%apply(make_spec('tau_land_w2', main%surface, W2, twod=.true.))
    call processor%apply(make_spec('tau_ssi_w2', main%surface, W2, twod=.true.))

    ! Field on soil levels and land tiles
    ! vector_space => function_space_collection%get_fs(twod_mesh, 0, 0, W3,
    !     get_ndata_val('soil_levels_and_tiles'))
    call processor%apply(make_spec('tile_water_extract', main%surface, W3,      &
        mult='soil_levels_and_tiles', twod=.true.))

    !========================================================================
    ! Fields owned by the soil hydrology scheme
    !========================================================================

    ! 2D fields, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('soil_albedo', main%soil,                    &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_roughness', main%soil,                 &
        ckp=(checkpoint_flag .and. l_vary_z0m_soil)))
    call processor%apply(make_spec('soil_thermal_cond', main%soil,              &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('mean_topog_index', main%soil,               &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('a_sat_frac', main%soil,                     &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('c_sat_frac', main%soil,                     &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('a_wet_frac', main%soil,                     &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('c_wet_frac', main%soil,                     &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_sat_frac', main%soil,                  &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('water_table', main%soil,                    &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('wetness_under_soil', main%soil,             &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_moist_wilt', main%soil,                &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_moist_crit', main%soil,                &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_moist_sat', main%soil,                 &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_cond_sat', main%soil,                  &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_thermal_cap', main%soil,               &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_suction_sat', main%soil,               &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('clapp_horn_b', main%soil,                   &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('surface_runoff', main%soil, W3, twod=.true.))
    call processor%apply(make_spec('sub_surface_runoff', main%soil, W3,         &
        twod=.true.))


    ! Fields on soil levels
    call processor%apply(make_spec('soil_temperature', main%soil,               &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_moisture', main%soil,                  &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('unfrozen_soil_moisture', main%soil,         &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('frozen_soil_moisture', main%soil,           &
          ckp=checkpoint_flag))

    ! 2D fields, don't need checkpointing
    call processor%apply(make_spec('soil_moist_avail', main%soil, W3, twod=.true.))
    call processor%apply(make_spec('thermal_cond_wet_soil', main%soil, W3,      &
        twod=.true.))

    !========================================================================
    ! Fields owned by the snow scheme
    !========================================================================

    ! Fields on surface tiles, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('tile_snow_mass', main%snow,                 &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('tile_snow_rgrain', main%snow,               &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('n_snow_layers', main%snow,                  &
          ckp=checkpoint_flag, is_int=.true.))
    call processor%apply(make_spec('snow_depth', main%snow,                     &
         ckp=checkpoint_flag))
    call processor%apply(make_spec('snowpack_density', main%snow,               &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('snow_under_canopy', main%snow,              &
          ckp=checkpoint_flag))

    ! Fields on snow layers
    call processor%apply(make_spec('snow_layer_thickness', main%snow,           &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('snow_layer_ice_mass', main%snow,            &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('snow_layer_liq_mass', main%snow,            &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('snow_layer_temp', main%snow,                &
          ckp=checkpoint_flag))
    call processor%apply(make_spec('snow_layer_rgrain', main%snow,              &
          ckp=checkpoint_flag))

    ! 2D fields
    call processor%apply(make_spec('snow_soot', main%snow,                      &
        ckp=checkpoint_flag))

    ! Fields which don't need checkpointing
    call processor%apply(make_spec('snow_unload_rate', main%snow, W3,           &
        mult='plant_func_types', twod=.true.))

    !========================================================================
    ! Fields owned by the chemistry scheme
    !========================================================================

    ! Chemistry emissions - only accessed if the UKCA (algorithm &) kernel is
    ! called. This occurs when 'glomap_mode' aerosol scheme is chosen,
    ! irrespective of the chemistry scheme since currently there is no option
    ! to use chemistry scheme without aerosols.
    ! The emissions are never advected but checkpointed for chem_scheme_strattrop
    if ( aerosol == aerosol_um .and.            &
         ( glomap_mode == glomap_mode_ukca .or. &
           glomap_mode == glomap_mode_dust_and_clim ) ) then

      if ( chem_scheme == chem_scheme_strattrop .or.                           &
         chem_scheme == chem_scheme_strat_test ) then
        checkpoint_flag = .true.
        is_empty = .false.
      else
        checkpoint_flag = .false.
        is_empty = .true.
      end if

      call processor%apply(make_spec('emiss_c2h6', main%chemistry,             &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_c3h8', main%chemistry,             &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_c5h8', main%chemistry,             &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_ch4', main%chemistry,              &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_co', main%chemistry,               &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_hcho', main%chemistry,             &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_me2co', main%chemistry,            &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_mecho', main%chemistry,            &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_nh3', main%chemistry,              &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_no', main%chemistry,               &
          ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('emiss_meoh', main%chemistry,             &
          ckp=checkpoint_flag, empty=is_empty))
      ! 3-D emissions
      call processor%apply(make_spec('emiss_no_aircrft', main%chemistry,       &
          ckp=checkpoint_flag, empty=is_empty))

    endif  ! glomap = ukca

    ! Chemistry tracers and prognostic fields - need to be always active since
    ! they are accessed in ukca and conv_gr kernels, but checkpointed or
    ! advected only when chemistry is active and created with empty data when
    ! not used.
    if ( chem_scheme == chem_scheme_strattrop .or.                             &
         chem_scheme == chem_scheme_strat_test ) then
      advection_flag  = .true.
      checkpoint_flag = .true.
      is_empty        = .false.
    else
      advection_flag  = .false.
      checkpoint_flag = .false.
      is_empty        = .true.
    end if

    call processor%apply(make_spec('o3p', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('n', main%chemistry, empty=is_empty,        &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('no', main%chemistry, empty=is_empty,       &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('lumped_n', main%chemistry, empty=is_empty, &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('n2o5', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('ho2no2', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('hono2', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    is_rad = ( (ch4_rad_opt == ch4_rad_opt_ancil) .or. &
               (ch4_rad_opt == ch4_rad_opt_prognostic) )
    call processor%apply(make_spec('ch4', main%chemistry, &
      empty=(is_empty .and. .not. is_rad), &
      adv_coll=if_adv((advection_flag .or. &
      (ch4_rad_opt == ch4_rad_opt_prognostic)), adv%last_con), &
      ckp=(checkpoint_flag .or. is_rad)))
    is_rad = ( (co_rad_opt == co_rad_opt_ancil) .or. &
               (co_rad_opt == co_rad_opt_prognostic) )
    call processor%apply(make_spec('co', main%chemistry, &
      empty=(is_empty .and. .not. is_rad), &
      adv_coll=if_adv((advection_flag .or. &
      (co_rad_opt == co_rad_opt_prognostic)), adv%last_con), &
      ckp=(checkpoint_flag .or. is_rad)))
    call processor%apply(make_spec('hcho', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('meooh', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('h', main%chemistry, empty=is_empty,        &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('cl', main%chemistry, empty=is_empty,       &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('cl2o2', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('clo', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('oclo', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('br', main%chemistry, empty=is_empty,       &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('lumped_br', main%chemistry, empty=is_empty,&
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('brcl', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('brono2', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    is_rad = ( (n2o_rad_opt == n2o_rad_opt_ancil) .or. &
               (n2o_rad_opt == n2o_rad_opt_prognostic) )
    call processor%apply(make_spec('n2o', main%chemistry, &
      empty=(is_empty .and. .not. is_rad), &
      adv_coll=if_adv((advection_flag .or. &
      (n2o_rad_opt == n2o_rad_opt_prognostic)), adv%last_con), &
      ckp=(checkpoint_flag .or. is_rad)))
    call processor%apply(make_spec('lumped_cl', main%chemistry, empty=is_empty,&
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('hocl', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('hbr', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('hobr', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('clono2', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('cfcl3', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('cf2cl2', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('mebr', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('hono', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('c2h6', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('etooh', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('mecho', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('pan', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('c3h8', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('n_prooh', main%chemistry, empty=is_empty,  &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('i_prooh', main%chemistry, empty=is_empty,  &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('etcho', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('me2co', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('mecoch2ooh', main%chemistry,empty=is_empty,&
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('ppan', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('meono2', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('c5h8', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('isooh', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('ison', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('macr', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('macrooh', main%chemistry, empty=is_empty,  &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('mpan', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('hacet', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('mgly', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('nald', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('hcooh', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('meco3h', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('meco2h', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    is_rad = ( (h2_rad_opt == h2_rad_opt_ancil) .or. &
               (h2_rad_opt == h2_rad_opt_prognostic) )
    call processor%apply(make_spec('h2', main%chemistry, &
      empty=(is_empty .and. .not. is_rad), &
      adv_coll=if_adv((advection_flag .or. &
      (h2_rad_opt == h2_rad_opt_prognostic)), adv%last_con), &
      ckp=(checkpoint_flag .or. is_rad)))
    call processor%apply(make_spec('meoh', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('msa', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    is_rad = ( (nh3_rad_opt == nh3_rad_opt_ancil) .or. &
               (nh3_rad_opt == nh3_rad_opt_prognostic) )
    call processor%apply(make_spec('nh3', main%chemistry, &
      empty=(is_empty .and. .not. is_rad), &
      adv_coll=if_adv((advection_flag .or. &
      (nh3_rad_opt == nh3_rad_opt_prognostic)), adv%last_con), &
      ckp=(checkpoint_flag .or. is_rad)))
    call processor%apply(make_spec('cs2', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('csul', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('h2s', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('so3', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('passive_o3', main%chemistry,               &
      empty=is_empty,                                                          &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('age_of_air', main%chemistry,               &
      empty=is_empty,                                                          &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Non-UKCA gases that can be radiatively active
    is_rad = ( (co2_rad_opt == co2_rad_opt_ancil) .or. &
               (co2_rad_opt == co2_rad_opt_prognostic) )
    call processor%apply(make_spec('co2', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((co2_rad_opt == co2_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = (h2o_rad_opt == h2o_rad_opt_ancil) .or. &
             ! If dry, this field is used instead of the standard mr field
             ( (h2o_rad_opt == h2o_rad_opt_prognostic) .and. &
               (moisture_formulation == moisture_formulation_dry) )
    call processor%apply(make_spec('h2o', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((h2o_rad_opt == h2o_rad_opt_prognostic) .and. &
                      (moisture_formulation == moisture_formulation_dry), &
                      adv%last_con), &
      ckp=is_rad))
    is_rad = ( (hcn_rad_opt == hcn_rad_opt_ancil) .or. &
               (hcn_rad_opt == hcn_rad_opt_prognostic) )
    call processor%apply(make_spec('hcn', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((hcn_rad_opt == hcn_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (he_rad_opt == he_rad_opt_ancil) .or. &
               (he_rad_opt == he_rad_opt_prognostic) )
    call processor%apply(make_spec('he', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((he_rad_opt == he_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (n2_rad_opt == n2_rad_opt_ancil) .or. &
               (n2_rad_opt == n2_rad_opt_prognostic) )
    call processor%apply(make_spec('n2', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((n2_rad_opt == n2_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (o2_rad_opt == o2_rad_opt_ancil) .or. &
               (o2_rad_opt == o2_rad_opt_prognostic) )
    call processor%apply(make_spec('o2', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((o2_rad_opt == o2_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (cs_rad_opt == cs_rad_opt_ancil) .or. &
               (cs_rad_opt == cs_rad_opt_prognostic) )
    call processor%apply(make_spec('cs', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((cs_rad_opt == cs_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (k_rad_opt == k_rad_opt_ancil) .or. &
               (k_rad_opt == k_rad_opt_prognostic) )
    call processor%apply(make_spec('k', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((k_rad_opt == k_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (li_rad_opt == li_rad_opt_ancil) .or. &
               (li_rad_opt == li_rad_opt_prognostic) )
    call processor%apply(make_spec('li', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((li_rad_opt == li_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (na_rad_opt == na_rad_opt_ancil) .or. &
               (na_rad_opt == na_rad_opt_prognostic) )
    call processor%apply(make_spec('na', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((na_rad_opt == na_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (rb_rad_opt == rb_rad_opt_ancil) .or. &
               (rb_rad_opt == rb_rad_opt_prognostic) )
    call processor%apply(make_spec('rb', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((rb_rad_opt == rb_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (tio_rad_opt == tio_rad_opt_ancil) .or. &
               (tio_rad_opt == tio_rad_opt_prognostic) )
    call processor%apply(make_spec('tio', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((tio_rad_opt == tio_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))
    is_rad = ( (vo_rad_opt == vo_rad_opt_ancil) .or. &
               (vo_rad_opt == vo_rad_opt_prognostic) )
    call processor%apply(make_spec('vo', main%chemistry, &
      empty=(.not. is_rad), &
      adv_coll=if_adv((vo_rad_opt == vo_rad_opt_prognostic), adv%last_con), &
      ckp=is_rad))


    ! RO2 class of species, not advected if l_ukca_ro2_ntp = true
    if ( chem_scheme == chem_scheme_strattrop  .or.  &
         chem_scheme == chem_scheme_strat_test ) then
      checkpoint_flag = .true.
      is_empty = .false.
      if ( l_ukca_ro2_ntp ) then
        advection_flag = .false.
      else
        advection_flag = .true.
      end if
    else
      checkpoint_flag = .false.
      advection_flag  = .false.
      is_empty        = .true.
    end if
    call processor%apply(make_spec('meoo', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('etoo', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('meco3', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('n_proo', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('i_proo', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('etco3', main%chemistry, empty=is_empty,    &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('mecoch2oo', main%chemistry, empty=is_empty,&
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))

    ! Fields that are never advected ('lumped' versions of no2,bro,hcl are)
    advection_flag = .false.
    call processor%apply(make_spec('o1d', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('no2', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('bro', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('hcl', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('iso2', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('macro2', main%chemistry, empty=is_empty,   &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))

    ! Species which are active in all chemistry schemes, but some not advected
    ! for Offline oxidants
    advection_flag = .false.
    checkpoint_flag = .false.
    is_empty = .true.
    if ( chem_scheme == chem_scheme_offline_ox .or.       &
         chem_scheme == chem_scheme_strattrop  .or.       &
         chem_scheme == chem_scheme_strat_test ) then
      ! Don't need advecting or checkpointing for dust only
      checkpoint_flag = glomap_mode == glomap_mode_ukca
      advection_flag  = glomap_mode == glomap_mode_ukca
      is_empty        = .false.
    end if
    ! H2O2 - advected under all schemes
    call processor%apply(make_spec('h2o2', main%chemistry, empty=is_empty,     &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))

    if ( chem_scheme == chem_scheme_offline_ox ) then
      advection_flag = .false.
    end if
    is_rad = (o3_rad_opt == o3_rad_opt_prognostic)
    ! Special case: ozone from radiation_fields is used instead
    call processor%apply(make_spec('o3', main%chemistry, &
      empty=(is_empty .and. .not. is_rad), &
      adv_coll=if_adv((advection_flag .or. is_rad), adv%last_con), &
      ckp=(checkpoint_flag .or. is_rad)))
    call processor%apply(make_spec('no3', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('oh', main%chemistry, empty=is_empty,       &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('ho2', main%chemistry, empty=is_empty,      &
      adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))

    ! Aerosol precursors - always active but checkpointed and advected only for
    ! glomap_mode_ukca aerosol scheme
    if ( aerosol == aerosol_um .and.           &
         glomap_mode == glomap_mode_ukca ) then
      ! Don't need advecting or checkpointing for dust only
      checkpoint_flag = .true.
      advection_flag = .true.
      is_empty       = .false.
    else
      checkpoint_flag = .false.
      advection_flag = .false.
      is_empty       = .true.
    end if
    call processor%apply(make_spec('dms', main%chemistry, empty=is_empty,      &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    is_rad = ( (so2_rad_opt == so2_rad_opt_ancil) .or. &
               (so2_rad_opt == so2_rad_opt_prognostic) )
    call processor%apply(make_spec('so2', main%chemistry, &
      empty=(is_empty .and. .not. is_rad), &
      adv_coll=if_adv((advection_flag .or. &
      (so2_rad_opt == so2_rad_opt_prognostic)), adv%last_con), &
      ckp=(checkpoint_flag .or. is_rad)))
    call processor%apply(make_spec('h2so4', main%chemistry, empty=is_empty,    &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('dmso', main%chemistry, empty=is_empty,     &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('monoterpene', main%chemistry,              &
        empty=is_empty,                                                        &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    call processor%apply(make_spec('secondary_organic', main%chemistry,        &
        empty=is_empty,                                                        &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Upper limit for H2O2 (ancillary field) only active for glomap_mode and
    !  checkpointed for offline oxidants, never advected
    call processor%apply(make_spec('h2o2_limit', main%chemistry,               &
         empty=is_empty, ckp=(chem_scheme == chem_scheme_offline_ox) ) )

    !========================================================================
    ! Fields owned by the aerosol scheme
    !========================================================================

    ! Flags for some 2D and 3D fields
    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then
      checkpoint_flag = .true.
      checkpoint_GC3 = (emissions == emissions_GC3)
      checkpoint_GC5 = (emissions == emissions_GC5)
      is_empty = .false.
    else
      checkpoint_flag = .false.
      checkpoint_GC3 = .false.
      checkpoint_GC5 = .false.
      is_empty = .true.
    end if
    ! 2D fields, might need checkpointing
    if ( aerosol == aerosol_um .and.            &
         ( glomap_mode == glomap_mode_ukca .or. &
           glomap_mode == glomap_mode_dust_and_clim ) ) then
      call processor%apply(make_spec('dms_conc_ocean', main%aerosol,           &
           empty=is_empty, ckp=checkpoint_flag))
      call processor%apply(make_spec('emiss_bc_biofuel', main%aerosol,         &
           empty=is_empty, ckp=checkpoint_flag))
      call processor%apply(make_spec('emiss_bc_fossil', main%aerosol,          &
           empty=is_empty, ckp=checkpoint_flag))
      call processor%apply(make_spec('emiss_dms_land', main%aerosol,           &
           empty=is_empty, ckp=checkpoint_flag))
      call processor%apply(make_spec('emiss_monoterp', main%aerosol,           &
           empty=is_empty, ckp=checkpoint_flag))
      call processor%apply(make_spec('emiss_om_biofuel', main%aerosol,         &
           empty=is_empty, ckp=checkpoint_flag))
      call processor%apply(make_spec('emiss_om_fossil', main%aerosol,          &
           empty=is_empty, ckp=checkpoint_flag))
      call processor%apply(make_spec('emiss_so2_low', main%aerosol,            &
           empty=is_empty, ckp=checkpoint_flag))
      call processor%apply(make_spec('emiss_so2_high', main%aerosol,           &
           empty=is_empty, ckp=checkpoint_GC3))
      call processor%apply(make_spec('emiss_bc_biomass_high', main%aerosol,    &
           empty=is_empty, ckp=checkpoint_GC5))
      call processor%apply(make_spec('emiss_bc_biomass_low', main%aerosol,     &
           empty=is_empty, ckp=checkpoint_GC5))
      call processor%apply(make_spec('emiss_om_biomass_high', main%aerosol,    &
           empty=is_empty, ckp=checkpoint_GC5))
      call processor%apply(make_spec('emiss_om_biomass_low', main%aerosol,     &
           empty=is_empty, ckp=checkpoint_GC5))
      call processor%apply(make_spec('surf_wetness', main%aerosol,             &
           empty=is_empty, ckp=checkpoint_flag))
    end if

    ! 3D fields, might need checkpointing - flags set above
    if ( aerosol == aerosol_um .and.                                           &
         ( glomap_mode == glomap_mode_ukca .or.                                &
           glomap_mode == glomap_mode_dust_and_clim ) ) then
      call processor%apply(make_spec('emiss_bc_biomass', main%aerosol,         &
           empty=is_empty, ckp=checkpoint_GC3))
      call processor%apply(make_spec('emiss_om_biomass', main%aerosol,         &
           empty=is_empty, ckp=checkpoint_GC3))
      call processor%apply(make_spec('emiss_so2_nat', main%aerosol,            &
           empty=is_empty, ckp=checkpoint_flag))
    end if

    ! Flags for more 2D fields
    if (aerosol == aerosol_um .and. (glomap_mode == glomap_mode_ukca .or.             &
                                     glomap_mode == glomap_mode_dust_and_clim )) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if

    call processor%apply(make_spec('soil_clay', main%aerosol,                  &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('soil_sand', main%aerosol,                  &
        ckp=checkpoint_flag))

    ! 3D fields, might need checkpointing and/or advecting
    ! Nucleation mode is only used with UKCA
    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then
      checkpoint_flag = .true.
      advection_flag = .true.
      is_empty = .false.
    else
      checkpoint_flag = .false.
      advection_flag = .false.
      is_empty = .true.
    end if
    ! Nucleation soluble mode number mixing ratio
    call processor%apply(make_spec('n_nuc_sol', main%aerosol, empty=is_empty,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Nucleation soluble H2SO4 aerosol mmr
    call processor%apply(make_spec('nuc_sol_su', main%aerosol, empty=is_empty, &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Nucleation soluble organic carbon aerosol mmr
    call processor%apply(make_spec('nuc_sol_om', main%aerosol, empty=is_empty, &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))

    ! Set flag defaults
    checkpoint_flag     = .false.
    advection_flag      = .false.
    advection_flag_dust = .false.

    ! 3D fields, might need checkpointing and/or advecting
    if ( aerosol == aerosol_um .and.                                           &
         ( glomap_mode == glomap_mode_climatology .or.                         &
           glomap_mode == glomap_mode_ukca        .or.                         &
           glomap_mode == glomap_mode_dust_and_clim ) ) then
      checkpoint_flag = .true.
    end if

    if ( aerosol == aerosol_um  ) then

      select case ( glomap_mode )
      case( glomap_mode_ukca )
        advection_flag = .true.
        advection_flag_dust = .true.

      case( glomap_mode_dust_and_clim )
        advection_flag_dust = .true.

      end select

    end if

    ! Aitken soluble mode number mixing ratio
    call processor%apply(make_spec('n_ait_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Aitken soluble H2SO4 aerosol mmr
    call processor%apply(make_spec('ait_sol_su', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Aitken soluble black carbon aerosol mmr
    call processor%apply(make_spec('ait_sol_bc', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Aitken soluble organic carbon aerosol mmr
    call processor%apply(make_spec('ait_sol_om', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Accumulation soluble mode number mixing ratio
    call processor%apply(make_spec('n_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Accumulation soluble H2SO4 aerosol mmr
    call processor%apply(make_spec('acc_sol_su', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Accumulation soluble black carbon aerosol mmr
    call processor%apply(make_spec('acc_sol_bc', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Accumulation soluble organic carbon aerosol mmr
    call processor%apply(make_spec('acc_sol_om', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Accumulation soluble sea salt aerosol mmr
    call processor%apply(make_spec('acc_sol_ss', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Coarse soluble mode number mixing ratio
    call processor%apply(make_spec('n_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol, &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Coarse soluble H2SO4 aerosol mmr
    call processor%apply(make_spec('cor_sol_su', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Coarse soluble black carbon aerosol mmr
    call processor%apply(make_spec('cor_sol_bc', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Coarse soluble organic carbon aerosol mmr
    call processor%apply(make_spec('cor_sol_om', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Coarse soluble sea salt aerosol mmr
    call processor%apply(make_spec('cor_sol_ss', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Aitken insoluble mode number mixing ratio
    call processor%apply(make_spec('n_ait_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Aitken insoluble black carbon aerosol mmr
    call processor%apply(make_spec('ait_ins_bc', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Aitken insoluble organic carbon aerosol mmr
    call processor%apply(make_spec('ait_ins_om', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,  &
        adv_coll=if_adv(advection_flag, adv%last_con), ckp=checkpoint_flag))
    ! Accumulation insoluble mode number mixing ratio
    call processor%apply(make_spec('n_acc_ins', main%aerosol, Wtheta, coarse=.false.,   &
        adv_coll=if_adv(advection_flag_dust, adv%last_con), ckp=checkpoint_flag))
    ! Accumulation insoluble dust aerosol mmr
    call processor%apply(make_spec('acc_ins_du', main%aerosol, Wtheta, coarse=.false.,  &
        adv_coll=if_adv(advection_flag_dust, adv%last_con), ckp=checkpoint_flag))
    ! Coarse insoluble mode number mixing ratio
    call processor%apply(make_spec('n_cor_ins', main%aerosol, Wtheta, coarse=.false.,   &
        adv_coll=if_adv(advection_flag_dust, adv%last_con), ckp=checkpoint_flag))
    ! Coarse insoluble dust aerosol mmr
    call processor%apply(make_spec('cor_ins_du', main%aerosol, Wtheta, coarse=.false.,  &
        adv_coll=if_adv(advection_flag_dust, adv%last_con), ckp=checkpoint_flag))

    ! 3D fields, might need checkpointing
    if (aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if

    ! Cloud droplet number concentration
    call processor%apply(make_spec('cloud_drop_no_conc', main%aerosol,          &
        ckp=checkpoint_flag))

    if (aerosol == aerosol_um) then
        ! Dry diameter Aitken mode (Soluble)
        call processor%apply(make_spec('drydp_ait_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Dry diameter Accumulation mode (Soluble)
        call processor%apply(make_spec('drydp_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Dry diameter Coarse mode (Soluble)
        call processor%apply(make_spec('drydp_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Dry diameter Aitken mode (Insoluble)
        call processor%apply(make_spec('drydp_ait_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Dry diameter Accumulation mode (Insoluble)
        call processor%apply(make_spec('drydp_acc_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Dry diameter Coarse mode (Insoluble)
        call processor%apply(make_spec('drydp_cor_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Wet diameter Aitken mode (Soluble)
        call processor%apply(make_spec('wetdp_ait_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Wet diameter Accumulation mode (Soluble)
        call processor%apply(make_spec('wetdp_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Wet diameter Coarse mode (Soluble)
        call processor%apply(make_spec('wetdp_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Particle density Aitken mode (Soluble)
        call processor%apply(make_spec('rhopar_ait_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Particle density Accumulation mode (Soluble)
        call processor%apply(make_spec('rhopar_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Particle density Coarse mode (Soluble)
        call processor%apply(make_spec('rhopar_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Particle density Aitken mode (Insoluble)
        call processor%apply(make_spec('rhopar_ait_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Particle density Accumulation mode (Insoluble)
        call processor%apply(make_spec('rhopar_acc_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Particle density Coarse mode (Insoluble)
        call processor%apply(make_spec('rhopar_cor_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume of water Aitken mode (Soluble)
        call processor%apply(make_spec('pvol_wat_ait_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume of water Accumulation mode (Soluble)
        call processor%apply(make_spec('pvol_wat_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume of water Coarse mode (Soluble)
        call processor%apply(make_spec('pvol_wat_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Sulphate Aitken mode (Soluble)
        call processor%apply(make_spec('pvol_su_ait_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Black Carbon Aitken mode (Soluble)
        call processor%apply(make_spec('pvol_bc_ait_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Organic Matter Aitken mode (Soluble)
        call processor%apply(make_spec('pvol_om_ait_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Sulphate Accumulation mode (Soluble)
        call processor%apply(make_spec('pvol_su_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Black Carbon Accumulation mode (Soluble)
        call processor%apply(make_spec('pvol_bc_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Organic Matter Accumulation mode (Soluble)
        call processor%apply(make_spec('pvol_om_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Sea Salt Accumulation mode (Soluble)
        call processor%apply(make_spec('pvol_ss_acc_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Sulphate Coarse mode (Soluble)
        call processor%apply(make_spec('pvol_su_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Black Carbon Coarse mode (Soluble)
        call processor%apply(make_spec('pvol_bc_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Organic Matter Coarse mode (Soluble)
        call processor%apply(make_spec('pvol_om_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Sea Salt Coarse mode (Soluble)
        call processor%apply(make_spec('pvol_ss_cor_sol', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Black Carbon Aitken mode (Insoluble)
        call processor%apply(make_spec('pvol_bc_ait_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Organic Matter Aitken mode (Insoluble)
        call processor%apply(make_spec('pvol_om_ait_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Dust Accumulation mode (Insoluble)
        call processor%apply(make_spec('pvol_du_acc_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
        ! Partial volume component Dust Coarse mode (Insoluble)
        call processor%apply(make_spec('pvol_du_cor_ins', main%aerosol, Wtheta, coarse=coarse_rad_aerosol,   &
            ckp=checkpoint_flag))
    end if

    ! Fields on dust space, might need checkpointing
    if (aerosol == aerosol_um .and. (glomap_mode == glomap_mode_ukca .or.             &
                                     glomap_mode == glomap_mode_dust_and_clim )) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if

    call processor%apply(make_spec('dust_mrel', main%aerosol,                   &
          ckp=checkpoint_flag))

    ! Fields on dust space, don't need checkpointing
    call processor%apply(make_spec('dust_flux', main%aerosol, W3,               &
        mult='dust_divisions', twod=.true.))

    ! 3D fields, don't need checkpointing
    ! Sulphuric Acid aerosol MMR
    call processor%apply(make_spec('sulphuric', main%aerosol, Wtheta))

    ! Murk field
    call processor%apply(make_spec('murk', main%aerosol, &
                   adv_coll=if_adv(murk_prognostic, adv%last_con), &
                   ckp=murk_prognostic, empty = (.not. murk) ))
    call processor%apply(make_spec('murk_source', main%aerosol, Wtheta, &
                         empty = (.not. murk_prognostic) ))

    ! EasyAerosol fields, might need checkpointing
    if ( easyaerosol_sw ) then
       checkpoint_flag = .true.

      call processor%apply(make_spec('easy_absorption_sw', main%aerosol,        &
         ckp=checkpoint_flag))

      call processor%apply(make_spec('easy_extinction_sw', main%aerosol,        &
         ckp=checkpoint_flag))

      call processor%apply(make_spec('easy_asymmetry_sw', main%aerosol,         &
         ckp=checkpoint_flag))

    end if

    if ( easyaerosol_lw ) then
       checkpoint_flag = .true.

      call processor%apply(make_spec('easy_absorption_lw', main%aerosol,        &
         ckp=checkpoint_flag))

      call processor%apply(make_spec('easy_extinction_lw', main%aerosol,        &
         ckp=checkpoint_flag))

      call processor%apply(make_spec('easy_asymmetry_lw', main%aerosol,         &
         ckp=checkpoint_flag))

    end if

    !========================================================================
    ! Fields owned by the stochastic physics scheme
    !========================================================================

    ! 3D fields, don't need checkpointing
    call processor%apply(make_spec('dtheta_stph', main%stph, Wtheta, ckp=.false.))
    call processor%apply(make_spec('dmv_stph', main%stph, Wtheta, ckp=.false.))
    call processor%apply(make_spec('du_stph', main%stph, W2, ckp=.false.))

    ! 2D fields, might need checkpointing
    if ( stochastic_physics == stochastic_physics_um ) then
      if (blpert_type /= blpert_type_off) then
        is_empty = .false.
        checkpoint_flag = .true.
      else
        is_empty = .true.
        checkpoint_flag = .false.
      end if
      call processor%apply(make_spec('blpert_rand_fld', main%stph, W3, &
          twod=.true., ckp=checkpoint_flag, empty=is_empty))
      call processor%apply(make_spec('blpert_flag', main%stph, W3, &
          twod=.true., is_int=.true., ckp=checkpoint_flag, empty=is_empty))
    end if

#endif


  end subroutine process_physics_prognostics



  !>@brief Routine to initialise the field objects required by the physics
  !> @param[in]    mesh                 The current 3d mesh
  !> @param[in]    twod_mesh            The current 2d mesh
  !> @param[in]    coarse_mesh          The coarse 3d mesh
  !> @param[in]    coarse_twod_mesh     The coarse 2d mesh
  !> @param[in]    field_mapper         Provides access to the field collections
  !> @param[in]    clock                Model clock
  subroutine create_physics_prognostics( mesh,                  &
                                         twod_mesh,             &
                                         coarse_mesh,           &
                                         coarse_twod_mesh,      &
                                         field_mapper,          &
                                         clock )
    implicit none

    type( mesh_type ), intent(in), pointer         :: mesh
    type( mesh_type ), intent(in), pointer         :: twod_mesh
    type( mesh_type ), intent(in), pointer         :: coarse_mesh
    type( mesh_type ), intent(in), pointer         :: coarse_twod_mesh
    type( field_mapper_type ), intent(in)          :: field_mapper
    class( clock_type ), intent(in)                :: clock

    type( field_maker_type ) :: creator

    call log_event( 'Create physics prognostics', LOG_LEVEL_INFO )

    call creator%init(mesh, twod_mesh, coarse_mesh, coarse_twod_mesh, field_mapper, clock)

    call field_mapper%sanity_check()

    call process_physics_prognostics(creator)

  end subroutine create_physics_prognostics

end module create_physics_prognostics_mod
