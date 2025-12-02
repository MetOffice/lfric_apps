!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module init_ancils_mod

  use constants_mod,                  only : i_def, l_def, str_def,   &
                                             r_def
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use field_mod,                      only : field_type
  use field_parent_mod,               only : read_interface, &
                                             write_interface
  use io_config_mod,                  only : use_xios_io, checkpoint_read
  use linked_list_mod,                only : linked_list_type
  use lfric_xios_read_mod,            only : read_field_generic
  use lfric_xios_write_mod,           only : write_field_generic
  use field_collection_mod,           only : field_collection_type
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3, WTheta
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use lfric_xios_time_axis_mod,       only : time_axis_type
  use jules_control_init_mod,         only : n_land_tile
  use jules_physics_init_mod,         only : snow_lev_tile
  use jules_surface_types_mod,        only : npft
  use dust_parameters_mod,            only : ndiv
  use initialization_config_mod,      only : ancil_option,          &
                                             ancil_option_idealised,&
                                             ancil_option_updating, &
                                             sst_source,            &
                                             sst_source_start_dump, &
                                             sst_source_surf,       &
                                             init_option,           &
                                             init_option_fd_start_dump, &
                                             snow_source,               &
                                             snow_source_surf,          &
                                             sea_ice_source,        &
                                             sea_ice_source_surf
  use aerosol_config_mod,             only : glomap_mode,               &
                                             glomap_mode_climatology,   &
                                             glomap_mode_dust_and_clim, &
                                             glomap_mode_ukca,          &
                                             emissions, emissions_GC3,  &
                                             emissions_GC5,             &
                                             easyaerosol_cdnc,          &
                                             easyaerosol_sw,            &
                                             easyaerosol_lw,            &
                                             murk_prognostic
  use socrates_init_mod,              only : n_sw_band, n_lw_band
  use jules_surface_config_mod,       only : l_vary_z0m_soil, l_urban2t
  use jules_sea_seaice_config_mod,    only : amip_ice_thick
  use jules_radiation_config_mod,     only : l_sea_alb_var_chl, l_albedo_obs
  use radiation_config_mod,           only : topography, topography_slope, &
                                             topography_horizon, &
                                             n_horiz_ang, n_horiz_layer
  use specified_surface_config_mod,   only : internal_flux_method, &
                                             internal_flux_method_non_uniform, &
                                             surf_temp_forcing, &
                                             surf_temp_forcing_int_flux
  use derived_config_mod,             only : l_couple_sea_ice
  use chemistry_config_mod,           only : chem_scheme,                      &
                                             chem_scheme_strattrop,            &
                                             chem_scheme_strat_test,           &
                                             chem_scheme_offline_ox
  use radiative_gases_config_mod,     only : &
    ch4_rad_opt, ch4_rad_opt_ancil, &
    cs_rad_opt, cs_rad_opt_ancil, &
    co_rad_opt, co_rad_opt_ancil, &
    co2_rad_opt, co2_rad_opt_ancil, &
    h2_rad_opt, h2_rad_opt_ancil, &
    h2o_rad_opt, h2o_rad_opt_ancil, &
    hcn_rad_opt, hcn_rad_opt_ancil, &
    he_rad_opt, he_rad_opt_ancil, &
    k_rad_opt, k_rad_opt_ancil, &
    li_rad_opt, li_rad_opt_ancil, &
    n2_rad_opt, n2_rad_opt_ancil, &
    na_rad_opt, na_rad_opt_ancil, &
    nh3_rad_opt, nh3_rad_opt_ancil, &
    o2_rad_opt, o2_rad_opt_ancil, &
    rb_rad_opt, rb_rad_opt_ancil, &
    so2_rad_opt, so2_rad_opt_ancil, &
    tio_rad_opt, tio_rad_opt_ancil, &
    vo_rad_opt, vo_rad_opt_ancil

  implicit none

  public   :: create_fd_ancils,           &
              create_fd_ancils_idealised, &
              setup_ancil_field

contains

  !> @details Organises fields to be read from ancils into ancil_fields
  !           collection then reads them.
  !> @param[in,out] depository The depository field collection
  !> @param[in,out] ancil_fields Collection for ancillary fields
  !> @param[in] mesh      The current 3d mesh
  !> @param[in] twod_mesh The current 2d mesh
  subroutine create_fd_ancils( depository, ancil_fields, mesh, &
                               twod_mesh, aerosol_mesh, aerosol_twod_mesh, ancil_times_list )

    implicit none

    type( field_collection_type ), intent( inout ) :: depository
    type( field_collection_type ), intent( inout )   :: ancil_fields

    type( mesh_type ), intent(in), pointer :: mesh
    type( mesh_type ), intent(in), pointer :: twod_mesh
    type( mesh_type ), intent(in), pointer :: aerosol_mesh
    type( mesh_type ), intent(in), pointer :: aerosol_twod_mesh

    type(linked_list_type), intent(out) :: ancil_times_list

    ! Time axis objects for different ancil groups - must be saved to be
    ! available after function call
    type(time_axis_type), save :: sea_time_axis
    type(time_axis_type), save :: sst_time_axis
    type(time_axis_type), save :: sea_ice_time_axis
    type(time_axis_type), save :: snow_time_axis
    type(time_axis_type), save :: aerosol_time_axis
    type(time_axis_type), save :: albedo_vis_time_axis
    type(time_axis_type), save :: albedo_nir_time_axis
    type(time_axis_type), save :: pft_time_axis
    type(time_axis_type), save :: ozone_time_axis
    type(time_axis_type), save :: murk_time_axis
    type(time_axis_type), save :: em_bc_bf_time_axis
    type(time_axis_type), save :: em_bc_ff_time_axis
    type(time_axis_type), save :: em_bc_bb_time_axis
    type(time_axis_type), save :: em_bc_bb_hi_time_axis
    type(time_axis_type), save :: em_bc_bb_lo_time_axis
    type(time_axis_type), save :: em_dms_lnd_time_axis
    type(time_axis_type), save :: dms_ocn_time_axis
    type(time_axis_type), save :: em_mterp_time_axis
    type(time_axis_type), save :: em_om_bf_time_axis
    type(time_axis_type), save :: em_om_ff_time_axis
    type(time_axis_type), save :: em_om_bb_time_axis
    type(time_axis_type), save :: em_om_bb_hi_time_axis
    type(time_axis_type), save :: em_om_bb_lo_time_axis
    type(time_axis_type), save :: em_so2_lo_time_axis
    type(time_axis_type), save :: em_so2_hi_time_axis
    type(time_axis_type), save :: em_c2h6_time_axis
    type(time_axis_type), save :: em_c3h8_time_axis
    type(time_axis_type), save :: em_c5h8_time_axis
    type(time_axis_type), save :: em_ch4_time_axis
    type(time_axis_type), save :: em_co_time_axis
    type(time_axis_type), save :: em_hcho_time_axis
    type(time_axis_type), save :: em_me2co_time_axis
    type(time_axis_type), save :: em_mecho_time_axis
    type(time_axis_type), save :: em_nh3_time_axis
    type(time_axis_type), save :: em_no_time_axis
    type(time_axis_type), save :: em_meoh_time_axis
    type(time_axis_type), save :: em_no_aircrft_time_axis
    type(time_axis_type), save :: h2o2_limit_time_axis
    type(time_axis_type), save :: ho2_time_axis
    type(time_axis_type), save :: no3_time_axis
    type(time_axis_type), save :: o3_time_axis
    type(time_axis_type), save :: oh_time_axis
    type(time_axis_type), save :: cloud_drop_no_conc_time_axis
    type(time_axis_type), save :: easy_asymmetry_sw_time_axis
    type(time_axis_type), save :: easy_asymmetry_lw_time_axis
    type(time_axis_type), save :: easy_absorption_sw_time_axis
    type(time_axis_type), save :: easy_absorption_lw_time_axis
    type(time_axis_type), save :: easy_extinction_sw_time_axis
    type(time_axis_type), save :: easy_extinction_lw_time_axis

    ! Time axis options
    logical(l_def),   parameter :: interp_flag=.true.

    ! Set up ancil_fields collection
    write(log_scratch_space,'(A,A)') "Create ancil fields: "// &
          "Setting up ancil field collection"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    ! Here ancil fields are set up with a call to setup_ancil_field. For ancils
    ! that are time-varying, the time-axis is passed to the setup_ancil_field
    ! subroutine.

    !=====  LAND ANCILS  =====
    if (init_option == init_option_fd_start_dump .and. &
         .not. checkpoint_read) then
      call setup_ancil_field("land_area_fraction", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("land_tile_fraction", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.,          &
                              ndata=n_land_tile)
    end if

    call pft_time_axis%initialise("plant_func_time",          &
                                  file_id="plant_func_ancil", &
                                  interp_flag=interp_flag, pop_freq="daily")
    call setup_ancil_field("canopy_height", depository, ancil_fields,         &
                              mesh, twod_mesh, twod=.true., ndata=npft, &
                              time_axis=pft_time_axis)
    call setup_ancil_field("leaf_area_index", depository, ancil_fields,       &
                              mesh, twod_mesh, twod=.true., ndata=npft, &
                              time_axis=pft_time_axis)
    call ancil_times_list%insert_item(pft_time_axis)

    if ( l_urban2t .and. (init_option == init_option_fd_start_dump .and. &
         .not. checkpoint_read) ) then
      call setup_ancil_field("urbwrr", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("urbhwr", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("urbhgt", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
    endif

    !=====  SEA ANCILS  =====
    if ( l_sea_alb_var_chl ) then
      call sea_time_axis%initialise("sea_time", file_id="sea_ancil", &
                                    interp_flag=interp_flag, pop_freq="daily")
      call setup_ancil_field("chloro_sea", depository, ancil_fields, mesh, &
                              twod_mesh, twod=.true.,                      &
                              time_axis=sea_time_axis)
      call ancil_times_list%insert_item(sea_time_axis)
    end if

    if (sst_source /= sst_source_start_dump) then
      if (sst_source == sst_source_surf) then
        call sst_time_axis%initialise("sst_time", file_id="sst_ancil", &
                                      interp_flag=.false., pop_freq="daily", &
                                      window_size=1)

      else !sst_source == 'ancil'
        call sst_time_axis%initialise("sst_time", file_id="sst_ancil", &
                                      interp_flag=interp_flag, pop_freq="daily")
      end if
      call setup_ancil_field("tstar_sea", depository, ancil_fields, mesh, &
                              twod_mesh, twod=.true.,                   &
                              time_axis=sst_time_axis)
      call ancil_times_list%insert_item(sst_time_axis)
    end if

    !=====  SEA ICE ANCILS  =====
    if (.not. l_couple_sea_ice) then
      if (sea_ice_source == sea_ice_source_surf) then
        call sea_ice_time_axis%initialise("sea_ice_time", file_id="sea_ice_ancil", &
                                        interp_flag=.false., pop_freq="daily", &
                                        window_size=1)
      else
        call sea_ice_time_axis%initialise("sea_ice_time", file_id="sea_ice_ancil", &
                                        interp_flag=interp_flag, pop_freq="daily")
      end if
      if (.not. amip_ice_thick) then
        call setup_ancil_field("sea_ice_thickness", depository, ancil_fields, &
                  mesh, twod_mesh, twod=.true., time_axis=sea_ice_time_axis)
      end if
      call setup_ancil_field("sea_ice_fraction", depository, ancil_fields, &
                mesh, twod_mesh, twod=.true., time_axis=sea_ice_time_axis)
      call ancil_times_list%insert_item(sea_ice_time_axis)
    endif

    !=====  SNOW ANCILS ====
    if ( snow_source == snow_source_surf ) then
      call snow_time_axis%initialise("snow_time", file_id="snow_analysis_ancil", &
                                      yearly=.false., interp_flag=.false., &
                                      pop_freq="daily", window_size=1)

      call setup_ancil_field("tile_snow_rgrain_in", depository, ancil_fields,  &
                              mesh, twod_mesh, twod=.true., ndata=n_land_tile, &
                              time_axis=snow_time_axis)
      call setup_ancil_field("tile_snow_mass_in", depository, ancil_fields,    &
                             mesh, twod_mesh, twod=.true., ndata=n_land_tile,  &
                             time_axis=snow_time_axis)
      call setup_ancil_field("snow_under_canopy_in", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true., ndata=n_land_tile, &
                              time_axis=snow_time_axis)
      call setup_ancil_field("snow_depth_in", depository, ancil_fields,        &
                              mesh, twod_mesh, twod=.true., ndata=n_land_tile, &
                              time_axis=snow_time_axis)
      call setup_ancil_field("snowpack_density_in", depository, ancil_fields,  &
                              mesh, twod_mesh, twod=.true., ndata=n_land_tile, &
                              time_axis=snow_time_axis)
      call setup_ancil_field("n_snow_layers_in", depository, ancil_fields,     &
                              mesh, twod_mesh, twod=.true., ndata=n_land_tile, &
                              time_axis=snow_time_axis)
      call setup_ancil_field("snow_layer_thickness", depository, ancil_fields, &
                            mesh, twod_mesh, twod=.true., ndata=snow_lev_tile, &
                            time_axis=snow_time_axis)
      call setup_ancil_field("snow_layer_ice_mass", depository, ancil_fields,  &
                            mesh, twod_mesh, twod=.true., ndata=snow_lev_tile, &
                            time_axis=snow_time_axis)
      call setup_ancil_field("snow_layer_liq_mass", depository, ancil_fields,  &
                            mesh, twod_mesh, twod=.true., ndata=snow_lev_tile, &
                            time_axis=snow_time_axis)
      call setup_ancil_field("snow_layer_temp", depository, ancil_fields,      &
                            mesh, twod_mesh, twod=.true., ndata=snow_lev_tile, &
                            time_axis=snow_time_axis)
      call setup_ancil_field("snow_layer_rgrain", depository, ancil_fields,    &
                            mesh, twod_mesh, twod=.true., ndata=snow_lev_tile, &
                            time_axis=snow_time_axis)
      call ancil_times_list%insert_item(snow_time_axis)
    end if

    !=====  RADIATION ANCILS  =====
    if ( l_albedo_obs ) then
      call albedo_vis_time_axis%initialise("albedo_vis_time",          &
                                           file_id="albedo_vis_ancil", &
                                           interp_flag=interp_flag,    &
                                           pop_freq="daily")
      call setup_ancil_field("albedo_obs_vis", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,         &
                             time_axis=albedo_vis_time_axis)
      call ancil_times_list%insert_item(albedo_vis_time_axis)

      call albedo_nir_time_axis%initialise("albedo_nir_time",          &
                                           file_id="albedo_nir_ancil", &
                                           interp_flag=interp_flag,    &
                                           pop_freq="daily")
      call setup_ancil_field("albedo_obs_nir", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.,        &
                              time_axis=albedo_nir_time_axis)
      call ancil_times_list%insert_item(albedo_nir_time_axis)
    end if

    !=====  SOIL ANCILS  =====
    if (init_option == init_option_fd_start_dump .and. &
         .not. checkpoint_read) then
      call setup_ancil_field("soil_albedo", depository, ancil_fields, mesh, &
                              twod_mesh, twod=.true.)
      if ( l_vary_z0m_soil ) then
        call setup_ancil_field("soil_roughness", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
      endif
      call setup_ancil_field("soil_thermal_cond", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("soil_moist_wilt", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("soil_moist_crit", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("soil_moist_sat", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("soil_cond_sat", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("soil_thermal_cap", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("soil_suction_sat", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("clapp_horn_b", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("mean_topog_index", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("stdev_topog_index", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)

    !=====  OROGRAPHY ANCILS  =====
      call setup_ancil_field("sd_orog", depository, ancil_fields, mesh, &
                              twod_mesh, twod=.true.)
      call setup_ancil_field("grad_xx_orog", depository, ancil_fields, mesh,&
                              twod_mesh, twod=.true.)
      call setup_ancil_field("grad_xy_orog", depository, ancil_fields, mesh,&
                              twod_mesh, twod=.true.)
      call setup_ancil_field("grad_yy_orog", depository, ancil_fields, mesh,&
                              twod_mesh, twod=.true.)
      call setup_ancil_field("peak_to_trough_orog", depository, ancil_fields,  &
                              mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("silhouette_area_orog", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
      if (topography == topography_slope .or. &
           topography == topography_horizon) then
        call setup_ancil_field("grad_x_orog", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
        call setup_ancil_field("grad_y_orog", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
      end if
      if (topography == topography_horizon) then
        call setup_ancil_field("horizon_angle", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true., &
                                ndata=n_horiz_ang*n_horiz_layer)
        call setup_ancil_field("horizon_aspect", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true., &
                                ndata=n_horiz_ang)
      end if
    end if

    !=====  OZONE ANCIL  =====
    call ozone_time_axis%initialise("ozone_time", file_id="ozone_ancil", &
                                    interp_flag=interp_flag, pop_freq="monthly")
    call setup_ancil_field("ozone", depository, ancil_fields, mesh, &
                           twod_mesh, time_axis=ozone_time_axis,    &
                           alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
    call ancil_times_list%insert_item(ozone_time_axis)

    !=====  AEROSOL ANCILS =====
    if (murk_prognostic) then
      call murk_time_axis%initialise("murk_time", file_id="emiss_murk_ancil", &
                                    interp_flag=interp_flag, pop_freq="monthly")
      call setup_ancil_field("murk_source", depository, ancil_fields, mesh, &
                             twod_mesh, time_axis=murk_time_axis)
      call ancil_times_list%insert_item(murk_time_axis)
    end if

    if ( ( glomap_mode == glomap_mode_climatology ) .or. &
         ( glomap_mode == glomap_mode_dust_and_clim ) ) then
      call aerosol_time_axis%initialise( "aerosols_time",          &
                                         file_id="aerosols_ancil", &
                                         interp_flag=interp_flag,  &
                                         pop_freq="daily" )
      call setup_ancil_field("acc_sol_bc", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("acc_sol_om", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("acc_sol_su", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("acc_sol_ss", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("n_acc_sol",  depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("ait_sol_bc", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("ait_sol_om", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("ait_sol_su", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("n_ait_sol",  depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("ait_ins_bc", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("ait_ins_om", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("n_ait_ins",  depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("cor_sol_bc", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("cor_sol_om", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("cor_sol_su", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("cor_sol_ss", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)
      call setup_ancil_field("n_cor_sol",  depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis,        &
                             alt_mesh=aerosol_mesh, alt_twod_mesh=aerosol_twod_mesh)

      ! The following fields will need including when dust is available in the
      ! ancillary file:
      !   acc_sol_du, cor_sol_du, n_acc_ins, acc_ins_du, n_cor_ins, cor_ins_du

      call ancil_times_list%insert_item(aerosol_time_axis)
    end if

    !=====  EMISSION ANCILS (dust only) =====
    if ( ( glomap_mode == glomap_mode_dust_and_clim .or. &
           glomap_mode == glomap_mode_ukca ) .and. &
           init_option == init_option_fd_start_dump .and. &
           .not. checkpoint_read) then

      ! -- Single level ancils
      call setup_ancil_field( "soil_clay", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true. )
      call setup_ancil_field( "soil_sand", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true. )
      call setup_ancil_field( "dust_mrel", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true., ndata=ndiv )

    endif  ! glomap_dust_and_clim

    !=====  EMISSION ANCILS  =====
    if ( glomap_mode == glomap_mode_ukca   .and.                         &
         ancil_option == ancil_option_updating )  then
      ! -- Single level ancils
      call em_bc_bf_time_axis%initialise("em_bc_bf_time",                &
                                       file_id="emiss_bc_biofuel_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_bc_biofuel", depository, ancil_fields,   &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_bc_bf_time_axis)
      call ancil_times_list%insert_item(em_bc_bf_time_axis)

      call em_bc_ff_time_axis%initialise("em_bc_ff_time",                &
                                       file_id="emiss_bc_fossil_ancil",  &
                                       interp_flag=interp_flag,          &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_bc_fossil", depository, ancil_fields,    &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_bc_ff_time_axis)
      call ancil_times_list%insert_item(em_bc_ff_time_axis)

      call em_dms_lnd_time_axis%initialise("em_dms_lnd_time",              &
                                         file_id="emiss_dms_land_ancil",   &
                                         interp_flag=interp_flag,          &
                                         pop_freq="daily")
      call setup_ancil_field("emiss_dms_land", depository, ancil_fields,     &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_dms_lnd_time_axis)
      call ancil_times_list%insert_item(em_dms_lnd_time_axis)

      call dms_ocn_time_axis%initialise("dms_ocn_time",                 &
                                      file_id="dms_conc_ocean_ancil",   &
                                      interp_flag=interp_flag,          &
                                      pop_freq="daily")
      call setup_ancil_field("dms_conc_ocean", depository, ancil_fields,     &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=dms_ocn_time_axis)
      call ancil_times_list%insert_item(dms_ocn_time_axis)

      call em_mterp_time_axis%initialise("em_mterp_time",               &
                                       file_id="emiss_monoterp_ancil",  &
                                       interp_flag=interp_flag,         &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_monoterp", depository, ancil_fields,     &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_mterp_time_axis)
      call ancil_times_list%insert_item(em_mterp_time_axis)

      call em_om_bf_time_axis%initialise("em_om_bf_time",                &
                                       file_id="emiss_om_biofuel_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_om_biofuel", depository, ancil_fields,   &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_om_bf_time_axis)
      call ancil_times_list%insert_item(em_om_bf_time_axis)

      call em_om_ff_time_axis%initialise("em_om_ff_time",                &
                                       file_id="emiss_om_fossil_ancil",  &
                                       interp_flag=interp_flag,          &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_om_fossil", depository, ancil_fields,   &
                           mesh, twod_mesh, twod=.true.,              &
                           time_axis=em_om_ff_time_axis)
      call ancil_times_list%insert_item(em_om_ff_time_axis)

      call em_so2_lo_time_axis%initialise("em_so2_lo_time",              &
                                        file_id="emiss_so2_low_ancil",   &
                                        interp_flag=interp_flag,         &
                                        pop_freq="daily")
      call setup_ancil_field("emiss_so2_low", depository, ancil_fields,     &
                           mesh, twod_mesh, twod=.true.,              &
                           time_axis=em_so2_lo_time_axis)
      call ancil_times_list%insert_item(em_so2_lo_time_axis)

      if (emissions == emissions_GC3) then
        call em_so2_hi_time_axis%initialise("em_so2_hi_time",              &
                                          file_id="emiss_so2_high_ancil",  &
                                          interp_flag=interp_flag,         &
                                          pop_freq="daily")
        call setup_ancil_field("emiss_so2_high", depository, ancil_fields,  &
                                 mesh, twod_mesh, twod=.true.,              &
                                 time_axis=em_so2_hi_time_axis)
        call ancil_times_list%insert_item(em_so2_hi_time_axis)
      else if (emissions == emissions_GC5) then
        call em_bc_bb_hi_time_axis%initialise("em_bc_bb_hi_time",             &
                                         file_id="emiss_bc_biomass_hi_ancil", &
                                         interp_flag=interp_flag,             &
                                         pop_freq="daily")
        call setup_ancil_field("emiss_bc_biomass_high", depository, ancil_fields,&
                             mesh, twod_mesh, twod=.true.,                       &
                             time_axis=em_bc_bb_hi_time_axis)
        call ancil_times_list%insert_item(em_bc_bb_hi_time_axis)

        call em_bc_bb_lo_time_axis%initialise("em_bc_bb_lo_time",             &
                                         file_id="emiss_bc_biomass_lo_ancil", &
                                         interp_flag=interp_flag,             &
                                         pop_freq="daily")
        call setup_ancil_field("emiss_bc_biomass_low", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,                       &
                             time_axis=em_bc_bb_lo_time_axis)
        call ancil_times_list%insert_item(em_bc_bb_lo_time_axis)

        call em_om_bb_hi_time_axis%initialise("em_om_bb_hi_time",             &
                                         file_id="emiss_om_biomass_hi_ancil", &
                                         interp_flag=interp_flag,             &
                                         pop_freq="daily")
        call setup_ancil_field("emiss_om_biomass_high", depository, ancil_fields,&
                             mesh, twod_mesh, twod=.true.,                       &
                             time_axis=em_om_bb_hi_time_axis)
        call ancil_times_list%insert_item(em_om_bb_hi_time_axis)

        call em_om_bb_lo_time_axis%initialise("em_om_bb_lo_time",             &
                                         file_id="emiss_om_biomass_lo_ancil", &
                                         interp_flag=interp_flag,             &
                                         pop_freq="daily")
        call setup_ancil_field("emiss_om_biomass_low", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,                       &
                             time_axis=em_om_bb_lo_time_axis)
        call ancil_times_list%insert_item(em_om_bb_lo_time_axis)
      end if ! GC3 or GC5

      ! -- 3-D ancils
      !-- natural SO2 emissions, currently single-time
      if (init_option == init_option_fd_start_dump .and. &
           .not. checkpoint_read) then
        call setup_ancil_field("emiss_so2_nat", depository, ancil_fields,     &
                               mesh, twod_mesh)
      end if

      if (emissions == emissions_GC3) then
        call em_bc_bb_time_axis%initialise("em_bc_bb_time",                &
                                         file_id="emiss_bc_biomass_ancil", &
                                         interp_flag=interp_flag,          &
                                         pop_freq="daily")
        call setup_ancil_field("emiss_bc_biomass", depository, ancil_fields,  &
                                mesh, twod_mesh,                              &
                                time_axis=em_bc_bb_time_axis)   ! 3-D
        call ancil_times_list%insert_item(em_bc_bb_time_axis)

        call em_om_bb_time_axis%initialise("em_om_bb_time",                &
                                         file_id="emiss_om_biomass_ancil", &
                                         interp_flag=interp_flag,          &
                                         pop_freq="daily")
        call setup_ancil_field("emiss_om_biomass", depository, ancil_fields,  &
                                mesh, twod_mesh,                              &
                                time_axis=em_om_bb_time_axis)
        call ancil_times_list%insert_item(em_om_bb_time_axis)
      end if

      !=====  OFFLINE OXIDANT ANCILS  =====
      if ( chem_scheme == chem_scheme_offline_ox ) then
      call h2o2_limit_time_axis%initialise("h2o2_limit_time",         &
                                         file_id="h2o2_limit_ancil",  &
                                         interp_flag=interp_flag,     &
                                         pop_freq="daily")
      call setup_ancil_field("h2o2_limit", depository, ancil_fields,        &
                           mesh, twod_mesh,                           &
                           time_axis=h2o2_limit_time_axis)
      call ancil_times_list%insert_item(h2o2_limit_time_axis)

      call ho2_time_axis%initialise("ho2_time", file_id="ho2_ancil", &
                                    interp_flag=interp_flag,         &
                                    pop_freq="daily")
      call setup_ancil_field("ho2", depository, ancil_fields,               &
                           mesh, twod_mesh, time_axis=ho2_time_axis)
      call ancil_times_list%insert_item(ho2_time_axis)

      call no3_time_axis%initialise("no3_time", file_id="no3_ancil", &
                                    interp_flag=interp_flag,         &
                                    pop_freq="daily")
      call setup_ancil_field("no3", depository, ancil_fields,               &
                           mesh, twod_mesh, time_axis=no3_time_axis)
      call ancil_times_list%insert_item(no3_time_axis)

      call o3_time_axis%initialise("o3_time", file_id="o3_ancil",    &
                                   interp_flag=interp_flag,          &
                                   pop_freq="daily")
      call setup_ancil_field("o3", depository, ancil_fields,                &
                           mesh, twod_mesh, time_axis=o3_time_axis)
      call ancil_times_list%insert_item(o3_time_axis)

      call oh_time_axis%initialise("oh_time", file_id="oh_ancil",    &
                                   interp_flag=interp_flag,          &
                                   pop_freq="daily")
      call setup_ancil_field("oh", depository, ancil_fields,                &
                           mesh, twod_mesh, time_axis=oh_time_axis)
      call ancil_times_list%insert_item(oh_time_axis)


      end if  ! Offline Oxidants
    endif  ! ancil_updating, glomap_ukca
    ! ====== Easy aerosols ======
    ! CDNC
    if ( easyaerosol_cdnc .and. ancil_option == ancil_option_updating )  then
       call cloud_drop_no_conc_time_axis%initialise("cloud_drop_no_conc_time", &
                                       file_id="cloud_drop_no_conc_ancil",     &
                                       interp_flag=interp_flag,                &
                                       pop_freq="five_days")
       call setup_ancil_field("cloud_drop_no_conc", depository, ancil_fields,  &
                           mesh, twod_mesh,                                    &
                           time_axis=cloud_drop_no_conc_time_axis)   ! 3-D
       call ancil_times_list%insert_item(cloud_drop_no_conc_time_axis)
    endif ! easyaerosol_cdnc

    ! SW ASYMMETRY
    if ( easyaerosol_sw .and. ancil_option == ancil_option_updating )  then
       call easy_asymmetry_sw_time_axis%initialise("easy_asymmetry_sw_time",   &
                                       file_id="easy_asymmetry_sw_ancil",      &
                                       interp_flag=interp_flag,                &
                                       pop_freq="five_days")
       call setup_ancil_field("easy_asymmetry_sw", depository, ancil_fields,   &
                           mesh, twod_mesh,                                    &
                           ndata=n_sw_band, ndata_first=.true.,&
                           time_axis=easy_asymmetry_sw_time_axis)   ! 3-D
       call ancil_times_list%insert_item(easy_asymmetry_sw_time_axis)

    ! SW ABSORPTION
       call easy_absorption_sw_time_axis%initialise("easy_absorption_sw_time", &
                                       file_id="easy_absorption_sw_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
       call setup_ancil_field("easy_absorption_sw", depository, ancil_fields,  &
                           mesh, twod_mesh,                                    &
                           ndata=n_sw_band, ndata_first=.true.,                &
                           time_axis=easy_absorption_sw_time_axis)   ! 3-D
       call ancil_times_list%insert_item(easy_absorption_sw_time_axis)

    ! SW EXTINCTION
       call easy_extinction_sw_time_axis%initialise("easy_extinction_sw_time", &
                                       file_id="easy_extinction_sw_ancil",     &
                                       interp_flag=interp_flag,                &
                                       pop_freq="five_days")
       call setup_ancil_field("easy_extinction_sw", depository, ancil_fields,  &
                           mesh, twod_mesh,                                    &
                           ndata=n_sw_band, ndata_first=.true.,&
                           time_axis=easy_extinction_sw_time_axis)   ! 3-D
       call ancil_times_list%insert_item(easy_extinction_sw_time_axis)
    endif ! easyaerosol_sw

    ! LW ASYMMETRY
    if ( easyaerosol_lw .and. ancil_option == ancil_option_updating )  then
       call easy_asymmetry_lw_time_axis%initialise("easy_asymmetry_lw_time", &
                                       file_id="easy_asymmetry_lw_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
       call setup_ancil_field("easy_asymmetry_lw", depository, ancil_fields,   &
                           mesh, twod_mesh,                                    &
                           ndata=n_lw_band, ndata_first=.true.,&
                           time_axis=easy_asymmetry_lw_time_axis)   ! 3-D
       call ancil_times_list%insert_item(easy_asymmetry_lw_time_axis)

    ! LW ABSORPTION
       call easy_absorption_lw_time_axis%initialise("easy_absorption_lw_time", &
                                       file_id="easy_absorption_lw_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
       call setup_ancil_field("easy_absorption_lw", depository, ancil_fields,   &
                           mesh, twod_mesh,                                    &
                           ndata=n_lw_band, ndata_first=.true.,&
                           time_axis=easy_absorption_lw_time_axis)   ! 3-D
       call ancil_times_list%insert_item(easy_absorption_lw_time_axis)
       call easy_extinction_lw_time_axis%initialise("easy_extinction_lw_time", &
                                       file_id="easy_extinction_lw_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")

    ! LW EXTINCTION
       call setup_ancil_field("easy_extinction_lw", depository, ancil_fields,   &
                           mesh, twod_mesh,                                    &
                           ndata=n_lw_band, ndata_first=.true.,&
                           time_axis=easy_extinction_lw_time_axis)   ! 3-D
       call ancil_times_list%insert_item(easy_extinction_lw_time_axis)
    endif ! easyaerosol_lw

    ! === Chemistry Ancils ======
    if ( (chem_scheme == chem_scheme_strattrop .or.                  &
          chem_scheme == chem_scheme_strat_test)  .and.             &
         ancil_option == ancil_option_updating) then

      call em_c2h6_time_axis%initialise("em_c2h6_time",              &
                                        file_id="emiss_c2h6_ancil",  &
                                        interp_flag=interp_flag,     &
                                        pop_freq="daily")
      call setup_ancil_field("emiss_c2h6", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,           &
                             time_axis=em_c2h6_time_axis)
      call ancil_times_list%insert_item(em_c2h6_time_axis)

      call em_c3h8_time_axis%initialise("em_c3h8_time",              &
                                        file_id="emiss_c3h8_ancil",  &
                                        interp_flag=interp_flag,     &
                                        pop_freq="daily")
      call setup_ancil_field("emiss_c3h8", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,           &
                             time_axis=em_c3h8_time_axis)
      call ancil_times_list%insert_item(em_c3h8_time_axis)

      call em_c5h8_time_axis%initialise("em_c5h8_time",              &
                                        file_id="emiss_c5h8_ancil",  &
                                        interp_flag=interp_flag,     &
                                        pop_freq="daily")
      call setup_ancil_field("emiss_c5h8", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,           &
                             time_axis=em_c5h8_time_axis)
      call ancil_times_list%insert_item(em_c5h8_time_axis)

      call em_ch4_time_axis%initialise("em_ch4_time",                &
                                       file_id="emiss_ch4_ancil",    &
                                       interp_flag=interp_flag,      &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_ch4", depository, ancil_fields,  &
                             mesh, twod_mesh, twod=.true.,           &
                             time_axis=em_ch4_time_axis)
      call ancil_times_list%insert_item(em_ch4_time_axis)

      call em_co_time_axis%initialise("em_co_time",                  &
                                      file_id="emiss_co_ancil",      &
                                      interp_flag=interp_flag,       &
                                      pop_freq="daily")
      call setup_ancil_field("emiss_co", depository, ancil_fields,   &
                             mesh, twod_mesh, twod=.true.,           &
                            time_axis=em_co_time_axis)
      call ancil_times_list%insert_item(em_co_time_axis)

      call em_hcho_time_axis%initialise("em_hcho_time",              &
                                       file_id="emiss_hcho_ancil",   &
                                       interp_flag=interp_flag,      &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_hcho", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,           &
                             time_axis=em_hcho_time_axis)
      call ancil_times_list%insert_item(em_hcho_time_axis)

      call em_me2co_time_axis%initialise("em_me2co_time",            &
                                        file_id="emiss_me2co_ancil", &
                                        interp_flag=interp_flag,     &
                                        pop_freq="daily")
      call setup_ancil_field("emiss_me2co", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,            &
                             time_axis=em_me2co_time_axis)
      call ancil_times_list%insert_item(em_me2co_time_axis)

      call em_mecho_time_axis%initialise("em_mecho_time",             &
                                         file_id="emiss_mecho_ancil", &
                                         interp_flag=interp_flag,     &
                                         pop_freq="daily")
      call setup_ancil_field("emiss_mecho", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,            &
                             time_axis=em_mecho_time_axis)
      call ancil_times_list%insert_item(em_mecho_time_axis)

      call em_nh3_time_axis%initialise("em_nh3_time",                 &
                                       file_id="emiss_nh3_ancil",     &
                                       interp_flag=interp_flag,       &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_nh3", depository, ancil_fields,   &
                             mesh, twod_mesh, twod=.true.,            &
                             time_axis=em_nh3_time_axis)
      call ancil_times_list%insert_item(em_nh3_time_axis)

      call em_no_time_axis%initialise("em_no_time", file_id="emiss_no_ancil",  &
                                     interp_flag=interp_flag,         &
                                     pop_freq="daily")
      call setup_ancil_field("emiss_no", depository, ancil_fields,    &
                             mesh, twod_mesh, twod=.true.,            &
                             time_axis=em_no_time_axis)
      call ancil_times_list%insert_item(em_no_time_axis)

      call em_meoh_time_axis%initialise("em_meoh_time",               &
                                       file_id="emiss_meoh_ancil",    &
                                       interp_flag=interp_flag,       &
                                       pop_freq="daily")
      call setup_ancil_field("emiss_meoh", depository, ancil_fields,  &
                             mesh, twod_mesh, twod=.true.,            &
                             time_axis=em_meoh_time_axis)
      call ancil_times_list%insert_item(em_meoh_time_axis)

      call em_no_aircrft_time_axis%initialise("em_no_aircrft_time",    &
                                              file_id="emiss_no_aircrft_ancil",&
                                              interp_flag=interp_flag, &
                                              pop_freq="daily")
      call setup_ancil_field("emiss_no_aircrft", depository, ancil_fields, &
                             mesh, twod_mesh,                              &
                             time_axis=em_no_aircrft_time_axis)
      call ancil_times_list%insert_item(em_no_aircrft_time_axis)

    endif  ! if chem_scheme_strattrop/ strat_test

  end subroutine create_fd_ancils


  !> @details Organises fields to be read from ancils into ancil_fields
  !           collection then reads them in an idealised setup.
  !> @param[in,out] depository The depository field collection
  !> @param[in,out] ancil_fields Collection for ancillary fields
  !> @param[in] mesh      The current 3d mesh
  !> @param[in] twod_mesh The current 2d mesh
  subroutine create_fd_ancils_idealised( depository, ancil_fields, &
                                         mesh, twod_mesh )

    implicit none

    type( field_collection_type ), intent( inout ) :: depository
    type( field_collection_type ), intent( inout ) :: ancil_fields

    type( mesh_type ), intent(in), pointer :: mesh
    type( mesh_type ), intent(in), pointer :: twod_mesh

    ! Set up ancil_fields collection
    write(log_scratch_space,'(A,A)') "Create ancil fields: "// &
          "Setting up ancil field collection"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    ! Here ancil fields are set up with a call to setup_ancil_field.
    if ( ancil_option == ancil_option_idealised ) then
      if (ch4_rad_opt == ch4_rad_opt_ancil) then
        call setup_ancil_field("ch4", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (co_rad_opt == co_rad_opt_ancil) then
        call setup_ancil_field("co", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (co2_rad_opt == co2_rad_opt_ancil) then
        call setup_ancil_field("co2", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (h2_rad_opt == h2_rad_opt_ancil) then
        call setup_ancil_field("h2", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (h2o_rad_opt == h2o_rad_opt_ancil) then
        call setup_ancil_field("h2o", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (hcn_rad_opt == hcn_rad_opt_ancil) then
        call setup_ancil_field("hcn", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (he_rad_opt == he_rad_opt_ancil) then
        call setup_ancil_field("he", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (n2_rad_opt == n2_rad_opt_ancil) then
        call setup_ancil_field("n2", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (nh3_rad_opt == nh3_rad_opt_ancil) then
        call setup_ancil_field("nh3", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (o2_rad_opt == o2_rad_opt_ancil) then
        call setup_ancil_field("o2", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (so2_rad_opt == so2_rad_opt_ancil) then
        call setup_ancil_field("so2", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (cs_rad_opt == cs_rad_opt_ancil) then
        call setup_ancil_field("cs", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (k_rad_opt == k_rad_opt_ancil) then
        call setup_ancil_field("k", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (li_rad_opt == li_rad_opt_ancil) then
        call setup_ancil_field("li", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (na_rad_opt == na_rad_opt_ancil) then
        call setup_ancil_field("na", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (rb_rad_opt == rb_rad_opt_ancil) then
        call setup_ancil_field("rb", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (tio_rad_opt == tio_rad_opt_ancil) then
        call setup_ancil_field("tio", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
      if (vo_rad_opt == vo_rad_opt_ancil) then
        call setup_ancil_field("vo", depository, ancil_fields, &
          mesh, twod_mesh)
      endif
    endif

    ! Now the field collection is set up, the fields will be initialised in
    ! gungho_model_data_mod

    ! Prescribed forcing of surface temperature e.g. via internal flux
    if ( surf_temp_forcing == surf_temp_forcing_int_flux .and. &
         internal_flux_method == internal_flux_method_non_uniform ) then
      call setup_ancil_field("internal_flux", depository, ancil_fields, &
        mesh, twod_mesh)
    end if

  end subroutine create_fd_ancils_idealised

  !> @details Adds fields to the ancil collection, sets up their read and write
  !>      behaviour and creates them in the depository if they do not yet exist
  !> @param[in] name The field name
  !> @param[in, out] depository The depository field collection
  !> @param[in, out] ancil_fields The ancil field collection
  !> @param[in] mesh                The current 3d mesh
  !> @param[in, optional] twod_mesh The current 2d mesh
  !> @param[in, optional] ndata Number of non-spatial dimensions for multi-data
  !>                            field
  !> @param[in, out, optional] time_axis Time axis associated with ancil field
  !> @param[in, optional] alt_mesh      Alternative 3d mesh for time axis fields
  !> @param[in, optional] alt_twod_mesh Alternative 2d mesh for time axis fields
  subroutine setup_ancil_field( name, depository, ancil_fields, mesh, &
                                twod_mesh, twod, ndata, ndata_first,  &
                                time_axis, alt_mesh, alt_twod_mesh  )

    implicit none

    character(*),                   intent(in)          :: name
    type( field_collection_type ),  intent(inout)       :: depository
    type( field_collection_type ),  intent(inout)       :: ancil_fields
    type( mesh_type ),    pointer,  intent(in)          :: mesh
    type( mesh_type ),    pointer,  intent(in)          :: twod_mesh
    logical(l_def),       optional, intent(in)          :: twod
    integer(i_def),       optional, intent(in)          :: ndata
    logical(l_def),       optional, intent(in)          :: ndata_first
    type(time_axis_type), optional, intent(inout)       :: time_axis
    type( mesh_type ), optional, pointer, intent(in)    :: alt_mesh
    type( mesh_type ), optional, pointer, intent(in)    :: alt_twod_mesh

    ! Local variables
    type(field_type)          :: new_field
    integer(i_def)            :: ndat, time_ndat
    logical(l_def)            :: twod_field
    logical(l_def)            :: ndat_first
    integer(i_def), parameter :: fs_order_h = 0
    integer(i_def), parameter :: fs_order_v = 0

    ! Pointers
    type(function_space_type),       pointer :: vec_space => null()
    procedure(read_interface),       pointer :: tmp_read_ptr => null()
    procedure(write_interface),      pointer :: tmp_write_ptr => null()
    type(field_type),                pointer :: fld_ptr => null()
    class(pure_abstract_field_type), pointer :: abs_fld_ptr => null()

    ! Set field ndata if argument is present, else leave as default value
    if (present(ndata)) then
      ndat = ndata
    else
      ndat = 1
    end if
    if (present(twod)) then
      twod_field = twod
    else
      twod_field = .false.
    end if
    if (present(ndata_first)) then
       ndat_first = ndata_first
    else
       ndat_first = .false.
    end if

    ! If field does not yet exist, then create it
    if ( .not. depository%field_exists( name ) ) then
      write(log_scratch_space,'(3A,I6)') &
           "Creating new field for ", trim(name)
      call log_event(log_scratch_space,LOG_LEVEL_INFO)
      tmp_write_ptr => write_field_generic
      if (twod_field) then
        vec_space => function_space_collection%get_fs( twod_mesh, fs_order_h, &
                                                       fs_order_v,            &
                                                       W3, ndat )
      else
        vec_space => function_space_collection%get_fs( mesh, fs_order_h, &
                                                       fs_order_v,       &
                                                       WTheta, ndat )
      end if
      call new_field%initialise( vec_space, name=trim(name))
      call new_field%set_write_behaviour(tmp_write_ptr)
      ! Add the new field to the field depository
      call depository%add_field(new_field)
    end if

    ! If field is time-varying, also create field storing raw data to be
    ! interpolated
    if ( present(time_axis) ) then
      write(log_scratch_space,'(3A,I6)') &
           "Creating time axis field for ", trim(name)
      call log_event(log_scratch_space,LOG_LEVEL_INFO)

      ! Multiply ndat by the number of time windows
      time_ndat = ndat * time_axis%get_window_size()
      if (twod_field) then
        if ( present(alt_twod_mesh) ) then
          vec_space => function_space_collection%get_fs( alt_twod_mesh, fs_order_h, &
                                                         fs_order_v,                &
                                                         W3, time_ndat,             &
                                                         ndata_first )
        else
          vec_space => function_space_collection%get_fs( twod_mesh, fs_order_h, &
                                                         fs_order_v,            &
                                                         W3, time_ndat,         &
                                                         ndata_first )
        end if
      else
        if ( present(alt_mesh) ) then
          vec_space => function_space_collection%get_fs( alt_mesh, fs_order_h, &
                                                         fs_order_v,           &
                                                         Wtheta, time_ndat,    &
                                                         ndata_first )
        else
          vec_space => function_space_collection%get_fs( mesh, fs_order_h,     &
                                                         fs_order_v,           &
                                                         Wtheta, time_ndat,    &
                                                         ndata_first )
        end if
      end if
      call new_field%initialise( vec_space, name=trim(name) )
      call time_axis%add_field(new_field)
    end if

    ! Get a field pointer from the depository
    call depository%get_field(name, fld_ptr)

    if (.not. present(time_axis)) then
      !Set up field read behaviour for 2D and 3D fields
      tmp_read_ptr => read_field_generic
      ! Set field read behaviour for target field
      call fld_ptr%set_read_behaviour(tmp_read_ptr)
    end if

    ! Add the field pointer to the target field collection
    call depository%get_field(name, fld_ptr)
    abs_fld_ptr => fld_ptr
    call ancil_fields%add_reference_to_field(abs_fld_ptr)

    ! Nullify pointers
    nullify(vec_space)
    nullify(tmp_read_ptr)
    nullify(tmp_write_ptr)
    nullify(fld_ptr)

  end subroutine setup_ancil_field

end module init_ancils_mod
