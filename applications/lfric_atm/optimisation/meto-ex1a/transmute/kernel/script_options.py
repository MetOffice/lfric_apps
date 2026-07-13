# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
This file lifts optional overrides, where possible, into a localised,
single location.
The goal is to allow the creation of a simpler local.py which controls how OMP
is added around loops.
SCRIPT_OPTIONS_DICT is pulled into the local.py script, and overrides some
local functionality.
Each key represents a different kernel file which uses the local script.
Each file key in turn contains a dictionary of overrides.
'''

# Needs to be lifted and likely set by the build system longer term.
# The filename passed to PSyclone, this is the pre-processed FTN source.
FILE_EXTEN = ".xu90"

# Basic initialisation
SCRIPT_OPTIONS_DICT = {}

# File keys
SCRIPT_OPTIONS_DICT["aerosol_ukca_kernel_mod"+str(FILE_EXTEN)] = {
        "ignore_dependencies_for": [
            "surf_wetness", "l_tile_active", "o3p", "o3",
            "n", "no", "no3", "lumped_n", "n2o5", "ho2no2",
            "hono2", "h2o2", "ch4", "co", "hcho", "meoo",
            "meooh", "h", "oh", "ho2", "cl", "cl2o2", "clo",
            "oclo", "br", "lumped_br", "brcl", "brono2",
            "n2o", "lumped_cl", "hocl", "hbr", "hobr", "hobr",
            "clono2", "cfcl3", "cf2cl2", "mebr", "hono", "c2h6",
            "etoo", "etooh", "mecho", "meco3", "pan", "c3h8",
            "n_proo", "i_proo", "i_proo", "i_prooh", "etcho",
            "etco3", "me2co", "mecoch2oo", "mecoch2ooh", "ppan",
            "meono2", "c5h8", "iso2", "isooh", "ison", "macr",
            "macrooh", "macro2", "mpan", "hacet", "mgly",
            "nald", "hcooh", "meco3h", "meco2h", "h2", "meoh",
            "n_prooh", "msa", "nh3", "cs2", "csul", "h2s", "so3",
            "passive_o3", "age_of_air", "dms", "so2", "h2so4",
            "dmso", "monoterpene", "secondary_organic", "n_nuc_sol",
            "nuc_sol_su", "nuc_sol_om", "n_ait_sol", "ait_sol_su",
            "nuc_sol_om", "n_ait_sol", "ait_sol_su", "ait_sol_om",
            "n_acc_sol", "acc_sol_su", "acc_sol_bc", "acc_sol_om",
            "acc_sol_ss", "n_cor_sol", "cor_sol_su", "cor_sol_bc",
            "cor_sol_om", "cor_sol_ss", "n_ait_ins", "ait_ins_bc",
            "ait_ins_om", "n_acc_ins", "acc_ins_du", "n_cor_ins",
            "cor_ins_du", "cloud_drop_no_conc", "drydp_ait_sol",
            "drydp_acc_sol", "drydp_cor_sol", "drydp_ait_ins",
            "drydp_acc_ins", "drydp_cor_ins", "wetdp_ait_sol",
            "wetdp_acc_sol", "wetdp_cor_sol", "rhopar_ait_sol",
            "rhopar_acc_sol", "rhopar_cor_sol", "rhopar_ait_ins",
            "rhopar_acc_ins", "rhopar_cor_ins", "pvol_su_ait_sol",
            "pvol_bc_ait_sol", "pvol_om_ait_sol", "pvol_wat_ait_sol",
            "pvol_su_acc_sol", "pvol_bc_acc_sol", "pvol_om_acc_sol",
            "pvol_ss_acc_sol", "pvol_wat_acc_sol", "pvol_su_cor_sol",
            "pvol_bc_cor_sol", "pvol_om_cor_sol", "pvol_ss_cor_sol",
            "pvol_wat_cor_sol", "pvol_bc_ait_ins", "pvol_om_ait_ins",
            "pvol_du_acc_ins", "pvol_du_cor_ins", "no2", "bro", "hcl",
            "o1d", "ait_sol_bc",
    ]
}

SCRIPT_OPTIONS_DICT["bl_exp_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "tnuc_nlcl", "dw_bl", "surf_interp", "bq_bl", "bt_bl",
        "dtrdz_tq_bl", "level_ent", "level_ent_dsc", "ent_we_lim",
        "ent_t_frac", "ent_zrzi", "ent_we_lim_dsc", "ent_t_frac_dsc",
        "ent_zrzi_dsc", "fd_taux", "fd_tauy", "gradrinr", "rhokm_bl",
        "rhokh_bl", "moist_flux_bl", "heat_flux_bl", "gradrinr", "lmix_bl",
        "ngstress_bl", "tke_bl", "bl_weight_1dbl", "zh_nonloc", "z_lcl",
        "inv_depth", "qcl_at_inv_top", "shallow_flag", "uw0_flux", "vw0_flux",
        "lcl_height", "parcel_top", "level_parcel_top", "wstar_2d", "thv_flux",
        "parcel_buoyancy", "qsat_at_lcl", "bl_type_ind", "visc_m_blend",
        "visc_h_blend", "zh_2d", "zhsc_2d", "ntml_2d", "cumulus_2d",
        "rh_crit", "mix_len_bm", "dsldzm", "wvar", "zht", "oblen",
        ]
}

SCRIPT_OPTIONS_DICT["bl_imp_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "dqw_wth", "dtl_wth", "dqw_nt_wth", "dtl_nt_wth",
        "ct_ctq_wth", "qw_wth", "tl_wth", "dqw1_2d",
        "dtl1_2d", "ct_ctq1_2d",
        ]
}

SCRIPT_OPTIONS_DICT["bl_imp2_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "fric_heating_blyr", "fric_heating_incv", "nblyr", "cca",
        "ccw", "z_lcl", "cf_area", "cf_bulk", "cf_ice", "cf_liq",
        "dtheta_bl", "m_v", "m_cl", "m_s", "fqw_star_w3", "ftl_star_w3",
        "heat_flux_bl", "moist_flux_bl", "rhokh_bl", "dqw_wth",
        "dtl_wth", "qw_wth", "tl_wth"
        ]
}

SCRIPT_OPTIONS_DICT["bm_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "theta_inc", "m_v", "m_cl", "cf_bulk", "cf_liq", "cf_ice", "cf_area",
        "sskew_bm", "svar_bm", "svar_tb",
        ]
}

SCRIPT_OPTIONS_DICT["bm_tau_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "tau_dec_bm", "tau_hom_bm", "tau_mph_bm"
        ]
}

SCRIPT_OPTIONS_DICT["conv_comorph_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "precfrac", "cca_2d", "cape_diluted", "lowest_cca_2d",
        "dt_conv", "dmv_conv", "dmcl_conv", "dms_conv", "m_ci",
        "m_r", "m_g", "massflux_up_half", "massflux_up", "massflux_down",
        "pressure_inc_env", "entrain_up", "entrain_down", "detrain_up",
        "du_conv", "dv_conv", "conv_prog_dtheta", "dmv_conv", "dbcf_conv",
        "dcff_conv", "dcfl_conv", "dt_conv", "dmv_conv", "dmcl_conv",
        "dms_conv", "dd_mf_cb", "cca", "ccw", "cv_top", "cv_base",
        "lowest_cv_top", "lowest_cv_base", "pres_cv_top", "pres_cv_base",
        "pres_lowest_cv_top", "pres_lowest_cv_base", "tke_bl", "o3p", "o1d",
        "o3", "nit", "no", "no3", "lumped_n", "n2o5", "ho2no2", "hono2",
        "h2o2", "ch4", "co", "hcho", "meoo", "meooh", "h", "oh", "ho2",
        "cl", "cl2o2", "clo", "oclo", "br", "lumped_br", "brcl", "brono2",
        "n2o", "lumped_cl", "hocl", "hbr", "hobr", "clono2", "cfcl3", "cf2cl2",
        "mebr", "hono", "c2h6", "etoo", "etooh", "mecho", "meco3", "pan",
        "c3h8", "n_proo", "i_proo", "n_prooh", "i_prooh", "etcho", "etcho",
        "me2co", "mecoch2oo", "mecoch2ooh", "ppan", "meono2", "c5h8", "iso2",
        "isooh", "ison", "macr", "macro2", "macrooh", "mpan", "hacet", "mgly",
        "nald", "hcooh", "hcooh", "meco2h", "h2", "meoh", "msa", "nh3", "cs2",
        "csul", "h2s", "so3", "passive_o3", "age_of_air", "dms", "so2",
        "h2so4", "dmso", "monoterpene", "secondary_organic", "n_nuc_sol",
        "nuc_sol_su", "nuc_sol_om", "n_ait_sol", "ait_sol_su", "ait_sol_bc",
        "ait_sol_om", "n_acc_sol", "acc_sol_su", "acc_sol_bc", "acc_sol_om",
        "acc_sol_ss", "n_cor_sol", "cor_sol_su", "cor_sol_bc", "cor_sol_bc",
        "cor_sol_om", "cor_sol_ss", "n_ait_ins", "ait_ins_bc", "ait_ins_om",
        "n_acc_ins", "acc_ins_du", "n_cor_ins", "cor_ins_du", "detrain_down",
        "conv_prog_dmv", "etco3", "etco3", "meco3h",
        ]
}

SCRIPT_OPTIONS_DICT["jules_exp_kernel_mod"+str(FILE_EXTEN)] = {
    "node_type_check": False,
    "ignore_dependencies_for": [
        "z0msea_2d", "tstar_land", "sea_ice_pensolar", "rhostar_2d",
        "recip_l_mo_sea_2d", "t1_sd_2d", "q1_sd_2d", "surf_interp",
        "rhokm_bl", "rhokh_bl", "moist_flux_bl", "heat_flux_bl",
        "gradrinr", "alpha1_tile", "fracaero_t_tile", "fracaero_s_tile",
        "z0h_tile", "z0m_tile", "chr1p5m_tile", "resfs_tile",
        "gc_tile", "canhc_tile", "ashtf_prime_tile", "dtstar_tile",
        "rhokh_tile", "blend_height_tq", "z0m_eff", "ustar",
        "soil_moist_avail", "snow_unload_rate", "tile_temperature",
        "tile_heat_flux", "tile_moisture_flux", "z0m_2d",
        "dust_div_flux", "tile_water_extract", "net_prim_prod",
        "surface_conductance", "thermal_cond_wet_soil", "soil_respiration",
        "gross_prim_prod", "z0h_eff", "fluxes%tstar_ij", "forcing%pstar_ij",
        "qs_star",
        ],
    "safe_pure_calls": ["qsat_mix"]
}

SCRIPT_OPTIONS_DICT["jules_extra_kernel_mod"+str(FILE_EXTEN)] = {
    "node_type_check": False,
    "ignore_dependencies_for": [
        "canopy_water", "tile_snow_mass", "n_snow_layers", "snow_depth",
        "tile_snow_rgrain", "snow_under_canopy", "snowpack_density",
        "snowice_melt", "soil_sat_frac", "water_table", "wetness_under_soil",
        "surface_runoff", "sub_surface_runoff", "soil_moisture_content",
        "grid_snow_mass", "throughfall", "snow_layer_thickness",
        "snow_layer_ice_mass", "snow_layer_liq_mass", "snow_layer_temp",
        "snow_layer_rgrain", "soil_temperature", "soil_moisture",
        "unfrozen_soil_moisture", "frozen_soil_moisture",
    ]
}

SCRIPT_OPTIONS_DICT["jules_imp_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "tstar_land", "sea_ice_pensolar", "ashtf_prime_sea", "dtstar_sea",
        "ashtf_prime", "dtstar_sice", "heat_flux_bl", "moist_flux_bl",
        "tile_heat_flux", "tile_moisture_flux", "tile_temperature",
        "screen_temperature", "tile_heat_flux", "tile_moisture_flux",
        "snowice_sublimation", "surf_heat_flux", "canopy_evap",
        "snowice_melt", "time_since_transition", "surf_ht_flux",
        "water_extraction", "lake_evap", "snomlt_surf_htf", "soil_evap",
        "soil_surf_ht_flux", "t1p5m", "q1p5m", "qcl1p5m", "rh1p5m",
        "t1p5m_ssi", "q1p5m_ssi", "qcl1p5m_ssi", "rh1p5m_ssi",
        "t1p5m_land_loc", "q1p5m_land_loc", "t1p5m_land",
        "q1p5m_land", "qcl1p5m_land", "rh1p5m_land", "t1p5m_surft",
        "q1p5m_surft", "latent_heat", "surf_sw_net", "surf_radnet",
        "surf_lw_up", "surf_lw_down", "sea_ice_temperature", "latent_heat",
    ]
}

SCRIPT_OPTIONS_DICT["lsp_prognostic_tnuc_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "tnuc",
    ]
}

SCRIPT_OPTIONS_DICT["lw_rad_tile_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "tile_lw_grey_albedo", "emis_nvg", "emis_pft",
    ],
    "safe_pure_calls": ["emis_nvg", "emis_pft"]
}

SCRIPT_OPTIONS_DICT["mphys_kernel_mod"+str(FILE_EXTEN)] = {
    "node_type_check": False,
    "ignore_dependencies_for": [
        "dml_wth", "dmv_wth", "dml_wth", "dtheta", "refl_1km", "dms_wth",
        "dmr_wth", "dmg_wth", "murk", "dbcf_wth", "dcff_wth", "dcfl_wth",
        "ls_rain_2d", "ls_snow_2d", "ls_graup_2d", "lsca_2d", "ls_rain_3d",
        "ls_snow_3d", "precfrac", "autoconv", "accretion", "rim_cry",
        "rim_agg", "refl_tot", "superc_liq_wth", "superc_rain_wth",
        "sfwater", "sfrain", "sfsnow",
    ]
}

SCRIPT_OPTIONS_DICT["mphys_turb_gen_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "dbcf", "dcfl", "dmr_cl", "dmr_v", "dtheta",
    ]
}

SCRIPT_OPTIONS_DICT["pc2_initiation_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "dtheta_inc_wth", "dmv_inc_wth", "dmcl_inc_wth", "dmci_inc_wth",
        "dms_inc_wth", "dcfl_inc_wth", "dcff_inc_wth", "dbcf_inc_wth",
        "sskew_bm", "svar_bm", "svar_tb",
    ]
}

SCRIPT_OPTIONS_DICT["photol_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "photol_rates",
    ]
}

SCRIPT_OPTIONS_DICT["sw_rad_tile_kernel_mod"+str(FILE_EXTEN)] = {
    "node-type-check": False,
    "ignore_dependencies_for": [
        "albedo_obs_scaling", "tile_sw_direct_albedo",
        "tile_sw_diffuse_albedo", "sea_ice_pensolar_frac_direct",
        "sea_ice_pensolar_frac_diffuse",
    ]
}

SCRIPT_OPTIONS_DICT["smith_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "theta_inc", "m_v", "m_cl", "cf_bulk", "cf_liq", "cf_ice", "cf_area",
    ]
}

SCRIPT_OPTIONS_DICT["spectral_gwd_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "du_spectral_gwd", "dv_spectral_gwd", "dtemp_spectral_gwd",
        "tau_east_spectral_gwd", "tau_south_spectral_gwd",
        "tau_west_spectral_gwd", "tau_north_spectral_gwd",
    ]
}

SCRIPT_OPTIONS_DICT["tracer_mix_kernel_mod"+str(FILE_EXTEN)] = {
    "ignore_dependencies_for": [
        "tracer",
    ]
}
