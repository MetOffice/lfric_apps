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
