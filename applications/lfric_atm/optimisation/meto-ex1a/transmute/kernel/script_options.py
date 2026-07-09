# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Lift options list and similar from each individual script up into this file.
Aim is to allow the creation of a simple global.py which adds OMP over all
loops. The option which matches the file being worked on can be pulled in
and referenced. This reduces the number of files needed
'''

# Needs to be lifted and likely set by the build system longer term.
# The filename passed to PSyclone, this is the pre-processed FTN source.
FILE_EXTEN = ".xu90"

# Basic initialisation
SCRIPT_OPTIONS_DICT = {}

# Kernels
SCRIPT_OPTIONS_DICT["bl_exp_kernel_mod"+str(FILE_EXTEN)] = {
    "false_dep_vars": [
        "tnuc_nlcl",
        "dw_bl",
        "surf_interp",
        "bq_bl",
        "bt_bl",
        "dtrdz_tq_bl",
        "level_ent",
        "level_ent_dsc",
        "ent_we_lim",
        "ent_t_frac",
        "ent_zrzi",
        "ent_we_lim_dsc",
        "ent_t_frac_dsc",
        "ent_zrzi_dsc",
        "fd_taux",
        "fd_tauy",
        "gradrinr",
        "rhokm_bl",
        "rhokh_bl",
        "moist_flux_bl",
        "heat_flux_bl",
        "gradrinr",
        "lmix_bl",
        "ngstress_bl",
        "tke_bl",
        "bl_weight_1dbl",
        "zh_nonloc",
        "z_lcl",
        "inv_depth",
        "qcl_at_inv_top",
        "shallow_flag",
        "uw0_flux",
        "vw0_flux",
        "lcl_height",
        "parcel_top",
        "level_parcel_top",
        "wstar_2d",
        "thv_flux",
        "parcel_buoyancy",
        "qsat_at_lcl",
        "bl_type_ind",
        "visc_m_blend",
        "visc_h_blend",
        "zh_2d",
        "zhsc_2d",
        "ntml_2d",
        "cumulus_2d",
        "rh_crit",
        "mix_len_bm",
        "dsldzm",
        "wvar",
        "zht",
        "oblen",
        ]
}

SCRIPT_OPTIONS_DICT["bl_imp_kernel_mod"+str(FILE_EXTEN)] = {
    "false_dep_vars": [
        "dqw_wth",
        "dtl_wth",
        "dqw_nt_wth",
        "dtl_nt_wth",
        "ct_ctq_wth",
        "qw_wth",
        "tl_wth",
        "dqw1_2d",
        "dtl1_2d",
        "ct_ctq1_2d",
        ]
}

SCRIPT_OPTIONS_DICT["bl_imp2_kernel_mod"+str(FILE_EXTEN)] = {
    "false_dep_vars": [
        "fric_heating_blyr",
        "fric_heating_incv",
        "nblyr",
        "cca",
        "ccw",
        "z_lcl",
        "cf_area",
        "cf_bulk",
        "cf_ice",
        "cf_liq",
        "dtheta_bl",
        "m_v",
        "m_cl",
        "m_s",
        "fqw_star_w3",
        "ftl_star_w3",
        "heat_flux_bl",
        "moist_flux_bl",
        "rhokh_bl",
        "dqw_wth",
        "dtl_wth",
        "qw_wth",
        "tl_wth"
        ]
}

SCRIPT_OPTIONS_DICT["mphys_kernel_mod"+str(FILE_EXTEN)] = {

    "options": {
        "node-type-check": False,
        "ignore_dependencies_for": [
            "dtheta",           # First and Second i, k loop
            "dmv_wth",          # First and Second i, k loop
            "dml_wth",          # First and Second i, k loop
            "dms_wth",          # First and Second i, k loop
            "dmr_wth",          # Third i, k loop
            "dmg_wth",          # Forth i, k loop
            "murk",             # First k, i loop
            "dbcf_wth",         # Fifth i, k loop
            "dcfl_wth",         # Fifth i, k loop
            "dcff_wth",         # Fifth i, k loop
            "ls_rain_2d",       # Fifth i, k loop
            "ls_snow_2d",       # Fifth i, k loop
            "ls_graup_2d",      # Fifth i, k loop
            "lsca_2d",          # Fifth i, k loop
            "ls_rain_3d",       # Fifth i, k loop
            "ls_snow_3d",       # Fifth i, k loop
            "precfrac",         # Fifth i, k loop
            "refl_tot",         # Fifth i, k loop
            "autoconv",         # Fifth i, k loop
            "accretion",        # Fifth i, k loop
            "rim_cry",          # Fifth i, k loop
            "rim_agg",          # Fifth i, k loop
            "refl_1km",         # Fifth i, k loop
            "superc_liq_wth",   # Sixth i, k loop
            "superc_rain_wth",  # Seventh i, k loop
            "sfwater",          # Eighth i, k loop
            "sfwater",          # Second k, i loop
            "sfrain",           # Third k, i loop
            "sfsnow",           # Fourth k, i loop
        ]
    }
}
