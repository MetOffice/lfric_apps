# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Lift options list and similar from each individual script up into this file.
Aim is to allow the creation of a simple global.py which adds OMP over all
loops. The option which matches the file being worked on can be pulled in
and referenced. This helps reduce the number of files needed.

The following files have overrides below:
* bl_lsp
* btq_int
* ex_flux_tq
* ex_flux_uv
* kmkh
* tr_mix
* imp_mix
* fm_drag
'''

FILE_EXTEN = ".xu90"

# Basic initialisation, will be used by the local script
SCRIPT_OPTIONS_DICT = {}

# ## Local.py options for boundary layer ##

SCRIPT_OPTIONS_DICT["bl_lsp"+str(FILE_EXTEN)] = {
    "first_private_list": [
        "newqcf"
        ]
}

SCRIPT_OPTIONS_DICT["ex_flux_tq"+str(FILE_EXTEN)] = {
    "options": {
        "ignore_dependencies_for": [
            "bl_diag%grad_t_adj"
            ]
        },
    "first_private_list": [
        "f2_fqw",
        "f2_ftl",
        "fsc_fqw",
        "fsc_ftl",
        "non_grad_fqw",
        "non_grad_ftl"
        ]
}

SCRIPT_OPTIONS_DICT["ex_flux_uv"+str(FILE_EXTEN)] = {
    "options": {
        "ignore_dependencies_for": [
            "tau_x_y"
            ]
        }
}

SCRIPT_OPTIONS_DICT["tr_mix"+str(FILE_EXTEN)] = {
    "options": {
        "ignore_dependencies_for": [
            "f_field",
            "rhok_dep",
            "surf_dep_flux",
            "gamma_rhokh_rdz"
            ]
        },
    "first_private_list": [
        "dfield_inv",
        "dz_disc",
        "dzlkp1",
        "f_field_ent",
        "km1"
        ]
}

SCRIPT_OPTIONS_DICT["bdy_impl3"+str(FILE_EXTEN)] = {
    "options": {
        "ignore_dependencies_for": [
            "dtl1_1",
            "ct_prod",
            "dqw1_1",
            "dqw",
            "dtl",
            "ct_ctq",
            "temp",
            "dqw1",
            "dtl1",
            "ctctq1"
            ]
        },
    "safe_pure_calls": [
        "oneover_v"
        ],
    "max_threads_parse": True,
    "loop_type_init": [
        "max_threads",
        "blm1",
        "tdims_omp_block",
        "tdims_seg_block"
        ]
}

SCRIPT_OPTIONS_DICT["bdy_impl4"+str(FILE_EXTEN)] = {
    "options": {
        "ignore_dependencies_for": [
            "tl",
            "qw",
            "dtl",
            "dqw"
            ]
        },
    "loop_type_init": [
        "tdims_omp_block",
        "tdims_seg_block"
        ]
}

SCRIPT_OPTIONS_DICT["imp_mix"+str(FILE_EXTEN)] = {
    "options": {
        "ignore_dependencies_for": [
            "d_field",
            "af",
            "field",
            "gamma_rhok_dep",
            "f_field",
            "surf_dep_flux"
            ]
        },
    "max_threads_parse": True,
    "loop_type_init": [
        "max_threads",
        "blm1",
        "pdims_omp_block",
        "pdims_seg_block"
        ]
}

SCRIPT_OPTIONS_DICT["fm_drag"+str(FILE_EXTEN)] = {
    "options": {
        "ignore_dependencies_for": [
            "k_for_buoy",
            "u_hm",
            "v_hm",
            "tl_hm",
            "qw_hm",
            "tau_fd_x",
            "tau_fd_y"
            ]
        },
    "first_private_list": [
        "fp_x_low",
        "fp_x_steep",
        "fp_y_low",
        "fp_y_steep",
        "rib_fn",
        "tausx",
        "tausy",
        "wta",
        "wtb"
        ]
}

# ## Other boundary layer options ##
