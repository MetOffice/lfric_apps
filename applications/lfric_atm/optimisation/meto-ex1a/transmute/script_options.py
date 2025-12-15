# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Lift options list and similar from each induvidial script up into this file.
Aim is to allow the creation of a simple global.py which adds OMP over all
loops. The option which matches the file being worked on can be pulled in
and referenced. This reduces the number of files needed
'''

# Needs to be lifted and likely set by the build system longer term.
# The filename passed to PSyclone, this is the pre-processed FTN source.
FILE_EXTEN = ".xu90"

# Basic initialisation, will be used by global script
SCRIPT_OPTIONS_DICT = {}

# Kernels
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
