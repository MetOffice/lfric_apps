##############################################################################
# Copyright (c) 2025,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Bespoke Opt script for mphys_kernel_mod to add OpenMP to loops present in
the Kernel.
As PSyclone is detecting a number of lhs rhs dependencies due to the
copies in and out of the kernel, we are needing to provide a list of
ignore_dependencies_for in the transformation options list.
'''

import logging
from psyclone.transformations import (
    OMPLoopTrans,
    TransformationError)
from psyclone.psyir.nodes import Loop


omp_transform_par_do = OMPLoopTrans(
    omp_schedule="static",
    omp_directive="paralleldo")


def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply paralleldo transformations
    to each loop in large scale precip.
    '''

    for loop in psyir.walk(Loop):
        if not loop.ancestor(Loop):
            options = {"ignore_dependencies_for": [
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
                ],
                "node-type-check": False}
            try:
                omp_transform_par_do.apply(loop, options)

            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)
