# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Bespoke PSyclone transformation script for sw_rad_tile_kernel_mod.
'''

import logging
from psyclone.transformations import (
    OMPLoopTrans,
    TransformationError)
from psyclone.psyir.nodes import Loop


omp_transform_par_do = OMPLoopTrans(
    omp_schedule="static",
    omp_directive="paralleldo")

options = {"ignore_dependencies_for": [
    "albedo_obs_scaling",
    "tile_sw_direct_albedo",
    "tile_sw_diffuse_albedo",
    "sea_ice_pensolar_frac_direct",
    "sea_ice_pensolar_frac_diffuse",
    ],
    "node-type-check": False}


def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply OMP paralleldo transformations
    to each loop in sw_rad_tile_kernel_mod.
    '''

    for loop in psyir.walk(Loop):
        if not loop.ancestor(Loop):
            try:
                omp_transform_par_do.apply(loop, options)

            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)

# Ignore loops setting these as order dependent:
#     land_field l ainfo%land_index sea_pts
#     ainfo%sea_index ainfo%sice_pts_ncat ainfo%sice_index_ncat
