# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Bespoke PSyclone transformation script for jules_extra_kernel_mod.
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
    each schedule (or subroutine) and apply OMP paralleldo transformations
    to each loop in jules_extra_kernel_mod.
    '''

    for loop in psyir.walk(Loop):
        if not loop.ancestor(Loop):
            options = {"ignore_dependencies_for": [
                "canopy_water",
                "tile_snow_mass",
                "n_snow_layers",
                "snow_depth",
                "tile_snow_rgrain",
                "snow_under_canopy",
                "snowpack_density",
                "snowice_melt",
                "soil_sat_frac",
                "water_table",
                "wetness_under_soil",
                "surface_runoff",
                "sub_surface_runoff",
                "soil_moisture_content",
                "grid_snow_mass",
                "throughfall",
                "snow_layer_thickness",
                "snow_layer_ice_mass",
                "snow_layer_liq_mass",
                "snow_layer_temp",
                "snow_layer_rgrain",
                "soil_temperature",
                "soil_moisture",
                "unfrozen_soil_moisture",
                "frozen_soil_moisture",
                ],
                "node-type-check": False}
            try:
                omp_transform_par_do.apply(loop, options)

            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)

# Ignore loops setting these as order dependent:
#    land_pts l ainfo%land_index soil_pts
#    ainfo%soil_index lice_pts ainfo%lice_index
