# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Bespoke PSyclone transformation script for jules_imp_kernel_mod.
'''

import logging
from psyclone.transformations import (
    OMPLoopTrans,
    TransformationError)
from psyclone.psyir.nodes import Loop, IfBlock


omp_transform_par_do = OMPLoopTrans(
    omp_schedule="static",
    omp_directive="paralleldo")

options = {"ignore_dependencies_for": [
        "tstar_land",
        "sea_ice_pensolar",
        "ashtf_prime_sea",
        "dtstar_sea",
        "ashtf_prime",
        "dtstar_sice",
        "heat_flux_bl",
        "moist_flux_bl",
        "tile_heat_flux",
        "tile_moisture_flux",
        "tile_temperature",
        "screen_temperature",
        "tile_heat_flux",
        "tile_moisture_flux",
        "snowice_sublimation",
        "surf_heat_flux",
        "canopy_evap",
        "snowice_melt",
        "time_since_transition",
        "surf_ht_flux",
        "water_extraction",
        "lake_evap",
        "snomlt_surf_htf",
        "soil_evap",
        "soil_surf_ht_flux",
        "t1p5m",
        "q1p5m",
        "qcl1p5m",
        "rh1p5m",
        "t1p5m_ssi",
        "q1p5m_ssi",
        "qcl1p5m_ssi",
        "rh1p5m_ssi",
        "t1p5m_land_loc",
        "q1p5m_land_loc",
        "t1p5m_land",
        "q1p5m_land",
        "qcl1p5m_land",
        "rh1p5m_land",
        "t1p5m_surft",
        "q1p5m_surft",
        "latent_heat",
        "surf_sw_net",
        "surf_radnet",
        "surf_lw_up",
        "surf_lw_down",
        "sea_ice_temperature",
        "latent_heat",
        ]}

# Set loop inference rules - used to find loops by index name
Loop.set_loop_type_inference_rules({
    "l": {"variable": "l"},
    "n": {"variable": "n"}})


def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply paralleldo transformations
    to each loop in jules_imp_kernel_mod.
    '''
    for loop in psyir.walk(Loop):
        if not loop.ancestor(Loop) and loop.loop_type == 'n':
            # Check if the loop over 'n' contains a loop over 'l'
            loop_descendents = [descendent for descendent in loop.walk(Loop)
                                if descendent is not loop]
            if loop_descendents[0].loop_type == 'l':
                # Now check if there are any if statements in this loop
                if_statements = [
                    descendent for descendent in
                    loop_descendents[0].walk(IfBlock, depth=None)
                    if descendent is not loop
                    ]
                # There is only one loop like this so we can just skip the
                # transformation for it
                if len(if_statements) > 0:
                    pass

        # Otherwise, transform all other loops
        elif not loop.ancestor(Loop):
            try:
                omp_transform_par_do.apply(loop, options)

            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform:\n %s", err)

# Ignore loops setting these as order dependent:
#   land_field l ainfo%land_index sice_pts ainfo%sice_index
#   sea_pts ainfo%sea_inde ainfo%sice_pts_ncat
