# -----------------------------------------------------------------------------
# (C) 2026 Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
A global script to add OpenMP to loops present in the file provided.
This script imports a SCRIPT_OPTIONS_DICT which can be used to override
small aspects of this script per file it is applied to.
Overrides currently include:
* ignore_dependencies_for
* node_type_check
'''

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import (
    Loop,
    OMPParallelDoDirective,
    OMPParallelDirective,
    OMPDoDirective,)
from transmute_psytrans.transmute_functions import (
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC
)
from script_options import (
    SCRIPT_OPTIONS_DICT
)

ignore_dependencies_for = [
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
        "ainfo%sice_pts_ncat", "ainfo%sice_index_ncat", "flandg",
    ]

def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply paralleldo transformations
    to each loop.
    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    '''

    # Work through each loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):
        # If there is an OMP ancestor skip.
        if (
            loop.ancestor(OMPParallelDoDirective) is not None
            or loop.ancestor(OMPDoDirective) is not None
            or loop.ancestor(OMPParallelDirective) is not None
        ):
            continue
        # If there is not an ancestor node which is a loop
        # Most ideal loops in this file are top loops
        if not loop.ancestor(Loop) or loop.variable.name == 'n':
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
                    loop, 
                    ignore_dependencies_for=ignore_dependencies_for)

            except (TransformationError, IndexError) as err:
                logging.warning(f"Could not transform because:{err}")
