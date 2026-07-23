# -----------------------------------------------------------------------------
# (C) 2026 Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Bespoke PSyclone transformation script for sw_rad_tile_kernel_mod.
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
    get_children,
    are_variables_present,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
)


ignore_dependencies_for = [
    "albedo_obs_scaling", "tile_sw_direct_albedo",
    "tile_sw_diffuse_albedo", "sea_ice_pensolar_frac_direct",
    "sea_ice_pensolar_frac_diffuse", "flandg", "weight_blue"
    ]

do_not_omp_over = ["ainfo"]


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

        # If there is an loop over n, which has a child loop of i or l, skip
        if loop.variable.name == 'n':
            children = get_children(loop, node_type=Loop)
            if children:
                if children[0].variable.name in ['i', 'l']:
                    continue

        # If there is not an ancestor node which is a loop
        # Most ideal loops in this file are top loops
        # However, loops over i, l and n are okay, where the attempt
        # to paralellise the outer loop has failed or was skipped.
        # Loops over the i_band are skipped due to low iterations
        if ((not loop.ancestor(Loop) or loop.variable.name in ['i', 'l', 'n'])
                and loop.variable.name not in ['i_band']):
            # A loop over ainfo%sice_pts_ncat should not be parallised
            if loop.variable.name == 'n':
                if are_variables_present(loop, do_not_omp_over):
                    continue
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
                    loop,
                    ignore_dependencies_for=ignore_dependencies_for)

            except (TransformationError, IndexError) as err:
                logging.warning(f"Could not transform because:{err}")
