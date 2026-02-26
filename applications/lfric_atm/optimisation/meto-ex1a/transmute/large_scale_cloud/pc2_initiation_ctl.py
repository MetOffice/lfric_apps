##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Optimisation script that replaces existing OpenMP parallelisation with
PSyclone-generated directives to parallelise additional loops.
"""

import logging
from psyclone.transformations import TransformationError
from psyclone.psyir.nodes import Loop
from transmute_psytrans.transmute_functions import (
    get_outer_loops,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
)


def trans(psyir):
    """
    Apply OpenMP Directives
    """

    # Identify outer loops in the subroutine
    outer_loops = [loop for loop in get_outer_loops(psyir)
                   if not loop.ancestor(Loop)]

    try:
        # Parallelise k-loops and i-loops (j-loops have a trip count of 1)
        for loop in outer_loops:
            if loop.variable.name == 'k':
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop)
            elif loop.variable.name == 'j':
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop.walk(Loop)[1])
    except (TransformationError, IndexError) as err:
        logging.warning("OMPParallelLoopTrans failed: %s", err)
