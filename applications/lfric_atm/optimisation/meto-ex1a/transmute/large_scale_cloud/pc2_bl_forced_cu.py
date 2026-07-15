# -----------------------------------------------------------------------------
# (C) 2025 Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
Optimisation script that replaces existing OpenMP parallelisation with
PSyclone-generated directives to target loops over index i instead of
index j. Trip count of j loops is 1 in LFRic, which prevents parallel
execution.
"""

import logging
from psyclone.transformations import TransformationError
from psyclone.psyir.nodes import Loop, Routine
from transmute_psytrans.transmute_functions import (
    loop_replacement_of,
    get_outer_loops,
    OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC
)


def trans(psyir):
    """
    Apply OpenMP Directives
    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    """

    # Remove any loops relating to specified loop type
    for node in psyir.walk(Routine):
        loop_replacement_of(node, "j")

    # Identify outer loops
    outer_loops = [loop for loop in get_outer_loops(psyir)
                   if not loop.ancestor(Loop)]

    # Apply OpenMP parallel do directives and use workaround for
    # firstprivate variable issue; replicate dynamic and static
    # schedules of the original implementation
    
    for idx, loop in enumerate(outer_loops):
        if idx == 0:
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC.apply(loop.walk(Loop)[0])
            except (TransformationError, IndexError) as err:
                logging.warning("OMPParallelLoopTrans failed: %s", err)
        else:
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop.walk(Loop)[0])
            except (TransformationError, IndexError) as err:
                logging.warning("OMPParallelLoopTrans failed: %s", err)
