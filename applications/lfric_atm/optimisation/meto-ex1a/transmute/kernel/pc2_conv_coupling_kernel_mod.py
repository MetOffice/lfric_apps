# -----------------------------------------------------------------------------
# (C) 2026 Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
Optimisation script that adds OpenMP parallel do worksharing-loop directives.
The main loop requires dynamic schedule to improve load balancing between
threads. Some PSyclone dependency errors and a subroutine thread safety
check need to be overridden; these assignments and subroutine call can be
safely parallelised. Multiple arrays need to be declared OpenMP-private.
"""

import logging
from psyclone.transformations import TransformationError
from psyclone.psyir.nodes import (
    Loop,
    OMPParallelDoDirective,
    OMPParallelDirective,
    OMPDoDirective,
)
from transmute_psytrans.transmute_functions import (
    set_pure_subroutines,
    OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
)

# Variables that appear on the left-hand side of assignments
# for which PSyclone dependency errors can be ignored
ignore_dependencies_for = [
    "bcf_incr", "bcf_work", "cff_work", "cfl_forcing", "cfl_incr",
    "cfl_work", "p_forcing", "p_work", "qcl_incr", "qcl_work",
    "qv_forcing", "qv_incr", "qv_work", "t_forcing", "t_incr",
    "t_work", "zeros", "dt_conv_wth", "dmv_conv_wth", "dmcl_conv_wth",
    "dcfl_conv_wth", "dbcf_conv_wth"
]

# Arrays that appear on the left-hand side of assignments
# which trigger automatic array privatisation in PSyclone
force_privates = [
    "p_forcing", "p_work", "t_forcing", "t_work", "qv_forcing",
    "cfl_forcing", "qv_work", "qcl_work", "cfl_work", "cff_work",
    "bcf_work", "t_incr", "qv_incr", "qcl_incr", "cfl_incr", "bcf_incr",
]


def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply paralleldo transformations
    to each loop.
    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    '''

    # Declare subroutine "pc2_hom_conv" as pure to enable parallelisation
    # of the encompassing loop
    set_pure_subroutines(psyir, "pc2_hom_conv")

    # Work through each loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):
        # If there is an OMP ancestor skip.
        if (
            loop.ancestor(OMPParallelDoDirective) is not None
            or loop.ancestor(OMPDoDirective) is not None
            or loop.ancestor(OMPParallelDirective) is not None
        ):
            continue

        if not loop.ancestor(Loop):
            # We wish to parallelise over the outmost loops.
            # The k loop, also over the call down to pc2_hom_conv, we wish to
            # parallelise. As each level is likely load balanced differently,
            # a dynamic OMP schedule is better.
            # All other loops in this file are single dimension i loops.
            if loop.variable.name == 'k':
                try:
                    OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC.apply(
                        loop,
                        force_private=force_privates,
                        ignore_dependencies_for=ignore_dependencies_for,
                        node_type_check=False)
                except (TransformationError, IndexError) as err:
                    logging.warning(f"Could not transform because:{err}")
            else:
                try:
                    OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
                        loop,
                        ignore_dependencies_for=ignore_dependencies_for)
                except (TransformationError, IndexError) as err:
                    logging.warning(f"Could not transform because:{err}")
