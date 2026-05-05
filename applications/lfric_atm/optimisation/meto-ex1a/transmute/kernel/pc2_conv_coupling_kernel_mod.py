# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
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
from psyclone.transformations import (OMPLoopTrans, TransformationError)
from psyclone.psyir.nodes import Loop
from transmute_psytrans.transmute_functions import (
    set_pure_subroutines,
    match_lhs_assignments,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC
)

# Variables that appear on the left-hand side of assignments
# for which PSyclone dependency errors can be ignored
false_dep_vars = [
    "dt_conv_wth",
    "dmv_conv_wth",
    "dmcl_conv_wth",
    "dcfl_conv_wth",
    "dbcf_conv_wth"
]

# Arrays that appear on the left-hand side of assignments
# which trigger automatic array privatisation in PSyclone
private_arrays = [
    "p_forcing",
    "p_work",
    "t_forcing",
    "t_work",
    "qv_forcing",
    "cfl_forcing",
    "qv_work",
    "qcl_work",
    "cfl_work",
    "cff_work",
    "bcf_work",
    "t_incr",
    "qv_incr",
    "qcl_incr",
    "cfl_incr",
    "bcf_incr",
]


def trans(psyir):
    """
    Apply OpenMP Directives
    """

    # Declare subroutine "pc2_hom_conv" as pure to enable parallelisation
    # of the encompassing loop
    set_pure_subroutines(psyir, "pc2_hom_conv")

    # OMPLoopTrans supports automatic array privatisation
    omp_pardo_dyn = OMPLoopTrans(omp_schedule="dynamic",
                                 omp_directive="paralleldo")

    # Add parallel do directives to outer loops
    for loop in psyir.walk(Loop):
        if not loop.ancestor(Loop):
            # Check if any variables appear on the LHS of assignment
            # expressions that require "privatisation" or lead to false
            # dependency errors in the parallel loop transformation.
            options = {"privatise_arrays":
                       len(match_lhs_assignments(loop, private_arrays)) > 0}
            ignore_deps_vars = match_lhs_assignments(loop, false_dep_vars)
            if len(ignore_deps_vars) > 0:
                options["ignore_dependencies_for"] = ignore_deps_vars

            # The loop with privatised arrays also needs dynamic schedule
            if options["privatise_arrays"]:
                try:
                    omp_pardo_dyn.apply(loop, options)
                except (TransformationError, IndexError) as err:
                    logging.warning("OMPLoopTrans failed: %s", err)
            else:
                try:
                    OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop, options)
                except (TransformationError, IndexError) as err:
                    logging.warning("OMPLoopTrans failed: %s", err)
