# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
Optimisation script that adds OpenMP worksharing-loop directives to speed up
various loops. Some PSyclone dependency errors need to be overridden; these
assignments can be safely parallelised.
"""

import logging
from psyclone.transformations import TransformationError
from psyclone.psyir.nodes import Loop
from transmute_psytrans.transmute_functions import (
    match_lhs_assignments,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC
)

# Variables that appear on the left-hand side of assignments
# for which PSyclone dependency errors can be ignored
false_dep_vars = [
    "dtheta_inc_wth",
    "dmv_inc_wth",
    "dmcl_inc_wth",
    "dmci_inc_wth",
    "dms_inc_wth",
    "dcfl_inc_wth",
    "dcff_inc_wth",
    "dbcf_inc_wth",
    "sskew_bm",
    "svar_bm",
    "svar_tb"
]


def trans(psyir):
    """
    Apply OpenMP Directives
    """

    # Add parallel do directives to outer loops
    for loop in psyir.walk(Loop):
        if not loop.ancestor(Loop):
            # Check if any eligible variables appear on the LHS of
            # assignment expressions; these lead to false dependency
            # errors in the parallel loop transformation that can be
            # ignored
            ignore_deps_vars = match_lhs_assignments(loop, false_dep_vars)
            options = {}
            if len(ignore_deps_vars) > 0:
                options["ignore_dependencies_for"] = ignore_deps_vars

            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop, options)
            except (TransformationError, IndexError) as err:
                logging.warning("OMPParallelLoopTrans failed: %s", err)
