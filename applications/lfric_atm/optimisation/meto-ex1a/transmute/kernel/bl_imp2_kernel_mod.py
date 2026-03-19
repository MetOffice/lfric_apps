# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Custom script for bl_imp2_kernel_mod, where we cannot effectively add OMP
around all of the outer loops as it causes KGO issues with full and fast debug.
Instead if it is placed around the i loop, over the seg len loop range, we
maximise the parallelism, whilst preserving KGOs.
'''

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import Loop
from transmute_psytrans.transmute_functions import (
    match_lhs_assignments,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC
)


def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply paralleldo transformations
    to each loop.
    '''

    # Variables that appear on the left-hand side of assignments
    # for which PSyclone dependency errors can be ignored
    false_dep_vars = [
        "fric_heating_blyr",
        "fric_heating_incv",
        "nblyr",
        "cca",
        "ccw",
        "z_lcl",
        "cf_area",
        "cf_bulk",
        "cf_ice",
        "cf_liq",
        "dtheta_bl",
        "m_v",
        "m_cl",
        "m_s",
        "fqw_star_w3",
        "ftl_star_w3",
        "heat_flux_bl",
        "moist_flux_bl",
        "rhokh_bl",
        "dqw_wth",
        "dtl_wth",
        "qw_wth",
        "tl_wth"
        ]

    # Work through each loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):
        if loop.variable.name == 'i':
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
                logging.warning(
                    "Could not transform because:\n %s", err)
