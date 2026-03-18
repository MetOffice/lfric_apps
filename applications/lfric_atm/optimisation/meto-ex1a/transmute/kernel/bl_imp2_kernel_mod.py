##############################################################################
# Copyright (c) 2025,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
A global script to add OpenMP to loops present in the file provided.
This script imports a SCRIPT_OPTIONS_DICT which can be used to override
small aspects of this script per file it is applied to.
Overrides currently include:
* Options list for transformations
* safe pure calls for loops over calls which can be parallelised
'''

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import Loop, Literal
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
        if (loop.variable.name == 'i'):
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
