# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Custom script for bl_exp_kernel_mod, where we cannot effectively add OMP
around all of the outer loops.
A k loop has a low number of iterations, requiring the parallelism to be
pushed down onto the i loop.
A custom function, check_literals assists this.
'''

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import Loop, Literal
from transmute_psytrans.transmute_functions import (
    match_lhs_assignments,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC
)


def check_literals(node):
    '''
    Check if the bounds, found in the literals are a value
    or are a variable.
    '''
    literals = node.walk(Literal)
    # Position 2 in the array relates to the range of the loop node
    # This range, if a variable is used is 1
    # Therefore if a number is used, it will be above 1
    # Relative to this file, 3 is the only range in question
    if int(literals[1].value) > 1:
        # Return False if there is one
        return False
    # Return True if there isn't one, indicating an appropriate loop
    return True


def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply paralleldo transformations
    to each loop, adjusting according to loop ranges.
    '''

    # Variables that appear on the left-hand side of assignments
    # for which PSyclone dependency errors can be ignored
    false_dep_vars = [
            "land_index",
            "tnuc_nlcl",
            "dw_bl",
            "surf_interp",
            "bq_bl",
            "bt_bl",
            "dtrdz_tq_bl",
            "level_ent",
            "level_ent_dsc",
            "ent_t_frac",
            "ent_t_frac_dsc",
            "ent_we_lim",
            "ent_we_lim_dsc",
            "ent_zrzi",
            "ent_zrzi_dsc",
            "fd_taux",
            "fd_tauy",
            "gradrinr",
            "heat_flux_bl",
            "lmix_bl",
            "moist_flux_bl",
            "ngstress_bl",
            "rhokh_bl",
            "rhokm_bl",
            "tke_bl",
            "bl_weight_1dbl",
            "zh_nonloc",
            "z_lcl",
            "inv_depth",
            "qcl_at_inv_top",
            "shallow_flag",
            "uw0_flux",
            "vw0_flux",
            "lcl_height",
            "parcel_top",
            "level_parcel_top",
            "wstar_2d",
            "thv_flux",
            "parcel_buoyancy",
            "qsat_at_lcl",
            "bl_type_ind",
            "visc_m_blend",
            "visc_h_blend",
            "zh_2d",
            "zhsc_2d",
            "ntml_2d",
            "cumulus_2d",
            "rh_crit",
            "mix_len_bm",
            "dsldzm",
            "wvar",
            "zht",
            "oblen"
        ]

    # Work through each loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):

        # For the checks below, break this out into a flag
        safe_to_transform = False

        # save a loop ancestor reference per loop step
        ancestor = loop.ancestor(Loop)

        # There are three loop variables in this file, l, k or i
        if loop.variable.name == 'k':
            # check if a k loop bound is not a variable and not ideal to thread
            # over due to the low range, and there are no loop ancestors
            if check_literals(loop) and not ancestor:
                safe_to_transform = True

        elif loop.variable.name == 'i':
            # If there is an ancestor, we generally don't want to
            # unless that ancestor was not ideal to thread over
            if ancestor:
                if not check_literals(ancestor):
                    safe_to_transform = True
            else:
                safe_to_transform = True

        else:  # (loop.variable.name == 'l')
            if not ancestor:
                safe_to_transform = True

        if safe_to_transform:
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
