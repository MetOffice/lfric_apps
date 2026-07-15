# -----------------------------------------------------------------------------
# (C) 2025 Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
Optimisation script that replaces existing OpenMP parallelisation with
PSyclone-generated directives to target loops over index i instead of
index j. Trip count of j loops is 1 in LFRic, which prevents parallel
execution. Private variables need to be declared explicitly as PSyclone
analysis currently misses a scalar variable that a subroutine modifies in
a parallel region. PSyclone thread safety checks need to be overridden;
the subroutines can be safely parallelised. Compiler directives used in
the original code are re-inserted for performance and consistency of output.
"""

import logging
from psyclone.transformations import TransformationError
from psyclone.psyir.nodes import (
    Loop,
    UnknownDirective,
    OMPParallelDirective
)
from transmute_psytrans.transmute_functions import (
    set_pure_subroutines,
    get_outer_loops,
    get_compiler,
    first_priv_red_init,
    match_lhs_assignments,
    match_call_args,
    OMP_PARALLEL_REGION_TRANS,
    OMP_DO_LOOP_TRANS_STATIC
)

# Variables in parallel region that need to be private
private_variables = [
    "alphal", "alx", "i", "j", "k", "km1", "kp1", "mux", "tmp",
    "frac_init", "kk", "kkm1", "kkp1", "qc", "qc_points", "qsl",
    "tlx", "qsi", "idx", "deltacl_c", "deltacf_c", "deltaql_c",
    "cf_c", "cfl_c", "cff_c"
]

private_variable_par_sec = [
    "qsl", "qsi", "qsl_lay", "qsi_lay", "sigx", "deltacl_c", "idx",
    "deltacf_c", "cf_c", "cfl_c", "cff_c", "deltaql_c", "qnx_min",
    "qnx_max", "tl_lay", "ql_lay", "qsl_lay", "qsi_lay", "tdc_lay",
    "inv_thm_lay", "inv_tmp_lay", "wvar_lay", "dtldz_lay", "dqtdz_lay",
    "l_set_modes"]

# Subroutines that need to be declared as "pure"
pure_subroutines = ["qsat", "qsat_mix", "qsat_wat", "qsat_wat_mix"]

# Variables that appear on the left-hand side of assignments
# or as call arguments for which PSyclone dependency errors
# can be ignored
ignore_dependencies_for = [
    "qc_points",
    "idx",
    "tl_in",
    "p_theta_levels",
    "qsi_lay",
    "qsl_lay",
]


def trans(psyir):
    """
    Apply OpenMP and Compiler Directives
    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    """

    # Declare subroutines as pure to enable parallelisation
    # of the encompassing loops
    set_pure_subroutines(psyir, pure_subroutines)

    # Identify outer loops for setting up parallel regions
    outer_loops = [loop for loop in get_outer_loops(psyir)
                   if not loop.ancestor(Loop)]

    # Check if first OpenMP region can be parallelised and
    # apply directives
    try:
        OMP_PARALLEL_REGION_TRANS.validate(outer_loops[0:2])
        OMP_PARALLEL_REGION_TRANS.apply(outer_loops[0:2])
        OMP_DO_LOOP_TRANS_STATIC.apply(outer_loops[0])
        OMP_DO_LOOP_TRANS_STATIC.apply(outer_loops[1].walk(Loop)[1])
    except (TransformationError, IndexError) as err:
        logging.warning("Parallelisation of the 1st region failed: %s", err)

    # Parallelise the second region and insert compiler directives
    # Add redundant variable initialisation to work around a known
    # PSyclone issue when using CCE
    try:
        OMP_PARALLEL_REGION_TRANS.validate(outer_loops[2:3])
        OMP_PARALLEL_REGION_TRANS.apply(
            [outer_loops[2]], 
            force_private=private_variable_par_sec)
    except (TransformationError, IndexError) as err:
                logging.warning("Transformation failed as:", err)

    # Insert before OpenMP directives to avoid PSyclone errors
    if get_compiler() == "cce":
        for loop in outer_loops[2].walk(Loop)[3:5]:
            dir = UnknownDirective(" NOFISSION", "DIR")
            insert_at = loop.parent.children.index(loop)
            loop.parent.children.insert(insert_at, dir)

    for loop in outer_loops[2].walk(Loop)[13:18]:
        dir = UnknownDirective(" IVDEP", "DIR")
        insert_at = loop.parent.children.index(loop)
        loop.parent.children.insert(insert_at, dir)

    for loop in outer_loops[2].walk(Loop)[2:7]:
        if loop.ancestor(OMPParallelDirective):
            try:
                OMP_DO_LOOP_TRANS_STATIC.apply(
                    loop,
                    ignore_dependencies_for=ignore_dependencies_for,
                    force_private=private_variables)
            except (TransformationError, IndexError) as err:
                logging.warning("Transformation failed as:", err)

    for loop in outer_loops[2].walk(Loop)[8:13:2]:
        if loop.ancestor(OMPParallelDirective):
            try:
                OMP_DO_LOOP_TRANS_STATIC.apply(
                    loop, 
                    ignore_dependencies_for=ignore_dependencies_for,
                    force_private=private_variables)
            except (TransformationError, IndexError) as err:
                logging.warning("Transformation failed as:", err)
