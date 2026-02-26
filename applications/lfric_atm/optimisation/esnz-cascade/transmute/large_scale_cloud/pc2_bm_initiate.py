##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
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
from psyclone.psyir.nodes import (Loop, CodeBlock)
from transmute_psytrans.transmute_functions import (
    set_pure_subroutines,
    get_outer_loops,
    mark_explicit_privates,
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

# Subroutines that need to be declared as "pure"
pure_subroutines = ["qsat", "qsat_mix", "qsat_wat", "qsat_wat_mix"]

# Variables that appear on the left-hand side of assignments
# or as call arguments for which PSyclone dependency errors
# can be ignored
false_dep_vars = [
    "qc_points",
    "idx",
    "tl_in",
    "p_theta_levels",
    "qsi_lay",
    "qsl_lay",
]


class CompilerDirective():
    """
    Custom compiler directive class to avoid an issue
    with fparser.two.Fortran2003.Directive that will
    be resolved in an upcoming fparser release.
    """
    def __init__(self, directive):
        self.directive = directive

    def tofortran(self):
        """
        Return directive with prefix
        """
        return "!DIR$ " + self.directive


def trans(psyir):
    """
    Apply OpenMP and Compiler Directives
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

    # Declare private symbols for the last loop nest explicitly,
    # PSyclone misses one
    mark_explicit_privates(outer_loops[2], private_variables)

    # Parallelise the second region and insert compiler directives
    # Add redundant variable initialisation to work around a known
    # PSyclone issue when using CCE
    try:
        if get_compiler() == 'cce':
            first_priv_red_init(outer_loops[2], ["i", "j", "k"])

        OMP_PARALLEL_REGION_TRANS.validate(outer_loops[2:3])
        OMP_PARALLEL_REGION_TRANS.apply(outer_loops[2])

        # Insert before OpenMP directives to avoid PSyclone errors
        if get_compiler() == "cce":
            for loop in outer_loops[2].walk(Loop)[3:5]:
                cblock = CodeBlock([CompilerDirective("NOFISSION")],
                                   CodeBlock.Structure.STATEMENT)
                insert_at = loop.parent.children.index(loop)
                loop.parent.children.insert(insert_at, cblock)

        for loop in outer_loops[2].walk(Loop)[13:16]:
            cblock = CodeBlock([CompilerDirective("IVDEP")],
                               CodeBlock.Structure.STATEMENT)
            insert_at = loop.parent.children.index(loop)
            loop.parent.children.insert(insert_at, cblock)

        for loop in outer_loops[2].walk(Loop)[2:7]:
            # Check if any eligible variables appear in subroutine
            # call arguments; these lead to false dependency errors
            # in the parallel loop transformation that can be
            # ignored
            ignore_deps_vars = match_call_args(loop, false_dep_vars)
            options = {}
            if len(ignore_deps_vars) > 0:
                options["ignore_dependencies_for"] = ignore_deps_vars
            OMP_DO_LOOP_TRANS_STATIC.apply(loop, options)

        for loop in outer_loops[2].walk(Loop)[8:13:2]:
            # Check if any eligible variables appear on the LHS of
            # assignment expressions to ignore false dependency errors
            ignore_deps_vars = match_lhs_assignments(loop, false_dep_vars)
            options = {}
            if len(ignore_deps_vars) > 0:
                options["ignore_dependencies_for"] = ignore_deps_vars

            OMP_DO_LOOP_TRANS_STATIC.apply(loop, options)

    except (TransformationError, IndexError) as err:
        logging.warning("Parallelisation of the 2nd region failed: %s", err)
