# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
A local.py script for all kernels, where instead of adding OMP across the
outermost loop, it is placed around the i loop, over the seg len loop range,
or across the l loop.
'''

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import Loop
from transmute_psytrans.transmute_functions import (
    match_lhs_assignments,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC
)
from script_options import (
    SCRIPT_OPTIONS_DICT
)


def trans(psyir):
    '''
    PSyclone function call, run through psyir object,
    each schedule (or subroutine) and apply paralleldo transformations
    to each loop.
    '''

    fortran_file_name = str(psyir.root.name)

    false_dep_vars = []

    if fortran_file_name in SCRIPT_OPTIONS_DICT:
        file_overrides = SCRIPT_OPTIONS_DICT[fortran_file_name]
        if "false_dep_vars" in file_overrides.keys():
            false_dep_vars = file_overrides["false_dep_vars"]

    # Work through each loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):
        if loop.variable.name in ['i', 'l']:
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
