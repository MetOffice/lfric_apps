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

    node_type_check = True
    ignore_dependencies_for = []

    if fortran_file_name in SCRIPT_OPTIONS_DICT:
        file_overrides = SCRIPT_OPTIONS_DICT[fortran_file_name]
        if "ignore_dependencies_for" in file_overrides.keys():
            ignore_dependencies_for = file_overrides[
                    "ignore_dependencies_for"]
        if "node_type_check" in file_overrides.keys():
            node_type_check = file_overrides[
                    "node_type_check"]

    # Work through each loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):
        if loop.variable.name in ['i', 'l']:

            try:
                #OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop, options)
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
                    loop, 
                    ignore_dependencies_for=ignore_dependencies_for,
                    node_type_check=node_type_check)

            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)
                print(f"Could not transform because:\n {err}")
