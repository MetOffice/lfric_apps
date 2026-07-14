# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
A local.py script for all kernels, where instead of adding OMP across the
outermost loop, it is placed around the i loop, or across the l loop.
This script imports a SCRIPT_OPTIONS_DICT which can be used to override
small aspects of this script per file it is applied to.
Overrides currently include:
* ignore_dependencies_for
* node_type_check
* safe_pure_calls
'''

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import (
    Loop, Call,
    OMPParallelDoDirective,
    OMPParallelDirective,
    OMPDoDirective,)
from transmute_psytrans.transmute_functions import (
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
    safe_pure_calls = []

    if fortran_file_name in SCRIPT_OPTIONS_DICT:
        file_overrides = SCRIPT_OPTIONS_DICT[fortran_file_name]
        if "ignore_dependencies_for" in file_overrides.keys():
            ignore_dependencies_for = file_overrides[
                    "ignore_dependencies_for"]
        if "node_type_check" in file_overrides.keys():
            node_type_check = file_overrides[
                    "node_type_check"]
        if "safe_pure_calls" in file_overrides.keys():
            safe_pure_calls = file_overrides[
                    "safe_pure_calls"]
        
    # Set the calls to 'pure', given the provided override.
    # pure allows PSyclone to parallelise over them with OMP.
    if safe_pure_calls:
        for call in psyir.walk(Call):
            if call.routine.symbol.name in safe_pure_calls:
                call.routine.symbol.is_pure = True

    # Work through each loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):
        # If there is an OMP ancestor skip.
        if (
            loop.ancestor(OMPParallelDoDirective) is not None
            or loop.ancestor(OMPDoDirective) is not None
            or loop.ancestor(OMPParallelDirective) is not None
        ):
            continue
        # Allow loops over 'i' and 'l' indexes to be parallelised.
        if loop.variable.name in ['i', 'l']:
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
                    loop, 
                    ignore_dependencies_for=ignore_dependencies_for,
                    node_type_check=node_type_check)
            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)
