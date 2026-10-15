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
from psyclone.psyir.transformations import (
    ArrayAssignment2LoopsTrans,
    OMPLoopTrans,
    OMPMinimiseSyncTrans,
    TransformationError,
    MaximalOMPParallelRegionTrans,
)
from psyclone.psyir.nodes import (
    Assignment,
    Routine,
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
    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    '''

    loop_trans = OMPLoopTrans()
    minsync_trans = OMPMinimiseSyncTrans()

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

    # First convert assignments to loops whenever possible
    for assignment in psyir.walk(Assignment):
        try:
            ArrayAssignment2LoopsTrans().apply(assignment)
        except TransformationError:
            pass

    # Apply loop_trans to all the loops possible.
    for loop in psyir.walk(Loop):
        if loop.ancestor(OMPDoDirective) is not None:
            continue
        if loop.variable.name in ['i', 'l']:
            try:
                loop_trans.apply(
                    loop,
                    ignore_dependencies_for=ignore_dependencies_for,
                    nowait=True)
            except (TransformationError, IndexError) as err:
                logging.warning(
                    f"{fortran_file_name} Could not transform \
                    because:\n {err}")

    # Apply the largest possible parallel regions and remove any barriers that
    # can be removed.
    for routine in psyir.walk(Routine):
        MaximalOMPParallelRegionTrans().apply(routine)
        minsync_trans.apply(routine)
