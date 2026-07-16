# -----------------------------------------------------------------------------
# (C) 2026 Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Local script for Boundary layer.
This by default removes the j loop(s) and replaces it with a init,
spans a PARALLEL region across the whole file,
and then adds OMP DO to the top most loop in each group of loops.
There are some small bespoke needs for the 6x files that this currently
affects which are captured below as they are initialised.

This is currently used by the following files:
* bl_lsp
* btq_int
* ex_flux_tq
* ex_flux_uv
* kmkh
* tr_mix
* imp_mix
* fm_drag
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
    Directive,
    Loop,
    Routine,
    OMPDoDirective,
)
from transmute_psytrans.transmute_functions import (
    loop_replacement_of,
    replace_n_threads,
)
from script_options import (
    SCRIPT_OPTIONS_DICT
)
def trans(psyir):
    """
    Local.py script for boundary layer.
    This spans a PARALLEL section across the whole file,
    and then adds OMP to either to top loop of a nest, or k

    This is currently used by the following files:
    * bdy_impl4
    * bl_lsp
    * btq_int
    * ex_flux_tq
    * ex_flux_uv
    * kmkh
    * tr_mix
    * imp_mix
    * fm_drag

    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """

    loop_trans = OMPLoopTrans()
    minsync_trans = OMPMinimiseSyncTrans()

    ignore_dependencies_for = []
    max_threads_parse = False

    # Get the file name to use with the SCRIPT_OPTIONS_DICT
    fortran_file_name = str(psyir.root.name)
    # Check if file is in the script_options_dict
    # Copy out anything that's needed
    # Only the options list is currently
    if fortran_file_name in SCRIPT_OPTIONS_DICT:
        file_overrides = SCRIPT_OPTIONS_DICT[fortran_file_name]
        # Update the respective lists if the filename override exists
        if "ignore_dependencies_for" in file_overrides.keys():
            ignore_dependencies_for = file_overrides["ignore_dependencies_for"]
        if "max_threads_parse" in file_overrides.keys():
            max_threads_parse = file_overrides["max_threads_parse"]

    # Replace max_threads = 1
    if max_threads_parse:
        replace_n_threads(psyir, "max_threads")

    # Remove any loops relating to specified loop type
    for node in psyir.walk(Routine):
        loop_replacement_of(node, "j")

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
        #if not loop.ancestor(Directive):
        if loop.variable.name in ['i', 'ii', 'l']:
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

