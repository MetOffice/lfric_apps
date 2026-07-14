# -----------------------------------------------------------------------------
# (C) 2026 Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
A global script to add OpenMP to loops present in the file provided.
This script imports a SCRIPT_OPTIONS_DICT which can be used to override
small aspects of this script per file it is applied to.
Overrides currently include:
* ignore_dependencies_for
* node_type_check
'''

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import Loop
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
    :returns:
    '''

    fortran_file_name = str(psyir.root.name)

    node_type_check = True
    ignore_dependencies_for = []

    fortran_file_name = str(psyir.root.name)
    # Check if file is in the script_options_dict
    # Copy out anything that's needed
    # options list and a pure calls override
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
        if not loop.ancestor(Loop):
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
                    loop, 
                    ignore_dependencies_for=ignore_dependencies_for,
                    node_type_check=node_type_check)

            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)
