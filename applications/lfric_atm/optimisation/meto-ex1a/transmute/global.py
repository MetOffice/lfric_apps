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
from psyclone.psyir.nodes import Loop
from transmute_psytrans.transmute_functions import (
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
    set_pure_subroutines,
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

    # options list for transformation.
    options = {}

    # Designate calls in regions as safe to parallelise over.
    safe_pure_calls = []

    fortran_file_name = str(psyir.root.name)
    # Check if file is in the script_options_dict
    # Copy out anything that's needed
    # options list and a pure calls override
    if fortran_file_name in SCRIPT_OPTIONS_DICT:
        file_overrides = SCRIPT_OPTIONS_DICT[fortran_file_name]
        if "options" in file_overrides.keys():
            options = file_overrides["options"]
        if "safe_pure_calls" in file_overrides.keys():
            safe_pure_calls = file_overrides["safe_pure_calls"]

    # Set the pure calls if needed
    if safe_pure_calls:
        set_pure_subroutines(psyir, safe_pure_calls)

    # Work through each loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):
        if not loop.ancestor(Loop):
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop, options)

            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)
