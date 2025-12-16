##############################################################################
# Copyright (c) 2025,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Bespoke Opt script for bl_diags_mod to add OpenMP to loops present in
the file. - To be replaced by a global in #900
'''

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import Loop
from transmute_psytrans.transmute_functions import (
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
)


def trans(psyir):
    '''
    Work through each loop in bl_diags_mod and add a parallel region
    '''
    for loop in psyir.walk(Loop):
        if not loop.ancestor(Loop):
            options = {
                "node-type-check": False}
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop, options)

            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)
