# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
PSyclone Script for applying OpenMP transformations to sw kernel mod
"""

import logging
from psyclone.transformations import (
    TransformationError)
from psyclone.psyir.nodes import (
    Call, Loop)
from transmute_psytrans.transmute_functions import (
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
    replace_n_threads,
    first_priv_red_init)

# Associate ll index with segmentation_ll tool type
Loop.set_loop_type_inference_rules({"segmentation_ll": {"variable": "ll"}})


def trans(psyir):
    """
    Entry point for OpenMP transformations for `sw_kernel_mod.F90`.
    """

    # Replace max_threads = 1
    replace_n_threads(psyir, "max_threads")

    # Enable pure for runes call to allow PSyclone to add OMP around it
    for call in psyir.walk(Call):
        if call.routine.name == "runes":
            call.routine.symbol.is_pure = True

    # Walk the loops of the psyir obj
    for loop in psyir.walk(Loop):
        # if the loop is of the set type above, segmentation_ll
        if loop.loop_type == "segmentation_ll":
            # first private workaround
            first_priv_red_init(loop, ["n_profile_list_seg", "seg_end"])
            # Apply the transformation
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop)
            except (TransformationError, IndexError) as err:
                logging.warning(
                    "Could not transform because:\n %s", err)
