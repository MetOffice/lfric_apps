# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Bespoke script for kmkhz_9c_wtrac.
Remove any j loops.
Place a OMP parallel do inside the i_wt loop.
'''

import logging
from psyclone.psyir.nodes import (
    Schedule,
    Routine,
    Loop,
    OMPParallelDirective,
)
from psyclone.transformations import TransformationError
from transmute_psytrans.transmute_functions import (
    loop_replacement_of,
    get_compiler,
    first_priv_red_init,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
)

Loop.set_loop_type_inference_rules({
    "i_wt": {"variable": "i_wt"},
    "j": {"variable": "j"}})


def trans(psyir):
    """
    Bespoke script for kmkhz_9c_wtrac
    """

    # For the PSyclone 3.1 bug with CCE. Certain files
    # are causing a first private to be generated in the parallel section.
    first_private_list = [
        "dz_disc",
        "qw_lapse",
        "k",
        "fa_tend",
        "inv_tend",
        "ml_tend",
        "totqf_efl"
        ]

    # Remove any loops relating to j loop type
    for node in psyir.walk(Routine):
        loop_replacement_of(node, "j")

    # Span a parallel section across the whole routine,
    # apart for a few exceptions provided
    for loop in psyir.walk(Loop):
        loop_ancestor_type = ""
        try:
            loop_ancestor_type = loop.ancestor(Loop).loop_type
        # pylint: disable=bare-except
        except:  # noqa: E722
            pass
        if loop_ancestor_type:
            if loop_ancestor_type == "i_wt":
                # Span the region
                try:
                    OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
                        loop)
                except (TransformationError, IndexError) as err:
                    logging.warning("OMPParallelLoopTrans failed: %s", err)

    # CCE first private issue with 3.1, to be removed longer term
    # pylint: disable=too-many-nested-blocks
    if get_compiler() == "cce":
        for routine in psyir.walk(Routine):
            for node in routine.children:
                for schedule in node.walk(Schedule):
                    for child in schedule.children:
                        if isinstance(child, OMPParallelDirective):
                            first_priv_red_init(child, first_private_list)
                            break
