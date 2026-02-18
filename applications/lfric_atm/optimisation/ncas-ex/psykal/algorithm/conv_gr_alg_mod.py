# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
PSyclone transformation script for the LFRic (Dynamo0p3) API to apply
colouring and redundant computation to the level-1 halo for
the initialisation built-ins generically.
"""

from psyclone_tools import redundant_computation_setval, colour_loops, view_transformed_schedule


def trans(psyir):
    """
    Applies PSyclone colouring and redundant computation transformations.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    redundant_computation_setval(psyir)
    colour_loops(psyir)
    view_transformed_schedule(psyir)
