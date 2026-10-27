##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
PSyclone transformation script for the LFRic (Dynamo0p3) API to apply
colouring and redundant computation to the level-1 halo for
the initialisation built-ins generically.

The inner product used by the checksum algorithm produces different
results on different numbers of threads if OpenMP is applied. It means
the checksums are different even if the model result is the
same. Therefore, we do not apply the OpenMP transform to the
algorithm.

'''
from psyclone_tools import (redundant_computation_setval, colour_loops,
                            view_transformed_schedule)


def trans(psyir):
    '''
    Applies PSyclone colouring and redundant computation transformations.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    colour_loops(psyir)
    view_transformed_schedule(psyir)
