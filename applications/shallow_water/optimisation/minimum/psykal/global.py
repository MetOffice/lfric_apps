##############################################################################
# (c) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''
PSyclone transformation script for the LFRic (Dynamo0p3) API to apply
the minumum set of transformations required to permit the application
to run safely in debug mode on all supported platforms.

Applying the 'redundant_computation_setval' transformation for the
initialisation built-ins permits developers to set safe values in halos.

'''
from psyclone_tools import (redundant_computation_setval,
                            view_transformed_schedule)


def trans(psyir):
    '''
    Applies PSyclone redundant computation transformations on
    initialisation built-ins only.

    :param psyir: the PSyIR of the PSyIR-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    view_transformed_schedule(psyir)
