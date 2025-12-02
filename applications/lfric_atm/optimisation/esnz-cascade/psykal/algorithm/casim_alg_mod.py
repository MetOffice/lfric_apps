##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
PSyclone transformation script for the LFRic (Dynamo0p3) API to apply
colouring and redundant computation to the level-1 halo for
the initialisation built-ins generically.

'''

from psyclone_tools import (redundant_computation_setval, colour_loops,
                            view_transformed_schedule)
from psyclone.transformations import (Dynamo0p3OMPLoopTrans,
                                      OMPParallelTrans)
from psyclone.psyGen import InvokeSchedule


def trans(psyir):
    '''
    Applies PSyclone colouring and redundant computation transformations.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    colour_loops(psyir)

    # Extracted from psyclone tools#
    # openmp_parallelise_loops#
    # To avoid adding OMP around loops in invoke_1_casim_kernel_type

    # But otherwise still add openmp accross the algortithm and other calls
    # Currently there are no means in the psyclone tools to
    # avoid applying a transformation around something specific
    # such as this where we want to avoid the subroutine
    # casim_kernel_type or as psyclone knows is in it's representaion
    # invoke_1_casim_kernel_type

    otrans = Dynamo0p3OMPLoopTrans()
    oregtrans = OMPParallelTrans()

    # Loop over all the InvokeSchedule in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        if subroutine.invoke != psyir.invokes.get("invoke_casim_kernel_type"):
            # Add OpenMP to loops unless they are over colours or are null
            for loop in subroutine.loops():
                if loop.loop_type not in ["colours", "null"]:
                    oregtrans.apply(loop)
                    otrans.apply(loop, options={"reprod": True})
    # Extracted from psyclone tools#

    view_transformed_schedule(psyir)
