##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
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


def trans(psy):
    '''
    Applies PSyclone colouring and redundant computation transformations.
    '''
    redundant_computation_setval(psy)
    colour_loops(psy)

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
    # Loop over all the Invokes in the PSy object
    for invoke in psy.invokes.invoke_list:
        if invoke != psy.invokes.get("invoke_casim_kernel_type"):
            schedule = invoke.schedule

            # Add OpenMP to loops unless they are over colours or are null
            for loop in schedule.loops():
                if loop.loop_type not in ["colours", "null"]:
                    oregtrans.apply(loop)
                    otrans.apply(loop, options={"reprod": True})
    # Extracted from psyclone tools#

    view_transformed_schedule(psy)
