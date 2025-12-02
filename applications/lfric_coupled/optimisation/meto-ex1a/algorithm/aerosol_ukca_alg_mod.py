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


def trans(psy):
    '''
    Applies PSyclone colouring and redundant computation transformations.
    Note that PSyclone currently cannot apply OpenMP transformations
    to UKCA because it is not thread-safe.

    '''
    redundant_computation_setval(psy)
    colour_loops(psy)
    view_transformed_schedule(psy)
