# -----------------------------------------------------------------------------
#  (C) Crown copyright 2025 Met Office. All rights reserved.
#  The file LICENCE, distributed with this code, contains details of the terms
#  under which the code may be used.
# -----------------------------------------------------------------------------

'''
PSyclone transformation script for adjoint_tests.
In particular, contains changes to allow for filtering of routines based on if
they should be executed in serial.
'''

from psyclone_tools import (redundant_computation_setval, colour_loops,
                            view_transformed_schedule)

from psyclone.transformations import (Dynamo0p3OMPLoopTrans,
                                      OMPParallelTrans)
from psyclone.psyGen import InvokeSchedule

# NOTE: Whilst gen_*_lookup_code kernels are called from psy-lite,
#       this optimisation script has no impact.
def kernel_requires_serial_execution(kernel_name: str) -> bool:
    '''
    Returns a boolean depending on if a kernel requires serial execution
    (no OpenMP), through comparing with a set list of kernel names.

    :param kernel_name: String name of the kernel code call.
    :type kernel_name: str
    '''

    # ====================================================
    # === Known kernels that require serial execution. ===
    # glob-like wildcard patterns to pattern match against kernel names.
    SERIAL_KERNELS = ["gen_*_lookup_code"]
    # ====================================================

    for serial in SERIAL_KERNELS:
        # Match a glob pattern
        if "*" in serial:
            glob_sep = serial.split("*")

            # Ensure `*` only occurs once in the pattern.
            assert len(glob_sep) == 2
            start, end = glob_sep[0], glob_sep[1]

            if kernel_name[:len(start)] == start and \
               kernel_name[-len(end):] == end:
                # Glob match found
                return True

        # Match an exact kernel name
        elif kernel_name == serial:
            return True

    return False


def openmp_parallelise_loops_adj(psyir):
    '''
    Custom routine for adding OpenMP loops: Additional filtering
    for any kernels that require serial execution.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    '''
    otrans = Dynamo0p3OMPLoopTrans()
    oregtrans = OMPParallelTrans()

    # Loop over all the InvokeSchedule in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        # Add OpenMP to loops unless they are over colours, should be run in
        # serial, or are null
        for loop in subroutine.loops():
            needs_serial_exec = any(kernel_requires_serial_execution(k.name)
                                    for k in loop.kernels())

            if not needs_serial_exec and \
               loop.loop_type not in ["colours", "null"]:
                oregtrans.apply(loop)
                otrans.apply(loop, options={"reprod": True})


def trans(psyir):
    '''
    Applies PSyclone colouring, OpenMP and redundant computation
    transformations.

    :param psyir: the PSyIR of the PSyIR-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    redundant_computation_setval(psyir)
    colour_loops(psyir)
    openmp_parallelise_loops_adj(psyir)
    view_transformed_schedule(psyir)
