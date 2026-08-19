# -----------------------------------------------------------------------------
#  (C) Crown copyright 2023 Met Office. All rights reserved.
#  The file LICENCE, distributed with this code, contains details of the terms
#  under which the code may be used.
# -----------------------------------------------------------------------------


"""
This file contains frequently used transformations to simplify
their application in PSyclone optimisations scripts.

"""

from psyclone.domain.lfric import LFRicConstants
from psyclone.psyGen import InvokeSchedule
from psyclone.psyir.nodes import (
    Loop,
    Routine,
    Directive,
    Container,
    OMPParallelDirective,
    OMPParallelDoDirective,
    OMPDoDirective,
    FileContainer,
    ProfileNode
)
try:
    from psyclone.psyir.transformations import OMPParallelTrans
except ImportError:
    # Support for psyclone < 3.3
    from psyclone.transformations import OMPParallelTrans
try:
    from psyclone.domain.lfric.transformations import (
        LFRicRedundantComputationTrans)
except ImportError:
    # Support for psyclone < 3.3
    from psyclone.transformations import LFRicRedundantComputationTrans

from psyclone.transformations import (
    LFRicColourTrans,
    LFRicOMPLoopTrans,
    TransformationError
)
from psyclone.psyir.transformations import ProfileTrans

# List of allowed 'setval_*' built-ins for redundant computation transformation
SETVAL_BUILTINS = ["setval_c"]


# -----------------------------------------------------------------------------
def redundant_computation_setval(psyir: FileContainer):
    """
    Applies the redundant computation transformation to loops over DoFs
    for the initialision built-ins, 'setval_*'.

    To reduce MPI communications, current PSyclone-LFRic strategy does not
    apply halo swaps on input arguments to kernels with increment
    operations on continuous fields such as 'GH_INC'. For such kernels,
    PSy-layer code needs to loop into the halo to correctly compute owned
    DoFs on the boundary between the halo and the domain. Therefore values
    of the remaining DoFs in the first halo cell need to be initialised to
    values that will not induce numerical errors.

    By default, the initialisation 'setval_*' built-ins do not initialise
    into the halos. This transform causes them to do so, and so permits
    developers to set safe values in halos.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    :raises Exception: if there is more than one built-in call per DoF loop.

    """
    # Import redundant computation transformation
    rtrans = LFRicRedundantComputationTrans()

    # Loop over all the InvokeSchedule in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        # Make setval_* built-ins compute redundantly to the level-1 halo
        # if they are in their own loop
        for loop in subroutine.loops():
            if loop.iteration_space == "dof":
                if len(loop.kernels()) != 1:
                    raise Exception(
                        f"Expecting loop to contain 1 call but found "
                        f"'{len(loop.kernels())}'"
                    )
                if loop.kernels()[0].name in SETVAL_BUILTINS:
                    rtrans.apply(loop, options={"depth": 1})


# -----------------------------------------------------------------------------
def colour_loops(psyir: FileContainer, enable_tiling=False):
    """
    Applies the colouring transformation to all applicable loops and optionally
    enables tiling.
    It creates the instance of `LFRicColourTrans` only once.

    :param psyir: the PSyIR of the PSy-layer.
    :param enable_tiling: a bool to enable tiling. Default False.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    const = LFRicConstants()
    ctrans = LFRicColourTrans()

    # Loop over all the subroutines in the PSyIR object
    for subroutine in psyir.walk(Routine):
        # Colour loops over cells unless they are on discontinuous
        # spaces or over DoFs
        for child in subroutine.children:
            # Check if the profiling calipers have been added before the 
            # colouring.
            if isinstance(child, ProfileNode):
                raise TransformationError(
                "Must apply colour_loops BEFORE profile_loops function "
                "in optimisation script.")
            if (
                isinstance(child, Loop)
                and child.iteration_space.endswith("cell_column")
                and child.field_space.orig_name
                not in const.VALID_DISCONTINUOUS_NAMES
            ):
                ctrans.apply(child, options={"tiling": enable_tiling})

# -----------------------------------------------------------------------------
def profile_loops(psyir: FileContainer, colours_only=True):
    """
    Applies timing calipers to kernels during the psyclone build. The default
    is to only profile coloured loops but colours_only can be set to False to
    profile every instance of a coded kernel.

    :param psyir: the PSyIR of the PSy-layer.
    :param colours_only: profile only the coloured kernels. Default True.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    profile_trans = ProfileTrans()
    leave_loops = ["cells_in_colour",
                   "tiles_in_colour",
                   "cells_in_tile"]

    # Loop over all the InvokeSchedule in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        # Add timing calipers to coloured loops. This should be done
        # before the application of the openmp transformation.
        count = 0
        for loop in subroutine.loops():
            if not loop.coded_kernels():
                continue
            # Insert profiler calls before loop over colours
            if ((loop.loop_type == "colours") or 
                (colours_only is False and loop.loop_type not in leave_loops)):
                # First check that the transformation is not being made inside 
                # an OMP region.
                if (loop.ancestor(OMPParallelDirective) 
                    or loop.ancestor(OMPParallelDoDirective)
                    or loop.ancestor(OMPDoDirective)):
                    raise TransformationError(
                        "Must apply profile_loops BEFORE "
                        "openmp_parallelise_loops function in optimisation "
                        "script.")
                # Constructing unique calliper name based on kernel name,
                # invoke name and kernel count
                k_object = loop.ancestor(InvokeSchedule).coded_kernels()[count]
                k_name = k_object.name
                invoke_name = loop.ancestor(InvokeSchedule).invoke.name
                file_name = loop.ancestor(Container).name
                # Make region name
                region_name = invoke_name + ":" + k_name + "_k"  + str(count)
                options = {"region_name": (file_name, region_name)}
                profile_trans.apply(loop, options=options)
                # Count here is to distinguish kernels of the same name
                # in the same invoke.
                count += 1

# -----------------------------------------------------------------------------
def openmp_parallelise_loops(psyir: FileContainer):
    """
    Applies OpenMP Loop transformation to each applicable loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    otrans = LFRicOMPLoopTrans()
    oregtrans = OMPParallelTrans()

    # Loop over all the InvokeSchedule in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        # Add OpenMP to loops unless they are over colours, are null,
        # or if an outer loop is already parallelised (OpenMP is applied
        # to loop over tiles instead of cells if tiling is enabled)
        for loop in subroutine.loops():
            if loop.loop_type not in ["colours", "null"] and \
               not loop.ancestor(Directive):
                oregtrans.apply(loop)
                otrans.apply(loop, options={"reprod": True})


# -----------------------------------------------------------------------------
def view_transformed_schedule(psyir: FileContainer):
    """
    Provides view of transformed Invoke schedule in the PSy-layer.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    setval_count = 0

    # Loop over all the Invokes in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        print(f"Transformed invoke '{subroutine.name}' ...")

        # Count instances of setval_* built-ins
        for loop in subroutine.loops():
            if loop.iteration_space == "dof":
                if loop.kernels()[0].name in SETVAL_BUILTINS:
                    setval_count += 1

        # Take a look at what we have done
        print(f"Found {setval_count} {SETVAL_BUILTINS} calls")
        print(subroutine.view())
