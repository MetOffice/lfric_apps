##############################################################################
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################


'''PSyclone transformation script for spt_main_alg_mod to apply colouring
and GPU offloading/CPU parallelization. Also adds redundant computation to
the level-1 halo for setval_* generically. This is based on
https://github.com/stfc/PSyclone/blob/master/examples/lfric/
scripts/gpu_offloading.py .

'''

import os
import sys
from psyclone.domain.lfric import LFRicConstants
from psyclone.psyir.nodes import Directive, Loop, Routine
from psyclone.psyir.transformations import (
    ACCKernelsTrans, TransformationError, OMPTargetTrans,
    OMPDeclareTargetTrans)
from psyclone.transformations import (
    Dynamo0p3ColourTrans, Dynamo0p3OMPLoopTrans,
    Dynamo0p3RedundantComputationTrans, OMPParallelTrans,
    ACCParallelTrans, ACCLoopTrans, ACCRoutineTrans,
    OMPLoopTrans)
from psyclone.domain.common.transformations import KernelModuleInlineTrans


# Names of any invoke that we won't add any GPU offloading
INVOKE_EXCLUSIONS = [
]

# Names of any kernel that we won't add parallelization
KERNEL_EXCLUSIONS = [
]

# Names of any kernels that we won't offload to GPU
GPU_KERNEL_EXCLUSIONS = ["spt_saturation_cap_code",]
# spt_saturation_cap_code: GPU transformation cannot be applied because of
# using qsat_wat_mix from qsat_mod. As qsat_mod is going to be modified in
# future, this falls out of the scope of the NGARCH project.
# Error message: Transformation Error: Kernel 'spt_saturation_cap_code'
# accesses the symbol 'qsat_wat_mix: RoutineSymbol<NoType, pure=unknown,'
# 'elemental=unknown>' which is imported.
# If this symbol represents data then it must first be converted to a
# Kernel argument using the KernelImportsToArguments transformation.

OFFLOAD_DIRECTIVES = os.getenv('LFRIC_OFFLOAD_DIRECTIVES', "none")


def trans(psyir):
    '''Applies PSyclone colouring and GPU offloading transformations. Any
    kernels that cannot be offloaded to GPU are parallelised using OpenMP
    on the CPU if they can be parallelised. Any setval_* kernels are
    transformed so as to compute into the L1 halos.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''
    inline_trans = KernelModuleInlineTrans()
    rtrans = Dynamo0p3RedundantComputationTrans()
    ctrans = Dynamo0p3ColourTrans()
    otrans = Dynamo0p3OMPLoopTrans()
    const = LFRicConstants()
    cpu_parallel = OMPParallelTrans()

    if OFFLOAD_DIRECTIVES == "omp":
        # Use OpenMP offloading
        loop_offloading_trans = OMPLoopTrans(
            omp_directive="teamsdistributeparalleldo",
            omp_schedule="none"
        )
        # OpenMP does not have a kernels parallelism directive equivalent
        # to OpenACC 'kernels'
        kernels_trans = None
        gpu_region_trans = OMPTargetTrans()
        gpu_annotation_trans = OMPDeclareTargetTrans()
    elif OFFLOAD_DIRECTIVES == "acc":
        # Use OpenACC offloading
        loop_offloading_trans = ACCLoopTrans()
        kernels_trans = ACCKernelsTrans()
        gpu_region_trans = ACCParallelTrans(default_present=False)
        gpu_annotation_trans = ACCRoutineTrans()
    elif OFFLOAD_DIRECTIVES == "none":
        pass
    else:
        print(
            f"The PSyclone transformation script expects the "
            f"LFRIC_OFFLOAD_DIRECTIVES to be set to 'omp' or 'acc' or "
            f"'none' but found '{OFFLOAD_DIRECTIVES}'."
        )
        sys.exit(-1)

    print(f"PSy name = '{psyir.name}'")

    for subroutine in psyir.walk(Routine):

        print("Transforming invoke '{0}' ...".format(subroutine.name))

        # Make setval_* compute redundantly to the level 1 halo if it
        # is in its own loop
        for loop in subroutine.loops():
            if loop.iteration_space == "dof":
                if len(loop.kernels()) == 1:
                    if loop.kernels()[0].name in ["setval_c"]:
                        rtrans.apply(loop, options={"depth": 1})

        if (
            psyir.name.lower() in INVOKE_EXCLUSIONS
            or OFFLOAD_DIRECTIVES == "none"
        ):
            print(
                f"Not adding GPU offloading to invoke '{subroutine.name}'"
            )
            offload = False
        else:
            offload = True

        # Keep a record of any kernels we fail and succeed to offload
        succeeded_offload = set()
        failed_to_offload = set()

        # Colour loops over cells unless they are on discontinuous spaces
        # (alternatively we could annotate the kernels with atomics)
        for loop in subroutine.loops():
            if loop.iteration_space.endswith("cell_column"):
                if (loop.field_space.orig_name not in
                        const.VALID_DISCONTINUOUS_NAMES):
                    ctrans.apply(loop)

        # Mark kernels inside the loops over cells as GPU-enabled
        # and inline them.
        for loop in subroutine.loops():
            if loop.iteration_space.endswith("cell_column"):
                if offload:
                    for kern in loop.kernels():
                        if kern.name.lower() in (
                            GPU_KERNEL_EXCLUSIONS + KERNEL_EXCLUSIONS +
                            list(succeeded_offload)
                        ):
                            continue

                        try:
                            gpu_annotation_trans.apply(
                                kern, options={'force': True}
                            )
                            print(f"GPU-annotated kernel '{kern.name}'")

                            try:
                                inline_trans.apply(kern)
                                print(f"Module-inlined kernel '{kern.name}'")
                                succeeded_offload.add(kern.name.lower())
                            except TransformationError as err:
                                print(
                                    f"Failed to module-inline '{kern.name}'"
                                    f" due to:\n{err.value}"
                                )
                        except TransformationError as err:
                            failed_to_offload.add(kern.name.lower())
                            print(
                                f"Failed to annotate '{kern.name}' with "
                                f"GPU-enabled directive due to:\n{err.value}"
                            )
                        # For annotated or inlined kernels we could attempt to
                        # provide compile-time dimensions for the temporary
                        # arrays and convert to code unsupported intrinsics.

        # Add GPU offloading to loops unless they are over colours or are null.
        for loop in subroutine.walk(Loop):
            kernel_names = [k.name.lower() for k in loop.kernels()]
            if offload and all(
                name not in (
                    list(failed_to_offload) + GPU_KERNEL_EXCLUSIONS +
                    KERNEL_EXCLUSIONS
                )
                for name in kernel_names
            ):
                try:
                    if loop.loop_type == "colours":
                        pass
                    if loop.loop_type == "colour":
                        loop_offloading_trans.apply(
                            loop, options={"independent": True})
                        gpu_region_trans.apply(loop.ancestor(Directive))
                    if loop.loop_type == "":
                        loop_offloading_trans.apply(
                            loop, options={"independent": True})
                        gpu_region_trans.apply(loop.ancestor(Directive))
                    if loop.loop_type == "dof":
                        # Loops over dofs can contains reductions
                        if kernels_trans:
                            # If kernel offloading is available it should
                            # manage them
                            kernels_trans.apply(loop)
                        else:
                            # Otherwise, if the reductions exists, they will
                            # be detected by the dependencyAnalysis and raise
                            # a TransformationError captured below
                            loop_offloading_trans.apply(
                                loop, options={"independent": True})
                            gpu_region_trans.apply(loop.ancestor(Directive))
                        # Alternatively we could use loop parallelism with
                        # reduction clauses
                    print(f"Successfully offloaded loop with {kernel_names}")
                except TransformationError as err:
                    print(f"Failed to offload loop with {kernel_names} "
                          f"because: {err}")

        # Apply OpenMP thread parallelism for any kernels we've not been able
        # to offload to GPU.
        for loop in subroutine.walk(Loop):
            if any(
                kern.name.lower() in KERNEL_EXCLUSIONS
                for kern in loop.kernels()
            ):
                continue

            if (
                not offload
                or any(
                    kern.name.lower() in (
                        list(failed_to_offload) + GPU_KERNEL_EXCLUSIONS
                    )
                    for kern in loop.kernels()
                )
            ):
                if loop.loop_type not in ["colours", "null"]:
                    cpu_parallel.apply(loop)
                    otrans.apply(loop, options={"reprod": True})

        print(subroutine.view())
