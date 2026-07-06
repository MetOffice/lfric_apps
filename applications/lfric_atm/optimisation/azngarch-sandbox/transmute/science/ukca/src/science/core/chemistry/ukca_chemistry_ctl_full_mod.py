##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# This transformation introduces OpenMP directives around loops inside
# UKCA full-domain mode.

import logging
import os

from psyclone.psyir.symbols import (
    ArrayType,
    DataSymbol,
)
from psyclone.psyir.nodes import (
    Loop,
    Reference,
    Routine
)
from psyclone.transformations import (
    OMPLoopTrans,
    TransformationError
)
from psyclone.version import __MAJOR__, __MICRO__, __MINOR__

# PSyclone version
psy_version = (__MAJOR__, __MINOR__, __MICRO__)

def match_loop(loop: Loop, var_name: str, stop_name: str) -> bool:
    # Return true only if loop's variable is named var_name
    # and loop's stop expression is a reference named stop_name.
    return (loop.variable.name == var_name and
                isinstance(loop.stop_expr, Reference) and
                loop.stop_expr.name == stop_name)

def trans(psyir):
    # All loops are dynamically scheduled
    omp_trans = OMPLoopTrans(omp_directive="paralleldo",
                             omp_schedule="dynamic")

    for routine in psyir.walk(Routine):
        if routine.name == "ukca_chemistry_ctl_full":
            for loop in routine.walk(Loop):
                try:
                    # Parallelise the "DO l = 1, dim_ntp" loops
                    if match_loop(loop, "l", "dim_ntp"):
                        omp_trans.apply(loop)

                    # Parallelise the "DO jspf = 1, jpcspf" loop
                    if match_loop(loop, "jspf", "jpcspf"):
                        omp_trans.apply(loop, force=True)

                    # Parallelise the 3D chunking loop
                    if match_loop(loop, "zi", "model_levels"):
                        # Find all "chunk_" arrays (to be marked as private)
                        privates = set()
                        for sym in loop.get_all_accessed_symbols():
                            if (sym.name.startswith("chunk_") and
                                   isinstance(sym, DataSymbol) and
                                   isinstance(sym.datatype, ArrayType)):
                                privates.add(sym)

                        # Apply the transformation
                        parent, position = loop.parent, loop.position
                        omp_trans.apply(loop, force=True, collapse=3)

                        # Mark explicitly private variables
                        if psy_version < (3, 3, 0):
                            loop.explicitly_private_symbols.update(privates)
                        else:
                            directive = parent.children[position]
                            directive.explicitly_private_symbols.update(
                                privates)

                except TransformationError as err:
                    err_msg = ("ukca_chemistry_ctl_full_mod.py: Error: "
                               "could not apply OMP transformation "
                               f"to loop '{loop.variable.name}': "
                               f"{err.message_text}")
                    raise TransformationError(err_msg) from err
