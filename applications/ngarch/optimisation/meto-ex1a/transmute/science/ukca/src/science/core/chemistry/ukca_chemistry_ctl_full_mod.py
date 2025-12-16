##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# Summary
# =======
#
# This transformation introduces a chunking loop around the call to the ASAD
# solver in `ukca_chemistry_ctl_full_mod.F90`.  A single call to the ASAD
# solver is replaced with multiple calls, each of which operates on a "chunk"
# of the full domain. The chunk size is taken at compile time from the
# environment variable UKCA_FULL_CHUNK_SIZE.  If this variable is not set then
# the source code is passed through unmodified.
#
# Chunking is mainly achieved by slicing the arguments to the ASAD solver.
# However, the ASAD solver is dependent not only on its arguments but
# also on global ASAD arrays, several of which are read and written directly by
# UKCA full-domain mode.  When chunking is enabled, any ASAD arrays that were
# originally full-domain sized become chunk sized. To resolve this size change,
# the transformer needs to be told about all of these ASAD arrays via
# the following parameter.
#
#   * asad_vars:              a dict mapping the ASAD arrays (originally
#                             full-domain sized), which are accessed by UKCA
#                             full-domain mode, to their associated ranks
#
# Unfortunately, this parameter cannot be inferred automatically because the
# ASAD arrays are dynamically allocated (and hence we don't know which ones
# were originally full-domain sized at compile time). Any changes to UKCA
# full-domain mode must therefore ensure that asad_vars is updated, if
# necessary.
#
# For each variable in asad_vars, the transformer introduces a
# full-domain-sized counterpart.  Outside the chunking loop, it renames each
# access of a narrow (chunk sized) ASAD array into an access of its wide
# (full-domain sized) counterpart.  Inside the chunking loop, slices of the
# newly introduced wide arrays are copied into narrow ASAD arrays, then the
# ASAD solver is called, and then narrow ASAD arrays are copied back into
# slices of their wide counterparts.
#
# In addition to asad_vars, the transformer uses the following parameters.
#
#   * fulldom_size_name:      name of variable holding full-domain size
#   * asad_call_name:         name of the top-level ASAD solver routine
#
# Example
# =======
#
# Given the program
#
#   subroutine main()
#     integer, parameter :: n = 256
#     integer :: arr1(n)
#     integer :: arr2(n, 2)
#     integer :: asad_arr1(32)
#     arr1(:) = foo
#     asad_arr1(:) = arr1(:) + 1
#     call asad_cdrive(arr1, arr2, n)
#     arr1(:) = arr1(:) + asad_arr1(:)
#     arr1(:) = arr1(:) + 1
#   end subroutine
#
# and the parameters
#
#   fulldom_size_name     = "n"
#   asad_vars             = {"asad_arr1" : 1}
#   asad_call_name        = "asad_cdrive"
#   UKCA_FULL_CHUNK_SIZE  = 32
#
# the following program is produced.
#
#   subroutine main()
#     integer, parameter :: n = 256
#     integer :: arr1(n)
#     integer :: arr2(n, 2)
#     integer :: asad_arr1(32)
#     integer :: full_asad_arr1(n)
#     integer :: chunk_begin, chunk_end, chunk_size
#
#     arr1(:) = foo
#     full_asad_arr1(:) = arr1(:) + 1
#     do chunk_begin = 1, n, 32
#       chunk_end = min(n, chunk_begin+chunk_size-1)
#       chunk_size = 1 + chunk_end - chunk_begin
#       asad_arr1(1:chunk_size) = full_asad_arr1(chunk_begin:chunk_end)
#       call asad_cdrive(arr1(chunk_begin:chunk_end),                       &
#                        arr2(chunk_begin:chunk_end, :),                    &
#                        chunk_size)
#       full_asad_arr1(chunk_begin:chunk_end) = asad_arr1(1:chunk_size)
#     end do
#     arr1(:) = arr1(:) + full_asad_arr1(:)
#     arr1(:) = arr1(:) + 1
#   end subroutine

# Imports
# =======

import os

from psyclone.version import (__MAJOR__, __MINOR__, __MICRO__)
from psyclone.psyir.nodes import (
    ArrayReference, Assignment, BinaryOperation, Call, IntrinsicCall,
    Literal, Loop, Routine, Reference, Schedule)
from psyclone.psyir.symbols import (
    INTEGER_TYPE, REAL_TYPE, CHARACTER_TYPE, Symbol, DataSymbol, ArrayType,
    RoutineSymbol)
from psyclone.psyir.transformations.reference2arrayrange_trans import (
    Reference2ArrayRangeTrans)

psy_version = (__MAJOR__, __MINOR__, __MICRO__)

# Transformation Parameters
# =========================

# Name of variable holding the full-domain size
fulldom_size_name = "tot_n_pnts"

# ASAD arrays in use (and their ranks)
asad_vars = {"rk":  2, "sph2o": 1, "sphno3": 1, "tnd": 1,
             "y":   2, "za":    1, "dpd":    2, "dpw": 2,
             "prk": 2, "fpsc1": 1, "fpsc2":  1}

# Name of the top-level ASAD call
asad_call_name = "asad_cdrive"

# Name of the routine in which to apply the transformation
routine_name = "ukca_chemistry_ctl_full"

# Transformation
# ==============


def trans(psyir):
    desired_chunk_size = os.getenv("UKCA_FULL_CHUNK_SIZE")
    if desired_chunk_size is None:
        # Do nothing if the chunk size is not set
        return
    elif desired_chunk_size == "FULL_DOMAIN":
        # Message to print (via umPrint) when chunking enabled
        message_text = ("UKCA full-domain chunking enabled with " +
                        "a chunk size equal to the size of the full " +
                        "domain")
        # We use None to represent the full-domain chunk size
        desired_chunk_size = None
    else:
        # Message to print (via umPrint) when chunking enabled
        message_text = ("UKCA full-domain chunking enabled with " +
                        "a chunk size of " + desired_chunk_size)

    # Locate correct routine within which to apply the transformation
    for routine in psyir.walk(Routine):
        if routine.name != routine_name:
            continue

        # Find variable holding the size of the full domain
        try:
            array_size_var = routine.symbol_table.lookup(fulldom_size_name)
        except Exception:
            continue

        # Find the call to the ASAD solver
        asad_call = None
        for call in routine.walk(Call):
            if call.routine.name == asad_call_name:
                if call.parent is routine:
                    asad_call = call
        if asad_call is None:
            continue

        # Find references to ASAD arrays before and after the call
        # --------------------------------------------------------

        refs_before = set()
        refs_after = set()
        for stmt in routine.children[:asad_call.position]:
            for ref in stmt.walk(Reference):
                if ref.name in asad_vars:
                    refs_before.add(ref.name)
        for stmt in routine.children[asad_call.position+1:]:
            for ref in stmt.walk(Reference):
                if ref.name in asad_vars:
                    refs_after.add(ref.name)

        # Introduce variable to hold the desired chunk size
        # -------------------------------------------------

        desired_chunk_size_var = routine.symbol_table.find_or_create_tag(
            "desired_chunk_size",
            symbol_type=DataSymbol,
            datatype=INTEGER_TYPE)

        if desired_chunk_size is None:
            assign_desired_chunk_size = Assignment.create(
                Reference(desired_chunk_size_var),
                Reference(array_size_var))
        else:
            assign_desired_chunk_size = Assignment.create(
                Reference(desired_chunk_size_var),
                Literal(str(desired_chunk_size), INTEGER_TYPE))

        # Introduce full-domain array for each ASAD array
        # -----------------------------------------------

        full_vars = {}
        for (var_name, var_rank) in asad_vars.items():
            # Create array bounds
            bounds = [Reference(array_size_var)]
            for i in range(2, var_rank+1):
                bounds.append(IntrinsicCall.create(
                    IntrinsicCall.Intrinsic.SIZE,
                    [Reference(Symbol(var_name)),
                     ("dim", Literal(str(i), INTEGER_TYPE))]))
            # Create variables
            new_var = routine.symbol_table.find_or_create_tag(
                "full_" + var_name,
                symbol_type=DataSymbol,
                datatype=ArrayType(REAL_TYPE, bounds))
            full_vars[var_name] = (bounds, new_var)
            # Add initialiser
            if var_name in refs_before:
                initialiser = Assignment.create(
                    ArrayReference.create(new_var, [":" for b in bounds]),
                    Literal("0.0", REAL_TYPE))
                routine.addchild(initialiser, index=0)

        # Replace each use of ASAD array with full-domain counterpart
        # -----------------------------------------------------------

        for stmt in routine.children:
            for ref in stmt.walk(Reference):
                if ref.name in full_vars:
                    (_var_bounds, var_sym) = full_vars[ref.name]
                    ref.symbol = var_sym

        # Identify array references in call
        # ---------------------------------

        for arg in asad_call.arguments:
            for ref in arg.walk(Reference):
                ref2arraytrans = Reference2ArrayRangeTrans()
                if not isinstance(ref, ArrayReference):
                    if psy_version <= (3, 1, 0):
                        if ref.is_array:
                            ref2arraytrans.apply(ref)
                    else:
                        if ref.symbol.is_array:
                            ref2arraytrans.apply(
                                ref, allow_call_arguments=True)

        # Create a new "chunking" loop
        # ----------------------------

        # Create loop variables
        chunk_begin_var = routine.symbol_table.find_or_create_tag(
            "chunk_begin",
            symbol_type=DataSymbol,
            datatype=INTEGER_TYPE)
        chunk_end_var = routine.symbol_table.find_or_create_tag(
            "chunk_end",
            symbol_type=DataSymbol,
            datatype=INTEGER_TYPE)
        chunk_size_var = routine.symbol_table.find_or_create_tag(
            "chunk_size",
            symbol_type=DataSymbol,
            datatype=INTEGER_TYPE)

        # Create assignment for chunk_end
        minop = IntrinsicCall.create(
            IntrinsicCall.Intrinsic.MIN,
            [Reference(array_size_var),
             BinaryOperation.create(
                 BinaryOperation.Operator.ADD,
                 Reference(chunk_begin_var),
                 BinaryOperation.create(
                     BinaryOperation.Operator.SUB,
                     Reference(desired_chunk_size_var),
                     Literal("1", INTEGER_TYPE)))])
        assign_chunk_end = Assignment.create(Reference(chunk_end_var), minop)

        # Create assignment for chunk_size
        chunk_size = BinaryOperation.create(
            BinaryOperation.Operator.SUB,
            BinaryOperation.create(
                BinaryOperation.Operator.ADD,
                Literal("1", INTEGER_TYPE),
                Reference(chunk_end_var)),
            Reference(chunk_begin_var))
        assign_chunk_size = Assignment.create(Reference(chunk_size_var),
                                              chunk_size)

        # Create chunking loop
        loop = Loop(variable=chunk_begin_var)
        asad_call.replace_with(loop)
        loop.children = [Literal("1", INTEGER_TYPE),
                         Reference(array_size_var),
                         Reference(desired_chunk_size_var),
                         Schedule(parent=loop, children=[asad_call])]

        # Copy full-size arrays into chunk-size arrays
        for var_name in refs_before:
            (bounds, full_var) = full_vars[var_name]
            var_sym = DataSymbol(
                var_name, datatype=ArrayType(REAL_TYPE, bounds))
            assign_full_var = Assignment.create(
                ArrayReference.create(var_sym, [":" for b in bounds]),
                ArrayReference.create(full_var, [":" for b in bounds]))
            loop.loop_body.addchild(assign_full_var, index=asad_call.position)

        # Copy chunk-size arrays into full-size arrays
        for var_name in refs_after:
            (bounds, full_var) = full_vars[var_name]
            var_sym = DataSymbol(
                var_name, datatype=ArrayType(REAL_TYPE, bounds))
            assign_full_var = Assignment.create(
                ArrayReference.create(full_var, [":" for b in bounds]),
                ArrayReference.create(var_sym, [":" for b in bounds]))
            loop.loop_body.addchild(assign_full_var)

        # Update references to fulldom_size_name
        for ref in loop.loop_body.walk(Reference):
            if ref.name == fulldom_size_name:
                ref.replace_with(Reference(chunk_size_var))

        # Update references to arrays
        for ref in loop.loop_body.walk(ArrayReference):
            if ref.name in asad_vars.keys():
                ref.indices[0].start = Literal("1", INTEGER_TYPE)
                ref.indices[0].stop = Reference(chunk_size_var)
            else:
                ref.indices[0].start = Reference(chunk_begin_var)
                ref.indices[0].stop = Reference(chunk_end_var)

        # Adding chunk assignments
        loop.loop_body.addchild(assign_chunk_end, index=0)
        loop.loop_body.addchild(assign_chunk_size, index=1)

        # Add print statement
        # -------------------

        print_call = Call()
        print_call.addchild(Reference(RoutineSymbol("umPrint")))
        print_call.addchild(Literal(message_text, CHARACTER_TYPE))
        loop.parent.addchild(print_call, index=loop.position)

        # Assign desired chunk size
        # -------------------------

        loop.parent.addchild(assign_desired_chunk_size, index=loop.position)
