##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# Summary
# =======
#
# This transformation is used to adjust the number of points passed to the ASAD
# initialisation routine to account for full-domain chunking.  The chunk size
# is taken at compile time from the environment variable UKCA_FULL_CHUNK_SIZE.
# If this variable is not set then the source code is passed through
# unmodified.
#
# Example
# =======
#
# Give the program
#
#   subroutine main()
#     if (ukca_config%l_ukca_asad_full) then
#       n_pnts = tot_n_pnts
#     endif
#   end subroutine
#
# and the parameters
#
#   match_if  = "ukca_config%l_ukca_asad_full"
#   match_lhs = "n_pnts"
#   match_rhs = "tot_n_pnts"
#   UKCA_FULL_CHUNK_SIZE = 256
#
# the following program is produced.
#
#   subroutine main()
#     if (ukca_config%l_ukca_asad_full) then
#       n_pnts = 256
#       umPrint("Initialising ASAD solver with 256 points")
#     end if
#   end subroutine

# Imports
# =======

import os

from psyclone.psyir.nodes import (
    Assignment, Reference, Literal, IfBlock, Call)
from psyclone.psyir.symbols import (
    INTEGER_TYPE, RoutineSymbol, CHARACTER_TYPE)

# Transformation Parameters
# =========================

# The if condition to match against
match_if = "ukca_config%l_ukca_asad_full"

# The left-hand-side to match against
match_lhs = "n_pnts"

# The right-hand-side to match against
match_rhs = "tot_n_pnts"

# Transformation
# ==============


def trans(psyir):
    chunk_size = os.getenv("UKCA_FULL_CHUNK_SIZE")
    if chunk_size is None:
        # Do nothing if the chunk size is not set
        return
    elif chunk_size == "FULL_DOMAIN":
        # Message to print (via umPrint) when chunking enabled
        message_text = "Initialising ASAD solver for full-domain"
        # We use None to represent the full-domain chunk size
        chunk_size = None
    else:
        # Message to print (via umPrint) when chunking enabled
        message_text = ("Initialising ASAD solver with " +
                        chunk_size +
                        " points due to full-domain chunking")

    # Search and replace
    found = None
    for ifblock in psyir.walk(IfBlock):
        if isinstance(ifblock.condition, Reference):
            (sig, _inds) = ifblock.condition.get_signature_and_indices()
            if str(sig) == match_if:
                for assign in ifblock.walk(Assignment):
                    if (
                            isinstance(assign.lhs, Reference) and
                            assign.lhs.name == match_lhs and
                            isinstance(assign.rhs, Reference)
                    ):
                        if assign.rhs.name == match_rhs:
                            if chunk_size is not None:
                                assign.rhs.replace_with(
                                    Literal(str(chunk_size), INTEGER_TYPE))
                            found = assign

    # Insert print call
    if found:
        print_call = Call()
        print_call.addchild(Reference(RoutineSymbol("umPrint")))
        print_call.addchild(Literal(message_text, CHARACTER_TYPE))
        found.parent.addchild(print_call, index=found.position+1)
