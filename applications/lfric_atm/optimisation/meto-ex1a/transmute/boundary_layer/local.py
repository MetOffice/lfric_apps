# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Local script for Boundary layer.
This by default removes the j loop(s) and replaces it with a init,
spans a PARALLEL region across the whole file,
and then adds OMP DO to the top most loop in each group of loops.
There are some small bespoke needs for the 6x files that this currently
affects which are captured below as they are initialised.

This is currently used by the following files:
* bl_lsp
* btq_int
* ex_flux_tq
* ex_flux_uv
* kmkh
* tr_mix
* imp_mix
* fm_drag
'''

import logging
from psyclone.psyir.nodes import (
    Routine,
    Loop,
    OMPParallelDirective,
)
from psyclone.transformations import (TransformationError)
from transmute_psytrans.transmute_functions import (
    loop_replacement_of,
    get_compiler,
    first_priv_red_init,
    remove_unspanable_nodes,
    set_pure_subroutines,
    replace_n_threads,
    OMP_PARALLEL_REGION_TRANS,
    OMP_DO_LOOP_TRANS_STATIC,
)
from boundary_layer.script_options import (
    SCRIPT_OPTIONS_DICT
)


Loop.set_loop_type_inference_rules({
    "i_wt": {"variable": "i_wt"},  # For ex_flux_tq.F90
    "ient": {"variable": "ient"},  # For tr_mix.F90
    "ii": {"variable": "ii"},  # For bdy_impl3.F90
    "k": {"variable": "k"},  # For all files which use script
    "j": {"variable": "j"},  # For all files which use script
    "i": {"variable": "i"}})  # For bdy_impl3.F90


# Longer term we will raise some of these into a override import, see Apps#900
# pylint: disable=too-many-locals
# pylint: disable=too-many-statements
# pylint: disable=too-many-branches
def trans(psyir):
    """
    Local.py script for boundary layer.
    This spans a PARALLEL section across the whole file,
    and then adds OMP to either to top loop of a nest, or k

    This is currently used by the following files:
    * bl_lsp
    * btq_int
    * ex_flux_tq
    * ex_flux_uv
    * kmkh
    * tr_mix
    * imp_mix
    * fm_drag
    """

    # options list for Transformation.
    options = {}
    # First privates created by CCE redundant inits.
    first_private_list = []
    # Designate calls in regions as safe to parallelise over.
    safe_pure_calls = []
    # Do we update the max_threads variable with a library call?
    max_threads_parse = False
    # Assignment nodes that we do not wish to parallelise over at the start
    loop_type_init = ["j"]

    # Lifted and extended from Global.py
    # Given the file, update the above lists with extra options

    # Get the file name to use with the SCRIPT_OPTIONS_DICT
    fortran_file_name = str(psyir.root.name)
    # Check if file is in the script_options_dict
    # Copy out anything that's needed
    # Only the options list is currently
    if fortran_file_name in SCRIPT_OPTIONS_DICT:
        file_overrides = SCRIPT_OPTIONS_DICT[fortran_file_name]
        # Update the respective lists if the filename override exists
        if "options" in file_overrides.keys():
            options = file_overrides["options"]
        if "first_private_list" in file_overrides.keys():
            first_private_list = file_overrides["first_private_list"]
        if "safe_pure_calls" in file_overrides.keys():
            safe_pure_calls = file_overrides["safe_pure_calls"]
        if "max_threads_parse" in file_overrides.keys():
            max_threads_parse = file_overrides["max_threads_parse"]
        if "loop_type_init" in file_overrides.keys():
            for name in file_overrides["loop_type_init"]:
                loop_type_init.append(name)

    # Set up some specifics to this local script

    # Remove any j loops and add an init for j
    remove_loop_type = ["j"]
    # Avoid any nodes related to the timers
    timer_routine_names = ["lhook"]

    # Set the pure calls if needed
    if safe_pure_calls:
        set_pure_subroutines(psyir, safe_pure_calls)

    # Replace max_threads = 1
    if max_threads_parse:
        replace_n_threads(psyir, "max_threads")

    # Remove any loops relating to specified loop type
    for node in psyir.walk(Routine):
        for removal_type in remove_loop_type:
            loop_replacement_of(node, removal_type)

    # Span a parallel section across the whole routine,
    # apart for a few exceptions provided
    for routine in psyir.walk(Routine):
        routine_children = remove_unspanable_nodes(
            routine,
            timer_routine_names,
            loop_type_init
        )
        # Span the region across filtered down node list
        try:
            OMP_PARALLEL_REGION_TRANS.apply(
                routine_children)
        except (TransformationError, IndexError) as err:
            logging.warning("OMPParallelTrans failed: %s", err)

    # CCE first private issue with 3.1, to be removed longer term
    if get_compiler() == "cce" and first_private_list:
        for routine in psyir.walk(Routine):
            for node in routine.children:
                if isinstance(node, OMPParallelDirective):
                    first_priv_red_init(node, first_private_list)
                    break

    # Loop ancestor type that a loop cannot have.
    avoid_loop_ancestor = ["ii", "k", "i"]
    # A loop type which a loop cannot have an OMP do section
    avoid_loop_type = ["i_wt", "ient"]
    # Work through the loops now in the spanned section and try transformations
    for loop in psyir.walk(Loop):
        loop_ancestor_type = ""
        try:
            loop_ancestor_type = loop.ancestor(Loop).loop_type
        # pylint: disable=bare-except
        except:  # noqa: E722
            pass
        # Default is there is no Loop ancestor
        if (((not loop.ancestor(Loop)) or
                # Or the rest, the loop ancestor is not to be avoided
                loop_ancestor_type not in avoid_loop_ancestor) and
                # And the loop is not of certain loop types
                str(loop.loop_type) not in avoid_loop_type):
            try:
                OMP_DO_LOOP_TRANS_STATIC.apply(
                    loop, options)
            except (TransformationError, IndexError) as err:
                logging.warning("OMPLoopTrans failed: %s", err)
