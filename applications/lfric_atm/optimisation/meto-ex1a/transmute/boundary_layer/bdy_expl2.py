# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Transmute script for bdy_expl2
'''
import logging
from psyclone.psyir.nodes import (
    OMPParallelDirective,
    OMPBarrierDirective,
    Routine,
    Loop,
    Assignment,
    Call
)
from psyclone.transformations import (TransformationError)
from transmute_psytrans.transmute_functions import (
    loop_replacement_of,
    replace_n_threads,
    get_compiler,
    first_priv_red_init,
    set_pure_subroutines,
    get_ancestors,
    get_descendents,
    OMP_PARALLEL_REGION_TRANS,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
    OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC,
    OMP_DO_LOOP_TRANS_DYNAMIC,
    OMP_DO_LOOP_TRANS_STATIC,
)

Loop.set_loop_type_inference_rules({
    "ient": {"variable": "ient"},
    "i_wt": {"variable": "i_wt"},
    "ii": {"variable": "ii"},
    "k": {"variable": "k"},
    "j": {"variable": "j"},
    "i": {"variable": "i"},
    "l": {"variable": "l"}})


# pylint: disable=too-many-locals
def trans(psyir):
    """
    Main section to apply OpenMP Directives to bdy_expl2
    """
    # A list of children nodes will be generated for the subroutine.
    # At these given nodes, start and end a parallel region.
    # Each time a region is spanned, under the hood this list is updated by
    # PSyclone and the next index changes as the new region becomes a node in
    # this list.
    parallel_section_descriptors = [
        # LHS  , # RHS
        ["ntml_save", "ntml"],          # start
        ["bl_diag", "fb_surf"],         # stop
        ["grad_t_adj", "min"],          # start
        ["bl_diag","dvdzm"],      # stop
        ["zh", "z_uv"],          # start
        # ["wtrac_bl", "zero"],           # stop 
        # ["cumulus", "false"],           # start
        ["cloud_base_found", "true"],  # stop
        ["weight1", "r_theta_levels"],  # start
        ["ntml", "ntml_local"],         # stop
        # ["mix_len_bm", "max"],          # start
        # ["mix_len_bm", "bm_tiny"],      # stop
        # ["visc_m", "visc_m"],           # start
        # ["rhokm", "visc_m"],            # stop
        ]
    # Note: These only go one level deep when spanning.
    # If blocks count as nodes, and this routine is full of them.
    # As such these are a little conservative with the regions,
    # so we may loose out a tiny bit of performance with
    # no waits, but this can be improved in the future if needed.

    # Allow spanning parallel and do clauses over the following calls
    safe_pure_calls = [
        "qsat",
        "qsat_mix",
        "qsat_wat",
        "qsat_wat_mix"]
    # A list of known variables to redundantly initialise with CCE
    first_private_list = [
        "frac_lev",
        "qc_tot",
        "slope",
        "zpr",
        "ntop",
        "z_scale",
        "f_log",
        "lambdah",
        "vkz",
        "weight1",
        "weight2",
        "weight3",
        "sh",
        "var_fac",
        "km",
        "kp",
        "var_fac"]

    ## Currently identical to avoid_loop_type ## 
    parrallel_avoid_loop_type = [
        "k",
        "i_wt",
        "ient"
        ]

    # Replace max_threads = 1
    replace_n_threads(psyir, "max_threads")

    # Set the pure
    set_pure_subroutines(psyir, safe_pure_calls)

    bdy_expl2_routine = None
    # Grab the routine reference to use throughout the script
    for routine in psyir.walk(Routine):
        if routine.name == "bdy_expl2":
            bdy_expl2_routine = routine

    # Strip out loops relating to cells_j loop type.
    # for routine in psyir.walk(Routine):
    #     if routine.name == "bdy_expl2":
    #         loop_replacement_of(routine, "j")
    loop_replacement_of(bdy_expl2_routine, "j")

    # for routine in psyir.walk(Routine):
    #     if routine.name == "bdy_expl2":
    routine_children = bdy_expl2_routine.children

    current_node = bdy_expl2_routine
    

    # # Grab the indexes given the provided list of nodes and targets
    parallel_section_index_locations = find_node_index(
        routine_children,
        parallel_section_descriptors)

    # # These are intially indexes in a unaltered tree.
    # # As we span sections, we need to account for the loss in nodes
    parallel_section_index_locations = correct_index_list_for_trans(
        parallel_section_index_locations)

    # span parallel regions
    for parallel_region_bounds in parallel_section_index_locations:
        start = parallel_region_bounds[0]
        end = parallel_region_bounds[1]
        try:
            OMP_PARALLEL_REGION_TRANS.apply(
                routine_children[start:end])
        except (TransformationError, IndexError) as err:
            logging.warning(
                "Could not transform because:\n %s", err)

    # Setup for both inside region and otherwise
    dynamic_loop_type = ["ii"]
    avoid_loop_direct_ancestor = [
        "i",
        "ii"
        ]
    avoid_loop_type = [
        "k",
        "i_wt",
        "ient"
        ]

    # options = {
    #     "ignore_dependencies_for": [
    #         "zh_local"
    #         ]}

    options = {
        "ignore_dependencies_for": [
            "visc_m",
            "visc_h",
            "sigma_h",
            "cumulus",
            "l_shallow",
            "bl_type_5",
            "bl_type_6",
            "cu_over_orog",
            "ntml",
            "topbl",
            "zh_local",
            "pstar",
            "qssurf",
            "tstar",
            "p_theta_levels",
            "qs_tl",
            "tl"],
        "nowait": True
        }

    for node in routine_children:

        # To be used in the parallel sections to help determine good nodes
        if isinstance(node, Call):
            for pure_call in safe_pure_calls:
                node_ast_string = str(node.ast).lower()
                if node_ast_string.find(pure_call) != -1:
                    print("Found a pure call")

        if isinstance(node,OMPParallelDirective):

            parallel_dir_loops = node.walk(Loop)
            parallel_dir_loops = get_descendents(node, Loop)
            
            for loop in parallel_dir_loops:

                if str(loop.loop_type) not in avoid_loop_type:

                    static_transformation = OMP_DO_LOOP_TRANS_STATIC

                    all_ancestors = get_ancestors(loop)

                    loop_check_ancestors = True
                    for loop_ancestor in all_ancestors:
                        if loop_ancestor.variable.name in avoid_loop_direct_ancestor:
                            print(loop_ancestor.variable.name)
                            print(avoid_loop_direct_ancestor)
                            loop_check_ancestors = False

                    if (not loop.ancestor(Loop) or loop_check_ancestors): 
                        try:
                            static_transformation.apply(
                                loop, options)
                        except (TransformationError, IndexError) as err:
                            # logging.warning(
                            #     "Could not transform because:\n %s", err)
                            print("Could not transform as")
                            print(err)

    removing_unsafe_barriers = True
    while removing_unsafe_barriers is True:
        removing_unsafe_barriers = False
        for node in routine_children:
            if isinstance(node, OMPBarrierDirective):
                print("Unsafe Barrier removal, should be fixed in PSyclone")
                tmp = node.detach()
                removing_unsafe_barriers = True
                break

    # # Setup options for the do inside a spanned region
    # options = {
    #     "ignore_dependencies_for": [
    #         "zh_local"
    #         ]}
    #     #     ,
    #     # "nowait": True
    #     # }

    # for node in routine_children:
    #     print(node)

    # # First pass to add the OMP DO
    # parallelise_loops(
    #     psyir,
    #     True,
    #     avoid_loop_type,
    #     avoid_loop_direct_ancestor,
    #     dynamic_loop_type,
    #     options)

    # options = {
    #     "ignore_dependencies_for": [
    #         "visc_m",
    #         "visc_h",
    #         "sigma_h",
    #         "cumulus",
    #         "l_shallow",
    #         "bl_type_5",
    #         "bl_type_6",
    #         "cu_over_orog",
    #         "ntml",
    #         "topbl",
    #         "zh_local"]}

    # # Second pass to add OMP Parallel do over remaining loops
    # parallelise_loops(
    #     psyir,
    #     False,
    #     avoid_loop_type,
    #     avoid_loop_direct_ancestor,
    #     dynamic_loop_type,
    #     options)

    # if get_compiler() == "cce":
    #     for routine in psyir.walk(Routine):
    #         for node in routine.children:
    #             if isinstance(node, OMPParallelDirective):
    #                 first_priv_red_init(node, first_private_list)
    #                 break


# Loop through loops and parallelise
# pylint: disable=dangerous-default-value
# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments
def parallelise_loops(
    psyir,
    inside_parallel=False,
    avoid_loop_type=[],
    avoid_loop_direct_ancestor=[],
    dynamic_loop_type=[],
    options={}
):

    '''
    For bdy_expl2, insert do directives given a specific set of circumstances.
    This is called twice, once inserting OMP DO inside PARALLEL regions and a
    second adding PARALLEL DOs around loops without existing regions.
    '''

    # Are we expecting to insert OMP DO or PARALLEL DOs
    # Grab a reference to pre configured transformations
    if inside_parallel is True:
        dynamic_transformation = OMP_DO_LOOP_TRANS_DYNAMIC
        static_transformation = OMP_DO_LOOP_TRANS_STATIC
    else:
        dynamic_transformation = OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC
        static_transformation = OMP_PARALLEL_LOOP_DO_TRANS_STATIC

    for loop in psyir.walk(Loop):

        # Check if we want to avoid this loop type or not
        if str(loop.loop_type) not in avoid_loop_type:
            # Setup some basic checks to be used below
            proceed_with_do = False

            # Check the OMP Parallel Ancestry against our expectation
            loop_omp_ancestor = loop.ancestor(OMPParallelDirective)
            # Expecting to be inside a parallel and there is an ancestor
            if loop_omp_ancestor and inside_parallel:
                proceed_with_do = True
            # Expecting to not be inside a parallel
            # and there is not an ancestor
            elif not loop_omp_ancestor and not inside_parallel:
                proceed_with_do = True

            # If we okay to proceed given the above
            if proceed_with_do:
                # Grab the loop ancestry
                # List of all to check number of nests
                all_ancestors = []
                all_ancestors = get_ancestors(loop)
                # A loop ancestor type
                # pylint: disable=bare-except
                loop_ancestor_type = ""
                try:
                    loop_ancestor_type = loop.ancestor(Loop).loop_type
                except:  # noqa: E722
                    pass

                # The main if for guarding adding a directive
                # It is strict on the ancestry, with given exceptions.
                # ( There is not an ancestor
                #   OR
                #   ( The ancestor is not in the avoid list
                #    AND there are not multiple ancestors )
                # AND (everything) the loop is safe given the above checks.
                if (not loop.ancestor(Loop) or
                        (((loop_ancestor_type not in avoid_loop_direct_ancestor) and
                            len(all_ancestors) <= 1))):
                    # Depending on whether we want the schedule to be
                    # OMP static or dynamic for this loop.
                    if loop.loop_type in dynamic_loop_type:
                        transformation_do_sched = dynamic_transformation
                    else:
                        transformation_do_sched = static_transformation
                    try:
                        transformation_do_sched.apply(
                            loop, options)
                    except (TransformationError, IndexError) as err:
                        logging.warning(
                            "Could not transform because:\n %s", err)


def find_node_index(routine_children, parallel_section_descriptors):
    '''
    Find a list of indexes of the first time a node contianing an assignment
    which matches the values in the inside array.
    If multiple nodes have similar lhs/rhs, it is good to pick
    unique references.
    '''

    # work inits
    next_region = []
    parallel_section_index_locations = []
    step_indexer = 0

    # Work through the children of the routine, hold the index for use later
    for index, node in enumerate(routine_children):
        # Reset variables per routine child node
        match_index = False
        list_of_attributes = []
        # Get the children of a node, if they exisit
        try:
            list_of_attributes = node.walk(Assignment)
        except AttributeError:
            pass

        # Walk through assignments of the child.
        for assignment in list_of_attributes:
            # Check each assignment
            match_index = check_parallel_section_list(
                parallel_section_descriptors[step_indexer],
                assignment
                )

            # If we found a matching node
            if match_index is True:

                # Is it the start node of a parallel section
                if step_indexer % 2 == 0:
                    next_region.append(index)
                # Or the end node of a parallel section
                else:
                    next_region.append(index+1)
                    parallel_section_index_locations.append(next_region)
                    next_region = []

                # Cycle whether it is start of stop, position in the passed
                # parallel_section_descriptors list will describe this
                # Note: If any fail to match, this will fall out of sync.
                step_indexer += 1

            # If we can found the number of matches against the list provided
            # Return the list of indexes
            if step_indexer >= len(parallel_section_descriptors):
                return parallel_section_index_locations

    # If there is no match, there is a default return of an empty list.
    return []


# pylint: disable=too-many-branches
def check_parallel_section_list(lhs_rhs_node, assignment):
    '''
    For this method, we are checking a LHS and RHS string which
    will match against an assignment contained within a node.
    This node is ideally a loop (though if blocks work), where
    we can span a parallel section.

    '''

    # Setup some defaults
    # LHS can always be grabbed if it is an assignment
    lhs_string = str(assignment.lhs.name)
    # The RHS is not known yet
    rhs_string = ""

    # if lhs_string.lower() == "cumulus":
    #     print(assignment.rhs)
    #     print(assignment.rhs.value)

    # If the LHS matches the LHS of the 'node representation' from our
    # list of a desired targeted nodes properties.
    if lhs_string.lower() == lhs_rhs_node[0].lower():
        # Now we check for the RHS
        # There are many RHS options which could be present
        # We try and capture them below.

        # Setup some defaults
        check_string = str(lhs_rhs_node[1].lower())
        rhs_string = ""
        rhs_name_found = False
        rhs_ast_found = False
        rhs_value_found = False

        # Try and get the RHS name
        # pylint: disable=bare-except
        try:
            rhs_name = assignment.rhs.name
            if rhs_name is not None:
                rhs_name_found = True
        except:  # noqa: E722
            pass
        # Try and get the ast (intrinsic function) reference
        # pylint: disable=bare-except
        try:
            rhs_ast = assignment.rhs.ast
            if rhs_ast is not None:
                rhs_ast_found = True
        except:  # noqa: E722
            pass

        # Try and get the assignment (intrinsic function) reference
        # pylint: disable=bare-except
        try:
            rhs_value = assignment.rhs.value
            if rhs_value is not None:
                rhs_value_found = True
        except:  # noqa: E722
            pass

        # This has to loop through the children too sadly
        if not rhs_name_found and not rhs_ast_found:
            for child in assignment.rhs.children:
                # pylint: disable=bare-except
                try:
                    rhs_child_name = str(child.name).lower()
                    if rhs_child_name == check_string:
                        return True
                        # break
                except:  # noqa: E722
                    pass

        # If there is a name string, and not an ast (intrinsic function)
        if rhs_name_found and not rhs_ast_found:
            rhs_string = str(assignment.rhs.name).lower()
            if rhs_string == check_string:
                return True

        # If there is not a name, but is a ast (intrinsic function)
        elif rhs_ast_found and not rhs_name_found:
            # Check if check_string is in rhs_string
            # in does not seem to work, this does.
            rhs_string = str(assignment.rhs.ast).lower()
            if rhs_string.find(check_string) != -1:
                return True

        # If there is a tre/false value
        elif rhs_value_found and not rhs_name_found and not rhs_ast_found:
            rhs_string = str(assignment.rhs.value).lower()
            if rhs_string == check_string:
                return True

        else:
            return False

    return False


def correct_index_list_for_trans(index_list):
    '''
    Each time a parallel section is spanned, all of the nodes contained inside
    are placed inside a the new node. This means that the index list is out of
    date. This function corrects the indexes by reducing them by the number
    of nodes spanned inside each new node.
    '''

    # Start the first indexes at 0 and do not effect them
    reduce_by = 0
    # For each Start / Stop parallel section pair inside the list
    for start_stop_pair in index_list:
        # Reduce the Start / Stop index by the difference between them
        start_stop_pair[0] = start_stop_pair[0] - reduce_by
        start_stop_pair[1] = start_stop_pair[1] - reduce_by

        # Increase the reduce_by variable derived by the difference
        # between the current pair ready for the next
        reduce_by = reduce_by + (start_stop_pair[1] - start_stop_pair[0] - 1)

    # Return the updated list
    return index_list
