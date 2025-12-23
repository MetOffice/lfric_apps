# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Transmute script for bdy_expl2
'''
import logging
from psyclone.psyir.nodes import OMPParallelDirective, Node, Routine, Loop, Assignment, IfBlock, Container, Literal, Schedule
from psyclone.transformations import (TransformationError)
from transmute_psytrans.transmute_functions import (
    loop_replacement_of,
    replace_n_threads,
    get_compiler,
    first_priv_red_init,
    set_pure_calls,
    get_ancestors,
    loop_init_under,
    ancestor_loop_init_over,
    OMP_PARALLEL_REGION_TRANS,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
    OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC,
    OMP_DO_LOOP_TRANS_DYNAMIC,
    OMP_DO_LOOP_TRANS_STATIC,
)
#from transmute_psytrans.tools import find_node_index

Loop.set_loop_type_inference_rules({
    "ient": {"variable": "ient"},
    "i_wt": {"variable": "i_wt"},
    "ii": {"variable": "ii"},
    "k": {"variable": "k"},
    "j": {"variable": "j"},
    "i": {"variable": "i"},
    "l": {"variable": "l"}})


def trans(psyir):
    """
    Main section to apply OpenMP Directives to bdy_expl2
    """
    # A list of children nodes will be generated for the subroutine.
    # At these given nodes, start and end a parallel region.
    # Each time a region is spanned, under the hood this list is updated by
    # PSyclone and the next index changes as the new region becomes a node in
    # this list.
    # It might be perceived as fragile, however capturing the right loop any
    # other way is currently impractical, and there is a simple tool to assist
    # finding the node added with this ticket, #1056.
    # Parallelise over these specific nodes:
    # parallel_section_index_locations = [
    #     [5,10],
    #     [6,19],
    #     [8,10],
    #     [11,13],
    #     [17,19],
    #     [22,23],
    #     [27,30],
    #     [30,36],
    # ]
    parallel_section_descriptors = [
        ["ntml_save","ntml"], # start
        ["BL_diag","fb_surf"], # stop
        ["grad_t_adj","min"],
        ["seg_slice_start","ii"],
        # ["grad_t_adj","max_t_grad"],
        # ["BL_diag","dvdzm"],
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

    # Replace max_threads = 1
    replace_n_threads(psyir, "max_threads")

    # Set the pure
    set_pure_calls(psyir, safe_pure_calls)

    # Strip out loops relating to cells_j loop type.
    for routine in psyir.walk(Routine):
        if routine.name == "bdy_expl2":
            loop_replacement_of(routine,"j")


    for routine in psyir.walk(Routine):
        if routine.name == "bdy_expl2":
            routine_children=routine.children

            # Grab the indexes given the provided list of nodes and targets
            parallel_section_index_locations = find_node_index(
                routine_children,
                parallel_section_descriptors)

            # These are intially indexes in a unaltered tree.
            # As we span sections, we need to account for the loss in nodes
            parallel_section_index_locations = correct_index_list_for_trans(
                parallel_section_index_locations)

            # span parallel regions
            for parallel_region_bounds in parallel_section_index_locations:
                start=parallel_region_bounds[0]
                end=parallel_region_bounds[1]
                try:
                    OMP_PARALLEL_REGION_TRANS.apply(routine_children[start:end])
                except (TransformationError, IndexError) as err:
                    print("Could not span parallel because:")
                    print(err)
            

    # Setup for both inside region and otherwise
    check_loop_value = {}
    check_loop_value["k"] = 2
    dynamic_loop_type = ["ii"]

    # Setup options for the do inside a spanned region
    avoid_loop_ancestor = ["k","i","ii"]
    avoid_loop_type = ["i_wt"]
    options = {"ignore_dependencies_for": ["sigma_h","bl_diag%dvdzm","topbl","zh_local"]}

    #First pass to add the OMP DO
    # parallelise_loops(
    #     psyir,
    #     True,
    #     avoid_loop_type,
    #     avoid_loop_ancestor,
    #     check_loop_value,
    #     dynamic_loop_type,
    #     options)   

    # A loop type which a loop cannot have an OMP parallel do section
    avoid_loop_type = ["i_wt","ient","k"]
    avoid_loop_ancestor = ["i","ii"]
    options = {"ignore_dependencies_for": ["visc_m","visc_h"]}

    # Second pass to add OMP Parallel do over remaining loops
    # parallelise_loops(
    #     psyir,
    #     False,
    #     avoid_loop_type,
    #     avoid_loop_ancestor,
    #     check_loop_value,
    #     dynamic_loop_type,
    #     options)

    if get_compiler() == "cce":
        for routine in psyir.walk(Routine):
            for node in routine.children:
                if isinstance(node, OMPParallelDirective):
                    first_priv_red_init(node, first_private_list)
                    break


# Loop through loops and parallelise
def parallelise_loops(
    psyir,
    inside_parallel=False,
    avoid_loop_type=[],
    avoid_loop_ancestor=[],
    check_loop_value={},
    dynamic_loop_type=[],
    options={}):

    '''
    For bdy_expl2, insert do directives given a specific set of circumstances.
    This is called twice, once inserting OMP DO inside PARALLEL regions and a
    second adding PARALLEL DOs around loops without existing regions.
    '''

    # Are we expecting to insert OMP DO or PARALLEL DOs
    # Grab a reference to pre configured transformations
    if inside_parallel == True:
        dynamic_transformation = OMP_DO_LOOP_TRANS_DYNAMIC
        static_transformation = OMP_DO_LOOP_TRANS_STATIC
    else:
        dynamic_transformation = OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC
        static_transformation = OMP_PARALLEL_LOOP_DO_TRANS_STATIC

    for loop in psyir.walk(Loop):

        # Check if we want to avoid this loop type or not
        if str(loop.loop_type) not in avoid_loop_type:
            # Setup some basic checks to be used below
            nest_loop = False
            is_safe_loop = True
            proceed_with_do = False
        
            # Check the OMP Parallel Ancestry against our expectation
            loop_OMP_ancestor = loop.ancestor(OMPParallelDirective)
            # Expecting to be inside a parallel and there is an ancestor
            if loop_OMP_ancestor and inside_parallel:
                proceed_with_do = True
            # Expecing to not be insidea a parallel and there is not an ancestor
            elif not loop_OMP_ancestor and not inside_parallel:
                proceed_with_do = True

            # If we okay to proceed given the above
            if proceed_with_do:
                # Grab the loop ancestry
                # List of all to check number of nests
                all_ancestors=[]
                all_ancestors=get_ancestors(loop)
                # A loop ancestor type 
                # (wraped in try/except in case there is not suitable ancestor)
                loop_ancestor_type=""
                try:
                    loop_ancestor_type=loop.ancestor(Loop).loop_type
                except:
                    pass

                # Check if this is a safe loop, default is it is
                is_safe_loop = loop_init_under(
                    loop, 
                    check_loop_value)

                # Check if the loop is under a node which we are not
                # parallelising over, making this nestable
                nest_loop = ancestor_loop_init_over(
                    loop,
                    loop_ancestor_type, 
                    check_loop_value)

                # The main if for guarding adding a directive
                # It is strict on the ancestry, with given exceptions.
                # ( There is not an ancestor
                #   OR 
                #   ( ( The ancestor is not in the avoid list
                #     OR
                #     Is it a nestable loop ) 
                #    AND there are not multiple ancestors ) ) 
                # AND (everything) the loop is safe given the above checks
                if (
                (not loop.ancestor(Loop) or
                (((loop_ancestor_type not in avoid_loop_ancestor) or 
                nest_loop) and
                len(all_ancestors) <=1 )) and
                is_safe_loop):
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
                        print("Could not transform because:")
                        print(err)


def find_node_index(routine_children, parallel_section_descriptors):
    '''
    Find a list of indexes of the first time a node contianing an assignment
    which matches the values in the inside array.
    If multiple nodes have similar lhs/rhs, it is good to pick unique references.
    '''

    # parallel_section_descriptors = [
    #     ["ntml_save","ntml"],["BL_diag","fb_surf"],
    #     ["grad_t_adj","max_t_grad"],["BL_diag","dvdzm"],
    #     ]
    
    # work inits
    next_region = []
    parallel_section_index_locations = []
    step_indexer = 0
    
    #print(a_parallel_section)
    for index, node in enumerate(routine_children):
        match_index = None
        list_of_attributes = []
        try:
            list_of_attributes = node.walk(Assignment)
        except AttributeError:
            pass

        for assignment in list_of_attributes:
            # insert new bits here
            
            #print(a_parallel_section[step_indexer])
            match_index = check_parallel_section_list(
                index, 
                parallel_section_descriptors[step_indexer],
                assignment)

            print(match_index)

            if match_index is not None:
            
                if step_indexer % 2 == 0:
                    next_region.append(match_index)
                else:
                    next_region.append(match_index+1)
                    parallel_section_index_locations.append(next_region)
                    next_region = []
                
                step_indexer+=1
                print(step_indexer)

            if step_indexer >= len(parallel_section_descriptors):
                print(parallel_section_index_locations)
                return parallel_section_index_locations

    return []

def check_parallel_section_list(index, lhs_rhs_node, assignment):

    lhs_string = str(assignment.lhs.name)
    rhs_string = ""

    if lhs_string.lower() == lhs_rhs_node[0].lower():
        print("Found LHS")
        rhs_name_found = False
        rhs_ast_found = False
        found_rhs=False
        check_string = str(lhs_rhs_node[1].lower())
        rhs_string = ""
        print(lhs_rhs_node)
        try:
            rhs_store1 = assignment.rhs.name
            if rhs_store1 is not None:
                rhs_name_found = True
            # if rhs_string == check_string:
            #     found_rhs =   True
        except:
            pass
        try:
            rhs_store2 = assignment.rhs.ast
            if rhs_store2 is not None:
                rhs_ast_found = True
            # # Check if check_string is in rhs_string
            # # in does not seem to work, this does.
            # if rhs_string.find(check_string) != -1:
            #     found_rhs = True
        except:
            pass

        if rhs_name_found and not rhs_ast_found:
            #print("1st")
            rhs_string = str(assignment.rhs.name).lower()
            if rhs_string == check_string:
                found_rhs = True

        elif not rhs_name_found and rhs_ast_found:
           # print("Second")
            # Check if check_string is in rhs_string
            # in does not seem to work, this does.
            rhs_string = str(assignment.rhs.ast).lower()
            if rhs_string.find(check_string) != -1:
                found_rhs = True

        # print(rhs_name_found)
        # print(rhs_ast_found)

        if found_rhs:
            print("Found RHS - Y")
            # print(index)
            # print(rhs_string)
            # #print(rhs_string1)
            # print(check_string)
            # we find the first match and return
            # if there are duplicates, pick a different one
            return index


def correct_index_list_for_trans(index_list):

    reduce_by = 0
    for start_stop_pair in index_list:
        start_stop_pair[0] = start_stop_pair[0] - reduce_by
        start_stop_pair[1] = start_stop_pair[1] - reduce_by

        reduce_by = start_stop_pair[1] - start_stop_pair[0] - 1

    print(index_list)
    return index_list