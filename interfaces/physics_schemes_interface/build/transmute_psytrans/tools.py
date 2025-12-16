# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Contains some functions which might assist in generating scripts, but otherwise
would not be used in the script actively.
'''


# Tool, not active part of script
def find_node_index(routine_children, find_lhs, find_rhs):
    '''
    A tool to assist developers finding the index of a node to
    use when spanning a parallel section. 
    You can use the variable name of the LHS and either the variable
    name of the RHS (if it exisists), or the name of the function used
    on the RHS.
    This is called with a list of routine_children inside a walk of the
    routines in the PSyIR object.
    '''
    for index, node in enumerate(routine_children):
        list_of_attributes = []
        try:
            list_of_attributes = node.walk(Assignment)
        except AttributeError:
            pass

        for assignment in list_of_attributes:
            if str(assignment.lhs.name) == find_lhs:
                found_rhs=False
                print("#### New find ####")
                print("found LHS")
                print(index)
                try:
                    print(assignment.rhs.ast)
                    if find_rhs in str(assignment.rhs.ast):
                        found_rhs = True
                except:
                    pass
                try:
                    if str(assignment.rhs.name) == find_rhs:
                        found_rhs = True
                except:
                    pass
                if found_rhs:
                    print("Found RHS - Y")
                    print(index)
                    break
