##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Function to parse the precision settings of the build config of a rose-stem
task. Will search for a default precision specified by NNbit and then any mixed
precision settings with format eg. rbl32. Based on the settings it will decide
what the default should be (most common) and return a string with the new
default and mixed precision settings in alphabetical order. If this new string
doesn't match the task name the rose-stem suite will raise an error.
Intended as a Jinja2 Custom Filter
"""

import re


def get_precision(build_string):
    """
    Main function for Jinja2Filter
    Parse a build config string from a task name
    Return a precision config string and dictionary of precisions
    """

    # Find a default with format 'NNbit'
    # Use 64 if not set
    try:
        default = int(re.search(r"(\d+)bit", build_string).group(1))
    except AttributeError:
        default = 64

    # Dictionary of type names and their default values
    rbl = ["rbl", default]
    rdef = ["rdef", default]
    rphys = ["rphys", default]
    rsolver = ["rsolver", 32]
    rtran = ["rtran", default]
    types = [rbl, rdef, rphys, rsolver, rtran]
    precisions = {}

    # Loop over all types and search for 'typeNN' string matches in build config
    # Extract the precision, NN, as an integer.
    # If type not in the build_string use the default
    for precision_type in types:
        precision_name = precision_type[0]
        try:
            num = int(
                re.search(rf"{precision_name}(\d+)", build_string).group(1)
            )
        except AttributeError:
            num = precision_type[1]

        precisions[precision_name] = num

    # Construct the output string defining the precision
    # Potentially change the default such that it is the most common
    # This ensures that all tasks with the same build get the same build string
    precision_string = f"{default}bit"
    for precision_type in types:
        num = precisions[precision_type[0]]
        if num != precision_type[1]:
            precision_string += f"-{precision_type[0]}{num}"

    return [precisions, precision_string]
