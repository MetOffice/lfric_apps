#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Prints a set of paths for use in Make scripts. These paths are derived from
the list of PSyAD files from psyad_files_list.txt.
'''

from pathlib import Path
from default_parser import default_parser


if __name__ == "__main__":

    PARSER = default_parser()
    PARSER.add_argument("-f",
                        nargs="+",
                        type=Path,
                        help="List of relative paths to PSyAD files",
                        required=True)
    PARSER.add_argument("-opt",
                        help="Path option select, "
                             "can be either 'base' or 'kernel'",
                        required=True)
    ARGS = PARSER.parse_args()

    # We want a set of the directories
    DIRS = set()
    for FILE_ in ARGS.f:
        if ARGS.opt == "base":
            # Adding to set everything before kernel
            DIRS.add(str(FILE_).rsplit("/kernel/", maxsplit=1)[0])
        elif ARGS.opt == "kernel":
            # Adding to set everything before last instance of "/",
            # and after source.
            DIRS.add(str(FILE_.parent).rsplit("/source/", maxsplit=1)[-1])
        else:
            pass

    # Outputting to Make variable
    for DIR_ in DIRS:
        print(DIR_)
