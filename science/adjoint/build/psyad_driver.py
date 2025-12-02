#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Add the adjoint tests to the driver module, based on the psyad_files in
psyad_vars.py. Invoked in psyad_driver.mk script.
'''
import subprocess
from datetime import datetime
from pathlib import Path
from default_parser import default_parser


def add_to_driver(psyad_file: Path, driver_target: Path) -> None:
    '''
    Adds the adjoint test to the driver.
    '''
    # Determining kernel name from psyad file path.
    name = str(psyad_file).rsplit("/", maxsplit=1)[-1]
    name = name.replace("_kernel_mod", "")
    name = name.split(".")[0]

    # Generating name of algorithm from kernel name.
    if name[0:3] == "tl_":
        alg_name = f"atlt_{name[3:]}"
    else:
        alg_name = f"adjt_{name}"

    now = datetime.now().strftime("%H:%M:%S")
    print(f"{now} *Adding* {name} to driver")
    module_line = f"use {alg_name}_alg_mod, only : {alg_name}_alg"
    call_line = f"call {alg_name}_alg( mesh, chi, panel_id )"
    # Inserts commands before comment in driver
    module_command = fr'sed -i "/! <adjoint_test_mod>/i \   \ {module_line}"'
    module_command = f'{module_command} {driver_target}'
    call_command = fr'sed -i "/! <adjoint_test_call>/i \   \ {call_line}"'
    call_command = f'{call_command} {driver_target}'
    subprocess.call(module_command, shell=True)
    subprocess.call(call_command, shell=True)


if __name__ == "__main__":

    PARSER = default_parser()
    PARSER.add_argument("-d",
                        type=Path,
                        help="Path to driver target file.",
                        required=True)
    PARSER.add_argument("-f",
                        nargs="+",
                        type=Path,
                        help="List of relative paths to PSyAD files",
                        required=True)
    ARGS = PARSER.parse_args()

    for FILE in ARGS.f:
        add_to_driver(FILE, ARGS.d)
