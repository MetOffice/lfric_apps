#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file LICENCE
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
# Some content in this file was generated or refactored with assistance from
# - Met Office Github Copilot Enterprise (Claude Opus 5.5)
'''
Wrapper script for a set of commands/script used to create weights files.

'''

import subprocess
import os


def run_command(command):
    '''
    Takes in a list of commands and, if n commands > 1, pipes the commands
    to the CLI and executes.

    param list command: A list of commands (incl files) to execute.

    '''
    cmd = subprocess.Popen(command, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    stdout, stderr = cmd.communicate()
    return_code = cmd.returncode
    print(' ')
    print('Output from running: ', ' '.join(command), ':')
    print(stdout.decode())
    print(stderr.decode())
    if return_code != 0:
        raise Exception(' '.join(command) + ' failed.')


def weight_gen(um_ptype, int_method):
    '''
    This function sets environment variables from passed arguments as well
    as setting the output filename dependent on the interpolation method
    (int_method).

    After setting env vars and outfile name, weight_gen utilises the
    run_command function declared in this module to create a um_grid object,
    generate weights for this object, and relocate the output to the outfile
    name.

    param str um_ptype: The grid type that dictates the coordinate shift
                        performed in 'create_um_grid.py'.
    param str int_method: The interpolation method that will dictate the
                          output filename and some weight generation settings.

    '''
    os.environ["GRID"] = um_ptype
    os.environ["INT_METHOD"] = int_method
    os.environ["GRID_PATH_UM"] = "UM_grid_" + um_ptype + ".nc"
    direct = os.environ.get('INT_DIRECT')
    if direct == 'um2lfric':
        outfilename = um_ptype + "_to_FACE_CENTRE_" + int_method + ".nc"
    elif direct == 'lfric2um':
        outfilename = "FACE_CENTRE_to_" + um_ptype + "_" + int_method + ".nc"
    else:
        print('Interpolation direction not supported')

    run_command(["create_um_grid.py"])
    run_command(["generate_weights.py"])
    run_command(["mv", "regrid_weights.nc", outfilename])


if __name__ == "__main__":

    option = os.environ.get('UM_GRID_INPUT')

    # Conditionals to check if namelist file needs to be generated and
    # calls appropriate script with environment variables set as
    # arguments.
    # Fixed resolution: generates namelist from 'NLIST_FILE'
    if option == 'dump-file':
        argument = os.environ.get('DUMPPATH')
        run_command(["create_um_grid_nml.py", "-d", argument])
    elif option == 'specified-resolution':
        argument = os.environ.get('RESOLUTION')
        run_command(["create_um_grid_nml.py", "-r", argument])

    int_methods = ["bilinear", "neareststod"]

    # Optionally add a mask field to a local copy of the LFRic mesh, used
    # only for the land_mask_bilinear weights
    mask_file = os.environ.get('MASK_FILE')
    if mask_file:
        if os.environ.get('INT_DIRECT') != 'lfric2um':
            raise Exception('MASK_FILE is only supported for lfric2um')
        masked_mesh = "lfric_mesh_masked.nc"
        run_command(["add_mask_to_mesh.py",
                     "--mesh-in", os.environ.get('GRID_PATH_LFRIC'),
                     "--mesh-name", os.environ.get('MESH_NAME_LFRIC'),
                     "--mask-file", mask_file,
                     "--mask-var", os.environ.get('MASK_VAR'),
                     "--mesh-out", masked_mesh,
                     "--out-var", os.environ.get('MASK_VAR')])
        os.environ["GRID_PATH_LFRIC_MASKED"] = masked_mesh
        os.environ["MASK_VAR_LFRIC"] = os.environ.get('MASK_VAR')
        int_methods.append("land_mask_bilinear")

    # Loop over grid types and generate weights for each
    for um_ptype in ["P", "U", "V"]:
        for int_method in int_methods:
            weight_gen(um_ptype, int_method)
