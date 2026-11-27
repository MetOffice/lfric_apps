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


def add_mask(kind, mask_file, mask_var, invert):
    '''
    Add a mask field to a local copy of the LFRic mesh, for generating
    masked weights.

    param str kind: Name of the mask type (e.g. 'land'), used to name the
                    masked mesh file.
    param str mask_file: NetCDF file holding the mask on the LFRic mesh faces.
    param str mask_var: Name of the mask variable in mask_file.
    param bool invert: Whether to invert the mask before adding it.
    return: The path of the masked mesh file.

    '''
    if os.environ.get('INT_DIRECT') != 'lfric2um':
        raise Exception('Masked weights are only supported for lfric2um')
    masked_mesh = "lfric_mesh_" + kind + "_masked.nc"
    command = ["add_mask_to_mesh.py",
               "--mesh-in", os.environ.get('GRID_PATH_LFRIC'),
               "--mesh-name", os.environ.get('MESH_NAME_LFRIC'),
               "--mask-file", mask_file,
               "--mask-var", mask_var,
               "--mesh-out", masked_mesh,
               "--out-var", mask_var]
    if invert:
        command.append("--invert")
    run_command(command)
    return masked_mesh


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

    int_methods = {"P": ["bilinear", "neareststod"],
                   "U": ["bilinear", "neareststod"],
                   "V": ["bilinear", "neareststod"]}

    # Optionally generate land- and sea-masked weights. Land and sea
    # fields are only on P points, so these are only generated for the
    # P grid. The sea mask file is a land mask in which cells containing
    # any sea are 0, so it is inverted to leave just those cells unmasked
    masks = {}
    for kind, file_env, var_env, invert in [
            ("land", "MASK_FILE", "MASK_VAR", False),
            ("sea", "SEA_MASK_FILE", "SEA_MASK_VAR", True)]:
        mask_file = os.environ.get(file_env)
        if mask_file:
            mask_var = os.environ.get(var_env)
            masks[kind] = (add_mask(kind, mask_file, mask_var, invert),
                           mask_var)
            int_methods["P"].extend([kind + "_mask_bilinear",
                                     kind + "_mask_neareststod"])

    # Loop over grid types and generate weights for each
    for um_ptype, ptype_int_methods in int_methods.items():
        for int_method in ptype_int_methods:
            kind = int_method.split("_mask_")[0]
            if kind in masks:
                os.environ["GRID_PATH_LFRIC_MASKED"] = masks[kind][0]
                os.environ["MASK_VAR_LFRIC"] = masks[kind][1]
            weight_gen(um_ptype, int_method)
