#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file LICENCE
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""
Script to generate um_grid namelist file

A namelist file for output must be provided. This is currently given in:
lfricinputs/rose-stem/app/generate_weights/weights_gen_wrapper.py

"""

import sys
import os


def write_nml_file(igrid_targ, rotated, lambda_pole, phi_pole, lambda_points,
                   phi_points, delta_lambda, delta_phi, lambda_origin,
                   phi_origin):
    '''
    Overwrites an existing namelist file with coordinate and grid information
    set within this script and read from a namelist file (NLIST_FILE) which
    is currently accessed as an environment variable that stores a path.

    '''

    fname = os.environ.get('NLIST_FILE')
    with open(fname, "w+") as inoutfile:
        inoutfile.write('&GRID\n'
                        f'  POINTS_LAMBDA_TARG = {lambda_points}\n'
                        f'  POINTS_PHI_TARG = {phi_points}\n'
                        f'  LAMBDA_ORIGIN_TARG = {lambda_origin}\n'
                        f'  PHI_ORIGIN_TARG = {phi_origin}\n'
                        f'  PHI_POLE = {phi_pole}\n'
                        f'  LAMBDA_POLE = {lambda_pole}\n'
                        f'  ROTATED = {rotated}\n'
                        f'  DELTA_LAMBDA_TARG = {delta_lambda}\n'
                        f'  DELTA_PHI_TARG = {delta_phi}\n'
                        f'  IGRID_TARG = {igrid_targ}\n'
                        '/\n')


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--resolution", help="resolution of UM grid",
                        type=int)
    parser.add_argument("-d", "--dump", help="path to UM dump file to use")
    args = parser.parse_args()

    # pylint: disable=invalid-name
    # TODO: Make constants conform to UPPER_CASE style
    igrid_targ = 6

    if args.resolution:

        res = args.resolution
        print('Creating GLOBAL UM NON-ROTATED grid namelist file from',
              'given resolution: ', res)

        # Basic, unrotated global mesh
        lambda_pole = 0.0
        phi_pole = 90.0
        rotated = 'F'                           # False

        # Start at south pole
        # N
        start_lon = 0.0
        start_lat = -90.0

        lambda_points = int(2 * res)            # Number of longitude lines
        phi_points = int(3 * res / 2)           # Number of latitude lines
        delta_lambda = 360.0 / lambda_points    # Line spacing
        delta_phi = 180.0 / phi_points          # Line spacing
        lambda_origin = start_lon + (delta_lambda / 2)
        phi_origin = start_lat + (delta_phi / 2)

        write_nml_file(igrid_targ, rotated, lambda_pole, phi_pole,
                       lambda_points, phi_points, delta_lambda, delta_phi,
                       lambda_origin, phi_origin)

    elif args.dump:

        # Mule documentation:
        # https://code.metoffice.gov.uk/doc/um/mule/2022.05.1/examples.html
        import mule
        df_path = args.dump
        print('Creating UM grid namelist file from dump: ', df_path)
        df = mule.DumpFile.from_file(df_path)   # decoding the dump file

        if df.fixed_length_header.grid_staggering != 6:
            print('ERROR: Only ENDGAME grids are currently valid input')
            sys.exit(1)

        lambda_pole = df.real_constants.north_pole_lon
        phi_pole = df.real_constants.north_pole_lat

        if lambda_pole != 0.0 or phi_pole != 90.0:
            rotated = 'T'   # True
        else:
            rotated = 'F'   # False

        # Unlike when resolution is specified, all of the information we
        # need is provided in the dump file.
        start_lon = df.real_constants.start_lon
        start_lat = df.real_constants.start_lat

        lambda_points = df.integer_constants.num_cols
        phi_points = df.integer_constants.num_rows
        delta_lambda = df.real_constants.col_spacing
        delta_phi = df.real_constants.row_spacing
        lambda_origin = start_lon + (delta_lambda / 2)
        phi_origin = start_lat + (delta_phi / 2)

        write_nml_file(igrid_targ, rotated, lambda_pole, phi_pole,
                       lambda_points, phi_points, delta_lambda, delta_phi,
                       lambda_origin, phi_origin)

    else:

        print('')
        print('Use the --help or -h flags to see commandline options')
        print('')
