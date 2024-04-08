#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file LICENCE
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''Run ESMF regrid program to generate weights'''
import subprocess
import sys
import os


class INTER:
    '''Class to hold ESMF regrid options'''
    def __init__(self):
        self.method = None
        self.pole = None
        self.points = None
        self.unmapped = None
        self.degenerate = None
        self.area = None
        self.line = None
        self.norm = None
        self.lfric_grid_type = None
        self.lfric_grid_path = None
        self.um_grid_type = None
        self.um_grid_path = None

    def set_regrid_info(self):
        '''Set information from environment variables'''
        self.method = os.environ.get('INT_METHOD')
        self.pole = os.environ.get('INT_POLE')
        self.unmapped = os.environ.get('INT_UNMAPPED')
        self.degenerate = os.environ.get('INT_DEGENERATE')
        self.area = os.environ.get('INT_AREA')
        self.line = os.environ.get('INT_LANE')
        self.norm = os.environ.get('INT_NORM')
        self.interpolation_direction = os.environ.get('INT_DIRECT')
        self.um_grid_type = os.environ.get('GRID_TYPE_UM')
        self.um_grid_path = os.environ.get('GRID_PATH_UM')
        self.lfric_grid_type = os.environ.get('GRID_TYPE_LFRIC')
        self.lfric_grid_path = os.environ.get('GRID_PATH_LFRIC')

    def set_arguments_um2lfric(self):
        '''Set ESMF regrid command line options'''
        options = ' -s ' + self.um_grid_path
        options = options + ' -d ' + self.lfric_grid_path
        options = options + ' --dst_loc center'
        options = options + ' --64bit_offset --check -m ' + self.method
        options = options + ' --extrap_method neareststod'
        if self.um_grid_type == 'regional':
            options = options + ' --src_regional '
        if self.lfric_grid_type == 'regional':
            options = options + ' --dst_regional '
        if self.unmapped == "yes":
            options = options + ' -i '
        if self.degenerate == "yes":
            options = options + ' --ignore_degenerate '
        if self.line is not None:
            options = options + ' -l ' + self.line
        if self.area == "yes" and self.method == 'conserve':
            options = options + ' --user_areas '
        options = options + ' -p ' + self.pole + ' '
        if self.norm != "ignore":
            options = options + ' --norm_type ' + self.norm
        return options

    def set_arguments_lfric2um(self):
        '''Set ESMF regrid command line options'''
        options = ' -s ' + self.lfric_grid_path
        options = options + ' --src_loc center'
        options = options + ' -d ' + self.um_grid_path
        options = options + ' --64bit_offset --check -m ' + self.method
        if self.um_grid_type == 'regional':
            options = options + ' --src_regional '
        if self.lfric_grid_type == 'regional':
            options = options + ' --dst_regional '
        if self.unmapped == "yes":
            options = options + ' -i '
        if self.degenerate == "yes":
            options = options + ' --ignore_degenerate '
        if self.line is not None:
            options = options + ' -l ' + self.line
        if self.area == "yes" and self.method == 'conserve':
            options = options + ' --user_areas '
        options = options + ' -p ' + self.pole + ' '
        if self.norm != "ignore":
            options = options + ' --norm_type ' + self.norm
        return options

    def print_info(self):
        '''Print options to the output'''
        print('INTERPOLATION Options:')
        print('method: ', self.method)
        print('pole: ', self.pole)
        print('pole points: ', self.points)
        print('unmapped: ', self.unmapped)
        print('degenerate: ', self.degenerate)
        print('area: ', self.area)
        print('line: ', self.line)
        print('norm: ', self.norm)
        print('UM grid type: ', self.um_grid_type)
        print('LFRic grid type: ', self.lfric_grid_type)
        print('Interpolation direction: ', self.interpolation_direction)


def run_exe(arg, outfile, style):
    '''Run ESMF regrid'''
    model = ' ESMF_RegridWeightGen ' + arg + ' -w ' + outfile
    nproc = os.environ.get('MPIRUN_N_JOBS')
    cmd = 'mpirun -n ' + nproc + model
    print(cmd)
    retcode = subprocess.call(cmd, shell=True)
    if retcode != 0:
        print("Execution terminated by the signal", retcode)
        sys.exit(1)

    if style == "SCRIP":
        cmd = "ncrename -h "
        cmd = cmd + "-d n_a,src_grid_size -d n_b,dst_grid_size "
        cmd = cmd + "-d n_s,num_links "
        cmd = cmd + "-d nv_a,src_grid_corners -d nv_b,dst_grid_corners "
        cmd = cmd + "-v yc_a,src_grid_center_lat -v yc_b,dst_grid_center_lat "
        cmd = cmd + "-v xc_a,src_grid_center_lon -v xc_b,dst_grid_center_lon "
        cmd = cmd + "-v yv_a,src_grid_corner_lat -v yv_b,dst_grid_corner_lat "
        cmd = cmd + "-v xv_a,src_grid_corner_lon -v xv_b,dst_grid_corner_lon "
        cmd = cmd + "-v mask_a,src_grid_imask -v mask_b,dst_grid_imask "
        cmd = cmd + "-v area_a,src_grid_area -v area_b,dst_grid_area "
        cmd = cmd + "-v frac_a,src_grid_frac -v frac_b,dst_grid_frac "
        cmd = cmd + "-v col,src_address -v row,dst_address "
        cmd = cmd + outfile
        print(cmd)
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            print("Execution terminated by the signal", retcode)
            sys.exit(1)
        cmd = 'ncap2 -h -s "remap_matrix[num_links,num_wgts]=S" ' + outfile
        print(cmd)
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            print("Execution terminated by the signal", retcode)
            sys.exit(1)
        cmd = "ncks --abc -O -h -x -v S " + outfile + " -o " + outfile
        print(cmd)
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            print("Execution terminated by the signal", retcode)
            sys.exit(1)
        cmd = 'ncatted -h -O -a '
        cmd = cmd + 'history,global,a,c,'
        cmd = cmd + '"Convert to SCRIP style weights file using NCO tools\n" '
        cmd = cmd + outfile
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            print("Execution terminated by the signal", retcode)
            sys.exit(1)


if __name__ == "__main__":
    MY_INTER = INTER()
    MY_INTER.set_regrid_info()
    MY_INTER.print_info()
    style = os.environ.get('FORMAT')
    direct = os.environ.get('INT_DIRECT')
    if direct == 'um2lfric':
        run_exe(MY_INTER.set_arguments_um2lfric(),
                "regrid_weights.nc", style)
    elif direct == 'lfric2um':
        run_exe(MY_INTER.set_arguments_lfric2um(),
                "regrid_weights.nc", style)
    else:
        print('Interpolation direction not supported')
