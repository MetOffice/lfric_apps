#!/usr/bin/env python
#-----------------------------------------------------------------------------
# (C) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
#-----------------------------------------------------------------------------
"""Run ESMF regrid program to generate weights"""
import subprocess
import sys
import os
from create_grid import GRID,get_env_info,get_data,validate_input,transform_and_write


class INTER:
    """Class to hold ESMF regrid options"""
    def __init__(self):
        self.method = None
        self.pole = None
        self.unmapped = None
        self.degenerate = None
        self.area = None
        self.line = None
        self.norm = None
        self.dst_grid_type = None
        self.dst_grid_name = None
        self.src_grid_type = None
        self.src_grid_name = None

    def set_regrid_info(self):
        """Set information from environment variables"""
        self.method = os.environ.get('INT_METHOD')
        self.pole = os.environ.get('INT_POLE')
        self.unmapped = os.environ.get('INT_UNMAPPED')
        self.degenerate = os.environ.get('INT_DEGENERATE')
        self.area = os.environ.get('INT_AREA')
        self.line = os.environ.get('INT_LANE')
        self.norm = os.environ.get('INT_NORM')
        self.src_grid_type = os.environ.get('SRC_GRID_TYPE')
        self.src_grid_name = os.environ.get('SRC_MESH_NAME')
        self.dst_grid_type = os.environ.get('DST_GRID_TYPE')
        self.dst_grid_name = os.environ.get('DST_MESH_NAME')
        
    def set_arguments(self):
        """Set ESMF regrid command line options"""
        options = ' -s ' + self.src_grid_name + ".nc"
        options = options + ' -d ' + self.dst_grid_name + ".nc"
        options = options + ' --dst_loc center'
        options = options + ' --src_loc center'
        options = options + ' --64bit_offset --check -m ' + self.method
        options = options + ' --extrap_method none'
        options = options + ' -p ' + self.pole + ' '
        if self.src_grid_type == 'regional' and self.dst_grid_type == 'regional':
          options = options + ' -r '
        else:
          if self.src_grid_type == 'regional':
            options = options + ' --src_regional '
          if self.dst_grid_type == 'regional':
            options = options + ' --dst_regional '
        if self.unmapped == "yes":
            options = options + ' -i '
        if self.degenerate == "yes":
            options = options + ' --ignore_degenerate '
        if self.line is not None:
            options = options + ' -l ' + self.line
        if self.area == "yes" and self.method == 'conserve':
            options = options + ' --user_areas '
        if self.norm is not None:
            options = options + ' --norm_type ' + self.norm
        return options

    def print_info(self):
        """Print options to the output"""
        print('INTERPOLATION Options:')
        print('method: ', self.method)
        print('pole: ', self.pole)
        print('unmapped: ', self.unmapped)
        print('degenerate: ', self.degenerate)
        print('area: ', self.area)
        print('line: ', self.line)
        print('norm: ', self.norm)
        print('SRC grid type: ', self.src_grid_type)
        print('DST grid type: ', self.dst_grid_type)
        print('SRC grid name: ', self.src_grid_name)
        print('DST grid name: ', self.dst_grid_name)

def run_exe(arg, outfile, style):
    """Run ESMF regrid"""
    model = ' ESMF_RegridWeightGen ' + arg + ' -w ' + outfile
    nproc = os.environ.get('MPIRUN_N_JOBS')
    cmd = 'mpiexec -n ' + nproc + model
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
    SRC_GRID = GRID()
    SRC_GRID = get_env_info(SRC_GRID, 'SRC')
    SRC_GRID = get_data(SRC_GRID)
    validate_input(SRC_GRID)
    # end processing source grid
    print('SOURCE grid ')
    DST_GRID = GRID()
    DST_GRID = get_env_info(DST_GRID, 'DST')
    DST_GRID = get_data(DST_GRID)
    # end processing destination grid
    validate_input(DST_GRID)
    print('DESTINATION grid ')

    print('Save SRC -> DST grid info')
    transform_and_write(SRC_GRID)
    transform_and_write(DST_GRID)

    MY_INTER = INTER()
    MY_INTER.set_regrid_info()
    MY_INTER.print_info()
    style = os.environ.get('FORMAT')
    run_exe(MY_INTER.set_arguments(),
                "regrid_weights.nc", style)
