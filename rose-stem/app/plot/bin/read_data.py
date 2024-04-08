#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Python script for reading LFRic diagnostic output data with separate functions
for nodal text file data and UGRID NetCDF
'''

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import pandas as pd

import glob
import sys

import iris

if iris.__version__ < "3.0.0":
    iris.FUTURE.netcdf_promote = True

'''
Read from per-processor nodal text files
# filestem - path to files
# ncomp - number of components (3 for vector field and 1 for scalar)
# comp - which component we are interested in (for filtering levels)
'''


def read_nodal_data(filestem, ncomp, comp):

    if ncomp == 3:  # vector fields with 3 components
        col_names = ['x', 'y', 'z', 'level', 'c1', 'c2', 'c3']
        col_types = {'x': np.float32, 'y': np.float32, 'z': np.float32,
                     'level': np.float32, 'c1': np.float64,
                     'c2': np.float64, 'c3': np.float64}
    else:         # scalar fields with 1 component
        col_names = ['x', 'y', 'z', 'level', 'c1']
        col_types = {'x': np.float32, 'y': np.float32, 'z': np.float32,
                     'level': np.float32, 'c1': np.float64}

    # get the list of files to stitch together
    dirlist = glob.glob(filestem)

    # Create an empty DataFrame
    all_data = pd.DataFrame()

    # If no files are found then don't try to process them and issue a warning
    if len(dirlist) < 1:
        print("Warning: No files found to plot")
    else:

        for f in dirlist:

            tmp_data = pd.read_csv(f, header=None, names=col_names, skiprows=1,
                                   comment=']', dtype=col_types, sep='\s+')

            all_data = all_data.append(tmp_data)

        # Select full or half levels for vector fields, if required
        if ncomp == 3:  # vector field
            if comp != '3':  # u or v component half levels
                all_data = all_data.loc[all_data['level'].mod(1) == 0.5]
            else:  # w component full levels
                all_data = all_data.loc[all_data['level'].mod(1) == 0.0]

    return all_data

def load_cube_by_varname(filename, var):
    variable_constraint = iris.Constraint(cube_func=(
                                lambda c: c.var_name == var))
    return iris.load_cube(filename, constraint=variable_constraint)

def read_ugrid_data(filestem, field):
    '''
    Read from UGRID NetCDF file
    # filestem - path to file
    # field - field name
    '''
    # Read the file and extract field of interest
    cube = load_cube_by_varname(filestem, field)

    return cube