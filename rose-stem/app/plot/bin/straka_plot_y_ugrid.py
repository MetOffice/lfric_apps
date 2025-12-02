#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Basic python script to plot the x-z profile minus a constant state of 300
from a Dynamo output file.

This version extracts data from UGRID output files using Iris and interpolates
onto a regulat grid for plotting

Levels are determined from the data.


'''

import sys
import numpy as np
from read_data import read_ugrid_data, load_cube_by_varname
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')  # noqa: E402

import matplotlib.pyplot as plt
import matplotlib.cm as cm

cube = None
n_levs = None


def make_figure(plotpath, nx, ny, field, timestep=-1):

    """  Create a figure for the input field and output to file. """

    # Get coordinates
    mesh_cube = load_cube_by_varname(datapath, cube.attributes['mesh'])
    if mesh_cube.attributes['geometry'] == 'planar':
        y_coord_name = "projection_y_coordinate"
    else:
        y_coord_name = "latitude"

    y = np.around(cube.coord(y_coord_name).points)

    time = np.around(cube.coord('time').points)[timestep]

    # Sort by y and retain indices for reshape
    sortidx = np.argsort(y)

    slice_fig = plt.figure(figsize=(15, 10))
    # get min and max of x,y data for plot axes
    ymin = min(y)
    ymax = max(y)
    zmin = 0.0
    zmax = 6400.0

    r2d = 1.0/1000.0

    nx = int(nx)
    ny = int(ny)
    nz = n_levs

    # create 2D plot

    vali = np.zeros([ny, nx, n_levs])
    yi = np.zeros([ny, nx, n_levs])

    for p in range(n_levs):

        # get the data for this level
        data = cube.data[timestep, p]

        vali[:, :, p] = data[sortidx].reshape((ny, nx))
        yi[:, :, p] = y[sortidx].reshape((ny, nx))

    # create meshgrid to get x_i and y_i for plotting
    y2d = np.linspace(ymin, ymax, ny)
    z2d = np.linspace(zmin, zmax, nz)
    y_i, x_i = np.meshgrid(z2d, y2d)

    dz = np.zeros([ny, nz])
    for i in range(ny):
        dz[i, :] = vali[i, 0, :] - 300.0

    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    c_map = cm.summer
    cc = np.linspace(-16, -1, 16)
    cf = plt.contourf(x_i * r2d, y_i * r2d, np.round(dz, 10), cc, cmap=c_map)
    plt.contour(x_i * r2d, y_i * r2d, np.round(dz, 10), cc,
                linewidths=1.0, colors='k', extend='min')
    plt.axis([0, 16, 0, 5])
    plt.xlabel("y (km)")
    plt.ylabel("z (km)")
    plt.title('max: %2.4e, min: %2.4e' % (np.max(dz), np.min(dz)))
    plt.colorbar(cf, cmap=c_map)

    out_file_name = plotpath + \
        "/" + 'straka_y_ugrid' + "_T{:06n}".format(time) + ".png"
    slice_fig.savefig(out_file_name, bbox_inches='tight')


if __name__ == "__main__":

    try:
        datapath, nx, ny, fields, plotpath = sys.argv[1:6]
    except ValueError:
        print("Usage: {0} <datapath> <nx> <ny> <field_names> <plotpath>"
              .format(sys.argv[0]))
        sys.exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    for field in field_list:

        cube = read_ugrid_data(datapath, field)

        n_levs = cube.data[0, :].shape[0]

        # Only try to plot if we found some data for this field
        if n_levs > 0:
            make_figure(plotpath, nx, ny, field)
