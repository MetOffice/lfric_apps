#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2019 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Basic Python script to plot the x-z profile minus a constant state of 300
from a GungHo output file.

This version extracts data from UGRID output files using Iris and interpolates
onto a regulat grid for plotting.

Levels are determined from the data.

'''

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib
# Need to set a non-interactive backend for suites
matplotlib.use('Agg')  # noqa: E402
from six.moves import range

from iris.pandas import as_data_frame

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys

from read_data import read_ugrid_data


cube = None
n_levs = None


def make_figure(plotpath, nx, ny, field, timestep):

    # Get coordinates

    x = np.around(cube.coord('longitude').points)
    y = np.around(cube.coord('latitude').points)
    time = np.around(cube.coord('time').points, decimals=1)
    print('times = ', time)
    # Sort by y and retain indices for reshape
    sortidx = np.argsort(x)

    # Get min and max of x,y data for plot axes
    xmin = min(x)
    xmax = max(x)
    ymin = min(y)
    ymax = max(y)
    zmin = 0.0
    zmax = 1500.0

    scale = 1.0/1000.0

    nx = int(nx)
    ny = int(ny)
    nz = n_levs

    # Create 2D plot

    vali = np.zeros([nx, ny, n_levs])

    for tt in timestep:
        slice_fig = plt.figure(figsize=(15, 10))
        t = np.int(tt)
        for p in range(n_levs):

            # Get the data for this level
            data = cube.data[t, p]

            vali[:, :, p] = data[sortidx].reshape((nx, ny))

            # Create meshgrid to get x_i and y_i for plotting
            x2d = np.linspace(xmin, xmax, nx)
            z2d = np.linspace(zmin, zmax, nz)
            y_i, x_i = np.meshgrid(z2d, x2d)

            dz = np.zeros([nx, nz])
            background = 0.0
            if field == 'air_potential_temperature' or field == 'theta':
                background = 303.05
                cc = np.linspace(-0.1, 0.6, 15)
            elif field == 'Xi2':
                cc = np.linspace(-0.22, 0.22, 12)
            elif field == 'm_v':
                background = 0.0
                cc = np.linspace(0.013, 0.0271, 12)
            else:
                cc = np.linspace(np.amin(zi), np.amax(zi), 11)

            for i in range(nx):
                dz[i, :] = vali[i, 0, :] - background

            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
            c_map = cm.summer
            cf = plt.contourf(x_i * scale, y_i * scale,
                              np.round(dz, 10), cc, cmap=c_map)
            cl = plt.contour(x_i * scale, y_i * scale,
                             np.round(dz, 10), cc, linewidths=1.0, colors='k',
                             linestyle="", extend='min')
            plt.axis([-.5, .5, 0, 1.5])
            plt.xlabel("x (km)")
            plt.ylabel("z (km)")
            plt.title('max: %2.4e, min: %2.4e' % (np.max(dz), np.min(dz)))
            plt.colorbar(cf,  cmap=c_map)

            out_file_name = plotpath + "/" + \
                'robert_' + field + "_t" + str(time[t]) + ".png"
            slice_fig.savefig(out_file_name, bbox_inches='tight')


if __name__ == "__main__":

    try:
        datapath, nx, ny, fields, timesteps, plotpath = sys.argv[1:7]
    except ValueError:
        print(
            "Usage: {0} <datapath> <nx> <ny> <field_names> <timesteps> "
            "<plotpath>".format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')
    # Split out the list of timeslices
    time_list = timesteps.split(':')

    for field in field_list:

        cube = read_ugrid_data(datapath, field)

        n_levs = cube.data[0, :].shape[0]

        # Only try to plot if we found some data for this field
        if n_levs > 0:
            make_figure(plotpath, nx, ny, field, time_list)
