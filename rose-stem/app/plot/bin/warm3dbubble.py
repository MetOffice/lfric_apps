#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2018 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
##############################################################################
'''
Basic python script to plot the x-z profile minus a constant state of 303
from a Dynamo output file.

This version takes nodal format output files and
interpolates onto a regular grid.

Levels are determined from the data

This version stitches together a directory of files
and extracts all levels so it can work in the serial
case where there is one file or the parallel case where
there is a file for each processor.

'''


from __future__ import absolute_import
from __future__ import print_function
import numpy as np
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys

from read_data import read_nodal_data

# Use viridis colormap
from python_maps import viridis_data
from matplotlib.colors import ListedColormap
viridis = ListedColormap(viridis_data, name='viridis')
plt.register_cmap(name='viridis', cmap=viridis)

levels = None
data = None


def make_figure(plotpath, nx, ny, field, component, timestep):

    val_col = 'c' + str(component)

    # get min and max of x,y data for plot axes
    min_lev = min(levels)

    xmin = data.loc[data['level'] == min_lev]['x'].min()
    xmax = data.loc[data['level'] == min_lev]['x'].max()
    ymin = data.loc[data['level'] == min_lev]['y'].min()
    ymax = data.loc[data['level'] == min_lev]['y'].max()

    zmin = 0.0
    zmax = 1500.0

    r2d = 1.0

    nx = int(nx)
    ny = int(ny)
    nz = len(levels)
    zi = np.zeros([ny, nx, len(levels)])

    c_map = viridis

    for p in range(len(levels)):
        p_data = data.loc[data['level'] == levels[p]]
        zi[:, :, p] = (p_data[val_col].values).reshape((ny, nx))

    background = 0.0
    if field == 'theta':
        background = 300.0

    cc = np.linspace(0.05, 0.5, 10)

    if timestep == 'T000000':
        plotlevel = 38
    elif timestep == 'T000160':
        plotlevel = 54
    else:
        plotlevel = 80

    f = plt.figure(figsize=(20, 10))

    # x-z plot

    ax1 = f.add_axes([.1, .1, .3, .8])

    # create meshgrid to get x_i and y_i for plotting
    x2d = np.linspace(xmin, xmax, nx)
    z2d = np.linspace(zmin, zmax, nz)
    y_i, x_i = np.meshgrid(z2d, x2d)

    zp = z2d[plotlevel]

    dz = np.zeros([nx, len(levels)])
    for i in range(nx):
        dz[i, :] = zi[ny//2, i, :] - background

    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cf = ax1.contourf(x_i * r2d, y_i * r2d, np.round(dz, 10),
                      cc, cmap=c_map, extend='max')
    cl = ax1.contour(x_i * r2d, y_i * r2d, np.round(dz, 10), cc,
                     linewidths=2.0, colors='k', linestyle="", extend='min')
    ax1.plot([x2d[0], x2d[-1]], [zp, zp], 'k--', linewidth=2)

    ax1.set(xlim=[-500, 500], ylim=[0, 1000], aspect=1)
    ax1.set_xlabel("x (m)", fontsize=32)
    ax1.set_ylabel("z (m)", fontsize=32)
    ax1.set_xticks(np.arange(-500, 750, 250))
    ax1.tick_params(axis='both', labelsize=32)

    # Plot one-dimensional slice
    slice_1d_fig = plt.figure(figsize=(10, 10))
    plt.plot(np.round(dz[nx//2, :], 10), y_i[nx//2, :] * r2d, 'k', linewidth=4)
    plt.ylim([0, 1000])
    plt.ylabel("z (m)", fontsize=24)
    plt.xlabel(r"$\Delta \theta$ (K)", fontsize=24)
    plt.tick_params(axis='both', labelsize=24)
    out_file_name = plotpath + "/" + field + "_1d_" + timestep + ".png"
    slice_1d_fig.savefig(out_file_name, bbox_inches='tight')

    # x-y plot

    ax2 = f.add_axes([.5, .1, .3, .8])

    # create meshgrid to get x_i and y_i for plotting
    x2d = np.linspace(xmin, xmax, nx)
    y2d = np.linspace(ymin, ymax, ny)
    y_i, x_i = np.meshgrid(y2d, x2d)

    dz = zi[:, :, plotlevel] - background

    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cf = ax2.contourf(x_i * r2d, y_i * r2d, np.round(dz, 10),
                      cc, cmap=c_map, extend='max')
    cl = ax2.contour(x_i * r2d, y_i * r2d, np.round(dz, 10), cc,
                     linewidths=2.0, colors='k', linestyle="", extend='min')

    ax2.set(xlim=[-500, 500], ylim=[-500, 500], aspect=1)
    ax2.set_xlabel("x (m)", fontsize=32)
    ax2.set_ylabel("y (m)", fontsize=32)
    ax2.set_xticks(np.arange(-500, 750, 250))
    ax2.set_yticks(np.arange(-500, 750, 250))
    ax2.tick_params(axis='both', labelsize=32)

    cax = f.add_axes([.825, .2, 0.0125, .6])
    cb = f.colorbar(cf, cax=cax)
    for l in cb.ax.yaxis.get_ticklabels():
        l.set_fontsize(24)

    out_file_name = plotpath + "/" + field + "_xz_xy_" + timestep + ".png"
    f.savefig(out_file_name, bbox_inches='tight')

if __name__ == "__main__":

    try:
        config, datapath, nx, ny, fields, timesteps, plotpath, = sys.argv[1:8]
    except ValueError:
        print("Usage: {0} <file_stem_name> <datapath> <nx> <ny> <field_names>"
              "<timestep_list> <plotpath>".format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    # Split out the list of timesteps
    ts_list = timesteps.split(':')

    for field in field_list:

        for ts in ts_list:

            filestem = datapath + "/" + config + "_nodal_" + \
                field + "_" + ts + "*"
            data = read_nodal_data(filestem, 1, 1)

            # Sort the data (needed to be able to reshape and not regrid)
            data = data.sort_values(['y', 'x', 'z'])

            levels = np.sort(data.level.unique())

            # Only try to plot if we found some files for this timestep
            if len(levels) > 0:
                make_figure(plotpath, nx, ny, field, 1, ts)
