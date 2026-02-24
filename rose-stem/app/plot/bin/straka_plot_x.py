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

import sys

from read_data import read_nodal_data

# Use viridis colormap
try:
    viridis = plt.get_cmap("viridis")
except ValueError:
    from matplotlib.colors import ListedColormap
    from python_maps import viridis_data
    viridis = ListedColormap(viridis_data, name='viridis')
    plt.register_cmap(name='viridis', cmap=viridis)

levels = None
data = None


def make_figure(plotpath, nx, ny, field, component, timestep):

    val_col = 'c' + str(component)

    slice_fig, ax = plt.subplots(figsize=(12, 5))

    # get min and max of x,y data for plot axes
    min_lev = min(levels)

    xmin = data.loc[data['level'] == min_lev]['x'].min()
    xmax = data.loc[data['level'] == min_lev]['x'].max()
    ymin = data.loc[data['level'] == min_lev]['y'].min()
    ymax = data.loc[data['level'] == min_lev]['y'].max()

    zmin = 0.0
    zmax = 6400.0

    r2d = 1.0/1000.0

    nx = int(nx)
    ny = int(ny)
    nz = len(levels)

    zi = np.zeros([ny, nx, len(levels)])

    for p in range(len(levels)):
        p_data = data.loc[data['level'] == levels[p]]
        zi[:, :, p] = (p_data[val_col].values).reshape((ny, nx))

    c_map = viridis

    # create meshgrid to get x_i and y_i for plotting
    x2d = np.linspace(xmin, xmax, nx)
    z2d = np.linspace(zmin, zmax, nz)
    y_i, x_i = np.meshgrid(z2d, x2d)

    dz = np.zeros([nx, len(levels)])
    back = 0.0
    if field == 'theta':
        back = 300.0
        cc = np.linspace(-16, -1, 16)
    elif field == 'w3projection_xi2':
        cc = np.linspace(-0.06, 0.06, 11)
    else:
        cc = np.linspace(np.amin(zi), np.amax(zi), 11)

    for i in range(nx):
        dz[i, :] = zi[0, i, :] - back

    cf = ax.contourf(x_i * r2d, y_i * r2d, np.round(dz, 10),
                     cc, cmap=c_map, extend='min')
    # Add contour lines -- solid for negative
    _ = ax.contour(x_i * r2d, y_i * r2d, np.round(dz, 10), cc,
                   colors='k', linewidths=0.5, linestyles='solid')

    ax.set_xlabel("x (km)", fontsize=28)
    ax.set_ylabel("z (km)", fontsize=28)
    ax.set_xlim([0.0, 16.0])
    ax.set_ylim([0.0, 5.0])
    ax.set_xticks(np.arange(0, 18, 2))
    ax.set_yticks(np.arange(0, 6, 1))
    # Set tick label size
    ax.tick_params(axis='both', which='major', labelsize=24)
    cb = plt.colorbar(cf,  cmap=c_map, pad=0.04)
    cb.ax.tick_params(labelsize=24)

    # Front position is symmetrized
    print(nx, np.min(dz), np.max(dz))

    out_file_name = plotpath + "/" + 'straka_x_' + field + "_" + timestep \
        + ".png"
    slice_fig.savefig(out_file_name, bbox_inches='tight')


if __name__ == "__main__":

    try:
        config, datapath, nx, ny, fields, timesteps, plotpath = sys.argv[1:8]
    except ValueError:
        print("Usage: {0} <file_stem_name> <datapath> <nx> <ny> "
              "<field_names> <timestep_list> "
              "<plotpath>".format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    # Split out the list of timesteps
    ts_list = timesteps.split(':')

    any_plots = False

    for field in field_list:

        if field in ['rho', 'theta', 'exner', 'buoyancy']:
            # Scalar fields
            ncomp = 1
            comp = 1
        else:
            # Vector fields
            ncomp = 3
            # W3 projected U, V, W and XI components
            if field in ['w3projection_u1', 'w3projection_u2',
                         'w3projection_u3', 'w3projection_xi1',
                         'w3projection_xi2', 'w3projection_xi3']:
                comp = 1
            elif (field == 'u' or field == 'xi'):
                comp = [1, 2, 3]

        for ts in ts_list:

            filestem = datapath + "/" + config + "_nodal_" + \
                field + "_" + ts + "*"

            if (field != 'u' and field != 'xi'):
                data = read_nodal_data(filestem, ncomp, comp)
                if (not data.empty):
                    # Sort the data
                    # (needed to be able to reshape and not regrid)
                    data = data.sort_values(['y', 'x', 'z'])
                    levels = np.sort(data.level.unique())
                    make_figure(plotpath, nx, ny, field, comp, ts)
                    any_plots = True
            else:
                for comp_u in comp:
                    data = read_nodal_data(filestem, ncomp, comp_u)
                    if (not data.empty):
                        # Sort the data
                        # (needed to be able to reshape and not regrid)
                        data = data.sort_values(['y', 'x', 'z'])
                        levels = np.sort(data.level.unique())
                        make_figure(plotpath, nx, ny, field, comp_u, ts)
                        any_plots = True

    if not any_plots:
        print("Error: No plots made.")
        exit(2)
