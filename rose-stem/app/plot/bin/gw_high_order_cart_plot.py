#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2018 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution
##############################################################################

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.interpolate import griddata
import math
import sys

from read_data import read_nodal_data

levels = None
data = None


def make_figure(datapath, plotpath, field, component, timestep, nx, ny):

    fig = plt.figure(figsize=(10, 5.25))
    # set fontsize
    fsize = 28

    val_col = 'c' + str(component)

    # Sort levels in asscending order, this is needed for high order spaces
    sorted_levels = sorted(levels)
    l2h = np.zeros(len(levels))
    for i in range(len(levels)):
        for j in range(len(levels)):
            if (sorted_levels[i] == levels[j]):
                l2h[i] = j

    # Get min and max of x,y data for plot axes
    min_lev = min(levels)
    max_lev = max(levels)
    xmin = data.loc[data['level'] == min_lev]['x'].min()
    xmax = data.loc[data['level'] == min_lev]['x'].max()
    ymin = data.loc[data['level'] == min_lev]['y'].min()
    ymax = data.loc[data['level'] == min_lev]['y'].max()
    zmin = data.loc[data['level'] == min_lev]['z'].min()
    zmax = data.loc[data['level'] == max_lev]['z'].max()

    r2d = 1.0/1000.0
    nz = len(levels)

    # Create 2D plot
    x2d = np.linspace(xmin, xmax, nx)
    y2d = np.linspace(ymin, ymax, ny)
    z2d = np.linspace(zmin, zmax, nz)
    wi = np.zeros([ny, nx, len(levels)])
    xi, yi = np.meshgrid(x2d, y2d)
    for p in range(len(levels)):
        pp = int(l2h[p])
        p_data = data.loc[data['level'] == levels[pp]]
        wi[:, :, p] = griddata((p_data['x'].values, p_data['z'].values),
                               p_data[val_col].values, (xi, yi),
                               method='nearest')

    dw = np.zeros([len(levels), nx])
    for i in range(nx):
        dw[:, i] = wi[0, i, :] - wi[0, 0, :]

    ax = fig.add_subplot(1, 1, 1)

    cc = np.linspace(-0.0030, 0.0030, 13)
    c_map = cm.summer
    xi, zi = np.meshgrid(x2d, z2d)
    cf = ax.contourf(xi * r2d, zi * r2d, dw, cc, cmap=c_map)
    cl = ax.contour(xi * r2d, zi*r2d, dw, cc, linewidths=0.5, colors='k')
    ax.tick_params(axis='both', which='major', labelsize=fsize)
    ax.set_xlabel('x (km)', fontsize=fsize)
    ax.set_ylabel('z (km)', fontsize=fsize)
    ax.set_xlim([-150.0, 150.0])
    ax.set_ylim([0.0, 10.0])

    plt.tight_layout()
    out_file_name = plotpath + "/" + 'gravity_wave' + "_" + field + "_" + \
        timestep + ".png"
    plt.savefig(out_file_name, bbox_inches='tight')


if __name__ == "__main__":

    try:
        config, datapath, nx, ny, fields, timesteps, plotpath = sys.argv[1:8]
    except ValueError:
        print("Usage: {0} <file_stem_name> <datapath>, <nx>, <ny>, "
              "<fields>, <timesteps>, <plotpath> ".format(sys.argv[0]))
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
            filestem = datapath + "/" + config + "_nodal_" + field + "_" + \
                ts + "*"

            if (field != 'u' and field != 'xi'):
                data = read_nodal_data(filestem, ncomp, comp)
                if (not data.empty):
                    levels = data.level.unique()
                    make_figure(datapath, plotpath, field, comp, ts,
                                int(nx), int(ny))
                    any_plots = True
            else:
                for comp_u in comp:
                    data = read_nodal_data(filestem, ncomp, comp_u)
                    if (not data.empty):
                        levels = data.level.unique()
                        make_figure(datapath, plotpath, field, comp_u, ts,
                                    int(nx), int(ny))
                        any_plots = True

    if not any_plots:
        print("Error: No plots made.")
        exit(2)
