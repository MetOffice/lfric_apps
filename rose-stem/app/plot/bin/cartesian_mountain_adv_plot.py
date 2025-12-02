#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2018 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Basic python script to plot the x-z profile of variables for the orography
test cases on subplots.

This version takes nodal format output files and
interpolates onto a regular grid.

Levels are determined from the data

This version stitches together a directory of files
and extracts all levels so it can work in the serial
case where there is one file or the parallel case where
there is a file for each processor.

'''

import numpy as np
import sys
from read_data import read_nodal_data
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm

levels = None
data = None


def make_figure(plotpath, nx, ny, field, component, timestep, cntrs):

    slice_fig = plt.figure(figsize=(10, 15))

    val_col = 'c' + str(component)

    # get min and max of x,y data for plot axes
    min_lev = min(levels)

    xmin = data.loc[data['level'] == min_lev]['x'].min()
    xmax = data.loc[data['level'] == min_lev]['x'].max()
    ymin = data.loc[data['level'] == min_lev]['y'].min()
    ymax = data.loc[data['level'] == min_lev]['y'].max()

    # zmin is min of bottom level
    zmin = data['z'].min()

    zmax = 25000.0

    val_col = 'c' + str(component)

    nx = int(nx)
    ny = int(ny)
    nl = len(levels)

    # create 2D plot
    val_i = np.zeros([ny, nx, nl])
    x_i = np.zeros([ny, nx, nl])
    height_i = np.zeros([ny, nx, nl])

    # interpolate field values and heights onto xy for each level
    for p in range(len(levels)):
        p_data = data.loc[data['level'] == levels[p]]

        # Using reshape of numpy array
        val_i[:, :, p] = (p_data[val_col].values).reshape((ny, nx))
        height_i[:, :, p] = (p_data['z'].values).reshape((ny, nx))
        x_i[:, :, p] = (p_data['x'].values).reshape((ny, nx))/1000.0

    # Take actual orography heights as z
    zi_adj = height_i[0, :, :]/1000.0

    # extract the slice for plotting
    dval = np.zeros([nx, len(levels)])
    dval = val_i[1, :, :]

    x_plt = x_i[0, :, :]

    ax1 = slice_fig.add_subplot(2, 1, 1)

    # Setting contour limits and intervals
    cmin_t = 0.0
    cmax_t = 1.0
    nc_t = 11
    if cntrs == 'lines':
        cc = np.linspace(cmin_t, cmax_t, nc_t)
        cf = ax1.contour(x_plt, zi_adj, np.round(dval, 10), cc,
                         cmap=cm.coolwarm, linewidths=3)
        plt.colorbar(cf, cmap=cm.coolwarm)
    elif cntrs == 'colours':
        cf = ax1.pcolormesh(x_plt, zi_adj, np.round(dval, 10), vmin=cmin_t,
                            vmax=cmax_t, cmap=cm.coolwarm)
        plt.colorbar(cf, cmap=cm.coolwarm)

    ax1.set_title('max: %2.4e, min: %2.4e' % (np.max(dval), np.min(dval)))
    ax1.set_xlabel("x(m)", fontsize=24)
    ax1.set_ylabel("z(m)", fontsize=24)

    xmin_lim = -150.0
    xmax_lim = 150.0
    zmin = 0.0
    zmax = 25.0
    ax1.set_xticks(np.arange(xmin_lim, xmax_lim, 50))
    ax1.set_yticks(np.arange(zmin, zmax, 5))
    ax1.set_ylim(zmin, zmax)

    # Plot errors to analytic solution
    ax2 = slice_fig.add_subplot(2, 1, 2)

    error = np.zeros([nx, len(levels)])
    nsteps = int(timestep[1:])
    u0 = 10.0
    dt = 50.0
    x0 = -50000.0 + u0*dt*nsteps
    if x0 > 150000.0:
        x0 = x0 - 300000.0
    z0 = 9000.0
    r1 = 25000.0
    r2 = 3000.0
    for i in range(nx):
        for k in range(len(levels)):
            l1 = np.sqrt(((x_plt[i, k]*1000.0 - x0)/r1)**2 +
                         ((zi_adj[i, k]*1000.0 - z0)/r2)**2)
            if l1 < 1.0:
                analytic_soln = 1.0*np.cos(0.5*l1*np.pi)**2
            else:
                analytic_soln = 0.0
            error[i, k] = dval[i, k] - analytic_soln

    # Setting contour limits and intervals
    cmin_t = -0.3
    cmax_t = 0.3
    nc_t = 11
    if cntrs == 'lines':
        cc = np.linspace(cmin_t, cmax_t, nc_t)
        cf = ax2.contour(x_plt, zi_adj, np.round(error, 10), cc,
                         cmap=cm.coolwarm, linewidths=3)
        plt.colorbar(cf, cmap=cm.coolwarm)
    elif cntrs == 'colours':
        cf = ax2.pcolormesh(x_plt, zi_adj, np.round(error, 10), vmin=cmin_t,
                            vmax=cmax_t, cmap=cm.coolwarm)
        plt.colorbar(cf, cmap=cm.coolwarm)

    ax2.set_title('max: %2.4e, min: %2.4e' % (np.max(error), np.min(error)))
    ax2.set_xlabel("x(m)", fontsize=24)
    ax2.set_ylabel("z(m)", fontsize=24)

    xmin_lim = -150.0
    xmax_lim = 150.0
    zmin = 0.0
    zmax = 25.0
    ax2.set_xticks(np.arange(xmin_lim, xmax_lim, 50))
    ax2.set_yticks(np.arange(zmin, zmax, 5))
    ax2.set_ylim(zmin, zmax)

    # Save figure
    out_file_name = plotpath + "/" "BiP_mountain_" + field + '_' + ts + ".png"
    plt.tight_layout()
    plt.savefig(out_file_name)


if __name__ == "__main__":

    try:
        datapath, nx, ny, fields, component, timesteps, cntrs, plotpath \
            = sys.argv[1:9]
    except ValueError:
        print("Usage: {0} <datapath> <nx> <ny> <field_names> <component>"
              "<timestep_list> <cntrs> <plotpath>".format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    # Split out the list of timesteps
    ts_list = timesteps.split(':')

    any_plots = False

    for ts in ts_list:

        for field in field_list:

            filestem = datapath + "/diagGungho_nodal_" + field + "_" + ts + "*"

            if field in ['u']:
                data = read_nodal_data(filestem, 3, component)
            else:
                data = read_nodal_data(filestem, 1, component)

            if (not data.empty):
                # Sort the data (needed to be able to reshape and not regrid)
                data = data.sort_values(['y', 'x', 'z'])

                levels = np.sort(data.level.unique())

                # Only try to plot if we found some files for this timestep
                if len(levels) > 0:
                    make_figure(plotpath, nx, ny, field, component, ts, cntrs)
                    any_plots = True

    if not any_plots:
        print("Error: No plots made.")
        exit(2)
