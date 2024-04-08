#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2018 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Python script to plot slize of the difference of a field
from the initial timelevel 0 value

This version takes nodal format output files and
interpolates onto a regular grid.

Filename hardcoded.

Levels are determined from the data

This version stitches together a directory of files
and extracts all levels so it can work in the serial
case where there is one file or the parallel case where
there is a file for each processor.

'''

import numpy as np
from scipy.interpolate import griddata
import math
import sys
from read_data import read_nodal_data
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm


levels = None
data = None
data0 = None


def make_figure(plotpath, field, component, timestep, plotlong, plotlat,
                plotlevel, err_range):
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

    zmin = min_lev*1000.0
    zmax = max_lev*1000.0

    r2d = 180.0/np.pi
    nx, ny, nz = 360, 180, len(levels)

    val_col = 'c' + str(component)

    # Create 2D plot
    x2d = np.linspace(xmin, xmax, nx)
    z2d = np.linspace(zmin, zmax, nz)
    y2d = np.linspace(ymin, ymax, ny)
    zi = np.zeros([ny, nx, len(levels)])
    zi = np.zeros([ny, nx, len(levels)])
    zp = np.zeros([ny, nx, len(levels)])
    xi, yi = np.meshgrid(x2d, y2d)
    for p in range(len(levels)):
        pp = int(l2h[p])
        p_data = data.loc[data['level'] == levels[p]]
        zi[:, :, p] = griddata((p_data['x'].values, p_data['y'].values),
                               p_data[val_col].values, (xi, yi),
                               method='linear')
        zp[:, :, p] = griddata((p_data['x'].values, p_data['y'].values),
                               p_data[val_col].values, (xi, yi),
                               method='linear')

        # Subtract off initial value
        p0_data = data0.loc[data0['level'] == levels[p]]
        zi[:, :, p] = zi[:, :, p] \
            - griddata((p0_data['x'].values, p0_data['y'].values),
                       p0_data[val_col].values, (xi, yi), method='linear')

    e_min = float(err_range)
    cc = np.linspace(-e_min, e_min, 13)
    ccp = np.linspace(-0.2, 1.2, 15)
    c_map = cm.summer

    # xz plot
    if int(plotlat) >= -90 and int(plotlat) <= 90:

        lat = int(plotlat)+90

        yi, xi = np.meshgrid(z2d, x2d)
        dz = np.zeros([nx, len(levels)])
        for i in range(nx):
            dz[i, :] = zi[lat, i, :]

        fig = plt.figure(figsize=(10, 5))
        cf = plt.contourf(xi*r2d, yi/1000.0, dz, cc, cmap=c_map)
        plt.colorbar(cf,  cmap=c_map)
        cl = plt.contour(xi * r2d, yi / 1000.0, dz, cc, linewidths=0.5,
                         colors='k')
        plt.title('max: %e, min: %e' % (np.max(dz), np.min(dz)))
        plt.xlabel('Longitude')
        plt.ylabel('z')
        if (field != 'u' and field != 'xi'):
            out_file_name = plotpath + "/" "slice_xz_" + field + "_" + timestep
            + ".png"
        else:
            out_file_name = plotpath + "/" "slice_xz_" + field + str(component)
            + "_" + timestep + ".png"
        plt.savefig(out_file_name, bbox_inches='tight')

    # yz plot
    if int(plotlong) >= 0 and int(plotlong) <= 360:

        yi, xi = np.meshgrid(z2d, y2d)
        dz = np.zeros([ny, len(levels)])
        for i in range(ny):
            dz[i, :] = zi[i, int(plotlong), :]

        fig = plt.figure(figsize=(10, 5))
        cf = plt.contourf(xi*r2d, yi/1000.0, dz, cc, cmap=c_map)
        plt.colorbar(cf,  cmap=c_map)
        cl = plt.contour(xi * r2d, yi / 1000.0, dz, cc, linewidths=0.5,
                         colors='k')
        plt.title('max: %e, min: %e' % (np.max(dz), np.min(dz)))
        plt.xlabel('Latitude')
        plt.ylabel('z')
        if (field != 'u' and field != 'xi'):
            out_file_name = plotpath + "/" "slice_yz_" + field + "_"
            + timestep + ".png"
        else:
            out_file_name = plotpath + "/" "slice_yz_" + field + str(component)
            + "_" + timestep + ".png"
        plt.savefig(out_file_name, bbox_inches='tight')

    # xy plot
    if int(plotlevel) >= 0 and int(plotlevel) < nz:
        slice_fig = plt.figure(figsize=(10, 12))
        ax1 = slice_fig.add_subplot(2, 1, 1)
        xi, yi = np.meshgrid(x2d, y2d)
        dz = zi[:, :, int(plotlevel)]
        cf = ax1.contourf(xi*r2d, yi*r2d, dz, cc, cmap=c_map)
        plt.colorbar(cf, cmap=c_map)
        cl = ax1.contour(xi * r2d, yi * r2d, dz, cc, linewidths=0.5,
                         colors='k')
        ax1.set_xlabel('Longitude')
        ax1.set_ylabel('Latitude')
        ax1.set_title('Error')

        ax2 = slice_fig.add_subplot(2, 1, 2)
        xi, yi = np.meshgrid(x2d, y2d)
        dz = zp[:, :, int(plotlevel)]
        cf = ax2.contourf(xi*r2d, yi*r2d, dz, ccp, cmap=c_map)
        plt.colorbar(cf, cmap=c_map)
        cl = ax2.contour(xi * r2d, yi * r2d, dz, ccp, linewidths=0.5,
                         colors='k')
        ax2.set_xlabel('Longitude')
        ax2.set_ylabel('Latitude')
        ax2.set_title('Field')
        if (field != 'u' and field != 'xi'):
            out_file_name = plotpath + "/" "slice_xy_" + field + "_" \
                + timestep + ".png"
        else:
            out_file_name = plotpath + "/" "slice_xy_" + field \
                + str(component) + "_" + timestep + ".png"
        plt.savefig(out_file_name, bbox_inches='tight')


if __name__ == "__main__":
    try:
        config, datapath, fields, timesteps, plotlong, plotlat, plotlevel, \
            err_range, plotpath = sys.argv[1:10]
    except ValueError:
        print("Usage: {0} <config> <datapath> <field_names> <timestep_list>"
              "<plotlong> <plotlat> <plotlevel> <err_range> <plotpath>"
              .format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    # Split out the list of timesteps
    ts_list = timesteps.split(':')

    any_plots = False

    for field in field_list:

        # Scalar fields only
        ncomp = 1
        comp = 1

        filestem = datapath + "/" + config + "_nodal_" + field + "_" \
            + "T000000" + "*"
        data0 = read_nodal_data(filestem, ncomp, comp)

        for ts in ts_list:

            filestem = datapath + "/" + config + "_nodal_" + field + "_" \
                + ts + "*"

            data = read_nodal_data(filestem, ncomp, comp)

            if (not data.empty):
                levels = data.level.unique()
                make_figure(plotpath, field, comp, ts, plotlong, plotlat,
                            plotlevel, err_range)
                any_plots = True

    if not any_plots:
        print("Error: No plots made.")
        exit(2)
