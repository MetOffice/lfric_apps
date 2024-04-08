#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Python script to plot xz slices along the y=0 and xy slices on a specified
level

This version takes nodal format output files and
interpolates onto a regular grid.

Filename hardcoded.

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
matplotlib.use('Agg')  # noqa: E402

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.interpolate import griddata
import math
import sys

from read_data import read_nodal_data

levels = None
data = None


def make_figure(plotpath, field, component, timestep,
                plotlong, plotlat, plotlevel):
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
    xi, yi = np.meshgrid(x2d, y2d)
    for p in range(len(levels)):
        pp = int(l2h[p])
        p_data = data.loc[data['level'] == levels[p]]
        zi[:, :, p] = griddata((p_data['x'].values, p_data['y'].values),
                               p_data[val_col].values, (xi, yi),
                               method='linear')

    if field == 'constant':
        cc = np.linspace(0.85,1.15, 13)
    else:
        cc = np.linspace(np.amin(data[val_col].values),
                         np.amax(data[val_col].values), 13)
    c_map = cm.summer

    # xz plot
    if int(plotlat) >= -90 and int(plotlat) <= 90:
        lat = int(plotlat)+90
        yi, xi = np.meshgrid(z2d, x2d)
        dz = np.zeros([nx, len(levels)])
        for i in range(nx):
            dz[i, :] = zi[lat, i, :]

        fig = plt.figure(figsize=(10, 5))
        cf = plt.contourf(xi * r2d, yi / 1000.0, dz, cc, cmap=c_map)
        plt.colorbar(cf, cmap=c_map)
        cl = plt.contour(xi * r2d, yi / 1000.0, dz, cc, linewidths=0.5,
                         colors='k')
        plt.title('max: %e, min: %e' % (np.max(dz), np.min(dz)))
        plt.xlabel('Longitude')
        plt.ylabel('z')
        if (field != 'u' and field != 'xi' and field != 'wind'):
            out_file_name = (plotpath + "/" "slice_xz_" + field + "_"
                             + timestep + ".png")
        else:
            out_file_name = (plotpath + "/" "slice_xz_" + field
                             + str(component) + "_" + timestep + ".png")
        plt.savefig(out_file_name, bbox_inches='tight')
        plt.close()

    # yz plot
    if int(plotlong) >= 0 and int(plotlong) <= 360:
        yi, xi = np.meshgrid(z2d, y2d)
        dz = np.zeros([ny, len(levels)])
        for i in range(ny):
            dz[i, :] = zi[i, int(plotlong), :]

        fig = plt.figure(figsize=(10, 5))
        cf = plt.contourf(xi * r2d, yi / 1000.0, dz, cc, cmap=c_map)
        plt.colorbar(cf,  cmap=c_map)
        cl = plt.contour(xi * r2d, yi / 1000.0, dz, cc, linewidths=0.5,
                         colors='k')
        plt.title('max: %e, min: %e' % (np.max(dz), np.min(dz)))
        plt.xlabel('Latitude')
        plt.ylabel('z')
        if (field != 'u' and field != 'xi' and field != 'wind'):
            out_file_name = (plotpath + "/" "slice_yz_" + field + "_"
                             + timestep + ".png")
        else:
            out_file_name = (plotpath + "/" "slice_yz_" + field
                             + str(component) + "_" + timestep + ".png")
        plt.savefig(out_file_name, bbox_inches='tight')
        plt.close()

    # xy plot
    if int(plotlevel) >= 0 and int(plotlevel) < nz:
        fig = plt.figure(figsize=(10, 5))
        xi, yi = np.meshgrid(x2d, y2d)
        dz = zi[:, :, int(plotlevel)]
        cf = plt.contourf(xi * r2d, yi * r2d, dz, cc, cmap=c_map)
        plt.colorbar(cf,  cmap=c_map)
        cl = plt.contour(xi * r2d, yi * r2d, dz, cc, linewidths=0.5,
                         colors='k')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        if (field != 'u' and field != 'xi' and field != 'wind'):
            out_file_name = (plotpath + "/" "slice_xy_" + field + "_"
                             + timestep + ".png")
        else:
            out_file_name = (plotpath + "/" "slice_xy_" + field
                             + str(component) + "_" + timestep + ".png")
        plt.savefig(out_file_name, bbox_inches='tight')
        plt.close()


if __name__ == "__main__":

    try:
        config, datapath, fields, timesteps, plotlong, plotlat, plotlevel, \
            plotpath = sys.argv[1:9]
    except ValueError:
        print("Usage: {0} <file_stem_name> <datapath> <field_names> \
            <timestep_list> <plotlong> <plotlat> <plotlevel> \
            <plotpath>".format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    # Split out the list of timesteps
    ts_list = timesteps.split(':')

    for field in field_list:

        if field in ['rho', 'theta', 'exner', 'buoyancy', 'density',
                     'pressure', 'm_v', 'tracer_con', 'tracer_adv',
                     'constant']:
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
            elif (field == 'u' or field == 'xi' or field == 'wind'):
                comp = [1, 2, 3]

        for ts in ts_list:

            filestem = (datapath + "/" + config + "_nodal_" + field + "_"
                        + ts + "*")

            if (field != 'u' and field != 'xi' and field != 'wind'):
                data = read_nodal_data(filestem, ncomp, comp)

                levels = data.level.unique()

                # Only try to plot if we found some files for this timestep
                if len(levels) > 0:
                    make_figure(plotpath, field, comp, ts, plotlong,
                                plotlat, plotlevel)

            else:
                for comp_u in comp:
                    data = read_nodal_data(filestem, ncomp, comp_u)

                    levels = data.level.unique()

                    # Only try to plot if we found some files for this timestep
                    if len(levels) > 0:
                        make_figure(plotpath, field, comp_u, ts, plotlong,
                                    plotlat, plotlevel)
