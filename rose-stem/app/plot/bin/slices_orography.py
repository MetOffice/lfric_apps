#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2019 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
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
from six.moves import range
matplotlib.use('Agg')  # noqa: E402

import matplotlib.pyplot as plt

from scipy.interpolate import griddata
import sys

from read_data import read_nodal_data

# Use viridis colormap
from python_maps import viridis_data
from matplotlib.colors import ListedColormap

viridis = ListedColormap(viridis_data, name='viridis')
plt.register_cmap(name='viridis', cmap=viridis)

levels = None
data = None


def make_figure(plotpath, field, component, timestep, plotlong, plotlat,
                plotlevel, test):
    # Sort levels in asscending order, this is needed for high order spaces
    sorted_levels = sorted(levels)
    l2h = np.zeros(len(levels))
    for i, l1 in enumerate(levels):
        for j, l2 in enumerate(levels):
            if sorted_levels[i] == levels[j]:
                l2h[i] = j

    # Get min and max of x,y data for plot axes
    min_lev = min(levels)

    xmin = data.loc[data['level'] == min_lev]['x'].min()
    xmax = data.loc[data['level'] == min_lev]['x'].max()
    ymin = data.loc[data['level'] == min_lev]['y'].min()
    ymax = data.loc[data['level'] == min_lev]['y'].max()

    zmin = 0.0
    zmax = 30000.0

    r2d = 180.0/np.pi
    nx, ny, nz = 360, 180, len(levels)

    val_col = 'c' + str(component)

    # Create 2D plot
    x2d = np.linspace(xmin, xmax, nx)
    z2d = np.linspace(zmin, zmax, nz)
    y2d = np.linspace(ymin, ymax, ny)
    zi = np.zeros([ny, nx, len(levels)])
    height_i = np.zeros([ny, nx, len(levels)])
    xi, yi = np.meshgrid(x2d, y2d)
    for p, l in enumerate(levels):
        p_data = data.loc[data['level'] == levels[p]]
        zi[:, :, p] = griddata((p_data['x'].values, p_data['y'].values),
                               p_data[val_col].values, (xi, yi),
                               method='linear')
        height_i[:, :, p] = griddata((p_data['x'].values, p_data['y'].values),
                                     p_data['z'].values, (xi, yi),
                                     method='linear')

    if test == '2-1':  # dcmip21 test
        cc = np.concatenate((np.linspace(-1.6, -0.2, 8),
                             np.linspace(0.2, 1.6, 8)), axis=0)
    elif test == '2-2':  # dcmip22 test
        cc = np.concatenate((np.linspace(-3.0, -0.2, 15),
                             np.linspace(0.2, 3.0, 15)), axis=0)
    elif test == '2-0-0':  # dcmip200 test
        if component == 3:
            cc = np.linspace(-0.0011, 0.0011, 12)
        else:
            cc = np.linspace(-0.11, 0.11, 12)
    else:
        cc = np.concatenate((np.linspace(-3.0, -0.2, 15),
                             np.linspace(0.2, 3.0, 15)), axis=0)

    c_map = viridis

    # xz plot
    if int(plotlat) >= -90 and int(plotlat) <= 90:

        lat = int(plotlat)+90

        yi, xi = np.meshgrid(z2d, x2d)
        dz = np.zeros([nx, len(levels)])
        for i in range(nx):
            dz[i, :] = zi[lat, i, :]

        # Take actual orography heights as z
        zi_adj = height_i[lat, :, :]

        fig = plt.figure(figsize=(10, 5))
        cf = plt.contourf(xi*r2d, zi_adj, dz, cc, cmap=c_map)
        plt.colorbar(cf, cmap=c_map)
        cl = plt.contour(xi*r2d, zi_adj, dz, cc, linewidths=0.5, colors='k')
        plt.title('max: %e, min: %e' % (np.max(dz), np.min(dz)))
        plt.xlabel('Longitude')
        plt.ylabel('z')
        if (field != 'u' and field != 'xi'):
            out_file_name = (plotpath + "/slice_xz_" + field +
                             "_" + timestep + ".png")
        else:
            out_file_name = (plotpath + "/slice_xz_" + field +
                             str(component) + "_" + timestep + ".png")
        plt.savefig(out_file_name, bbox_inches='tight')

    # xy plot
    if int(plotlevel) >= 0 and int(plotlevel) < nz:
        fig = plt.figure(figsize=(10, 5))
        xi, yi = np.meshgrid(x2d, y2d)
        dz = zi[:, :, int(plotlevel)]
        cf = plt.contourf(xi*r2d, yi*r2d, dz, cc, cmap=c_map)
        plt.colorbar(cf, cmap=c_map)
        cl = plt.contour(xi * r2d, yi * r2d, dz, cc, linewidths=0.5,
                         colors='k')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        if (field != 'u' and field != 'xi'):
            out_file_name = (plotpath + "/slice_xy_" + field + "_" +
                             timestep + ".png")
        else:
            out_file_name = (plotpath + "/slice_xy_" + field +
                             str(component) + "_" + timestep + ".png")
        plt.savefig(out_file_name, bbox_inches='tight')


if __name__ == "__main__":

    try:
        (config, datapath, fields, timesteps, plotlong, plotlat, plotlevel,
         plotpath, test) = sys.argv[1:10]
    except ValueError:
        print("Usage: {0} <file_stem_name>  <datapath> <field_names> "
              "<timestep_list> <plotlong> <plotlat> <plotlevel> "
              "<plotpath> <test>".format(sys.argv[0]))

        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    # Split out the list of timesteps
    ts_list = timesteps.split(':')

    for field in field_list:

        if field in ['rho', 'theta', 'exner', 'buoyancy']:
            # Scalar fields
            ncomp = 1
            comp = 1
        else:
            # Vector fields
            ncomp = 3
            # W3 projected U, V, W and XI components
            if (field == 'u' or field == 'xi'):
                comp = [1, 2, 3]
            else:
                comp = 1

        for ts in ts_list:

            filestem = (datapath + "/" + config + "_nodal_" + field + "_" +
                        ts + "*")
            if(field != 'u' and field != 'xi'):
                data = read_nodal_data(filestem, ncomp, comp)
                levels = data.level.unique()

                # Only try to plot if we found some files for this timestep
                if len(levels) > 0:
                    make_figure(plotpath, field, comp, ts, plotlong, plotlat,
                                plotlevel, test)

            else:
                for comp_u in comp:
                    data = read_nodal_data(filestem, ncomp, comp_u)
                    levels = data.level.unique()

                    # Only try to plot if we found some files for this timestep
                    if len(levels) > 0:
                        make_figure(plotpath, field, comp_u, ts, plotlong,
                                    plotlat, plotlevel, test)
