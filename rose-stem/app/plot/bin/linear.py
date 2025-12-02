#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Plot horizontal and vertical cross sections of desired fields with
fixed contour intervals
'''
import sys

import numpy as np

# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')  # noqa: E402
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from read_data import read_ugrid_data
from scipy.interpolate import griddata

# Use summer colormap
C_MAP = cm.summer

# Size of regular grid
NY, NX = 360, 720


def interpolate_data(cube, n_levs, time_value):
    '''
    Interpolate the data
    '''

    # Compute the horizontal grid
    x_coord = np.around(cube.coord('longitude').points, decimals=5)
    y_coord = np.around(cube.coord('latitude').points, decimals=5)

    xmin = np.amin(x_coord)
    xmax = np.amax(x_coord)
    ymin = np.amin(y_coord)
    ymax = np.amax(y_coord)

    # Generate a regular grid to interpolate the data.
    xi_grid = np.linspace(xmin, xmax, NX)
    yi_grid = np.linspace(ymin, ymax, NY)

    x_mesh, y_mesh = np.meshgrid(xi_grid, yi_grid)

    # Interpolate using delaunay triangularization
    plot_data = np.zeros((NY, NX, n_levs))

    for _, lev in enumerate(range(n_levs)):

        data = cube.data[time_value, lev]
        linear_data = griddata((x_coord, y_coord), data,
                               (x_mesh, y_mesh), method='linear')
        nearest_data = griddata((x_coord, y_coord), data,
                                (x_mesh, y_mesh), method='nearest')
        linear_data[np.isnan(linear_data)] = (nearest_data[
            np.isnan(linear_data)])

        plot_data[:, :, lev] = linear_data

    return(plot_data, xi_grid, yi_grid)


def make_figures(filein, plotpath, fields, lid, figname, idx_list):
    '''
    Read, interpolate and plot the data
    '''

    # assume uniform grid
    zmin = 0.0
    zmax = lid
    zi_full = np.linspace(zmin, zmax, lid+1)
    zi_half = 0.5*(zi_full[1:] + zi_full[0:lid])

    directions = ['xy', 'yz', 'xz']

    for time_value in [-1]:

        nyplots = len(directions)
        nxplots = 1

        for _, field in enumerate(fields):

            interp_fig = plt.figure(figsize=(20, 10))
            fields_name = field
            cube = read_ugrid_data(filein, field)

            # Set some levels for contours:
            levels = None
            if figname == 'linear_model-semi-implicit':
                if field == 'rho':
                    levels = np.linspace(-0.001, 0.0001, 11)
            if figname == 'linear_model-runge-kutta':
                if field == 'rho':
                    levels = np.linspace(-0.002, 0.0001, 11)
            if figname == 'linear_model-dcmip':
                if field == 'theta':
                    levels = np.linspace(-0.12, 0.12, 13)
                if field == 'u_in_w2h':
                    levels = np.linspace(-0.4, 0.4, 11)
                if field == 'v_in_w2h':
                    levels = np.linspace(-0.2, 0.2, 11)
                if field == 'w_in_wth':
                    levels = np.linspace(-0.2, 0.2, 11)
                if field == 'exner':
                    levels = np.linspace(-0.00003, 0.00003, 11)
                if field == 'rho':
                    levels = np.linspace(-0.0007, 0.0002, 11)
            if figname == 'linear_model-nwp-gal9_C12':
                if field == 'theta':
                    levels = np.linspace(-5.0, 5.0, 10)
                if field == 'u_in_w2h':
                    levels = np.linspace(-6.0, 6.0, 10)
                if field == 'v_in_w2h':
                    levels = np.linspace(-6.0, 6.0, 10)
                if field == 'w_in_wth':
                    levels = np.linspace(-0.05, 0.05, 10)
                if field == 'exner':
                    levels = np.linspace(-0.003, 0.003, 10)
                if field == 'rho':
                    levels = np.linspace(-0.005, 0.005, 10)
                if field == 'm_v':
                    levels = np.linspace(-0.0005, 0.0005, 10)
                if field == 'm_r':
                    levels = np.linspace(-0.0001, 0.0001, 10)
                if field == 'm_cl':
                    levels = np.linspace(-0.0001, 0.0001, 10)
                if field == 'm_s':
                    levels = np.linspace(-0.0001, 0.0001, 10)
            if figname == 'linear_model-nwp-gal9_C224':
                if field == 'theta':
                    levels = np.linspace(-6.0, 6.0, 6)
                if field == 'u_in_w2h':
                    levels = np.linspace(-12.0, 12.0, 6)
                if field == 'v_in_w2h':
                    levels = np.linspace(-12.0, 12.0, 6)
                if field == 'w_in_wth':
                    levels = np.linspace(-0.1, 0.1, 6)
                if field == 'exner':
                    levels = np.linspace(-0.003, 0.003, 10)
                if field == 'rho':
                    levels = np.linspace(-0.05, 0.05, 6)
                if field == 'm_v':
                    levels = np.linspace(-0.0005, 0.0005, 6)
                if field == 'm_r':
                    levels = np.linspace(-0.0001, 0.0001, 6)
                if field == 'm_cl':
                    levels = np.linspace(-0.0001, 0.0001, 6)
                if field == 'm_s':
                    levels = np.linspace(-0.0001, 0.0001, 6)

            # Vertical levels will be last entry in dimension coords
            levels_name = cube.dim_coords[-1].name()
            n_levs = len(cube.coord(levels_name).points)

            time = np.around(cube.coord('time').points, decimals=1)

            (plot_data, xi, yi) = interpolate_data(cube, n_levs, time_value)

            # Choose the correct vertical level set
            if lid == n_levs-1:
                zi_plot = zi_full
            else:
                zi_plot = zi_half

            # idx_list contains degrees/level to plot
            plot_long = (int(idx_list[0])+180)*int(NX/360)
            plot_lat = (int(idx_list[1])+90)*int(NY/180)
            plot_level = int(idx_list[2])

            for d_iterate, direction in enumerate(directions):

                plotnum = d_iterate + 1

                interp_fig.add_subplot(nyplots, nxplots, plotnum)

                if direction == 'xz':
                    h_x1, h_x2 = np.meshgrid(xi, zi_plot)
                    h_x3 = plot_data[plot_lat, :, :].T
                    plt.title([field, direction, ' lat = ',
                               yi[plot_lat]*360./np.real(NX)])
                if direction == 'yz':
                    h_x1, h_x2 = np.meshgrid(yi, zi_plot)
                    h_x3 = plot_data[:, plot_long, :].T
                    plt.title([field, direction, 'long = ',
                               xi[plot_long]*360./np.real(NX)])
                if direction == 'xy':
                    h_x2, h_x1 = np.meshgrid(yi, xi)
                    h_x3 = plot_data[:, :, plot_level].T
                    plt.title([field, direction,
                               'Height = ', zi_plot[plot_level]])

                plt.contourf(h_x1, h_x2, h_x3, levels=levels, cmap=C_MAP)
                plt.colorbar(cmap=C_MAP)
                plt.contour(h_x1, h_x2, h_x3,
                            levels=levels, linewidths=0.5, colors='k')

            pngfile = '%s/%s-%s-time%s.png' % (plotpath,
                                               figname,
                                               fields_name,
                                               time[time_value])

            plt.tight_layout()
            plt.savefig(pngfile)

        plt.close()


if __name__ == "__main__":

    try:
        ARGS = sys.argv[:]
        FILES, PLOTPATH, LID, FIGNAME = ARGS[1:5]
        FIELD_LIST = ARGS[5].split(':')
        IDX_LIST = ARGS[6].split(':')
    except ValueError:
        print("Usage: {0} <filein> <plotpath> <lid>"
              " <figname> <fields_list> <idx_list>"
              .format(sys.argv[0]))
        sys.exit(1)

    FILE_LIST = FILES.split(':')
    for FILEIN in FILE_LIST:
        make_figures(FILEIN, PLOTPATH, FIELD_LIST,
                     int(LID), FIGNAME, IDX_LIST)
