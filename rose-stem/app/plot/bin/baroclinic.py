#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2019 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Need to set a non-interactive backend for suites
'''
Plots horizontal slices and vertical cross sections specified at particular
longitude and latitude values.
'''
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
# Need to set a non-interactive backend for suites
matplotlib.use('Agg')  # noqa: E402
from six.moves import range
import sys
import iris

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from read_data import read_ugrid_data
# Use magma colormap
try:
    magma = plt.get_cmap("magma")
except ValueError:
    from matplotlib.colors import ListedColormap
    from magma import magma_data
    magma = ListedColormap(magma_data, name='magma')
    plt.register_cmap(name='magma', cmap=cmap)

# Size of regular grid
ny, nx = 360, 720

plot_lon = 360
plot_lat = 280


def make_figures(filein, plotpath, fields, vertical_spacing, formulation):

    # Compute the vertical grid
    lid = 30.
    n_full = 31

    if vertical_spacing == 'um':
        # UM L38 set
        lid = 40.
        n_full = 39
        zi_f = np.array([
            .0, .0005095,  .0020380,  .0045854,  .0081519,  .0127373,
            .0183417,  .0249651,  .0326074,  .0412688,  .0509491,
            .0616485,  .0733668,  .0861040,  .0998603,  .1146356,
            .1304298,  .1472430,  .1650752,  .1839264,  .2037966,
            .2246857,  .2465938,  .2695209,  .2934670,  .3184321,
            .3444162,  .3714396,  .3998142,  .4298913,  .4620737,
            .4968308,  .5347160,  .5763897,  .6230643,  .6772068,
            .7443435,  .8383348, 1.000000])*lid

    elif vertical_spacing == 'dcmip':
        # DCMIP Stretched grid
        mu = 15.
        zi_f = np.zeros([n_full])
        for k in range(n_full):
            zi_f[k] = lid * \
                      (np.sqrt(mu * (float(k)/float(n_full))**2 + 1.) - 1.) \
                      / (np.sqrt(mu+1.) - 1)
    else:
        # Assume uniform grid
        zmin = 0.0
        zmax = lid
        zi_f = np.linspace(zmin, zmax, n_full)

    zi_h = 0.5*(zi_f[1:] + zi_f[0:n_full-1])

    direction = 'xy'  # , 'yz', 'xz'

    for t in [-1]:
        if fields is None:
            fields = ['surface_pressure_temperature']

        for field in fields:
            if field == 'surface_pressure_temperature':
                combined_fields = ['exner', 'theta']
            else:
                combined_fields = [field]

            interp_fig, ax = plt.subplots(figsize=(20, 10))

            for cfield in combined_fields:

                cube = read_ugrid_data(filein, cfield)
                # Vertical levels will be last entry in dimension coords
                levels_name = cube.dim_coords[-1].name()
                # Set some levels for contours:
                levels = None
                if cfield == 'air_potential_temperature' or cfield == 'theta':
                    levels = np.linspace(220, 330, 12)
                elif cfield == 'eastward_wind' or cfield == 'u_in_w2h':
                    if (formulation == 'moist'):
                        levels = np.arange(-54, 72, 6.)
                    else:
                        levels = np.arange(-20, 36, 4.)
                elif cfield == 'northward_wind' or cfield == 'v_in_w2h':
                    if (formulation == 'moist'):
                        levels = np.linspace(-50, 75, 14)
                    else:
                        levels = np.linspace(-25, 40, 14)
                elif cfield == 'upward_air_velocity' or cfield == 'w_in_wth':
                    if (formulation == 'moist'):
                        levels = np.linspace(-0.4, 0.4, 11)
                    else:
                        levels = np.linspace(-0.3, 0.3, 11)
                elif cfield == 'exner_pressure' or cfield == 'exner':
                    # Exner will be converted to hPa
                    levels = np.linspace(916, 1020, 14)
                elif cfield == 'air_density':
                    levels = np.arange(0, 1.4, 0.05)
                elif cfield == 'm_cl':
                    levels = np.arange(0.0, 7.5e-3, 0.5e-3)
                elif cfield == 'm_v':
                    levels = np.arange( 0.0, 5.0e-3, 2.0e-4)

                n_levs = len(cube.coord(levels_name).points)

                plot_data = np.zeros((ny, nx, n_levs))

                time = np.around(cube.coord('time').points, decimals=1)

                # Compute the horizontal grid
                x = np.around(cube.coord('longitude').points, decimals=5)
                y = np.around(cube.coord('latitude').points, decimals=5)

                xmin = np.amin(x)
                xmax = np.amax(x)
                ymin = np.amin(y)
                ymax = np.amax(y)

                # Generate a regular grid to interpolate the data.
                xi = np.linspace(xmin, xmax, nx)
                yi = np.linspace(ymin, ymax, ny)

                xf, yf = np.meshgrid(xi, yi)

                # Choose the correct vertical level set
                if n_full == n_levs:
                    zi = zi_f
                else:
                    zi = zi_h

                # Interpolate using delaunay triangularization
                for p, level in enumerate(range(n_levs)):
                    data = cube.data[t, level]
                    fi = griddata((x, y), data, (xf, yf), method='linear')
                    fi_n = griddata((x, y), data, (xf, yf), method='nearest')
                    fi[np.isnan(fi)] = fi_n[np.isnan(fi)]

                    plot_data[:, :, level] = fi

                    if cfield == 'exner':
                        # Convert to hPa
                        rd = 287.05
                        p0 = 100000.0
                        kappa = rd/1005.0
                        plot_data[:, :, level] = 0.01*fi**(1.0/kappa) * p0

                if (cfield == 'm_v' or  cfield == 'm_cl'):
                    # Plot level 10 for moisture fields
                    level = 10
                    cmap = magma.reversed()
                else:
                    level = 0
                    cmap = magma

                if direction == 'xz':
                    lon, height = np.meshgrid(xi, zi)
                    CS = ax.contourf(lon, height,
                                        plot_data[:, plot_lat, :].T,
                                        levels=levels, cmap=cmap)
                    plt.colorbar(cmap=cmap)
                    CL = ax.contour(lat, height,
                                        plot_data[:, plot_lat, :].T,
                                        levels=levels, linewidths=0.5,
                                        colors='k')
                    ax.title(['lat = ', yi[plot_lat]*360./np.real(nx)])
                if direction == 'yz':
                    lat, height = np.meshgrid(yi, zi)
                    CS = ax.contourf(lat, height,
                                        plot_data[:, plot_long, :].T,
                                        levels=levels, cmap=cmap)
                    ax.colorbar(cmap=cmap)
                    CL = ax.contour(lat, height,
                                        plot_data[:, plot_long, :].T,
                                        levels=levels, linewidths=0.5,
                                        colors='k')
                    ax.title(['long = ', xi[plot_long]*360./np.real(nx)])
                if direction == 'xy':
                    lat, lon = np.meshgrid(yi, xi)
                    if cfield == 'exner':
                        # Extrapolate data to the surface
                        dz = plot_data[:, :, 0] + (zi_f[0] - zi_h[0]) * \
                            (plot_data[:, :, 0] - plot_data[:, :, level]) \
                            / (zi_h[0] - zi_h[1])
                    else:
                        dz = plot_data[:, :, level]
                    if cfield != 'exner':
                        CS = ax.contourf(lon, lat,
                                            plot_data[:, :, level].T,
                                            levels=levels, cmap=cmap)
                        _ = interp_fig.colorbar(CS, ax=ax)
                    if cfield != 'theta':
                        CL = ax.contour(lon, lat, dz.T, levels=levels,
                                            linewidths=1.0, colors='k')
                        ax.clabel(CL, CL.levels[1::2], fontsize=15,
                                    inline=1, fmt='%3.1f')

            pngfile = '%s/baroclinic_plot-%s-time%s-%s.png' % \
                (plotpath, cfield, time[t], direction)
            plt.savefig(pngfile)
            plt.close()


if __name__ == "__main__":

    try:
        args = sys.argv[:]
        filein, plotpath, vertical_grid, formulation = args[1:5]
        field_list = None
        if len(args[:]) > 5:
            field_list = args[5].split(':')
    except ValueError:
        print("Usage: {0} <filein> <plotpath> <vertical_grid> <formulation> [<fields_list>] "
              .format(sys.argv[0]))
        exit(1)

    make_figures(filein, plotpath, field_list, vertical_grid, formulation)
