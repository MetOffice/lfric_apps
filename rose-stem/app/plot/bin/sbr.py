#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2019 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Need to set a non-interactive backend for suites
'''
Plot horizontal and vertical cross sections of desired fields with
fixed contour intervals
'''
import matplotlib
matplotlib.use('Agg')

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
    plt.register_cmap(name='magma', cmap=magma)

if iris.__version__ < "3.0.0":
    iris.FUTURE.netcdf_promote = True

# Size of regular grid
ny, nx = 360, 720


def make_figures(filein, plotpath, fields, vertical_spacing, lid, n_full,
                 figname, idx_list):

    if vertical_spacing == 'um38':
        # um L38 set
        zi_f = np.array([.0, .0005095, .0020380, .0045854, .0081519, .0127373,
                         .0183417, .0249651, .0326074, .0412688, .0509491,
                         .0616485, .0733668, .0861040, .0998603, .1146356,
                         .1304298, .1472430, .1650752, .1839264, .2037966,
                         .2246857, .2465938, .2695209, .2934670, .3184321,
                         .3444162, .3714396, .3998142, .4298913, .4620737,
                         .4968308, .5347160, .5763897, .6230643, .6772068,
                         .7443435, .8383348, 1.000000])*lid

    elif vertical_spacing == 'um70':
        # um L70 set
        zi_f = np.array([.0000000, .0002500, .0006667, .0012500, .0020000,
                         .0029167, .0040000, .0052500, .0066667, .0082500,
                         .0100000, .0119167, .0140000, .0162500, .0186667,
                         .0212500, .0240000, .0269167, .0300000, .0332500,
                         .0366667, .0402500, .0440000, .0479167, .0520000,
                         .0562500, .0606667, .0652500, .0700000, .0749167,
                         .0800000, .0852500, .0906668, .0962505, .1020017,
                         .1079213, .1140113, .1202745, .1267154, .1333406,
                         .1401592, .1471838, .1544313, .1619238, .1696895,
                         .1777643, .1861929, .1950307, .2043451, .2142178,
                         .2247466, .2360480, .2482597, .2615432, .2760868,
                         .2921094, .3098631, .3296378, .3517651, .3766222,
                         .4046373, .4362943, .4721379, .5127798, .5589045,
                         .6112759, .6707432, .7382500, .8148403, .9016668,
                         1.0000000])*lid

    elif vertical_spacing == 'um85':
        # um L85 set
        zi_f = np.array([
            0.0000000E+00, 0.2352941E-03, 0.6274510E-03, 0.1176471E-02, 0.1882353E-02, 
            0.2745098E-02, 0.3764706E-02, 0.4941176E-02, 0.6274510E-02, 0.7764705E-02, 
            0.9411764E-02, 0.1121569E-01, 0.1317647E-01, 0.1529412E-01, 0.1756863E-01, 
            0.2000000E-01, 0.2258823E-01, 0.2533333E-01, 0.2823529E-01, 0.3129411E-01, 
            0.3450980E-01, 0.3788235E-01, 0.4141176E-01, 0.4509804E-01, 0.4894118E-01, 
            0.5294117E-01, 0.5709804E-01, 0.6141176E-01, 0.6588235E-01, 0.7050980E-01, 
            0.7529411E-01, 0.8023529E-01, 0.8533333E-01, 0.9058823E-01, 0.9600001E-01, 
            0.1015687E+00, 0.1072942E+00, 0.1131767E+00, 0.1192161E+00, 0.1254127E+00, 
            0.1317666E+00, 0.1382781E+00, 0.1449476E+00, 0.1517757E+00, 0.1587633E+00, 
            0.1659115E+00, 0.1732221E+00, 0.1806969E+00, 0.1883390E+00, 0.1961518E+00, 
            0.2041400E+00, 0.2123093E+00, 0.2206671E+00, 0.2292222E+00, 0.2379856E+00, 
            0.2469709E+00, 0.2561942E+00, 0.2656752E+00, 0.2754372E+00, 0.2855080E+00, 
            0.2959203E+00, 0.3067128E+00, 0.3179307E+00, 0.3296266E+00, 0.3418615E+00, 
            0.3547061E+00, 0.3682416E+00, 0.3825613E+00, 0.3977717E+00, 0.4139944E+00, 
            0.4313675E+00, 0.4500474E+00, 0.4702109E+00, 0.4920571E+00, 0.5158098E+00, 
            0.5417201E+00, 0.5700686E+00, 0.6011688E+00, 0.6353697E+00, 0.6730590E+00, 
            0.7146671E+00, 0.7606701E+00, 0.8115944E+00, 0.8680208E+00, 0.9305884E+00, 
            0.1000000E+01])*lid


    elif vertical_spacing == 'dcmip':
        # dcmip Stretched grid
        mu = 15.
        zi_f = np.zeros([n_full])
        for k in range(n_full):
            eta = float(k)/float(n_full)
            zi_f[k] = ((np.sqrt(mu*eta**2 + 1.) - 1.)
                      /(np.sqrt(mu + 1.) - 1) * lid)
    else:
        # assume uniform grid
        zmin = 0.0
        zmax = lid
        zi_f = np.linspace(zmin, zmax, n_full)

    zi_h = 0.5*(zi_f[1:] + zi_f[0:n_full-1])

    directions = ['xy', 'yz', 'xz']

    for t in [-1]:
        if fields is None:
            fields_name = 'winds'
            fields = ['u_in_w2h', 'v_in_w2h', 'w_in_wth']
            separate_plots = False
        elif len(fields) == 1:
            fields_name = fields[0]
            separate_plots = False
        else:
            fields_name = 'fields'
            separate_plots = True

        nxplots = len(fields)
        nyplots = len(directions)

        if separate_plots:
            nxplots = 1
        else:
            interp_fig = plt.figure(figsize=(20, 10))
        nplots = nxplots*nyplots
        for f, field in enumerate(fields):
            if separate_plots:
                interp_fig = plt.figure(figsize=(20, 10))
                fields_name = field
            cube = read_ugrid_data(filein, field)
            # Vertical levels will be last entry in dimension coords
            levels_name = cube.dim_coords[-1].name()

            # Set some levels for contours:
            levels = None
            if field == 'theta':
                levels = np.linspace(220, 330, 12)
            if field == 'u_in_w2h':
                levels = np.linspace(-5, 45, 11)
            if field == 'v_in_w2h':
                levels = np.linspace(-1.5, 1.5, 11)
            if field == 'w_in_wth':
                levels = np.linspace(-0.25, 0.25, 11)
            if field == 'exner':
                # exner will be converted to hPa
                levels = np.linspace(916, 1020, 14)
            if field == 'density':
                levels = np.arange(0, 1.4, 0.05)

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

            if idx_list is None:
                plot_long = int(nx/2)
                plot_lat = int(ny/2)
                plot_level = n_levs-2
            else:
                # idx_list contains degrees/level to plot
                plot_long = (int(idx_list[0])+180)*int(nx/360)
                plot_lat = (int(idx_list[1])+90)*int(ny/180)
                plot_level = int(idx_list[2])


            # Interpolate using delaunay triangularization
            for p, l in enumerate(range(n_levs)):
                data = cube.data[t, l]
                fi = griddata((x, y), data, (xf, yf), method='linear')
                fi_n = griddata((x, y), data, (xf, yf), method='nearest')
                fi[np.isnan(fi)] = fi_n[np.isnan(fi)]

                plot_data[:, :, l] = fi

                if field == 'exner_pressure':
                    # Convert to hPa
                    rd = 287.05
                    p0 = 100000.0
                    kappa = rd/1005.0
                    plot_data[:, :, l] = 0.01*fi**(1.0/kappa)*p0

            for d, direction in enumerate(directions):

                if separate_plots:
                    plotnum = d + 1
                else:
                    plotnum = d*nxplots + f + 1

                ax = interp_fig.add_subplot(nyplots, nxplots, plotnum)

                if direction == 'xz':
                    x1, x2 = np.meshgrid(xi, zi)
                    x3 = plot_data[plot_lat, :, :].T
                    plt.title([field, direction, ' lat = ',
                               yi[plot_lat]*360./np.real(nx)])
                if direction == 'yz':
                    x1, x2 = np.meshgrid(yi, zi)
                    x3 = plot_data[:, plot_long, :].T
                    plt.title([field, direction, 'long = ',
                               xi[plot_long]*360./np.real(nx)])
                if direction == 'xy':
                    x2, x1 = np.meshgrid(yi, xi)
                    x3 = plot_data[:, :, plot_level].T
                    plt.title([field, direction, 'Height = ', zi[plot_level]])

                CS = plt.contourf(x1, x2, x3, levels=levels, cmap=magma)
                plt.colorbar(cmap=magma)
                CL = plt.contour(x1, x2, x3,
                                 levels=levels, linewidths=0.5, colors='k')

            pngfile = '%s/%s-%s-time%s.png' % (plotpath, figname, fields_name, time[t])

            if separate_plots:
                plt.tight_layout()        
                plt.savefig(pngfile)

        if ~separate_plots:
            plt.tight_layout()       
            plt.savefig(pngfile)

        plt.close()

if __name__ == "__main__":

    try:
        args = sys.argv[:]
        files, plotpath, vertical_grid, lid, n_full, figname = args[1:7]
        field_list = None
        if len(args[:]) > 7:
            field_list = args[7].split(':')
        idx_list = None
        if len(args[:]) > 8:
            idx_list = args[8].split(':')
    except ValueError:
        print("Usage: {0} <filein> <plotpath> <vertical_grid> <lid>"
              " <n_full> <figname> [<fields_list>]"
              .format(sys.argv[0]))
        exit(1)

    file_list = files.split(':')
    for filein in file_list:
        make_figures(filein, plotpath, field_list, vertical_grid,
                     int(lid), int(n_full), figname, idx_list)
