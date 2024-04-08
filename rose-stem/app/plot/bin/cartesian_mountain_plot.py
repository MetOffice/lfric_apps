#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
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
levels0 = None
data = None
data0 = None


def make_figure(plotpath, field, component, timestep, plotlevel_x,
                plotlevel_y, case, zoom, cntrs):

    val_col = 'c' + str(component)

    # get min and max of x,y data for plot axes

    min_lev = min(levels)

    xmin = data.loc[data['level'] == min_lev]['x'].min()
    xmax = data.loc[data['level'] == min_lev]['x'].max()
    ymin = data.loc[data['level'] == min_lev]['y'].min()
    ymax = data.loc[data['level'] == min_lev]['y'].max()

    # zmin is min of bottom level

    zmin = data['z'].min()/1000.0

    # zmax is a test-dependent value
    if zoom == 'zoom_1':  # zooming as in paper
        if (case != 'cosine'):
            zmax = 12.0
        else:
            print('Error: cosine orography does not have the zoom option')
    elif zoom == 'zoom_2':  # zooming up to bottom of damping layer
        if case == 'nhmw':
            zmax = 25.0
        elif case == 'hmw':
            zmax = 30.0
        elif (case == 'schar') or (case == 'schar3d'):
            zmax = 20.0
        elif case == 'bell':
            zmax = 10.0
        else:
            print('Error: cosine orography does not have the zoom option')
    elif zoom == 'no_zoom':  # plotting the entire domain depth
        if case == 'nhmw':
            zmax = 35.0
        elif case == 'hmw':
            zmax = 50.0
        elif (case == 'schar') or (case == 'schar3d'):
            zmax = 30.0
        elif case == 'cosine':
            zmax = 6.0
        elif case == 'bell':
            zmax = 16.0

    val_col = 'c' + str(component)
    if(case == 'schar3d'):
        # Hard-coded for vertical-slice-like runs
        if(field == 'u' and (component == '1' or component == '2')):
            nx = len(data.loc[data['level'] == min_lev])//400
            ny = 400
        else:
            nx = len(data.loc[data['level'] == min_lev])//200
            ny = 200
    elif(case == 'bell'):
        # Hard-coded for vertical-slice-like runs
        if(field == 'u' and (component == '1' or component == '2')):
            nx = len(data.loc[data['level'] == min_lev])//600
            ny = 600
        else:
            nx = len(data.loc[data['level'] == min_lev])//200
            ny = 200
    else:
        # Hard-coded for vertical-slice-like runs
        if(field == 'u' and (component == '1' or component == '2')):
            nx = len(data.loc[data['level'] == min_lev])//8
            ny = 8
        else:
            nx = len(data.loc[data['level'] == min_lev])//4
            ny = 4

    nl = len(levels)

    # Create 2D plot
    val_i = np.zeros([ny, nx, nl])
    val_i0 = np.zeros([ny, nx, nl])
    x_i = np.zeros([ny, nx, nl])
    y_i = np.zeros([ny, nx, nl])
    height_i = np.zeros([ny, nx, nl])

    # Interpolate field values and heights onto xy for each level
    for p in range(len(levels)):
        p_data = data.loc[data['level'] == levels[p]]

        # Using reshape of numpy array
        val_i[:, :, p] = (p_data[val_col].values).reshape((ny, nx))
        height_i[:, :, p] = (p_data['z'].values).reshape((ny, nx))/1000.0
        x_i[:, :, p] = (p_data['x'].values).reshape((ny, nx))/1000.0
        y_i[:, :, p] = (p_data['y'].values).reshape((ny, nx))/1000.0

    if field in ('theta', 'm_cl', 'rho'):
        for p in range(len(levels0)):
            p_data = data0.loc[data0['level'] == levels0[p]]
            # Using reshape of numpy array
            val_i0[:, :, p] = (p_data[val_col].values).reshape((ny, nx))

    if field in ('theta', 'm_cl', 'rho'):
        # subtracting entire background profile (profile at initial time)
        val_i -= val_i0

    # Setting contour limits and intervals for vertical velocity

    if case == 'nhmw':
        cmin_w = -0.0048
        cmax_w = 0.0048
        nc_w = 17
        cmin_t = -0.002
        cmax_t = 0.002
        nc_t = 9
    elif case == 'hmw':
        cmin_w = -0.004
        cmax_w = 0.004
        nc_w = 17
    elif (case == 'schar') or (case == 'schar3d'):
        cmin_w = -0.5
        cmax_w = 0.5
        nc_w = 21
    elif case == 'cosine':
        cmin_w = -3.0
        cmax_w = 3.0
        nc_w = 13
    elif case == 'bell':
        cmin_w = -2.5
        cmax_w = 2.5
        nc_w = 21

    c_map = viridis

    # xz Plot

    # Take actual orography heights as z
    zi_adj = height_i[int(plotlevel_y), :, :]

    # Extract the slice for plotting
    dval = np.zeros([nx, len(levels)])

    if case == 'bell':
        dval = 0.5*(val_i[int(plotlevel_y), :, :]
                    + val_i[int(plotlevel_y)+1, :, :])
    else:
        dval = val_i[int(plotlevel_y), :, :]

    x_plt = x_i[int(plotlevel_y), :, :]

    fig_height = 7.5

    fig = plt.figure(figsize=(15, fig_height))

    # Plots
    if (field == 'w3projection_u1' or (field == 'u' and component == '1')):
        if case == 'nhmw':
            if cntrs == 'lines':
                cc = np.linspace(9.99, 10.01, 21)
                cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                                 cc, cmap=cm.Spectral, linewidths=3)
            elif cntrs == 'colours':
                cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10),
                                    vmin=9.99, vmax=10.01, cmap=cm.coolwarm)
        elif case == 'cosine':
            if cntrs == 'lines':
                cc = np.linspace(7.0, 12.6, 15)
                cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                                 cc, colors='k', linewidths=3)
                plt.clabel(cf, fontsize=10)
            elif cntrs == 'colours':
                cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10),
                                    vmin=7.0, vmax=12.6, cmap=cm.coolwarm)
        else:
            if cntrs == 'lines':
                cc = np.linspace(np.min(dval), np.max(dval), 13)
                cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                                 cc, cmap=cm.Spectral, linewidths=3)
            else:
                cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10))
                plt.colorbar(cf,  cmap=cm.coolwarm)
    elif (field == 'w3projection_u3' or (field == 'u' and component == '3')):
        if cntrs == 'lines':
            if case == 'bell':
                cc = np.concatenate((np.linspace(-2.5, -0.25, 10),
                                     np.linspace(0.25, 2.5, 10)), axis=0)
            else:
                cc = np.linspace(cmin_w, cmax_w, nc_w)
            if case == 'cosine':
                cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                                 cc, colors='k', linewidths=3)
                plt.clabel(cf, fontsize=10)
            else:
                cf = plt.contourf(x_plt, zi_adj, np.round(dval, 10),
                                  cc, cmap=c_map, linewidths=3, extend='both')
                cl = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                                 cc, linewidths=2, colors='k')
                if case == 'schar':
                    cb = plt.colorbar(cf)
                elif case == 'hmw':
                    cb = plt.colorbar(cf)
                else:
                    cb = plt.colorbar(cf)
                cb.ax.tick_params(labelsize=24)
        elif cntrs == 'colours':
            if timestep == 'T000000':
                cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10),
                                    cmap=cm.coolwarm)
            else:
                cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10),
                                    vmin=cmin_w, vmax=cmax_w,
                                    cmap=cm.coolwarm)
    elif (field == 'w3projection_u2' or (field == 'u' and component == '2')):
        if cntrs == 'lines':
            cc = np.linspace(np.min(dval), np.max(dval), 13)
            if case == 'cosine':
                cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                                 cc, colors='k', linewidths=3)
                plt.clabel(cf, fontsize=10)
            else:
                cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                                 cc, cmap=cm.Spectral, linewidths=3)
        elif cntrs == 'colours':
            cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10),
                                vmin=np.min(dval), vmax=np.max(dval),
                                cmap=cm.coolwarm)
    elif field == 'theta':
        if cntrs == 'lines':
            cc = np.linspace(cmin_t, cmax_t, nc_t)
            cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                             cc, cmap=cm.Spectral,
                             linewidths=3)
            plt.colorbar(cf,  cmap=cm.Spectral)
        elif cntrs == 'colours':
            cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10),
                                vmin=cmin_t, vmax=cmax_t,
                                cmap=cm.coolwarm)
            plt.colorbar(cf,  cmap=cm.coolwarm)
    elif field == 'rho':
        if cntrs == 'lines':
            cc = np.linspace(-0.00001, 0.00001, 21)
            cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                             cc, cmap=cm.Spectral, linewidths=3)
            plt.colorbar(cf,  cmap=cm.Spectral)
        elif cntrs == 'colours':
            cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10),
                                vmin=-0.00001, vmax=-0.00001,
                                cmap=cm.coolwarm)
            plt.colorbar(cf,  cmap=cm.coolwarm)
    else:
        if cntrs == 'lines':
            cc = np.linspace(-0.0005, 0.0005, 21)
            cf = plt.contour(x_plt, zi_adj, np.round(dval, 10),
                             cc, cmap=cm.Spectral, linewidths=3)
            plt.colorbar(cf,  cmap=cm.Spectral)
        elif cntrs == 'colours':
            cf = plt.pcolormesh(x_plt, zi_adj, np.round(dval, 10),
                                cmap=cm.coolwarm)
            plt.colorbar(cf,  cmap=cm.coolwarm)

    plt.title('max: %2.4e, min: %2.4e' % (np.max(dval), np.min(dval)))
    plt.xlabel("x(km)", fontsize=32)
    plt.ylabel("z(km)", fontsize=32)

    if zoom == 'no_zoom':  # plotting the entire domain horizontal length
        if case == 'nhmw':
            xmin_lim = -72
            xmax_lim = 72
        elif case == 'hmw':
            xmin_lim = -120
            xmax_lim = 120
        elif case == 'schar':
            xmin_lim = -50
            xmax_lim = 50
        elif case == 'bell':
            xmin_lim = -30
            xmax_lim = 30
    else:  # plotting only a portion of horizontal domain length
        if case == 'nhmw':
            xmin_lim = -12.5
            xmax_lim = 32.5
            plt.xticks(np.arange(-10, 35, 5))
        elif case == 'hmw':
            xmin_lim = -40
            xmax_lim = 40
        elif case == 'bell':
            xmin_lim = -8
            xmax_lim = 24
            plt.xticks(np.arange(-8, 28, 4))
        elif (case == 'schar') or (case == 'schar3d'):
            xmin_lim = -20
            xmax_lim = 20
            plt.xticks(np.arange(-20, 25, 5))
        else:
            print('Error: cosine orography does not have the zoom option.')

    if zoom == 'zoom_1':  # Ticks as in paper
        plt.yticks(np.arange(0, 14, 2))
    elif zoom == 'zoom_2':
        if case == 'bell':
            plt.yticks(np.arange(0, 20, 5))
    else:
        if case == 'cosine':
            plt.yticks(np.arange(0, 7, 1))

    plt.xlim(xmin_lim, xmax_lim)
    plt.tick_params(axis='both', labelsize=32)
    plt.ylim(zmin, zmax)

    if field in ['u', 'xi']:
        out_file_name = plotpath + "/" "nodal_slices_xz_" \
            + field + str(component) + '_' + ts + ".png"
    else:
        out_file_name = plotpath + "/" "nodal_slices_xz_" \
            + field + '_' + ts + ".png"
    plt.tight_layout()
    plt.savefig(out_file_name)

    # Plot yz  and xy slices if 3d
    if(case == 'schar3d' or case == 'bell'):

        # yz Plot

        # Take actual orography heights as z
        zi_adj = height_i[:, int(plotlevel_x), :]

        # Extract the slice for plotting
        dval = np.zeros([ny, len(levels)])
        if case == 'bell':
            dval = 0.5*(val_i[:, int(plotlevel_x), :]
                        + val_i[:, int(plotlevel_x)+1, :])
        else:
            dval = val_i[:, int(plotlevel_x), :]

        y_plt = y_i[:, int(plotlevel_x), :]

        fig = plt.figure(figsize=(15, fig_height))

        # Plots
        if(field == 'w3projection_u1' or (field == 'u' and component == '1')):
            cf = plt.pcolormesh(y_plt, zi_adj, np.round(dval, 10))
            plt.colorbar(cf,  cmap=cm.coolwarm)
        elif(field == 'w3projection_u3' or
             (field == 'u' and component == '3')):
            if cntrs == 'lines':
                if case == 'schar3d':
                    cc = np.linspace(cmin_w, cmax_w, nc_w)
                elif case == 'bell':
                    cc = np.concatenate((np.linspace(-0.6, -0.1, 6),
                                         np.linspace(0.1, 0.6, 6)), axis=0)

                cf = plt.contourf(y_plt, zi_adj, np.round(dval, 10),
                                  cc, cmap=c_map, linewidths=3, extend='both')
                cl = plt.contour(y_plt, zi_adj, np.round(dval, 10),
                                 cc, linewidths=2, colors='k')
                cb = plt.colorbar(cf)
                cb.ax.tick_params(labelsize=24)
            elif cntrs == 'colours':
                cf = plt.pcolormesh(y_plt, zi_adj, np.round(dval, 10),
                                    vmin=cmin_w, vmax=cmax_w,
                                    cmap=cm.coolwarm)
        elif (field == 'w3projection_u2' or
              (field == 'u' and component == '2')):
            if cntrs == 'lines':
                cc = np.linspace(np.min(dval), np.max(dval), 13)
                cf = plt.contour(y_plt, zi_adj, np.round(dval, 10), cc,
                                 cmap=cm.Spectral, linewidths=3)
            elif cntrs == 'colours':
                cf = plt.pcolormesh(y_plt, zi_adj, np.round(dval, 10),
                                    vmin=np.min(dval), vmax=np.max(dval),
                                    cmap=cm.coolwarm)
        elif field == 'theta':
            if cntrs == 'lines':
                cc = np.linspace(cmin_t, cmax_t, nc_t)
                cf = plt.contour(y_plt, zi_adj, np.round(dval, 10), cc,
                                 cmap=cm.Spectral, linewidths=3)
                plt.colorbar(cf,  cmap=cm.Spectral)
            elif cntrs == 'colours':
                cf = plt.pcolormesh(y_plt, zi_adj, np.round(dval, 10),
                                    vmin=cmin_t, vmax=cmax_t,
                                    cmap=cm.coolwarm)
                plt.colorbar(cf,  cmap=cm.coolwarm)
        elif field == 'rho':
            if cntrs == 'lines':
                cc = np.linspace(-0.00001, 0.00001, 21)
                cf = plt.contour(y_plt, zi_adj, np.round(dval, 10), cc,
                                 cmap=cm.Spectral, linewidths=3)
                plt.colorbar(cf,  cmap=cm.Spectral)
            elif cntrs == 'colours':
                cf = plt.pcolormesh(y_plt, zi_adj, np.round(dval, 10),
                                    vmin=-0.00001, vmax=-0.00001,
                                    cmap=cm.coolwarm)
                plt.colorbar(cf,  cmap=cm.coolwarm)
        else:
            if cntrs == 'lines':
                cc = np.linspace(-0.0005, 0.0005, 21)
                cf = plt.contour(y_plt, zi_adj, np.round(dval, 10), cc,
                                 cmap=cm.Spectral, linewidths=3)
                plt.colorbar(cf,  cmap=cm.Spectral)
            elif cntrs == 'colours':
                cf = plt.pcolormesh(y_plt, zi_adj, np.round(dval, 10),
                                    cmap=cm.coolwarm)
                plt.colorbar(cf,  cmap=cm.coolwarm)

        plt.title('max: %2.4e, min: %2.4e' % (np.max(dval), np.min(dval)))
        plt.xlabel("y(km)", fontsize=32)
        plt.ylabel("z(km)", fontsize=32)

        if zoom == 'no_zoom':  # plotting the entire domain horizontal length
            if case == 'schar3d':
                xmin_lim = -50
                xmax_lim = 50
            elif case == 'bell':
                xmin_lim = -20
                xmax_lim = 20
            else:
                xmin_lim = -2
                xmax_lim = 2
        else:  # plotting only a portion of horizontal domain length
            if case == 'bell':
                xmin_lim = -10
                xmax_lim = 10
            elif case == 'schar3d':
                xmin_lim = -20
                xmax_lim = 20
                plt.xticks(np.arange(-20, 25, 5))
            else:
                print('Error: please select one of schar3d, bell orography.')

        if zoom == 'zoom_1':  # Ticks as in paper
            plt.yticks(np.arange(0, 14, 2))
        elif zoom == 'zoom_2':  # Ticks as in paper
            if case == 'bell':
                plt.yticks(np.arange(0, 15, 5))

        plt.xlim(xmin_lim, xmax_lim)
        plt.tick_params(axis='both', labelsize=32)
        plt.ylim(zmin, zmax)

        if field in ['u', 'xi']:
            out_file_name = plotpath + "/" "nodal_slices_yz_" + field + \
                str(component) + '_' + ts + ".png"
        else:
            out_file_name = plotpath + "/" "nodal_slices_yz_" + field + \
                '_' + ts + ".png"
        plt.tight_layout()
        plt.savefig(out_file_name)

        # xy plot
        # extract the slice for plotting, hard-coded level 10/100
        dval = np.zeros([ny, len(levels)])

        dval = val_i[:, :, 10]

        x_plt = x_i[:, :, 10]
        y_plt = y_i[:, :, 10]

        fig = plt.figure(figsize=(15, fig_height))

        if case == 'schar3d':
            cc = np.linspace(np.min(dval), np.max(dval), 13)
            cf = plt.contourf(x_plt, y_plt, np.round(dval, 10), cc,
                             cmap=c_map, linewidths=3)
            cl = plt.contour(x_plt, y_plt, np.round(dval, 10),
                             cc, linewidths=2, colors='k')


            plt.xlim(-20, 20)
            plt.ylim(-20, 20)
            plt.tick_params(axis='both', labelsize=24)
            plt.xticks(np.arange(-20, 25, 5))
            plt.yticks(np.arange(-20, 25, 5))

        elif case == 'bell':
            cc = np.concatenate((np.linspace(-1.0, -0.1, 10),
                                 np.linspace(0.1, 1.0, 10)), axis=0)
            cf = plt.contourf(x_plt, y_plt, np.round(dval, 10), cc,
                              cmap=c_map, linewidths=3, extend='both')
            cl = plt.contour(x_plt, y_plt, np.round(dval, 10),
                             cc, linewidths=2, colors='k')
            cb = plt.colorbar(cf)
            cb.ax.tick_params(labelsize=24)

            plt.xlim(-8, 24)
            plt.ylim(-8, 8)
            plt.xticks(np.arange(-8, 28, 4))
            plt.yticks(np.arange(-8, 16, 8))

        plt.title('max: %2.4e, min: %2.4e' % (np.max(dval), np.min(dval)))
        plt.xlabel("x(km)", fontsize=32)
        plt.ylabel("y(km)", fontsize=32)
        plt.tick_params(axis='both', labelsize=32)

        if field in ['u', 'xi']:
            out_file_name = plotpath + "/" "nodal_slices_xy_" + field + \
                str(component) + '_' + ts + ".png"
        else:
            out_file_name = plotpath + "/" "nodal_slices_xy_" + field + \
                '_' + ts + ".png"
        plt.tight_layout()
        plt.savefig(out_file_name)

        if case == 'bell':
            # Plot also the xy slice at level 4 for the bell case
            dval = np.zeros([ny, len(levels)])
            dval = val_i[:, :, 4]

            x_plt = x_i[:, :, 4]
            y_plt = y_i[:, :, 4]

            fig = plt.figure(figsize=(15, fig_height))

            cc = np.concatenate((np.linspace(-1.5, -0.1, 15),
                                 np.linspace(0.1, 1.5, 15)), axis=0)
            cf = plt.contourf(x_plt, y_plt, np.round(dval, 10), cc,
                              cmap=c_map, linewidths=3, extend='both')
            cl = plt.contour(x_plt, y_plt, np.round(dval, 10),
                             cc, linewidths=2, colors='k')
            cb = plt.colorbar(cf)
            cb.ax.tick_params(labelsize=24)

            plt.xlim(-8, 24)
            plt.ylim(-8, 8)
            plt.tick_params(axis='both', labelsize=32)
            plt.xticks(np.arange(-8, 28, 4))
            plt.yticks(np.arange(-8, 16, 8))

            plt.xlabel("x(km)", fontsize=32)
            plt.ylabel("y(km)", fontsize=32)
            plt.tick_params(axis='both', labelsize=32)

            out_file_name = plotpath + "/" "nodal_slices_xy_" + field + \
                component + '_' + ts + "_zlev4.png"
            plt.tight_layout()
            plt.savefig(out_file_name)

if __name__ == "__main__":

    try:
         config, datapath, fields, component, timesteps, plotlevel_x, plotlevel_y, \
            case, zoom, cntrs, plotpath = sys.argv[1:12]
    except ValueError:
        print("Usage: {0} <file_stem_name> <datapath> <field_names> <component> "
              "<timestep_list> <plotlevel_x> <plotlevel_y> <case> "
              "<zoom> <cntrs> <plotpath>".format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    # Split out the list of timesteps
    ts_list = timesteps.split(':')

    for field in field_list:

        if((field != 'u') & (field != 'xi') & (field != 'w3projection_u1') &
           (field != 'w3projection_u2') & (field != 'w3projection_u3')):
            if (component != '1'):
                print("Scalars have only one component!")
                exit(1)
        if field in ('theta', 'm_cl', 'rho'):
            # Create initial data for theta
            filestem = datapath + "/" + config + "_nodal_" + \
                field + "_T000000" + "*"
            data0 = read_nodal_data(filestem, 1, component)

            # Sort the data (needed to be able to reshape and not regrid)
            data0 = data0.sort(['y', 'x', 'z'])
            levels0 = np.sort(data0.level.unique())

    for ts in ts_list:

        # Making space for more fields in the plot.

        for i, field in enumerate(field_list):

            filestem = datapath + "/" + config + "_nodal_" + \
                field + "_" + ts + "*"

            if field in ['u', 'w3projection_u1',
                         'w3projection_u2', 'w3projection_u3']:
                data = read_nodal_data(filestem, 3, component)
            else:
                data = read_nodal_data(filestem, 1, component)

            if data is not None:

                levels = np.sort(data.level.unique())

                # Sort the data (needed to be able to reshape and not regrid)
                data = data.sort_values(['y', 'x', 'z'])

                # Only try to plot if we found some files for this timestep
                if len(levels) > 0:
                    make_figure(plotpath, field, component, ts,
                                plotlevel_x, plotlevel_y, case, zoom, cntrs)
