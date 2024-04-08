#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2017 Met Office. All rights reserved.
# For further details please refer to the file LICENCE which you should have
# received as part of this distribution.
##############################################################################
'''
Python script to plot SWE slices.

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

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import math

from scipy.interpolate import griddata

import sys
# Use magma colormap
try:
    magma = plt.get_cmap("magma")
except ValueError:
    from matplotlib.colors import ListedColormap
    from magma import magma_data
    magma = ListedColormap(magma_data, name='magma')
    plt.register_cmap(name='magma', cmap=magma)

from read_data import read_nodal_data

def make_figures(datapath, plotpath, nx, field_list, ts_list):
    """Function to plot shallow water miniapp field output.
    :arg datapath: path to data files, including configuration part of file names
    :arg plotpath: path to folder where plots should be saved
    :arg nx: either 'C' for spherical mesh, or integer for planar square resolution
    :arg field_list: list of fields to be plotted, as given in dataset
    :arg ts_list: list of timestep strings, e.g. 'T000000' for initial time"""

    # Assuming equal resolution in x and y directions for planar domain
    if nx[0] == 'C':
        nnxx = 'C'
        nx, ny, r2d = 720, 360, 180.0/np.pi
        xax, yax = 'Longitude', 'Latitude'
        xmin, xmax = -np.pi, np.pi
        ymin, ymax = -np.pi/2., np.pi/2.
    elif nx[0] == 'P':
        nnxx = 'P'
        _, res, length = nx.split(':')
        nx, lx = int(res), int(length)
        ny, r2d = nx, 1
        xax, yax = 'x', 'y'
        xmin, xmax = -lx/2., lx/2.
        ymin, ymax = -lx/2., lx/2.

    # Create meshgrid for data to be interpolated onto
    xi, yi = np.meshgrid(np.linspace(xmin, xmax, nx),
                         np.linspace(ymin, ymax, ny))
    zi = np.empty([len(ts_list), ny, nx])
    zi[:] = np.nan
    # ze is normalised error
    ze = np.empty([len(ts_list), ny, nx])
    ze[:] = np.nan
    # za is absolute error
    za = np.empty([len(ts_list), ny, nx])
    za[:] = np.nan
    maxerror = np.empty([len(ts_list)])
    maxerror[:] = np.nan
    l2error = np.empty([len(ts_list)])
    l2error[:] = np.nan

    # Use flag to figure out if buoyancy is in a vector or scalar space
    bspace_is_v0, find_bspace = False, False
    # Same for vorticity space
    qspace_is_v0, find_qspace = True, True

    # Loop through fields and create frames for each field
    for field in field_list:
        if field[-1] in ['1', '2', '3']:
            fieldname = field[:-1]
        elif field == 'depth':
            fieldname = 'geopot'
        else:
            fieldname = field

        nodata = []
        for frame, ts in enumerate(ts_list):
            # Find data and interpolate it onto meshgrid for each frame
            filestem = datapath + fieldname + "_" + ts + "*"
            if field == 'geopot' or field == 'depth':
                data = read_nodal_data(filestem, 1, 1)
                val_col = 'c1'
            elif field == 'buoyancy':
                #if find_bspace:
                    # Figure out in first iteration which bspace was used, using data structure
                    #data = read_nodal_data(filestem, 3, 3)
                    #if data.empty:
                    #    data = read_nodal_data(filestem, 1, 1)
                    #    if data.empty:
                    #        # If buoyancy was a field to plot but has no data, skip it
                    #        print('No buoyancy data found. Did you run non-thermal SWE?')
                    #        nodata.extend(ts_list)
                    #        break
                    #    else:
                    #        bspace_is_v0 = False
                    #find_bspace= False

                #if bspace_is_v0:
                #    data = read_nodal_data(filestem, 3, 3)
                #    val_col = 'c3'
                #else:
                data = read_nodal_data(filestem, 1, 1)
                val_col = 'c1'
            elif field == 'q':
                #if find_qspace:
                #    # Figure out in first iteration which qspace was used, using data structure
                #    data = read_nodal_data(filestem, 3, 3)
                #    if 'c3' not in data.keys():
                #        raise IOError('Required column not in data (i.e. either wrong data or data missing)')
                #    if all(np.isnan(data['c3'])):
                #        data = read_nodal_data(filestem, 1, 1)
                #        if all(np.isnan(data['c1'])):
                #            # If q was a field to plot but has no data, skip it
                #            print('No vorticity data found.')
                #            nodata.extend(ts_list)
                #            break
                #        else:
                #            qspace_is_v0 = False
                #    find_qspace= False
                #
                #if qspace_is_v0:
                #    data = read_nodal_data(filestem, 3, 3)
                #    val_col = 'c3'
                #else:
                data = read_nodal_data(filestem, 1, 1)
                val_col = 'c1'
            else:
                data = read_nodal_data(filestem, 3, int(field[-1]))
                val_col = 'c' + field[-1]

            if (not data.empty):
                levels = data.level.unique() # levels = 0.5
                p_data = data.loc[data['level'] == levels[0]]
                zi[frame, :, :] = griddata((p_data['x'].values, p_data['y'].values),
                                            p_data[val_col].values, (xi, yi),
                                            method='linear')
                if field == 'depth':
                    gravity = 9.80616
                    zi[frame, :, :] /= gravity

                # jk
                if frame==0:
                    p_0 = p_data
            else:
                print('Data empty for field {0}, frame {1}. '\
                      ' Moving on to next frame.'.format(field, frame))
                nodata.append(frame)
                continue

            # Normalised Error

            ze[frame, :, :] = (zi[0, :, :] - zi[frame, :, :])/np.nanmax(abs(zi[0, :, :]))

            # Absolute Error

            za[frame, :, :] = (zi[0, :, :] - zi[frame, :, :])

            # calculate the L_inf error

            maxerror[frame] = np.nanmax(ze[frame, :, :])

            # calculate the L_2 error

            if nnxx == 'C':

                l2error[frame] = math.sqrt(np.nansum((abs(za[frame, :, :])**2)*np.cos(yi[:,:]))/np.nansum((abs(zi[0, :, :])**2)*np.cos(yi[:,:])))

            elif nnxx == 'P':

                l2error[frame] = math.sqrt(np.nansum((abs(za[frame, :, :])**2))/np.nansum((abs(zi[0, :, :])**2)))

        # Skip fields that did not come with any data
        if len(nodata) == len(ts_list):
            continue

        # Split ts loop in two to get global contour limits

        ze_min, ze_max = np.nanmin(ze), np.nanmax(ze)

        if abs(ze_max - ze_min) < 1e-13:
            cc_max, cc_min = ze_max + 1e-13, ze_min - 1e-13
            cc_nr = 8
        else:
            off_limit = (ze_max - ze_min)/10.
            cc_max = ze_max + off_limit
            cc_min = ze_min - off_limit
            cc_nr = 21
            # Rounding setup does not work in special cases; can turn it off
            max_shift, min_shift = cc_max*10**13, cc_min*10**13
            nrdgt_max = len(str(abs(int(max_shift))))
            nrdgt_min = len(str(abs(int(min_shift))))
            # Round to second digit of limit with larger magnitude
            if nrdgt_max > nrdgt_min:
                cc_max = (int(str(int(max_shift))[:2]) + 1)*10**(nrdgt_max - 2)
                cc_min = round(min_shift, 2 - nrdgt_max)
            elif nrdgt_max < nrdgt_min:
                cc_max = round(max_shift, 2 - nrdgt_min)
                cc_min = (int(str(int(min_shift))[:2]) + 1)*10**(nrdgt_min - 2)
            else:
                cc_max = round(max_shift, 2 - nrdgt_min)
                cc_min = round(min_shift, 2 - nrdgt_max)

            # Sometimes this rounding setup fails, in which case we
            # go back to the unrounded setup
            if cc_max <= cc_min:
                cc_max = (ze_max + off_limit)*10**13
                cc_min = (ze_min - off_limit)*10**13

            # Set up rounded contour values before shifting back
            diff_string = str(int(cc_max-cc_min))[:3]
            if diff_string[-1] == '0':
                diff_string = diff_string[:-1]
            # Make sure number of contour lines is not too small or large
            cc_nr = int(diff_string) + 1
            if cc_nr > 100:
                cc_nr = (cc_nr - 1)/4 + 1
            elif cc_nr > 50:
                cc_nr = (cc_nr - 1)/2 + 1
            if cc_nr < 10:
                cc_nr = (cc_nr - 1)*4 + 1
            elif cc_nr < 20:
                cc_nr = (cc_nr - 1)*2 + 1

            cc_max /= 10**13
            cc_min /= 10**13
            # Sometimes this rounding setup fails, in which case we
            # go back to the unrounded setup
            if cc_max < ze_max or cc_min > ze_min:
                cc_max = ze_max + off_limit
                cc_min = ze_min - off_limit
                cc_nr = 21
            if cc_nr > 80:
                cc_nr = 21

        cc = np.linspace(cc_min, cc_max, int(cc_nr))

        reset_cc = False
        if reset_cc:
            cc = np.linspace(ze_min, ze_max, 11)

        # Save plots for each frame
        for frame, ts in enumerate(ts_list):
            if frame in nodata:
                continue
            plt.clf()
            plt.figure(figsize=(10, 5))
            cf = plt.contourf(xi*r2d, yi*r2d, ze[frame,:,:], cc, cmap=magma)
            plt.colorbar(cf, cmap=magma)
            cl = plt.contour(xi*r2d, yi*r2d, ze[frame,:,:], cc, linewidths=0.5, colors='k')
            plt.xlabel(xax)
            plt.ylabel(yax)
            #plt.title('Maximum Normalised Error: %1.3f' % maxerror[frame])
            plt.title('Error: L_2: '  +"{:.2e}".format(l2error[frame])+ ', L_inf: ' +"{:.2e}".format(maxerror[frame]))
            out_file_name = plotpath + "/" "error_" + field + "_" + ts + ".png"
            #out_file_name = plotpath + "/" "swe_vortex_64x64_" + field + "_" + ts + ".png"
            #out_file_name = plotpath + "/" "mountain_C24_" + field + "_" + ts + ".png"
            plt.savefig(out_file_name, bbox_inches='tight', dpi=150)


if __name__ == "__main__":
    # Use plot_swe with commands such as
    # python plot_error.py shallow_water $PWD P:32:1 buoyancy:geopot:wind:q 20:100 $PWD
    # For files of name form shallow_water_nodal_buoyancy_T000080.m.m
    # With files stored in folder $PWD,
    # Planar resolution of 32x32, side length 1 (for sphere simply use 'C' instead)
    # Plotting fields buoyancy, geopot, wind, q
    # For time steps at multiples of 20 up to time step 100.

    try:
        config, datapath, nx, fields, timesteps, plotpath = sys.argv[1:8]
    except ValueError:
        print("Usage: {0} <file_stem_name> <datapath>"
              "<'P' if plane, with resolution and extent, 'C' if sphere> <fields_list>"
              "<timestep_list or dumpfreq:tmax> <plotpath>".format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    extend_components_q, extend_components_wind = None, None
    field_list = fields.split(':')
    for field in field_list:
        if field not in ['depth', 'geopot', 'buoyancy', 'wind', 'q']:
            raise IOError('Must choose fields from geopot, buoyancy, wind, and q')
        if field == 'wind':
            extend_components_wind = True

    if extend_components_wind:
        field_list.remove('wind')
        field_list.extend(['wind1', 'wind2'])

    # Split out the list of timesteps
    if timesteps[0] == 'T':
        ts_list = timesteps.split(':')
    else:
        dfr, tmax = [int(strg) for strg in timesteps.split(':')]
        ts_list, c_ = [], 0
        while dfr*c_ <= tmax:
            entry = str(dfr*c_)
            while len(entry) < 6:
                entry = '0' + entry
            entry = 'T' + entry
            ts_list.append(entry)
            c_ += 1

    # Update datapath to include configuration
    datapath += "/" + config + "_nodal_"

    make_figures(datapath, plotpath, nx, field_list, ts_list)
