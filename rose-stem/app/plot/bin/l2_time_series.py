#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# example plot string:
# python l2_time_series.py data/control/ data/pert/  plots T000072:T000144:T000216:T000288:T000360:T000432 theta:exner
'''
Plots l2 norm of the difference field between two run files
'''
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
# Need to set a non-interactive backend for suites
matplotlib.use('Agg')  # noqa: E402
import sys

import numpy as np
import matplotlib.pyplot as plt

from read_data import read_ugrid_data

# file name variables:
data_init  = 'lfric_initial.nc'
data_run   = 'lfric_diag.nc'

def make_figures(data_loc1, data_loc2, plotpath, fields, tsteps):

    fig, ax = plt.subplots(figsize=(20, 10))
    fig.suptitle(r'Normalised $L_2$ norm of the difference field between ' +
                 'two run files')

    theta_min = np.inf
    min_time  = 0
    theta_max = -np.inf
    max_time  = 0

    if fields is None:
            fields = ['theta', 'exner']

    for field in fields:
        if 0.0 in tsteps:
            cube1_init = read_ugrid_data(f'{data_loc1}/{data_init}', field)
            cube2_init = read_ugrid_data(f'{data_loc2}/{data_init}', field)
        cube1 = read_ugrid_data(f'{data_loc1}/{data_run}', field)
        cube2 = read_ugrid_data(f'{data_loc2}/{data_run}', field)

        l2_values = []
        for t in tsteps:
            if t == 0:
                plot_data = np.zeros(cube1_init.data[:].shape)
                plot_data[:, :] = cube1_init.data[:] - \
                                  cube2_init.data[:]
                control_data = np.zeros(cube1_init.data[:].shape)
                control_data[:, :] = cube1_init.data[:]
            else:
                time_points = cube1.coord('time').points
                tstep_index = np.where(time_points == t)[0][0]

                plot_data = np.zeros(cube1.data[tstep_index, :].shape)
                plot_data[:, :] = cube1.data[tstep_index, :] - \
                                cube2.data[tstep_index, :]
                control_data = np.zeros(cube1.data[tstep_index, :].shape)
                control_data[:, :] = cube1.data[tstep_index, :]

            l2_value = np.linalg.norm(plot_data) /\
                        np.linalg.norm(control_data)
            l2_values.append(l2_value)

            if field == 'theta':
                if np.amin(abs(plot_data)) < theta_min:
                    theta_min = np.amin(abs(plot_data))
                    min_time  = t
                if np.amax(abs(plot_data)) > theta_max:
                    theta_max = np.amax(abs(plot_data))
                    max_time  = t

        ax.plot(tsteps, l2_values, label=field)

    if 'theta' in fields:
        ax.set_title(r'$min(\theta_{\mathrm{pert}} - \theta)$ = ' +
                     f'{theta_min:.3e} at time {min_time}\n' +
                     r'$max(\theta_{\mathrm{pert}} - \theta)$ = ' +
                     f'{theta_max:.3e} at time {max_time}')

    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$L_2(f_{\mathrm{pert}} - f)/L_2(f)$')
    pngfile = '%s/l2error.png' % plotpath
    plt.savefig(pngfile)
    plt.close()


if __name__ == "__main__":

    try:
        args = sys.argv[:]
        data_loc1, data_loc2, plotpath, tsteps = args[1:5]
        field_list = None
        tsteps_list = tsteps.split(':')
        tsteps_list = [float(tstep[1:])*10**3 for tstep in tsteps_list]
        if len(args[:]) > 5:
            field_list = args[5].split(':')
    except ValueError:
        print("Usage: {0} <data_loc1> <data_loc2> <plotpath> <tsteps> [<fields_list>] "
              .format(sys.argv[0]))
        exit(1)

    make_figures(data_loc1, data_loc2, plotpath, field_list, tsteps_list)
