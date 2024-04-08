#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2018 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution
##############################################################################
'''
Script plotting time series of maxima of absolute value of velocities
'''

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # noqa: E402

import matplotlib.pyplot as plt

'''

Usage example:

python time_series_max.py datapath plotpath plot_ref num_proc


'''

if __name__ == "__main__":

    try:
        config, datapath, plotpath, plot_ref, num_proc = sys.argv[1:6]
    except ValueError:
        print(f"Usage: {sys.argv[0]} <file_stem_name> <datapath> <plotpath> "
              f"<plot_ref> <num_proc>")
        sys.exit(1)

    fig = plt.figure(figsize=(30, 8))

    if num_proc == '1':
        data_u = np.genfromtxt(datapath + "/" + config + "_nodal_minmax_u.m")
    else:
        data_u = np.genfromtxt(datapath +
                               "/" + config + "_nodal_minmax_u.Rank000000.m")
        for n in range(1, int(num_proc)):  # works up to 10 cores only!
            data_u = \
                np.concatenate((data_u, np.genfromtxt(
                    datapath + "/" + config +
                    "_nodal_minmax_u.Rank00000" + str(n) + ".m")),
                    axis=1)

    JJ = len(data_u)//3  # Number of time steps

    zvel_min = data_u[::3, -1]
    zvel_max = data_u[::3, -2]

    mvel_min = data_u[1::3, -1]
    mvel_max = data_u[1::3, -2]

    max_abs_zvel = np.zeros((JJ))
    max_abs_mvel = np.zeros((JJ))

    for i in range(len(zvel_min)):
        max_abs_zvel[i] = abs(zvel_min[i]) \
            if abs(zvel_min[i]) > abs(zvel_max[i]) \
            else abs(zvel_max[i])

    for i in range(len(mvel_min)):
        max_abs_mvel[i] = abs(mvel_min[i]) \
            if abs(mvel_min[i]) > abs(mvel_max[i]) \
            else abs(mvel_max[i])

    vvel_min = data_u[2::3, -1]
    vvel_max = data_u[2::3, -2]

    max_abs_vvel = np.zeros((JJ))

    for i in range(len(vvel_min)):
        max_abs_vvel[i] = abs(vvel_min[i]) \
            if abs(vvel_min[i]) > abs(vvel_max[i]) \
            else abs(vvel_max[i])

    x_coords = list(range(0, JJ))

    if plot_ref == '1':
        # Plotting reference data

        if num_proc == '1':
            data_u = np.genfromtxt(datapath + "/" + config + "_nodal_minmax_u"
                                              "_ref.m")
        else:
            data_u = np.genfromtxt(datapath + "/" + config + "_nodal_minmax_u"
                                              ".Rank000000_ref.m")

        zvel_min_ref = data_u[::3, -1]
        zvel_max_ref = data_u[::3, -2]

        mvel_min_ref = data_u[1::3, -1]
        mvel_max_ref = data_u[1::3, -2]

        max_abs_zvel_ref = np.zeros((JJ))
        max_abs_mvel_ref = np.zeros((JJ))

        for i in range(len(zvel_min_ref)):
            max_abs_zvel_ref[i] = abs(zvel_min_ref[i]) \
                if abs(zvel_min_ref[i]) > abs(zvel_max_ref[i]) \
                else abs(zvel_max_ref[i])

        for i in range(len(mvel_min_ref)):
            max_abs_mvel_ref[i] = abs(mvel_min_ref[i]) \
                if abs(mvel_min_ref[i]) > abs(mvel_max_ref[i]) \
                else abs(mvel_max_ref[i])

        vvel_min_ref = data_u[2::3, -1]
        vvel_max_ref = data_u[2::3, -2]

        max_abs_vvel_ref = np.zeros((JJ))

        for i in range(len(vvel_min_ref)):
            max_abs_vvel_ref[i] = abs(vvel_min_ref[i]) \
                if abs(vvel_min_ref[i]) > abs(vvel_max_ref[i]) \
                else abs(vvel_max_ref[i])

    # Plotting all components of velocity
    ax_max_zvel = fig.add_subplot(1, 3, 1)
    s_zvel = ax_max_zvel.plot(x_coords, max_abs_zvel, 'b')
    if plot_ref == '1':
        s_zvel = ax_max_zvel.plot(x_coords[::100], max_abs_zvel_ref[::100],
                                  'k', linewidth=2)
    ax_max_zvel.get_xaxis().get_major_formatter().set_scientific(False)
    ax_max_zvel.set_xlabel('Time step #')
    ax_max_zvel.set_ylabel(r'$\max||u||_\infty [m/s]$')
    plt.title('Max |u| = %e' % np.max(max_abs_zvel))
    ax_max_zvel.legend(loc='center right', fancybox=True, shadow=True)
    plt.tight_layout()

    ax_max_mvel = fig.add_subplot(1, 3, 2)
    s_mvel = ax_max_mvel.plot(x_coords, max_abs_mvel, 'b')
    if plot_ref == '1':
        s_mvel = ax_max_mvel.plot(x_coords[::100], max_abs_mvel_ref[::100],
                                  'k', linewidth=2)
    ax_max_mvel.get_xaxis().get_major_formatter().set_scientific(False)
    ax_max_mvel.set_xlabel('Time step #')
    ax_max_mvel.set_ylabel(r'$\max||v||_\infty [m/s]$')
    ax_max_mvel.set_yscale("log")
    plt.title('Max |v| = %e' % np.max(max_abs_mvel))
    ax_max_mvel.legend(loc='center right', fancybox=True, shadow=True)
    plt.tight_layout()

    ax_max_vvel = fig.add_subplot(1, 3, 3)
    s_vvel = ax_max_vvel.plot(x_coords, max_abs_vvel, 'b')
    if plot_ref == '1':
        s_vvel = ax_max_vvel.plot(x_coords[::100], max_abs_vvel_ref[::100],
                                  'k', linewidth=2)
    ax_max_vvel.set_xlabel('Time step #')
    ax_max_vvel.set_ylabel(r'$\max||w||_\infty [m/s]$')
    plt.title('Max |w| = %e' % np.max(max_abs_vvel))
    ax_max_vvel.set_yscale("log")
    ax_max_vvel.legend(loc='center right', fancybox=True, shadow=True)

    plt.tight_layout()

    plt.tight_layout()
    png_file_name = plotpath + "/vel_time_series_0_" + str(JJ-1) + "ts.png"
    plt.savefig(png_file_name)
