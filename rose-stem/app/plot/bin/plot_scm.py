#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
''' Quick plot of for lfric_atm scm output '''

# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

# Note non-PEP8 collecting of imports as the backend needs to be
# set before we import iris.
import iris
import matplotlib.pyplot as plt

def load_cube_by_varname(filename, var):
   variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var))
   return iris.load_cube(filename, constraint=variable_constraint)

def do_plot(datapath, plotfield, plotpath='.'):
    ''' Do the plotting using data from datapath. Send output to plotpath '''

    print(f"Plotting {plotfield}")
    lfric = load_cube_by_varname(datapath, plotfield)
    lfric = lfric[:, :, 0]

    plt.figure(figsize=(15, 10))
    for n, time in enumerate([0, 2, 4, 8, 12, 14]):
        plt.subplot(2, 3, n+1)
        try:
           # first try wtheta fields
           plt.plot(lfric.data[time, :],
                    lfric.coord('full_levels').points[:],
                    linewidth=2)
        except:
           # then w3 fields
           plt.plot(lfric.data[time, :],
                    lfric.coord('half_levels').points[:],
                    linewidth=2)

        plt.xlabel(plotfield)
        plt.ylabel('Model Level Number')
        plt.title('Timestep = '+str(time+1))

    plt.savefig(plotpath+'/'+plotfield+'.png',
                bbox_inches='tight')
    plt.close()


def do_time_plot(datapath, plotfield, plotpath='.'):
    ''' Do the plotting using data from datapath. Send output to plotpath '''

    lfric = load_cube_by_varname(datapath, plotfield)
    lfric = lfric[:, 0]

    plt.figure(figsize=(15, 10))
    plt.plot(lfric.coord('time').points[:]/3600.0,
             lfric.data[:],
             linewidth=2)
    plt.xlabel('Time (hours)')
    plt.ylabel(plotfield)
    plt.savefig(plotpath+'/'+plotfield+'.png',
                bbox_inches='tight')
    plt.close()


if __name__ == "__main__":

    import sys
    try:
        opts = [opt for opt in sys.argv[1:] if opt.startswith('-')]
        args = [arg for arg in sys.argv[1:] if not arg.startswith('-')]
        datapath, plotpath = args[0:2]
        extra_plots_ukca = '-ukca' in opts
    except ValueError:
        print("Usage: {0} [-ukca] <datapath> <plotpath>".format(sys.argv[0]))
        exit(1)
    do_plot(datapath, 'theta', plotpath)
    do_plot(datapath, 'm_v',   plotpath)
    do_plot(datapath, 'm_cl',  plotpath)
    do_plot(datapath, 'm_cf',  plotpath)
    do_plot(datapath, 'u_in_w3', plotpath)
    do_plot(datapath, 'combined_cloud_amount', plotpath)
    do_time_plot(datapath, 'sw_down_surf', plotpath)
    do_time_plot(datapath, 'lw_up_toa', plotpath)
    do_time_plot(datapath, 'ls_prec', plotpath)
    do_time_plot(datapath, 'conv_rain', plotpath)
    if extra_plots_ukca:
        do_plot(datapath, 'h2o2', plotpath)
        do_plot(datapath, 'dms', plotpath)
        do_plot(datapath, 'so2', plotpath)
        do_plot(datapath, 'h2so4', plotpath)
        do_plot(datapath, 'dmso', plotpath)
        do_plot(datapath, 'monoterpene', plotpath)
        do_plot(datapath, 'secondary_organic', plotpath)
        do_plot(datapath, 'n_nuc_sol', plotpath)
        do_plot(datapath, 'n_ait_sol', plotpath)
        do_plot(datapath, 'n_acc_sol', plotpath)
        do_plot(datapath, 'n_cor_sol', plotpath)
        do_plot(datapath, 'n_ait_ins', plotpath)
        do_plot(datapath, 'n_acc_ins', plotpath)
        do_plot(datapath, 'n_cor_ins', plotpath)
        do_plot(datapath, 'nuc_sol_su', plotpath)
        do_plot(datapath, 'nuc_sol_om', plotpath)
        do_plot(datapath, 'ait_sol_su', plotpath)
        do_plot(datapath, 'ait_sol_bc', plotpath)
        do_plot(datapath, 'ait_sol_om', plotpath)
        do_plot(datapath, 'acc_sol_su', plotpath)
        do_plot(datapath, 'acc_sol_bc', plotpath)
        do_plot(datapath, 'acc_sol_om', plotpath)
        do_plot(datapath, 'acc_sol_ss', plotpath)
        do_plot(datapath, 'acc_sol_du', plotpath)
        do_plot(datapath, 'cor_sol_su', plotpath)
        do_plot(datapath, 'cor_sol_bc', plotpath)
        do_plot(datapath, 'cor_sol_om', plotpath)
        do_plot(datapath, 'cor_sol_ss', plotpath)
        do_plot(datapath, 'cor_sol_du', plotpath)
        do_plot(datapath, 'ait_ins_bc', plotpath)
        do_plot(datapath, 'ait_ins_om', plotpath)
        do_plot(datapath, 'acc_ins_du', plotpath)
        do_plot(datapath, 'cor_ins_du', plotpath)
        do_plot(datapath, 'cloud_drop_no_conc', plotpath)
        do_plot(datapath, 'sw_aer_optical_depth_rts', plotpath)
        do_plot(datapath, 'lw_aer_optical_depth_rts', plotpath)
        do_plot(datapath, 'aod_ukca_ait_sol', plotpath)
        do_plot(datapath, 'aaod_ukca_ait_sol', plotpath)
        do_plot(datapath, 'aod_ukca_acc_sol', plotpath)
        do_plot(datapath, 'aaod_ukca_acc_sol', plotpath)
        do_plot(datapath, 'aod_ukca_cor_sol', plotpath)
        do_plot(datapath, 'aaod_ukca_cor_sol', plotpath)
        do_plot(datapath, 'aod_ukca_ait_ins', plotpath)
        do_plot(datapath, 'aaod_ukca_ait_ins', plotpath)
        do_plot(datapath, 'aod_ukca_acc_ins', plotpath)
        do_plot(datapath, 'aaod_ukca_acc_ins', plotpath)
        do_plot(datapath, 'aod_ukca_cor_ins', plotpath)
        do_plot(datapath, 'aaod_ukca_cor_ins', plotpath)
