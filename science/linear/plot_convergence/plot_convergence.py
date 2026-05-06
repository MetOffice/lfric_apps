#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Using the output derived from the test_timesteps integration test with
10 values of gamma, plot a graph of the relative error (linearisation error)
against the size of the perturbation for different prognostic variables
(gamma).
'''
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plot_data(filename, axes, variable, color, shape):
    '''
    Create a data frame with columns gamma (size of the perturbation)
    and norm (the relative error). Plot this data with a log-log axis.
    '''

    datafile = open(filename, 'r')

    norm_line = []
    line = datafile.readline

    while line:
        line = datafile.readline()

        if variable in line:
            split_line = line.replace('norm', '').replace('\n', '').split('=')
            norm_line.append([float(split_line[1]), float(split_line[2])])

    norm_df = pd.DataFrame(norm_line, columns=['gamma', 'norm'])

    datafile.close()
    
    # Square root (as the read data is only the inner product and
    # has not included the square root to give the norm)
    norm_df['norm'] = np.sqrt(norm_df['norm'])

    # Normalise - to give a relative error
    normalise = norm_df['norm'].iloc[4]
    norm_df['norm'] = norm_df['norm'] / normalise

    # Plot the data
    if (CONFIG == 'nwp_gal9'):
        ymin = 10**-2
        ymax = 10**2
        xmin = 10**-3
        xmax = 10**1
    elif (CONFIG == 'semi_implicit'):
        ymin = 10**-2
        ymax = 10**2
        xmin = 10**-1
        xmax = 10**4
    elif (CONFIG == 'runge_kutta'):
        ymin = 10**-3
        ymax = 10**2
        xmin = 10**-1
        xmax = 10**4
    else:
        print(CONFIG+' not listed')

    norm_df.plot.scatter(x='gamma', y='norm', loglog=True,
                         xlim=(xmin, xmax),
                         ylim=(ymin, ymax),
                         ax=axes, color=color,
                         marker=shape, label = variable)

    # Check extremes
    norm_min = norm_df['norm'].min()
    norm_max = norm_df['norm'].max()
    if (norm_min < ymin) :
        print('Warning: Min value is'+ str(norm_min))
    if (norm_max > ymax) :
        print('Warning: Max value is'+ str(norm_max))
    
    # Plot the expected line
    centre = norm_df['gamma'].iloc[4] 
    expected_x = [10**-2 *centre, centre, 100* centre]
    expected_y = [10**-2, 10**0, 100]
    plt.plot(expected_x, expected_y)

    
def make_plot(directory, filename):
    '''
    Plot the data for the different prognostic variables, together with the
    expected gradient ( which is linear ), on the same plot.
    '''

    axs = plt.axes()

    plot_data(directory + filename, axs, 'gamma_rho', 'c', '<')
    plot_data(directory + filename, axs, 'gamma_u', 'r', 'o')
    plot_data(directory + filename, axs, 'gamma_exner', 'b', 's')
    plot_data(directory + filename, axs, 'gamma_theta', 'g', '^')
    plot_data(directory + filename, axs, 'gamma_mr', 'black', 'x')

    plt.legend(loc='lower right')
    plt.xlabel('Gamma')
    plt.ylabel('Relative error')
    plt.title('Validity of the tangent linear model')

    # To show the plot to the screen, uncommment plt.show()
    #plt.show()

    # Save the plot to a file
    plt.savefig(directory + filename + "convergence_plot.png") 


if __name__ == "__main__":

    DATA_DIRECTORY = os.getcwd()+'/'
    DATA_FILENAME = 'outfile'
    CONFIG = os.getenv('CONFIG')
    print(CONFIG)

    make_plot(DATA_DIRECTORY, DATA_FILENAME)
