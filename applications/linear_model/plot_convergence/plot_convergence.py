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

    norm_df.plot.scatter(x='gamma', y='norm', loglog=True, xlim=(10**0, 10**5),
                         ylim=(10**-5, 10**0), ax=axes, color=color,
                         marker=shape)


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
    plot_data(directory + filename, axs, 'gamma_total', 'black', 'x')
    plot_data(directory + filename, axs, 'gamma_mr', 'm', '*')

    expected_x = [10**0, 10**5]
    expected_y = [10**-5, 10**0]
    plt.plot(expected_x, expected_y)

    plt.legend(['Expected gradient (linear)', 'density', 'momentum',
                'exner pressure', 'potential temperature', 'total',
                'moisture mixing ratios'],
               loc='lower right')
    plt.xlabel('Gamma')
    plt.ylabel('Relative error')
    plt.title('Validity of the tangent linear model')

    plt.show()


if __name__ == "__main__":

    DATA_DIRECTORY = os.getcwd()+'/'
    DATA_FILENAME = 'outfile'

    make_plot(DATA_DIRECTORY, DATA_FILENAME)
