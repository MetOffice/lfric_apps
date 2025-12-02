#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2017 Met Office. All rights reserved.
# For further details please refer to the file LICENCE which you should have
# received as part of this distribution.
##############################################################################

'''
Python script to plot all levels from a Dynamo output file. Levels are determined
from the data as they are different for different fields.

This version takes nodal format output files and interpolates onto a regular
grid.

This version stitches together a directory of files and extracts all levels
so it can work in the serial case where there is one file or the parallel
case where there is a file for each processor.

This version is for plotting under suites and accepts command line args
for the field and timestep to plot. It also plots to file rather than to screen 
'''
from __future__ import absolute_import
from __future__ import print_function
import numpy as np

# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys

from read_data import read_nodal_data

levels = None
data = None


def make_figure(field, nx, ny, component, timestep):

  val_col = 'c' + str(component)

  fig = plt.figure(figsize=(15,10))
  plt.plot()

#  for p in range(len(levels)):
  for p in range(1):

    p_data = data.loc[data['level'] == levels[p]]

    # get min and max of x,y data for plot axes
    xmin = p_data['x'].min()
    xmax = p_data['x'].max()
    ymin = p_data['y'].min()
    ymax = p_data['y'].max()

    # Size of regular grid
    nx = int(nx)
    ny = int(ny)

    # Using reshape of numpy array
    w3fieldi = (p_data[val_col].values).reshape((ny, nx))
    xi = (p_data['x'].values).reshape((ny, nx))
    yi = (p_data['y'].values).reshape((ny, nx))

    cf = plt.pcolormesh(xi, yi, w3fieldi,cmap=cm.hsv,vmin=-1.0,vmax=10.0) # use this for biperiodic grid
    plt.colorbar(cf,cmap=cm.hsv)

    out_file_name = plotpath + "/" "biperiodic-cosmic_" + field + "_" + timestep +  ".png"
    plt.savefig(out_file_name , bbox_inches='tight')

  print(out_file_name)


if __name__ == "__main__":

  try:
    config, datapath, nx, ny, fields, timesteps, plotpath = sys.argv[1:8]
  except ValueError:
    print("Usage: {0} <file_stem_name> <datapath> <nx> <ny> <field_names> <timestep_list> <plotpath>".format(sys.argv[0]))
    exit(1)

  # Split out the list of fields
  field_list = fields.split(':')

  # Split out the list of timesteps
  ts_list = timesteps.split(':')

  for field in field_list:

    for ts in ts_list:

      filestem =  datapath + "/" + config + "_nodal_" + field + "_" + ts + "*"

      if field in ['u','xi']:
        # Vector field - plot w component by default
        component = 3
        data = read_nodal_data(filestem, 3, component)
      else:
        # Scalar field - only one component to plot
        component = 1
        data = read_nodal_data(filestem, 1, component)

      # Sort the data (needed to be able to reshape and not regrid)
      data = data.sort_values(['y','x','z'])

      levels = np.sort(data.level.unique())

      make_figure(field, nx, ny, component, ts)


