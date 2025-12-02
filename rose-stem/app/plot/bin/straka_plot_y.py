#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Basic python script to plot the x-z profile minus a constant state of 300
from a Dynamo output file.

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

levels = None
data = None


def make_figure(plotpath, nx, ny, field, component, timestep):

  val_col = 'c' + str(component)

  slice_fig = plt.figure(figsize=(15,10))

  # get min and max of x,y data for plot axes
  min_lev = min(levels)

  xmin = data.loc[data['level'] == min_lev]['x'].min()
  xmax = data.loc[data['level'] == min_lev]['x'].max()
  ymin = data.loc[data['level'] == min_lev]['y'].min()
  ymax = data.loc[data['level'] == min_lev]['y'].max()

  zmin = 0.0
  zmax = 6400.0

  r2d = 1.0/1000.0;

  nx = int(nx)
  ny = int(ny)
  nz = len(levels)


  zi = np.zeros([ny,nx,len(levels)])

  for p in range(len(levels)):
    p_data = data.loc[data['level'] == levels[p]]
    zi[:,:,p] = (p_data[val_col].values).reshape((ny, nx))

  # create meshgrid to get x_i and y_i for plotting
  y2d = np.linspace(ymin, ymax, ny)
  z2d = np.linspace(zmin, zmax, nz)
  y_i, x_i = np.meshgrid(z2d, y2d) 

  dz = np.zeros([ny,len(levels)])
  for i in range(ny):
    dz[i,:] = zi[i,0,:] - 300.0

  matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
  c_map = cm.summer
  cc = np.linspace(-16,-1,16)
  cf = plt.contourf(x_i *r2d, y_i * r2d, np.round(dz,10), cc, cmap=c_map)
  cl = plt.contour(x_i * r2d, y_i*r2d, np.round(dz,10), cc, linewidths=1.0,colors='k', linestyle="", extend='min')
  plt.axis([0, 16, 0, 5])
  plt.xlabel("y (km)")
  plt.ylabel("z (km)")
  plt.title('max: %2.4e, min: %2.4e'%(np.max(dz),np.min(dz)))
  plt.colorbar(cf,  cmap=c_map)

  out_file_name = plotpath + "/" + 'straka_y' + "_" + timestep +  ".png"
  slice_fig.savefig(out_file_name , bbox_inches='tight')

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

      data = read_nodal_data(filestem, 1, 1)

      # Sort the data (needed to be able to reshape and not regrid)
      data = data.sort_values(['y','x','z'])

      levels = np.sort(data.level.unique())

      # Only try to plot if we found some files for this timestep
      if len(levels) > 0:
        make_figure(plotpath, nx, ny, field, 1, ts)

