#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.interpolate import griddata

import math

import sys

from read_data import read_nodal_data

levels = None
data = None
data0 = None

def make_figure(plotpath, nx, ny, field, component, timestep):

  val_col = 'c' + str(component)

  fig = plt.figure(figsize=(7,7))

  # Sort levels in ascending order, this is needed for high order spaces
  sorted_levels = sorted(levels)
  l2h = np.zeros(len(levels))
  for i in range(len(levels)):
    for j in range(len(levels)):
      if ( sorted_levels[i] == levels[j] ):
        l2h[i] = j

   # get min and max of x,y data for plot axes

  min_lev = min(levels)

  xmin = data.loc[data['level'] == min_lev]['x'].min()
  xmax = data.loc[data['level'] == min_lev]['x'].max()
  ymin = data.loc[data['level'] == min_lev]['y'].min()
  ymax = data.loc[data['level'] == min_lev]['y'].max()

  zmin = 0.0
  zmax = 10000.0

  r2d = 1.0/1000.0;

  nx = int(nx)
  ny = int(ny)
  nz = len(levels)


  mid_z = math.floor(len(levels)/2)
  wi = np.zeros([ny,nx,len(levels)])
  wi0= np.zeros([ny,nx,len(levels)])
  wi_bg = np.zeros(len(levels))

  for p in range(len(levels)):
    pp = int(l2h[p])
    p_data = data.loc[data['level'] == levels[pp]]
    p_data0 = data0.loc[data0['level'] == levels[pp]]
    wi[:,:,p] = (p_data[val_col].values).reshape((ny, nx))
    wi0[:,:,p] = (p_data0[val_col].values).reshape((ny, nx))


  # subtracting background profile
  for l in range(len(levels)):
    wi_bg[l] = wi0[0,0,l]
    wi[:,:,l] -= wi_bg[l]

  dw = np.zeros([nx,len(levels)])
  for i in range(nx):
    dw[i,:] = wi[0,i,:]


  # create meshgrid to get z_i and x_i for plotting
  x2d = np.linspace(xmin, xmax, nx)
  z2d = np.linspace(zmin, zmax, nz)
  z_i, x_i = np.meshgrid(z2d, x2d)

  ax = fig.add_subplot(2,1,1)

  matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
  c_map = cm.summer
  cc = np.linspace(-0.003,0.003,13)
  cf = ax.contourf(x_i *r2d, z_i * r2d, dw, cc, cmap=c_map)
  plt.colorbar(cf,  cmap=c_map)
  cl = ax.contour(x_i * r2d, z_i*r2d, dw, cc, linewidths=0.5,colors='k')
  ax.set_title('max: %2.4e, min: %2.4e'%(np.max(dw),np.min(dw)))
  ax.set_xlabel('x (km)')
  ax.set_ylabel('z (km)')

  bx = fig.add_subplot(2,1,2)
  bx.plot(x_i*r2d, dw[:,int(round(mid_z))], color='k')
  bx.set_xlabel('x (km)')
  bx.set_ylabel('theta (K)')
  bx.set_ylim([-0.002, 0.003])
 
  out_file_name = plotpath + "/" + 'gravity_wave' + "_" + timestep + ".png"
  plt.savefig(out_file_name, bbox_inches='tight')

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
    # Create initial data
    filestem =  datapath + "/" + config + "_nodal_" + field + "_" + "T000000" + "*"

    data0 = read_nodal_data(filestem, 1, 1)

    for ts in ts_list:

      filestem =  datapath + "/" + config + "_nodal_" + field + "_" + ts + "*"

      data = read_nodal_data(filestem, 1, 1)
      levels = data.level.unique()

      # Only try to plot if we found some files for this timestep
      if len(levels) > 0:
        make_figure(plotpath, nx, ny, field, 1, ts)

