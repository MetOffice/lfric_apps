#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Python script to plot xz slices along the equator and
a hovmoller plot along the equator at height = 5 km of the dcmip test case

This version takes nodal format output files and
interpolates onto a regular grid.

Filename hardcoded.

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

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.interpolate import griddata

import sys

from read_data import read_nodal_data

levels = None
data = None
       
def make_xz_figure(plotpath, field, component, timestep):

  val_col = 'c' + str(component)

  plt.figure()

   # get min and max of x,y data for plot axes

  min_lev = min(levels)

  deltaz = data.loc[data['level'] == (min_lev+1)][val_col].min() - data.loc[data['level'] == (min_lev)][val_col].min()
  xmin = data.loc[data['level'] == min_lev]['x'].min()
  xmax = data.loc[data['level'] == min_lev]['x'].max()
  ymin = data.loc[data['level'] == min_lev]['y'].min()
  ymax = data.loc[data['level'] == min_lev]['y'].max()

  zmin = 0.0
  zmax = 10000.0
  
  r2d = 180/np.pi;
  nx,ny,nz = 80*4,2,11

  #create 2D plot
  x2d = np.linspace(xmin, xmax, nx)
  z2d = np.linspace(zmin, zmax, nz)
  y2d = 0.0
  xi, yi = np.meshgrid(x2d, y2d)  
  zi = np.zeros([1,nx,len(levels)])

  for p in range(len(levels)):
    p_data = data.loc[data['level'] == levels[p]]
    zi[:,:,p] = griddata((p_data['x'].values, p_data['y'].values), p_data[val_col].values, (xi, yi), method='linear')

  yi, xi = np.meshgrid(z2d, x2d) 
  dz = np.zeros([nx,len(levels)])
  for i in range(nx):
    dz[i,:] = zi[0,i,:] - zi[0,0,:]

  c_map = cm.summer
  cc = np.linspace(-0.12,0.12,13)
  cf = plt.contourf(xi *r2d, yi , dz, cc, cmap=c_map)
  plt.colorbar(cf,  cmap=c_map)
  cl = plt.contour(xi * r2d, yi, dz, cc, linewidths=0.5,colors='k')
  plt.title('max: %e, min: %e'%(np.max(dz),np.min(dz)))
  out_file_name = plotpath + "/" "dcmip301_xz_" + field + "_" + timestep +  ".png"
  plt.savefig(out_file_name , bbox_inches='tight')

def make_hovmoller_figure(datapath, plotpath, field, component, step):

  val_col = 'c' + str(component)

  fig = plt.figure(figsize=(15,10))

  r2d = 180/np.pi;
  nx,ny,nz,nt = 90*4,2,1,11

  t2d = np.zeros([nx,nt])


  for t in range(0,nt):

    timestep = "T"+str(step*t).zfill(6)
    filestem =  datapath + "/diagGungho_nodal_" + field + "_" + timestep + "*"

    data = read_nodal_data(filestem, 1, 1)
    levels = data.level.unique()
 
    # get min and max of x,y data for plot axes

    min_lev = min(levels)

    deltaz = data.loc[data['level'] == (min_lev+1)][val_col].min() - data.loc[data['level'] == (min_lev)][val_col].min()
    xmin = data.loc[data['level'] == min_lev]['x'].min()
    xmax = data.loc[data['level'] == min_lev]['x'].max()
    ymin = data.loc[data['level'] == min_lev]['y'].min()
    ymax = data.loc[data['level'] == min_lev]['y'].max()
    zmin = 0.0
    zmax = 10000.0
  
    #create 2D plot
    x2d = np.linspace(xmin, xmax, nx)
    z2d = 5
    y2d = 0.0
    xi, yi = np.meshgrid(x2d, y2d)   
    zi = np.zeros([1,nx])
    p = 5
    p_data = data.loc[data['level'] == levels[p]]
    zi[:,:] = griddata((p_data['x'].values, p_data['y'].values), p_data[val_col].values, (xi, yi), method='linear')

    for i in range(nx):
      t2d[i,t] = zi[0,i] - zi[0,0]

  time =  np.linspace(0, 3600, nt)
  ti, xi = np.meshgrid(time, x2d)  
  c_map = cm.summer 
  cc = np.linspace(-0.12,0.12,21)
  cf = plt.contourf(xi *r2d, ti , t2d, cc, cmap=c_map)
  plt.colorbar(cf,  cmap=c_map)
  cl = plt.contour(xi * r2d, ti, t2d, cc, linewidths=0.5,colors='k')
  plt.title('max: %e, min: %e'%(np.max(t2d),np.min(t2d)))

  # Add in constant speeds
  Lz = 20000.0
  N =  0.01
  u1 = (20.0 - N*Lz/(2.0*np.pi))
  u2 = (20.0 + N*Lz/(2.0*np.pi))

  #Convert to degrees
  a=6371229.0/125.0
  u1 = u1/a * r2d
  u2 = u2/a * r2d

  x0 = -1.0/3.0*np.pi * r2d

  plt.plot([x0,x0+u1*3600.0],[0,3600.0],'k',linewidth=4)
  plt.plot([x0,x0+u2*3600.0],[0,3600.0],'k',linewidth=4)

  out_file_name = plotpath + "/" "dcmip301_hovmoller_" + field + "_" + timestep +  ".png"
  plt.savefig(out_file_name , bbox_inches='tight')


if __name__ == "__main__":

  try:
    config, datapath, fields, timesteps, step, plotpath = sys.argv[1:7]
  except ValueError:
    print("Usage: {0} <file_stem_name> <datapath> <field_names> <timestep_list> <step> <plotpath>".format(sys.argv[0]))
    exit(1)

  # Split out the list of fields
  field_list = fields.split(':')

  # Split out the list of timesteps
  ts_list = timesteps.split(':')

  for field in field_list:

    for ts in ts_list:


      filestem =  datapath + "/" + config + "_nodal_" + field + "_" + ts + "*"

      data = read_nodal_data(filestem, 1, 1)
      levels = data.level.unique()

      # Only try to plot if we found some files for this timestep
      if len(levels) > 0:
        make_xz_figure(plotpath,field, 1, ts)


    make_hovmoller_figure(datapath, plotpath, field, 1, int(step))


