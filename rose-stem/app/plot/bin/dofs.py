#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Python script to plot xz slices along the y=0 and xy slices on a specified level 

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
import math
import sys

from read_data import read_nodal_data

levels = None
data = None

       
def make_figure(plotpath, field, component, timestep):

  val_col = 'c' + str(component)

  fig = plt.figure(figsize=(10,5))
  nz = len(levels)
  for p in range(len(levels)):
    p_data = data.loc[data['level'] == levels[p]]
    g = p/(2.0*nz) + 0.5
    plt.plot(p_data[val_col].values,color=(g,g,g),linewidth=2)

  plt.title('max: %e, min: %e'%(np.amax(data[val_col].values),np.amin(data[val_col].values)))
  out_file_name = plotpath + "/" "dofs_" + field + "_" + timestep +  ".png"
  plt.savefig(out_file_name , bbox_inches='tight')

 
if __name__ == "__main__":

  try:
    config, datapath, fields, timesteps, plotpath = sys.argv[1:9]
  except ValueError:
    print("Usage: {0} <datapath> <field_names> <timestep_list> <plotpath>".format(sys.argv[0]))
    exit(1)

  # Split out the list of fields
  field_list = fields.split(':')


  # Split out the list of timesteps
  ts_list = timesteps.split(':')

  for field in field_list: 

    for ts in ts_list:

      filestem =  datapath + "/" + config + "_nodal_" + field + "_" + ts + "*"     

      data = read_nodal_data(filestem, 3, 1)
      levels = data.level.unique()

      # Only try to plot if we found some files for this timestep
      if len(levels) > 0:
        make_figure(plotpath,field, 1, ts)

