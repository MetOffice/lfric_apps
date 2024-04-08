#!/usr/bin/env python
''' Quick plot for lfric_atm global output '''

# Need to set a non-interactive backend for suites
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

# Note non-PEP8 collecting of imports as the backend needs to be
# set before we import iris.
import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Index to the fields
varname=0
colbar_min=1
colbar_max=2

# Fields which are available to plot
theta           = ['theta',           292,  295]
m_v             = ['m_v',             6e-3, 8e-3]
m_cl            = ['m_cl',            0,    1e-3]
m_ci            = ['m_ci',            0,    2e-3]
ls_prec         = ['ls_prec',         1,    25]
w_in_wth        = ['w_in_wth',        -1,    3]
u_in_w3         = ['u_in_w3',         -5,    5]
v_in_w3         = ['v_in_w3',         -5,    5]

def load_cube_by_varname(filename, var):
   variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var))
   return iris.load_cube(filename, constraint=variable_constraint)

def do_plot(datapath, plotfield, plotpath='.', plotlevel=0):
    ''' Do the plotting using data from datapath. Send output to plotpath '''

    lfric = load_cube_by_varname(datapath, plotfield[varname])
    if lfric.ndim == 2:
        lfric = lfric[-1]
    else:
        lfric = lfric[-1, plotlevel]
    
    x_coord_name = 'longitude'
    y_coord_name = 'latitude'
    mesh_cube = load_cube_by_varname(datapath, lfric.attributes['mesh'])
    if mesh_cube.attributes['geometry'] == 'planar':
       x_coord_name = 'projection_x_coordinate'
       y_coord_name = 'projection_y_coordinate'
 
    if plotfield[varname] == 'ls_prec':
       import cf_units
       lfric.units = cf_units.Unit('mm s-1')
       lfric.convert_units('mm h-1')

    # Get the x and y co-ordinates
    x_coord = (np.around(lfric.coord(x_coord_name).points, decimals=5))/1000.0
    y_coord = (np.around(lfric.coord(y_coord_name).points, decimals=5))/1000.0
    t_coord = (np.around(lfric.coord('time').points, decimals=2))/3600.0

    # Save the min and max of the data
    field_min = np.around(np.min(lfric.data), decimals=7)
    field_max = np.around(np.max(lfric.data), decimals=7)

    # Reshape data (assuming square mesh)
    npts = len(x_coord)
    nx = np.sqrt(npts)
    nx_int = nx.astype(int)
    cube = np.reshape(lfric.data, (nx_int,nx_int))
    xx = np.reshape(x_coord, (nx_int,nx_int))
    yy = np.reshape(y_coord, (nx_int,nx_int))

    if plotfield[varname] == 'ls_prec':
       levels = np.array([0.125,0.25,0.5,1,2,4,8,16,32,64])
    else:
       spacing = (plotfield[colbar_max] - plotfield[colbar_min]) / 8.0
       levels = np.arange(plotfield[colbar_min],plotfield[colbar_max]+spacing,spacing)

    my_cmap = plt.cm.viridis
    my_norm = mpl.colors.BoundaryNorm(levels, my_cmap.N)

    plt.figure(figsize=(6, 5))
    plot = plt.pcolormesh(xx, yy, cube, norm=my_norm, rasterized=True, shading='auto')
    plt.colorbar(plot,orientation='vertical')

    plt.title(plotfield[varname]+', level = '+str(plotlevel)
                                +', time = '+str(np.around(t_coord[0], decimals=2))+' hr' 
                                +'\n min = '+str(field_min)
                                +', max = '+str(field_max) )
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')

    plt.savefig(plotpath+'/'+plotfield[varname]+'_level'+str(plotlevel)+'.png', bbox_inches='tight')


if __name__ == "__main__":

    import sys
    try:
        datapath, plotpath = sys.argv[1:3]
    except ValueError:
        print("Usage: {0} <datapath> <plotpath>".format(sys.argv[0]))
        exit(1)
    do_plot(datapath, theta,           plotpath, plotlevel=20)
    do_plot(datapath, m_v,             plotpath, plotlevel=20)
    do_plot(datapath, m_cl,            plotpath, plotlevel=20)
    do_plot(datapath, m_ci,            plotpath, plotlevel=20)
    do_plot(datapath, ls_prec,         plotpath)
    do_plot(datapath, w_in_wth,        plotpath, plotlevel=20)
    do_plot(datapath, u_in_w3,         plotpath, plotlevel=20)
    do_plot(datapath, v_in_w3,         plotpath, plotlevel=20)
