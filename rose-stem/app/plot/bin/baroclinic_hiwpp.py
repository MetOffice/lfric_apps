#!/usr/bin/env python3

# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import sys
import iris

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

from matplotlib import cm

# Use magma colormap
from magma import magma_data
from matplotlib.colors import ListedColormap
#magma = ListedColormap(magma_data, name='magma')
#plt.register_cmap(name='magma', cmap=magma)

#iris.FUTURE.netcdf_promote = True

# Size of regular grid
ny, nx = 768, 1536

plot_lon = 50#360
plot_lat = 125#280
plot_level = 1

def read_ugrid_data(filename, field):
    variable_constraint = iris.Constraint(cube_func=(
                                lambda c: c.var_name == field))
    return iris.load_cube(filename, constraint=variable_constraint)


def make_figures(filein, plotpath, fields, figname, vertical_spacing, direction):

    # Compute the vertical grid
    lid = 30.
    n_full = 31

    if vertical_spacing=='um':
        # um L38 set
        zi_f = np.array([.0, .0005095,  .0020380,  .0045854,  .0081519,  .0127373, 
                         .0183417,  .0249651,  .0326074,  .0412688,  .0509491,
                         .0616485,  .0733668,  .0861040,  .0998603,  .1146356, 
                         .1304298,  .1472430,  .1650752,  .1839264,  .2037966, 
                         .2246857,  .2465938,  .2695209,  .2934670,  .3184321, 
                         .3444162,  .3714396,  .3998142,  .4298913,  .4620737, 
                         .4968308,  .5347160,  .5763897,  .6230643,  .6772068, 
                         .7443435,  .8383348, 1.000000])*lid

    elif vertical_spacing=='dcmip':
        # dcmip Stretched grid
        mu = 15.
        zi_f = np.zeros([n_full])
        for k in range(n_full):
            zi_f[k] = lid * (np.sqrt(mu*(float(k)/float(n_full-1))**2 + 1.) - 1.)/(np.sqrt(mu+1.) - 1)
    else:
        # assume uniform grid
        zmin = 0.0
        zmax = lid
        zi_f = np.linspace(zmin, zmax, n_full)

  
    zi_h = 0.5*(zi_f[1:] + zi_f[0:n_full-1])

    #direction = 'xz'#, 'yz', 'xz'

    #t = -1
    if fields is None:
        fields = ['u_in_w2h', 'v_in_w2h', 'w_in_wth']

    for field in fields:
        interp_fig = plt.figure(figsize=(20,10))
        #cube = iris.load_cube(filein, field)
        cube = read_ugrid_data(filein, field)
        levels_name = cube.dim_coords[-1].name()
        #Set some levels for contours:
        levels=None
        if field=='air_potential_temperature':
            levels = np.linspace(220, 330, 12)
        if field=='eastward_wind' or field=='u1':
            levels = np.arange(-20,36,4.)
        if field == 'northward_wind' or field=='u2':
            levels = np.linspace(-25, 40, 14)
        if field == 'upward_air_velocity' or field=='w_in_wth':
            levels = 5*np.linspace(-0.0000011, 0.0000011, 12) 
        if field == 'exner_pressure':
            levels = np.linspace(916, 1020, 14) # exner will be converted to hPa
        if field == 'air_density':
            levels = np.arange(0,1.4,0.05)    
        if field == 'divergence_of_wind':
            levels = np.linspace(-0.00000011, 0.000000011, 12) 

        n_levs = len(cube.coord(levels_name).points)

        plot_data=np.zeros((ny,nx,n_levs))

        time = np.around(cube.coord('time').points, decimals=1)
        print('times = ',time)
        t = np.where(time == 86400.)[0][0]
        print( 't = ',t)
        # Compute the horizontal grid
        x = np.around(cube.coord('longitude').points, decimals=5)
        y = np.around(cube.coord('latitude').points, decimals=5)

        xmin = np.amin(x)
        xmax = np.amax(x)
        ymin = np.amin(y)
        ymax = np.amax(y)

            
        xmin = -180#130.0
        xmax = 180#140.0
        ymin = -90#-47.5
        ymax = 90#-22.5

        # Generate a regular grid to interpolate the data.
        xi = np.linspace(xmin, xmax, nx)
        yi = np.linspace(ymin, ymax, ny)

        xf, yf = np.meshgrid(xi, yi)

        # Choose the correct vertical level set
        if n_full == n_levs:
            zi = zi_f
        else:
            zi = zi_h

        # Interpolate using delaunay triangularization 
        for p,l in enumerate(range(n_levs)):
            data = cube.data[t,l]
            #data = cube.data[l]
            fi = griddata((x, y), data, (xf, yf), method='nearest')
            #fi_n = griddata((x, y), data, (xf, yf), method='nearest')
            #fi[np.isnan(fi)] = fi_n[np.isnan(fi)] 
         
            plot_data[:,:,l]=fi

            if field == 'exner_pressure':
                # Convert to hPa
                rd = 287.05
                p0 = 100000.0
                kappa = rd/1005.0
                plot_data[:,:,l] = 0.01*fi**(1.0/kappa) * p0

        ax = interp_fig.add_subplot(1,1,1)


        lat, lon = np.meshgrid(yi, xi)
        if field == 'exner_pressure' and plot_level == 0:
            # Extrapolate data to the surface
            dz = plot_data[:,:,0] + (zi_f[0] - zi_h[0])*(plot_data[:,:,0] - plot_data[:,:,1])/(zi_h[0] - zi_h[1])
        else:
            dz = plot_data[:,:,plot_level]

        if field == 'upward_air_velocity' or field == 'w_in_wth':
            # Remove zonal mean
            for j in range(ny):
                zonal_mean = np.mean(dz[j,:])
                dz[j,:] = dz[j,:] - zonal_mean

        #levels = np.linspace(-0.075,0.075,21)
        CS=plt.contourf(lon, lat, dz.T, levels=levels, cmap=cm.bwr, extend='both')
        cbar=plt.colorbar(cmap=cm.bwr, format="%.2e")
        cbar.ax.tick_params(labelsize=22)
        #CL=plt.contour(lon, lat, dz.T, levels=levels, linewidths=1.0, colors='k')
        #plt.clabel(CL, CL.levels[1::2], fontsize=15, inline=1, fmt='%3.1f')
        print('max/min = ',np.max(dz),np.min(dz))
        #print('long,lat,W at corner = ',lon[576,235],lat[576,235],dz[576,235])
        #plt.axis([130, 140, -47.5,-22.5])
        plt.xlabel("Longitude", fontsize=28)
        plt.ylabel("Latitude", fontsize=28)
        ax.tick_params(labelsize=24)
        pngfile='%s/%s-Zoom-%s-time%s-%s.png' % (plotpath, figname, field, time[t], direction)
        plt.savefig(pngfile)
        plt.close()


if __name__ == "__main__":

  try:
     args=sys.argv[:]
     filein, plotpath, figname, vertical_grid, direction, indices = args[1:7]
     field_list=None
     if len(args[:])>7: field_list=args[7].split(':')
  except ValueError:
     print("Usage: {0} <filein> <plotpath> <figname> <vertical_grid> <direction> [<fields_list>]".format(sys.argv[0]))
     exit(1)
 
  index_list = indices.split(':')
  make_figures(filein, plotpath, field_list, figname, vertical_grid, direction)
