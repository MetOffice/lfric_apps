#!/usr/bin/env python

# Need to set a non-interactive backend for suites
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import sys
import iris

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from read_data import read_ugrid_data

import matplotlib.cm as cm

#iris.FUTURE.netcdf_promote = True

# Size of regular grid
ny, nx = 360, 720

def make_figures(filein, plotpath, fields, vertical_spacing, lid, n_full, figname):

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

    directions = ['yz']

    if fields is None:
      fields = ['theta', 'u_in_w2h']
    for field in fields:

        cube = read_ugrid_data(filein, field)
        time = np.around(cube.coord('time').points, decimals=1)

        levels_name = cube.dim_coords[-1].name()
        #Set some levels for contours:
        levels=None
        if field=='theta':
            levels = np.linspace(200, 800, 31)
        if field=='u1':
            levels = np.arange(-20,36,4.)
        if field == 'u2':
            levels = np.linspace(-1.0, 1.0, 11)
        if field == 'u3':
            levels = np.linspace(-0.003, 0.003, 11) 
        if field == 'exner':
            levels = np.linspace(916, 1020, 14) # exner will be converted to hPa
        if field == 'rho':
            levels = np.arange(0,1.4,0.05)    

        n_levs = len(cube.coord(levels_name).points)
        plot_data=np.zeros((ny,nx,n_levs))

        # Compute the horizontal grid
        x = np.around(cube.coord('longitude').points, decimals=5)
        y = np.around(cube.coord('latitude').points, decimals=5)

        xmin = np.amin(x)
        xmax = np.amax(x)
        ymin = np.amin(y)
        ymax = np.amax(y)

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
            data = cube.data[0,l]
            fi = griddata((x, y), data, (xf, yf), method='linear')
            fi_n = griddata((x, y), data, (xf, yf), method='nearest')
            fi[np.isnan(fi)] = fi_n[np.isnan(fi)] 

            if field == 'exner':
                # Convert to hPa
                rd = 287.05
                p0 = 100000.0
                kappa = rd/1005.0
                inc = 0.01*fi**(1.0/kappa) * p0
            else:
                inc = fi

                plot_data[:,:,l] = plot_data[:,:,l] + inc

        for direction in directions:
            interp_fig = plt.figure(figsize=(20,10))

            if field=='theta' and direction == 'xy':
                levels = np.linspace(200, 320, 13)

            nxplots=1
            nyplots=1
            nplots=nxplots*nyplots

            # Color map
            c_map = cm.summer

            for l in range(nplots):
                ax = interp_fig.add_subplot(nxplots,nyplots,l+1)

                if direction == 'xz':
                    mean_plot_data = np.mean(plot_data,axis=0)
                    lon, height = np.meshgrid(xi, zi)
                    CS=plt.contourf(lon, height, mean_plot_data.T, levels=levels, cmap=c_map)
                    cbar=plt.colorbar()
                    cbar.ax.tick_params(labelsize=28)
                    CL=plt.contour(lon, height, mean_plot_data.T, levels=levels, linewidths=0.5, colors='k')
                    plt.xlabel("Longitude", fontsize=28)
                    plt.ylabel("Height", fontsize=28)
                if direction == 'yz':
                    mean_plot_data = np.mean(plot_data,axis=1)
                    lat, height = np.meshgrid(yi, zi)
                    CS=plt.contourf(lat, height, mean_plot_data.T, levels=levels, cmap=c_map)
                    cbar=plt.colorbar()
                    cbar.ax.tick_params(labelsize=28)
                    CL=plt.contour(lat, height, mean_plot_data.T, levels=levels, linewidths=0.5, colors='k')
                    plt.xlabel("Latitude", fontsize=28)
                    plt.ylabel("Height", fontsize=28)
                if direction == 'xy':
                    mean_plot_data = plot_data[:,:,0]#np.mean(plot_data,axis=2)
                    lat, lon = np.meshgrid(yi, xi)
                    CS=plt.contourf(lon, lat, mean_plot_data.T, levels=levels, cmap=c_map)
                    cbar=plt.colorbar()
                    cbar.ax.tick_params(labelsize=28)
                    CL=plt.contour(lon, lat, mean_plot_data.T, levels=levels, linewidths=0.5, colors='k')
                    plt.clabel(CL, CL.levels[1::2], fontsize=15, inline=1, fmt='%3.1f')
                    plt.xlabel("Longitude", fontsize=28)
                    plt.ylabel("Latitude", fontsize=28)

            plt.tight_layout()

            pngfile='%s/%s-%s-time-mean-%s.png' % (plotpath, figname, field, direction)
            plt.savefig(pngfile)
            plt.close()


if __name__ == "__main__":

  try:
     args=sys.argv[:]
     filein, plotpath, figname, vertical_grid, lid, n_full = args[1:7]
     field_list=None
     if len(args[:])>7: field_list=args[6].split(':')
  except ValueError:
     print("Usage: {0} <filein> <plotpath> <figname> <vertical_grid> <lid> <n_full> [<fields_list>]".format(sys.argv[0]))
     exit(1)

  make_figures(filein, plotpath, field_list, vertical_grid, int(lid), int(n_full), figname)
