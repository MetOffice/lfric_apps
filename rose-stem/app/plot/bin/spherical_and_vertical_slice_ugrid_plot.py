#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
Plots horizontal slices and vertical cross sections from ugrid NetCDF data.
Can loop over fields, slices and time steps.

Plotting is designed for vertical slice transport tests but should be
compatible with planar and spherical domains.
"""
import sys
import math
import iris
import matplotlib
# Need to set a non-interactive backend for suites
matplotlib.use('Agg')  # noqa: E402
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata
from read_data import read_ugrid_data
from six.moves import range


# --------------------------------------------------------------------------- #
# Rounding functions for neatly determining contour levels
# --------------------------------------------------------------------------- #


def roundup(number, digits=0):
    n = 10**-digits
    return round(math.ceil(number / n) * n, digits)


def rounddown(number, digits=0):
    n = 10**-digits
    return round(math.floor(number / n) * n, digits)

# --------------------------------------------------------------------------- #
# Main plotting script
# --------------------------------------------------------------------------- #


def make_figures(filein, plotpath, field_list, slice_list,
                 time_list, extrusion, testname):

    # If specific model levels / slices are desired, they can be put here
    plot_lat = None
    plot_lon = None
    plot_level = 5 if testname == 'lam_gw' else 0

# --------------------------------------------------------------------------- #
# Set some plotting style details
# --------------------------------------------------------------------------- #

    plt.rc('font', family='serif', size=36)

    # Preset conditions for some specific tests
    # Hard wire in top lid values
    if testname in ['curl_free', 'cylinder', 'div_free',
                    'eternal_fountain', 'rotational', 'translational']:
        spherical = False

        zmin = 0.0
        zmax = 2000.

    elif testname == 'bryan_fritsch':
        spherical = False
        zmin = 0.0
        zmax = 10000.

    elif testname == 'grabowski_clark':
        spherical = False
        zmin = 0.0
        zmax = 2400.

    elif testname in ['baroclinic', 'aquaplanet', 'spherical', 'lam_gw',
                      'sbr', 'dcmip101', 'vert_def', 'hadley', 'cos_phi',
                      'four_part_sbr', 'deformation', 'divergent']:
        spherical = True

        if testname in ['spherical', 'cos_phi', 'four_part_sbr', 
                        'deformation', 'divergent']:
            # This is a special 2D spherical shell
            zmin = 0.0
            zmax = 1.0
        elif testname == 'lam_gw':
            zmin = 0.0
            zmax = 10000.0  # A 10 km lid
        elif testname in ['sbr', 'dcmip101', 'vert_def', 'hadley']:
            zmin = 0.0
            zmax = 12000.0  # A 12 km lid
        elif spherical:
            zmin = 0.0
            zmax = 30000.   # Assume 30 km lid

    else:
        # Assume a spherical set-up if known test not specified
        spherical = True

    plot_zmin = zmin
    plot_zmax = zmax

    # Some specific field labels
    title_dict = {'theta': r'$\theta \ / $ K',
                  'theta_e': r'$\theta_e \ / $ K',
                  'theta_e_moist': r'$\theta_e \ / $ K',
                  'theta_vd_pert': r"$\theta'_{vd} \ / $ K",
                  'rho': r'$\rho \ / $ kg m$^{-3}$',
                  'density': r'$\rho \ / $ kg m$^{-3}$',
                  'm_v': r'$m_v \ / $ kg kg$^{-1}$',
                  'm_cl': r'$m_{cl} \ / $ kg kg$^{-1}$',
                  'tracer': r'$m_v \ / $ kg kg$^{-1}$',
                  'tracer_con': r'$q_c \ / $ kg kg$^{-1}$',
                  'tracer_adv': r'$q_a \ / $ kg kg$^{-1}$',
                  'constant': r'$q_c \ / $ kg kg$^{-1}$',
                  'buoyancy': r'$b \ / $ m s$^{-2}$',
                  'w2_vector1': r'$F_1 \ / $ m s$^{-1}$',
                  'w2_vector2': r'$F_2 \ / $ m s$^{-1}$',
                  'w2_vector3': r'$F_3 \ / $ m s$^{-1}$',
                  'w3_aerosol': r'$a_\rho \ / $ kg kg$^{-1}$',
                  'wt_aerosol': r'$a_\theta \ / $ kg kg$^{-1}$'}

    # Find number of full levels by asking for theta
    try:
        cube = read_ugrid_data(filein, 'theta')
    except iris.exceptions.ConstraintMismatchError:
        cube = read_ugrid_data(filein, 'buoyancy')
    levels_name = cube.dim_coords[-1].name()
    nz_full = len(cube.coord(levels_name).points)

    # Set plot_level for certain tests
    if testname in ['sbr', 'dcmip101', 'vert_def']:
        plot_level = int(np.floor(nz_full/2))

    # Find max horizontal extent for non-spherical domains
    try:
        cube = read_ugrid_data(filein, 'u_in_w2h')
    except iris.exceptions.ConstraintMismatchError:
        # If u_in_w2h is not in cube, then try wind1
        # (e.g. for transport miniapp)
        cube = read_ugrid_data(filein, 'wind1')

    x_coord_name = "longitude"
    y_coord_name = "latitude"
    try:
        lon_data = np.around(cube.coord(x_coord_name).points, decimals=5)
        lat_data = np.around(cube.coord(y_coord_name).points, decimals=5)
    except iris.exceptions.CoordinateNotFoundError:
        x_coord_name = 'projection_x_coordinate'
        y_coord_name = 'projection_y_coordinate'
        lon_data = np.around(cube.coord(x_coord_name).points, decimals=5)
        lat_data = np.around(cube.coord(y_coord_name).points, decimals=5)
    except Exception as err:
        msg = "error extracting planar mesh: " + str(err)
        print(msg)
        exit(1)

    plot_lonmax = np.amax(np.abs(lon_data))
    plot_latmax = np.amax(np.abs(lat_data))

    epsilon = 1e-12
    if testname == 'lam_gw':
        # This is a domain not centred on zero
        plot_lonmin = np.amin(lon_data)
    # Assume a sensible domain here
    elif np.amin(lon_data) < - epsilon:
        # Assume this means lon=0 in the middle of the domain
        plot_lonmin = -plot_lonmax
    else:
        # Assume this means lon=0 is the min of the domain
        plot_lonmin = 0.0
    if np.amin(lat_data) < - epsilon:
        # Assume this means lat=0 in the middle of the domain
        plot_latmin = -plot_latmax
    else:
        # Assume this means lat=0 is the min of the domain
        plot_latmin = 0.0

# --------------------------------------------------------------------------- #
# Extract data
# --------------------------------------------------------------------------- #

    for field in field_list:

        # This is a diagnostic field so we extract several fields
        # Remove this if we actually have theta_e as a diagnostic
        if field in ['theta_e', 'theta_e_moist']:
            cube = get_theta_e_cube(filein, zmin, zmax, field)
        elif field == 'theta_vd_pert':
            cube = get_theta_vd_pert_cube(filein)
        else:
            cube = read_ugrid_data(filein, field)

        # Vertical levels will be last entry in dimension coords
        levels_name = cube.dim_coords[-1].name()

        try:
            time = np.around(cube.coord('time').points, decimals=1)
        except iris.exceptions.CoordinateNotFoundError:
            # there is probably only one time level
            time = [0.0]

# --------------------------------------------------------------------------- #
# Create mesh-grid for plotting on
# --------------------------------------------------------------------------- #

        # Compute the horizontal grid
        x_coord_name = "longitude"
        y_coord_name = "latitude"
        if spherical:
            try:
                lon_data = np.around(
                    cube.coord(x_coord_name).points, decimals=5)
                lat_data = np.around(
                    cube.coord(y_coord_name).points, decimals=5)
            except iris.exceptions.CoordinateNotFoundError:
                x_coord_name = 'projection_x_coordinate'
                y_coord_name = 'projection_y_coordinate'
                lon_data = np.around(
                    cube.coord(x_coord_name).points, decimals=5)
                lat_data = np.around(
                    cube.coord(y_coord_name).points, decimals=5)
            except Exception as err:
                msg = "error extracting planar mesh: " + str(err)
                print(msg)
                exit(1)

            # Make a latitude-longitude grid with spacing of 0.5 degrees
            nlat, nlon = 360, 720
            if plot_lon is None:
                plot_lon = 360
            if plot_lat is None:
                plot_lat = 180
            if plot_level is None:
                plot_level = 0

            length_scaling = 1.0
            height_scaling = 1.0 / 1000.0  # Height in km

            if testname != "lam_gw":
                plot_latmin = -90.
                plot_latmax = 90.
                plot_lonmin = -180.
                plot_lonmax = 180.

        else:
            try:
                lon_data = np.around(
                    cube.coord(x_coord_name).points, decimals=5)
                lat_data = np.around(
                    cube.coord(y_coord_name).points, decimals=5)
            except iris.exceptions.CoordinateNotFoundError:
                x_coord_name = 'projection_x_coordinate'
                y_coord_name = 'projection_y_coordinate'
                lon_data = np.around(
                    cube.coord(x_coord_name).points, decimals=5)
                lat_data = np.around(
                    cube.coord(y_coord_name).points, decimals=5)
            except Exception as err:
                msg = "error extracting planar mesh: " + str(err)
                print(msg)
                exit(1)

            nlat, nlon = len(np.unique(lat_data)), len(np.unique(lon_data))
            if plot_lon is None:
                plot_lon = 0
            if plot_lat is None:
                plot_lat = 0
            if plot_level is None:
                plot_level = 0

            length_scaling = 1.0 / 1000.0  # Lengths in km
            height_scaling = 1.0 / 1000.0  # Height in km

        nz = len(cube.coord(levels_name).points)

        lonmin = np.amin(lon_data)
        lonmax = np.amax(lon_data)
        latmin = np.amin(lat_data)
        latmax = np.amax(lat_data)

        # Generate a regular grid to interpolate the data.
        lon1d = np.linspace(lonmin, lonmax, nlon)
        lat1d = np.linspace(latmin, latmax, nlat)

        lat2d, lon2d = np.meshgrid(lat1d, lon1d)

        z1d_full = make_extrusion(extrusion, nz_full, zmin, zmax)
        z1d_half = 0.5*(z1d_full[1:] + z1d_full[0:nz_full-1])

        # Choose whether to use heights from full level set or half levels
        if nz == nz_full:
            z1d = z1d_full
        else:
            z1d = z1d_half

# --------------------------------------------------------------------------- #
# Interpolate or read in data
# --------------------------------------------------------------------------- #

        # Shortcut for plotting all dumped timesteps if 'all' is in time list
        if 'all' in time_list:
            time_list = range(len(time))

        for time_point in time_list:

            t = int(time_point)

            plot_data = np.zeros((nlon, nlat, nz))

            # Interpolate using delaunay triangularization
            for level in range(nz):
                if len(time) > 1:
                    data = cube.data[t, level]
                else:
                    data = cube.data[level]  # There is no time coordinate
                fi = griddata((lon_data, lat_data), data,
                              (lon2d, lat2d), method='linear')
                fi_n = griddata((lon_data, lat_data), data,
                                (lon2d, lat2d), method='nearest')
                fi[np.isnan(fi)] = fi_n[np.isnan(fi)]

                plot_data[:, :, level] = fi

# --------------------------------------------------------------------------- #
# Loop through data and pull out meshgrid and data for each slice
# --------------------------------------------------------------------------- #

            for slice in slice_list:
                slice_fig = plt.figure(figsize=(15, 10))

                if slice == 'xy':
                    Y_meshgrid, X_meshgrid = \
                        np.meshgrid(lat1d*length_scaling, lon1d*length_scaling)
                    slice_data = plot_data[:, :, plot_level]
                    slice_label = r'$z$: %1.1e km,' % (z1d[plot_level] / 1000.)
                    if spherical:
                        plot_xlabel = r'$\lambda \ / $ deg'
                        plot_ylabel = r'$\phi \ / $ deg'
                    else:
                        plot_xlabel = r'$x \ / $ km'
                        plot_ylabel = r'$y \ / $ km'
                    plot_xmin = plot_lonmin * length_scaling
                    plot_xmax = plot_lonmax * length_scaling
                    plot_ymin = plot_latmin * length_scaling
                    plot_ymax = plot_latmax * length_scaling

                elif slice == 'xz':
                    if nz < 2:
                        message = 'Cannot plot xz slice for field %s' % field
                        message += 'which is only 1 layer thick'
                        raise ValueError(message)
                    Y_meshgrid, X_meshgrid = \
                        np.meshgrid(z1d*height_scaling, lon1d*length_scaling)
                    slice_data = plot_data[:, plot_lat, :]
                    if spherical:
                        slice_label = r'$\phi$: %1.1e deg,' \
                            % (lat1d[plot_lat] * length_scaling)
                        plot_xlabel = r'$\lambda \ / $ deg'
                    else:
                        slice_label = r'$y$: %1.1e km,' % \
                            (lat1d[plot_lat] * length_scaling)
                        plot_xlabel = r'$x \ / $ km'
                    plot_ylabel = r'$z \ / $ km'
                    plot_xmin = plot_lonmin * length_scaling
                    plot_xmax = plot_lonmax * length_scaling
                    plot_ymin = plot_zmin * height_scaling
                    plot_ymax = plot_zmax * height_scaling

                elif slice == 'yz':
                    if nz < 2:
                        message = 'Cannot plot yz slice for field %s' % field
                        message += ' which is only 1 layer thick'
                        raise ValueError(message)
                    Y_meshgrid, X_meshgrid = \
                        np.meshgrid(z1d*height_scaling, lat1d*length_scaling)
                    slice_data = plot_data[plot_lon, :, :]
                    if spherical:
                        slice_label = r'$\lambda$: %1.1e deg,' % \
                            (lon1d[plot_lon] * length_scaling)
                        plot_xlabel = r'$\phi \ / $ deg'
                    else:
                        slice_label = r'$x$: %1.1e km,' % \
                            (lon1d[plot_lon] * length_scaling)
                        plot_xlabel = r'$y \ / $ km'
                    plot_ylabel = r'$z \ / $ km'
                    plot_xmin = plot_latmin * length_scaling
                    plot_xmax = plot_latmax * length_scaling
                    plot_ymin = plot_zmin * height_scaling
                    plot_ymax = plot_zmax * height_scaling

                else:
                    raise ValueError('Cannot plot slice %s' % slice)

# --------------------------------------------------------------------------- #
# Set up contour details
# --------------------------------------------------------------------------- #

                # Special contours for our known tests
                if ((testname in ['cylinder', 'div_free',
                                  'eternal_fountain', 'rotational',
                                  'translational', 'curl_free'] and
                        field in ['theta', 'density', 'rho', 'm_v',
                                  'tracer', 'tracer_con', 'tracer_adv',
                                  'w2_vector1', 'w2_vector2',
                                  'wt_aerosol', 'w3_aerosol'])
                    or (testname in ['sbr', 'dcmip101', 'vert_def'] and
                        field in ['theta', 'm_v', 'tracer', 'tracer_con', 'tracer_adv'])):

                    # Hardwire contour details
                    # 2.0 is the background value (usually the minimum)
                    # 5.0 is the expected maximum value
                    # Take steps of 0.5 from 0.0 to 6.0
                    # but omit 2.0 for contour lines
                    step = 0.5
                    tracer_background = 2.0
                    tracer_max = 5.0
                    max_field = tracer_max + 2*step
                    min_field = tracer_background - 4*step
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = []
                    epsilon = 1e-14
                    for contour in contour_colours:
                        if abs(contour - tracer_background) > epsilon:
                            contour_lines.append(contour)

                elif (testname in ['cylinder', 'div_free',
                                    'eternal_fountain', 'rotational',
                                    'translational'] and
                        field == 'w2_vector3'):

                    # Hardwire contour details
                    # 0.0 is the background value (usually the minimum)
                    # 3.0 is the expected maximum value
                    # Take steps of 0.4 from -1.0 to 4.0
                    # but omit 2.0 for contour lines
                    step = 0.5
                    tracer_background = 0.0
                    tracer_max = 3.0
                    max_field = tracer_max + 2*step
                    min_field = tracer_background - 2*step
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = []
                    epsilon = 1e-14
                    for contour in contour_colours:
                        if abs(contour - tracer_background) > epsilon:
                            contour_lines.append(contour)

                elif (testname in ['hadley', 'cos_phi', 'deformation','divergent']
                      or (testname in ['sbr', 'dcmip101', 'vert_def'] and
                          field in ['density', 'rho'])):
                    step = 0.1
                    tracer_background = 0.0
                    tracer_max = 1.0
                    max_field = tracer_max + 2*step
                    min_field = tracer_background - 2*step
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = []
                    epsilon = 1e-14
                    for contour in contour_colours:
                        if abs(contour - tracer_background) > epsilon:
                            contour_lines.append(contour)
                elif (field == 'constant'):
                    step = 0.05
                    min_field = 0.8
                    max_field = 1.2
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = np.copy(contour_colours)
                elif (testname == 'four_part_sbr'):
                    step = 0.25
                    min_field = -3.0
                    max_field = 3.0
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = np.copy(contour_colours)

                elif (testname == 'lam_gw' and field == 'buoyancy'):
                    step = 0.002
                    min_field = -0.014
                    max_field = 0.014
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = np.copy(contour_colours)
                elif (testname == 'bryan_fritsch' and field in ['theta_e', 'theta_e_moist']):
                    step = 0.5
                    min_field = 318.
                    max_field = 330.
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = np.copy(contour_colours)
                elif (testname == 'bryan_fritsch' and field == 'theta'):
                    step = 0.5
                    min_field = 299.
                    max_field = 305.
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = np.copy(contour_colours)
                elif (testname == 'grabowski_clark' and field == 'm_cl'):
                    step = 2e-4
                    min_field = -4e-4
                    max_field = 2.4e-3
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = np.copy(contour_colours)
                elif (testname == 'grabowski_clark'
                      and field == 'theta_vd_pert'):
                    step = 0.5
                    min_field = -6.
                    max_field = 6.
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = np.copy(contour_colours)
                else:
                    # Set some nice contours based on the min and max fields
                    # Values are specific for this slice at this time
                    if np.amax(plot_data) != 0.0 and np.amin(plot_data) != 0.0:
                        digits = math.floor(-np.log10(np.amax(plot_data)
                                                      - np.amin(plot_data)))
                    else:
                        digits = 1
                    max_field = roundup(np.amax(plot_data), digits=digits)
                    min_field = rounddown(np.amin(plot_data), digits=digits)
                    if (max_field - min_field) < 1e-15:
                        max_field += 1e-14
                        min_field -= 1e-14
                    step = (max_field - min_field) / 20
                    contour_colours = np.arange(min_field, max_field+step,
                                                step=step)
                    contour_lines = np.copy(contour_colours)

# --------------------------------------------------------------------------- #
# Set up colour bar details
# --------------------------------------------------------------------------- #

                # This looks good also when converted to grey scale
                colour_code = 'Blues'
                # We scale colour map to avoid extreme colours
                colour_levels_scaling = 1.2

                # Make a colour map for more contours than we'll use
                actual_num_colour_levels = len(contour_colours)
                pure_num_colour_levels = \
                    np.ceil(actual_num_colour_levels*colour_levels_scaling)

                pure_cmap = cm.get_cmap(colour_code, pure_num_colour_levels)
                new_colours = pure_cmap(np.linspace(0,
                                                    1/colour_levels_scaling),
                                        actual_num_colour_levels)
                cmap = ListedColormap(new_colours)

# --------------------------------------------------------------------------- #
# Actually plot fields and contours
# --------------------------------------------------------------------------- #

                # Plot field with colours
                cf = plt.contourf(X_meshgrid, Y_meshgrid, slice_data,
                                  contour_colours, cmap=cmap, extend='both')

                # Set colours for over and under shoots
                cf.cmap.set_under('magenta')
                cf.cmap.set_over('yellow')

                # Plot contour lines
                cl = plt.contour(X_meshgrid, Y_meshgrid, slice_data,
                                 contour_lines, linewidths=1.0, colors='k',
                                 linestyle="", extend='min')

# --------------------------------------------------------------------------- #
# Plot labels
# --------------------------------------------------------------------------- #

                plt.xlim([plot_xmin, plot_xmax])
                plt.ylim([plot_ymin, plot_ymax])
                plt.xticks([plot_xmin, plot_xmax])
                plt.yticks([plot_ymin, plot_ymax])
                plt.xlabel(plot_xlabel, labelpad=-20)
                plt.ylabel(plot_ylabel, labelpad=-20)
                plt.title(slice_label+' max: %2.2e, min: %2.2e, time: %d s'
                          % (np.max(slice_data), np.min(slice_data), time[t]),
                          pad=40, size=30, loc='left')

                try:
                    colour_label = title_dict[field]
                except KeyError:
                    colour_label = field
                    # make strings latex friendly
                    colour_label = colour_label.replace('_', r'\_')
                plt.colorbar(cf, cmap=cmap, label=colour_label)

                # Save file
                pngfile = '%s/%s-%s-time%s-%s.png' % \
                    (plotpath, testname, field, time[t], slice)
                plt.savefig(pngfile, bbox_inches='tight')
                plt.close()

# --------------------------------------------------------------------------- #
# Set up heights with options for different extrusion types
# --------------------------------------------------------------------------- #


def make_extrusion(extrusion, nz, zmin, zmax):

    # Get level heights based on what the stretching is
    if extrusion == 'linear':
        z1d = np.linspace(zmin, zmax, nz)
    elif extrusion == 'quadratic':
        z1d = zmin + (zmax - zmin) * np.array([(float(k) / float(nz)) ** 2
                                               for k in range(nz)])
    elif extrusion == 'geometric':
        stretching = 1.03
        z1d = np.zeros([nz])
        for k in range(nz):
            z1d[k] = zmin + (zmax - zmin) * ((stretching**k - 1)
                                             / (stretching**nz - 1))

    elif extrusion == 'dcmip':
        mu = 15.
        z1d = np.zeros([nz])
        for k in range(nz):
            z1d[k] = zmin + ((zmax - zmin) / (np.sqrt(mu+1) - 1)
                             * (np.sqrt(mu * (float(k)/float(nz))**2 + 1) - 1))
    elif extrusion == 'um':
        if nz != 39:
            raise ValueError('For UM extrusion need 39 levels not %d' % nz)

        z1d = zmin + (zmax - zmin) * \
            np.array([.0000000, .0005095, .0020380, .0045854, .0081519,
                      .0127373, .0183417, .0249651, .0326074, .0412688,
                      .0509491, .0616485, .0733668, .0861040, .0998603,
                      .1146356, .1304298, .1472430, .1650752, .1839264,
                      .2037966, .2246857, .2465938, .2695209, .2934670,
                      .3184321, .3444162, .3714396, .3998142, .4298913,
                      .4620737, .4968308, .5347160, .5763897, .6230643,
                      .6772068, .7443435, .8383348, 1.000000])

    else:
        raise NotImplementedError(
            'Extrusion %s is not implemented' % extrusion)

    return z1d

# --------------------------------------------------------------------------- #
# Function for extracting theta_e cube
# --------------------------------------------------------------------------- #

def get_theta_e_cube(filein, zmin, zmax, field_name):
    """
    Takes name of file and returns cube for theta_e variable. If the field name
    is 'theta_e_moist', then the heat capacities of water species are accounted
    for rather than zeroed.
    """

    p_zero = 100000.0  # Reference pressure in Pa
    Rd = 287.05
    Rv = 461.51
    cpd = 1005.0
    cpv = 1850.0 if field_name == 'theta_e_moist' else 0.0
    cl = 4180.0 if field_name == 'theta_e_moist' else 0.0
    Lv0 = 2.501e6
    T_ref = 273.15

    cube_theta = read_ugrid_data(filein, 'theta')
    cube_mr_v = read_ugrid_data(filein, 'm_v')
    cube_mr_cl = read_ugrid_data(filein, 'm_cl')
    cube_exner = read_ugrid_data(filein, 'exner')
    cube_exner_wtheta = cube_theta.copy()
    cube_shape = np.shape(cube_exner.data)

    # Get height data
    levels_name = cube_theta.dim_coords[-1].name()
    nz_full = len(cube_theta.coord(levels_name).points)
    levels_name = cube_exner.dim_coords[-1].name()
    z1d_full = make_extrusion(extrusion, nz_full, zmin, zmax)
    z1d_half = 0.5*(z1d_full[1:] + z1d_full[0:nz_full-1])

    for j in range(0,cube_shape[1]):   # Loop through columns
        # Map from W3 to Wtheta
        weight_denom = z1d_half[j] - z1d_half[j-1]
        weight_upper = z1d_full[j] - z1d_half[j-1]
        weight_lower = z1d_half[j] - z1d_full[j]
        cube_exner_wtheta.data[:, j, :] = (
            weight_upper / weight_denom *cube_exner.data[:, j, :]
            + weight_lower / weight_denom *cube_exner.data[:, j, :]
        )

    # Bottom
    weight_denom = z1d_half[1] - z1d_half[0]
    weight_upper = z1d_full[0] - z1d_half[0]
    weight_lower = z1d_half[1] - z1d_full[0]
    cube_exner_wtheta.data[:, 0, :] = (
        weight_upper / weight_denom *cube_exner.data[:, 1, :]
        + weight_lower / weight_denom *cube_exner.data[:, 0, :]
    )
    # Top
    weight_denom = z1d_half[cube_shape[1]-1] - z1d_half[cube_shape[1]-2]
    weight_upper = z1d_full[cube_shape[1]] - z1d_half[cube_shape[1]-2]
    weight_lower = z1d_half[cube_shape[1]-1] - z1d_full[cube_shape[1]]
    cube_exner_wtheta.data[:, cube_shape[1], :] = (
        weight_upper / weight_denom * cube_exner.data[:, cube_shape[1]-1, :]
        + weight_lower / weight_denom *cube_exner.data[:, cube_shape[1]-2, :]
    )

    # Can't take powers of cubes, so manipulate data directly
    cube_exner_wtheta_pow = cube_exner_wtheta.copy()
    cube_exner_wtheta_pow.data = cube_exner_wtheta.data ** (cpd / Rd)

    # Theta_e
    cube_T = cube_theta * cube_exner_wtheta
    cube_p = p_zero * cube_exner_wtheta_pow / (1.0 + Rv / Rd * cube_mr_v)
    cube_Lv = cube_T.copy()
    cube_Lv.data = Lv0 - (cl - cpv)*(cube_T.data - T_ref)
    kappa = Rd / (cpd + cl * (cube_mr_v.data + cube_mr_cl.data))

    # Can't take powers of cubes, so manipulate data directly
    cube_p_pow = cube_p.copy()
    cube_p_pow.data = (cube_p.data / p_zero) ** (-kappa)

    exp_arg = cube_Lv * cube_mr_v / ((cpd + cl * (cube_mr_v + cube_mr_cl)) * cube_T)
    cube = cube_T * cube_p_pow * np.exp(exp_arg.data)

    return cube

# --------------------------------------------------------------------------- #
# Function for extracting theta_vd cube
# --------------------------------------------------------------------------- #


def get_theta_vd_pert_cube(filein):
    """
    Takes name of file and returns cube for perturbed theta_vd variable.
    """
    epsilon = 287. / 461.

    cube_theta = read_ugrid_data(filein, 'theta')
    cube_mr_v = read_ugrid_data(filein, 'm_v')
    cube_theta_vd = cube_theta * (1.0 + cube_mr_v / epsilon)

    # Adjust cube so that it is a perturbation from the background state
    cube_background = cube_theta_vd.copy()
    cube_shape = np.shape(cube_background.data)
    for i in range(cube_shape[0]):   # Loop through time steps
        for j in range(cube_shape[2]):   # Loop through columns
            # Get background state from first column on first time step
            cube_background.data[i, :, j] = cube_theta_vd.data[0, :, 0]

    cube = cube_theta_vd - cube_background

    return cube

# --------------------------------------------------------------------------- #
# Main function for passing fields to plot
# --------------------------------------------------------------------------- #


if __name__ == "__main__":

    try:
        args = sys.argv[:]
        filein, plotpath, fields, slices, \
            times, extrusion, testname = args[1:8]

        # Split out lists of fields, slices and times
        field_list = fields.split(':')
        slice_list = slices.split(':')
        # Should be indices of dumped times or 'all'
        time_list = times.split(':')

    except ValueError:
        print("Usage: {0} <filein> <plotpath> <field_list> \
              <slice_list> <time_list> <extrusion> <testname>"
              .format(sys.argv[0]))
        exit(1)

    make_figures(filein, plotpath, field_list, slice_list,
                 time_list, extrusion, testname)
