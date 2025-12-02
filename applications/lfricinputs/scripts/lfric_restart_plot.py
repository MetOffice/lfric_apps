#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file LICENCE.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

import os
import argparse
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from netCDF4 import Dataset

UM_RMDI = -32768.0*32768.0

# Set a non-interactive backend.
plt.switch_backend('Agg')

# Set default colormap to viridis
if hasattr(cm, "viridis"):
    plt.set_cmap(cm.viridis)
elif hasattr(cm, "winter"):
    plt.set_cmap(cm.winter)


def load_data_and_extract_field(datapath, field_name, level, time_idx=0):
    """
    Returns a single level's worth of data and corresponding longitude
    and latitude values.
    """

    # Load the checkpoint file
    if os.path.isfile(datapath):
        # Load as NetCDF
        dataset = Dataset(datapath)

        # Find the variable...
        field = None
        if field_name in dataset.variables:
            field = dataset.variables[field_name]
        else:
            # Maybe it's the long_name?
            for variable in dataset.variables.values():
                if (hasattr(variable, "long_name") and
                        variable.long_name == field_name):
                    field = variable
                    break

        if field is None:
            raise NameError("Field name {0} not found".format(field_name))

        # Get corresponding longitude and latitude arrays
        lat_name, lon_name = field.coordinates.split()
        longitudes = dataset.variables[lon_name]
        latitudes = dataset.variables[lat_name]

        # Get dimensions objects
        dims = field.get_dims()
        if "time" in dataset.variables:
            is_timeseries = True
            if len(dims) == 2:
                num_levels = 1
            else:
                num_levels= len(dims[1])
        else:
            is_timeseries = False
            if len(dims) == 1:
                num_levels = 1 # 2d field only have xy dims
            else:
                num_levels = len(dims[0]) #3d fields levels is first dim

        print("Number of levels is {0}".format(num_levels))

        # Sanity check
        if level > num_levels:
            raise ValueError(
                "Requested level is greater than number of levels of field")

        # Read in field data and clip to UM_RMDI to cater for fields containing
        # LFRic RMDI values
        if len(dims) == 1:
            level_data = field[:]
        else:
            if is_timeseries:
                level_data = field[time_idx, level, :]
            else:
                level_data = field[level, :]
        level_data[level_data < UM_RMDI] = UM_RMDI

        # Set lattitude and longitude and number of data points        
        level_lon = longitudes[:]
        level_lat = latitudes[:]
        num_points = level_data.shape[0]

        # Mask data using UM_RMDI value
        level_data = ma.masked_values(level_data, UM_RMDI, copy=False)

        # Copy masks to lat/lon arrays
        level_lon.mask = level_data.mask
        level_lat.mask = level_data.mask

        # Remove masked points from arrays
        level_data = level_data.compressed()
        level_lon = level_lon.compressed()
        level_lat = level_lat.compressed()
    else:
        raise ValueError("Can't find input file/directory: {0}"
                         .format(datapath))

    return level_lon, level_lat, level_data, num_points, field


def export_plot_data(file_name, lon, lat, data, num_points):
    """
    Used to export raw plotting data to file
    """
    f = open(file_name, mode='w')
    for (latitude, longitude, datum) in zip(lat, lon, data):
        f.write("{0},{1},{2}\n".format(latitude, longitude, datum))
    f.close()


def interpolate_to_regular_grid(level_lon, level_lat, level_data, num_points):
    """
    To smooth over the points near the poles when viewed in a more
    traditional lat-lon projection.
    """

    # Calculate the grid "C" value and a roughly equivalent "N" value...
    c_value = int((num_points/6)**0.5)
    n_value = np.round(c_value*(12.0/7.0)**0.5)

    # Round this up to a sensible(ish) looking N value
    # (just a multiple of 48 for now)
    n_value = int(np.ceil(n_value / 48.0) * 48)

    # Size (same as original dump here)
    nx, ny = 2*n_value, 1.5*n_value

    # Generate the regular grid
    xi = np.linspace(min(level_lon), max(level_lon), nx)
    yi = np.linspace(min(level_lat), max(level_lat), ny)
    xf, yf = np.meshgrid(xi, yi)

    # Interpolate using delaunay triangularization
    zi = griddata((level_lon, level_lat),
                  level_data, (xf, yf), method='linear')

    return xi, yi, zi, c_value, n_value


def gen_markersize(n_points):
    """
    The "points" based plot of the LFRic grid needs the size of the
    dots to be set sensibly or plotting artifacts are introduced.
    This method uses a rough exponential fit to pick a size.
    """
    c_value = int((n_points/6)**0.5)

    def exponential(x, a, b, c):
        return a*np.exp(-b*x)+c

    # C Value sample points
    x = [48, 108, 168, 192]
    # Markersize values
    y = [40, 12, 6, 1]

    # Create the fit
    pcoeffs, _ = curve_fit(exponential, x, y, p0=(1, 1e-6, 1))

    # Now derive the value (don't return a value less than 1)
    return max(round(exponential(c_value, *pcoeffs)), 1)


def make_single_figure(datapath, field_name, level, time_idx=0, interp=False,
                       plot_path=None):
    """
    Covers the 2 "single" figure plot types (either interpolated or just
    the point values).
    """

    # Get the data
    level_lon, level_lat, level_data, num_points, field = load_data_and_extract_field(
        datapath, field_name, level, time_idx)

    # Make the figure and axis
    fig = plt.figure(figsize=(20, 15))
    ax = fig.add_subplot(111)

    if interp:
        # Do the interpolation
        xi, yi, zi, c_value, n_value = interpolate_to_regular_grid(
            level_lon, level_lat, level_data)

        # Plot the interpolated data using contourf
        plot = ax.contourf(xi, yi, zi, 500)

        # For filename/title
        interp_fname = "_interp"
        interp_title = ", Interpolated to regular grid (N{0})".format(n_value)

    else:
        # Get a suitable size for the plot dots
        markersize = gen_markersize(num_points)
        plot = ax.scatter(
            level_lon, level_lat, c=level_data,
            edgecolor="none", s=markersize*2, alpha=0.8)

        # For filename/title
        interp_fname = "_points"
        interp_title = ""

    # Get units for title
    unit = ""
    if hasattr(field, "units") and field.units != "1":
        unit = " ({0})".format(field.units)

    # and name
    if hasattr(field, "long_name"):
        name = field.long_name.replace("_", " ").capitalize()
        file_name = field.long_name
    else:
        name = field_name
        file_name = field_name

    # Tidy up the axes
    ax.set_xlim([-180, 180])
    ax.set_ylim([-90, 90])
    ax.set_yticks(range(-90, 91, 30))
    ax.set_xticks(range(-180, 181, 60))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    fig.colorbar(plot, aspect=50, shrink=0.8, pad=0.01, format='%.4e')
    ax.set_title("{0}{1}, Layer: {2}{3}".format(
        name, unit, level, interp_title), fontsize=15)
    ax.patch.set_facecolor("#DDDDDD")
    ax.grid(color="white", linewidth=2)
    fig.subplots_adjust(right=0.999)
    # Save file as .png
    if plot_path is None:
        plot_path = "."
    out_file_name = os.path.join(plot_path,
                                 "{0}_{1}-{2}{3}.png"
                                 .format(file_name, level, time_idx, interp_fname))
    fig.savefig(out_file_name, bbox_inches='tight')

    data_file_name = os.path.join(plot_path,
                                 "{0}_{1}-{2}{3}.dat"
                                 .format(file_name, level, time_idx, interp_fname))
    export_plot_data(data_file_name, level_lon, level_lat, level_data, num_points)


def make_double_figure(datapath, field_name,  level, time_idx=0, plot_path=None):
    """
    Covers the "double" figure plot types.
    """

    # Get the data
    level_lon, level_lat, level_data, num_points, field = load_data_and_extract_field(
        datapath, field_name, level, time_idx)

    # Make the figure and axes
    fig = plt.figure(figsize=(24, 12))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax3 = fig.add_axes([0.2, 0.1, 0.6, 0.01])

    # Do the interpolation
    xi, yi, zi, c_value, n_value = interpolate_to_regular_grid(
        level_lon, level_lat, level_data, num_points)

    # Get the min/max color range to set it equal
    vmin = np.min([np.nanmin(zi), np.nanmin(level_data)])
    vmax = np.max([np.nanmax(zi), np.nanmax(level_data)])

    # Plot the interpolated data using contourf
    plot1 = ax2.contourf(xi, yi, zi, 500, vmin=vmin, vmax=vmax)

    # Get a suitable size for the plot dots
    markersize = gen_markersize(num_points)

    # Plot the non interpolated data
    plot2 = ax1.scatter(
        level_lon, level_lat, c=level_data,
        edgecolor="none", s=markersize, alpha=0.8, vmin=vmin, vmax=vmax)

    # Get units for title
    unit = ""
    if hasattr(field, "units") and field.units != "1":
        unit = " ({0})".format(field.units)

    # and name
    if hasattr(field, "long_name"):
        name = field.long_name.replace("_", " ").capitalize()
        file_name = field.long_name
    else:
        name = field_name
        file_name = field_name

    ax1.set_title("Point Values on LFRic grid (C{0})".format(c_value),
                  fontsize=15)

    ax2.set_title("LFRic points Interpolated onto regular grid (N{0})"
                  .format(n_value), fontsize=15)

    # Tidy up the axes - these properties are the same in both
    for ax in [ax1, ax2]:
        ax.set_xlim([-180, 180])
        ax.set_ylim([-90, 90])
        ax.set_yticks(range(-90, 91, 30))
        ax.set_xlabel("Longitude")
        ax.patch.set_facecolor("#DDDDDD")
        ax.grid(color="white", linewidth=2)

    # Plus some minor tweaks to account for the fact that the plots
    # are right next to each other
    ax1.set_xticks(range(-180, 180, 60))
    ax2.set_xticks(range(-180, 181, 60))
    ax1.set_ylabel("Latitude")
    ax2.set_yticklabels([])

    # Title and colorbar
    fig.suptitle("{0}{1}, Layer: {2}".format(
        name, unit, level), fontsize=20, y=0.95)
    fig.colorbar(plot2, cax=ax3, orientation="horizontal", )
    fig.subplots_adjust(wspace=0, bottom=0.175)

    # Save file as .png
    if plot_path is None:
        plot_path = "."
    out_file_name = os.path.join(plot_path,
                                 "{0}_{1}-{2}.png"
                                 .format(file_name, level, time_idx))
    fig.savefig(out_file_name, bbox_inches='tight')

    data_file_name = os.path.join(plot_path,
                                 "{0}_{1}-{2}{3}.dat"
                                 .format(file_name, level, time_idx, interp_fname))
    export_plot_data(data_file_name, level_lon, level_lat, level_data, num_points)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Checkpoint Restart Plotter for LFRic",
        )

    parser.add_argument("input_file",
                        help="path to an LFRic checkpoint restart file, "
                        "or a directory containing UM2LFRic output text files")
    parser.add_argument("field_name",
                        help="name of field to plot (can be the NetCDF "
                        "variable name OR the long_name) must match the name "
                        "in the filenames if using on UM2LFRic text files")
    parser.add_argument("level",
                        help="index of the level to be plotted",
                        type=int)
    parser.add_argument("time_idx",
                        help="index of the time to be plotted",
                        type=int)                  
    parser.add_argument("--plot_path",
                        help="destination directory for plot "
                        "(default is current dir)")
    parser.add_argument("--plot_type",
                        help="what sort of plot? 'points', "
                        "'interpolated', or 'both' (default is both)")

    args = parser.parse_args()

    plot_type = args.plot_type
    if plot_type is None:
        plot_type = "points"
    elif plot_type not in ("points", "interpolated", "both"):
        raise ValueError("Unknown plot type: {0}".format(plot_type))

    if plot_type == "both":
        make_double_figure(args.input_file,
                           args.field_name,
                           args.level,
                           time_idx=args.time_idx,
                           plot_path=args.plot_path)
    else:
        if plot_type == "interpolated":
            interpolate = True
        else:
            interpolate = False

        make_single_figure(args.input_file,
                           args.field_name,
                           args.level,
                           time_idx=args.time_idx,
                           plot_path=args.plot_path,
                           interp=interpolate)
