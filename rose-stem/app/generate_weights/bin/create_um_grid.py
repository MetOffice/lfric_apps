#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file LICENCE
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

# Modified: O. Brunt, Met Office

'''
This module generates a UM grid file in SCRIP format from a namelist file.
The process can be broken down into steps, like so:

    1. Instantiate a GRID object.
    2. Set GRID instance attributes by reading from namelist and variable
       resolution mesh files.
       - Namelist file supplies rotation, poles, delta and lambda targets and
         origins, and x and y origins.
       - Variable resolution mesh file supplies row- and column-dependent
         constants lamba_p, lambda_u, phi_p, and phi_v.
       - If the variable resolution mesh file can be parsed and sets all
         the aforementioned GRID atributes to a value other than None,
         `GRID.variable_grid` is set to True.
    3. Create an `iris.cube.Cube` instance from the GRID instance.
       - If GRID.variable_grid == True, x and y coordinates are set from the
         constants provided by the variable mesh file.
       - If GRID.variable_grid is False, x and y coordinates are generated
         from x and y origins at an interval defined by delta_x and delta_y.
         This creates a fixed resolution mesh.
    4. If the coordinates are rotated, they are unrotated.
    5. Coordinates of the corners for each cell are extracted and unrotated.
    6. The area of each grid cell is calculated.
    7. The grid shape is defined from the latitude of cell centres.
    8. The grid is written to a NetCDF file at GRID_PATH_UM.

    TODO: Refactor to provide sufficient exception handling.

'''
# pylint: disable=import-error
import os
from netCDF4 import Dataset
from um_utils.cutout import CoordRotator
import numpy
import iris
import iris.fileformats
import iris.analysis
import iris.coord_systems
import f90nml


class GRID:
    '''
    A UM-specific GRID class. Most attributes held by this class are read from
    a namelist file or generated from the namelist information.

    Row- and column-dependent constants are read from a file with the
    `accessa_grid` extension.

    '''
    def __init__(self):
        self.vname = None     # Title to use in netCDF file
        self.grid = None      # UM grid type: P/U/V points at grid cell centres
        self.l_area = None    # If not 'no', output grid cell surface areas
        self.dlon = None      # Longitude of grid cell centres
        self.dlat = None      # Lattitude of grid cell centres
        self.dclo = None      # Longitude of grid cell corners
        self.dcla = None      # Lattitude of cell corners
        self.darea = None     # Area of each grid cell
        self.ucoor = None     # Coordinate units
        self.uarea = None     # Area units of each grid cell
        self.shape = None     # Shape of grid

        # Information read from namelist
        self.rotated_grid = None  # Flag to indicated grid is rotated
        self.pole_lon = None  # Lambda pole
        self.pole_lat = None  # Phi pole
        self.delx = None  # Delta lambda
        self.dely = None  # Delta phi
        self.npx = None  # Points lambda
        self.npy = None  # Points phi
        self.xorigin = None  # lambda origin
        self.yorigin = None  # Phi origin

        # Information read from accessa_grid file
        self.variable_grid = False
        self.lambda_p = None
        self.lambda_u = None
        self.phi_p = None
        self.phi_v = None

    def read_namelist(self, grid_namelist):
        '''
        Reads in a namelist file by converting to a dictionary using f90nml.
        Class attributes are then set using the lookup dictionary defined here.

        :param str grid_namelist: A path to a .nml file containing data under
                                  the headers defined in the setatttr_dict.

        '''
        # f90nml docs: https://f90nml.readthedocs.io/en/latest/
        # f90nml creates a dictionary from the namelist file
        grid_def = f90nml.read(grid_namelist)

        print(grid_def['grid'])

        setattr_dict = {'rotated_grid': 'rotated',
                        'pole_lon':     'lambda_pole',
                        'pole_lat':     'phi_pole',
                        'delx':         'delta_lambda_targ',
                        'dely':         'delta_phi_targ',
                        'npx':          'points_lambda_targ',
                        'npy':          'points_phi_targ',
                        'xorigin':      'lambda_origin_targ',
                        'yorigin':      'phi_origin_targ'}

        for key, value in setattr_dict.items():
            try:
                setattr(self, key, grid_def['grid'][value])
            except KeyError:
                pass

    def read_accessa_grid(self, accessa_grid):
        '''
        Reads in a variable.accessa_grid file by converting to a
        dictionary using f90nml.
        Class attributes are then set using the lookup dictionary
        defined here.

        :param str accessa_grid: A path to a .accessa_grid file or a
                                 file containing lambda_p, lambda_u,
                                 phi_p, and/or phi_v coordinate data
                                 under the headers defined in the
                                 setatttr_dict.

        '''
        # Convert input file to a dictionary for parsing
        grid_def = f90nml.read(accessa_grid)

        # Lookup table for setattr.
        setattr_dict = {'lambda_p': 'lambda_input_p',
                        'lambda_u': 'lambda_input_u',
                        'phi_p':    'phi_input_p',
                        'phi_v':    'phi_input_v'}

        for key, value in setattr_dict.items():
            try:
                # For each key-value pair, function sets self.`key`
                # with `value` section of the input file dict.
                setattr(self, key, grid_def['horizgrid'][value])
            except KeyError:
                # Sets varaible grid to false if any attributes not
                # set from file correctly. Used in um_create_cube
                # to determine whether to create a variable resolution
                # grid.
                setattr(self, 'variable_grid', False)
            else:
                # If all are set correctly, set self.variable_grid to
                # True.
                setattr(self, 'variable_grid', True)


def um_guess_bounds(cube):
    '''
    Estimate area of grid cells by guessing bounds of 'latitude'
    and 'longitude' coordinates and return areas as an array.

    :param cube: the cube from which bounds are to be guessed.
    :type cube: :py:class:`iris.cube.Cube`

    :returns cube: an Iris cube with bounds on coordinates.
    :rtype cube: :py:class:`iris.cube.Cube`

    '''

    # Guess latitude and longitude bounds
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    return cube


def um_sort_cube_bounds(cube):
    '''
    This function modifies the input `iris.cube.Cube` by filling missing
    coordinate bounds, setting boundaries on latitude coordinates to be within
    -90 to 90 degrees, and performs a check and update for the cube being
    circular (global mesh).

    :param cube: the cube to be checked and updated.
    :type cube: :py:class:`iris.cube.Cube`

    '''

    # Fill any missing coordinate bounds
    for coord in cube.dim_coords:
        if not coord.has_bounds():
            coord.guess_bounds()

    # Replace the cube's latitude bounds with a copy that has been clipped
    # to sit within -90 to 90 degrees.
    if len(cube.coord_dims('latitude')) > 0:
        latbounds = cube.coord('latitude').bounds.copy()
        latbounds = numpy.clip(latbounds, -90., 90.)
        cube.coord('latitude').bounds = latbounds

    # Check whether the longitude minimum is exactly 360 degrees less than
    # the maximum, indicating a circular/global mesh.
    # Update the boolean attribute cube.coord('longitude').circular as True.
    if len(cube.coord_dims('longitude')) > 0:
        if cube.coord('longitude').bounds.max() == \
           cube.coord('longitude').bounds.min()+360:
            cube.coord('longitude').circular = True


def um_get_grid(cube):
    '''
    Creates ndarrays lat and lon which contain coordinate information for each
    space in the `iris.cube` grid.
    Returns the latitude point for each longitude as an array with the shape
    [longitude, latitude], as given by cube.shape.
    Also returns the longitude points for each latitudinal line with the
    shape [longitude, latitude], as given by cube.shape.

    :param cube: The cube which provides the latitude, longitude and shape
                 information.
    :type cube: :py:class:`iris.cube.Cube`

    :returns lon: An array of arrays, where each array contains every longitude
                  point.
    :rtype lon: ndarray

    :returns lat: An array of arrays, where each array contains one latitude
                  point for every value.
    :rtype lon: ndarray

    '''

    # `Cube.shape` is part of the `iris.cube.Cube.data` attribute and provides
    # the dimensions of the grid (nx, ny). [1]
    # [1]: https://scitools-iris.readthedocs.io/en/v3.0.0/generated/api/iris/cube.html # noqa

    # Sets lon as a reshaped version of the `iris.cube.Cube` longitude points
    # where cube.shape equals [x, y]. This gives the longitude points at each
    # latitude.
    lon = numpy.resize(cube.coord('longitude').points, (cube.shape))

    # Sets lat as an empty version of the `iris.cube.Cube` shape equal to
    # [longitude, latitude]. The loop gives the latitude points at each
    # latitude.
    lat = numpy.zeros((cube.shape))
    for i in range(cube.shape[1]):
        lat[:, i] = cube.coord('latitude').points

    return lon, lat


def um_gen_area(cube):
    # pylint: disable=line-too-long
    '''
    Area weights are calculated for each lat/lon cell using the built-in
    Iris analysis tool `area_weights`. Each cell is weighted by its own
    area, not as a fraction of the whole.

    See the documentation for `iris.analysis.cartography.area_weights`:
    https://scitools-iris.readthedocs.io/en/v3.0.1/_modules/iris/analysis/cartography.html#area_weights

    :param cube: The Iris cube containing latitude and longitude points to
                be weighted.
    :type cube: :py:class:`iris.cube.Cube`

    :returns: An  array of area weights, with the same dimensions as the cube.
    :rtype: array

    '''
    return iris.analysis.cartography.area_weights(cube)


def um_create_cube(my_grid):
    '''
    Creates a horizontal UM cube from a Grid object.
    Depending on the type of grid points (U/V/P), coordinates
    are shifted to account for staggering.

     --------[V]--------
     |                 |
     |                 |
     |                 |
   [U]       [P]      [U]
     |                 |
     |                 |
     |                 |
     --------[V]--------

    Either creates a non-variable resolution grid or reads coords
    of a variable resolution grid from `my_grid` attributes.

    Latitude coordinate bounds are checked and stop the function executing
    if they are beyond the bounds -90 or 90.

    Longitude coordinate bounds are checked and stop the function
    executing if they are beyond the bounds 0 or 360. This assumes the
    longitude range is given from 0 to 360 degrees.

    :param my_grid: The grid object from which a UM cube will be generated.
    :type my_grid: :py:class:`Grid`

    :returns cube: An Iris cube containing coordinate information for
    :rtype cube: :py:class:`iris.cube.Cube`

    '''
    tiny = 1.0e-10

    if my_grid.variable_grid:
        # Gets cell centre points
        xcoord = my_grid.lambda_p
        ycoord = my_grid.phi_p

        # Get u-points if grid = U
        if my_grid.grid == 'U':
            xcoord = my_grid.lambda_u

        # Get v-points if grid = V
        if my_grid.grid == 'V':
            ycoord = my_grid.phi_v

    else:
        # Calculates cell centre points
        xcoord = numpy.arange(my_grid.xorigin, my_grid.xorigin +
                              (my_grid.npx * my_grid.delx - tiny),
                              my_grid.delx)
        ycoord = numpy.arange(my_grid.yorigin, my_grid.yorigin +
                              (my_grid.npy * my_grid.dely - tiny),
                              my_grid.dely)

        # Shift in x direction for U-staggering
        if my_grid.grid == 'U':
            xcoord -= .5*my_grid.delx
            print(my_grid.npx, xcoord.shape, 'U-grid dims')

        # Shift in y direction for V-staggering
        if my_grid.grid == 'V':
            ycoord = numpy.arange(my_grid.yorigin, my_grid.yorigin +
                                  ((my_grid.npy + 1) * my_grid.dely - tiny),
                                  my_grid.dely)
            print(my_grid.npy, ycoord.shape, 'V-grid dims')
            ycoord -= .5*my_grid.dely
            ycoord[0] = numpy.float32(ycoord[0])
            ycoord[-1] = numpy.float32(ycoord[-1])

    # Assumes longitude range from 0 to 360 degrees
    if min(xcoord) < 0.:
        raise Exception('lowest longitude less than zero')

    if min(ycoord) < -90. or min(ycoord) > 90.:
        print(min(ycoord), min(ycoord))
        print('lowest latitude less than -90 or greater than 90')

    # Creates a 1D, numeric Coord with read-only points and bounds
    mlong = iris.coords.DimCoord(xcoord, standard_name='longitude',
                                 circular=False, units='degrees',
                                 coord_system=iris.coord_systems.GeogCS(
                                     iris.fileformats.pp.EARTH_RADIUS))

    # Creates a 1D, numeric Coord with read-only points and bounds
    mlat = iris.coords.DimCoord(ycoord, standard_name='latitude',
                                units='degrees',
                                coord_system=iris.coord_systems.GeogCS(
                                    iris.fileformats.pp.EARTH_RADIUS))

    # Define the shape of the cube
    data = numpy.zeros((len(ycoord), len(xcoord)))

    # Create the cube
    cube = iris.cube.Cube(data, dim_coords_and_dims=[(mlat, 0), (mlong, 1)],
                          var_name='cube')

    # Add bounds and sort (see um_sort_cube_bounds)
    um_guess_bounds(cube)
    um_sort_cube_bounds(cube)

    return cube


def um_lat_corners(cube):
    '''
    Creates a 3-dimensional array with the shape (lon, lat, 4),
    meaning there are `lon` number of array sets with `lat` number
    of arrays, all of length 4. Each array represents a cell and
    each value represents a corner of the cell.

    [2]----------------[3]
     |                  |
     |                  |
     |                  |
     |                  |
     |                  |
     |                  |
     |                  |
    [0]----------------[1]


    The first and second values in the every 4-point array correspond
    to the minimum latitude of each cell (points 0 and 1), and the
    final two correspond to the maximum (points 2 and 3 on the diagram).

    :param cube: The Iris cube containing latitude boundaries.
    :type cube: :py:class:`iris.cube.Cube`

    :returns corners: A 3-dimensional array where each array set
                      contains the latitude of a cell corner.
    :rtype corners: ndarray

    '''

    # Create and empty grid with matching shape to cube.
    lat_bounds = numpy.zeros((cube.shape))
    # Creates a 3D array with shape (lon, lat, 4)
    corners = numpy.zeros((cube.shape[0], cube.shape[1], 4))

    # Takes the minimum latitude for each cell and stores
    # it in a numpy array which can be broadcast to the shape
    # of `corners`.
    for i in range(cube.shape[1]):
        lat_bounds[:, i] = cube.coord('latitude').bounds[:, 0]
    corners[:, :, 0] = lat_bounds   # Set min bound as corner 0
    corners[:, :, 1] = lat_bounds   # Set min bound as corner 1

    # Takes the maximum latitude for each cell and stores
    # it in a numpy array which can be broadcast to the shape
    # of `corners`.
    for i in range(cube.shape[1]):
        lat_bounds[:, i] = cube.coord('latitude').bounds[:, 1]
    corners[:, :, 2] = lat_bounds   # Set max bound as corner 2
    corners[:, :, 3] = lat_bounds   # Set max bound as corner 3

    return corners


def um_lon_corners(cube):
    '''
    Creates a 3-dimensional array with the shape (lon, lat, 4),
    meaning there are `lon` number of array sets with `lat` number
    of arrays, all of length 4. Each array represents a cell and
    each value represents a corner of the cell.

    [2]----------------[3]
     |                  |
     |                  |
     |                  |
     |                  |
     |                  |
     |                  |
     |                  |
    [0]----------------[1]


    The first and second values in the every 4-point array correspond
    to the minimum longitude of each cell (points 0 and 2), and the
    final two correspond to the maximum (points 1 and 3 on the diagram).

    :param cube: The Iris cube containing longitude boundaries.
    :type cube: :py:class:`iris.cube.Cube`

    :returns corners: A 3-dimensional array where each array set
                      contains the longitude of a cell corner.
    :rtype corners: ndarray

    '''

    # Create and empty grid with matching shape to cube.
    lon_bounds = numpy.zeros((cube.shape))
    # Creates a 3D array with shape (lon, lat, 4)
    corners = numpy.zeros((cube.shape[0], cube.shape[1], 4))

    # Takes the minimum longitude for each cell and stores
    # it in a numpy array which can be broadcast to the shape
    # of `corners`.
    for i in range(cube.shape[0]):
        lon_bounds[i, :] = cube.coord('longitude').bounds[:, 0]
    corners[:, :, 0] = lon_bounds
    corners[:, :, 3] = lon_bounds

    # Takes the maximum longitude for each cell and stores
    # it in a numpy array which can be broadcast to the shape
    # of `corners`.
    for i in range(cube.shape[0]):
        lon_bounds[i, :] = cube.coord('longitude').bounds[:, 1]
    corners[:, :, 1] = lon_bounds
    corners[:, :, 2] = lon_bounds

    return corners


def transform(fin, name):
    '''
    Identifies the number of dimensions held by the input array and then
    transposes the array into either a one-dimensional array or a
    two-dimensional array with the shape (x, 4) if one of the shape's
    dimensions equals 4. This is a proxy check to identify arrays
    containing lat/lon coordinates for cell corners.
    If the array is already one-dimensional, the input array is returned
    without being transformed.

    :param fin: An array of either 1, 2, or 3-dimensions.
    :type fin: ndarray

    :param str name: The name to which the array will be referred.

    :returns fout: A one- or two-dimensional array.
    :rtype fout: ndarray

    '''
    # Get shape of array
    shape = fin.shape

    # Check if array is 2-dimensional
    if len(shape) == 2:
        # If there are 4 array sets of length x, returns x arrays each
        # of length 4.
        if shape[0] == 4:
            # from shape (x, y, z) to (y, z, x)
            fout = numpy.transpose(fin)
        # If there are x array sets each of length 4, no tranformations
        # are made.
        elif shape[1] == 4:
            fout = fin
        # If not in shape (x, 4), assumed to be a cell centres array.
        # Fill 1D array with values of `fin`.
        else:
            fout = numpy.zeros((shape[0]*shape[1]), dtype=numpy.float64)
            fout = fin.flatten()
    # Check if array is 3-dimensional
    elif len(shape) == 3:
        # If shape i (4, x, y), transpose into shape (x, y, 4) then
        # mutate fin to 2D shape (yx, 4).
        if shape[0] == 4:
            fin = fin.transpose((1, 2, 0))
            shape_o = shape
            shape = fin.shape
            print('transforming corner data for ', name, ' from ', shape_o,
                  ' to', shape)
        fout = numpy.zeros((shape[0]*shape[1], 4), dtype=numpy.float64)
        fout = numpy.reshape(fin, (shape[0]*shape[1], 4))
    # Check if array is 1-dimensional
    elif len(shape) == 1:
        # make no changes.
        fout = fin
    else:
        raise Exception('problem in transform')
    return fout


def transform_and_write(filename, my_grid):
    '''
    Writes attributes from a Grid object to a specified file in NetCDF format.
    For each variable, a section is created in the output file with the data-
    type and dimensions required, a name is set, and an array is written.
    For variables containing coordinate data, the function applies the
    `transform` function to cell centre and cell corner arrays to reduce the
    number of array dimensions by 1 (unless the array is 1D in which case it
    applies no transformations).

    :param str filename: A path to the desired output file, with extension.

    :param my_grid: A Grid object with attributes to be written to file.
    :type my_grid: :py:class:`Grid`

    '''

    # Set up variables
    rank = len(my_grid.shape)
    mrank = my_grid.shape
    corners = 4

    # Check grid is 2D
    if rank == 2:
        grid_size = my_grid.shape[0]*my_grid.shape[1]
    else:
        raise Exception('transform_and_write: ' +
                        'wrong number of dimensions in input data')

    # Set up NetCDF file
    fileout = filename
    outfile = Dataset(fileout, 'w')

    # Create NetCDF4 dimensions to be used by variables
    outfile.createDimension('grid_rank', numpy.int32(rank))
    outfile.createDimension('grid_size', numpy.int32(grid_size))
    outfile.createDimension('grid_corners', numpy.int32(corners))

    # Create `grid_dims` with dimensions: `grid_rank`
    rankout = outfile.createVariable('grid_dims',
                                     numpy.dtype('int32').char,
                                     ('grid_rank', ))
    rankout.long_name = "grid_dims"
    rankout[:] = numpy.asarray(mrank).astype(numpy.int32)

    # Create `grid_center_lat` with dimensions: `grid_size`
    latout = outfile.createVariable('grid_center_lat',
                                    numpy.dtype('float64').char,
                                    ('grid_size', ))
    latout.long_name = "grid_center_lat"
    latout.units = str(my_grid.ucoor)
    latout[:] = transform(my_grid.dlat, 'lat')

    # Create `grid_center_lon` with dimensions: `grid_size`
    lonout = outfile.createVariable('grid_center_lon',
                                    numpy.dtype('float64').char,
                                    ('grid_size', ))
    lonout.long_name = "grid_center_lon"
    lonout.units = str(my_grid.ucoor)
    lonout[:] = transform(my_grid.dlon, 'lon')

    # Create `grid_imask` with dimensions: `grid_size`
    maskout = outfile.createVariable('grid_imask',
                                     numpy.dtype('int32').char,
                                     ('grid_size', ))
    maskout.long_name = "grid_imask"
    maskout.units = "unitless"
    maskout[:] = numpy.int32(1)

    # Create `grid_area` with dimensions: `grid_size`
    if my_grid.l_area != 'no':
        srfout = outfile.createVariable('grid_area',
                                        numpy.dtype('float64').char,
                                        ('grid_size', ))
        srfout.long_name = "grid surface m^2"
        srfout.units = str(my_grid.uarea)
        srfout[:] = transform(my_grid.darea, 'area')

    # Create `grid_corner_lat` with dimensions: `grid_size` and `grid_corners`
    claout = outfile.createVariable('grid_corner_lat',
                                    numpy.dtype('float64').char,
                                    ('grid_size', 'grid_corners'))
    claout.long_name = "grid_corner_lat"
    claout.units = str(my_grid.ucoor)
    claout[:, :] = transform(my_grid.dcla, 'cla')

    # Create `grid_corner_lon` with dimensions: `grid_size` and `grid_corners`
    cloout = outfile.createVariable('grid_corner_lon',
                                    numpy.dtype('float64').char,
                                    ('grid_size', 'grid_corners'))
    cloout.long_name = "grid_corner_lon"
    cloout.units = str(my_grid.ucoor)
    cloout[:, :] = transform(my_grid.dclo, 'clo')

    outfile.title = my_grid.vname
    outfile.close()


def unrotate_coords(lon_rot, lat_rot, coords):
    '''
    Checks that latitude and longitude coordinates provided are the same
    shape, creates empty numpy arrays of their shape and then applies the
    mule utility `CoordRotator` to populate the arrays with unrotated
    coordinates.

    This function is used in `get_data` to unrotate cell center and cell
    corner coordinates stored in a Grid object.

    :param lon_rot: An array containing longitude coordinates to be
                    (un)rotated onto a regular lat-lon grid.
    :type lon_rot:  ndarray

    :param lat_rot: An array containing latitude coordinates to be
                    (un)rotated onto a regular lat-lon grid.
    :type lat_rot:  ndarray

    :param coords: A CoordRotator object from the um_utils library. This
                   is initialised in the `get_data` function.
    :type coords: :py:class:`mule.um_utils.CoordRotator`

    :returns lon: An array of unrotated coordinates.
    :rtype lon: ndarray

    :returns lat: An array of unrotated coordinates.
    :rtype lat: ndarray

    '''
    # Check grids are equal size and shape
    if lon_rot.shape == lat_rot.shape:
        lon = numpy.empty(lon_rot.shape)
        lat = numpy.empty(lat_rot.shape)
    else:
        raise Exception('unrotate_coords: lon/lat arrays have different sizes')

    for index, _ in numpy.ndenumerate(lon_rot):
        x = lon_rot[index]
        y = lat_rot[index]
        # Rotate a pair of coordinates from the grid defined by this object
        # onto a regular lat-lon grid.
        lon[index], lat[index] = coords.unrotate(x, y)

        if lon[index] < 0.0:
            lon[index] = lon[index] + 360.0

    return lon, lat


def get_data(my_grid):
    '''
    Utilises all of the functions defined here to generate a UM grid to:

    1. Instantiate a GRID object.
    2. Set GRID instance attributes by reading from namelist and variable
       resolution mesh files.
       - Namelist file supplies rotation, poles, delta and lambda targets and
         origins, and x and y origins.
       - Variable resolution mesh file supplies row and column dependent
         constants lamba_p, lambda_u, phi_p, and phi_v.
       - If the variable resolution mesh file can be parsed and sets all
         the aforementioned GRID atributes to a value other than None,
         `GRID.variable_grid` is set to True.
    3. Create an `iris.cube.Cube` instance from the GRID instance.
       - If GRID.variable_grid == True, x and y coordinates are set from the
         constants provided by the variable mesh file.
       - If GRID.variable_grid is False, x and y coordinates are generated
         from x and y origins at an interval defined by delta_x and delta_y.
         This creates a fixed resolution mesh.
    4. If the coordinates are rotated, they are unrotated.
    5. Coordinates of the corners for each cell are extracted and unrotated.
    6. The area of each grid cell is calculated.
    7. The grid shape is defined from the latitude of cell centers.
    8. The grid is written to a NetCDF file at GRID_PATH_UM.

    '''
    # Set my_grid attributes
    my_grid.vname = "um_grid"
    # .grid is set in `weights_gen_wrapper.py`: will loop over ['P', 'U', 'V']
    my_grid.grid = os.environ.get('GRID')
    my_grid.l_area = os.environ.get('LAREA')

    # Sets attributes from class method that parses namelist and
    # variable res files. If `VARRES_FILE` is None, no attributes are
    # set and the my_grid attribute `variable_grid` defaults to `False`.
    my_grid.read_namelist(os.environ.get('NLIST_FILE'))

    # Must check if VARRES_FILE exists since .read_accessa_grid method relies
    # on f90nml.read() which cannot accept None in filename argument position.
    varres_file = os.environ.get('VARRES_FILE')
    if varres_file:
        my_grid.read_accessa_grid(varres_file)

    attrs = vars(my_grid)
    print(', '.join("%s: %s" % item for item in attrs.items()))

    cube = um_create_cube(my_grid)

    rotated_coords = CoordRotator(my_grid.pole_lon, my_grid.pole_lat)

    # Get longitude and latitude of cell centres
    x, y = um_get_grid(cube)
    if my_grid.rotated_grid:
        my_grid.dlon, my_grid.dlat = unrotate_coords(x, y, rotated_coords)
    else:
        my_grid.dlon, my_grid.dlat = x, y

    # Get longitude and latitude of cell corners
    x = um_lon_corners(cube)
    y = um_lat_corners(cube)
    if my_grid.rotated_grid:
        my_grid.dclo, my_grid.dcla = unrotate_coords(x, y, rotated_coords)
    else:
        my_grid.dclo, my_grid.dcla = x, y

    # Get area
    my_grid.darea = um_gen_area(cube)

    # Set units
    my_grid.ucoor = "degrees"
    my_grid.uarea = "m2"

    # Set shape of grid
    my_grid.shape = numpy.asarray(my_grid.dlat.shape[::-1])

    return my_grid


if __name__ == "__main__":
    um_grid = GRID()
    um_grid = get_data(um_grid)
    transform_and_write(os.environ["GRID_PATH_UM"], um_grid)
