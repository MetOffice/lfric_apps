#!/usr/bin/env /data/users/support/ants/_dev/environments/ukcp18/bin/python2.7
#-----------------------------------------------------------------------------
# (C) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
#-----------------------------------------------------------------------------
"""Create grid information for models for use with ESMF regridding program"""
import os
from netCDF4 import Dataset
import numpy
import iris
import iris.fileformats
import iris.analysis
import iris.coord_systems
import f90nml
from common import get_var, set_onames

class GRID:
    """Holds information from suite"""
    def __init__(self):
        self.lgrid = None
        self.vname = None
        self.vname_lo = None
        self.vname_la = None
        self.fname = None
        self.vname_mask = None
        self.l_corners = None
        self.vname_corners_lo = None
        self.vname_corners_la = None
        self.vname_area = None
        self.dlon = None
        self.ulon = None
        self.dlat = None
        self.ulat = None
        self.dclo = None
        self.uclo = None
        self.dcla = None
        self.ucla = None
        self.uarea = None
        self.dmask = None
        self.umask = None
        self.shape = None
        self.nfold = False

def LFRIC_get_corners(in_co_file, my_grid):
    #sorners info
    my_grid.dclo, my_grid.uclo = get_var(in_co_file, my_grid.vname_corners_lo, False)
    my_grid.dcla, my_grid.ucla = get_var(in_co_file, my_grid.vname_corners_la, False)
    node_ids, _ = get_var(in_co_file, my_grid.id_corners, False)
    node_lon, _ = get_var(in_co_file, my_grid.vname + '_node_x', False)
    node_lat, _ = get_var(in_co_file, my_grid.vname + '_node_y', False)
    lon_tc = node_lon[node_ids[:, :] - 1]
    lat_tc = node_lat[node_ids[:, :] - 1]
    return lon_tc, lat_tc, my_grid.uclo, my_grid.ucla


def corners_transformation(field):
    """Transform corners to correct order"""
    return field.transpose((1, 2, 0))


def transform(fin, name):
    """Remove one dimension from the array"""
    shape = fin.shape
    if len(shape) == 2:
        if shape[0] == 4:
            fout = numpy.transpose(fin)
        elif shape[1] == 4:
            fout = fin
        else:
            fout = numpy.zeros((shape[0]*shape[1]), dtype=numpy.float64)
            fout = fin.flatten()
    elif len(shape) == 3:
        if shape[0] == 4:
            fin = corners_transformation(fin)
            shape_o = shape
            shape = fin.shape
            print('transforming corner data for ', name, ' from ', shape_o, \
                  ' to ', shape)
        fout = numpy.zeros((shape[0]*shape[1], 4), dtype=numpy.float64)
        fout = numpy.reshape(fin, (shape[0]*shape[1], 4))
    elif len(shape) == 1:
        fout = fin
    else:
        raise Exception('problem in transform')
    return fout


def transform_and_write(my_grid):
    """Transform arrays and write netcdf output"""
    rank = len(my_grid.shape)
    mrank = my_grid.shape
    print(rank, ' RANK')
    if rank == 2:
        grid_size = my_grid.shape[0]*my_grid.shape[1]
    elif rank == 1:
        grid_size = my_grid.shape[0]
    else:
        raise Exception('transform_and_write: Wrong number of dimensions in\
 input data')

    corners = 4
    fileout = my_grid.vname+'.nc'
    outfile = Dataset(fileout, 'w')
    outfile.createDimension('grid_rank', numpy.int32(rank))
    outfile.createDimension('grid_size', numpy.int32(grid_size))
    outfile.createDimension('grid_corners', numpy.int32(corners))
    rankout = outfile.createVariable('grid_dims',
                                     numpy.dtype('int32').char,
                                     ('grid_rank', ))
    rankout.long_name = "grid_dims"
    rankout.units = "none"
    rankout[:] = numpy.asarray(mrank).astype(numpy.int32)
    latout = outfile.createVariable('grid_center_lat',
                                    numpy.dtype('float64').char,
                                    ('grid_size', ))
    latout.long_name = "grid_center_lat"
    latout.units = my_grid.ulat
    latout[:] = transform(my_grid.dlat, 'lat')
    lonout = outfile.createVariable('grid_center_lon',
                                    numpy.dtype('float64').char,
                                    ('grid_size', ))
    lonout.long_name = "grid_center_lon"
    lonout.units = my_grid.ulon
    lonout[:] = transform(my_grid.dlon, 'lon')
    maskout = outfile.createVariable('grid_imask',
                                     numpy.dtype('int32').char,
                                     ('grid_size', ))
    maskout.long_name = "grid_imask"
    maskout.units = "unitless"
    print("Set mask to 1 everywhere.")
    maskout[:] = numpy.int32(1)

    print('Area not generated for model: ', my_grid.vname)

    claout = outfile.createVariable('grid_corner_lat',
                                    numpy.dtype('float64').char,
                                    ('grid_size', 'grid_corners'))
    claout.long_name = "grid_corner_lat"
    claout.units = my_grid.ucla
    claout[:, :] = transform(my_grid.dcla, 'cla')
    cloout = outfile.createVariable('grid_corner_lon',
                                    numpy.dtype('float64').char,
                                    ('grid_size', 'grid_corners'))
    cloout.long_name = "grid_corner_lon"
    cloout.units = my_grid.uclo
    cloout[:, :] = transform(my_grid.dclo, 'clo')
    outfile.title = my_grid.vname
    outfile.close()
    return


def get_env_info(my_grid, src):
    """Get information set in suite"""
    my_grid.fname = os.environ.get(src + '_GRID_PATH')

    my_grid.vname = os.environ.get(src + '_MESH_NAME')
    my_grid.vname_lo = my_grid.vname + '_face_x'
    my_grid.vname_la = my_grid.vname + '_face_y'
    my_grid.vname_corners_lo = my_grid.vname + '_node_x'
    my_grid.vname_corners_la = my_grid.vname +'_node_y'
    my_grid.id_corners = my_grid.vname + '_face_nodes'

    return my_grid

def validate_input(my_grid):
    """check sizes of arrays !!!!!!
       check units"""
    print('--------------validate---------------------')
    if my_grid.ulon is None:
        print(' Lon/Lat units not found')
        if max(numpy.abs(my_grid.dlon).max(),
               numpy.abs(my_grid.dlon).min()) > 6.5:
            my_grid.ulon = 'degrees'
            my_grid.ulat = 'degrees'
            print('  Based on values setting to degrees')
        else:
            my_grid.ulon = 'radians'
            my_grid.ulat = 'radians'
            print(' Based on values setting to radians. MAY cause problems\
                    with ESMF')
    else:
        if "deg" in (my_grid.ulon).lower():
            my_grid.ulon = 'degrees'
            my_grid.ulat = 'degrees'
        elif "rad" in (my_grid.ulon).lower():
            my_grid.ulon = 'radians'
            my_grid.ulat = 'radians'
        else:
            raise Exception('Unrecognised lon/lat units: ', my_grid.ulon)
    if my_grid.uclo is None:
        print(' Lon/Lat corners units not found')
        if max(numpy.abs(my_grid.dclo).max(),
               numpy.abs(my_grid.dclo).min()) > 6.5:
            my_grid.uclo = 'degrees'
            my_grid.ucla = 'degrees'
            print('  Based on values setting corners units to degrees')
        else:
            my_grid.uclo = 'radians'
            my_grid.ucla = 'radians'
            print('  Based on values setting corners units to radians')
    else:
        if "deg" in (my_grid.uclo).lower():
            my_grid.uclo = 'degrees'
            my_grid.ucla = 'degrees'
        elif "rad" in (my_grid.uclo).lower():
            my_grid.uclo = 'radians'
            my_grid.ucla = 'radians'
        else:
            raise Exception('Unrecognised lon/lat corners unit: ',
                            my_grid.ulon)
    if my_grid.uclo != my_grid.ulon:
        raise Exception('Problem with input data: longitude (',
                        my_grid.ulon, ') and longitude corners (',
                        my_grid.uclo, ') units differ')
    if my_grid.vname_area is not None:
        if my_grid.uarea is None:
            my_grid.uarea = "m^2"
            print('Setting area units to m2. If wrong add \
                   units to Netcdf input file!!')
    print('--------------validate info---------------------')
    print('vname', my_grid.vname)
    print('fname', my_grid.fname)
    print('vname_lo', my_grid.vname_lo)
    print('vname_la', my_grid.vname_la)
    print('vname_mask', my_grid.vname_mask)
    print('l_corners', my_grid.l_corners)
    print('vname_corners_lo', my_grid.vname_corners_lo)
    print('vname_corners_la', my_grid.vname_corners_la)
    print('vname_area', my_grid.vname_area)
    print('dlon', (my_grid.dlon).shape, my_grid.ulon)
    print('dlat', (my_grid.dlat).shape, my_grid.ulat)
    print('dclo', (my_grid.dclo).shape, my_grid.uclo)
    print('dcla', (my_grid.dcla).shape, my_grid.ucla)
    print('dmask', (my_grid.dmask).shape, my_grid.umask)
    print('ingore n-fold (NEMO only) ', my_grid.nfold)
    print('----------end validate---------------------')
    return


def get_data(my_grid):
    """Get data as requested in suite"""
    co_file = Dataset(my_grid.fname, 'r')
    # lon lat
    my_grid.dlon, my_grid.ulon = get_var(co_file, my_grid.vname_lo, False)
    my_grid.dlat, my_grid.ulat = get_var(co_file, my_grid.vname_la, False)
    my_grid.dclo, my_grid.dcla, my_grid.uclo, \
    my_grid.ucla = LFRIC_get_corners(co_file, my_grid)
    my_grid.vname_corners_lo = "atm.clo"
    my_grid.vname_corners_la = "atm.cla"
    co_file.close()

    my_grid.shape = numpy.asarray(my_grid.dlat.shape[::-1])

    # mask always from file or ignore
    my_grid.dmask = 1 + numpy.zeros(my_grid.dlon.shape)
    return my_grid

if __name__ == "__main__":
    SRC_GRID = GRID()
    SRC_GRID = get_env_info(SRC_GRID, 'SRC')
    SRC_GRID = get_data(SRC_GRID)
    validate_input(SRC_GRID)
    # end processing source grid
    print('SOURCE grid ')
    DST_GRID = GRID()
    DST_GRID = get_env_info(DST_GRID, 'DST')
    DST_GRID = get_data(DST_GRID)
    # end processing destination grid
    validate_input(DST_GRID)
    print('DESTINATION grid ')

    print('Save SRC -> DST grid info')
    transform_and_write(SRC_GRID)
    transform_and_write(DST_GRID)
