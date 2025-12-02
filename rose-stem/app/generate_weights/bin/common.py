#!/usr/bin/env /data/users/support/ants/_dev/environments/ukcp18/bin/python2.7
#-----------------------------------------------------------------------------
# (C) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
#-----------------------------------------------------------------------------
"""module with common functions"""
import os
import sys
import numpy as np
from netCDF4 import Dataset

def get_ncdata(in_file, sgrid):
    """Read data from file"""
    try:
        grid = in_file.variables[sgrid]
    except (IOError) as err:
        print("I/O error({0}): {1}".format(err.errno, err.strerror))
    except:
        raise Exception("Unexpected error:")
    return grid

def get_var(co_file, vname, lmask):
    """get values fron netcdf file"""
    tmp = get_ncdata(co_file, vname)
    dims = len(tmp.dimensions)
    if 'units' in tmp.__dict__:
        outu = tmp.units
    else:
        outu = None
    if dims == 1:
        out = tmp[:]
    elif dims == 2:
        out = tmp[:, :]
    elif dims == 3:
        out = tmp[:, :, :]
    elif dims == 4:
        out = tmp[:, :, :, :]
    else:
        raise Exception('get_var: Too many dimensions')
    if lmask and dims > 2:
        print('Too many dimensions for ', vname, ': ', out.shape, ' leaving only last two !!!')
        if dims == 3:
            out = tmp[0, :, :]
        else:
            out = tmp[0, 0, :, :]
    return np.squeeze(out), outu

def create_file_name(lflag, lfull):
    """Create input file name"""
    g1grid = os.environ.get('GRID1_GRID')
    g2grid = os.environ.get('GRID2_GRID')
    g1model = os.environ.get('GRID1_MODEL')
    g2model = os.environ.get('GRID2_MODEL')

    if lflag:
        src_vname = set_onames(g1grid.upper(), g1model)
        dst_vname = set_onames(g2grid.upper(), g2model)
        lignm_src = os.environ.get('GRID1_IGMASK') == 'yes'
        lignm_dst = os.environ.get('GRID2_IGMASK') == 'yes'
    else:
        dst_vname = set_onames(g1grid.upper(), g1model)
        src_vname = set_onames(g2grid.upper(), g2model)
        lignm_src = os.environ.get('GRID2_IGMASK') == 'yes'
        lignm_dst = os.environ.get('GRID1_IGMASK') == 'yes'

    int_method = os.environ.get('INT_METHOD')
    int_norm = os.environ.get('INT_NORM')
    pwd = os.environ.get('OUT_DIR')

    if not (pwd and pwd.strip()):
        shrd = os.environ.get('CYLC_SUITE_SHARE_DIR')
        cycle = os.environ.get('CYLC_TASK_CYCLE_TIME')
        pwd = os.path.join(shrd, 'data', cycle)
    mname = pwd + '/rmp_' + src_vname + '_to_' + dst_vname  \

    mname = mname + '_' + set_method(int_method, int_norm)

    if lignm_src and lignm_dst:
        mname = mname + '_nomask'
    elif lignm_src:
        mname = mname + '_nomasksrc'
    elif lignm_dst:
        mname = mname + '_nomaskdst'

    if lfull:
        mname = mname + '.nc'

    return mname, src_vname, dst_vname

def set_method(method, norm):
    """Set name based on interpolation method"""

    if 'CONSERVE2ND' in method.upper():
        intn = "conservative2nd"
    elif "CONSERVE" in method.upper():
        intn = "conservative1st"
    else:
        intn = method.lower()

    if ('CONSERVE' in method.upper()) and (norm.lower() != 'fracarea')\
                                      and (norm.lower() != 'dstarea'):
        sys.exit("Wrong normalization for conservative remapping.\
 Only fracarea and dstarea are valid")

    if 'fracarea' in norm.lower():
        intn = intn + "_fracarea"
    elif "dstarea" in norm.lower():
        intn = intn + "_dstarea"
    else:
        print("Ignore normalization")

    return intn


def set_onames(gtp, mdl):
    """Set variable names for output files"""
    if mdl == 'LFRIC':
        if gtp == 'T':
            sname = 'lfricT'
        else:
            raise Exception('LFRic can handle Grid T only')
    elif mdl == 'NEMO':
        if gtp == 'T':
            sname = 'oceT'
        elif gtp == 'U':
            sname = 'oceU'
        elif gtp == 'V':
            sname = 'oceV'
        elif gtp == 'F':
            sname = 'oceF'
        else:
            raise Exception('Unrecognised Grid: ', gtp, 'for model ', mdl)
        print('NEMO GRID: ', gtp)
    else:
        raise Exception('Unrecognised Model: ', mdl)
    return sname

def ord2_corr(src_field, src_mask, rank_src):
    """calculate lon and lat gradients to correct
       destination values for conservative
       weights """
#we need
    src_mask_r = np.reshape(src_mask, (rank_src[1], rank_src[0]))
    grad_lon = np.zeros((rank_src[0]*rank_src[1], 1), dtype=np.float64)
    grad_lat = np.zeros((rank_src[0]*rank_src[1], 1), dtype=np.float64)
#find unique source grid points
    for j in range(rank_src[1]):
        for i in range(rank_src[0]):
            delew = .5
            delsn = .5
            ipt = j*rank_src[0]+i
            if(src_mask_r[j, i]):
#find N,E,S,W points for col[i], use fortran indexing (starts with 1)
#follow scrip approach. Incorrect for NEMO grid for global model!
#find i-1 and i+1
#i-1
                ipt_m1 = i - 1
                if ipt_m1 < 0:
                    ipt_m1 = rank_src[0] -1 + ipt_m1
#account for mask i-1
                if(src_mask_r[j, ipt_m1] == 0):
                    ipt_m1 = i
                    delew = 1.
#i+1
                ipt_p1 = i + 1
                if ipt_p1 > rank_src[0] -1:
                    ipt_p1 = ipt_p1 - rank_src[1] - 1
#account for mask i+1
                if(src_mask_r[j, ipt_p1] == 0):
                    ipt_p1 = i
                    delew = 1.
#j+1
                jpt_p1 = j + 1
                if(jpt_p1 > rank_src[1]) - 1:
                    jpt_p1 = j
                    delsn = 1.
#account for mask j+1
                if src_mask_r[jpt_p1, i] == 0:
                    jpt_p1 = j
                    delsn = 1.
#j-1
                jpt_m1 = j - 1
                if jpt_m1 < 0:
                    jpt_m1 = j
                    delsn = 1.
#account for mask j-1
                if src_mask_r[jpt_m1, i] == 0:
                    jpt_m1 = j
                    delsn = 1.
#shift index python (0.0) fortran (1,1)
                ipt_w = j*rank_src[0]+ipt_m1
                ipt_e = j*rank_src[0]+ipt_p1
                jpt_n = jpt_p1*rank_src[0]+i
                jpt_s = jpt_m1*rank_src[0]+i
                grad_lon[ipt] = delew*(src_field[ipt_e] - src_field[ipt_w])
                grad_lat[ipt] = delsn*(src_field[jpt_n] - src_field[jpt_s])
    return  grad_lon, grad_lat


def get_remapping_info_oasis(fname):
    """Get remapping info from file"""
    in_file = Dataset(fname+'.nc', 'r')
    col = get_ncdata(in_file, 'src_address')[:]
    row = get_ncdata(in_file, 'dst_address')[:]
    swgh = get_ncdata(in_file, 'remap_matrix')[:, :]
    x_src = get_ncdata(in_file, 'src_grid_center_lon')[:]
    y_src = get_ncdata(in_file, 'src_grid_center_lat')[:]
    x_dst = get_ncdata(in_file, 'dst_grid_center_lon')[:]
    y_dst = get_ncdata(in_file, 'dst_grid_center_lat')[:]
    rank_dst = get_ncdata(in_file, 'dst_grid_dims')[:]
    rank_src = get_ncdata(in_file, 'src_grid_dims')[:]
    mask_a = get_ncdata(in_file, 'src_grid_imask')[:]
    mask_b = get_ncdata(in_file, 'dst_grid_imask')[:]
    frac_b = get_ncdata(in_file, 'dst_grid_frac')[:]
    in_file.close()

    return col, row, swgh, x_src, y_src, x_dst, \
        y_dst, rank_src, rank_dst, mask_a, mask_b, frac_b

def interpolate(row, swgh, col, src_field, fshape, frac_b, src_mask,\
                                               rank_src, fname, linv):
    """Use weights to interpolate test function between grids"""

    n_s = col.shape[0]
    dst_field = np.zeros(fshape)

    for i in range(n_s):
        dst_field[row[i]-1] = dst_field[row[i]-1]+swgh[i, 0]*src_field[col[i]-1]
#normalization
    if 'DSTAREA' in fname and 'CONSERV' in fname:
        for i in range(frac_b.shape[0]):
            if frac_b[i] != 0.0:
                dst_field[i] = dst_field[i]/frac_b[i]
#swap values for mask
    if linv:
        dst_field = 1. - dst_field

    return dst_field
