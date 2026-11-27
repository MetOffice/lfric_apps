#!/usr/bin/env python
##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Some content in this file was generated or refactored with assistance from
# - Met Office Github Copilot Enterprise (Claude Opus 5.5)
'''
Copy an LFRic UGRID mesh file and add a face-located mask field to it,
taken from another NetCDF file. The field is written on the mesh's
n<mesh>_face dimension with _FillValue=0. Optionally the mask can be
inverted first (e.g. to turn a land mask into a sea mask).
'''
import argparse
import shutil

import numpy as np
from netCDF4 import Dataset

# Set the value that will represent masked cells (i.e. _FillValue in netcdf)
FILL_VALUE = 0.0
# Metadata carried over from the source mask variable, where present
COPY_ATTRS = ['units', 'standard_name', 'long_name', 'um_stash_source']


def add_mask_to_mesh(mesh_in, mesh_name, mask_file, mask_variable_name,
                     mesh_out, out_var, invert=False):
    '''
    Write mesh_out: a copy of mesh_in with out_var added on the faces.

    param str mesh_in: Path to the input LFRic UGRID mesh file.
    param str mesh_name: Name of the mesh topology within mesh_in.
    param str mask_file: Path to the NetCDF file containing the mask.
    param str mask_variable_name: Name of the mask variable in mask_file.
    param str mesh_out: Path of the masked mesh file to write.
    param str out_var: Name to give the mask variable in mesh_out.
    param bool invert: If True, swap masked (0) and unmasked (1) cells.
    '''
    # Work on a copy so the original mesh is left untouched
    shutil.copyfile(mesh_in, mesh_out)

    with Dataset(mask_file, 'r') as src:
        if mask_variable_name not in src.variables:
            raise ValueError(f"Variable '{mask_variable_name}' "
                             f"not found in {mask_file}")
        src_var = src.variables[mask_variable_name]
        # Read raw values as a plain array rather than a numpy masked array
        src_var.set_auto_mask(False)
        # Drop singleton dimensions (e.g. time) to leave a 1D face array
        data = np.squeeze(src_var[...])
        attrs = {name: src_var.getncattr(name) for name in COPY_ATTRS
                 if name in src_var.ncattrs()}

    if invert:
        data = np.where(data == FILL_VALUE, 1.0, FILL_VALUE)
        # The source metadata describes the un-inverted mask
        attrs = {'units': '1',
                 'long_name': f'inverse of {mask_variable_name} '
                              f'from {mask_file}'}

    with Dataset(mesh_out, 'a') as mesh:
        # LFRic UGRID naming convention for the face dimension
        face_dim = 'n' + mesh_name + '_face'
        if face_dim not in mesh.dimensions:
            raise ValueError(f"Dimension '{face_dim}' not found in {mesh_in}")
        n_faces = len(mesh.dimensions[face_dim])
        if data.shape != (n_faces,):
            raise ValueError(
                f"Shape of '{mask_variable_name}' {data.shape} does not match "
                f"'{face_dim}' ({n_faces})")

        var = mesh.createVariable(out_var, 'f8', (face_dim,),
                                  fill_value=FILL_VALUE)
        var.setncatts(attrs)
        # UGRID attributes linking the field to the mesh faces
        var.mesh = mesh_name
        var.location = 'face'
        var[:] = data

    print(f"Added '{out_var}' to {mesh_out}: "
          f"{int(np.count_nonzero(data != FILL_VALUE))} of {n_faces} "
          f"faces unmasked")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh-in', required=True)
    parser.add_argument('--mesh-name', required=True)
    parser.add_argument('--mask-file', required=True)
    parser.add_argument('--mask-variable-name', required=True)
    parser.add_argument('--mesh-out', required=True)
    parser.add_argument('--out-var', required=True)
    parser.add_argument('--invert', action='store_true',
                        help='Swap masked (0) and unmasked (non-zero) cells')
    args = parser.parse_args()

    add_mask_to_mesh(args.mesh_in, args.mesh_name, args.mask_file,
                     args.mask_variable_name, args.mesh_out, args.out_var,
                     args.invert)
