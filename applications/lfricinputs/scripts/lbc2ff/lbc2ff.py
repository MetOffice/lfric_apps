#!/usr/bin/env python3
#
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file LICENCE.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Converts an LBC file into standard fieldsfile.

Takes a single positional argument - input filename. Script will output a
fieldsfile in the working directory with the same name as the input file with
'.ff' appended.

'''
# pylint: disable=import-error
import argparse
import mule
from mule.lbc import LBCToMaskedArrayOperator
from lbc_stash_map import LBC_STASH_MAP


def main():
    '''
    This function parses an LBC file from a path specified on the command line
    and creates a FieldsFile from the parsed LBC by calling the
    `create_ff_from_lbc` function.
    The LBC fields are then each transposed from a 1-dimensional LBC array into
    a 3-dimensional arrray of levels, rows, and columns and the array is
    sliced to remove the halo regions. The FieldsFile created in this function
    then has each transposed and trimmed field appended to it and is finally
    written to an outfile of the same name as the LBC but with the `.ff`
    extension appended.

    '''
    parser = argparse.ArgumentParser(usage=str('Convert an LBC file into a' +
                                     ' fieldsfile'))
    parser.add_argument("input_filename", help=argparse.SUPPRESS)
    args = parser.parse_args()
    input_filename = args.input_filename

    # Open file
    lbc = mule.lbc.LBCFile.from_file(input_filename)
    fldfle = create_ff_from_lbc(lbc)

    # Use the mule provided LBC operator. This converts the 1d LBC array
    # (which contains all levels) into a standard 3d array of levels, rows
    # and columns.
    lbc_to_masked = LBCToMaskedArrayOperator()
    for field in lbc.fields:
        field = lbc_to_masked(field)
        ncols = field.lbnpt
        nrows = field.lbrow
        # Rim and halo widths
        halo_code = field.lbuser3
        rimwidth = int(halo_code // 10000)
        halo_ns = int(halo_code - rimwidth * 10000) // 100
        halo_ew = int(halo_code - rimwidth * 10000 - halo_ns * 100)
        for level_num, single_level in enumerate(field.get_data(), 1):
            # Copy whole field to get metadata
            field_2d = field.copy()
            # Data provider is mule method for telling new field
            # where to get data from. Slice the array to remove halo regions
            # as not needed in LFRic LBC file
            array_provider = mule.ArrayDataProvider(
                single_level.filled(mule._REAL_MDI)[
                    halo_ns:nrows + halo_ns,
                    halo_ew:ncols + halo_ew])
            field_2d.set_data_provider(array_provider)
            # Update the stash code to the standard prognostic version
            field_2d.lbuser4 = LBC_STASH_MAP[field_2d.lbuser4]
            # Update level number
            field_2d.lblev = level_num
            # Update lbhem variable to be what is expected by fieldsfile
            field_2d.lbhem = fldfle.fixed_length_header.horiz_grid_type % 100
            fldfle.fields.append(field_2d)
    # Set dataset version
    fldfle.fixed_length_header.data_set_format_version = 20
    # Update dataset type to be fieldsfile
    fldfle.fixed_length_header.dataset_type = 3
    fldfle.to_file(input_filename + ".ff")


def create_ff_from_lbc(lbc):
    '''
    This function takes in a pre-parsed LBC file that contains all the
    information required to describe a Unified Model (UM) mesh. The LBC object
    has its data accessed and copied into a FieldsFile which is initialised
    within the scope of this function.

    param lbc: A pre-parsed LBC file.
    type lbc: :py:class:`mule.lbc.LBCFile`

    return fldfle: A FieldsFile containing identical information to the input
               LBC file.
    rtype fldfle: :py:class:`mule.FieldsFile`

    '''
    # Assign Classes to variables so they do not need line breaks (pycodestyle)
    FF_LDC = mule.ff.FF_LevelDependentConstants     # Level
    FF_RDC = mule.ff.FF_RowDependentConstants       # Row
    FF_CDC = mule.ff.FF_ColumnDependentConstants    # Column

    # Create new file to copy into
    fldfle = mule.FieldsFile()

    # Copy across all standard headers found in an LBC
    fldfle.fixed_length_header = mule.FixedLengthHeader(
        lbc.fixed_length_header.raw[1:])
    fldfle.integer_constants = mule.ff.FF_IntegerConstants(
        lbc.integer_constants.raw[1:])
    fldfle.real_constants = mule.ff.FF_RealConstants(
        lbc.real_constants.raw[1:])
    fldfle.level_dependent_constants = FF_LDC.empty(
        lbc.level_dependent_constants.raw.shape[0])
    fldfle.level_dependent_constants.raw[:, 1:5] = (
        lbc.level_dependent_constants.raw[:, 1:])

    # For variable resolution files both row_dependent_constants (phi_p and
    # phi_v), and column_dependent_constants (lamba_u and lambda_p) will be
    # present and must be copied across.

    # Check Row Dependent constants
    if lbc.row_dependent_constants is not None:
        fldfle.row_dependent_constants = FF_RDC.empty(
            lbc.row_dependent_constants.raw.shape[0]
            )
        fldfle.row_dependent_constants = FF_RDC(
            lbc.row_dependent_constants.raw[:, 1:3]
            )
    # Check Column Dependent Constants
    if lbc.column_dependent_constants is not None:
        fldfle.column_dependent_constants = FF_CDC.empty(
            lbc.column_dependent_constants.raw.shape[0]
            )
        fldfle.column_dependent_constants = FF_CDC(
            lbc.column_dependent_constants.raw[:, 1:3]
            )
    return fldfle


if __name__ == "__main__":
    main()
