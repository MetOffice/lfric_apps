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
import numpy as np
from operator import itemgetter

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
    sliced to remove the halo regions. All fields created in this function
    are written to an outfile of the same name as the LBC but with the `.ff`
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

    with open(input_filename + ".ff", "wb+") as outfile:

        # Counter for number of fields processed
        num_fields = 0

        # Required for navigating open file (will be set later)
        data_loc = 0
        lookups_loc = 0

        # The Fixed Length Header
        header = fldfle.fixed_length_header

        # Navigate to the COMPONENTS section of the open file
        outfile.seek(header._NUM_WORDS * fldfle.WORD_SIZE)
        fldfle._write_singular_headers(outfile)

        # Record where we left off, we will write LOOKUPs here
        single_headers_end = (outfile.tell() // fldfle.WORD_SIZE)
        # Store this in our FLH
        header.lookup_start = single_headers_end + 1

        # First we need to calculate how many fields we will generate.
        # This is crucial for setting fixed_length_headers values for where
        # data start.
        for field in lbc.fields:
            num_fields += 1                         # each main field is used
            num_values = field.num_values()         # record n_vals
            levels = enumerate(field.get_data())    # number each level in data
            num_fields += max(levels, key=itemgetter(0))[0]     # add max lev

        # calculate where to write the first field's data
        word_number = ((single_headers_end + 1) + num_values *
                        num_fields)
        offset = word_number - 1
        offset -= offset % -fldfle._DATA_START_ALIGNMENT
        header.data_start = offset + 1

        # navigate to that location
        outfile.seek((header.data_start - 1) * fldfle.WORD_SIZE)

        # Now we can loop over and write the data to file
        # first, reset num_fields
        num_fields = 0

        # Loop over the fields in the LBC, unpack, write to file
        for field in lbc.fields:
            field = lbc_to_masked(field)
            ncols = field.lbnpt
            nrows = field.lbrow
            # Rim and halo widths
            halo_code = field.lbuser3
            rimwidth = int(halo_code // 10000)
            halo_ns = int(halo_code - rimwidth * 10000) // 100
            halo_ew = int(halo_code - rimwidth * 10000 - halo_ns * 100)

            # Update the stash code to the standard prognostic version
            field.lbuser4 = LBC_STASH_MAP[field.lbuser4]
            # Update lbhem variable to be what is expected by fieldsfile
            field.lbhem = fldfle.fixed_length_header.horiz_grid_type % 100

            for level_num, single_level in enumerate(field.get_data(), 1):
                # Copy field to get metadata
                field_2d = field.copy()
                # Counter for number of fields processed - this lives in here
                # because we are creating a new fieldsfile field for every
                # level in each LBC field
                num_fields += 1
                # Data provider is mule method for telling new field
                # where to get data from. Slice the 1D-array to remove halo
                # regions as not needed in LFRic LBC file
                array_provider = mule.ArrayDataProvider(
                    single_level.filled(mule._REAL_MDI)[
                        halo_ns:nrows + halo_ns,
                        halo_ew:ncols + halo_ew])
                field_2d.set_data_provider(array_provider)

                # Update level number
                field_2d.lblev = level_num

                # Set lookup dimensions in our header
                header.lookup_dim1 = field_2d.num_values()
                header.lookup_dim2 = num_fields

                # Write to file
                # Both provides and returns data and lookup locations because
                # they change after each iteration
                data_loc, lookups_loc = write_field(outfile, field_2d, fldfle,
                                                    num_fields,
                                                    data_loc, lookups_loc)

        # Return to the beginning to write out fixed_length_header
        outfile.seek(0)
        # Set dataset version
        fldfle.fixed_length_header.data_set_format_version = 20
        # Update dataset type to be fieldsfile
        fldfle.fixed_length_header.dataset_type = 3
        fldfle.fixed_length_header.to_file(outfile)


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


def write_field(outfile, field, fieldsfile, num_fields, data_loc, lookups_loc):
    '''
    This function writes a single :py:class:`mule.Field3` instance to an open
    binary file.

    This function is a modified version of the `write_to_file` method of the
    :py:class:`mule.UMFile` class. The original method is designed to work on
    a list of fields owned by the :py:class:`mule.UMFile` instance.

    == Parameters
    param outfile: An open binary file in write mode.
    type lbc: :py:class:`_io.BufferedRandom`

    param field: The field to be written to the open file
    type field: :py:class:`mule.Field3`

    param fieldsfile: The fieldsfile instance that contains the appropriate
                      fixed length header
    type fieldsfile: :py:class:`mule.ff.FieldsFile`

    param int data_loc: The location that this field's data will be written

    param in lookups_loc: The location that this field's LOOKUPs will be written

    == Returns
    return int data_loc: The location in the binary file that is next available
                         (unwritten) in the file's data array

    return int lookup_loc: The location in the binary file that is next
                           available (unwritten) in the file's LOOKUP array

    '''
    # A reference to the header
    flh = fieldsfile.fixed_length_header

    # sense check for invalid release number
    if field.lbrel == -99.0:
        return

    # =========================================================
    # 0.0 Record where the this field's data start
    # =========================================================
    # If first field, data starts at current location
    if num_fields == 1:
        data_loc = outfile.tell()

    # record where field starts and navigate to that location
    field.lbegin = data_loc // fieldsfile.WORD_SIZE
    outfile.seek(data_loc)

    # =====================================================
    # 1.0 Check and set packing code
    # =====================================================
    # WGDOS packed fields can be tagged with an accuracy of
    # -99.0; this indicates that they should not be packed,
    # so reset the packing code here accordingly
    if field.lbpack % 10 == 1 and int(field.bacc) == -99:
        field.lbpack = 10*(field.lbpack//10)

    # =====================================================
    # 2.0 Get and write data
    # =====================================================
    # Strip just the n1-n3 digits from the lbpack value
    # since the later digits are not relevant
    lbpack321 = "{0:03d}".format(
                        field.lbpack - ((field.lbpack//1000) % 10)*1000)

    # Select an appropriate operator for writing the data
    # (if one is available for the given packing code)
    if lbpack321 in fieldsfile.WRITE_OPERATORS:
        write_operator = fieldsfile._write_operators[lbpack321]
    else:
        raise ValueError(f'Cannot save data with lbpack={field.lbpack} : '
                'packing not supported.')

    # Use the write operator to prepare the field data for
    # writing to disk
    data_bytes, data_size = write_operator.to_bytes(field)

    # # The bytes returned by the operator are in the exact
    # # format to be written
    outfile.write(data_bytes)
    # and the operator also returns the exact number of
    # words/records taken up by the data; this is exactly
    # what needs to go in the Field's lblrec
    field.lblrec = data_size

    # The other record header, lbnrec, is the number of
    # words/records used to store the data; this may be
    # different to the above in the case of packed data;
    # if the packing method has a different word size.
    # Calculate the actual on-disk word size here
    size_on_disk = ((write_operator.WORD_SIZE*data_size) //
                    fieldsfile.WORD_SIZE)

    # =========================================================
    # 2.1 Add any padding to create a uniform data section
    # =========================================================
    # Padding will also be applied to ensure that the next
    # block of data is aligned with a sector boundary
    field.lbnrec = (
        size_on_disk -
        (size_on_disk % -fieldsfile._WORDS_PER_SECTOR))

    # Check how big each sector should be (used to check if padding needed)
    sector_size = fieldsfile._WORDS_PER_SECTOR * fieldsfile.WORD_SIZE

    # Pad out the data section to a whole number of sectors.
    overrun = outfile.tell() % sector_size
    if overrun != 0:
        padding = np.zeros(sector_size - overrun, 'i1')
        outfile.write(padding)

    # =========================================================
    # 2.2 Update data location
    # =========================================================
    # The current position is returned and can be used on a subsequent
    # function call to locate the next unused location after the last
    # data sector
    data_loc = outfile.tell()

    # =========================================================
    # 3.0 Write LOOKUP component for field
    # =========================================================
    # update fixed_length header with latest data size
    flh.data_dim1 = ((outfile.tell() // fieldsfile.WORD_SIZE) -
                      flh.data_start + 1
                      )
    # Go back and write the LOOKUP component.
    if num_fields == 1:
        outfile.seek((flh.lookup_start - 1) * fieldsfile.WORD_SIZE)
    else:
        outfile.seek(lookups_loc)
    # write field lookups to file
    field.to_file(outfile)

    # Store where this lookup ended (actually this minus 1)
    lookups_loc = outfile.tell()

    return data_loc, lookups_loc


if __name__ == "__main__":
    main()
