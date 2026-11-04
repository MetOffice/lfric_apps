#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

"""
Create the namcouple file that is used by the Oasis coupler in LFric2LFRic.
The details of the entries are taken from the iodef file. The iodef file is
defined in the global variables, along with the weights files and multidata
axis sizes.
"""

import xml.etree.ElementTree as ET

# ---------------------------------------------------------------
# Global variables - these can be edited for different configurations
# ---------------------------------------------------------------

# Name of the namcouple file (output)
NAMCOUPLE_FILE = "namcouple"

# Name of the iodef file (input)
IODEF_FILE = "iodef_FULL_DUMP_CF.xml"

# Prefix name for the output fields in the iodef file
PREFIX_OUT = "checkpoint"

# Prefix name for the input fields in the iodef file
PREFIX_IN = "restart"

# Arbitrary large number governing the frequency of the coupling
# for the special fields
MAX_FREQ = 5000

# Dictionary for the multidata axes, with the associated number of
# multidata levels, and the associated weights file.
AXIS_MAP = {
    # Fields on land/sea are masked
    "surface_tiles":         [11, "weight_file_masked.nc"],
    "land_tiles":            [8, "weight_file_masked.nc"],
    "sea_ice_categories":    [1, "weight_file_masked.nc"],
    # Multilayer snow fields use nearest neighbour, but NOT masked
    "snow_layers_and_tiles": [27, "weight_file_nn.nc"],
    # Fields on land/sea are masked
    "soil_levels":           [4, "weight_file_masked.nc"],
    "plant_func_types":      [5, "weight_file_masked.nc"],
    "dust_divisions":        [6, "weight_file_masked.nc"],
    # General regridding - e.g. for atmospheric fields
    "general":               [1, "weight_file.nc"],
}

# Dictionary for special fields, with the associated number of
# multidata levels, and the associated weights file. These values
# overide the values in the AXIS_MAP.
SPECIAL_FIELDS = {
    # Snow fields on land use nearest neighbour AND are masked
    "tile_snow_mass":       [11, "weight_file_nn_masked.nc"],
    "tile_snow_grain":      [11, "weight_file_nn_masked.nc"],
    "n_snow_layers":        [11, "weight_file_nn_masked.nc"],
    "snow_depth":           [11, "weight_file_nn_masked.nc"],
    "snowpack_density":     [11, "weight_file_nn_masked.nc"],
    "snow_under_canopy":    [11, "weight_file_nn_masked.nc"],
    # Single layer fields on land are masked
    "z0msea":               [1, "weight_file_masked.nc"],
    "surface_conductance":  [1, "weight_file_masked.nc"],
    "soil_sat_frac":        [1, "weight_file_masked.nc"],
    "wetness_under_soil":   [1, "weight_file_masked.nc"],
    "snow_soot":            [1, "weight_file_masked.nc"],
}

# ---------------------------------------------------------------
# Functions
# ---------------------------------------------------------------

def write_string(namcouple, count, src, dst, freq, weightfile):
    """
    Write the individual mapping entry information for variable being mapped
    from src to dst, with frequency freq and using the specified weightfile.
    """
    namcouple.write(f"# Mapping entry {count} \n")
    namcouple.write(f"{src} {dst} 1 {freq} 1 "
                    "lfric2lfric_restart_oa.nc EXPORTED \n")
    namcouple.write("512 1 1024 1 lfric2lfric lfric2lfric SEQ=+1 \n")
    namcouple.write("R 0 R 0 \n")
    namcouple.write("MAPPING \n")
    namcouple.write(f"{weightfile} \n")
    namcouple.write("############################################ \n")


def write_header(namcouple, count):
    """
    Write out the initial part of the namcouple file. This includes
    $NFIELDS (specified by count).
    """
    namcouple.write("# Control file for OASIS3-MCT \n")
    namcouple.write("# Oasis reads this file at run-time \n")
    namcouple.write("# \n")
    namcouple.write("# This file was generated automatically by "
                    "create_namcouple.py \n")
    namcouple.write("# \n")
    namcouple.write("############################################# \n")
    namcouple.write("# Total number of fields being exchanged = "
                    "number MAPPING entries below \n")
    namcouple.write(" $NFIELDS \n")
    namcouple.write(f"{count} \n")
    namcouple.write(" $END \n")
    namcouple.write("############################################# \n")
    namcouple.write(" $RUNTIME \n")
    namcouple.write("# This gives you the total simulated time for "
                    "this run in seconds \n")
    namcouple.write(" 50000 \n")
    namcouple.write(" $END \n")
    namcouple.write("############################################# \n")
    namcouple.write(" $NLOGPRT \n")
    namcouple.write("# Debugging output level: 0 -> 30 \n")
    namcouple.write("# Timer activation level: 0 -> 3  \n")
    namcouple.write(" 0 0 \n")
    namcouple.write(" $END \n")
    namcouple.write("############################################# \n")
    namcouple.write(" $STRINGS \n")
    namcouple.write("# Definitions of the particular fields being "
                    "exchanged \n")
    namcouple.write("# \n")
    namcouple.write("############################################# \n")


def write_footer(namcouple):
    """
    Write out the last part of the namcouple file. This is required to
    close the $STRINGS that was opened in write_header.
    """
    namcouple.write("$END \n")
    namcouple.write("############################################# \n")


def process_twod_fields(namcouple, field, count):
    """
    For given special field, determine the field name and axis.
    If this is a special field (a 2D field) then
    create the src and dst name, get the weights file from AXIS_MAP
    or SPECIAL_FIELDS, and then write the entry to the namcouple file.
    """
    var = field.get("id")
    domain = field.get("domain_ref")
    axis = field.get("axis_ref")
    input_var = var.replace(PREFIX_OUT + '_', '')

    # 2D fields have domain = face
    if domain == "face":
        # Get the associated multidata levels and weightfile
        # from the axis
        if axis is None:
            axis = 'general'
        try:
            axis_data = AXIS_MAP[axis]
        except KeyError as exc:
            raise ValueError(axis + " is not in AXIS_MAP") from exc

        # Override the associated multidata levels and weightfile
        # if a special field
        axis_data = SPECIAL_FIELDS.get(input_var, axis_data)

        # axis_data[0] gives the number of levels in the multidata axis
        # axis_data[1] gives the name of the associated weights file
        for axis_level in range(1, axis_data[0]+1):
            src = "src_" + input_var + "_level_" + str(axis_level)
            dst = "dst_" + input_var + "_level_" + str(axis_level)
            count += 1

            # Write out the entry to the namcouple
            write_string(namcouple, count, src, dst, MAX_FREQ, axis_data[1])
    return count


def which_fields(root, field_list):
    """
    Determine which fields will be regridded. These fields are in the output
    (denoted by PREFIX_OUT) and the input (denoted by PREFIX_IN).
    """
    for field in root.iter("field"):
        output_var = field.get("id")

        if output_var.startswith(PREFIX_OUT):
            input_var = output_var.replace(PREFIX_OUT, PREFIX_IN)

            if root.findall('.//field[@id="' + input_var + '"]'):
                field_list.append(field)
    return field_list


def count_twod_fields(field_list):
    """
    Count the total number of 2D fields. This includes counting
    the number of levels for multidata fields.
    """
    count = 0
    for field in field_list:
        domain = field.get("domain_ref")
        if domain == "face":
            axis = field.get("axis_ref")
            if axis is None:
                axis = 'general'
            try:
                axis_data = AXIS_MAP[axis]
            except KeyError as exc:
                raise ValueError(axis + " is not in AXIS_MAP") from exc
            count = count + axis_data[0]
    print('Total number of special entries in the namcouple is ', count)
    return count


def main():
    """
    Use the iodef file to create a namcouple file. The namcouple file
    is used by Oasis in LFRic2LFRic.

    All 3D fields use the 'general' weights file for regridding and
    reuse the same namcouple entry.
    All 2D fields use the specific weights file associated with that
    field or multidata axis, and there is a separate namcouple entry
    for each field AND multidata-level.

    In lfric2lfric, the 2D fields are all interpolated at the same
    pseudo time, whilst the 3D fields are interpolated at successive
    pseudo times. This is implemented here by using a frequency of 1
    for the 3D fields, and a MAX_FREQ for the 2D fields.
    """
    tree = ET.parse(IODEF_FILE)
    root = tree.getroot()

    with open(NAMCOUPLE_FILE, "w", encoding="utf-8") as namcouple:

        # Determine which fields will be regridded
        field_list = []
        field_list = which_fields(root, field_list)

        # The total number of entries needs to be written in the
        # namcouple header. So first count the fields that need a
        # separate entry in the namcouple

        # Count the number of 2D fields (with separate multidata levels)
        # to be regridded
        total_count = count_twod_fields(field_list)

        # Add 1 to the count (entry for 3D fields)
        total_count = total_count + 1

        # Write the header, including the total count of entries
        write_header(namcouple, total_count)

        # Write the string for the 3D fields (Entry 1)
        # using frequency of 1, and weight_file defined in AXIS_MAP.
        count = 1
        frequency = 1
        axis = 'general'
        try:
            axis_data = AXIS_MAP[axis]
        except KeyError as exc:
            raise ValueError(axis + " is not in AXIS_MAP") from exc
        write_string(namcouple, count, "src_face", "dst_face",
                     frequency, axis_data[1])

        # Write the strings for the special fields (Entries 2 to total_count)
        for field in field_list:
            count = process_twod_fields(namcouple, field, count)

        # Write the footer of the file
        write_footer(namcouple)


if __name__ == "__main__":
    main()
