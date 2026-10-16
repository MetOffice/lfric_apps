#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

"""
Create the namcouple file that is used by the Oasis coupler in LFric2LFRic.
The details of the entries are taken from the iodef file and the axis_map.
"""

import xml.etree.ElementTree as ET

NAMCOUPLE_FILE = "new_namcouple"
IODEF_FILE = "iodef_ral3.xml"
PREFIX = "checkpoint"
MAX_FREQ = 5000

AXIS_MAP = {
    "land_tiles":            [8, "weight_file_masked.nc"],
    "snow_layers_and_tiles": [11, "weight_file_nearest.nc"],
    "soil_levels":           [5, "weight_file_masked.nc"],
    "general":               [1, "weight_file_bilinear.nc"],
}


def write_string(f, count, src, dst, freq, weightfile):
    """
    Write the individual mapping entry information for variable being mapped
    from src to dst, with frequency freq and using the specified weightfile.
    """
    f.write("# Mapping entry %s \n" % (count))
    f.write("%s %s 1 %s 1 lfric2lfric_restart_oa.nc EXPORTED \n" % (src, dst, freq))
    f.write("512 1 1024 1 lfric2lfric lfric2lfric SEQ=+1 \n")
    f.write("R 0 R 0 \n")
    f.write("MAPPING \n")
    f.write("%s \n" % (weightfile))
    f.write("########################################################## \n")


def write_header(f, count):
    """
    Write out the initial part of the namcouple file. This includes $NFIELDS (specified
    by count).
    """
    f.write("# Control file for OASIS3-MCT \n")
    f.write("# Oasis reads this file at run-time \n")
    f.write("# \n")
    f.write("# This file was generated automatically by create_namcouple.py \n")
    f.write("# \n")
    f.write("########################################################## \n")
    f.write("# Total number of fields being exchanged = number MAPPING entries below \n")
    f.write(" $NFIELDS \n")
    f.write(" %s \n" % (count))
    f.write(" $END \n")
    f.write("########################################################## \n")
    f.write(" $RUNTIME \n")
    f.write("# This gives you the total simulated time for this run in seconds \n")
    f.write(" 50000 \n")
    f.write(" $END \n")
    f.write("########################################################## \n")
    f.write(" $NLOGPRT \n")
    f.write("# Debugging output level: 0 -> 30 \n")
    f.write("# Timer activation level: 0 -> 3  \n")
    f.write(" 0 0 \n")
    f.write(" $END \n")
    f.write("########################################################## \n")
    f.write(" $STRINGS \n")
    f.write("# Definitions of the particular fields being exchanged \n")
    f.write("# \n")
    f.write("########################################################## \n")


def write_footer(f):
    """
    Write out the last part of the namcouple file. This is required to close the $STRINGS
    that was opened in write_header.
    """
    f.write("$END \n")
    f.write("########################################################## \n")


def process_field(f, field, count, write):
    """
    For given field, determine the field name and axis.
    If this is a special field (a multidata field with an axis) then
    create the src and dst name, get the weights file from AXIS_MAP, and
    then write the entry to the namcouple file.
    """
    var = field.get("id")
    axis = field.get("axis_ref")

    # Only process the fields that start with PREFIX. There may be both input and output
    # fields in the iodef file. We only want to process the output fields.
    if var.startswith(PREFIX):
        if axis is not None:
            try:
                axis_data = AXIS_MAP[axis]
            except KeyError:
                raise ValueError(axis + " is not in AXIS_MAP")
            # axis_data[0] gives the number of levels in the multidata axis
            # axis_data[1] gives the name of the associated weights file
            for n in range(1, axis_data[0]+1):
                src = "src_" + var + "_level_" + str(n)
                dst = "dst_" + var + "_level_" + str(n)
                count += 1
                if (write):
                    write_string(f, count, src, dst, MAX_FREQ, axis_data[1])
    return count


def main():
    """
    Use the iodef file to create a namcouple file. The namcouple file is used by Oasis
    in LFRic2LFRic.
    All non-multidata fields use the 'general' weights file for regridding and reuse the
    same namcouple entry.
    All multidata fields use the specific weights file associated with that axis, and there
    is a separate namcouple entry for each field AND level.
    In lfric2lfric, the multidata fields are all read at the same pseudo time, whilst the
    non-multidata fields (on all levels) are read at successive pseduo times. This is
    implemented here by using a frequency of 1 for the non-multidata fields, and a MAX_FREQ
    for the multidata fields.
    Note: If there are any non-multidata fields that require a different weights file to the
    general file, this code will need to be modified to add an exception list.
    """
    tree = ET.parse(IODEF_FILE)
    root = tree.getroot()

    with open(NAMCOUPLE_FILE, "w", encoding="utf-8") as f:

        # Count how many fields will be written (without writing)
        count = 1
        write = False
        for field in root.iter("field"):
            count = process_field(f, field, count, write)

        # Write the header, including the total count of fields
        write_header(f, count)

        # Write the string for the general regridding (Entry 1)
        count = 1
        try:
            axis_data = AXIS_MAP['general']
        except KeyError:
            raise ValueError(axis + " is not in AXIS_MAP")
        write_string(f, count, "src_face", "dst_face", 1, axis_data[1])

        # Write the strings for the special fields
        write = True
        for field in root.iter("field"):
            count = process_field(f, field, count, write)

        # Write the footer of the file
        write_footer(f)


if __name__ == "__main__":
    main()
