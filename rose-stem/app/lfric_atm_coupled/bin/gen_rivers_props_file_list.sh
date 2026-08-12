#! /usr/bin/env bash

# *****************************COPYRIGHT*******************************
# (C) Crown copyright 2024 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

###############################################################################
# This script writes the rivers props file list compatible with the
# version of JULES in which it is present
#
# The path to the data is given by the environment variable
# $RIV_NUMBER_ANCILLARY
#
# The Rivers ancillary file list is written to file_list.txt
# in the current directory
###############################################################################

tee file_list.txt > /dev/null <<EOF
'$RIV_NUMBER_ANCILLARY/qrparm.rivseq.nc'
'$RIV_NUMBER_ANCILLARY/qrparm.rivseq.nc'
'$RIV_NUMBER_ANCILLARY/river_number_trip.nc'
'$RIV_NUMBER_ANCILLARY/qrclim.rivstor.nc'
EOF
