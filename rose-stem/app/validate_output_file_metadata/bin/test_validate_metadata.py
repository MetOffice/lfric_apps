##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Module to analyse metadata structures in output files from LFRic.
Rules for metadata output are encoded.

"""

import netCDF4
import pytest
@pytest.fixture
def load_rootgrp(infile):
    rootgrp = netCDF4.Dataset(infile, format="NETCDF4")
    print(f"loading ... ... {infile}")
    return rootgrp

# def test_time_instant_exists(load_rootgrp):
#     "validate that the output file has a variable named 'time_instant'"
#     assert "time_instant" in load_rootgrp.variables

def test_time_centered_exists(load_rootgrp):
    "validate that the output file has a variable named 'time_instant'"
    assert "time_centered" in load_rootgrp.variables

def test_time_counter_not_exists(load_rootgrp):
    "validate that the output file does not have a variable named 'time_counter'"
    assert "time_counter" not in load_rootgrp.variables
