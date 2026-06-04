##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Module to analyse metadata structures in lfric_diagnostics output files
from LFRic.
Rules for metadata output are encoded.

"""

import netCDF4
import os
import pytest

@pytest.fixture
def load_rootgrp(diag_infile):
    rootgrp = None
    if os.path.exists(diag_infile):
        rootgrp = netCDF4.Dataset(diag_infile, format="NETCDF4")

    return rootgrp

class TestTemporalMetadata:
    def test_time_centered_exists(self, load_rootgrp):
        """
        Validate that the output file has a variable named 'time_instant'
        """
        assert "time_centered" in load_rootgrp.variables

    def test_time_counter_not_exists(self, load_rootgrp):
        """
        Validate that the output file does not have a variable named 'time_counter'
        """
        assert "time_counter" not in load_rootgrp.variables

class TestSpatialMetadata:
    def test_Mesh2d_node_x_exists(self, load_rootgrp):
        """
        Validate that the output file has a variable named 'Mesh2d_node_x'
        """
        assert "Mesh2d_node_x" in load_rootgrp.variables

    def test_Mesh2d_node_x_stdname(self, load_rootgrp):
        """
        Validate that the output file x variable has longitude standard_name
        """
        assert load_rootgrp.variables.get("Mesh2d_node_x").standard_name == 'longitude'

    def test_Mesh2d_node_x_units(self, load_rootgrp):
        """
        Validate that the output file x variable has degrees_east units
        """
        assert load_rootgrp.variables.get("Mesh2d_node_x").units == 'degrees_east'

    def test_Mesh2d_node_y_exists(self, load_rootgrp):
        """
        Validate that the output file has a variable named 'Mesh2d_node_y'
        """
        assert "Mesh2d_node_y" in load_rootgrp.variables

    def test_Mesh2d_node_y_stdname(self, load_rootgrp):
        """
        Validate that the output file y variable has latitude standard_name
        """
        assert load_rootgrp.variables.get("Mesh2d_node_y").standard_name == 'latitude'

    def test_Mesh2d_node_y_units(self, load_rootgrp):
        """
        Validate that the output file y variable has degrees_north units
        """
        assert load_rootgrp.variables.get("Mesh2d_node_y").units == 'degrees_north'

class TestPackingMetadata:
    """
    We use the temperature value as our reference, as this is a floating point
    value that should be packed to a 16 bit integer in our packed case
    """
    def test_temp_variable_unpacked(self, load_rootgrp):
        """
        Test that if the variable is unpacked, add_offset and scale_factor
        are not present in the metadata. Packed dtype is int16
        """
        temp = load_rootgrp.variables.get('temperature')
        if temp.dtype != 'int16':
            assert not hasattr(temp, 'add_offset')
            assert not hasattr(temp, 'scale_factor')

    def test_temp_variable_packed(self, load_rootgrp):
        """
        Test that if the variable is packed, add_offset and scale_factor
        are present in the metadata. Packed dtype is int16. We won't check
        the actual values as these may change between configurations
        """
        temp = load_rootgrp.variables.get('temperature')
        if temp.dtype == 'int16':
            assert hasattr(temp, 'add_offset')
            assert hasattr(temp, 'scale_factor')
