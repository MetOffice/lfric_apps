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
        Validate that the output file has a variable named 'Mesh2d_node_x'
        """
        assert load_rootgrp.variables.get("Mesh2d_node_x").standard_name == 'longitude'

    def test_Mesh2d_node_x_units(self, load_rootgrp):
        """
        Validate that the output file has a variable named 'Mesh2d_node_x'
        """
        assert load_rootgrp.variables.get("Mesh2d_node_x").units == 'degrees_east'

    def test_Mesh2d_node_y_exists(self, load_rootgrp):
        """
        Validate that the output file has a variable named 'Mesh2d_node_y'
        """
        assert "Mesh2d_node_y" in load_rootgrp.variables

    def test_Mesh2d_node_y_stdname(self, load_rootgrp):
        """
        Validate that the output file has a variable named 'Mesh2d_node_y'
        """
        assert load_rootgrp.variables.get("Mesh2d_node_y").standard_name == 'latitude'

    def test_Mesh2d_node_y_units(self, load_rootgrp):
        """
        Validate that the output file has a variable named 'Mesh2d_node_y'
        """
        assert load_rootgrp.variables.get("Mesh2d_node_y").units == 'degrees_north'
