##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--lfric_diag_file", action="store", required=True,
        help="Input lfric_diagnostics.nc netCDF file"
    )

@pytest.fixture
def diag_infile(request):
    return request.config.getoption("--lfric_diag_file")
