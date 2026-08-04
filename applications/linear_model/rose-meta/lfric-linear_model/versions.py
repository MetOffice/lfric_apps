import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version31_32 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__



class vn32_t683(MacroUpgrade):
    """
    Upgrade macro for ticket #683 by Alan J Hewitt.
    Users can now select GLOMAP setting via namelist.
    NWP option glomap_mode_dust_and_clim is hard coded to i_mode_setup == 8
    Simple option glomap_mode_climatology is hard coded to i_mode_setup == 8
    """

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t683"

    def upgrade(self, config, meta_config=None):
        # Add settings
        self.add_setting( config, ["namelist:aerosol",
                                   "i_mode_setup"], "8" )
        
        return config, self.reports

