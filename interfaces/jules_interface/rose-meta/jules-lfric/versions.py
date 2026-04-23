import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version30_31 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""


class vn31_t400(MacroUpgrade):
    """Upgrade macro for ticket #400 by Dan Copsey."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t400"

    def upgrade(self, config, meta_config=None):

        # Is this a coupled model? The easiest way to get this is from the number of sea ice categories.
        nice = int(self.get_setting_value(config, ["namelist:jules_sea_seaice", "nice"]))
        is_coupled = nice > 1
 
        # Add the new namelist settings
        self.add_setting(config, ["namelist:jules_sea_seaice", "l_zenith_albedo"], ".true.")
        if is_coupled:
            self.add_setting(config, ["namelist:jules_sea_seaice", "meltpond_alb_vn"], "'malinka'")
        else:
            self.add_setting(config, ["namelist:jules_sea_seaice", "meltpond_alb_vn"], "'none'")
        self.add_setting(config, ["namelist:jules_sea_seaice", "snow_grain_size_max"], "100.0")
        self.add_setting(config, ["namelist:jules_sea_seaice", "snow_grain_size_min"], "70.0")
        self.add_setting(config, ["namelist:jules_sea_seaice", "snowpatch"], "0.02")

        return config, self.reports
