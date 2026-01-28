import sys

from metomi.rose.upgrade import MacroUpgrade

from .version22_30 import *


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


class vn30_t99(MacroUpgrade):
    """Upgrade macro for ticket #99 by Fred Wobus."""

    BEFORE_TAG = "vn3.0"
    AFTER_TAG = "vn3.0_t99"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-lfric_atm
        """Set segmentation size for Gregory-Rowntree convection kernel"""
        self.add_setting(config, ["namelist:physics", "conv_gr_segment"], "16")
        return config, self.reports


class vn30_t182(MacroUpgrade):
    # Upgrade macro for #182 by Tom Hill

    BEFORE_TAG = "vn3.0_t99"
    AFTER_TAG = "vn3.0_t182"

    def upgrade(self, config, meta_config=None):
        """Add linear boundary layer physics scheme"""
        self.add_setting(config, ["namelist:linear", "Blevs_m"], "15")
        self.add_setting(config, ["namelist:linear", "e_folding_levs_m"], "10")
        self.add_setting(config, ["namelist:linear", "l_0_m"], "80.0")
        self.add_setting(config, ["namelist:linear", "l_boundary_layer"], ".true.")
        self.add_setting(config, ["namelist:linear", "log_layer"], "2")
        self.add_setting(config, ["namelist:linear", "u_land_m"], "0.4")
        self.add_setting(config, ["namelist:linear", "u_sea_m"], "0.4")
        self.add_setting(config, ["namelist:linear", "z_land_m"], "0.05")
        self.add_setting(config, ["namelist:linear", "z_sea_m"], "0.0005")
        return config, self.reports
