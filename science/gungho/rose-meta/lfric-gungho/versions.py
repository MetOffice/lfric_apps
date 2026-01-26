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


class vn30_t238(MacroUpgrade):
    """Upgrade macro for ticket #238 by Thomas Bendall."""

    BEFORE_TAG = "vn3.0"
    AFTER_TAG = "vn3.0_t238"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(
            config, ["namelist:finite_element", "coord_space"], "'Wchi'"
        )
        coord_order = self.get_setting_value(
            config, ["namelist:finite_element", "coord_order"]
        )
        self.add_setting(
            config,
            ["namelist:finite_element", "coord_order_multigrid"],
            coord_order,
        )

        return config, self.reports


class vn30_t180(MacroUpgrade):
    """Upgrade macro for ticket #180 by Thomas Bendall."""

    BEFORE_TAG = "vn3.0_t238"
    AFTER_TAG = "vn3.0_t180"

    def upgrade(self, config, meta_config=None):
        # Get values
        n_orog_smooth = self.get_setting_value(
            config, ["namelist:initialization", "n_orog_smooth"]
        )
        w0_mapping = self.get_setting_value(
            config, ["namelist:initialization", "w0_orography_mapping"]
        )
        coord_order = self.get_setting_value(
            config, ["namelist:finite_element", "coord_order"]
        )

        # Add new settings
        self.add_setting(
            config, ["namelist:orography", "n_orog_smooth"], n_orog_smooth
        )
        self.add_setting(
            config, ["namelist:orography", "w0_multigrid_mapping"], w0_mapping
        )
        self.add_setting(
            config, ["namelist:orography", "orography_order"], coord_order
        )

        # Remove old settings
        self.remove_setting(
            config, ["namelist:initialization", "n_orog_smooth"]
        )
        self.remove_setting(
            config, ["namelist:initialization", "w0_orography_mapping"]
        )

        return config, self.reports
