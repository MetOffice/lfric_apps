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


class vn30_t192(MacroUpgrade):
    """Upgrade macro for ticket #192 by Juan M. Castillo."""

    BEFORE_TAG = "vn3.0_t99"
    AFTER_TAG = "vn3.0_t192"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-lfric2lfric
        self.rename_setting(
            config,
            ["namelist:partitioning(source)", "mesh_type"],
            ["namelist:partitioning(source)", "mesh_target"],
        )
        self.rename_setting(
            config,
            ["namelist:partitioning(destination)", "mesh_type"],
            ["namelist:partitioning(destination)", "mesh_target"],
        )

        return config, self.reports
