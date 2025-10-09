import sys

from metomi.rose.upgrade import MacroUpgrade

from .version21_22 import *


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


class vn22_t887(MacroUpgrade):
    """Upgrade macro for ticket #887 by Mike Whitall."""

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t887"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        nml = "namelist:cloud"
        self.add_setting(config, [nml, "dbsdtbs_turb_0"], "1.5E-4")
        self.add_setting(config, [nml, "i_pc2_erosion_numerics"], "'implicit'")

        return config, self.reports
