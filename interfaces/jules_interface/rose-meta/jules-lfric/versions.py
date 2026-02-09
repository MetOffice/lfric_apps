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

class vn21_t148(MacroUpgrade):
    # Upgrade macro for #148 by Dan Copsey

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t148"

    def upgrade(self, config, meta_config=None):
        # Add settings
        self.add_setting(config, ["namelist:jules_hydrology", "l_inland"], ".false.")
        return config, self.reports

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
