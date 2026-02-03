import re
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


class vn30_t135(MacroUpgrade):
    """Upgrade macro for ticket #135 by James Manners."""

    BEFORE_TAG = "vn3.0"
    AFTER_TAG = "vn3.0_t135"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/socrates-radiation
        self.add_setting(config, ["namelist:cosp", "n_cosp_step"], "1")

        return config, self.reports
