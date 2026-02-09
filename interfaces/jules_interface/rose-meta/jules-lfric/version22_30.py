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


class vn22_t202(MacroUpgrade):
    """Upgrade macro for ticket #202 by Katty Huang."""

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t202"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        self.add_setting(
            config, ["namelist:jules_surface", "anthrop_heat_mean"], "20.0"
        )
        self.add_setting(
            config, ["namelist:jules_surface", "anthrop_heat_option"], "'dukes'"
        )
        return config, self.reports


class vn22_t1012(MacroUpgrade):
    """Upgrade macro for ticket #1012 by Maggie Hendry."""

    BEFORE_TAG = "vn2.2_t202"
    AFTER_TAG = "vn2.2_t1012"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Blank upgrade macro to bump tag
        return config, self.reports


class vn22_t34(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.2_t1012"
    AFTER_TAG = "vn3.0"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Blank Upgrade Macro
        return config, self.reports
