import re
import sys

from metomi.rose.upgrade import MacroUpgrade


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


class vn20_t334(MacroUpgrade):
    """Upgrade macro for ticket #334 by Ian Boutle."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t334"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-aerosol
        self.add_setting(config, ["namelist:aerosol", "horiz_d"], "2.25")
        self.add_setting(config, ["namelist:aerosol", "us_am"], "1.45")
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t334"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-aerosol
        # Blank Upgrade Macro
        return config, self.reports
