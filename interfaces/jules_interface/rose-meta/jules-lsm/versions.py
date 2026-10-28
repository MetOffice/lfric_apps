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


class vn32_tXXX(MacroUpgrade):
    """Upgrade macro for ticket #XXX by Ian Boutle."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_tXXX"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lsm
        """Add the initial water table and wetness under the soil column"""
        self.add_setting(
            config, ["namelist:ideal_surface", "water_table"], "1.0"
        )
        self.add_setting(
            config, ["namelist:ideal_surface", "wetness_under_soil"], "1.0"
        )
        return config, self.reports
