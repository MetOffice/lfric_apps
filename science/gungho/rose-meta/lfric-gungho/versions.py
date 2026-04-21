import re
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


class vn31_t118(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t118"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Blank Upgrade Macro
        return config, self.reports


class vn31_t363(MacroUpgrade):
    """Upgrade macro for ticket #363 by Jaffery Irudayasamy."""

    BEFORE_TAG = "vn3.1_t118"
    AFTER_TAG = "vn3.1_t363"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Set segmentation size limit for short and long wave radiation kernels"""
        self.add_setting(config, ["namelist:physics", "sw_segment_limit"], "32")
        self.add_setting(config, ["namelist:physics", "lw_segment_limit"], "32")

        return config, self.reports
