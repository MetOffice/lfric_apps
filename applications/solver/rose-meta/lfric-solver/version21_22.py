import re
import sys

from metomi.rose.upgrade import MacroUpgrade

from .version20_21 import *


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


class vn21_t708(MacroUpgrade):
    """Upgrade macro for ticket #708 by mark Hedley."""

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t708"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-solver
        """Add new log reduction switch to namelist logging"""
        self.add_setting(
            config, ["namelist:logging", "log_to_rank_zero_only"], ".false."
        )
        return config, self.reports


class vn21_t4604(MacroUpgrade):
    """Upgrade macro for ticket #4604 by Mike Hobson."""

    BEFORE_TAG = "vn2.1_t708"
    AFTER_TAG = "vn2.1_t4604"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-solver
        """Rename generate_inner_haloes to generate_inner_halos"""
        self.rename_setting(
            config,
            ["namelist:partitioning", "generate_inner_haloes"],
            ["namelist:partitioning", "generate_inner_halos"],
        )
        return config, self.reports


class vn21_t871(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.1_t4604"
    AFTER_TAG = "vn2.2"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-solver
        # Blank Upgrade Macro
        return config, self.reports
