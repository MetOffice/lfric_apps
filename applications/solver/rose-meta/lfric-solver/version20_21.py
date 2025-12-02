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


class vn20_t358(MacroUpgrade):
    """Upgrade macro for ticket #358 by Joshua Dendy."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t358"

    def upgrade(self, config, meta_config=None):
        # Commands From: applications/solver/rose-meta/lfric-solver
        """
        Add element_order_h and element_order_v to namelist finite_element,
        replacing element_order
        """
        self.add_setting(
            config, ["namelist:finite_element", "element_order_h"], "0"
        )
        self.add_setting(
            config, ["namelist:finite_element", "element_order_v"], "0"
        )
        self.remove_setting(
            config, ["namelist:finite_element", "element_order"]
        )
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t358"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-solver
        # Blank Upgrade Macro
        return config, self.reports
