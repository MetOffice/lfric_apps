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


class vn31_t249(MacroUpgrade):
    # Upgrade macro for #249 by Mike Whitall

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t249"

    def upgrade(self, config, meta_config=None):
        # Add settings

        # Blank macro needed just to update meta-data version
        # (apps using the new option 'smooth_fix' under the existing
        #  multi-option switch 'pc2_init_logic' fail checks against
        #  the existing meta-data).

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
