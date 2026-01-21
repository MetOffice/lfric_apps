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

class vn30_t171(MacroUpgrade):
    # Upgrade macro for #171 by James Kent

    BEFORE_TAG = "vn3.0"
    AFTER_TAG = "vn3.0_t171"

    def upgrade(self, config, meta_config=None):
        # Add adjust_tracer_equation to transport namelist
        self.add_setting(
            config, ["namelist:transport", "adjust_tracer_equation"], ".false."
        )

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
