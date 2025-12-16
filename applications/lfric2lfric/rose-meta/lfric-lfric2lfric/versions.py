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

class vn30_t48(MacroUpgrade):
    """Upgrade macro for ticket #48 by Juan M Castillo."""

    BEFORE_TAG = "vn3.0"
    AFTER_TAG = "vn3.0_t48"

    def upgrade(self, config, meta_config=None):
        self.add_setting(config, ["namelist:lfric2lfric", "mode"], "'ics'")
        self.add_setting(
            config,
            ["namelist:lfric2lfric", "source_file_lbc"],
            "'source_file_lbc'",
        )
        self.add_setting(
            config,
            ["namelist:lfric2lfric", "weight_file_lbc"],
            "'weight_file_lbc'",
        )

        return config, self.reports
