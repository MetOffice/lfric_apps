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
class vn31_t504(MacroUpgrade):
    """Upgrade macro for 504 by Adrian Lock"""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t504"

    def upgrade(self, config, meta_config=None):
        # Add settings
        self.add_setting(config,
                  ["namelist:initial_temperature", "theta_pert_start"],
                   "0.0")
        self.add_setting(config,
                  ["namelist:initial_temperature", "theta_pert_end"],
                   "0.0")
        self.add_setting(config,
                  ["namelist:initial_temperature", "theta_pert_size"],
                   "0.0")
        return config, self.reports
