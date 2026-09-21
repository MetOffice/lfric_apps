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


class vn32_t802(MacroUpgrade):
    """ Upgrade macro for #802 by thomas.melvin"""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t802"

    def upgrade(self, config, meta_config=None):
        # Add lmin and lmax settings for the mixed solver
        self.add_setting(
            config,
            ["namelist:mixed_solver", "chebyshev_lmax"],
            "2.5",
        )
        self.add_setting(
            config,
            ["namelist:mixed_solver", "chebyshev_lmin"],
            "1.0",
        )
        return config, self.reports
