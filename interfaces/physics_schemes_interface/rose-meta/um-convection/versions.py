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

class vn31_t368(MacroUpgrade):
    """Upgrade macro for ticket #368 by Ian Boutle."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t368"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-convection
        self.add_setting(
            config, ["namelist:convection", "llcs_first_outer"], ".false."
        )

        return config, self.reports
