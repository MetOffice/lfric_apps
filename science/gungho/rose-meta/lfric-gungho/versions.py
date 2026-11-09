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

class vn32_tNNN(MacroUpgrade):
    """Upgrade macro for ticket #NNN by Thomas Bendall."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_tNNN"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Replace horizontal_transport_predictor with wind_transport_predictor."""
        old = ["namelist:formulation", "horizontal_transport_predictor"]
        new = ["namelist:formulation", "wind_transport_predictor"]
        value = self.get_setting_value(config, old)
        self.remove_setting(config, old)
        if value == ".true.":
            self.add_setting(config, new, "'no_etadot'")
        else:
            self.add_setting(config, new, "'3d'")
        return config, self.reports


