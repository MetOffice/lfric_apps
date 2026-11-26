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


class vn32_txxx(MacroUpgrade):
    """Upgrade macro for ticket #xxx by Thomas Bendall."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_txxx"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Convert theta_dispersion_correction from logical to enumeration"""
        section = ["namelist:transport", "theta_dispersion_correction"]
        value = self.get_setting_value(config, section)
        if value == ".false.":
            self.change_setting_value(config, section, "'none'")
        elif value == ".true.":
            self.change_setting_value(config, section, "'original'")
        return config, self.reports
