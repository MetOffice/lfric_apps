714import sys

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


class vn32_t379(MacroUpgrade):
    """Upgrade macro for ticket #714 by Thomas Bendall."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t714"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Semi-Implicit setting changes ----------------------------------------
        # Get number of outer/inner iterations from the existing namelist
        outer = self.get_setting_value(
            config, ["namelist:timestepping", "outer_iterations"]
        )
        inner = self.get_setting_value(
            config, ["namelist:timestepping", "inner_iterations"]
        )
        # Determine inner iterations settings
        if outer in [None, "''", ""]:
            inner_array = "''"
        else:
            inner_array = f"{inner}"
            for _ in range(int(outer) - 1):
                inner_array += f",{inner}"
        self.add_setting(
            config,
            ["namelist:timestepping", "inner_iterations_si"],
            inner_array,
        )
        # Remove old setting
        self.remove_setting(
            config, ["namelist:timestepping", "inner_iterations"]
        )
        # TR-BDF2 settings -----------------------------------------------------
        # This is easier because it's new and we can prescribe what to do
        self.add_setting(
            config, ["namelist:timestepping", "inner_iterations_tr"], "2,1"
        )
        self.add_setting(
            config, ["namelist:timestepping", "inner_iterations_bdf2"], "1,1"
        )

        return config, self.reports
