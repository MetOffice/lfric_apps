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
class vn32_t634(MacroUpgrade):
    # Upgrade macro for #634 by Ian Boutle

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t634"

    def upgrade(self, config, meta_config=None):
        # Add settings
        nml="namelist:boundaries"
        self.add_setting(config,[nml,"lbc_bal_meth"], "'keep_rho'")
        self.add_setting(config,[nml,"lbc_sort_theta"], ".true.")
        nml="namelist:initialization"
        eos_height=self.get_setting_value(config,[nml, "model_eos_height"])
        self.remove_setting(config,[nml,"model_eos_height"])
        self.add_setting(config,[nml,"init_eos_height"], eos_height)
        self.add_setting(config,[nml,"init_exner_method"], "'hydrostatic'")
        self.add_setting(config,[nml,"init_sort_theta"], ".true.")
        return config, self.reports
