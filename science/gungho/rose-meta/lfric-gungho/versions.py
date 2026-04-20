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
class vn31_t443(MacroUpgrade):
    # Upgrade macro for PR #443 by Samantha Pullen

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t443"

    def upgrade(self, config, meta_config=None):
        # Add name entry to iau_addinf_io namelist
        self.add_setting(
            config, ["namelist:iau_addinf_io(addinf1)", "name"], "''"
        )
        self.add_setting(
            config, ["namelist:iau_addinf_io(addinf2)", "name"], "''"
        )
        # Add name entry to iau_ainc_io namelist
        self.add_setting(
            config, ["namelist:iau_ainc_io(ainc1)", "name"], "''"
        )
        self.add_setting(
            config, ["namelist:iau_ainc_io(ainc2)", "name"], "''"
        )
        # Add name entry to iau_bcorr_io namelist
        self.add_setting(
            config, ["namelist:iau_bcorr_io(bcorr1)", "name"], "''"
        )
        return config, self.reports