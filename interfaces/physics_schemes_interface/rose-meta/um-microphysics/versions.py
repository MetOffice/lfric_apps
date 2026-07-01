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


class vn31_t243(MacroUpgrade):
    """Upgrade macro for ticket #243 by Mike Whitall."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t243"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-microphysics
        nml = "namelist:microphysics"
        self.add_setting(config, [nml, "l_improve_precfrac_checks"], ".false.")
        return config, self.reports


class vn31_t360(MacroUpgrade):
    """Upgrade macro for ticket #360 by Ian Boutle."""

    BEFORE_TAG = "vn3.1_t243"
    AFTER_TAG = "vn3.1_t360"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-microphysics
        self.add_setting(config, ["namelist:microphysics", "aut_qc"], "2.47")
        self.add_setting(config, ["namelist:microphysics", "ai"], "2.57e-2")
        upd_precfrac_opt = self.get_setting_value(
            config, ["namelist:microphysics", "i_update_precfrac"]
        )
        self.remove_setting(
            config, ["namelist:microphysics", "i_update_precfrac"]
        )
        self.add_setting(
            config,
            ["namelist:microphysics", "update_precfrac_opt"],
            upd_precfrac_opt,
        )

        return config, self.reports
