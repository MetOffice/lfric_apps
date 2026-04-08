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


class vn31_t360(MacroUpgrade):
    # Upgrade macro for #360 by Ian Boutle

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t360"

    def upgrade(self, config, meta_config=None):
        # Add settings
        self.add_setting(config, ["namelist:microphysics","aut_qc"],"2.47")
        self.add_setting(config, ["namelist:microphysics","ai"],"2.57e-2")
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
