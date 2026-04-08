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
        self.add_setting(config, ["namelist:aerosol","ukca_scale_marine_pom_ems"],".false.")
        self.add_setting(config, ["namelist:aerosol","marine_pom_ems_scaling"],"1.0")
        self.add_setting(config, ["namelist:aerosol","ukca_scale_sea_salt_ems"],".false.")
        self.add_setting(config, ["namelist:aerosol","sea_salt_ems_scaling"],"1.0")
        return config, self.reports

