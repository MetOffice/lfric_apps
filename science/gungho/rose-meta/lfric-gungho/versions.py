import sys
import re

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


class vn32_t800(MacroUpgrade):
    """Upgrade macro for PR #800 by Chris Smith."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t800"

    def upgrade(self, config, meta_config=None):
        """Add vapour_relax namelist to configuration source list"""
        source = self.get_setting_value(config, ["file:configuration.nml","source"])
        source = re.sub(r'\(namelist:vapour_forcing\)',
                        r'(namelist:vapour_forcing)' + '\n' + ' (namelist:vapour_relax)',
                        source)
        self.change_setting_value(config, ["file:configuration.nml","source"], source)
        """Add vapour_relaxation setting to external_forcing namelist"""
        self.add_setting(config, ["namelist:external_forcing", "vapour_relaxation"], ".false.")
        """Data for vapour_relax namelist"""
        self.add_setting(config, ["namelist:vapour_relax"])
        self.add_setting(config, ["namelist:vapour_relax", "coordinate"], "'height'")
        self.add_setting(config, ["namelist:vapour_relax", "heights"], "0.0")
        self.add_setting(config, ["namelist:vapour_relax", "number_heights"], "1")
        self.add_setting(config, ["namelist:vapour_relax", "number_times"], "1")
        self.add_setting(config, ["namelist:vapour_relax", "profile_data"], "0.0")
        self.add_setting(config, ["namelist:vapour_relax", "times"], "0.0")
        self.add_setting(config, ["namelist:vapour_relax", "timescale"], "1.0")
        self.add_setting(config, ["namelist:vapour_relax", "variable"], "'mr'")
        return config, self.reports
