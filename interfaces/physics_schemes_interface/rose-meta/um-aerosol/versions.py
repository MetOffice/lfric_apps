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


class vn32_t698(MacroUpgrade):
    """Upgrade macro for ticket #698 by Alan J Hewitt."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t698"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-aerosol
        # Add new settings with the default option SUBCOCSSDU_7mode
        self.add_setting(
            config, ["namelist:aerosol", "mode_setup"], "'SUBCOCSSDU_7mode'"
        )
        # Default to false since this is the setting in all existing tests
        self.add_setting(
            config, ["namelist:aerosol", "l_dust_mp_ageing"], ".false."
        )
        # Default to true since this is the setting in all existing tests
        self.add_setting(
            config, ["namelist:aerosol", "l_ukca_radaer_sustrat"], ".true."
        )

        return config, self.reports
