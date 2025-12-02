import re
import sys

from metomi.rose.upgrade import MacroUpgrade

from .version20_21 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn21_t756(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t756"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/socrates-radiation
        self.add_setting(
            config, ["namelist:radiative_gases", "o3_profile_data"], "0"
        )
        self.add_setting(
            config, ["namelist:radiative_gases", "o3_profile_heights"], "0.0"
        )
        self.add_setting(
            config, ["namelist:radiative_gases", "o3_profile_size"], "0"
        )
        return config, self.reports


class vn21_t871(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.1_t756"
    AFTER_TAG = "vn2.2"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/socrates-radiation
        # Blank Upgrade Macro
        return config, self.reports
