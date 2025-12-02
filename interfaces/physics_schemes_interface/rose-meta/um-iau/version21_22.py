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


class vn21_t207(MacroUpgrade):
    """Upgrade macro for ticket #207 by Mike Thurlow."""

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t207"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-iau
        """Add new variables and default data to iau namelist"""
        self.add_setting(
            config, ["namelist:iau", "iau_ainc_multifile"], ".false."
        )
        self.add_setting(
            config, ["namelist:iau", "iau_tendency_addinf"], ".true."
        )
        self.add_setting(
            config, ["namelist:iau", "iau_tendency_ainc"], ".false."
        )
        self.add_setting(
            config, ["namelist:iau", "iau_tendency_bcorr"], ".true."
        )
        self.add_setting(
            config, ["namelist:iau", "iau_tendency_pertinc"], ".false."
        )
        self.add_setting(config, ["namelist:iau", "iau_use_addinf"], ".false.")
        self.add_setting(config, ["namelist:iau", "iau_use_bcorr"], ".false.")
        self.add_setting(config, ["namelist:iau", "iau_use_pertinc"], ".false.")
        return config, self.reports


class vn21_t871(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.1_t207"
    AFTER_TAG = "vn2.2"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-iau
        # Blank Upgrade Macro
        return config, self.reports
