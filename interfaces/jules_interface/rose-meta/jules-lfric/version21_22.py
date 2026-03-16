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


class vn21_t596(MacroUpgrade):
    """Upgrade macro for ticket #596 by Maggie Hendry."""

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t596"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Only used in LFRic apps
        self.rename_setting(
            config,
            ["namelist:jules_surface", "check_soilm_negatives"],
            ["namelist:surface", "check_soilm_negatives"],
        )
        self.rename_setting(
            config,
            ["namelist:jules_surface", "lake_water_conservation"],
            ["namelist:surface", "lake_water_conservation"],
        )
        # jules-sea-seaice from jules-lfric being shared rather than um-atmos
        self.rename_setting(
            config,
            ["namelist:surface", "amip_ice_thick"],
            ["namelist:jules_sea_seaice", "amip_ice_thick"],
        )
        self.rename_setting(
            config,
            ["namelist:surface", "z0h_specified"],
            ["namelist:jules_sea_seaice", "z0h_specified"],
        )
        self.rename_setting(
            config,
            ["namelist:surface", "z0m_specified"],
            ["namelist:jules_sea_seaice", "z0m_specified"],
        )
        return config, self.reports


class vn21_t871(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.1_t596"
    AFTER_TAG = "vn2.2"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Blank Upgrade Macro
        return config, self.reports
