import sys

from metomi.rose.upgrade import MacroUpgrade

from .version21_22 import *


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


class vn22_t202(MacroUpgrade):
    """Upgrade macro for ticket #202 by Katty Huang."""

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t202"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        self.add_setting(
            config, ["namelist:jules_surface", "anthrop_heat_mean"], "20.0"
        )
        self.add_setting(
            config, ["namelist:jules_surface", "anthrop_heat_option"], "'dukes'"
        )
        return config, self.reports


class vn22_t1012(MacroUpgrade):
    """Upgrade macro for ticket #1012 by Maggie Hendry."""

    BEFORE_TAG = "vn2.2_t202"
    AFTER_TAG = "vn2.2_t1012"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Blank upgrade macro to bump tag

        return config, self.reports


class vn22_t814(MacroUpgrade):
    # Upgrade macro for #814 by Maggie Hendry

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t814"

    def upgrade(self, config, meta_config=None):
        # Add jules_model_environment_lfric namelist
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:jules_hydrology",
            r"namelist:jules_hydrology)"
            + "\n"
            + " (namelist:jules_model_environment_lfric",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        self.add_setting(
            config,
            ["namelist:jules_model_environment_lfric", "l_jules_parent"],
            "'lfric'",
        )
        # Add jules_surface namelist items
        self.add_setting(
            config,
            ["namelist:jules_surface", "all_tiles"],
            "'off'",
        )
        self.add_setting(config, ["namelist:jules_surface", "beta1"], "0.83")
        self.add_setting(config, ["namelist:jules_surface", "beta2"], "0.93")
        self.add_setting(
            config, ["namelist:jules_surface", "beta_cnv_bl"], "0.04"
        )
        self.add_setting(
            config,
            ["namelist:jules_surface", "fd_hill_option"],
            "'capped_lowhill'",
        )
        self.add_setting(config, ["namelist:jules_surface", "fwe_c3"], "0.5")
        self.add_setting(
            config, ["namelist:jules_surface", "fwe_c4"], "20000.0"
        )
        self.add_setting(config, ["namelist:jules_surface", "hleaf"], "5.7e4")
        self.add_setting(config, ["namelist:jules_surface", "hwood"], "1.1e4")
        self.add_setting(
            config, ["namelist:jules_surface", "i_modiscopt"], "'on'"
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_epot_corr"], ".true."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_land_ice_imp"], ".true."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_mo_buoyancy_calc"], ".true."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "orog_drag_param"], "0.15"
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_flake_model"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_elev_land_ice"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_elev_lw_down"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_point_data"], ".false."
        )

        return config, self.reports
