import re
import sys

from metomi.rose.upgrade import MacroUpgrade

from .version22_30 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn30_t99(MacroUpgrade):
    """Upgrade macro for ticket #99 by Fred Wobus."""

    BEFORE_TAG = "vn3.0"
    AFTER_TAG = "vn3.0_t99"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-lfric_atm
        """Set segmentation size for Gregory-Rowntree convection kernel"""
        self.add_setting(config, ["namelist:physics", "conv_gr_segment"], "16")
        return config, self.reports


class vn30_t146(MacroUpgrade):
    """Upgrade macro for ticket #146 by Maggie Hendry."""

    BEFORE_TAG = "vn3.0_t99"
    AFTER_TAG = "vn3.0_t146"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
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


class vn30_t135(MacroUpgrade):
    """Upgrade macro for ticket #135 by James Manners."""

    BEFORE_TAG = "vn3.0_t146"
    AFTER_TAG = "vn3.0_t135"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/socrates-radiation
        self.add_setting(config, ["namelist:cosp", "n_cosp_step"], "1")
        return config, self.reports


class vn30_t171(MacroUpgrade):
    """Upgrade macro for ticket #171 by James Kent."""

    BEFORE_TAG = "vn3.0_t135"
    AFTER_TAG = "vn3.0_t171"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Add adjust_tracer_equation to transport namelist
        self.add_setting(
            config, ["namelist:transport", "adjust_tracer_equation"], ".false."
        )
        return config, self.reports


class vn30_t214(MacroUpgrade):
    """Upgrade macro for ticket #214 by mark Hedley."""

    BEFORE_TAG = "vn3.0_t171"
    AFTER_TAG = "vn3.0_t214"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Set segments configuration to true."""
        self.change_setting_value(
            config, ["namelist:physics", "configure_segments"], ".true."
        )

        return config, self.reports
