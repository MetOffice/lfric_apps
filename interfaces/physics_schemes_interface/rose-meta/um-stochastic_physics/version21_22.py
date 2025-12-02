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


class vn21_t164(MacroUpgrade):
    """Upgrade macro for ticket #164 by annemccabe."""

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t164"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-stochastic_physics
        """Add new namelist settings to stochastic_physics"""
        alnir = self.get_setting_value(
            config, ["namelist:jules_pftparm", "alnir_io"]
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alnir"], alnir
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alnir_min"], alnir
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alnir_max"], alnir
        )
        alpar = self.get_setting_value(
            config, ["namelist:jules_pftparm", "alpar_io"]
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alpar"], alpar
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alpar_min"], alpar
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alpar_max"], alpar
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_mp_fxd_cld_num"],
            "150.0e+06,150.0e+06,150.0e+06",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_mp_ice_fspd"],
            "6000000.0,6000000.0,6000000.0",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_mp_mp_czero"],
            "10.0,10.0,10.0",
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_mp_mpof"], "0.5,0.5,0.5"
        )
        omega = self.get_setting_value(
            config, ["namelist:jules_pftparm", "omega_io"]
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omega"], omega
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omega_min"], omega
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omega_max"], omega
        )
        omnir = self.get_setting_value(
            config, ["namelist:jules_pftparm", "omnir_io"]
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omnir"], omnir
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omnir_min"], omnir
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omnir_max"], omnir
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_lsfc_orog_drag_param"],
            "0.15,0.15,0.15",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_mp_snow_fspd"],
            "12.0,12.0,12.0",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_lsfc_z0_urban_mult"],
            "1.0,1.0,1.0",
        )
        soil = self.get_setting_value(
            config, ["namelist:jules_surface_types", "soil"]
        )
        npft = self.get_setting_value(
            config, ["namelist:jules_surface_types", "npft"]
        )
        z0_nvg_io = self.get_setting_value(
            config, ["namelist:jules_nvegparm", "z0_nvg_io"]
        )
        z0_soil = z0_nvg_io.split(",")[int(soil) - int(npft) - 1]
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_lsfc_z0_soil"],
            z0_soil + "," + z0_soil + "," + z0_soil,
        )
        z0hm_nvg_io = self.get_setting_value(
            config, ["namelist:jules_nvegparm", "z0hm_nvg_io"]
        )
        z0hm_soil = z0hm_nvg_io.split(",")[int(soil) - int(npft) - 1]
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_lsfc_z0hm_soil"],
            z0hm_soil + "," + z0hm_soil + "," + z0hm_soil,
        )
        z0v_io = self.get_setting_value(config, ["namelist:surface", "z0v"])
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_z0v"], z0v_io
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_z0v_min"], z0v_io
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_z0v_max"], z0v_io
        )
        return config, self.reports


class vn21_t871(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.1_t164"
    AFTER_TAG = "vn2.2"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-stochastic_physics
        # Blank Upgrade Macro
        return config, self.reports
