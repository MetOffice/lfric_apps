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
        # 0.66 and 1.2 are tuned GC6 values (0.5 and 0.8 originally)
        self.add_setting(config, ["namelist:convection","r_det"],"0.5")
        self.add_setting(config, ["namelist:convection","cca_md_scaling"],"0.8")
        # These settings are unchanged
        self.add_setting(config, ["namelist:convection","prog_ent_grad"],"-1.1")
        self.add_setting(config, ["namelist:convection","prog_ent_int"],"-2.9")
        self.add_setting(config, ["namelist:convection","prog_ent_max"],"2.5")
        self.add_setting(config, ["namelist:convection","cpress_term"],"0.3")
        self.add_setting(config, ["namelist:convection","ent_fac_sh"],"1.0")
        self.add_setting(config, ["namelist:convection","mparwtr"],"1.0e-3")
        self.add_setting(config, ["namelist:convection","thpixs_mid"],"0.5")
        self.add_setting(config, ["namelist:convection","c_mass_sh"],"0.03")
        cv_scheme = self.get_setting_value(config, ["namelist:convection", "cv_scheme"])
        if cv_scheme == "'comorph'":
            self.add_setting(config, ["namelist:convection","conv_prog_dtheta"],".false.")
            self.add_setting(config, ["namelist:convection","conv_prog_dq"],".false.")
        else:
            self.add_setting(config, ["namelist:convection","conv_prog_dtheta"],".true.")
            self.add_setting(config, ["namelist:convection","conv_prog_dq"],".true.")
        return config, self.reports

