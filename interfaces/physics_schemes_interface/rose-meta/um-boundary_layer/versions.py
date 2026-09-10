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

class vn32_t46(MacroUpgrade):
    """Upgrade macro for ticket #46 by Shusuke Nishimoto."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t46"

    def upgrade(self, config, meta_config=None):
        mixing_method = self.get_setting_value(
            config, ["namelist:mixing", "method"]
        )
        self.add_setting(config, ["namelist:blayer", "bl_scheme"], "'Kprof'")
        self.add_setting(
            config, ["namelist:blayer", "blending_kprof"], mixing_method
        )
        self.add_setting(
            config, ["namelist:blayer", "blending_hoc"], "'3dte_mk1'"
        )
        self.add_setting(
            config, ["namelist:blayer", "adv_turb_field"], ".true."
        )
        self.add_setting(config, ["namelist:blayer", "bdy_tke"], "'my3'")
        self.add_setting(
            config, ["namelist:blayer", "local_above_tkelvs"], ".true."
        )
        self.add_setting(config, ["namelist:blayer", "my_condense"], ".true.")
        self.add_setting(
            config, ["namelist:blayer", "my_force_initialise"], ".false."
        )
        self.add_setting(
            config, ["namelist:blayer", "my_ini_dbdz_min"], "1.0e-5"
        )
        self.add_setting(
            config, ["namelist:blayer", "my_lowest_pd_surf"], "'bh91'"
        )
        self.add_setting(config, ["namelist:blayer", "my_prod_adj"], ".true.")
        self.add_setting(
            config, ["namelist:blayer", "my_simeq_solver"], "'gauss'"
        )
        self.add_setting(config, ["namelist:blayer", "shcu_buoy"], ".false.")
        self.add_setting(config, ["namelist:blayer", "shcu_levels"], "-1")
        self.add_setting(config, ["namelist:blayer", "tke_levels"], "-1")
        self.add_setting(config, ["namelist:blayer", "use_l_sq"], ".true.")

        return config, self.reports
