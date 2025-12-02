import re
import sys

from metomi.rose.upgrade import MacroUpgrade


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


class vn20_t481(MacroUpgrade):
    """Upgrade macro for ticket #481 by Mike Whitall."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t481"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-cloud
        # All altered settings are in the cloud namelist
        nml = "namelist:cloud"
        # Logical ez_subcrit_only replaced by a 3-way integer switch:
        ez_subcrit = self.get_setting_value(config, [nml, "ez_subcrit"])
        self.remove_setting(config, [nml, "ez_subcrit"])
        if ez_subcrit == ".true.":
            i_bm_ez = "'subcrit'"
        else:
            i_bm_ez = "'orig'"
        # end if
        self.add_setting(config, [nml, "i_bm_ez_opt"], i_bm_ez)
        # New settings added
        self.add_setting(config, [nml, "ent_coef_bm"], "0.2")
        self.add_setting(config, [nml, "l_bm_sigma_s_grad"], ".false.")
        self.add_setting(config, [nml, "l_bm_tweaks"], ".false.")
        self.add_setting(config, [nml, "max_sigmas"], "3.0")
        self.add_setting(config, [nml, "min_sigx_ft"], "0.0")
        self.add_setting(config, [nml, "turb_var_fac_bm"], "1.0")
        return config, self.reports


class vn20_t334(MacroUpgrade):
    """Upgrade macro for ticket #334 by Ian Boutle."""

    BEFORE_TAG = "vn2.0_t481"
    AFTER_TAG = "vn2.0_t334"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-cloud
        cvscheme = self.get_setting_value(
            config, ["namelist:convection", "cv_scheme"]
        )
        if cvscheme == "'comorph'":
            self.add_setting(
                config, ["namelist:cloud", "fsd_conv_const"], "3.0"
            )
            self.add_setting(
                config, ["namelist:cloud", "fsd_min_conv_frac"], "0.02"
            )
            self.add_setting(
                config, ["namelist:cloud", "fsd_nonconv_const"], "0.8"
            )
        else:
            self.add_setting(
                config, ["namelist:cloud", "fsd_conv_const"], "2.81"
            )
            self.add_setting(
                config, ["namelist:cloud", "fsd_min_conv_frac"], "0.0"
            )
            self.add_setting(
                config, ["namelist:cloud", "fsd_nonconv_const"], "1.14"
            )
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t334"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        # Blank Upgrade Macro
        return config, self.reports
