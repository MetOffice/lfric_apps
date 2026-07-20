import re
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


class vn31_t464(MacroUpgrade):
    """Upgrade macro for ticket #464 by Ian Boutle."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t464"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        self.add_setting(
            config, ["namelist:cloud", "pc2_turb_horiz"], ".false."
        )
        return config, self.reports


class vn31_t243(MacroUpgrade):
    """Upgrade macro for ticket #243 by Mike Whitall."""

    BEFORE_TAG = "vn3.1_t464"
    AFTER_TAG = "vn3.1_t243"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        nml = "namelist:cloud"
        self.add_setting(config, [nml, "l_ensure_max_in_cloud_pc2"], ".false.")
        return config, self.reports


class vn31_t249(MacroUpgrade):
    """Upgrade macro for ticket #249 by Mike Whitall."""

    BEFORE_TAG = "vn3.1_t243"
    AFTER_TAG = "vn3.1_t249"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        # Blank macro needed just to update meta-data version
        # (apps using the new option 'smooth_fix' under the existing
        #  multi-option switch 'pc2_init_logic' fail checks against
        #  the existing meta-data).
        return config, self.reports


class vn31_t247(MacroUpgrade):
    """Upgrade macro for ticket #247 by Mike Whitall."""

    BEFORE_TAG = "vn3.1_t249"
    AFTER_TAG = "vn3.1_t247"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        # Add new switch controlling PC2 homogeneous forcing option.
        # Previously, this was hardwired in um_physics_init to use one
        # option if using the comorph convection scheme and another if not.
        # So need to implement the same logic here:
        # Load "cv_scheme" from the convection namelist
        nml = "namelist:convection"
        cv_scheme = self.get_setting_value(config, [nml, "cv_scheme"])
        if cv_scheme == "'comorph'":
            # Use the "weight by PDF width" option if using comorph
            pc2_homog_g = "'width'"
        else:
            # Use the "weight as a function of cloud-fraction" option otherwise
            pc2_homog_g = "'cf'"
        # Add new settings with the specified option
        nml = "namelist:cloud"
        self.add_setting(config, [nml, "pc2_homog_g_method"], pc2_homog_g)
        # Rename some other cloud-scheme namelist inputs for clarity...
        nml = "namelist:cloud"
        bm_ez_opt = self.get_setting_value(config, [nml, "i_bm_ez_opt"])
        pc2_erosion_num = self.get_setting_value(
            config, [nml, "i_pc2_erosion_numerics"]
        )
        pc2_init_method = self.get_setting_value(config, [nml, "pc2ini"])
        self.remove_setting(config, [nml, "i_bm_ez_opt"])
        self.remove_setting(config, [nml, "i_pc2_erosion_numerics"])
        self.remove_setting(config, [nml, "pc2ini"])
        self.add_setting(config, [nml, "bm_ez_opt"], bm_ez_opt)
        self.add_setting(config, [nml, "pc2_erosion_numerics"], pc2_erosion_num)
        self.add_setting(config, [nml, "pc2_init_method"], pc2_init_method)
        return config, self.reports


class vn31_t611(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn3.1_t247"
    AFTER_TAG = "vn3.2"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        # Blank Upgrade Macro
        return config, self.reports
