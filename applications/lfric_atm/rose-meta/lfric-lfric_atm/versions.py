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


class vn31_t118(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t118"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Blank Upgrade Macro
        return config, self.reports


class vn31_t363(MacroUpgrade):
    """Upgrade macro for ticket #363 by Jaffery Irudayasamy."""

    BEFORE_TAG = "vn3.1_t118"
    AFTER_TAG = "vn3.1_t363"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Set segmentation size limit for short and long wave radiation kernels"""
        self.add_setting(config, ["namelist:physics", "sw_segment_limit"], "32")
        self.add_setting(config, ["namelist:physics", "lw_segment_limit"], "32")

        return config, self.reports


class vn31_t348(MacroUpgrade):
    """Upgrade macro for ticket #348 by Ian Boutle."""

    BEFORE_TAG = "vn3.1_t363"
    AFTER_TAG = "vn3.1_t348"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Use PMSL halo calculations by default
        self.add_setting(
            config, ["namelist:physics", "pmsl_halo_calcs"], ".true."
        )

        return config, self.reports


class vn31_t368(MacroUpgrade):
    """Upgrade macro for ticket #368 by Ian Boutle."""

    BEFORE_TAG = "vn3.1_t348"
    AFTER_TAG = "vn3.1_t368"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-convection
        self.add_setting(
            config, ["namelist:convection", "llcs_first_outer"], ".false."
        )

        return config, self.reports


class vn31_t238(MacroUpgrade):
    """Upgrade macro for ticket #238 by Thomas Bendall."""

    BEFORE_TAG = "vn3.1_t368"
    AFTER_TAG = "vn3.1_t238"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(
            config, ["namelist:finite_element", "coord_space"], "'Wchi'"
        )
        coord_order = self.get_setting_value(
            config, ["namelist:finite_element", "coord_order"]
        )
        self.add_setting(
            config,
            ["namelist:finite_element", "coord_order_nonprime"],
            coord_order,
        )

        return config, self.reports


class vn31_t443(MacroUpgrade):
    """Upgrade macro for ticket #443 by Samantha Pullen."""

    BEFORE_TAG = "vn3.1_t238"
    AFTER_TAG = "vn3.1_t443"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Add name entry to iau_addinf_io namelist
        self.add_setting(
            config, ["namelist:iau_addinf_io(addinf1)", "name"], "''"
        )
        self.add_setting(
            config, ["namelist:iau_addinf_io(addinf2)", "name"], "''"
        )
        # Add name entry to iau_ainc_io namelist
        self.add_setting(config, ["namelist:iau_ainc_io(ainc1)", "name"], "''")
        self.add_setting(config, ["namelist:iau_ainc_io(ainc2)", "name"], "''")
        # Add name entry to iau_bcorr_io namelist
        self.add_setting(
            config, ["namelist:iau_bcorr_io(bcorr1)", "name"], "''"
        )

        return config, self.reports


class vn31_t464(MacroUpgrade):
    """Upgrade macro for ticket #464 by Ian Boutle."""

    BEFORE_TAG = "vn3.1_t443"
    AFTER_TAG = "vn3.1_t464"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        self.add_setting(
            config, ["namelist:cloud", "pc2_turb_horiz"], ".false."
        )

        return config, self.reports


class vn31_t382(MacroUpgrade):
    """Upgrade macro for ticket #382 by Benjamin Went."""

    BEFORE_TAG = "vn3.1_t464"
    AFTER_TAG = "vn3.1_t382"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-lfric_atm
        """Set segmentation size for the Boundary Layer"""
        self.change_setting_value(
            config, ["namelist:physics", "bl_segment"], "16"
        )

        return config, self.reports


class vn31_t243(MacroUpgrade):
    """Upgrade macro for ticket #243 by Mike Whitall."""

    BEFORE_TAG = "vn3.1_t382"
    AFTER_TAG = "vn3.1_t243"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-microphysics
        nml = "namelist:microphysics"
        self.add_setting(config, [nml, "l_improve_precfrac_checks"], ".false.")

        # Commands From: rose-meta/um-cloud
        nml = "namelist:cloud"
        self.add_setting(config, [nml, "l_ensure_max_in_cloud_pc2"], ".false.")

        return config, self.reports
