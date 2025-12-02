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


class vn20_t334(MacroUpgrade):
    """Upgrade macro for ticket #334 by Ian Boutle."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t334"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-boundary_layer
        self.add_setting(config, ["namelist:blayer", "dec_thres_cu"], "0.05")
        return config, self.reports


class vn20_t562(MacroUpgrade):
    """Upgrade macro for ticket #562 by Paul Burns."""

    BEFORE_TAG = "vn2.0_t334"
    AFTER_TAG = "vn2.0_t562"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-boundary_layer
        """Add ng_stress to namelist blayer"""
        self.add_setting(
            config, ["namelist:blayer", "ng_stress"], "'BG97_limited'"
        )
        return config, self.reports


class vn20_t472(MacroUpgrade):
    """Upgrade macro for ticket #472 by Mike Whitall."""

    BEFORE_TAG = "vn2.0_t562"
    AFTER_TAG = "vn2.0_t472"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-boundary_layer
        nml = "namelist:blayer"
        self.add_setting(config, [nml, "dzrad_disc_opt"], "'level_ntm1'")
        self.add_setting(config, [nml, "entr_smooth_dec"], "'on'")
        self.add_setting(config, [nml, "l_converge_ga"], ".false.")
        self.add_setting(config, [nml, "l_use_sml_dsc_fixes"], ".false.")
        self.add_setting(config, [nml, "num_sweeps_bflux"], "3")
        # New setting sc_diag_opt replaces relax_sc_over_cu
        l_relax = self.get_setting_value(config, [nml, "relax_sc_over_cu"])
        self.remove_setting(config, [nml, "relax_sc_over_cu"])
        if l_relax == ".true.":
            sc_diag_opt = "'cu_relax'"
        else:
            sc_diag_opt = "'orig'"
        self.add_setting(config, [nml, "sc_diag_opt"], sc_diag_opt)
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t472"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-boundary_layer
        # Blank Upgrade Macro
        return config, self.reports
