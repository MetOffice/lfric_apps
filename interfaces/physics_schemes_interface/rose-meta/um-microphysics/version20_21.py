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


class vn20_t467(MacroUpgrade):
    """Upgrade macro for ticket #467 by Mike Whitall."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t467"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-microphysics
        # Set l_mcr_precfrac (switch for prognostic precip fraction);
        # So-far this has been hardwired in the code to be true if using
        # the comorph convection scheme, and false if not.
        # Replicate this when upgrading existing workflows from vn2.0,
        # to ensure they preserve answers
        nml = "namelist:convection"
        cv_scheme = self.get_setting_value(config, [nml, "cv_scheme"])
        if cv_scheme == "'comorph'":
            l_mcr_precfrac = ".true."
        else:
            l_mcr_precfrac = ".false."
        # end if
        # Add new microphysics settings to namelist
        nml = "namelist:microphysics"
        self.add_setting(config, [nml, "i_update_precfrac"], "'homog'")
        self.add_setting(config, [nml, "l_mcr_precfrac"], l_mcr_precfrac)
        self.add_setting(config, [nml, "l_proc_fluxes"], ".false.")
        return config, self.reports


class vn20_t334(MacroUpgrade):
    """Upgrade macro for ticket #334 by Ian Boutle."""

    BEFORE_TAG = "vn2.0_t467"
    AFTER_TAG = "vn2.0_t334"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-microphysics
        self.add_setting(config, ["namelist:microphysics", "mp_dz_scal"], "2.0")
        return config, self.reports


class vn20_t547(MacroUpgrade):
    """Upgrade macro for ticket #547 by Adrian Lock."""

    BEFORE_TAG = "vn2.0_t334"
    AFTER_TAG = "vn2.0_t547"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-microphysics
        nml = "namelist:microphysics"
        self.add_setting(config, [nml, "heavy_rain_evap_fac"], "0.0")
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t547"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-microphysics
        # Blank Upgrade Macro
        return config, self.reports
