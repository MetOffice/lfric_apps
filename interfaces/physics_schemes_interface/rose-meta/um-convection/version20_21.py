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


class vn20_t334(MacroUpgrade):
    """Upgrade macro for ticket #334 by Ian Boutle."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t334"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-convection
        self.add_setting(config, ["namelist:convection", "efrac"], "1.0")
        self.add_setting(
            config, ["namelist:convection", "orig_mdet_fac"], "1.0"
        )
        self.add_setting(config, ["namelist:convection", "prog_ent_min"], "0.5")
        self.add_setting(config, ["namelist:convection", "qlmin"], "4.0e-4")
        return config, self.reports


class vn20_t552(MacroUpgrade):
    """Upgrade macro for ticket #552 by Adrian Lock."""

    BEFORE_TAG = "vn2.0_t334"
    AFTER_TAG = "vn2.0_t552"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-convection
        nml = "namelist:convection"
        self.add_setting(config, [nml, "resdep_precipramp"], ".false.")
        self.add_setting(config, [nml, "dx_ref"], "50000.0")
        return config, self.reports


class vn20_t472(MacroUpgrade):
    """Upgrade macro for ticket #472 by Mike Whitall."""

    BEFORE_TAG = "vn2.0_t552"
    AFTER_TAG = "vn2.0_t472"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-convection
        nml = "namelist:convection"
        self.add_setting(config, [nml, "l_cvdiag_ctop_qmax"], ".false.")
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t472"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-convection
        # Blank Upgrade Macro
        return config, self.reports
