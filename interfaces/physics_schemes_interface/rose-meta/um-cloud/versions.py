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

class vn31_t247(MacroUpgrade):
    # Upgrade macro for Issue #247 by Mike Whitall

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t247"

    def upgrade(self, config, meta_config=None):
        # Add settings

        # Add new switch controlling PC2 homogeneous forcing option.
        # Previously, this was hardwired in um_physics_init to use one
        # option if using the comorph convection scheme and another if not.
        # So need to implement the same logic here:

        # Load "cv_scheme" from the convection namelist
        nml = "namelist:convection"
        cv_scheme = self.get_setting_value(config, [nml, "cv_scheme"])

        if (cv_scheme == "'comorph'"):
            # Use the "weight by PDF width" option if using comorph
            i_pc2_homog_g = "'width'"
        else:
            # Use the "weight as a function of cloud-fraction" option otherwise
            i_pc2_homog_g = "'cf'"

        # Add new settings with the specified option
        nml = "namelist:cloud"
        self.add_setting(config, [nml, "i_pc2_homog_g_method"], i_pc2_homog_g)

        return config, self.reports
