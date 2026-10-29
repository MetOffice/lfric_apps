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


class vn32_tXXX(MacroUpgrade):
    """Upgrade macro for ticket #XXX by Ian Boutle."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_tXXX"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lsm
        """Add the soil properties, topographic index, initial water table and
        wetness under the soil column, with the values previously hard-coded"""
        for name, value in [
            ("clapp_horn_b", "11.2"),
            ("mean_topog_index", "1.0"),
            ("soil_albedo", "0.11"),
            ("soil_cond_sat", "0.0015"),
            ("soil_moist_crit", "0.37"),
            ("soil_moist_sat", "0.456"),
            ("soil_moist_wilt", "0.263"),
            ("soil_suction_sat", "0.324"),
            ("soil_thermal_cap", "1230000.0"),
            ("soil_thermal_cond", "0.22"),
            ("stdev_topog_index", "1.0"),
            ("water_table", "1.0"),
            ("wetness_under_soil", "1.0"),
        ]:
            self.add_setting(config, ["namelist:ideal_surface", name], value)
        return config, self.reports
