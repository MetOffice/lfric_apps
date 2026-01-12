import sys

from metomi.rose.upgrade import MacroUpgrade

from .version21_22 import *


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


class vn22_t850(MacroUpgrade):
    """Upgrade macro for ticket #850 by Shusuke Nishimoto."""

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t850"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-stochastic_physics
        idealised_test_name = self.get_setting_value(
            config, ["namelist:idealised", "test"]
        )
        l_multigrid = self.get_setting_value(
            config, ["namelist:formulation", "l_multigrid"]
        )
        limited_area = self.get_setting_value(
            config, ["namelist:boundaries", "limited_area"]
        )
        if (
            idealised_test_name == "'none'"
            and limited_area == ".true."
            and l_multigrid == ".true."
        ):
            self.change_setting_value(
                config,
                ["namelist:section_choice", "stochastic_physics"],
                "'um'",
            )
            self.add_setting(
                config,
                ["namelist:physics", "stochastic_physics_placement"],
                "'fast'",
            )
            blpert_type = "'theta_and_moist'"
            mesh_names = self.get_setting_value(
                config, ["namelist:multigrid", "chain_mesh_tags"]
            )
            coarsest_mesh_name = mesh_names.split(",")[-1]
        else:
            blpert_type = "'off'"
            coarsest_mesh_name = "''"
        self.add_setting(
            config, ["namelist:stochastic_physics", "blpert_type"], blpert_type
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_mesh_name"],
            coarsest_mesh_name,
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_time_correlation"],
            ".true.",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_decorrelation_time"],
            "600.0",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_only_near_edge"],
            ".true.",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_npts_from_edge"],
            "24",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_noncumulus_points"],
            ".false.",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_height_bottom"],
            "0.0",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_height_top"],
            "1500.0",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_add_vertical_shape"],
            ".true.",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "blpert_max_magnitude"],
            "1.0",
        )
        return config, self.reports


class vn22_t34(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.2_t850"
    AFTER_TAG = "vn3.0"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-stochastic_physics
        # Blank Upgrade Macro
        return config, self.reports
