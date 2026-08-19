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


class vn32_t670(MacroUpgrade):
    """Upgrade macro for ticket #670 by Thomas Bendall."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t670"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Add new nudging namelist options
        self.add_setting(
            config, ["namelist:nudging", "nudging_method"], "'convolution'"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_relax_time_theta"], "6.0"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_relax_time_u"], "6.0"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_relax_time_v"], "6.0"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_spinup_start"], "12.0"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_spinup_end"], "24.0"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_stop_time"], "10000.0"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_level_taper_bottom"], "5"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_level_taper_top"], "52"
        )
        self.add_setting(
            config, ["namelist:nudging", "nudging_min_tropopause_level"], "48"
        )
        self.change_setting_value(
            config, ["namelist:nudging", "nudging_level_bottom"], "6"
        )
        self.add_setting(
            config, ["namelist:nudging", "num_ref_data_levels"], "137"
        )
        self.add_setting(config, ["namelist:nudging", "spectral_kmax"], "20")
        self.add_setting(config, ["namelist:nudging", "spectral_kmin"], "2")
        self.add_setting(
            config, ["namelist:nudging", "spectral_stencil_extent"], "12"
        )
        self.add_setting(
            config, ["namelist:nudging", "spectral_envelope_width"], "0.1"
        )
        # Remove retired settings
        self.remove_setting(config, ["namelist:nudging", "nudging_source"])
        self.remove_setting(config, ["namelist:nudging", "nudging_width_bottom"])
        self.remove_setting(config, ["namelist:nudging", "nudging_width_top"])
        self.remove_setting(config, ["namelist:nudging", "nudge_data_levels"])
        # If nudging_mesh_name is still the default '' value, set it to
        # match dynamics_mesh_name
        nudging_mesh_name = self.get_setting_value(
            config, ["namelist:multires_coupling", "nudging_mesh_name"]
        )
        if nudging_mesh_name == "''":
            dynamics_mesh_name = self.get_setting_value(
                config, ["namelist:multires_coupling", "dynamics_mesh_name"]
            )
            if dynamics_mesh_name is not None:
                self.change_setting_value(
                    config,
                    ["namelist:multires_coupling", "nudging_mesh_name"],
                    dynamics_mesh_name,
                )

        return config, self.reports
