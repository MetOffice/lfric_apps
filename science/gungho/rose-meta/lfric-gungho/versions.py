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


class vn32_t634(MacroUpgrade):
    """Upgrade macro for ticket #634 by Ian Boutle."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t634"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        nml = "namelist:boundaries"
        self.add_setting(config, [nml, "lbc_bal_meth"], "'keep_rho'")
        self.add_setting(config, [nml, "lbc_sort_theta"], ".true.")
        nml = "namelist:initialization"
        eos_height = self.get_setting_value(config, [nml, "model_eos_height"])
        self.remove_setting(config, [nml, "model_eos_height"])
        self.add_setting(config, [nml, "init_eos_height"], eos_height)
        self.add_setting(config, [nml, "init_exner_method"], "'hydrostatic'")
        self.add_setting(config, [nml, "init_sort_theta"], ".true.")
        return config, self.reports


class vn32_t479(MacroUpgrade):
    """Upgrade macro for ticket #479 by Shusuke Nishimoto."""

    BEFORE_TAG = "vn3.2_t634"
    AFTER_TAG = "vn3.2_t479"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        self.add_setting(config, ["namelist:mixing", "fullstress"], ".false.")
        return config, self.reports


class vn32_t744(MacroUpgrade):
    """Upgrade macro for ticket #744 by Maggie Hendry."""

    BEFORE_TAG = "vn3.2_t479"
    AFTER_TAG = "vn3.2_t744"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lsm
        # Bump tag to pick up metadata changes
        return config, self.reports


class vn32_t698(MacroUpgrade):
    """Upgrade macro for ticket #698 by Alan J Hewitt."""

    BEFORE_TAG = "vn3.2_t744"
    AFTER_TAG = "vn3.2_t698"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-aerosol
        # Add new settings with the default option SUBCOCSSDU_7mode
        self.add_setting(
            config, ["namelist:aerosol", "mode_setup"], "'SUBCOCSSDU_7mode'"
        )
        # Default to false since this is the setting in all existing tests
        self.add_setting(
            config, ["namelist:aerosol", "l_dust_mp_ageing"], ".false."
        )
        # Default to true since this is the setting in all existing tests
        self.add_setting(
            config, ["namelist:aerosol", "l_ukca_radaer_sustrat"], ".true."
        )
        return config, self.reports


class vn32_t725(MacroUpgrade):
    """Upgrade macro for ticket #725 by Ian Boutle."""

    BEFORE_TAG = "vn3.2_t698"
    AFTER_TAG = "vn3.2_t725"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        self.add_setting(
            config, ["namelist:mixing", "leonard_inc_ice"], ".false."
        )
        self.add_setting(
            config, ["namelist:mixing", "leonard_inc_with_bl"], ".false."
        )
        return config, self.reports


class vn32_t699(MacroUpgrade):
    """Upgrade macro for ticket #699 by thomas.melvin."""

    BEFORE_TAG = "vn3.2_t725"
    AFTER_TAG = "vn3.2_t699"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Add native_w2_wind_transport to namelist transport"""
        self.add_setting(
            config,
            ["namelist:transport", "native_w2_wind_transport"],
            ".false.",
        )
        return config, self.reports


class vn32_t760(MacroUpgrade):
    """Upgrade macro for ticket #760 by Chris Smith."""

    BEFORE_TAG = "vn3.2_t699"
    AFTER_TAG = "vn3.2_t760"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        self.add_setting(
            config,
            ["namelist:initial_temperature", "profile_variable"],
            "'potential'",
        )
        self.add_setting(
            config, ["namelist:initial_vapour", "profile_variable"], "'mr'"
        )

        return config, self.reports

class vn32_t670(MacroUpgrade):
    """Upgrade macro for ticket #670 by Thomas Bendall."""

    BEFORE_TAG = "vn3.2_t760"
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