import re
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


class vn32_t717(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn3.2_t760"
    AFTER_TAG = "vn3.2_t717"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-convection
        # Add new comorph namelist
        # Append after 'convection' in configuration.nml
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"(\(?)namelist:convection(\)?)(\n)",
            r"\1namelist:convection\2\3 (namelist:comorph)\n",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        # Move existing comorph namelist entries from the "convection"
        # namelist to the new "comorph" namelist.
        nml1 = "namelist:convection"
        nml2 = "namelist:comorph"
        for entry in [
            "par_gen_mass_fac",
            "par_gen_rhpert",
            "par_radius_ppn_max",
            "resdep_precipramp",
            "dx_ref",
        ]:
            source = self.get_setting_value(config, [nml1, entry])
            self.remove_setting(config, [nml1, entry])
            self.add_setting(config, [nml2, entry], source)
        # Add new namelist entries with previous hardwired default values
        nml = "namelist:comorph"
        # Top-level settings
        self.add_setting(config, [nml, "n_dndraft_types"], "1")
        self.add_setting(config, [nml, "cf_conv_fac"], "2.0")
        self.add_setting(config, [nml, "wind_w_buoy_fac"], "1.0")
        self.add_setting(config, [nml, "overlap_power"], "0.5")
        self.add_setting(config, [nml, "rain_area_min"], "0.05")
        # Conv triggering and parcel initialisation
        self.add_setting(config, [nml, "par_gen_pert_fac"], "0.333")
        self.add_setting(config, [nml, "par_gen_core_fac"], "3.0")
        self.add_setting(
            config, [nml, "par_radius_init_method"], "'linear_p_q'"
        )
        self.add_setting(config, [nml, "par_radius_knob"], "0.45")
        self.add_setting(config, [nml, "par_radius_knob_max"], "2.0")
        self.add_setting(config, [nml, "ass_min_radius"], "500.0")
        # Plume model
        self.add_setting(config, [nml, "ent_coef"], "0.2")
        self.add_setting(config, [nml, "core_ent_cmr"], ".true.")
        self.add_setting(config, [nml, "core_ent_fac"], "1.0")
        self.add_setting(config, [nml, "min_cmr"], "2.0")
        self.add_setting(config, [nml, "max_cmr"], "6.0")
        self.add_setting(config, [nml, "drag_coef_par"], "0.5")
        self.add_setting(
            config, [nml, "par_radius_evol_method"], "'no_detrain'"
        )
        # In-plume microphysics
        self.add_setting(config, [nml, "autoc_opt"], "'quadratic'")
        self.add_setting(config, [nml, "coef_auto"], "0.025")
        self.add_setting(config, [nml, "hetnuc_temp"], "263.0")
        self.add_setting(config, [nml, "cf_area_coef"], "10.0")
        self.add_setting(config, [nml, "wind_w_fac"], "1.0")
        self.add_setting(config, [nml, "col_eff_coef"], "1.0")
        self.add_setting(config, [nml, "drag_coef_cond"], "0.5")
        self.add_setting(config, [nml, "vent_factor"], "0.25")
        self.add_setting(config, [nml, "rho_rim"], "600.0")
        self.add_setting(config, [nml, "nconc_cl"], "1.0E8")
        self.add_setting(config, [nml, "nconc_rain"], "1000.0")
        self.add_setting(config, [nml, "nconc_cf"], "300.0")
        self.add_setting(config, [nml, "nconc_snow"], "300.0")
        self.add_setting(config, [nml, "nconc_graup"], "100.0")
        self.add_setting(config, [nml, "tdep_n_cl"], "0.0")
        self.add_setting(config, [nml, "tdep_n_cf"], "8.18")
        return config, self.reports


class vn32_t754(MacroUpgrade):
    """Upgrade macro for ticket #754 by Chris Smith."""

    BEFORE_TAG = "vn3.2_t717"
    AFTER_TAG = "vn3.2_t754"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        self.add_setting(
            config, ["namelist:initialization", "regrav_interp"], ".false."
        )
        return config, self.reports


class vn32_t800(MacroUpgrade):
    """Upgrade macro for ticket #800 by Chris Smith."""

    BEFORE_TAG = "vn3.2_t754"
    AFTER_TAG = "vn3.2_t800"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Add vapour_relax namelist to configuration source list"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"\(namelist:vapour_forcing\)",
            r"(namelist:vapour_forcing)" + "\n" + " (namelist:vapour_relax)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        # Add vapour_relaxation setting to external_forcing namelist
        self.add_setting(
            config,
            ["namelist:external_forcing", "vapour_relaxation"],
            ".false.",
        )
        # Data for vapour_relax namelist
        self.add_setting(config, ["namelist:vapour_relax"])
        self.add_setting(
            config, ["namelist:vapour_relax", "coordinate"], "'height'"
        )
        self.add_setting(config, ["namelist:vapour_relax", "heights"], "0.0")
        self.add_setting(
            config, ["namelist:vapour_relax", "number_heights"], "1"
        )
        self.add_setting(config, ["namelist:vapour_relax", "number_times"], "1")
        self.add_setting(
            config, ["namelist:vapour_relax", "profile_data"], "0.0"
        )
        self.add_setting(config, ["namelist:vapour_relax", "times"], "0.0")
        self.add_setting(config, ["namelist:vapour_relax", "timescale"], "1.0")
        self.add_setting(config, ["namelist:vapour_relax", "variable"], "'mr'")

        return config, self.reports
