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


class vn20_t249(MacroUpgrade):
    """Upgrade macro for ticket #249 by Denis Sergeev."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t249"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/jules_interface/rose-meta/jules-lfric
        nml = "namelist:specified_surface"
        self.add_setting(config, [nml, "surf_temp_forcing"], "'none'")
        self.add_setting(config, [nml, "internal_flux_method"], "'uniform'")
        self.add_setting(config, [nml, "internal_flux_value"], "0.0")
        return config, self.reports


class vn20_t588(MacroUpgrade):
    """Upgrade macro for ticket #588 by Maggie Hendry."""

    BEFORE_TAG = "vn2.0_t249"
    AFTER_TAG = "vn2.0_t588"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Rename surface -> jules_radiation
        items_surface = [
            "sea_alb_method",
            "albedo_obs",
            "non_iso_scatter",
            "sea_alb_var_chl",
            "blue_sky_alb",
        ]
        items_jules = [
            "i_sea_alb_method",
            "l_albedo_obs",
            "l_niso_direct",
            "l_sea_alb_var_chl",
            "l_spec_alb_bs",
        ]
        for surface, jules in zip(items_surface, items_jules):
            self.rename_setting(
                config,
                ["namelist:surface", surface],
                ["namelist:jules_radiation", jules],
            )
        # Rename surface -> jules_snow
        items_jules = [
            "i_basal_melting_opt",
            "i_grain_growth_opt",
            "i_relayer_opt",
        ]
        items_surface = [
            "basal_melting",
            "grain_growth",
            "relayer_opt",
        ]
        for surface, jules in zip(items_surface, items_jules):
            self.rename_setting(
                config,
                ["namelist:surface", surface],
                ["namelist:jules_snow", jules],
            )
        # jules_soil
        """Add jules_soil namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:jules_radiation",
            r"namelist:jules_radiation)" + "\n" + " (namelist:jules_soil",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        # Rename surface -> jules_soil
        items_surface = [
            "dpsids_dsdz",
            "soil_sat_down",
            "l_vg_soil",
        ]
        items_jules = [
            "l_dpsids_dsdz",
            "l_soil_sat_down",
            "l_vg_soil",
        ]
        for surface, jules in zip(items_surface, items_jules):
            self.rename_setting(
                config,
                ["namelist:surface", surface],
                ["namelist:jules_soil", jules],
            )
        # jules_pftparm
        """Add jules_pftparm namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:jules_nvegparm",
            r"namelist:jules_nvegparm)" + "\n" + " (namelist:jules_pftparm",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        # Rename surface -> jules_pftparm
        items_surface = [
            "alb_snocov_max",
            "alb_leaf_nir",
            "alb_leaf_vis",
            "catch0",
            "dcatch_dlai",
            "fsmc_p0",
            "light_extinct",
            "knl",
            "scat_coef_vis",
            "scat_coef_nir",
            "z0v",
            "z0hm_ratio_pft",
        ]
        items_jules = [
            "albsnc_max_io",
            "alnir_io",
            "alpar_io",
            "catch0_io",
            "dcatch_dlai_io",
            "fsmc_p0_io",
            "kext_io",
            "knl_io",
            "omega_io",
            "omnir_io",
            "z0v_io",
            "z0hm_pft_io",
        ]
        for surface, jules in zip(items_surface, items_jules):
            self.rename_setting(
                config,
                ["namelist:surface", surface],
                ["namelist:jules_pftparm", jules],
            )
        # jules_sea_seaice
        """Add jules_sea_seaice namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:jules_radiation",
            r"namelist:jules_radiation)" + "\n" + " (namelist:jules_sea_seaice",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        # Rename surface -> jules_sea_seaice
        items_surface = [
            "alb_sice_melt",
            "buddy_sea",
            "cdn_hw_sea",
            "cdn_max_sea",
            "dt_ice_albedo",
            "emis_sea",
            "emis_sice",
            "evap_scale_sea",
            "heat_cap_sea",
            "i_high_wind_drag",
            "iceformdrag_lupkes",
            "l_10m_neut",
            "n_sea_ice_tile",
            "sea_surf_alg",
            "sice_heatflux",
            "stability_lupkes",
            "therm_cond_sea",
            "therm_cond_sice",
            "therm_cond_sice_snow",
            "u_cdn_hw",
            "u_cdn_max",
            "use_variable_sst",
        ]
        items_jules = [
            "alpham",
            "buddy_sea",
            "cdn_hw_sea",
            "cdn_max_sea",
            "dtice",
            "emis_sea",
            "emis_sice",
            "beta_evap",
            "hcap_sea",
            "i_high_wind_drag",
            "l_iceformdrag_lupkes",
            "l_10m_neut",
            "nice",
            "iseasurfalg",
            "l_sice_heatflux",
            "l_stability_lupkes",
            "kappa_seasurf",
            "kappai",
            "kappai_snow",
            "u_cdn_hw",
            "u_cdn_max",
            "l_use_dtstar_sea",
        ]
        for surface, jules in zip(items_surface, items_jules):
            self.rename_setting(
                config,
                ["namelist:surface", surface],
                ["namelist:jules_sea_seaice", jules],
            )
        # jules_hydrology
        """Add jules_hydrology namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:io",
            r"namelist:io" + "\n" + " (namelist:jules_hydrology)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        # Rename surface -> jules_hydrology
        items_surface = [
            "use_hydrology",
            "l_var_rainfrac",
        ]
        items_jules = [
            "l_hydrology",
            "l_var_rainfrac",
        ]
        for surface, jules in zip(items_surface, items_jules):
            self.rename_setting(
                config,
                ["namelist:surface", surface],
                ["namelist:jules_hydrology", jules],
            )
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t588"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Blank Upgrade Macro
        return config, self.reports
