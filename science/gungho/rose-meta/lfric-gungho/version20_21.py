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


class vn20_t85(MacroUpgrade):
    """Upgrade macro for ticket #85 by Chris Smith."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t85"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/gungho/rose-meta/lfric-gungho
        """Add new geostrophic_forcing namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:formulation",
            r"namelist:formulation" + "\n" + " (namelist:geostrophic_forcing)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add geostrophic_forcing setting to external_forcing namelist"""
        self.add_setting(
            config,
            ["namelist:external_forcing", "geostrophic_forcing"],
            ".false.",
        )
        """Add default data for geostrophic_forcing namelist"""
        self.add_setting(config, ["namelist:geostrophic_forcing"])
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "coordinate"], "'height'"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "heights"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "number_heights"], "1"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "number_times"], "1"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "profile_data_u"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "profile_data_v"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "times"], "0.0"
        )
        return config, self.reports


class vn20_t358(MacroUpgrade):
    """Upgrade macro for ticket #358 by Joshua Dendy."""

    BEFORE_TAG = "vn2.0_t85"
    AFTER_TAG = "vn2.0_t358"

    def upgrade(self, config, meta_config=None):
        # Commands From: components/driver/rose-meta/lfric-driver
        """
        Add element_order_h and element_order_v to namelist finite_element,
        replacing element_order
        """
        self.add_setting(
            config, ["namelist:finite_element", "element_order_h"], "0"
        )
        self.add_setting(
            config, ["namelist:finite_element", "element_order_v"], "0"
        )
        self.remove_setting(
            config, ["namelist:finite_element", "element_order"]
        )
        return config, self.reports


class vn20_t82(MacroUpgrade):
    """Upgrade macro for ticket #82 by Chris Smith."""

    BEFORE_TAG = "vn2.0_t358"
    AFTER_TAG = "vn2.0_t82"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/gungho/rose-meta/lfric-gungho
        # Commands From: science/gungho/rose-meta/lfric-gungho
        """Add vapour_forcing namelist to configuration source list"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = source + "\n" + " (namelist:vapour_forcing)"
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add vapour_forcing setting to external_forcing namelist"""
        self.add_setting(
            config, ["namelist:external_forcing", "vapour_forcing"], "'none'"
        )
        """Data for vapour_forcing namelist"""
        self.add_setting(config, ["namelist:vapour_forcing"])
        self.add_setting(
            config, ["namelist:vapour_forcing", "coordinate"], "'height'"
        )
        self.add_setting(config, ["namelist:vapour_forcing", "heights"], "0.0")
        self.add_setting(
            config, ["namelist:vapour_forcing", "number_heights"], "1"
        )
        self.add_setting(
            config, ["namelist:vapour_forcing", "number_times"], "1"
        )
        self.add_setting(
            config, ["namelist:vapour_forcing", "profile_data"], "0.0"
        )
        self.add_setting(config, ["namelist:vapour_forcing", "times"], "0.0")
        return config, self.reports


class vn20_t467(MacroUpgrade):
    """Upgrade macro for ticket #467 by Mike Whitall."""

    BEFORE_TAG = "vn2.0_t82"
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


class vn20_t339(MacroUpgrade):
    """Upgrade macro for ticket #339 by Joshua Dendy."""

    BEFORE_TAG = "vn2.0_t467"
    AFTER_TAG = "vn2.0_t339"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/gungho/rose-meta/lfric-gungho
        self.add_setting(
            config, ["namelist:idealised", "perturb_init"], ".false."
        )
        self.add_setting(
            config, ["namelist:idealised", "perturb_magnitude"], "0"
        )
        self.add_setting(config, ["namelist:idealised", "perturb_seed"], "0")
        return config, self.reports


class vn20_t541(MacroUpgrade):
    """Upgrade macro for ticket #541 by James Manners."""

    BEFORE_TAG = "vn2.0_t339"
    AFTER_TAG = "vn2.0_t541"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/socrates_interface/rose-meta/socrates-radiation
        self.add_setting(
            config, ["namelist:radiation", "cloud_entrapment"], "'zero'"
        )
        return config, self.reports


class vn20_t481(MacroUpgrade):
    """Upgrade macro for ticket #481 by Mike Whitall."""

    BEFORE_TAG = "vn2.0_t541"
    AFTER_TAG = "vn2.0_t481"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-cloud
        # All altered settings are in the cloud namelist
        nml = "namelist:cloud"
        # Logical ez_subcrit_only replaced by a 3-way integer switch:
        ez_subcrit = self.get_setting_value(config, [nml, "ez_subcrit"])
        self.remove_setting(config, [nml, "ez_subcrit"])
        if ez_subcrit == ".true.":
            i_bm_ez = "'subcrit'"
        else:
            i_bm_ez = "'orig'"
        # end if
        self.add_setting(config, [nml, "i_bm_ez_opt"], i_bm_ez)
        # New settings added
        self.add_setting(config, [nml, "ent_coef_bm"], "0.2")
        self.add_setting(config, [nml, "l_bm_sigma_s_grad"], ".false.")
        self.add_setting(config, [nml, "l_bm_tweaks"], ".false.")
        self.add_setting(config, [nml, "max_sigmas"], "3.0")
        self.add_setting(config, [nml, "min_sigx_ft"], "0.0")
        self.add_setting(config, [nml, "turb_var_fac_bm"], "1.0")
        return config, self.reports


class vn20_t334(MacroUpgrade):
    """Upgrade macro for ticket #334 by Ian Boutle."""

    BEFORE_TAG = "vn2.0_t481"
    AFTER_TAG = "vn2.0_t334"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-microphysics
        self.add_setting(config, ["namelist:microphysics", "mp_dz_scal"], "2.0")
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-convection
        self.add_setting(config, ["namelist:convection", "efrac"], "1.0")
        self.add_setting(
            config, ["namelist:convection", "orig_mdet_fac"], "1.0"
        )
        self.add_setting(config, ["namelist:convection", "prog_ent_min"], "0.5")
        self.add_setting(config, ["namelist:convection", "qlmin"], "4.0e-4")
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-cloud
        cvscheme = self.get_setting_value(
            config, ["namelist:convection", "cv_scheme"]
        )
        if cvscheme == "'comorph'":
            self.add_setting(
                config, ["namelist:cloud", "fsd_conv_const"], "3.0"
            )
            self.add_setting(
                config, ["namelist:cloud", "fsd_min_conv_frac"], "0.02"
            )
            self.add_setting(
                config, ["namelist:cloud", "fsd_nonconv_const"], "0.8"
            )
        else:
            self.add_setting(
                config, ["namelist:cloud", "fsd_conv_const"], "2.81"
            )
            self.add_setting(
                config, ["namelist:cloud", "fsd_min_conv_frac"], "0.0"
            )
            self.add_setting(
                config, ["namelist:cloud", "fsd_nonconv_const"], "1.14"
            )
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-boundary_layer
        self.add_setting(config, ["namelist:blayer", "dec_thres_cu"], "0.05")
        # Commands From: interfaces/physics_schemes_interface/rose-meta/um-aerosol
        self.add_setting(config, ["namelist:aerosol", "horiz_d"], "2.25")
        self.add_setting(config, ["namelist:aerosol", "us_am"], "1.45")
        return config, self.reports


class vn20_t249(MacroUpgrade):
    """Upgrade macro for ticket #249 by Denis Sergeev."""

    BEFORE_TAG = "vn2.0_t334"
    AFTER_TAG = "vn2.0_t249"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/jules_interface/rose-meta/jules-lfric
        nml = "namelist:specified_surface"
        self.add_setting(config, [nml, "surf_temp_forcing"], "'none'")
        self.add_setting(config, [nml, "internal_flux_method"], "'uniform'")
        self.add_setting(config, [nml, "internal_flux_value"], "0.0")
        # Commands From: science/gungho/rose-meta/lfric-gungho
        nml = "namelist:files"
        self.add_setting(config, [nml, "internal_flux_ancil_path"], "''")
        return config, self.reports


class vn20_t547(MacroUpgrade):
    """Upgrade macro for ticket #547 by Adrian Lock."""

    BEFORE_TAG = "vn2.0_t249"
    AFTER_TAG = "vn2.0_t547"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-microphysics
        nml = "namelist:microphysics"
        self.add_setting(config, [nml, "heavy_rain_evap_fac"], "0.0")
        return config, self.reports


class vn20_t429(MacroUpgrade):
    """Upgrade macro for ticket #429 by Denis Sergeev."""

    BEFORE_TAG = "vn2.0_t547"
    AFTER_TAG = "vn2.0_t429"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        # Add the Burrows & Sharp (1999) chemistry scheme option
        # (enabled if chem_scheme='flexchem')
        nml = "namelist:chemistry"
        self.add_setting(config, [nml, "flexchem_opt"], "'bs1999'")
        # Commands From: rose-meta/socrates-radiation
        # Add more radiatively active gases: Cs, K, Li, Na, Rb, TiO, VO
        nml = "namelist:radiative_gases"
        self.add_setting(config, [nml, "cs_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "cs_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "cs_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "cs_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "cs_mix_ratio"], "0")
        self.add_setting(config, [nml, "cs_rad_opt"], "'off'")
        self.add_setting(config, [nml, "k_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "k_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "k_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "k_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "k_mix_ratio"], "0")
        self.add_setting(config, [nml, "k_rad_opt"], "'off'")
        self.add_setting(config, [nml, "li_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "li_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "li_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "li_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "li_mix_ratio"], "0")
        self.add_setting(config, [nml, "li_rad_opt"], "'off'")
        self.add_setting(config, [nml, "na_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "na_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "na_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "na_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "na_mix_ratio"], "0")
        self.add_setting(config, [nml, "na_rad_opt"], "'off'")
        self.add_setting(config, [nml, "rb_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "rb_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "rb_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "rb_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "rb_mix_ratio"], "0")
        self.add_setting(config, [nml, "rb_rad_opt"], "'off'")
        self.add_setting(config, [nml, "tio_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "tio_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "tio_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "tio_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "tio_mix_ratio"], "0")
        self.add_setting(config, [nml, "tio_rad_opt"], "'off'")
        self.add_setting(config, [nml, "vo_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "vo_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "vo_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "vo_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "vo_mix_ratio"], "0")
        self.add_setting(config, [nml, "vo_rad_opt"], "'off'")
        return config, self.reports


class vn20_t552(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t429"
    AFTER_TAG = "vn2.0_t552"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-convection
        nml = "namelist:convection"
        self.add_setting(config, [nml, "resdep_precipramp"], ".false.")
        self.add_setting(config, [nml, "dx_ref"], "50000.0")
        return config, self.reports


class vn20_t562(MacroUpgrade):
    """Upgrade macro for ticket #562 by Paul Burns."""

    BEFORE_TAG = "vn2.0_t552"
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
        # Commands From: rose-meta/um-convection
        nml = "namelist:convection"
        self.add_setting(config, [nml, "l_cvdiag_ctop_qmax"], ".false.")
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


class vn20_t84(MacroUpgrade):
    """Upgrade macro for ticket #84 by Chris Smith."""

    BEFORE_TAG = "vn2.0_t472"
    AFTER_TAG = "vn2.0_t84"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Add new wind_forcing namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = source + "\n" + " (namelist:wind_forcing)"
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add default data for wind_forcing namelist"""
        self.add_setting(config, ["namelist:wind_forcing"])
        self.add_setting(
            config, ["namelist:wind_forcing", "coordinate"], "'height'"
        )
        self.add_setting(config, ["namelist:wind_forcing", "heights"], "0.0")
        self.add_setting(
            config, ["namelist:wind_forcing", "number_heights"], "1"
        )
        self.add_setting(config, ["namelist:wind_forcing", "number_times"], "1")
        self.add_setting(
            config, ["namelist:wind_forcing", "profile_data_u"], "0.0"
        )
        self.add_setting(
            config, ["namelist:wind_forcing", "profile_data_v"], "0.0"
        )
        self.add_setting(config, ["namelist:wind_forcing", "times"], "0.0")
        return config, self.reports


class vn20_t588(MacroUpgrade):
    """Upgrade macro for ticket #588 by Maggie Hendry."""

    BEFORE_TAG = "vn2.0_t84"
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
        # Commands From: rose-meta/um-stochastic_physics
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-spectral_gwd
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-orographic_drag
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-microphysics
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-iau
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-convection
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-cloud
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-chemistry
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-boundary_layer
        # Blank Upgrade Macro
        # Commands From: rose-meta/um-aerosol
        # Blank Upgrade Macro
        # Commands From: rose-meta/socrates-radiation
        # Blank Upgrade Macro
        # Commands From: rose-meta/jules-lfric
        # Blank Upgrade Macro
        # Commands From: rose-meta/lfric-driver
        # Blank Upgrade Macro
        # Commands From: rose-meta/lfric-gungho
        # Blank Upgrade Macro
        return config, self.reports
