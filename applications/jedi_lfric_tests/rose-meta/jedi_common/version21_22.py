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


class vn21_t83(MacroUpgrade):
    """Upgrade macro for ticket #83 by Chris Smith."""

    BEFORE_TAG = "vn2.1"
    AFTER_TAG = "vn2.1_t83"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Add theta_relax namelist to configuration source list"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"\(namelist:temp_tend_data\)",
            r"(namelist:temp_tend_data)" + "\n" + " (namelist:theta_relax)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add theta_relaxation setting to external_forcing namelist"""
        self.add_setting(
            config, ["namelist:external_forcing", "theta_relaxation"], ".false."
        )
        """Data for theta_relax namelist"""
        self.add_setting(config, ["namelist:theta_relax"])
        self.add_setting(
            config, ["namelist:theta_relax", "coordinate"], "'height'"
        )
        self.add_setting(config, ["namelist:theta_relax", "heights"], "0.0")
        self.add_setting(
            config, ["namelist:theta_relax", "number_heights"], "1"
        )
        self.add_setting(config, ["namelist:theta_relax", "number_times"], "1")
        self.add_setting(
            config, ["namelist:theta_relax", "profile_data"], "0.0"
        )
        self.add_setting(config, ["namelist:theta_relax", "times"], "0.0")
        self.add_setting(config, ["namelist:theta_relax", "timescale"], "1.0")
        return config, self.reports


class vn21_t476(MacroUpgrade):
    """Upgrade macro for ticket #476 by Dave Case."""

    BEFORE_TAG = "vn2.1_t83"
    AFTER_TAG = "vn2.1_t476"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Set configure_segments in physics namelist to be false to ignore segmentation"""
        self.add_setting(config, ["namelist:physics", "bl_segment"], "0")
        self.add_setting(config, ["namelist:physics", "gw_segment"], "0")
        self.add_setting(config, ["namelist:physics", "ussp_segment"], "0")
        return config, self.reports


class vn21_t663(MacroUpgrade):
    """Upgrade macro for ticket #663 by Ian Boutle."""

    BEFORE_TAG = "vn2.1_t476"
    AFTER_TAG = "vn2.1_t663"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-aerosol
        # Commands From: rose-meta/um-aerosol
        self.add_setting(config, ["namelist:aerosol", "murk_lbc"], ".false.")
        return config, self.reports


class vn21_t708(MacroUpgrade):
    """Upgrade macro for ticket #708 by mark Hedley."""

    BEFORE_TAG = "vn2.1_t663"
    AFTER_TAG = "vn2.1_t708"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(
            config, ["namelist:logging", "log_to_rank_zero_only"], ".false."
        )
        return config, self.reports


class vn21_t164(MacroUpgrade):
    """Upgrade macro for ticket #164 by annemccabe."""

    BEFORE_TAG = "vn2.1_t708"
    AFTER_TAG = "vn2.1_t164"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-stochastic_physics
        """Add new namelist settings to stochastic_physics"""
        alnir = self.get_setting_value(
            config, ["namelist:jules_pftparm", "alnir_io"]
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alnir"], alnir
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alnir_min"], alnir
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alnir_max"], alnir
        )
        alpar = self.get_setting_value(
            config, ["namelist:jules_pftparm", "alpar_io"]
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alpar"], alpar
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alpar_min"], alpar
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_alpar_max"], alpar
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_mp_fxd_cld_num"],
            "150.0e+06,150.0e+06,150.0e+06",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_mp_ice_fspd"],
            "6000000.0,6000000.0,6000000.0",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_mp_mp_czero"],
            "10.0,10.0,10.0",
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_mp_mpof"], "0.5,0.5,0.5"
        )
        omega = self.get_setting_value(
            config, ["namelist:jules_pftparm", "omega_io"]
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omega"], omega
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omega_min"], omega
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omega_max"], omega
        )
        omnir = self.get_setting_value(
            config, ["namelist:jules_pftparm", "omnir_io"]
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omnir"], omnir
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omnir_min"], omnir
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_omnir_max"], omnir
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_lsfc_orog_drag_param"],
            "0.15,0.15,0.15",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_mp_snow_fspd"],
            "12.0,12.0,12.0",
        )
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_lsfc_z0_urban_mult"],
            "1.0,1.0,1.0",
        )
        soil = self.get_setting_value(
            config, ["namelist:jules_surface_types", "soil"]
        )
        npft = self.get_setting_value(
            config, ["namelist:jules_surface_types", "npft"]
        )
        z0_nvg_io = self.get_setting_value(
            config, ["namelist:jules_nvegparm", "z0_nvg_io"]
        )
        z0_soil = z0_nvg_io.split(",")[int(soil) - int(npft) - 1]
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_lsfc_z0_soil"],
            z0_soil + "," + z0_soil + "," + z0_soil,
        )
        z0hm_nvg_io = self.get_setting_value(
            config, ["namelist:jules_nvegparm", "z0hm_nvg_io"]
        )
        z0hm_soil = z0hm_nvg_io.split(",")[int(soil) - int(npft) - 1]
        self.add_setting(
            config,
            ["namelist:stochastic_physics", "rp_lsfc_z0hm_soil"],
            z0hm_soil + "," + z0hm_soil + "," + z0hm_soil,
        )
        z0v_io = self.get_setting_value(config, ["namelist:surface", "z0v"])
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_z0v"], z0v_io
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_z0v_min"], z0v_io
        )
        self.add_setting(
            config, ["namelist:stochastic_physics", "rp_lsfc_z0v_max"], z0v_io
        )
        return config, self.reports


class vn21_t657(MacroUpgrade):
    """Upgrade macro for ticket #657 by Ian Boutle."""

    BEFORE_TAG = "vn2.1_t164"
    AFTER_TAG = "vn2.1_t657"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        nml = "namelist:cloud"
        fsd_nonconv = self.get_setting_value(config, [nml, "fsd_nonconv_const"])
        self.remove_setting(config, [nml, "fsd_nonconv_const"])
        self.add_setting(config, [nml, "fsd_nonconv_ice_const"], fsd_nonconv)
        self.add_setting(config, [nml, "fsd_nonconv_liq_const"], fsd_nonconv)
        return config, self.reports


class vn21_t596(MacroUpgrade):
    """Upgrade macro for ticket #596 by Maggie Hendry."""

    BEFORE_TAG = "vn2.1_t657"
    AFTER_TAG = "vn2.1_t596"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Only used in LFRic apps
        self.rename_setting(
            config,
            ["namelist:jules_surface", "check_soilm_negatives"],
            ["namelist:surface", "check_soilm_negatives"],
        )
        self.rename_setting(
            config,
            ["namelist:jules_surface", "lake_water_conservation"],
            ["namelist:surface", "lake_water_conservation"],
        )
        # jules-sea-seaice from jules-lfric being shared rather than um-atmos
        self.rename_setting(
            config,
            ["namelist:surface", "amip_ice_thick"],
            ["namelist:jules_sea_seaice", "amip_ice_thick"],
        )
        self.rename_setting(
            config,
            ["namelist:surface", "z0h_specified"],
            ["namelist:jules_sea_seaice", "z0h_specified"],
        )
        self.rename_setting(
            config,
            ["namelist:surface", "z0m_specified"],
            ["namelist:jules_sea_seaice", "z0m_specified"],
        )
        return config, self.reports


class vn21_t742(MacroUpgrade):
    """Upgrade macro for ticket #742 by Shusuke Nishimoto."""

    BEFORE_TAG = "vn2.1_t596"
    AFTER_TAG = "vn2.1_t742"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-boundary_layer
        cv_scheme = self.get_setting_value(
            config, ["namelist:convection", "cv_scheme"]
        )
        if cv_scheme == "'lambert_lewis'":
            bl_res_inv = "'off'"
        else:
            bl_res_inv = "'cosine_inv_flux'"
        self.add_setting(config, ["namelist:blayer", "bl_res_inv"], bl_res_inv)
        return config, self.reports


class vn21_t4604(MacroUpgrade):
    """Upgrade macro for ticket #4604 by Mike Hobson."""

    BEFORE_TAG = "vn2.1_t742"
    AFTER_TAG = "vn2.1_t4604"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        """Rename generate_inner_haloes to generate_inner_halos"""
        self.rename_setting(
            config,
            ["namelist:partitioning", "generate_inner_haloes"],
            ["namelist:partitioning", "generate_inner_halos"],
        )
        return config, self.reports


class vn21_t208(MacroUpgrade):
    """Upgrade macro for ticket #208 by Thomas Bendall."""

    BEFORE_TAG = "vn2.1_t4604"
    AFTER_TAG = "vn2.1_t208"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """
        Add transport_overwrite_freq to namelist boundaries. This takes the
        default value of 'final' unless the Method-of-Lines transport scheme is
        being used for a limited area model, in which case the value is set to
        'split_step'.
        """
        mol_transport = self.get_setting_value(
            config, ["namelist:transport", "horizontal_method"]
        )
        if "1" in mol_transport:
            self.add_setting(
                config,
                ["namelist:boundaries", "transport_overwrite_freq"],
                "'split_step'",
            )
        else:
            self.add_setting(
                config,
                ["namelist:boundaries", "transport_overwrite_freq"],
                "'final'",
            )
        return config, self.reports


class vn21_t255(MacroUpgrade):
    """Upgrade macro for ticket #255 by Christine Johnson."""

    BEFORE_TAG = "vn2.1_t208"
    AFTER_TAG = "vn2.1_t255"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Blank Upgrade Macro
        return config, self.reports


class vn21_t744(MacroUpgrade):
    """Upgrade macro for ticket #744 by Ian Boutle."""

    BEFORE_TAG = "vn2.1_t255"
    AFTER_TAG = "vn2.1_t744"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        nml = "namelist:cloud"
        liq_fsd = self.get_setting_value(config, [nml, "cloud_horizontal_fsd"])
        self.remove_setting(config, [nml, "cloud_horizontal_fsd"])
        self.add_setting(config, [nml, "cloud_horizontal_liq_fsd"], liq_fsd)
        self.add_setting(config, [nml, "cloud_horizontal_ice_fsd"], "0.0")
        return config, self.reports


class vn21_t590(MacroUpgrade):
    """Upgrade macro for ticket #590 by Thomas Bendall."""

    BEFORE_TAG = "vn2.1_t744"
    AFTER_TAG = "vn2.1_t590"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """
        Add share_stencil_extent to namelist departure_points, and set the
        default value to be true.
        """
        self.add_setting(
            config,
            ["namelist:departure_points", "share_stencil_extent"],
            ".true.",
        )
        return config, self.reports


class vn21_t638(MacroUpgrade):
    """Upgrade macro for ticket #638 by Tom Hill."""

    BEFORE_TAG = "vn2.1_t590"
    AFTER_TAG = "vn2.1_t638"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jedi_common
        # Blank Upgrade Macro
        return config, self.reports


class vn21_t207(MacroUpgrade):
    """Upgrade macro for ticket #207 by Samantha Pullen."""

    BEFORE_TAG = "vn2.1_t638"
    AFTER_TAG = "vn2.1_t207"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-iau
        """Add new variables and default data to iau namelist"""
        self.add_setting(
            config, ["namelist:iau", "iau_ainc_multifile"], ".false."
        )
        self.add_setting(
            config, ["namelist:iau", "iau_tendency_addinf"], ".true."
        )
        self.add_setting(
            config, ["namelist:iau", "iau_tendency_ainc"], ".false."
        )
        self.add_setting(
            config, ["namelist:iau", "iau_tendency_bcorr"], ".true."
        )
        self.add_setting(
            config, ["namelist:iau", "iau_tendency_pertinc"], ".false."
        )
        self.add_setting(config, ["namelist:iau", "iau_use_addinf"], ".false.")
        self.add_setting(config, ["namelist:iau", "iau_use_bcorr"], ".false.")
        self.add_setting(config, ["namelist:iau", "iau_use_pertinc"], ".false.")
        # Commands From: rose-meta/lfric-gungho
        """Add new iau_addinf_io, iau_ainc_io, and iau_bcorr_io namelists"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:helmholtz_solver",
            r"namelist:helmholtz_solver"
            + "\n"
            + " (namelist:iau_addinf_io(:))"
            + "\n"
            + " (namelist:iau_ainc_io(:))"
            + "\n"
            + " (namelist:iau_bcorr_io(:))",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add iau_addinf_path setting to files namelist"""
        self.add_setting(config, ["namelist:files", "iau_addinf_path"], "''")
        """Add iau_bcorr_path setting to files namelist"""
        self.add_setting(config, ["namelist:files", "iau_bcorr_path"], "''")
        """Add iau_pert_path setting to files namelist"""
        self.add_setting(config, ["namelist:files", "iau_pert_path"], "''")
        """Add default data for iau_addinf_io namelist"""
        self.add_setting(
            config, ["namelist:iau_addinf_io(addinf1)", "filename"], "''"
        )
        self.add_setting(
            config, ["namelist:iau_addinf_io(addinf1)", "start_time"], "0"
        )
        self.add_setting(
            config, ["namelist:iau_addinf_io(addinf2)", "filename"], "''"
        )
        self.add_setting(
            config, ["namelist:iau_addinf_io(addinf2)", "start_time"], "0"
        )
        """Add default data for iau_ainc_io namelist"""
        self.add_setting(
            config, ["namelist:iau_ainc_io(ainc1)", "filename"], "''"
        )
        self.add_setting(
            config, ["namelist:iau_ainc_io(ainc1)", "start_time"], "0"
        )
        self.add_setting(
            config, ["namelist:iau_ainc_io(ainc2)", "filename"], "''"
        )
        self.add_setting(
            config, ["namelist:iau_ainc_io(ainc2)", "start_time"], "0"
        )
        """Add default data for iau_bcorr_io namelist"""
        self.add_setting(
            config, ["namelist:iau_bcorr_io(bcorr1)", "filename"], "''"
        )
        self.add_setting(
            config, ["namelist:iau_bcorr_io(bcorr1)", "start_time"], "0"
        )
        """Add multifile_io setting to io namelist"""
        self.add_setting(config, ["namelist:io", "multifile_io"], ".false.")
        return config, self.reports


class vn21_t756(MacroUpgrade):
    """Upgrade macro for ticket #756 by Adrian Lock."""

    BEFORE_TAG = "vn2.1_t207"
    AFTER_TAG = "vn2.1_t756"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/socrates-radiation
        self.add_setting(
            config, ["namelist:radiative_gases", "o3_profile_data"], "0"
        )
        self.add_setting(
            config, ["namelist:radiative_gases", "o3_profile_heights"], "0.0"
        )
        self.add_setting(
            config, ["namelist:radiative_gases", "o3_profile_size"], "0"
        )
        # Commands From: rose-meta/lfric-gungho
        # Commands From: rose-meta/lfric-gungho
        """Add vertadvect namelist to configuration source list"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"\(namelist:vapour_forcing\)",
            r"(namelist:vapour_forcing)" + "\n" + " (namelist:vertadvect)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add vertadvect_forcing setting to external_forcing namelist"""
        self.add_setting(
            config,
            ["namelist:external_forcing", "vertadvect_forcing"],
            ".false.",
        )
        """Data for vertadvect namelist"""
        self.add_setting(config, ["namelist:vertadvect"])
        self.add_setting(config, ["namelist:vertadvect", "heights"], "0.0")
        self.add_setting(config, ["namelist:vertadvect", "number_heights"], "1")
        self.add_setting(config, ["namelist:vertadvect", "number_times"], "1")
        self.add_setting(
            config, ["namelist:vertadvect", "profile_data_w"], "0.0"
        )
        self.add_setting(config, ["namelist:vertadvect", "times"], "0.0")
        return config, self.reports


class vn21_t871(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.1_t756"
    AFTER_TAG = "vn2.2"

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
        # Commands From: rose-meta/lfric-linear
        # Blank Upgrade Macro
        # Commands From: rose-meta/lfric-adjoint
        # Blank Upgrade Macro
        # Commands From: rose-meta/jedi_common
        # Blank Upgrade Macro
        return config, self.reports
