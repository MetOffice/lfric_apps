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


class vn22_t885(MacroUpgrade):
    """Upgrade macro for ticket #885 by Samantha Pullen."""

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t885"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Add iau_sst_path to files namelist"""
        self.add_setting(config, ["namelist:files", "iau_sst_path"], "''")
        """Add iau_sst to section_choice namelist"""
        self.add_setting(
            config, ["namelist:section_choice", "iau_sst"], ".false."
        )
        return config, self.reports


class vn22_t4661(MacroUpgrade):
    """Upgrade macro for ticket #4661 by Denis Sergeev."""

    BEFORE_TAG = "vn2.2_t885"
    AFTER_TAG = "vn2.2_t4661"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(config, ["namelist:extrusion", "eta_values"], "''")
        return config, self.reports


class vn22_t771(MacroUpgrade):
    """Upgrade macro for ticket #771 by josephwallwork."""

    BEFORE_TAG = "vn2.2_t4661"
    AFTER_TAG = "vn2.2_t771"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        """Add new namelist options for chemistry timestep halving"""
        self.add_setting(
            config,
            ["namelist:chemistry", "i_chem_timestep_halvings"],
            value="0",
        )
        return config, self.reports


class vn22_t887(MacroUpgrade):
    """Upgrade macro for ticket #887 by Mike Whitall."""

    BEFORE_TAG = "vn2.2_t771"
    AFTER_TAG = "vn2.2_t887"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        nml = "namelist:cloud"
        self.add_setting(config, [nml, "dbsdtbs_turb_0"], "1.5E-4")
        self.add_setting(config, [nml, "i_pc2_erosion_numerics"], "'implicit'")
        return config, self.reports


class vn22_t886(MacroUpgrade):
    """Upgrade macro for ticket #886 by Samantha Pullen."""

    BEFORE_TAG = "vn2.2_t887"
    AFTER_TAG = "vn2.2_t886"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-iau
        # Blank Upgrade Macro
        return config, self.reports


class vn22_t850(MacroUpgrade):
    """Upgrade macro for ticket #850 by Shusuke Nishimoto."""

    BEFORE_TAG = "vn2.2_t886"
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


class vn22_t36(MacroUpgrade):
    """Upgrade macro for ticket #36 by Thomas Bendall."""

    BEFORE_TAG = "vn2.2_t850"
    AFTER_TAG = "vn2.2_t36"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """
        Reorganises the transport options that describe treatment of
        cubed-sphere panel edges.
        Replace all "special edges" options with "remapping".
        """
        # Get values of old options
        nml = "namelist:transport"
        special_edges = self.get_setting_value(
            config, [nml, "special_edges_treatment"]
        )
        extended = self.get_setting_value(config, [nml, "extended_mesh"])
        # Work out the new option for "panel_edge_treatment"
        if special_edges == ".true.":
            panel_edge_treatment = "'special_edges'"
            self.add_setting(config, [nml, "panel_edge_high_order"], ".true.")
        elif extended == ".true.":
            panel_edge_treatment = "'extended_mesh'"
            self.add_setting(config, [nml, "panel_edge_high_order"], ".false.")
        else:
            panel_edge_treatment = "'none'"
            self.add_setting(config, [nml, "panel_edge_high_order"], ".true.")
        # Add the new option and remove the old ones
        self.remove_setting(config, [nml, "extended_mesh"])
        self.remove_setting(config, [nml, "special_edges_treatment"])
        self.remove_setting(config, [nml, "special_edges_high_order"])
        self.add_setting(
            config, [nml, "panel_edge_treatment"], panel_edge_treatment
        )
        return config, self.reports


class vn22_t797(MacroUpgrade):
    """Upgrade macro for ticket #797 by Charlotte Norris."""

    BEFORE_TAG = "vn2.2_t36"
    AFTER_TAG = "vn2.2_t797"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        """
        Add fastjx_numwavel, fjx_solcyc_type, fjx_scat_file, fjx_solar_file,
        fastjx_dir, fjx_spec_file to namelist chemistry
        """
        self.add_setting(
            config, ["namelist:chemistry", "fastjx_numwavel"], "18"
        )
        self.add_setting(config, ["namelist:chemistry", "fjx_solcyc_type"], "0")
        self.add_setting(
            config, ["namelist:chemistry", "fjx_scat_file"], "'FJX_scat.dat'"
        )
        self.add_setting(
            config,
            ["namelist:chemistry", "fjx_solar_file"],
            "'FJX_solcyc_May17.dat'",
        )
        self.add_setting(
            config,
            ["namelist:chemistry", "fjx_spec_file"],
            "'FJX_spec_Nov11.dat'",
        )
        self.add_setting(
            config,
            ["namelist:chemistry", "fastjx_dir"],
            "'$UMDIR/vn13.9/ctldata/UKCA/fastj'",
        )
        return config, self.reports


class vn22_t995(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.2_t797"
    AFTER_TAG = "vn2.2_t995"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        self.add_setting(config, ["namelist:mixing", "smag_l_calc"], "'UseDx'")
        return config, self.reports


class vn22_t202(MacroUpgrade):
    """Upgrade macro for ticket #202 by Katty Huang."""

    BEFORE_TAG = "vn2.2_t995"
    AFTER_TAG = "vn2.2_t202"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        self.add_setting(
            config, ["namelist:jules_surface", "anthrop_heat_mean"], "20.0"
        )
        self.add_setting(
            config, ["namelist:jules_surface", "anthrop_heat_option"], "'dukes'"
        )
        return config, self.reports


class vn22_t827(MacroUpgrade):
    """Upgrade macro for ticket #827 by Thomas Bendall."""

    BEFORE_TAG = "vn2.2_t202"
    AFTER_TAG = "vn2.2_t827"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """
        Adds the "theta_moist_source" variable to the formulation namelist,
        which controls the option to add the missing term to the potential
        temperature equation, relating to the different heat capacities of moist
        phases. This is set to be trig-ignored when moisture settings are not
        'traditional'. In both cases the default value is .false.
        """
        nml = "namelist:formulation"
        moisture = self.get_setting_value(config, [nml, "moisture_formulation"])
        default_setting = ".false."
        if moisture == "'traditional'":
            self.add_setting(
                config, [nml, "theta_moist_source"], default_setting
            )
        else:
            # Trig-ignored as moisture not being used
            self.add_setting(
                config, [nml, "!!theta_moist_source"], default_setting
            )
        return config, self.reports


class vn22_t938(MacroUpgrade):
    """Upgrade macro for ticket #938 by Jon Elsey."""

    BEFORE_TAG = "vn2.2_t827"
    AFTER_TAG = "vn2.2_t938"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-aerosol
        # Add setting for ukca_mode_segment_size
        self.add_setting(
            config, ["namelist:aerosol", "ukca_mode_seg_size"], value="4"
        )
        return config, self.reports


class vn22_t903(MacroUpgrade):
    """Upgrade macro for ticket #903 by James Kent."""

    BEFORE_TAG = "vn2.2_t938"
    AFTER_TAG = "vn2.2_t903"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        # Add horizontal predictor options to formulation namelist
        self.add_setting(
            config,
            ["namelist:formulation", "horizontal_physics_predictor"],
            ".false.",
        )
        self.add_setting(
            config,
            ["namelist:formulation", "horizontal_transport_predictor"],
            ".false.",
        )
        return config, self.reports


class vn22_t1012(MacroUpgrade):
    """Upgrade macro for ticket #1012 by Maggie Hendry."""

    BEFORE_TAG = "vn2.2_t903"
    AFTER_TAG = "vn2.2_t1012"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Blank upgrade macro to bump tag
        return config, self.reports


class vn22_t953(MacroUpgrade):
    """Upgrade macro for ticket #593 by Thomas Bendall."""

    BEFORE_TAG = "vn2.2_t1012"
    AFTER_TAG = "vn2.2_t953"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """
        Add "wind_mono_top" and "wind_mono_top_depth" to the transport namelist.
        The "wind_mono_top" option is set to false by default.
        """
        num_layers = self.get_setting_value(
            config, ["namelist:extrusion", "number_of_layers"]
        )
        wind_mono_top = ".false."
        # Add wind_mono_top setting
        self.add_setting(
            config, ["namelist:transport", "wind_mono_top"], wind_mono_top
        )
        # Add wind_mono_top_depth setting
        # If the number of layers is greater than 10, set the depth to 5
        # otherwise set it to 0
        if int(num_layers) > 10:
            wind_mono_top_depth = "5"
        else:
            wind_mono_top_depth = "0"
        self.add_setting(
            config,
            ["namelist:transport", "wind_mono_top_depth"],
            wind_mono_top_depth,
        )
        return config, self.reports


class vn22_t4020(MacroUpgrade):
    """Upgrade macro for ticket #4020 by Andrew Coughtrie."""

    BEFORE_TAG = "vn2.2_t953"
    AFTER_TAG = "vn2.2_t4020"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(
            config, ["namelist:io", "end_of_run_checkpoint"], ".true."
        )
        self.add_setting(config, ["namelist:io", "checkpoint_times"], "")
        return config, self.reports


class vn22_t618(MacroUpgrade):
    """Upgrade macro for ticket #618 by Mohit Dalvi."""

    BEFORE_TAG = "vn2.2_t4020"
    AFTER_TAG = "vn2.2_t618"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        """
        Add Photolysis scheme and FastJX namelist items to chemistry namelist.
        Value of photol_scheme is based on chem_scheme choice.
        Change the value of fastjx_dir to replace hardwired path with
        link from central ctldata extract
        """
        chem_scheme = self.get_setting_value(
            config, ["namelist:chemistry", "chem_scheme"]
        )
        if chem_scheme == "'strattrop'":
            phot_scheme = "'fastjx'"
        else:
            phot_scheme = "'off'"  # Default
        self.add_setting(
            config, ["namelist:chemistry", "photol_scheme"], phot_scheme
        )
        self.add_setting(
            config, ["namelist:chemistry", "chem_timestep"], "3600"
        )
        self.add_setting(config, ["namelist:chemistry", "fastjx_mode"], "1")
        self.add_setting(
            config, ["namelist:chemistry", "fastjx_prescutoff"], "20.0"
        )
        self.add_setting(
            config, ["namelist:chemistry", "fjx_solcyc_months"], "0"
        )
        self.add_setting(config, ["namelist:chemistry", "fjx_solcyc_type"], "0")
        self.change_setting_value(
            config, ["namelist:chemistry", "fastjx_dir"], "'fastj'"
        )
        return config, self.reports


class vn22_t34(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.2_t618"
    AFTER_TAG = "vn3.0"

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
        # Commands From: rose-meta/jules-lsm
        # Blank Upgrade Macro
        # Commands From: rose-meta/lfric-driver
        # Blank Upgrade Macro
        # Commands From: rose-meta/lfric-gungho
        # Blank Upgrade Macro
        # Commands From: rose-meta/lfric-lfric_coupled
        # Blank Upgrade Macro
        return config, self.reports
