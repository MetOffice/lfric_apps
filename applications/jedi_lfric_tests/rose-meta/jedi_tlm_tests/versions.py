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


class vn22_t987(MacroUpgrade):
    """Upgrade macro for ticket #987 by Christine Johnson."""

    BEFORE_TAG = "vn2.2_t887"
    AFTER_TAG = "vn2.2_t987"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-linear
        scaling = self.get_setting_value(
            config, ["namelist:planet", "scaling_factor"]
        )
        if "125.0" in scaling:
            self.add_setting(config, ["namelist:linear", "fixed_ls"], ".false.")
        else:
            self.add_setting(config, ["namelist:linear", "fixed_ls"], ".true.")

        return config, self.reports


class vn22_t886(MacroUpgrade):
    """Upgrade macro for ticket #886 by Samantha Pullen."""

    BEFORE_TAG = "vn2.2_t987"
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
