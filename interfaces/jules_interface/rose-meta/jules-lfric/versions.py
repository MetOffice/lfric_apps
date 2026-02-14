import re
import sys

from metomi.rose.upgrade import MacroUpgrade

from .version22_30 import *


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


class vn30_t146(MacroUpgrade):
    """Upgrade macro for ticket #146 by Maggie Hendry."""

    BEFORE_TAG = "vn3.0"
    AFTER_TAG = "vn3.0_t146"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/jules-lfric
        # Add jules_model_environment_lfric namelist
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:jules_hydrology",
            r"namelist:jules_hydrology)"
            + "\n"
            + " (namelist:jules_model_environment_lfric",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        self.add_setting(
            config,
            ["namelist:jules_model_environment_lfric", "l_jules_parent"],
            "'lfric'",
        )
        # Add jules_surface namelist items
        self.add_setting(
            config,
            ["namelist:jules_surface", "all_tiles"],
            "'off'",
        )
        self.add_setting(config, ["namelist:jules_surface", "beta1"], "0.83")
        self.add_setting(config, ["namelist:jules_surface", "beta2"], "0.93")
        self.add_setting(
            config, ["namelist:jules_surface", "beta_cnv_bl"], "0.04"
        )
        self.add_setting(
            config,
            ["namelist:jules_surface", "fd_hill_option"],
            "'capped_lowhill'",
        )
        self.add_setting(config, ["namelist:jules_surface", "fwe_c3"], "0.5")
        self.add_setting(
            config, ["namelist:jules_surface", "fwe_c4"], "20000.0"
        )
        self.add_setting(config, ["namelist:jules_surface", "hleaf"], "5.7e4")
        self.add_setting(config, ["namelist:jules_surface", "hwood"], "1.1e4")
        self.add_setting(
            config, ["namelist:jules_surface", "i_modiscopt"], "'on'"
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_epot_corr"], ".true."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_land_ice_imp"], ".true."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_mo_buoyancy_calc"], ".true."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "orog_drag_param"], "0.15"
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_flake_model"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_elev_land_ice"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_elev_lw_down"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_surface", "l_point_data"], ".false."
        )

        return config, self.reports


class vn30_t205(MacroUpgrade):
    # Upgrade macro for #205 by Maggie Hendry

    BEFORE_TAG = "vn3.0_t146"
    AFTER_TAG = "vn3.0_t205"

    def upgrade(self, config, meta_config=None):
        # Add jules_pftparm items hard-wired in jules_physics_init
        jules_pftparm = {}
        jules_pftparm["a_wl_io"] = "0.65, 0.65, 0.005, 0.005, 0.10"
        jules_pftparm["a_ws_io"] = "10.0, 10.0, 1.0, 1.0, 10.0"
        jules_pftparm["albsnc_min_io"] = (
            "3.0e-1, 3.0e-1, 8.0e-1, 8.0e-1, 8.0e-1"
        )
        jules_pftparm["albsnf_maxl_io"] = "0.095, 0.059, 0.128, 0.106, 0.077"
        jules_pftparm["albsnf_maxu_io"] = "0.215, 0.132, 0.288, 0.239, 0.173"
        jules_pftparm["alnirl_io"] = "0.30, 0.23, 0.39, 0.39, 0.39"
        jules_pftparm["alniru_io"] = "0.75, 0.65, 0.95, 0.95, 0.87"
        jules_pftparm["alparl_io"] = "0.06, 0.04, 0.06, 0.06, 0.06"
        jules_pftparm["alparu_io"] = "0.15, 0.11, 0.25, 0.25, 0.25"
        jules_pftparm["alpha_io"] = "0.08, 0.08, 0.08, 0.04, 0.08"
        jules_pftparm["b_wl_io"] = "1.667, 1.667, 1.667, 1.667, 1.667"
        jules_pftparm["c3_io"] = "3*1, 0, 1"
        jules_pftparm["can_struct_a_io"] = "1.0, 1.0, 1.0, 1.0, 1.0"
        jules_pftparm["dgl_dm_io"] = "0.0, 0.0, 0.0, 0.0, 0.0"
        jules_pftparm["dgl_dt_io"] = "9.0, 9.0, 0.0, 0.0, 9.0"
        jules_pftparm["dqcrit_io"] = "0.090, 0.060, 0.100, 0.075, 0.100"
        jules_pftparm["dust_veg_scj_io"] = "0.0, 0.0, 1.0, 1.0, 0.5"
        jules_pftparm["dz0v_dh_io"] = "5.0e-2, 5.0e-2, 1.0e-1, 1.0e-1, 1.0e-1"
        jules_pftparm["emis_pft_io"] = "0.98, 0.99, 0.98, 0.98, 0.98"
        jules_pftparm["eta_sl_io"] = "0.01, 0.01, 0.01, 0.01, 0.01"
        jules_pftparm["f0_io"] = "0.875, 0.875, 0.900, 0.800, 0.900"
        jules_pftparm["fd_io"] = "0.015, 0.015, 0.015, 0.025, 0.015"
        jules_pftparm["fsmc_of_io"] = "0.0, 0.0, 0.0, 0.0, 0.0"
        jules_pftparm["g_leaf_0_io"] = "0.25, 0.25, 0.25, 0.25, 0.25"
        jules_pftparm["glmin_io"] = "1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6"
        jules_pftparm["gsoil_f_io"] = "1.0, 1.0, 1.0, 1.0, 1.0"
        jules_pftparm["hw_sw_io"] = "0.5, 0.5, 0.5, 0.5, 0.5"
        jules_pftparm["infil_f_io"] = "4.0, 4.0, 2.0, 2.0, 2.0"
        jules_pftparm["kn_io"] = "0.78, 0.78, 0.78, 0.78, 0.78"
        jules_pftparm["kpar_io"] = "0.5, 0.5, 0.5, 0.5, 0.5"
        jules_pftparm["lai_alb_lim_io"] = "0.005, 0.005, 0.005, 0.005, 0.005"
        jules_pftparm["lma_io"] = "0.0824, 0.2263, 0.0498, 0.1370, 0.0695"
        jules_pftparm["neff_io"] = "0.8e-3, 0.8e-3, 0.8e-3, 0.4e-3, 0.8e-3"
        jules_pftparm["nl0_io"] = "0.040, 0.030, 0.060, 0.030, 0.030"
        jules_pftparm["nmass_io"] = "0.0210, 0.0115, 0.0219, 0.0131, 0.0219"
        jules_pftparm["nr_io"] = "0.01726, 0.00784, 0.0162, 0.0084, 0.01726"
        jules_pftparm["nr_nl_io"] = "1.0, 1.0, 1.0, 1.0, 1.0"
        jules_pftparm["ns_nl_io"] = "0.1, 0.1, 1.0, 1.0, 0.1"
        jules_pftparm["nsw_io"] = "0.0072, 0.0083, 0.01604, 0.0202, 0.0072"
        jules_pftparm["omegal_io"] = "0.10, 0.10, 0.10, 0.12, 0.10"
        jules_pftparm["omegau_io"] = "0.23, 0.23, 0.35, 0.35, 0.35"
        jules_pftparm["omnirl_io"] = "0.50, 0.30, 0.53, 0.53, 0.53"
        jules_pftparm["omniru_io"] = "0.90, 0.65, 0.98, 0.98, 0.98"
        jules_pftparm["orient_io"] = "5*0"
        jules_pftparm["q10_leaf_io"] = "2.0, 2.0, 2.0, 2.0, 2.0"
        jules_pftparm["r_grow_io"] = "0.25, 0.25, 0.25, 0.25, 0.25"
        jules_pftparm["rootd_ft_io"] = "3.0, 1.0, 0.5, 0.5, 0.5"
        jules_pftparm["sigl_io"] = "0.0375, 0.1000, 0.0250, 0.0500, 0.0500"
        jules_pftparm["tleaf_of_io"] = "273.15, 243.15, 258.15, 258.15, 243.15"
        jules_pftparm["tlow_io"] = "0.0, -5.0, 0.0, 13.0, 0.0"
        jules_pftparm["tupp_io"] = "36.0, 31.0, 36.0, 45.0, 36.0"
        jules_pftparm["vint_io"] = "5.73, 6.32, 6.42, 0.00, 14.71"
        jules_pftparm["vsl_io"] = "29.81, 18.15, 40.96, 10.24, 23.15"
        for item, values in jules_pftparm.items():
            self.add_setting(config, ["namelist:jules_pftparm", item], values)

        return config, self.reports
