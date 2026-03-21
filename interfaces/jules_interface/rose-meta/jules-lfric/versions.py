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
        RMDI = str(-(2 ** 30))
        npft = int(
            self.get_setting_value(
                config, ["namelist:jules_surface_types", "npft"]
            )
        )

        # Add jules_pftparm items hard-wired in jules_physics_init
        jules_pftparm = {}
        jules_pftparm["a_wl_io"] = "0.65,0.65,0.005,0.005,0.10"
        jules_pftparm["a_ws_io"] = "10.0,10.0,1.0,1.0,10.0"
        jules_pftparm["albsnc_min_io"] = "3.0e-1,3.0e-1,8.0e-1,8.0e-1,8.0e-1"
        jules_pftparm["albsnf_maxl_io"] = "0.095,0.059,0.128,0.106,0.077"
        jules_pftparm["albsnf_maxu_io"] = "0.215,0.132,0.288,0.239,0.173"
        jules_pftparm["alnirl_io"] = "0.30,0.23,0.39,0.39,0.39"
        jules_pftparm["alniru_io"] = "0.75,0.65,0.95,0.95,0.87"
        jules_pftparm["alparl_io"] = "0.06,0.04,0.06,0.06,0.06"
        jules_pftparm["alparu_io"] = "0.15,0.11,0.25,0.25,0.25"
        jules_pftparm["alpha_io"] = "0.08,0.08,0.08,0.04,0.08"
        jules_pftparm["b_wl_io"] = "1.667,1.667,1.667,1.667,1.667"
        jules_pftparm["c3_io"] = "3*1,0,1"
        jules_pftparm["can_struct_a_io"] = "1.0,1.0,1.0,1.0,1.0"
        jules_pftparm["dgl_dm_io"] = "0.0,0.0,0.0,0.0,0.0"
        jules_pftparm["dgl_dt_io"] = "9.0,9.0,0.0,0.0,9.0"
        jules_pftparm["dqcrit_io"] = "0.090,0.060,0.100,0.075,0.100"
        jules_pftparm["dust_veg_scj_io"] = "0.0,0.0,1.0,1.0,0.5"
        jules_pftparm["dz0v_dh_io"] = "5.0e-2,5.0e-2,1.0e-1,1.0e-1,1.0e-1"
        jules_pftparm["emis_pft_io"] = "0.98,0.99,0.98,0.98,0.98"
        jules_pftparm["eta_sl_io"] = "0.01,0.01,0.01,0.01,0.01"
        jules_pftparm["f0_io"] = "0.875,0.875,0.900,0.800,0.900"
        jules_pftparm["fd_io"] = "0.015,0.015,0.015,0.025,0.015"
        jules_pftparm["fsmc_of_io"] = "0.0,0.0,0.0,0.0,0.0"
        jules_pftparm["g_leaf_0_io"] = "0.25,0.25,0.25,0.25,0.25"
        jules_pftparm["glmin_io"] = "1.0e-6,1.0e-6,1.0e-6,1.0e-6,1.0e-6"
        jules_pftparm["gsoil_f_io"] = "1.0,1.0,1.0,1.0,1.0"
        jules_pftparm["hw_sw_io"] = "0.5,0.5,0.5,0.5,0.5"
        jules_pftparm["infil_f_io"] = "4.0,4.0,2.0,2.0,2.0"
        jules_pftparm["kn_io"] = "0.78,0.78,0.78,0.78,0.78"
        jules_pftparm["kpar_io"] = "0.5,0.5,0.5,0.5,0.5"
        jules_pftparm["lai_alb_lim_io"] = "0.005,0.005,0.005,0.005,0.005"
        jules_pftparm["lma_io"] = "0.0824,0.2263,0.0498,0.1370,0.0695"
        jules_pftparm["neff_io"] = "0.8e-3,0.8e-3,0.8e-3,0.4e-3,0.8e-3"
        jules_pftparm["nl0_io"] = "0.040,0.030,0.060,0.030,0.030"
        jules_pftparm["nmass_io"] = "0.0210,0.0115,0.0219,0.0131,0.0219"
        jules_pftparm["nr_io"] = "0.01726,0.00784,0.0162,0.0084,0.01726"
        jules_pftparm["nr_nl_io"] = "1.0,1.0,1.0,1.0,1.0"
        jules_pftparm["ns_nl_io"] = "0.1,0.1,1.0,1.0,0.1"
        jules_pftparm["nsw_io"] = "0.0072,0.0083,0.01604,0.0202,0.0072"
        jules_pftparm["omegal_io"] = "0.10,0.10,0.10,0.12,0.10"
        jules_pftparm["omegau_io"] = "0.23,0.23,0.35,0.35,0.35"
        jules_pftparm["omnirl_io"] = "0.50,0.30,0.53,0.53,0.53"
        jules_pftparm["omniru_io"] = "0.90,0.65,0.98,0.98,0.98"
        jules_pftparm["orient_io"] = "5*0"
        jules_pftparm["q10_leaf_io"] = "2.0,2.0,2.0,2.0,2.0"
        jules_pftparm["r_grow_io"] = "0.25,0.25,0.25,0.25,0.25"
        jules_pftparm["rootd_ft_io"] = "3.0,1.0,0.5,0.5,0.5"
        jules_pftparm["sigl_io"] = "0.0375,0.1000,0.0250,0.0500,0.0500"
        jules_pftparm["tleaf_of_io"] = "273.15,243.15,258.15,258.15,243.15"
        jules_pftparm["tlow_io"] = "0.0,-5.0,0.0,13.0,0.0"
        jules_pftparm["tupp_io"] = "36.0,31.0,36.0,45.0,36.0"
        jules_pftparm["vint_io"] = "5.73,6.32,6.42,0.00,14.71"
        jules_pftparm["vsl_io"] = "29.81,18.15,40.96,10.24,23.15"

        # Add new switch and related parameters for photosynthesis model
        # (from JULES vn5.5_t864)
        jules_pftparm["act_jmax_io"] = ",".join(["50.0e3"] * npft)
        jules_pftparm["act_vcmax_io"] = ",".join(["72.0e3"] * npft)
        jules_pftparm["alpha_elec_io"] = ",".join(["0.4"] * npft)
        jules_pftparm["deact_jmax_io"] = ",".join(["200.0e3"] * npft)
        jules_pftparm["deact_vcmax_io"] = ",".join(["200.0e3"] * npft)
        jules_pftparm["ds_jmax_io"] = ",".join(["646.0"] * npft)
        jules_pftparm["ds_vcmax_io"] = ",".join(["649.0"] * npft)
        jules_pftparm["jv25_ratio_io"] = ",".join(["1.97"] * npft)

        # Parameters related to l_bvoc_emis (from UM vn10.1_t605)
        jules_pftparm["ief_io"] = "25.0,8.00,16.00,24.00,20.00"
        jules_pftparm["tef_io"] = "1.2,2.4,0.8,1.2,0.8"
        jules_pftparm["mef_io"] = "0.9,1.8,0.6,0.9,0.57"
        jules_pftparm["aef_io"] = "0.43,0.87,0.29,0.43,0.20"
        jules_pftparm["ci_st_io"] = "33.46,33.46,34.26,29.98,34.26"
        jules_pftparm["gpp_st_io"] = (
            "1.29E-07,2.58E-08,2.07E-07,3.42E-07,1.68E-007"
        )
        # Parameters related to l_inferno (from JULES vn4.4_t136)
        jules_pftparm["fef_co2"] = "1631,1576,1576,1654,1576"
        jules_pftparm["fef_co"] = "100,106,106,64,106"
        jules_pftparm["fef_ch4"] = "6.8,4.8,4.8,2.4,4.8"
        jules_pftparm["fef_nox"] = "2.55,3.24,3.24,2.49,3.24"
        jules_pftparm["fef_so2"] = "0.40,0.40,0.40,0.48,0.40"
        jules_pftparm["fef_oc"] = "4.3,9.1,9.1,3.2,9.1"
        jules_pftparm["fef_bc"] = "0.56,0.56,0.56,0.47,0.56"
        jules_pftparm["ccleaf_min"] = "0.8,0.8,0.8,0.8,0.8"
        jules_pftparm["ccleaf_max"] = "1.0,1.0,1.0,1.0,1.0"
        jules_pftparm["ccwood_min"] = "0.0,0.0,0.0,0.0,0.0"
        jules_pftparm["ccwood_max"] = "0.4,0.4,0.4,0.4,0.4"
        jules_pftparm["avg_ba"] = "0.6E6,0.6E6,1.4E6,1.4E6,1.2E6"
        # Parameters related to l_inferno (from JULES vn7.8_t1579)
        jules_pftparm["fef_c2h4_io"] = (
            "1.11E+00,1.54E+00,8.30E-01,1.99E+00,8.30E-01"
        )
        jules_pftparm["fef_c2h6_io"] = (
            "8.80E-01,9.70E-01,4.20E-01,1.01E+00,4.20E-01"
        )
        jules_pftparm["fef_c3h8_io"] = (
            "5.30E-01,2.90E-01,1.30E-01,3.12E-01,1.30E-01"
        )
        jules_pftparm["fef_hcho_io"] = (
            "2.40E+00,1.75E+00,1.23E+00,2.95E+00,1.23E+00"
        )
        jules_pftparm["fef_mecho_io"] = (
            "2.26E+00,8.10E-01,8.40E-01,2.02E+00,8.40E-01"
        )
        jules_pftparm["fef_nh3_io"] = (
            "1.33E+00,2.50E+00,8.90E-01,2.14E+00,8.90E-01"
        )
        jules_pftparm["fef_dms_io"] = (
            "2.00E-03,2.00E-03,8.00E-03,1.92E-02,8.00E-03"
        )

        # Parameters related to l_o3_damage (value from JULES configurations)
        jules_pftparm["dfp_dcuo_io"] = "0.04,0.02,0.25,0.13,0.03"
        jules_pftparm["fl_o3_ct_io"] = "1.6,1.6,5.0,5.0,1.6"

        # Parameters related to l_trif_fire (from JULES vn5.3_t872)
        jules_pftparm["fire_mort_io"] = ",".join(["1.0"] * npft)

        # Parameters related to Medlyn stomata model (from JULES vn5.3_t766)
        jules_pftparm["g1_stomata_io"] = ",".join(["2.0"] * npft)

        # Parameters related to SOX stomata model (from JULES vn7.4_t1491)
        jules_pftparm["sox_a_io"] = ",".join(["0.0"] * npft)
        jules_pftparm["sox_p50_io"] = ",".join(["0.0"] * npft)
        jules_pftparm["sox_rp_min_io"] = ",".join(["0.0"] * npft)

        # Parameters related to l_use_pft_psi (from JULES vn4.8_t541)
        jules_pftparm["psi_close_io"] = ",".join(["-1.5E6"]*npft)
        jules_pftparm["psi_open_io"] = ",".join(["-0.033E6"]*npft)

        # Parameters related to l_sugar (from JULES vn7.3_t1344)
        jules_pftparm["sug_grec_io"] = ",".join(["1.0"]*npft)
        jules_pftparm["sug_g0_io"] = ",".join(["1.0"]*npft)
        jules_pftparm["sug_yg_io"] = ",".join(["1.0"]*npft)

        # Remaining parameters added with missing data as no other information
        jules_pftparm["albsnf_max_io"] = ",".join([RMDI]*npft)
        jules_pftparm["fsmc_mod_io"]   = ",".join([RMDI]*npft)

        for item, values in jules_pftparm.items():
            self.add_setting(config, ["namelist:jules_pftparm", item], values)

        # Add jules_radiation switches related to added pftparms
        self.add_setting(
            config, ["namelist:jules_radiation", "l_spec_albedo"], ".true."
        )

        # Add jules_vegetation switches related to added pftparms
        self.add_setting(
            config, ["namelist:jules_vegetation", "photo_model"], "'collatz'"
        )
        self.add_setting(
            config, ["namelist:jules_vegetation", "stomata_model"], "'jacobs'"
        )
        self.add_setting(
            config, ["namelist:jules_vegetation", "l_bvoc_emis"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_vegetation", "l_inferno"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_vegetation", "l_o3_damage"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_vegetation", "l_sugar"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_vegetation", "l_trif_fire"], ".false."
        )
        self.add_setting(
            config, ["namelist:jules_vegetation", "l_use_pft_psi"], ".false."
        )


return config, self.reports
