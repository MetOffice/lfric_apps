# -----------------------------------------------------------------------------
# (C) 2026 Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
PSyclone script for applying OpenMP transformations specific to the
Gregory-Rowntree convection kernel.
"""
import logging
from psyclone.psyir.nodes import (
    Call,
    Loop,
    Node,
    OMPParallelDoDirective,
    Reference,
    Routine,
    OMPParallelDirective,
    OMPDoDirective,
)
from psyclone.transformations import (
    OMPParallelLoopTrans,
    TransformationError,
)
from transmute_psytrans.transmute_functions import (
    match_lhs_assignments,
    OMP_PARALLEL_LOOP_DO_TRANS_STATIC,
    OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC,
)


logger = logging.getLogger(__name__)

false_dep_vars_all = [
    "conv_rain", "conv_snow", "cca_2d", "cape_diluted", "lowest_cca_2d",
    "deep_in_col", "shallow_in_col", "mid_in_col", "freeze_level",
    "deep_prec", "shallow_prec", "mid_prec", "deep_term", "cape_timescale",
    "deep_cfl_limited", "mid_cfl_limited", "deep_tops", "dt_conv",
    "dmv_conv", "dmcl_conv", "dms_conv", "massflux_up", "massflux_down",
    "conv_rain_3d", "conv_snow_3d", "entrain_up", "entrain_down",
    "detrain_up", "detrain_down", "dd_dt", "dd_dq", "deep_massflux",
    "deep_dt", "deep_dq", "shallow_massflux", "shallow_dt", "shallow_dq",
    "mid_massflux", "mid_dt", "mid_dq", "cca_unadjusted",
    "massflux_up_half", "du_conv", "dv_conv", "dmv_conv",
    "conv_prog_precip", "conv_prog_precip", "dt_conv", "conv_prog_dtheta",
    "conv_prog_dmv", "dcfl_conv", "dcff_conv", "dbcf_conv", "dd_mf_cb",
    "cca", "ccw", "cv_top", "cv_base", "lowest_cv_top", "lowest_cv_base",
    "pres_cv_top", "pres_cv_base", "pres_lowest_cv_top",
    "pres_lowest_cv_base", "massflux_up_cmpta", "dth_conv_noshal",
    "dmv_conv_noshal", "tke_bl", "o3p", "o1d", "o3", "nit", "no", "no3",
    "lumped_n", "n2o5", "ho2no2", "hono2", "h2o2", "ch4", "co", "hcho",
    "meoo", "meooh", "h", "oh", "ho2", "cl", "cl2o2", "clo", "oclo",
    "br", "lumped_br", "brcl", "brono2", "n2o", "lumped_cl", "hocl",
    "hbr", "hobr", "clono2", "cfcl3", "cf2cl2", "mebr", "hono", "c2h6",
    "etoo", "etooh", "mecho", "meco3", "pan", "c3h8", "n_proo", "i_proo",
    "n_prooh", "i_prooh", "etcho", "etco3", "me2co", "mecoch2oo",
    "mecoch2ooh", "ppan", "meono2", "c5h8", "iso2", "isooh", "ison",
    "macr", "macro2", "macrooh", "mpan", "hacet", "mgly", "nald",
    "hcooh", "meco3h", "meco2h", "h2", "meoh", "msa", "nh3", "cs2",
    "csul", "h2s", "so3", "passive_o3", "age_of_air", "dms", "so2",
    "h2so4", "dmso", "monoterpene", "secondary_organic", "n_nuc_sol",
    "nuc_sol_su", "nuc_sol_om", "n_ait_sol", "ait_sol_su", "ait_sol_bc",
    "ait_sol_om", "n_acc_sol", "acc_sol_su", "acc_sol_bc", "acc_sol_om",
    "acc_sol_ss", "n_cor_sol", "cor_sol_su", "cor_sol_bc", "cor_sol_om",
    "cor_sol_ss", "n_ait_ins", "ait_ins_bc", "ait_ins_om", "n_acc_ins",
    "acc_ins_du", "n_cor_ins", "cor_ins_du"]

false_dep_vars_seg = [
    "ncells", "cumulus", "l_shallow", "l_mid", "bulk_cf_conv", "cape_ts_used",
    "ccp_strength", "cf_frozen_conv", "cf_liquid_conv", "conv_prog_flx",
    "conv_prog_precip_conv", "conv_type", "dbcfbydt", "dcffbydt", "dcflbydt",
    "deep_flag", "delta_smag", "delthvu", "dqbydt", "dqbydt_wtrac", "dqcfbydt",
    "dqcfbydt_wtrac", "dqclbydt", "dqclbydt_wtrac", "dthbydt", "dubydt_p",
    "dvbydt_p", "entrain_coef", "exner_rho_levels", "exner_rho_minus_one",
    "exner_theta_levels", "fqw", "freeze_lev", "ftl", "g_ccp", "h_ccp",
    "ind_cape_reduced", "it_area_dd", "it_area_ud", "it_cape_diluted",
    "it_cca", "it_cca0", "it_cca0_dp", "it_cca0_md", "it_cca0_sh", "it_cca_2d",
    "it_ccb", "it_ccb0", "it_cclwp", "it_cclwp0", "it_cct", "it_cct0",
    "it_ccw", "it_ccw0", "it_cg_term", "it_conv_rain", "it_conv_rain_3d",
    "it_conv_snow", "it_conv_snow_3d", "it_detrain_dwn", "it_detrain_up",
    "it_dp_cfl_limited", "it_dq_congest", "it_dq_dd", "it_dq_deep",
    "it_dq_midlev", "it_dq_shall", "it_dt_congest", "it_dt_dd", "it_dt_deep",
    "it_dt_midlev", "it_dt_shall", "it_du_congest", "it_du_dd", "it_du_deep",
    "it_du_midlev", "it_du_shall", "it_dv_congest", "it_dv_dd", "it_dv_deep",
    "it_dv_midlev", "it_dv_shall", "it_dwn_flux", "it_entrain_dwn",
    "it_entrain_up", "it_ind_deep", "it_ind_shall", "it_kterm_deep",
    "it_kterm_shall", "it_lcbase", "it_lcbase0", "it_lcca", "it_lctop",
    "it_mb1", "it_mb2", "it_md_cfl_limited", "it_mf_congest", "it_mf_deep",
    "it_mf_midlev", "it_mf_shall", "it_mid_level", "it_precip_cg",
    "it_precip_dp", "it_precip_md", "it_precip_sh", "it_up_flux",
    "it_up_flux_half", "it_uw_dp", "it_uw_mid", "it_uw_shall", "it_vw_dp",
    "it_vw_mid", "it_vw_mid", "it_vw_shall", "it_w2p", "it_wql_flux",
    "it_wqt_flux", "it_wstar_dn", "it_wstar_up", "it_wthetal_flux",
    "it_wthetav_flux", "l_congestus", "land_sea_mask", "ntml", "ntpar",
    "p_rho_minus_one", "p_star", "p_theta_levels", "past_conv_ht", "q1_sd",
    "q_conv", "q_wtrac", "qcf_conv", "qcf_wtrac", "qcl_conv", "qcl_wtrac",
    "ql_ad", "qsat_lcl", "r_rho_levels", "r_theta_levels", "rain_wtrac",
    "rho_dry", "rho_dry_theta", "rho_wet", "rho_wet_tq", "scm_convss_dg",
    "snow_wtrac",
    "t1_sd", "theta_conv", "tnuc_new", "tnuc_nlcl_um", "tot_tracer", "u_conv",
    "uw0", "v_conv", "vw0", "w", "w_max", "wstar", "wthvs", "z_rho", "z_theta",
    "zhpar", "zlcl", "zlcl_uv",
]

def trans(psyir: Routine):
    """
    Apply PSyClone OpenMP directives

    Transformation takes the following steps:
    - Add OMPParallelDirective around marked loops with OMPParallelTrans
    - Add OMPDoDirective to loops within the parallel regions with OMPLoopTrans
    - Add OMPParallelDoDirective to loops outside the parallel regions with OMPParallelLoopTrans

    Special treatment required for the loop containing glue_conv_6a.

    :param psyir: the PSyIR of the provided file.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`
    """


    # Identify extra parallel regions (loops inside callnumber loop, up to numseg)
    numseg_loop = None
    try:
        numseg_loop = next(filter(is_numseg_loop, psyir.walk(Loop)))
    except StopIteration:
        logger.error("Numseg loop not found in call-number loop.")

    # Special treatment of numseg loop for parallel do
    if numseg_loop:
        # Add "pure" property to specific symbols
        for call in numseg_loop.walk(Call):
            if call.routine.symbol.name in ["glue_conv_6a", "log_event"]:
                call.routine.symbol.is_pure = True
        try:
            OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC.apply(
                numseg_loop,
                ignore_dependencies_for=false_dep_vars_seg,
                node_type_check=False)
        except TransformationError as e:
            logger.warning(e)
            print(f"Trying loop but{e}")

    # Work through each other loop in the file and OMP PARALLEL DO
    for loop in psyir.walk(Loop):
        # Don't attempt to nest parallel directives
        if (
            loop.ancestor(OMPParallelDoDirective) is not None
            or loop.ancestor(OMPDoDirective) is not None
            or loop.ancestor(OMPParallelDirective) is not None
        ):
            continue
        if loop.variable.name in ['i']:
            try:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(
                    loop, 
                    ignore_dependencies_for=false_dep_vars_all,
                    node_type_check=False)

            except (TransformationError, IndexError) as err:
                logging.warning(f"Could not transform because:\n {err}")


def is_numseg_loop(node: Node):
    """
    Check if node is the num_seg loop: "do i = 1, num_seg"
    """
    return (
        isinstance(node, Loop)
        and node.variable.name == "i"
        and any(ref.name == "num_seg" for ref in node.stop_expr.walk(Reference))
    )

