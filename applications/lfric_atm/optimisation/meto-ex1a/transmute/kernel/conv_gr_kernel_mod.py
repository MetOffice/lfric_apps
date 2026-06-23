# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
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
    first_priv_red_init,
    get_compiler,
    add_omp_parallel_region,
)

logger = logging.getLogger(__name__)


def trans(psyir: Routine):
    """
    Apply PSyClone OpenMP directives

    Transformation takes the following steps:
    - Add OMPParallelDirective around marked loops with OMPParallelTrans
    - Add OMPDoDirective to loops within the parallel regions with OMPLoopTrans
    - Add OMPParallelDoDirective to loops outside the parallel regions with OMPParallelLoopTrans

    Special treatment required for the loop containing glue_conv_6a.
    """

    # Transformations
    omp_trans_parallel_loop_static = OMPParallelLoopTrans(
        omp_schedule="static", omp_directive="paralleldo"
    )
    omp_trans_parallel_loop_dynamic = OMPParallelLoopTrans(
        omp_schedule="dynamic", omp_directive="paralleldo"
    )

    default_trans_options = {"collapse": True}
    if get_compiler() == "cce":
        default_trans_options["collapse"] = False
    loop_trans_options = default_trans_options | {"node-type-check": False}

    for routine in psyir.walk(Routine):
        # Identify large parallel region
        try:
            callnumber_loop = next(filter(is_callnumber_loop, routine.children))
        except StopIteration:
            logger.error("Call-number loop not found in routine.")
            return
        first_loop = next(
            filter(lambda node: isinstance(node, Loop), routine.children)
        )

        case_loop = next(
            filter(
                lambda loop: "ntra_fld" in str(loop.stop_expr),
                routine.walk(Loop),
            ),
            None,
        )

        loops_to_ignore = [case_loop]

        # Ignore nested wtrac loops except inner-most (i) loops to agree with old script
        for loop in filter(is_wtrac_loop, routine.walk(Loop)):
            loops_to_ignore.extend(
                inner_loop
                for inner_loop in loop.walk(Loop)
                if inner_loop.variable.name != "i"
            )

        add_omp_parallel_region(
            first_loop,
            callnumber_loop,
            end_offset=-2,
            ignore_loops=loops_to_ignore,
            loop_trans_options=loop_trans_options,
        )

        # Identify extra parallel regions (loops inside callnumber loop, up to numseg)
        try:
            numseg_loop = next(filter(is_numseg_loop, routine.walk(Loop)))
        except StopIteration:
            logger.error("Numseg loop not found in call-number loop.")
            return

        callnumber_nested_loop = callnumber_loop.loop_body.children[0]

        add_omp_parallel_region(
            callnumber_nested_loop,
            numseg_loop,
            end_offset=-2,
            loop_trans_options=loop_trans_options,
        )

        # Special treatment of numseg loop for parallel do

        # Add "pure" property to specific symbols
        for call in numseg_loop.walk(Call):
            if call.routine.symbol.name in ["glue_conv_6a", "log_event"]:
                call.routine.symbol.is_pure = True
        try:
            omp_trans_parallel_loop_dynamic.apply(
                numseg_loop,
                options=loop_trans_options,
            )
        except TransformationError as e:
            logger.warning(e)

        for loop in routine.walk(Loop):
            # Identify each loop in the routine and add OMPParallelDoDirective to outer loops

            # Ignore specific outer loops (their inner loops will still be transformed)
            if loop.variable.name in ["call_number"]:
                continue

            # Don't attempt to nest parallel directives
            if (
                loop.ancestor(OMPParallelDoDirective) is not None
                or loop.ancestor(OMPDoDirective) is not None
                or loop.ancestor(OMPParallelDirective) is not None
            ):
                continue

            # Main transformations:
            # Add OMPParallelDoDirective outside OMPParallelDirective regions
            try:
                omp_trans_parallel_loop_static.apply(
                    loop,
                    options=default_trans_options
                    | {
                        "ignore_dependencies_for": (
                            ignore_dependencies_block2_after_numseg
                            + ignore_dependencies_block3
                        ),
                    },
                )
            except TransformationError as e:
                logger.warning(e)

        # Fix for firstprivate initialisation bug with CCE
        if get_compiler() == "cce":
            first_private_list = ["k", "conv_active", "orig_value"]
            for node in routine.children:
                if isinstance(node, OMPParallelDirective):
                    first_priv_red_init(node, first_private_list)
                    break


def is_numseg_loop(node: Node):
    """
    Check if node is the num_seg loop: "do i = 1, num_seg"
    """
    return (
        isinstance(node, Loop)
        and node.variable.name == "i"
        and any(ref.name == "num_seg" for ref in node.stop_expr.walk(Reference))
    )


def is_callnumber_loop(node: Node):
    """
    Check if node is the call_number loop: "do call_number = 1, n_conv_calls"
    """
    return (
        isinstance(node, Loop)
        and node.variable.name == "call_number"
        and any(
            ref.name == "n_conv_calls" for ref in node.stop_expr.walk(Reference)
        )
    )


def is_wtrac_loop(node: Node):
    """
    Check if node is a n_wtrac loop: "do _ = 1, n_wtrac, 1"
    """
    return isinstance(node, Loop) and any(
        ref.name == "n_wtrac" for ref in node.stop_expr.walk(Reference)
    )


ignore_dependencies_block2_after_numseg = [
    "conv_rain",
    "conv_snow",
    "cca_2d",
    "cape_diluted",
    "lowest_cca_2d",
    "deep_in_col",
    "shallow_in_col",
    "mid_in_col",
    "freeze_level",
    "deep_prec",
    "shallow_prec",
    "mid_prec",
    "deep_term",
    "cape_timescale",
    "deep_cfl_limited",
    "mid_cfl_limited",
    "deep_tops",
    "dt_conv",
    "dmv_conv",
    "dmcl_conv",
    "dms_conv",
    "massflux_up",
    "massflux_down",
    "conv_rain_3d",
    "conv_snow_3d",
    "entrain_up",
    "entrain_down",
    "detrain_up",
    "detrain_down",
    "dd_dt",
    "dd_dq",
    "deep_massflux",
    "deep_dt",
    "deep_dq",
    "shallow_massflux",
    "shallow_dt",
    "shallow_dq",
    "mid_massflux",
    "mid_dt",
    "mid_dq",
    "cca_unadjusted",
    "massflux_up_half",
    "du_conv",
    "dv_conv",
]

ignore_dependencies_block3 = [
    "o3p",
    "o1d",
    "o3",
    "nit",
    "no",
    "no3",
    "lumped_n",
    "n2o5",
    "ho2no2",
    "hono2",
    "h2o2",
    "ch4",
    "co",
    "hcho",
    "meoo",
    "meooh",
    "h",
    "oh",
    "ho2",
    "cl",
    "cl2o2",
    "clo",
    "oclo",
    "br",
    "lumped_br",
    "brcl",
    "brono2",
    "n2o",
    "lumped_cl",
    "hocl",
    "hbr",
    "hobr",
    "clono2",
    "cfcl3",
    "cf2cl2",
    "mebr",
    "hono",
    "c2h6",
    "etoo",
    "etooh",
    "mecho",
    "meco3",
    "pan",
    "c3h8",
    "n_proo",
    "i_proo",
    "n_prooh",
    "i_prooh",
    "etcho",
    "etco3",
    "me2co",
    "mecoch2oo",
    "mecoch2ooh",
    "ppan",
    "meono2",
    "c5h8",
    "iso2",
    "isooh",
    "ison",
    "macr",
    "macro2",
    "macrooh",
    "mpan",
    "hacet",
    "mgly",
    "nald",
    "hcooh",
    "meco3h",
    "meco2h",
    "h2",
    "meoh",
    "msa",
    "nh3",
    "cs2",
    "csul",
    "h2s",
    "so3",
    "passive_o3",
    "age_of_air",
    "dms",
    "so2",
    "h2so4",
    "dmso",
    "monoterpene",
    "secondary_organic",
    "n_nuc_sol",
    "nuc_sol_su",
    "nuc_sol_om",
    "n_ait_sol",
    "ait_sol_su",
    "ait_sol_bc",
    "ait_sol_om",
    "n_acc_sol",
    "acc_sol_su",
    "acc_sol_bc",
    "acc_sol_om",
    "acc_sol_ss",
    "n_cor_sol",
    "cor_sol_su",
    "cor_sol_bc",
    "cor_sol_om",
    "cor_sol_ss",
    "n_ait_ins",
    "ait_ins_bc",
    "ait_ins_om",
    "n_acc_ins",
    "acc_ins_du",
    "n_cor_ins",
    "cor_ins_du",
    "dmv_conv",
    "conv_prog_precip",
    "conv_prog_precip",
    "dt_conv",
    "conv_prog_dtheta",
    "dmv_conv",
    "conv_prog_dmv",
    "dcfl_conv",
    "dcff_conv",
    "dbcf_conv",
    "dt_conv",
    "dmv_conv",
    "dmcl_conv",
    "dms_conv",
    "dd_mf_cb",
    "cca",
    "ccw",
    "cv_top",
    "cv_base",
    "lowest_cv_top",
    "lowest_cv_base",
    "pres_cv_top",
    "pres_cv_base",
    "pres_lowest_cv_top",
    "pres_lowest_cv_base",
    "massflux_up_cmpta",
    "dth_conv_noshal",
    "dmv_conv_noshal",
    "tke_bl",
]
