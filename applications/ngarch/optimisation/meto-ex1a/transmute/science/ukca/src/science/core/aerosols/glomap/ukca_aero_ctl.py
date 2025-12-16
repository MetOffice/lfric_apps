# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
OMP CPU transform script for UKCA GLOMAP.

Provides:
    `trans(psyir)`: performs the OpenMP transformation of the main
      `ukca_aero_ctl` loop.

Rationale:
    GLOMAP processes are gridbox-independent. In the top-level GLOMAP
    subroutine (`ukca_aero_ctl`), the domain is chunked into segments. In order
    to maximise OpenMP coverage, the most efficient way to parallelise is by
    doing so over these chunks. Performance will vary depending on the size of
    the domain, and the user-specified chunk size (determined by
    `ukca_mode_seg_size` in the aerosol namelist).

Future work:
    Any changes to `ukca_aero_ctl` will possibly necessitate changes to the
    transformer.

    PSyclone as of version 3.1 is not capable of fully
    determining which variables need to be set PRIVATE, since many of these
    are arrays (and using `privatise_arrays` has not been found to work). The
    suggested fix from STFC is to use `explicitly_private_symbols`, which is
    captured in `mark_explicit_privates` from `transmute_functions`.
    The workaround is to set anything beginning in `seg_` to PRIVATE, since
    these variables are the ones used in the chunking and therefore need to
    be private to each thread.

    An issue has been identified with PSyclone version 3.2:
    https://github.com/stfc/PSyclone/issues/3173
    which causes the transformer to fail due to a bug where CodeBlocks did
    not have an `enters_scope` attribute. A fix is in progress at
    https://github.com/stfc/PSyclone/tree/3173_quickfix.
"""
from psyclone.psyir.nodes import Loop
from psyclone.transformations import OMPParallelLoopTrans, TransformationError
from psyclone.psyir.symbols import DataSymbol
from transmute_psytrans.transmute_functions import (
    first_priv_red_init,
    mark_explicit_privates,
)

OMP_TRANS = OMPParallelLoopTrans()
RESOLVE_IMPORTS = True


def get_private_symbols_from_name(node, search_str):
    """
    Find every symbol in the symbol table containing a user-defined string.
    A better alternative might be to add wildcard support to
    symbol_table.lookup(). This is necessary since PSyclone (as of version 3.1)
    does not make arrays PRIVATE, which is required for parallelising the
    loop over segments in GLOMAP.
    """
    symbols_to_add = []
    # Search through the symbol table of the node; if any entries match the
    # search string, add them to the list of symbols to be made PRIVATE.
    for symbol in node.parent.symbol_table.symbols:
        if search_str in symbol.name and isinstance(symbol, DataSymbol):
            symbols_to_add.append(symbol.name)
    return symbols_to_add


def trans(psyir):
    """
    Apply the PSyclone transformation. If it fails, raise an error.
    """
    opts = {
        # some non-PURE subroutines called within this loop
        "force": True,
        # several WRITE statements used for diagnostics
        "node-type-check": False,
    }
    # For the coarse-grained approach, we have *one* loop we want to work on
    # - the loop over segments. This gives almost complete coverage for GLOMAP,
    # with no race conditions (as each gridbox is independent).
    loops = psyir.walk(Loop)
    for loop in loops:
        # identify the loop in question - the loop over segments
        if hasattr(loop.stop_expr, "name") and loop.stop_expr.name in ["nseg"]:
            try:
                symbols_to_add = get_private_symbols_from_name(loop, "seg_")
                # add a few more symbols that don't fit this syntax
                symbols_to_add.extend(
                    ["i_end", "i_end_cp", "j", "nbs_index", "y"]
                )
                # set some symbols to be PRIVATE
                mark_explicit_privates(loop, symbols_to_add)
                # CCE compiler fix - initialise some FIRSTPRIVATE variables
                first_priv_red_init(
                    loop,
                    ["cmessage", "errcode", "i_end", "i_end_cp",
                     "i_start_cp", "iaer", "k", "nbs_index", "seg_aird",
                     "seg_airdm3", "seg_autoconv1d", "seg_bud_aer_mas",
                     "seg_ccn_1", "seg_ccn_2", "seg_ccn_3", "seg_ccn_4",
                     "seg_ccn_5", "seg_cdn", "seg_clf", "seg_clwc",
                     "seg_cn_3nm", "seg_craing", "seg_craing_up",
                     "seg_csnowg", "seg_delso2", "seg_delso2_2",
                     "seg_draing", "seg_drydp", "seg_dsnowg", "seg_dvisc",
                     "seg_dvol", "seg_erf_arg", "seg_erfterm", "seg_fac",
                     "seg_fconv_conv", "seg_height", "seg_het_rates",
                     "seg_htpblg", "seg_iarr", "seg_ilscat", "seg_jlabove",
                     "seg_karr", "seg_land_frac", "seg_larr", "seg_lday",
                     "seg_lowcloud", "seg_lwc", "seg_md", "seg_mdt",
                     "seg_mdtfixflag", "seg_mdtfixsink", "seg_mdwat",
                     "seg_mfpa", "seg_n_merge_1d", "seg_nbadmdt", "seg_nd",
                     "seg_plower", "seg_pmid", "seg_pupper", "seg_pvol",
                     "seg_pvol_wat", "seg_rh", "seg_rh_clr", "seg_rhoa",
                     "seg_rhopar", "seg_s", "seg_s0", "seg_s0_dot_condensable",
                     "seg_sarea", "seg_seaice", "seg_sm", "seg_surtp",
                     "seg_t", "seg_tr_rs", "seg_tsqrt", "seg_ustr",
                     "seg_v1d_tmp", "seg_vba", "seg_vconc", "seg_vfac",
                     "seg_wetdp", "seg_wvol", "seg_zh2o2", "seg_zho2",
                     "seg_znotg", "seg_zo3", "y"],
                    insert_at_start=True,
                )
                OMPParallelLoopTrans(omp_schedule="dynamic").apply(
                    loop, options=opts
                )

            except TransformationError as err:
                raise TransformationError(
                    f"ukca_aero_ctl.py: Error: "
                    f"could not apply OMP transformation "
                    f"to loop: {str(err)}"
                ) from err
