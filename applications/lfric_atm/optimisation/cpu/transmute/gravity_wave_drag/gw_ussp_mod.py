# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
PSyclone Script for applying OpenMP transformations to Spectral Gravity Wave.

Provides:
  - `trans(psyir)`: orchestrates all OpenMP transformations.
    * Forces a dynamic-scheduled OMP **PARALLEL DO** over
      `do i = 1, meta_segments%num_segments` with an explicit PRIVATE list.
    * Clusters adjacent top-level loops into **PARALLEL regions** to amortize
      thread overheads.
    * Applies **STATIC-scheduled** OMP DO / PARALLEL DO to identified heavy
      k-loops (profiling hotspots).

Context (specific to gw_ussp_mod.F90):
  - Default scheduling is **static** everywhere **except** for the
    `meta_segments%num_segments` loop, which is explicitly **dynamic** to
    accommodate variable per-segment workloads.
  - Certain scalars are initialised prior to OpenMP transformation so that
    backends that emit FIRSTPRIVATE have a well-defined initial value.

Rationale:
  - Profiling indicates k-index loops writing to heavy variables dominate
    runtime. Static scheduling generally performs best for these loops, while
    the meta-segmentation loop benefits from dynamic scheduling due to load
    imbalance.
  - Grouping adjacent outer loops into a single parallel region reduces
    enter/exit overhead and improves locality.

Future Work:
  - As PSyclone evolves, revisit nowait/lastprivate usage and refine region
    boundaries.
"""
import logging
from psyclone.psyir.nodes import Routine
from transmute_psytrans.transmute_functions import (
    add_parallel_do_over_meta_segments,
    parallel_regions_for_clustered_loops,
    omp_do_for_heavy_loops,
    get_compiler,
)

# ------------------------------------------------------------------------------
# OpenMP transformation objects
#
# Policy:
#   - STATIC schedule by default for heavy k-loops (best throughput observed).
#   - DYNAMIC schedule **only** for the PARALLEL DO
#     over meta_segments%num_segments.
# ------------------------------------------------------------------------------
# NOTE: OMP_* transform instances now live in the generic module; kept comment
# here to preserve file style and intent.

# Variables considered "heavy" when written inside loops:
#   - HEAVY_VARS_K → for k-loops
#   - HEAVY_VARS_I → for i-loops
# These guide where OMP DO / PARALLEL DO (static) should be applied.
HEAVY_VARS_K = {
    "rhont_smallhalo", "rho_th", "udotk", "g_x", "g_y", "t_inc",
    "r_u", "r_v", "g_xp_smallhalo", "g_yp_smallhalo",
}
HEAVY_VARS_I = {"rho_th", "nbv", "udotk"}

# --- meta_segments PARALLEL DO (forced, dynamic), with explicit privates ----
# Explicit PRIVATE list for the meta-segmentation PARALLEL DO
_META_PRIVATES = [
    # scalars
    "ii", "jj", "ss", "i", "k", "jdir",
    # allocatables used per-segment
    "s_sin_theta_lat", "s_totalppn", "s_rho_th", "s_nbv", "s_udotk", "s_fptot",
]


# ------------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------------
def trans(psyir):
    """
    Entry point for OpenMP transformations for `gw_ussp_mod.F90`.

    Passes:
      1) Force PARALLEL DO (dynamic) for meta-segmentation loop.
      2) Cluster adjacent top-level loops into PARALLEL regions.
      3) Apply OMP DO / PARALLEL DO (static) to heavy k-loops.
    """
    logging.info(
        "gw_ussp_mod.F90 Psyclone Optimisation for the CPU is running."
    )

    compiler = get_compiler()

    # GCC currently has issues with firstprivated indexes,
    # which is soon to be resolved in PSyclone.
    # However for now we will avoid using OpenMP
    # around meta_segments with with GCC.
    if compiler == "gnu":
        logging.info(
            "Skipping gw_ussp_mod meta_segment optimisations for GCC."
            )
    else:
        # 1) Force a PARALLEL DO around meta_segments%num_segments (dynamic)
        for routine in psyir.walk(Routine):
            add_parallel_do_over_meta_segments(
              routine,
              container_name="meta_segments",
              member_name="num_segments",
              privates=_META_PRIVATES,
            )

    # 2) Cluster adjacent outer loops into PARALLEL regions (no schedule here)
    for routine in psyir.walk(Routine):
        parallel_regions_for_clustered_loops(routine)

    # 3) Add OMP DO / PARALLEL DO on heavy k- and i-loops (static)
    #    (Skip the special i-loop handled above.)
    for routine in psyir.walk(Routine):
        omp_do_for_heavy_loops(routine, "k", HEAVY_VARS_K)
        omp_do_for_heavy_loops(
          routine, "i", HEAVY_VARS_I,
          skip_member_count=("i", "meta_segments", "num_segments")
        )
