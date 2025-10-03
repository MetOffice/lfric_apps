# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
import logging
import os
from typing import Sequence, Optional, Tuple, Set

from psyclone.psyir.nodes import (
    Loop, Assignment, Reference,
    OMPParallelDirective,
    OMPDoDirective, OMPParallelDoDirective,
    StructureReference, Member, Literal
)
from psyclone.psyir.symbols import DataSymbol
from psyclone.transformations import (
    OMPLoopTrans, TransformationError, OMPParallelTrans, OMPParallelLoopTrans
)

# ------------------------------------------------------------------------------
# OpenMP transformation objects
#
# Policy:
#   - STATIC schedule by default for heavy loops (best throughput observed).
#   - DYNAMIC schedule **only** for the special PARALLEL DO case that the app
#     identifies (e.g. a member-count loop like meta_segments%num_segments).
# ------------------------------------------------------------------------------
OMP_PARALLEL_REGION_TRANS = OMPParallelTrans()

# Default: static schedule for heavy loops
OMP_DO_LOOP_TRANS_STATIC = OMPLoopTrans(omp_schedule="static")
OMP_PARALLEL_LOOP_DO_TRANS_STATIC = OMPParallelLoopTrans(
    omp_schedule="static", omp_directive="paralleldo"
)

# Exception: dynamic schedule for the app-selected special loop
OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC = OMPParallelLoopTrans(
    omp_schedule="dynamic", omp_directive="paralleldo"
)
OMP_DO_LOOP_TRANS_DYNAMIC = OMPLoopTrans(omp_schedule="dynamic")
# ------------------------------------------------------------------------------


def is_heavy_loop(loop, heavy_vars: Set[str]):
    """
    Determine whether `loop` performs significant work on key variables.

    Parameters
    ----------
    loop : psyclone.psyir.nodes.Loop
        Candidate loop to inspect.
    heavy_vars : set[str]
        Names of variables (profiling hotspots) that mark a loop as heavy.

    Returns
    -------
    bool
        True if any Assignment within the loop writes to a Reference whose
        name is in `heavy_vars`; False otherwise.
    """
    for assign in loop.walk(Assignment):
        lhs = assign.lhs
        if isinstance(lhs, Reference) and lhs.name in heavy_vars:
            return True
    return False


def get_outer_loops(node):
    """
    Retrieve non-nested top-level Loop nodes under a PSyIR node.

    Returns only loops without Loop ancestors already collected, enabling
    parallel-region clustering at a consistent nesting level.
    """
    outer_loops = []
    for loop in node.walk(Loop):
        if loop.ancestor(Loop) not in outer_loops:
            outer_loops.append(loop)
    return outer_loops


def parallel_regions_for_clustered_loops(routine):
    """
    Enclose clusters of adjacent top-level loops in a single PARALLEL region.

    Notes
    -----
    - No schedule is specified at region level
      (loop-level directives handle it).
    """
    logging.info("Processing Routine for regions: '%s'", routine.name)

    outer_loops = get_outer_loops(routine)
    if not outer_loops:
        logging.info("No loops to regionize.")
        return

    # Build sortable (parent, child-index, loop) tuples and sort once.
    items = [
        (lp.parent, lp.parent.children.index(lp), lp) for lp in outer_loops
    ]
    items.sort(key=lambda t: (id(t[0]), t[1]))

    current_parent = None
    cluster = []
    prev_idx = None

    for parent, idx, loop in items:
        if parent is not current_parent:
            # Flush previous parent's tail cluster (if any)
            if (
                len(cluster) > 1
                and not any(
                    lp.ancestor(OMPParallelDirective) for lp in cluster
                )
            ):
                positions = f"{cluster[0].position}-{cluster[-1].position}"
                logging.info(
                    "Inserting region over loops at positions %s", positions
                )
                try:
                    OMP_PARALLEL_REGION_TRANS.apply(cluster)
                    logging.info("Region inserted.")
                except TransformationError as err:
                    logging.info("Region failed: %s", err)
            current_parent = parent
            cluster = [loop]
            prev_idx = idx
            continue

        if idx == prev_idx + 1:
            cluster.append(loop)
            prev_idx = idx
            continue

        # Non-adjacent: flush current cluster, start a new one.
        if (
            len(cluster) > 1
            and not any(lp.ancestor(OMPParallelDirective) for lp in cluster)
        ):
            positions = f"{cluster[0].position}-{cluster[-1].position}"
            logging.info(
                "Inserting region over loops at positions %s", positions
            )
            try:
                OMP_PARALLEL_REGION_TRANS.apply(cluster)
                logging.info("Region inserted.")
            except TransformationError as err:
                logging.info("Region failed: %s", err)

        cluster = [loop]
        prev_idx = idx

    # Final tail cluster.
    if (
        len(cluster) > 1
        and not any(lp.ancestor(OMPParallelDirective) for lp in cluster)
    ):
        positions = f"{cluster[0].position}-{cluster[-1].position}"
        logging.info("Inserting region over loops at positions %s", positions)
        try:
            OMP_PARALLEL_REGION_TRANS.apply(cluster)
            logging.info("Region inserted.")
        except TransformationError as err:
            logging.info("Region failed: %s", err)


def expr_contains_member(expr, container_name: str, member_name: str) -> bool:
    """
    Return True iff `expr` contains `<container_name>%<member_name>` anywhere
    within a StructureReference tree.

    Parameters
    ----------
    expr : PSyIR node (typically an expression)
    container_name : str
        The symbol name of the container (e.g. "meta_segments").
    member_name : str
        The member name (e.g. "num_segments").
    """
    for sref in expr.walk(StructureReference):
        # Defensive: not all StructureReference nodes guarantee a symbol
        try:
            if sref.symbol.name != container_name:
                continue
        except AttributeError:
            continue
        for mem in sref.walk(Member):
            if mem.name == member_name:
                return True
    return False


def omp_do_for_heavy_loops(
    routine,
    loop_var: str,
    heavy_vars: Set[str],
    skip_member_count: Optional[Tuple[str, str, str]] = None,
):
    """
    Insert OMP DO / PARALLEL DO (STATIC) for heavy loops over `loop_var`.

    Behavior
    --------
    - For each Loop where loop.variable.name == `loop_var` AND the loop writes
      to any name in `heavy_vars`:
        * If inside an OMP PARALLEL region -> apply OMP DO (static).
        * Otherwise                        -> apply PARALLEL DO (static).

    Special case
    ------------
    - If `skip_member_count=(loop_var, container, member)` is provided and
      matches the loop, SKIP that loop here so it can be handled by
      add_parallel_do_over_meta_segments() with a DYNAMIC schedule.
    """
    logging.info(
        "Processing Routine for heavy '%s'-loops: '%s'",
        loop_var,
        routine.name,
    )

    for loop in routine.walk(Loop):
        if not (loop.variable and loop.variable.name == loop_var):
            continue

        # Optional: skip an app-selected member-count loop (handled elsewhere)
        if skip_member_count is not None:
            lv, cont, mem = skip_member_count
            if loop_var == lv:
                stop_expr = getattr(loop, "stop_expr", None)
                if stop_expr and expr_contains_member(stop_expr, cont, mem):
                    continue

        # Only consider loops that write to heavy variables
        if not is_heavy_loop(loop, heavy_vars):
            continue

        # Avoid double annotation if already under an OMP DO/PARALLEL DO
        already_omp_do = bool(
            loop.ancestor((OMPDoDirective, OMPParallelDoDirective))
        )
        if already_omp_do:
            logging.info(
                "%s-loop at %s already inside OMP DO; skipping.",
                loop_var,
                loop.position,
            )
            continue

        in_parallel_region = bool(loop.ancestor(OMPParallelDirective))
        logging.info(
            "  %s-loop at %s: schedule=static (%s)",
            loop_var,
            loop.position,
            "in-region" if in_parallel_region else "standalone",
        )
        try:
            if in_parallel_region:
                OMP_DO_LOOP_TRANS_STATIC.apply(loop)
            else:
                OMP_PARALLEL_LOOP_DO_TRANS_STATIC.apply(loop)
            logging.warning("OMP applied to %s-loop (static).", loop_var)
        except TransformationError as err:
            logging.warning("Failed OMP on %s-loop: %s", loop_var, err)


def mark_explicit_privates(node, names):
    """
    Add symbols named in `names` to `node.explicitly_private_symbols`.

    Generic version of the original helper. Works with any PSyIR node that:
      - has a `scope.symbol_table`, and
      - provides an `explicitly_private_symbols` set-like attribute.

    Warns if a requested symbol cannot be found or is not a DataSymbol.
    """
    # Be forgiving so this helper can be used beyond Loop nodes
    scope = getattr(node, "scope", None)
    symtab = getattr(scope, "symbol_table", None)
    if symtab is None:
        logging.warning(
            "[warn] cannot set explicit privates:"
            "node has no scope.symbol_table."
        )
        return

    if not hasattr(
        node, "explicitly_private_symbols"
    ):
        logging.warning(
            "[warn] cannot set explicit privates:"
            " node has no 'explicitly_private_symbols'."
        )
        return

    for name in names:
        try:
            sym = symtab.lookup(name)
            if isinstance(sym, DataSymbol):
                node.explicitly_private_symbols.add(sym)
            else:
                logging.warning(
                    " [warn] private symbol '%s' is not a DataSymbol.",
                    name,
                )
        except KeyError:
            logging.warning(
                "[warn] private symbol '%s' not found in symbol table.",
                name,
            )


def get_compiler():
    """
    Best-effort compiler family from env.
    Prefers FC, then CC, then module hints.
    Returns: 'gnu' | 'intel' | 'cce' | 'nvhpc' | None
    """
    keys = ("COMPILER", "FC", "CC", "LOADEDMODULES", "_LMFILES_")
    for key in keys:
        val = os.environ.get(key)
        if not val:
            continue
        s = val.strip().lower()

        # GNU / GCC
        if ("gfortran" in s or "gcc" in s or "gnu" in s):
            return "gnu"

        # Intel (classic/oneAPI)
        if ("ifx" in s or "ifort" in s or "intel" in s
                or "icx" in s or "icc" in s):
            return "intel"

        # Cray CCE
        if ("cce" in s or "crayftn" in s or "cray" in s):
            return "cce"

        # NVIDIA HPC / PGI
        if ("nvfortran" in s or "nvc" in s or "nvhpc" in s
                or "pgfortran" in s or "pgi" in s):
            return "nvhpc"
    return None


def add_parallel_do_over_meta_segments(
    routine,
    container_name: str,
    member_name: str,
    privates: Sequence[str],
    init_scalars: Sequence[str] = ("jdir", "k"),
):
    """
    Force an OMP PARALLEL DO with **dynamic** schedule over meta-segments.

    Search
    ------
    Find `do i = 1, <container_name>%<member_name>` and:
      1) Insert initialisations for scalars that may be emitted as FIRSTPRIVATE
         (e.g., jdir, k) so they have a defined value.
      2) Mark explicit PRIVATE variables per `privates`.
      3) Apply:
         - **OMP DO (dynamic)** if the loop is
           **already inside** an OMP region.
         - **PARALLEL DO (dynamic)** otherwise.

    Notes
    -----
    - This transformation is forced (`options={"force": True}`) to ensure that
      scheduling is **dynamic** regardless of the default static policy.
    """
    logging.info(
        "Processing Routine for meta_segments loop: '%s'",
        routine.name,
    )

    # Locate the target loop: do i = 1, <container_name>%<member_name>
    target = None
    for loop in routine.walk(Loop):
        if not loop.variable or loop.variable.name != "i":
            continue
        stop_expr = getattr(loop, "stop_expr", None)
        if stop_expr and expr_contains_member(
            stop_expr, container_name, member_name
        ):
            target = loop
            break

    if not target:
        logging.info("  meta-segments style member-count loop not found.")
        return

    # Determine OpenMP context precisely:
    # - If already under an OMP DO/PARALLEL DO,
    #   skip to avoid double annotation.
    # - If inside an OMP PARALLEL region,
    #   we should apply only OMP DO (dynamic).
    # - Otherwise, apply OMP PARALLEL DO (dynamic).
    already_omp_do = bool(
        target.ancestor((OMPDoDirective, OMPParallelDoDirective))
    )
    if already_omp_do:
        logging.info(
            "Target loop already has an OMP DO/Parallel DO ancestor; skipping."
        )
        return
    in_parallel_region = bool(target.ancestor(OMPParallelDirective))

    # Ensure scalars that may be emitted as FIRSTPRIVATE have a value
    parent = target.parent
    insert_at = parent.children.index(target)
    for nm in init_scalars:  # e.g., ("jdir", "k")
        try:
            sym = target.scope.symbol_table.lookup(nm)
        except KeyError:
            continue
        init = Assignment.create(Reference(sym), Literal("0", sym.datatype))
        parent.children.insert(insert_at, init)
        insert_at += 1

    # Explicit privates per policy
    mark_explicit_privates(target, privates)

    # Apply the dynamic-scheduled directive (forced)
    try:
        if in_parallel_region:
            logging.info(
                "Found target loop at %s inside OMP parallel region: "
                "applying OMP DO (forced, dynamic).",
                target.position,
            )
            OMP_DO_LOOP_TRANS_DYNAMIC.apply(
                target, options={"force": True}
            )
        else:
            logging.info(
                "Found target loop at %s:"
                " applying OMP PARALLEL DO (forced, dynamic).",
                target.position,
            )
            OMP_PARALLEL_LOOP_DO_TRANS_DYNAMIC.apply(
                target, options={"force": True}
            )

        logging.info("Member-count PARALLEL DO inserted (dynamic).")
    except TransformationError:
        logging.warning(
            "Failed to insert dynamic PARALLEL DO", exc_info=True
        )
