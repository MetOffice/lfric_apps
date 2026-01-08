# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Transmute functions for optimisation scripts. Intension is for these to
eventually live in the PSyTran repo. See Ticket #906.
'''
import logging
import os
from itertools import dropwhile, takewhile
from typing import Sequence, Optional, Tuple, Set

from psyclone.psyir.nodes import (
    Loop,
    Call,
    Assignment,
    Reference,
    OMPParallelDirective,
    OMPDoDirective,
    OMPParallelDoDirective,
    StructureReference,
    Member,
    Literal,
    Schedule,
)
from psyclone.psyir.symbols import (
    DataSymbol,
    ContainerSymbol,
    RoutineSymbol,
    ImportInterface,
    UnsupportedFortranType,
    CHARACTER_TYPE,
)
from psyclone.transformations import (
    OMPLoopTrans,
    TransformationError,
    OMPParallelTrans,
    OMPParallelLoopTrans,
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
            if len(cluster) > 1 and not any(
                lp.ancestor(OMPParallelDirective) for lp in cluster
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
        if len(cluster) > 1 and not any(
            lp.ancestor(OMPParallelDirective) for lp in cluster
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
    if len(cluster) > 1 and not any(
        lp.ancestor(OMPParallelDirective) for lp in cluster
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

    if not hasattr(node, "explicitly_private_symbols"):
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
        if "gfortran" in s or "gcc" in s or "gnu" in s:
            return "gnu"

        # Intel (classic/oneAPI)
        if (
            "ifx" in s
            or "ifort" in s
            or "intel" in s
            or "icx" in s
            or "icc" in s
        ):
            return "intel"

        # Cray CCE
        if "cce" in s or "crayftn" in s or "cray" in s:
            return "cce"

        # NVIDIA HPC / PGI
        if (
            "nvfortran" in s
            or "nvc" in s
            or "nvhpc" in s
            or "pgfortran" in s
            or "pgi" in s
        ):
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
    first_priv_red_init(target, init_scalars)

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
            OMP_DO_LOOP_TRANS_DYNAMIC.apply(target, options={"force": True})
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
        logging.warning("Failed to insert dynamic PARALLEL DO", exc_info=True)


def first_priv_red_init(node_target, init_scalars, insert_at_start=False):
    """
    Add redundant initialisation before a Node, generally a Loop, where
    a OMP clause has FIRSTPRIVATE added by PSyclone.
    In PSyclone 3.1, many variables are made FIRSTPRIVATE; particularly
    variables that have definitions before an OpenMP parallel region, or
    that are in CodeBlocks. If these variables are uninitialised, then this
    causes the build to fail with some compilers (notably CCE).
    A fix for unnecessary FIRSTPRIVATE variables can be found in
    https://github.com/stfc/PSyclone/issues/2851, which is merged into
    PSyclone 3.2.
    Technical debt relating to this function is captured in lfric_apps:#906.

    Parameters
    ----------
    node_target : psyclone.psyir.nodes.Node
        Target Node to reference from, adds redundant initialisation before.
    init_scalars : List[str]
        List of str variable indexes to reference against.
    insert_at_start : bool, optional
        Toggles whether to insert references at the parent node
        (typically the start of the OpenMP region), or at the beginning of the
        Routine. This may be useful if e.g. your loop is initialised inside
        an if block so may or may not have an actual value at the start
        of the OpenMP region.

    Returns
    ----------
    None : Note the tree has been modified
    """
    # Ensure scalars that may be emitted as FIRSTPRIVATE have a value
    parent = node_target.parent
    # If True, add variables directly after the variable assignments
    # rather than later on in the tree, in case we have values that are
    # initialised inside conditionals later on
    if insert_at_start:
        insert_at = 0
    # Otherwise, put assignments directly before the parallel region
    else:
        insert_at = parent.children.index(node_target)
    for nm in init_scalars:  # e.g., ("jdir", "k")
        # Try and find the variable in the
        # continue rather than exit. This pattern should prevent stale
        # values from being used, since if we hit KeyError the Assignment is
        # not created.
        try:
            sym = node_target.scope.symbol_table.lookup(nm)
            # ensure character variables are initialised with CHARACTER_TYPE
            # rather than UnsupportedFortranType
            if isinstance(sym.datatype, UnsupportedFortranType):
                init = Assignment.create(
                    Reference(sym), Literal("", CHARACTER_TYPE)
                )
            else:
                init = Assignment.create(
                    Reference(sym), Literal("0", sym.datatype)
                )
            parent.children.insert(insert_at, init)
            insert_at += 1
        except KeyError:
            continue


def replace_n_threads(psyir, n_threads_var_name):
    '''
    If a scheme would use omp_get_max_threads() to determine
    how work is divided, it will often be done so through
    an omp clause. PSyclone will remove this.
    We will therefore need to be able to add it to the source.
    With a given variable name for n_threads know by the developer
    in the source, replace its initialisation (often to 1) with
    omp_get_max_threads().

    Parameters
    ----------
    psyir object : Uses whole psyir representation
    n_threads_var_name str : The name of the variable in the
                             Scheme. Set relative to the scheme.

    Returns
    ----------
    None : Note the tree has been modified
    '''

    imported_lib = False
    # Walk the schedules
    for sched in psyir.walk(Schedule):
        # Walk the Assignments
        for assign in sched.walk(Assignment):
            # If the assignment has a lhs and is a reference...
            if isinstance(assign.lhs, Reference):
                # and that lhs name is n_threads_var_name
                if assign.lhs.name == n_threads_var_name:
                    # Do this once, but only if needed
                    if imported_lib is False:
                        # Get the symbol table of the current schedule
                        symtab = sched.symbol_table
                        # Set up the omp_library symbol
                        omp_lib = symtab.find_or_create(
                            "omp_lib",
                            symbol_type=ContainerSymbol)
                        # Set up the omp_get_max_threads symbol
                        omp_get_max_threads = symtab.find_or_create(
                            "omp_get_max_threads",
                            symbol_type=RoutineSymbol,
                            # Import the reference
                            interface=ImportInterface(omp_lib))
                        imported_lib = True
                    # Replace the rhs of the reference with the
                    # omp_get_max_threads symbol
                    assign.rhs.replace_with(
                        # pylint: disable=possibly-used-before-assignment
                        Call.create(omp_get_max_threads))


def set_pure_subroutines(node, names):
    """
    Declare subroutines under `node` that are  named in
    `names` as pure, to enable parallelisation of the
    encompassing loops.
    """
    for call in node.walk(Call):
        if call.routine.name in names:
            call.routine.symbol.is_pure = True


def match_lhs_assignments(node, names):
    """
    Check if any symbol names in list `names` appear on the
    left-hand side of any assignments under `node` and
    return those names. Useful for handling, e.g., false
    dependencies and explicit variable privatisation.
    """
    lhs_names = []
    for assignment in node.walk(Assignment):
        if assignment.lhs.name in names:
            lhs_names.append(assignment.lhs.name)
    return lhs_names


def add_omp_parallel_region(
    start_node,
    end_node,
    *,
    end_offset=0,
    include_end=False,
    ignore_loops=None,
    loop_trans_options=None,
):
    """Add OMPParallelDirective around a span of nodes and OMPDoDirective around loops

    A parallel region will be created from siblings of start_node up to
    end_node.

    An end_offset may be used to add or remove nodes from the end of the region,
    and loops to be ignored can be supplied through ignore_loops.

    If end_node is not a sibling of start_node, loops and directives should be
    added up to the absolute position of end_node, although the interaction with
    offsets and include_end may be unexpected.
    """
    if ignore_loops in (None, [None]):
        ignore_loops = []

    schedule = start_node.siblings

    if end_offset != 0:
        local_idx = end_node.siblings.index(end_node)
        end_node = end_node.siblings[local_idx + end_offset]

    start_pos = start_node.abs_position
    end_pos = end_node.abs_position
    if include_end:
        end_pos += 1

    nodes_from_start = dropwhile(lambda node: node.abs_position < start_pos, schedule)
    all_nodes = list(takewhile(lambda node: node.abs_position < end_pos, nodes_from_start))

    OMP_PARALLEL_REGION_TRANS.apply(
        all_nodes,
        options={
            "node-type-check": False,
        },
    )

    for loop in start_node.parent.walk(Loop):
        # Identify each loop in the OMPParallelDirective and add OMPDoDirective to outer loops

        # Don't attempt to nest parallel directives or loops outside the parallel region
        if loop.ancestor(OMPDoDirective) is not None or loop.ancestor(OMPParallelDirective) is None:
            continue

        if loop in ignore_loops:
            continue

        # OMPDoDirective for outer loops inside OMPParallelDirective region
        try:
            OMP_DO_LOOP_TRANS_STATIC.apply(
                loop,
                options=loop_trans_options,
            )
        except TransformationError as e:
            logging.warning(e)

