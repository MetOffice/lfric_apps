"""Configuration parsing and validation for colocation.

The configuration is generic: any number of named executables, each sized by
its own environment variables. Nothing here knows what the executables are.
"""

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from interleaver.config import ExecutableConfig, ConfigError
from interleaver.constants import PLACEMENT_PATTERNS


# ---------------------------------------------------------------------------
# InterleaverConfig
# ---------------------------------------------------------------------------
# The full job-wide configuration: node size, every executable (in the order
# they appear on the mpiexec MPMD command line), which placement strategy to
# use, and where to read the PBS allocation from / write results to.
#
# Executables fall into two groups:
#   - "shared":   run side by side on the same nodes, on disjoint core ranges.
#   - "separate": listed in SEPARATE_EXECUTABLES; get dedicated nodes of their
#                 own (only supported by the asymmetric_two_block pattern).
@dataclass
class InterleaverConfig:
    """Multi-executable colocation configuration.

    Attributes:
        cores_per_node: Total cores available per node
        executables: Executables by name, in MPMD order
        placement_pattern: Placement strategy (see PLACEMENT_PATTERNS)
        pbs_nodefile: Path to PBS node allocation file
        cylc_task_work_dir: Work directory for output files
        separate_executables: Executables placed on their own nodes
    """
    cores_per_node: int
    executables: dict[str, ExecutableConfig]
    placement_pattern: str
    pbs_nodefile: Path
    cylc_task_work_dir: Path
    separate_executables: list[str] = field(default_factory=list)

    @property
    def mpmd_order(self) -> list[str]:
        """Executable names in MPMD order."""
        # dicts preserve insertion order (Python 3.7+), and from_env() inserts
        # executables in the order they appear in EXECUTABLES, so the keys are
        # already in MPMD order. This order MUST match the order of the
        # segments in the mpiexec command, because MPI assigns global ranks
        # sequentially across segments and the hostfile is read in that order.
        return list(self.executables)

    @property
    def shared_executables(self) -> list[str]:
        """Executables that share nodes with each other."""
        # Everything not explicitly marked separate is shared. Preserves MPMD
        # order because we iterate the ordered executables dict.
        return [n for n in self.executables if n not in self.separate_executables]

    @property
    def interleaved_nodes(self) -> int:
        """Nodes in the shared block."""
        # The shared block must be big enough for the hungriest shared
        # executable. (validate() additionally insists they all need the SAME
        # number of nodes, so in a valid config max() == each of them.)
        # default=0 avoids a ValueError from max() on an empty sequence; an
        # empty shared list is rejected later by validate() with a clear error.
        return max((self.executables[n].needed_nodes for n in self.shared_executables), default=0)

    @property
    def required_nodes(self) -> int:
        """Nodes the placement needs in total."""
        # Total = one shared block, plus a dedicated block for every separate
        # executable laid end to end. __main__ compares this against the
        # number of unique nodes PBS actually gave us and refuses to continue
        # on a mismatch, catching resource-request vs. config drift early.
        return self.interleaved_nodes + sum(
            self.executables[n].needed_nodes for n in self.separate_executables
        )

    @staticmethod
    def get_int(name, environ, default=None, allow_zero=False):
        '''
        Accept plain decimal digits only, matching the workflow contract.
        We deliberately do NOT just call int(v): int() happily accepts
        " 4", "+4", "0004", "4_000" and other forms that almost certainly
        indicate a templating mistake upstream. Being strict here means a
        broken Jinja2 substitution in the suite is caught immediately.
        '''
        v = environ.get(name)
        # Variable absent and a default exists -> use the default. Note
        # that an empty-string value is NOT treated as absent; it falls
        # through and fails the regex, because an empty value usually
        # means the workflow tried (and failed) to set it.
        if v is None and default is not None:
            return default
        # Up to 9 digits keeps values comfortably within sane bounds and
        # rejects absurd numbers. Leading zeros are disallowed (except a
        # lone "0" when zero is permitted, e.g. a start core of 0).
        pattern = r"(0|[1-9][0-9]{0,8})" if allow_zero else r"[1-9][0-9]{0,8}"
        if v is None or not re.fullmatch(pattern, v):
            kind = "non-negative" if allow_zero else "positive"
            raise ConfigError(f"{name} must be supplied as a {kind} integer by the workflow")
        return int(v)

    @staticmethod
    def get_path(name, environ):
        '''
        Fetch a path-valued variable and insist it both is set and exists.
        The explicit "not environ.get(name)" check matters: Path("") is ".",
        which always exists, so an unset variable would otherwise silently
        resolve to the current directory.
        '''
        p = Path(environ.get(name, ""))
        if not environ.get(name) or not p.exists():
            raise ConfigError(f"{name} must name an existing path")
        return p


    @classmethod
    def from_env(cls, environ=None) -> "InterleaverConfig":
        """Parse configuration from environment variables.

        Global variables:
          - EXECUTABLES: Space-separated executable names, in MPMD order
          - CORES_PER_NODE: Total cores per node
          - PLACEMENT_PATTERN: (optional) one of PLACEMENT_PATTERNS
          - SEPARATE_EXECUTABLES: (optional) executables on their own nodes
          - PBS_NODEFILE, CYLC_TASK_WORK_DIR

        For each executable NAME (variables use NAME in upper case):
          - {NAME}_TOTAL_RANKS: Total ranks
          - {NAME}_RANKS_PER_NODE: Ranks on each node
          - {NAME}_RANK_DEPTH: Cores per rank
          - {NAME}_OMP_NUM_THREADS: (optional) defaults to {NAME}_RANK_DEPTH
          - {NAME}_BINDING_START_CORE: (optional) defaults to the first core
            after the preceding shared executable, or 0 on separate nodes

        Raises:
            ConfigError: If a variable is missing or invalid
        """
        # ======================================================================
        # 0.0: set up some shortcuts
        # ======================================================================
        # static methods
        get_int = InterleaverConfig.get_int
        get_path = InterleaverConfig.get_path

        # Allow a plain dict to be injected (handy for unit tests); otherwise
        # read the real process environment set up by Cylc/Rose.
        environ = os.environ if environ is None else environ

        # ======================================================================
        # 1.0: Parse environment variables
        # ======================================================================
        # ----------------------------------------------------------------------
        # 1.1: Split EXECUTABLES list
        # ----------------------------------------------------------------------
        # EXECUTABLES is the master list and defines MPMD order. split() with
        # no argument handles any amount/kind of whitespace and drops empties.
        names = environ.get("EXECUTABLES", "").split()
        if not names:
            raise ConfigError("EXECUTABLES must list at least one executable")
        # Each name becomes part of environment variable names ({NAME}_...),
        # so it must be a valid identifier: start with a letter, then letters,
        # digits or underscores. Hyphens, dots etc. would produce variable
        # names the shell cannot export.
        for name in names:
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
                raise ConfigError(f"Executable name {name!r} must be usable in "
                                  f"environment variable names")
        # Check uniqueness on the UPPER-CASED form: "xios" and "XIOS" are
        # different names but would both read XIOS_* variables, so they would
        # silently share configuration. Treat that as a duplicate.
        prefixes = [n.upper() for n in names]
        if len(set(prefixes)) != len(prefixes):
            raise ConfigError("EXECUTABLES must not repeat a name")

        # ----------------------------------------------------------------------
        # 1.1: Get executables not using interleaving
        # ----------------------------------------------------------------------
        # Optional list of executables that should get their own nodes rather
        # than sharing. Validity (known names, compatible pattern) is checked
        # later in validate(), once the full config object exists.
        separate = environ.get("SEPARATE_EXECUTABLES", "").split()

        # ======================================================================
        # 2.0: Create a dictionary of executable configurations
        #      - Name, depth, ppn, total ranks, OMP, core to start binding
        # ======================================================================
        # Build one ExecutableConfig per name, preserving MPMD order.
        executables = {}
        # Running cursor for automatic core packing of shared executables:
        # the first shared executable starts at core 0, the next starts where
        # the previous one ended, and so on. This gives a sensible,
        # non-overlapping default without the user specifying start cores.
        next_shared_core = 0
        for name in names:
            prefix = name.upper()
            # Read depth first because it is also the default thread count.
            depth = get_int(f"{prefix}_RANK_DEPTH", environ)
            ranks_per_node = get_int(f"{prefix}_RANKS_PER_NODE", environ)
            # Separate executables own their nodes outright, so they naturally
            # start at core 0. Shared ones are packed after their predecessors.
            default_start = 0 if name in separate else next_shared_core
            exe = ExecutableConfig(
                name=name,
                total_ranks=get_int(f"{prefix}_TOTAL_RANKS", environ),
                ranks_per_node=ranks_per_node,
                rank_depth=depth,
                omp_num_threads=get_int(f"{prefix}_OMP_NUM_THREADS", environ, default=depth),
                binding_start_core=get_int(
                    f"{prefix}_BINDING_START_CORE", environ, default=default_start, allow_zero=True
                ),
            )
            # Advance the packing cursor only for shared executables; separate
            # ones live on other nodes and don't consume shared-node cores.
            # If the user overrode BINDING_START_CORE, packing continues from
            # wherever their override ended.
            if name not in separate:
                next_shared_core = exe.binding_end_core
            executables[name] = exe

        # ======================================================================
        # 3.0: Create an InterleaverConfig instance with all the information
        # ======================================================================
        # Assemble the final object. Note this only PARSES; semantic checks
        # (overlaps, node counts, pattern compatibility) happen in validate(),
        # which the caller must invoke explicitly.
        return cls(
            cores_per_node=get_int("CORES_PER_NODE", environ),
            executables=executables,
            placement_pattern=environ.get("PLACEMENT_PATTERN", "standard"),
            pbs_nodefile=get_path("PBS_NODEFILE", environ),
            cylc_task_work_dir=get_path("CYLC_TASK_WORK_DIR", environ),
            separate_executables=separate,
        )

    def validate(self) -> None:
        """Validate all configuration constraints.

        Checks:
        - Placement pattern is recognised
        - SEPARATE_EXECUTABLES names known executables and suits the pattern
        - Each executable fits within the node
        - Executables sharing nodes do not share cores
        - Executables sharing nodes span every node of their block, so each
          node holds the same mix of executables

        Raises:
            ConfigError: If validation fails
        """
        # ----------------------------------------------------------------------
        # 1. Unknown pattern names are rejected up front; otherwise __main__
        #    would silently fall back to StandardPlacement for a typo.
        # ----------------------------------------------------------------------
        if self.placement_pattern not in PLACEMENT_PATTERNS:
            raise ConfigError(
                f"PLACEMENT_PATTERN must be one of {', '.join(PLACEMENT_PATTERNS)}; "
                f"got {self.placement_pattern}"
            )

        # ----------------------------------------------------------------------
        # 2. Every SEPARATE_EXECUTABLES entry must refer to something listed
        #    in EXECUTABLES; a typo here would otherwise be ignored and the
        #    executable would quietly end up shared.
        # ----------------------------------------------------------------------
        for name in self.separate_executables:
            if name not in self.executables:
                raise ConfigError(f"SEPARATE_EXECUTABLES references unknown executable: {name}")

        # ----------------------------------------------------------------------
        # 3. Only the asymmetric pattern knows how to give executables their
        #    own nodes. The standard pattern puts everything everywhere, so a
        #    separate list would be meaningless (and misleading) there.
        # ----------------------------------------------------------------------
        if self.separate_executables and self.placement_pattern != "asymmetric_two_block":
            raise ConfigError("SEPARATE_EXECUTABLES requires PLACEMENT_PATTERN=asymmetric_two_block")

        # ----------------------------------------------------------------------
        # 4. There must be something in the shared block; a layout made
        #    entirely of separate executables has no "first block" to anchor.
        # ----------------------------------------------------------------------
        if not self.shared_executables:
            raise ConfigError("At least one executable must share nodes")

        # ----------------------------------------------------------------------
        # 5. Per-executable checks: threads fit in rank depth, and the core
        #    range fits on the node. See ExecutableConfig.validate.
        # ----------------------------------------------------------------------
        for exe in self.executables.values():
            exe.validate(self.cores_per_node)

        # ----------------------------------------------------------------------
        # 6. Shared executables live on the same nodes, so their core ranges
        #    must not overlap or two ranks would be pinned to the same cores.
        #    Compare every unordered pair (i < j) exactly once. Two half-open
        #    intervals [a, b) and [c, d) overlap iff a < d and c < b.
        #    Separate executables are excluded: they're on different nodes.
        # ----------------------------------------------------------------------
        shared = [self.executables[n] for n in self.shared_executables]
        for i, first in enumerate(shared):
            for second in shared[i + 1:]:
                if (first.binding_start_core < second.binding_end_core and
                        second.binding_start_core < first.binding_end_core):
                    raise ConfigError(
                        f"{first.name} (cores {first.binding_start_core}-{first.binding_end_core - 1}) "
                        f"overlaps {second.name} (cores {second.binding_start_core}-"
                        f"{second.binding_end_core - 1})"
                    )

        # ----------------------------------------------------------------------
        # 7. Every shared executable must need exactly the same number of
        #    nodes. If one needed fewer, some nodes in the shared block would
        #    carry only part of the mix, making per-node load uneven and
        #    breaking the "every node looks the same" assumption that the
        #    core-range layout relies on. The error lists each executable's
        #    need and the env vars driving it, so the fix is obvious.
        # ----------------------------------------------------------------------
        nodes = self.interleaved_nodes
        if any(exe.needed_nodes != nodes for exe in shared):
            needs = ", ".join(
                f"{exe.name} {exe.needed_nodes} "
                f"({exe.env_prefix}_TOTAL_RANKS={exe.total_ranks}, "
                f"{exe.env_prefix}_RANKS_PER_NODE={exe.ranks_per_node})"
                for exe in shared
            )
            raise ConfigError(
                f"Executables sharing nodes must span the same number of nodes; needs: {needs}"
            )
