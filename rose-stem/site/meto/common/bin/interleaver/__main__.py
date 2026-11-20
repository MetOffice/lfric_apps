#!/usr/bin/env python3
"""Generate colocation placement files from the environment.

Reads a generic multi-executable configuration (see InterleaverConfig.from_env),
loads the PBS node allocation and writes, in CYLC_TASK_WORK_DIR:
  nodes_uniq.txt   - Unique allocated nodes, in first-seen order
  nodes_place.txt  - MPMD hostfile (one hostname per rank, MPMD order)
  cpu_bind.txt     - CPU range per rank, for inspection
  mpmd_layout.txt  - One line per executable for building mpiexec:
                        name first_rank total_ranks ranks_per_node
                        binding_start_core rank_depth omp_num_threads
"""

# sys for the exit code and for writing errors to stderr.
import sys
from pathlib import Path

# Parsing/validation of the environment-driven configuration.
from interleaver.config import InterleaverConfig, ConfigError
# The two layout strategies and their error type.
from interleaver.placements import PlacementError
# Dict containing all the available placement patterns
from interleaver.constants import PLACEMENT_PATTERNS

def main() -> int:
    """Parse config, select placement strategy, generate and write files.

    Returns:
        0 on success, 1 on error
    """
    try:
        # Step 1: read every setting from the environment, then check that
        # the settings are mutually consistent (no core overlaps, etc.).
        config = InterleaverConfig.from_env()
        config.validate()

        # Step 2: work out which physical nodes we actually have.
        # Preserve first-seen order; PBS repeats nodes and may include blanks.
        # PBS_NODEFILE lists a host once per allocated slot (so the same node
        # can appear many times). We want each node exactly once, and we keep
        # PBS's order so node 0 is consistently the "first" node of the job.
        # A set alone would lose ordering, hence the list + set pair.
        nodes = []
        seen = set()
        for line in config.pbs_nodefile.read_text().splitlines():
            node = line.strip()
            if node and node not in seen:
                nodes.append(node)
                seen.add(node)

        # An empty node file means PBS gave us nothing usable - bail early.
        if not nodes:
            raise ConfigError("PBS_NODEFILE contains no nodes")
        # Require an EXACT match between the nodes the layout needs and the
        # nodes PBS allocated. Too few and ranks can't be placed; too many
        # means the job's resource request and the rank config have drifted
        # apart and we'd be wasting (and paying for) idle nodes.
        if len(nodes) != config.required_nodes:
            raise ConfigError(
                f"Requested placement requires {config.required_nodes} nodes, "
                f"but PBS_NODEFILE contains {len(nodes)}"
            )

        # Step 3: pick the placement strategy. validate() has already rejected
        # unknown pattern names, so anything other than the asymmetric pattern
        # here is guaranteed to be "standard".
        placement = PLACEMENT_PATTERNS.get(config.placement_pattern, default=None)

        if not placement:
            raise PlacementError(f"Selected placement not found as supported: "
                                 f"{config.placement_pattern} not in "
                                 f"{PLACEMENT_PATTERNS.keys()}")

        # Confirm the chosen strategy gives every executable enough nodes.
        placement.validate()

        # Step 4: write the outputs into the Cylc task work directory, so the
        # launcher that builds the mpiexec command can find them.
        work_dir = config.cylc_task_work_dir
        uniq_path = work_dir / "nodes_uniq.txt"
        hostfile_path = work_dir / "nodes_place.txt"
        binding_path = work_dir / "cpu_bind.txt"
        layout_path = work_dir / "mpmd_layout.txt"

        # Unique node list is written directly; the other three files are
        # produced using the output from the placement.
        uniq_path.write_text('\n'.join(nodes) + '\n')

        for path, lines in (
            (hostfile_path, placement.get_total_hostfile()),
            (binding_path, placement.get_total_binding()),
            (layout_path, placement.generate_layout()),
        ):
            Path(path).write_text('\n'.join(lines) + '\n')

        return 0

    # Expected failure modes are reported as a single clean line on stderr
    # (shown in job.err) with exit status 1, so the Cylc task fails clearly
    # without a noisy traceback. Anything else is a genuine bug and is left
    # to propagate with its full traceback.
    except (ConfigError, PlacementError, OSError) as e:
        print(f"Colocation error: {e}", file=sys.stderr)
        return 1


# Run via `python -m interleaver`; propagate main()'s return as the exit code.
if __name__ == "__main__":
    sys.exit(main())
