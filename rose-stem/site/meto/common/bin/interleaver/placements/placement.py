"""Placement strategies that map ranks to nodes and CPU cores."""

# ABC/abstractmethod let us define a common interface (PlacementStrategy) that
# concrete strategies must fill in, while sharing all the rank/core maths.
from abc import ABC, abstractmethod
# Only needed for type hints; strategies read sizing from the parsed config.
from interleaver.config import InterleaverConfig

# ---------------------------------------------------------------------------
# PlacementStrategy (abstract base)
# ---------------------------------------------------------------------------
# Splits the problem in two:
#   1. WHICH NODES each executable uses  -> get_node_blocks(), per strategy.
#   2. HOW RANKS FILL THOSE NODES/CORES  -> shared logic below, identical for
#      every strategy: ranks fill nodes in order, ranks_per_node at a time,
#      and within a node each rank takes the next rank_depth-wide core slot.
# Subclasses therefore only ever need to implement get_node_blocks().
class PlacementStrategy(ABC):
    """Abstract base for rank placement strategies.

    A strategy only decides which nodes each executable runs on. Ranks fill
    those nodes in order, ranks_per_node at a time, and each rank is bound to
    rank_depth cores starting at the executable's binding_start_core.

    Attributes:
        config: ColocationConfig instance
        nodes: List of unique node names from PBS allocation
    """
    def __init__(self, config: InterleaverConfig, nodes: list[str]) -> None:
        # Keep the config (sizes, order) and the de-duplicated, order-preserved
        # list of node hostnames from PBS. Node order matters: it determines
        # which physical nodes form which block.
        self.config = config
        self.nodes = nodes

    @abstractmethod
    def get_node_blocks(self) -> dict[str, list[str]]:
        """Return the nodes each executable runs on.

        Returns:
            {executable_name: [hostname1, hostname2, ...]}
        """
        pass

    def get_executable_order(self) -> list[str]:
        """Return executables in MPMD order."""
        # Delegates to the config so every component agrees on one ordering.
        return self.config.mpmd_order

    def get_hostfile_entries(self) -> dict[str, list[str]]:
        """Return one hostname per rank for each executable."""
        result = {}
        for name, block in self.get_node_blocks().items():
            exe = self.config.executables[name]
            # Local rank r lives on node 'r // ranks_per_node' of this
            # executable's block. With 4 ranks per node, ranks 0-3 -> block[0],
            # 4-7 -> block[1], etc. Integer division gives the node index.
            # The final node may receive fewer than ranks_per_node ranks.
            result[name] = [block[rank // exe.ranks_per_node] for rank in range(exe.total_ranks)]
        return result

    def get_binding_entries(self) -> dict[str, list[str]]:
        """Return one CPU range per rank for each executable."""
        result = {}
        for name, exe in self.config.executables.items():
            entries = []
            for rank in range(exe.total_ranks):
                # rank % ranks_per_node is the rank's slot index WITHIN its
                # node (0..ranks_per_node-1), mirroring the node assignment in
                # get_hostfile_entries. Each slot is rank_depth cores wide and
                # slots start at binding_start_core. Note: % and * share
                # precedence and associate left-to-right, so this reads as
                # start + ((rank % rpn) * depth), which is what we want.
                cpu_start = exe.binding_start_core + rank % exe.ranks_per_node * exe.rank_depth
                # Inclusive range string, e.g. "8-15", as consumed by CPU
                # binding tools; hence the -1 on the end core.
                entries.append(f"{cpu_start}-{cpu_start + exe.rank_depth - 1}")
            result[name] = entries
        return result

    def get_total_hostfile(self) -> list[str]:
        """Return concatenated hostfile in MPMD order."""
        entries = self.get_hostfile_entries()
        # Flatten executable-by-executable in MPMD order. MPI numbers global
        # ranks consecutively across MPMD segments, so line N of the hostfile
        # must be the host for global rank N.
        return [host for name in self.get_executable_order() for host in entries[name]]

    def get_total_binding(self) -> list[str]:
        """Return concatenated binding in MPMD order."""
        entries = self.get_binding_entries()
        # Same flattening as the hostfile so line N lines up with global rank
        # N in both files.
        return [cpus for name in self.get_executable_order() for cpus in entries[name]]

    def generate_layout(self) -> list[str]:
        """Return one line per executable describing its MPMD segment.

        Fields: name, first global rank, total ranks, ranks per node,
        binding start core, rank depth, OpenMP threads.
        """
        lines = []
        # Global rank of the first rank in the current MPMD segment. MPI
        # numbers ranks consecutively across segments, so each executable's
        # first rank is the running total of all ranks before it.
        first_rank = 0
        for name in self.get_executable_order():
            exe = self.config.executables[name]
            # Space-separated so a shell launcher can parse it with a plain
            # `read name first total rpn start depth threads` loop. Treat the
            # field order as a contract with any launcher - do not reorder.
            lines.append(
                f"{name} {first_rank} {exe.total_ranks} {exe.ranks_per_node} "
                f"{exe.binding_start_core} {exe.rank_depth} {exe.omp_num_threads}"
            )
            # Advance to where the next segment's ranks begin.
            first_rank += exe.total_ranks
        return lines

    def validate(self) -> None:
        """Verify every executable has enough nodes for its ranks."""
        for name, block in self.get_node_blocks().items():
            exe = self.config.executables[name]
            # If a strategy hands out too few nodes, get_hostfile_entries would
            # raise an opaque IndexError. Catch it here with a clear message.
            # Having MORE nodes than needed is fine; extras simply go unused.
            if len(block) < exe.needed_nodes:
                raise PlacementError(
                    f"{name} needs {exe.needed_nodes} nodes but was given {len(block)}"
                )
