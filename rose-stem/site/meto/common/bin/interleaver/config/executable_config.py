"""Configuration parsing and validation for interleaving.

The executable configuration is for a single executable though this class does
not have any knowledge of what the executables are or do.
"""

from dataclasses import dataclass
from math import ceil


from interleaver.config import ConfigError


# ---------------------------------------------------------------------------
# ExecutableConfig
# ---------------------------------------------------------------------------
# Describes how ONE executable (one MPMD segment of the mpiexec command line)
# should be laid out. All values come from {NAME}_* environment variables,
# parsed in ColocationConfig.from_env. The key mental model is:
#
#   On every node this executable uses, it occupies a contiguous run of cores
#   starting at binding_start_core. That run is split into ranks_per_node
#   slots, each rank_depth cores wide. Each MPI rank is pinned to one slot and
#   may run up to omp_num_threads OpenMP threads inside it.
#
#   e.g. binding_start_core=0, ranks_per_node=4, rank_depth=8 gives
#        rank slots 0-7, 8-15, 16-23, 24-31 on each node (32 cores total).
@dataclass
class ExecutableConfig:
    """Configuration for a single executable.

    Attributes:
        name:               Executable name (e.g., "lfric_atm", "xios", "nemo")
        total_ranks:        Total MPI ranks for this executable
        ranks_per_node:     MPI ranks placed on each node
        rank_depth:         Cores reserved for each rank (binding width)
        omp_num_threads:    OpenMP threads per rank
        binding_start_core: First core this executable uses on each node
    """
    name: str
    total_ranks: int
    ranks_per_node: int
    rank_depth: int
    omp_num_threads: int
    binding_start_core: int = 0

    @property
    def env_prefix(self) -> str:
        """Prefix of this executable's environment variables."""
        # Environment variables are conventionally upper case, so an
        # executable called "lfric_atm" is configured via LFRIC_ATM_* vars.
        # This is also used in error messages so users see the exact variable
        # name they need to change.
        return self.name.upper()

    @property
    def cores_per_node(self) -> int:
        """Cores this executable occupies on each node."""
        # Each rank reserves rank_depth cores, and there are ranks_per_node of
        # them on a node, so the footprint on a single node is the product.
        # Note this is the footprint on a FULL node; the last node may be
        # partially filled if total_ranks isn't a multiple of ranks_per_node,
        # but we still reserve the full width so the layout is uniform.
        return self.ranks_per_node * self.rank_depth

    @property
    def binding_end_core(self) -> int:
        """One past the last core this executable uses on each node."""
        # Exclusive upper bound (Python-range style). The last core actually
        # used is binding_end_core - 1. Using an exclusive bound makes the
        # overlap test in ColocationConfig.validate a simple interval check
        # and lets the next shared executable start exactly here.
        return self.binding_start_core + self.cores_per_node

    @property
    def effective_end_core(self) -> int:
        return self.effective_start_core + self.cores_per_node

    @property
    def needed_nodes(self) -> int:
        """Nodes needed to hold all ranks; the final node may be partial."""
        # Round UP: any leftover ranks still need a whole node to live on.
        # e.g. 10 ranks at 4 per node -> nodes hold 4, 4, 2 -> 3 nodes.
        return ceil(self.total_ranks / self.ranks_per_node)

    def validate(self, node_cores: int) -> None:
        """Validate this executable against the node size.

        Raises:
            ConfigError: If validation fails
        """
        # Cache the prefix so error messages can name the exact env vars.
        prefix = self.env_prefix

        # A rank is pinned to rank_depth cores. Asking for more OpenMP threads
        # than that would oversubscribe the cores (threads fighting over the
        # same core), which silently destroys performance. Fewer threads than
        # cores is allowed - e.g. deliberately leaving cores idle for memory
        # bandwidth reasons - so only the "greater than" case is an error.
        if self.omp_num_threads > self.rank_depth:
            raise ConfigError(f"{prefix}_OMP_NUM_THREADS must not exceed {prefix}_RANK_DEPTH")

        # The executable's core range must physically exist on the node. If
        # binding_end_core (exclusive) is beyond the node's core count, some
        # ranks would be bound to non-existent CPUs and mpiexec would fail or,
        # worse, wrap around onto cores used by something else.
        # This also catches cases where cores per node > node_cores
        if self.binding_end_core > node_cores:
            raise ConfigError(
                f"{self.name} requires cores {self.binding_start_core}-"
                f"{self.binding_end_core - 1}, but CORES_PER_NODE is {node_cores}"
            )
