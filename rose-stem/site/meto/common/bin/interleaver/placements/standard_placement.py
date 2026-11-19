from interleaver.placements import PlacementStrategy


class StandardPlacement(PlacementStrategy):
    """Standard colocation: all executables on all nodes.

    Different binding_start_core values keep executables on separate cores.
    """

    def get_node_blocks(self) -> dict[str, list[str]]:
        # Every executable gets the full node list. They never collide because
        # ColocationConfig.validate has already ensured their core ranges are
        # disjoint, so on each node they occupy different cores.
        return {name: self.nodes for name in self.config.executables}