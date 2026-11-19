from interleaver.placements import PlacementStrategy


class AsymmetricStripedPlacement(PlacementStrategy):
    """Asymmetric layout: partition nodes into contiguous blocks.

    Shared executables run together on the first block. Each executable in
    SEPARATE_EXECUTABLES then gets its own block, in MPMD order.

    For example, with executables A, B and C, and C separate:
    - Shared block (nodes 0-103): A + B
    - Separate block (nodes 104-109): C
    """

    def get_node_blocks(self) -> dict[str, list[str]]:
        # Size of the shared block = nodes needed by the shared executables
        # (validated to be equal for all of them).
        shared_nodes = self.config.interleaved_nodes
        # All shared executables get the SAME leading slice of nodes.
        blocks = {name: self.nodes[:shared_nodes] for name in self.config.shared_executables}
        # Cursor into the node list: separate blocks begin right after the
        # shared block and are handed out contiguously.
        start = shared_nodes
        # Walk in MPMD order so separate blocks appear in the node list in the
        # same order as their segments on the mpiexec command line.
        for name in self.get_executable_order():
            if name in self.config.separate_executables:
                needed = self.config.executables[name].needed_nodes
                # Give this executable exactly the nodes it needs, then move
                # the cursor past them for the next separate executable.
                blocks[name] = self.nodes[start:start + needed]
                start += needed
        return blocks
