from interleaver.placements import StandardPlacement, AsymmetricStripedPlacement

# The complete list of placement strategies that PLACEMENT_PATTERN may name.
#   "standard"             - every executable runs on every node, and they are
#                            kept apart purely by binding them to different
#                            core ranges on each node (see StandardPlacement).
#   "asymmetric_two_block" - the "shared" executables run together on a first
#                            block of nodes, and any SEPARATE_EXECUTABLES each
#                            get their own contiguous block of nodes after it
#                            (see AsymmetricStripedPlacement).
PLACEMENT_PATTERNS = {"standard": StandardPlacement,
                      "asymmetric_two_block": AsymmetricStripedPlacement}

# ==============================================================================
# Custom Error types
# ==============================================================================
# A dedicated exception type for anything wrong with the configuration that
# the workflow supplied. It subclasses ValueError because, semantically, these
# are "bad value" problems. __main__ catches this specifically so that config
# problems are reported as a clean one-line message (and exit code 1) instead
# of a full Python traceback in the Cylc job log.
class ConfigError(ValueError):
    """Configuration validation error."""
    pass

# Raised when a strategy cannot fit the requested ranks onto the nodes it has
# been given. Distinct from ConfigError so logs make clear whether the input
# was malformed (ConfigError) or the layout itself is impossible (this).
class PlacementError(RuntimeError):
    """Placement validation error."""
    pass
