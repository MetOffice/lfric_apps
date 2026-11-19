# XIOS Colocation Module

A Python 3.12 module for managing XIOS server colocation with LFRIC atmosphere ranks. Replaces shell scripts for hostfile and CPU binding generation with a well-structured, testable, and extensible implementation.

## Overview

This module handles the placement and binding of LFRIC atmosphere and XIOS server ranks across HPC nodes. It supports:

- **Standard interleaving**: All nodes interleaved uniformly
- **Asymmetric layout**: Two-block strategy with different node configurations
  - Interleaved block: Reduced atmosphere cores + XIOS ranks
  - Full-population block: All cores for atmosphere only

## Architecture

```
bin/colocation/
├── __init__.py           # Package metadata
├── __main__.py           # CLI entry point (python -m colocation)
├── config.py             # Configuration parsing & validation
├── placement.py          # Placement strategy classes
└── hostfile.py           # Hostfile & binding file generation
```

## Entry Points

### Colocation Module (Main)

Generate hostfile and binding files for colocated execution.

```bash
python -m colocation
```

**Exit codes:**
- 0: Success
- 1: Configuration or validation error
- 2: File I/O error
- 3: Runtime error

**Environment variables (required):**
- `EXECUTABLES` - An ordered list of executables to be included in the MPMD line.
                  e.g. 'lfric_atm xios nemo'
For each name in EXECUTABLES:
- `<EXEC>_RANKS_PER_NODE` - Mumber of processes per node
- `<EXEC>_RANK_DEPTH` - Number of cores between each process
- `<EXEC>_NUM_THREADS` - Number of OpenMP threads (must not exceed `<EXEC>_RANK_DEPTH`)

- `CORES_PER_NODE` - Total cores per node (architecture)
- `TOTAL_RANKS` - Total atmosphere ranks
- `OMP_NUM_THREADS` - Atmosphere thread depth (MPI process depth)
- `XIOS_SERVER_RANKS` - Total XIOS ranks (validated against configuration)
- `PBS_NODEFILE` - Path to PBS node allocation file
- `CYLC_TASK_WORK_DIR` - Work directory for output files

**Environment variables (optional, generic mode only):**
- `PLACEMENT_PATTERN` - 'standard' or 'asymmetric_two_block' (default: standard)
- `SEPARATE_EXECUTABLES` - Executables placed on their own node block

**Output files** (in `CYLC_TASK_WORK_DIR`):
- `my_nodes_uniq.txt` - Unique allocated nodes (first-seen order)
- `my_nodes_place.txt` - Hostfile (one hostname per rank, MPMD order)
- `my_cpu_bind.txt` - CPU binding ranges (one range per rank, format: "0-3")

## Supported Hardware Configurations

### Milan (128-core)
Standard colocation on all nodes:
- 112 cores/node for atmosphere
- 4 XIOS ranks/node × 4 cores deep = 16 cores XIOS
- 28 atmosphere ppn (112 / 4 threads)
- 4 XIOS ppn

**Example:** 40 nodes, 1120 atmosphere ranks, 160 XIOS ranks

### Genoa (192-core)
Asymmetric striped layout:
- **Colocated block (104 nodes):** 168 cores atmosphere + 24 cores XIOS
  - 42 atmosphere ppn (168 / 4 threads)
  - 6 XIOS ppn × 4 cores = 24 cores
  - 4368 atmosphere ranks + 624 XIOS ranks

- **Full-population block (105 nodes):** 192 cores atmosphere only
  - 48 atmosphere ppn (192 / 4 threads)
  - 5040 atmosphere ranks

**Total:** 209 nodes, 9408 atmosphere ranks, 624 XIOS ranks

## Testing

Tests live in `bin/colocation/tests/`, one file per module and one test class
per class under test:

| File | Tests |
|------|-------|
| `test_config.py` | `ExecutableConfig`, `ColocationConfig` |
| `test_placement.py` | `PlacementStrategy`, `StandardPlacement`, `AsymmetricStripedPlacement` |
| `test_hostfile.py` | `HostfileGenerator` |
| `test_main.py` | `main()` (`python3 -m colocation`) |
| `test_launcher.py` | `launch_colocated.sh` and `bind_rank.sh`, with a mock `mpiexec` |

Most tests share one small example of three made-up executables on three
16-core nodes, drawn out in `tests/conftest.py`. Run them with:

```bash
python3 -m pytest bin/colocation/tests
```

## Key Design Decisions

### 1. Placement Strategies (Abstract Base Class)
- Decouples placement logic from hostfile/binding generation
- Allows adding new patterns (e.g., alternating striped) without modifying core code
- Each strategy encapsulates node-to-rank mapping

### 2. Environment-Driven Configuration
- Maintains compatibility with current workflow
- All config read from environment variables at startup
- Validation happens before placement to catch errors early

### 3. Two-Block Asymmetric Striped
- Colocated and full-population nodes are contiguous blocks
- Simplifies hostfile generation and matches SMP allocation ordering
- Validates mathematically (rank counts must add to total)

### 4. Asymmetric Striped XIOS Validation Skip
- Standard formula (nodes × ppn) doesn't apply to asymmetric layout
- Pattern-aware validation: only validates when applicable
- Prevents false positives for two-block strategy

### 5. MPMD Hostfile Order
- Atmosphere ranks first, XIOS ranks second
- Matches mpiexec MPMD invocation: ` : <xios_block>`
- Binding file follows same order

## Implementation Notes

### Configuration Validation Flow
1. Load all environment variables (type checking)
2. Validate constraints:
   - OMP_NUM_THREADS divides ATMOS_CORES_PER_NODE
   - Total cores ≤ CORES_PER_NODE
   - XIOS_SERVER_RANKS matches expected (pattern-dependent)
   - Required files exist
3. Select placement strategy based on PLACEMENT_PATTERN

### Placement Strategy Selection
```python
if config.placement_pattern == "asymmetric_two_block":
    placement = AsymmetricStripedPlacement(config, nodes)
else:
    placement = StandardPlacement(config, nodes)
```

This logic is in `__main__.py`.

### Error Messages
All errors include:
- What failed (clear description)
- Why it failed (constraint violated, file not found, etc.)
- How to fix it (suggested values, required env vars, etc.)

Example:
```
ERROR: OMP_NUM_THREADS (4) must divide ATMOS_CORES_PER_NODE (112) exactly
  Valid options: 1, 2, 4, 7, 8, 14, 16, 28, 56, 112
  Set OMP_NUM_THREADS to one of these values
```

## Future Enhancements

1. **Alternating Striped Pattern**: Every Nth node colocated, configurable stripe width
2. **Dynamic XIOS Sizing**: Calculate optimal XIOS ppn based on I/O characteristics
3. **Node Heterogeneity**: Support mixed node types in single allocation
4. **Performance Profiling**: Built-in hooks for measuring placement effectiveness
5. **Validation Tool**: Standalone checker for comparing actual vs. expected rank placement

## Backwards Compatibility

- Replaces `mixed_hostfile.sh` and `colocation_config.sh`; `bin/launch_lfric.sh`
  calls `python3 -m colocation`, which validates the launch settings and writes
  the placement files
- `bind_rank.sh` (per-rank `taskset` binding) is still used by the launcher
- The allocation must contain exactly the required number of unique nodes
- All environment variable names unchanged
- Identical output file format and locations

## Development

### Adding a New Placement Strategy

1. Create subclass of `PlacementStrategy` in `placement.py`
2. Implement the two abstract methods, each returning a mapping of executable
   name to one entry per rank:
   - `get_hostfile_entries()`
   - `get_binding_entries()`
3. Update the strategy selection in `main()` in `__main__.py` to instantiate the new strategy
4. Document in this README

### Code Quality

- Type hints on all functions and classes
- Docstrings for all public methods
- Clear error messages with diagnostics
- 100% test coverage for non-trivial logic
- No external dependencies (Python stdlib only)

## References

- LFRIC model documentation
- XIOS I/O server documentation
- HPC job scheduling and MPI (PALS/Cray MPICH)
- Python 3.12 dataclasses and typing
