#!/usr/bin/env bash
# (c) Crown copyright Met Office. All rights reserved.

set -euo pipefail

# 0. Check that the rank calling this script has an associated PMI_RANK var
# -------------------------------------------------------------------------
# PMI_RANK is the global MPMD rank supplied by Cray MPICH/PALS.
if [[ ! "${PMI_RANK:-}" =~ ^(0|[1-9][0-9]{0,8})$ ]]; then
  echo "LFRIC LAUNCHER ERROR: A non-negative PMI_RANK from the MPI launcher was not found, this is required for per-rank CPU binding  " >&2
  exit 1
fi

# 1. Parse variables from the command line
# ----------------------------------------
FIRST_RANK="$1"
RANK_COUNT="$2"
PPN="$3"
CPU_START="$4"
DEPTH="$5"

# Move the script's view to the next item after the required args, this should
# be the executable.
shift 5

# 2. Derive the placement of the rank
# -----------------------------------
RANK=$((PMI_RANK - FIRST_RANK))

if ((RANK < 0 || RANK >= RANK_COUNT)); then
  echo "LFRIC LAUNCHER ERROR: PMI_RANK is outside the expected executable rank range ($PMI_RANK > $RANK_COUNT)" >&2
  exit 1
fi

# Ranks are placed into CPU blocks of size `depth`. `rank % ppn` gives
# this rank's position on the current node, which is multiplied by `depth`
# to find its CPU block. `cpu_start` offsets that block past any CPUs
# reserved for executables placed earlier on the node.
BASE=$((CPU_START + RANK % PPN * DEPTH))

# Restrict this rank to `depth` consecutive CPUs, starting at `base`.
# For example, if base=8 and depth=4, the process may run on CPUs 8-11.
# `exec` then replaces this wrapper script with the requested executable.
exec taskset --cpu-list "$BASE-$((BASE + DEPTH - 1))" "$@"