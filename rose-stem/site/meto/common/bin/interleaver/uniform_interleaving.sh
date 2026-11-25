#!/usr/bin/env bash
# (c) Crown copyright Met Office. All rights reserved.

set -euo pipefail

function validate_interleaving_settings () {

    local USED_CORES=0
    local name

    # Validate the uniform interleaving configuration
    for name in ${EXEC_LABELS}; do
        # Construct the actual environment-variable names
        local RANKS_PER_NODE_VAR="${name^^}_RANKS_PER_NODE"
        local OMP_NUM_THREADS_VAR="${name^^}_OMP_NUM_THREADS"
        local RANK_DEPTH_VAR="${name^^}_RANK_DEPTH"

        # Create references to those environment variables
        local -n EXEC_RANKS_PER_NODE="${RANKS_PER_NODE_VAR}"
        local -n EXEC_OMP_NUM_THREADS="${OMP_NUM_THREADS_VAR}"
        local -n EXEC_RANK_DEPTH="${RANK_DEPTH_VAR}"

        # Check that OMP_NUM_THREADS divides RANKS_PER_NODE exactly
        if (( EXEC_RANKS_PER_NODE % EXEC_OMP_NUM_THREADS != 0 )); then
            echo "LFRIC LAUNCHER ERROR: ${OMP_NUM_THREADS_VAR} must divide ${RANKS_PER_NODE_VAR} exactly when using \"uniform\" interleaving" >&2
            exit 1
        fi

        # Check that rank depth is not smaller than the thread count
        if (( EXEC_OMP_NUM_THREADS > EXEC_RANK_DEPTH )); then
            echo "LFRIC LAUNCHER ERROR: ${OMP_NUM_THREADS_VAR} must not exceed ${RANK_DEPTH_VAR} when using \"uniform\" interleaving" >&2
            exit 1
        fi

        # Update the total number of used cores
        (( USED_CORES += EXEC_RANKS_PER_NODE * EXEC_RANK_DEPTH ))
    done

    # Check that required cores do not exceed available cores
    if (( USED_CORES > CORES_PER_NODE )); then
        echo "LFRIC LAUNCHER ERROR: The total number of requested cores for interleaving exceeds the number of available cores (${USED_CORES} > ${CORES_PER_NODE})" >&2
        exit 1
    fi
}


function generate_interleaving_vars () {

    local SUM_TOTAL_RANKS=0
    local SUM_TOTAL_CORES=0
    local USED_CORES=0
    local name
    local first_exec="True"

    for name in ${EXEC_LABELS}; do
        # Extract info for each executable
        # 1. Pre-existing variables
        # -------------------------
        local RANKS_PER_NODE_VAR="${name^^}_RANKS_PER_NODE"
        local RANK_DEPTH_VAR="${name^^}_RANK_DEPTH"
        local TOTAL_RANKS_VAR="${name^^}_TOTAL_RANKS"

        local -n EXEC_RANKS_PER_NODE="${RANKS_PER_NODE_VAR}"
        local -n EXEC_RANK_DEPTH="${RANK_DEPTH_VAR}"
        local -n EXEC_TOTAL_RANKS="${TOTAL_RANKS_VAR}"

        # 2. Variables to be created
        # --------------------------
        local PPN_VAR="${name^^}_PPN"
        local START_RANK_VAR="${name^^}_START_RANK"
        local CPU_START_VAR="${name^^}_CPU_START"

        local -n EXEC_PPN="${PPN_VAR}"
        local -n EXEC_START_RANK="${START_RANK_VAR}"
        local -n EXEC_CPU_START="${CPU_START_VAR}"

        # 3. Create vars
        # --------------
        # Export the target variable's PPN by name, not its value
        EXEC_PPN=$(( EXEC_RANKS_PER_NODE * EXEC_RANK_DEPTH ))
        export "${PPN_VAR}"

        # The first exec will:
        # - start at rank 0 in the PMI list
        # - start from rank 0 on each node
        EXEC_START_RANK=$SUM_TOTAL_RANKS
        EXEC_CPU_START=$USED_CORES

        # Update the number of ranks used -> this will be used to set the
        # EXEC_START_RANK for the next exec in the list
        (( SUM_TOTAL_RANKS += EXEC_TOTAL_RANKS ))

        # Update the total number of used cores -> this will be used to set the
        # EXEC_CPU_START for the next exec in the list
        (( USED_CORES += EXEC_RANKS_PER_NODE * EXEC_RANK_DEPTH ))

        # Add to the number of required nodes
        (( SUM_TOTAL_CORES += EXEC_TOTAL_RANKS * EXEC_RANK_DEPTH ))

    done

    # Essentially a ceil calculation
    export REQ_INTERLEAVING_NODES=$(( (SUM_TOTAL_CORES + USED_CORES - 1) / USED_CORES ))

    # For uniform interleaving, all ranks have the same required nodes
    for name in ${EXEC_LABELS}; do
        local REQUIRED_NODES_VAR="${name^^}_REQUIRED_NODES"
        local -n EXEC_REQUIRED_NODES="${REQUIRED_NODES_VAR}"

        EXEC_REQUIRED_NODES=${REQ_INTERLEAVING_NODES}
        # Export the target variable by name, not its value
        export "${REQUIRED_NODES_VAR}"
    done

}