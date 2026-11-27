#!/usr/bin/env bash
# (c) Crown copyright Met Office. All rights reserved.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Wrapper to run relevant interleaving method
function utils::interleave () {

    if [ "${INTERLEAVING_METHOD}" = "uniform" ]; then
        source "${SCRIPT_DIR}/uniform_interleaving.sh"
    fi

    validate_interleaving_settings
    generate_interleaving_vars
    gen_nodefile
    gen_exec_command

}

# This function generates a list of node ids that each executable will be
# placed on.
function gen_nodefile () {
    # Set up variables
    # ----------------
    local nodes
    local name
    local rank

    # Safety checks
    # --------------------------------------------------------------------------
    # Read the unique, non-empty nodes allocated by PBS into an array
    mapfile -t nodes < <(awk 'NF && !seen[$0]++' "$PBS_NODEFILE")

    # Check PBS nodefile is not empty
    if (( ${#nodes[@]} == 0 )); then
        echo "PBS_NODEFILE contains no nodes" >&2
        exit 1
    fi

    # Check that the number of unique nodes is equal to the number of required
    # nodes
    if (( ${#nodes[@]} != REQ_INTERLEAVING_NODES )); then
        echo "LFRIC LAUNCHER ERROR: Interleaving config requires $REQ_INTERLEAVING_NODES nodes, but PBS_NODEFILE contains ${#nodes[@]}" >&2
        exit 1
    fi
    # --------------------------------------------------------------------------

    # Write the rank pattern to the MPMD hostfile
    # -------------------------------------------
    {
        for name in ${EXEC_LABELS}; do
            local TOTAL_RANKS_VAR="${name^^}_TOTAL_RANKS"
            local RANKS_PER_NODE_VAR="${name^^}_RANKS_PER_NODE"

            local -n EXEC_TOTAL_RANKS="${TOTAL_RANKS_VAR}"
            local -n EXEC_RANKS_PER_NODE="${RANKS_PER_NODE_VAR}"

            if [ "${INTERLEAVING_METHOD}" = "uniform" ]; then
                # Each executable occupies a consecutive block of hostfile entries,
                # with RANKS_PER_NODE ranks assigned to each node.
                for ((rank = 0; rank < EXEC_TOTAL_RANKS; rank++)); do
                    printf '%s\n' "${nodes[rank / EXEC_RANKS_PER_NODE]}"
                done
            fi
        done
    } > "${HOSTFILE}"


}

function gen_exec_command () {
    # This function requires the following variables to be set by the selected
    # interleaving method's functions:
    # <EXEC>_START_RANK             The rank number in the global list of ranks (PMI)
    #                               that marks the start of the executable's
    #                               assigned ranks. e.g.:
    #                    rank no.   0   1   ...   1999 |  2000   2001   2002
    #                    exec       a   a   a     a    |  b      b      b
    #                               ^ a start = 0      |  ^ b start = 2000
    #
    # <EXEC>_CPU_START              The rank number in the local (node-level)
    #                               list of ranks that marks the start of the
    #                               EXEC_LABELS assigned ranks. e.g.:
    #                    rank no.   0   1   ...   167  |  168    ...    191  (Genoa max)
    #                    exec       a   a   a     a    |  b      b      b
    #                               ^ a start = 0      |  ^ b start = 168

    EXEC_COMMAND=""
    local RANK_BINDER_SCRIPT="${SCRIPT_DIR}/rank_binder.sh"
    local first_exec="True"

    if [ "${RUN_METHOD}" = "mpiexec" ]; then
        # 0. Start the exec command
        # -------------------------
        EXEC_COMMAND+="mpiexec "

        # 0.1 If MEMORY_PROFILE is selected, add the initial parts here
        # ------------------------------------------------------------
        if [ "${MEMORY_PROFILE}" = "True" ]; then
            if [ "${EXTRA_LAUNCHER_OPTS}" != "" ]; then
                echo "LFRIC LAUNCHER INFO:  Overriding EXTRA_LAUNCHER_OPTS in favour of Memory Profiling (MEMORY_PROFILE = "${MEMORY_PROFILE}")" >&2
            fi
            EXEC_COMMAND+="--line-buffer --label "
        fi

        # 1. Add universal section (hostfile and cpu binding method)
        # ----------------------------------------------------------
        EXEC_COMMAND+="-hostfile="${HOSTFILE}" --cpu-bind=none "

        # 2. Generate launcher opts for each executable
        # ---------------------------------------------
        for name in ${EXEC_LABELS}; do
            local TOTAL_RANKS_VAR="${name^^}_TOTAL_RANKS"
            local RANKS_PER_NODE_VAR="${name^^}_RANKS_PER_NODE"
            local RANK_DEPTH_VAR="${name^^}_RANK_DEPTH"
            local PATH_VAR="${name^^}_PATH"
            local OMP_NUM_THREADS_VAR="${name^^}_OMP_NUM_THREADS"
            local PPN_VAR="${name^^}_PPN"
            local ARGS_VAR="${name^^}_ARGS"
            local START_RANK_VAR="${name^^}_START_RANK"
            local CPU_START_VAR="${name^^}_CPU_START"

            local -n EXEC_TOTAL_RANKS="${TOTAL_RANKS_VAR}"
            local -n EXEC_RANKS_PER_NODE="${RANKS_PER_NODE_VAR}"
            local -n EXEC_RANK_DEPTH="${RANK_DEPTH_VAR}"
            local -n EXEC_PATH="${PATH_VAR}"
            local -n EXEC_OMP_NUM_THREADS="${OMP_NUM_THREADS_VAR}"
            local -n EXEC_PPN="${PPN_VAR}"
            local -n EXEC_ARGS="${ARGS_VAR}"
            local -n EXEC_START_RANK="${START_RANK_VAR}"
            local -n EXEC_CPU_START="${CPU_START_VAR}"

            # Split executable options
            if [ "${first_exec}" = "False" ]; then
                EXEC_COMMAND+=" : "
            fi

            # Add opts:
            EXEC_COMMAND+="--depth=$EXEC_RANK_DEPTH "
            EXEC_COMMAND+="-n "${EXEC_TOTAL_RANKS}" "
            # Rank binding opts
            EXEC_COMMAND+="bash $RANK_BINDER_SCRIPT $EXEC_START_RANK $EXEC_TOTAL_RANKS $EXEC_PPN $EXEC_CPU_START $EXEC_RANK_DEPTH "
            EXEC_COMMAND+="env OMP_NUM_THREADS=$EXEC_OMP_NUM_THREADS "

            if [ "${MEMORY_PROFILE}" = "True" ]; then
                EXEC_COMMAND+="/usr/bin/time -f max_mem_${name}_%Mkb "
            else
                if [ "${EXTRA_LAUNCHER_OPTS}" != "" ]; then
                    EXEC_COMMAND+="${EXTRA_LAUNCHER_OPTS} "
                fi
            fi

            # Add Executable
            EXEC_COMMAND+="${EXEC_PATH}"

            if [ "${EXEC_ARGS}" != "" ]; then
                EXEC_COMMAND+=" ${EXEC_ARGS}"
            fi

            # Will be false from here until the end
            # means that " : " will be added before next set of exec stuff
            first_exec="False"
        done
    else
        echo "LFRIC LAUNCHER ERROR: Interleaving currently is supported for mpiexec launcher only." >&2
    fi

    export EXEC_COMMAND
}
