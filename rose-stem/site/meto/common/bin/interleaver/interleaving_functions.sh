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

    # Safety checks
    # --------------------------------------------------------------------------
    local nodes
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

    gen_nodefile
    gen_binding_file
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

    # Write the rank pattern to the MPMD hostfile
    # -------------------------------------------
    {
        for name in ${EXECUTABLES}; do
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


function gen_binding_file () {
    # Set up variables
    # ----------------
    local name
    local rank
    local core_offset=0

    # Generates the list of rank ids that belong to each MPI rank
    # -----------------------------------------------------------
    {
        for name in ${EXECUTABLES}; do
            # Common work
            local TOTAL_RANKS_VAR="${name^^}_TOTAL_RANKS"
            local RANKS_PER_NODE_VAR="${name^^}_RANKS_PER_NODE"
            local RANK_DEPTH_VAR="${name^^}_RANK_DEPTH"

            local -n EXEC_TOTAL_RANKS="${TOTAL_RANKS_VAR}"
            local -n EXEC_RANKS_PER_NODE="${RANKS_PER_NODE_VAR}"
            local -n EXEC_RANK_DEPTH="${RANK_DEPTH_VAR}"

            if [ "${INTERLEAVING_METHOD}" = "uniform" ]; then
                # Generate the CPU range for each rank.
                for ((rank = 0; rank < EXEC_TOTAL_RANKS; rank++)); do
                    base=$((core_offset + (rank % EXEC_RANKS_PER_NODE) * EXEC_RANK_DEPTH))
                    printf '%d-%d\n' "$base" "$((base + EXEC_RANK_DEPTH - 1))"
                done
                # Start the next executable after this executable's cores.
                (( core_offset += EXEC_RANKS_PER_NODE * EXEC_RANK_DEPTH ))
            fi
        done
    } > "${BINDINGS_FILE}"
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
    #                               executables assigned ranks. e.g.:
    #                    rank no.   0   1   ...   167  |  168    ...    191  (Genoa max)
    #                    exec       a   a   a     a    |  b      b      b
    #                               ^ a start = 0      |  ^ b start = 168

    local MODEL_LAUNCHER_OPTS
    local first_exec="True"

    if [ "${LAUNCHER}" = "mpiexec" ]; then
        # 0. If MEMORY_PROFILE is selected, add the initial parts here
        # ------------------------------------------------------------
        if [ "${MEMORY_PROFILE}" = "True" ]; then
            if [ "${EXTRA_LAUNCHER_OPTS}" != "" ]; then
                echo "LFRIC LAUNCHER INFO:  Overriding EXTRA_LAUNCHER_OPTS in favour of Memory Profiling (MEMORY_PROFILE = "${MEMORY_PROFILE}")" >&2
            fi
            MODEL_LAUNCHER_OPTS+="--line-buffer --label "
        fi

        # 1. Add universal section (hostfile and cpu binding method)
        # ----------------------------------------------------------
        MODEL_LAUNCHER_OPTS+="-hostfile="${HOSTFILE}" --cpu-bind=none "

        # 2. Generate launcher opts for each executable
        # ---------------------------------------------
        for name in ${EXECUTABLES}; do
            local TOTAL_RANKS_VAR="${name^^}_TOTAL_RANKS"
            local RANKS_PER_NODE_VAR="${name^^}_RANKS_PER_NODE"
            local RANK_DEPTH_VAR="${name^^}_RANK_DEPTH"
            local PATH_VAR="${name^^}_PATH"
            local START_RANK_VAR="${name^^}_START_RANK"
            local OMP_NUM_THREADS_VAR="${name^^}_OMP_NUM_THREADS"

            local -n EXEC_TOTAL_RANKS="${TOTAL_RANKS_VAR}"
            local -n EXEC_RANKS_PER_NODE="${RANKS_PER_NODE_VAR}"
            local -n EXEC_RANK_DEPTH="${RANK_DEPTH_VAR}"
            local -n EXEC_PATH="${PATH_VAR}"
            local -n EXEC_START_RANK="${START_RANK_VAR}"
            local -n EXEC_OMP_NUM_THREADS="${OMP_NUM_THREADS_VAR}"

            # Split executable options
            if [ "${first_exec}" = "False" ]; then
                MODEL_LAUNCHER_OPTS+=" : "
            fi

            # Add opts:
            MODEL_LAUNCHER_OPTS+="--depth=$EXEC_RANK_DEPTH "
            MODEL_LAUNCHER_OPTS+="-n "${EXEC_TOTAL_RANKS}" "
            # Rank binding opts
            MODEL_LAUNCHER_OPTS+="bash $RANK_BINDER_SCRIPT $EXEC_START_RANK $EXEC_TOTAL_RANKS $EXEC_PPN $EXEC_CPU_START $EXEC_RANK_DEPTH "
            MODEL_LAUNCHER_OPTS+="env OMP_NUM_THREADS=$EXEC_OMP_NUM_THREADS "

            if [ "${MEMORY_PROFILE}" = "True" ]; then
                MODEL_LAUNCHER_OPTS+="/usr/bin/time -f max_mem_${name}_%Mkb "
            else
                if [ "${EXTRA_LAUNCHER_OPTS}" != "" ]; then
                    MODEL_LAUNCHER_OPTS+="${EXTRA_LAUNCHER_OPTS} "
                fi
            fi

            # Add Executable
            MODEL_LAUNCHER_OPTS+="${EXEC_PATH} "

            # Add executable args
            if [ "${EXEC_ARGS}" != "" ]; then
                MODEL_LAUNCHER_OPTS+="${EXEC_ARGS} "
            fi

            # Will be false from here until the end
            # means that " : " will be added before next set of exec stuff
            first_exec="False"
        done
    else
        echo "LFRIC LAUNCHER ERROR: Interleaving currently is supported for mpiexec launcher only." >&2
    fi
}
