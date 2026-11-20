
# Helper functions
# Check if item exists in a space-separated list
function exists_in_list () {
    local ITEM=$1
    local LIST=$2

    for x in $LIST_ITEMS; do
        if [ "$x" = "$ITEM" ]; then
            return 0
        fi
    done

    echo "LFRIC LAUNCHER ERROR: ${ITEM}" not found in a supplied list of valid options: "${LIST}" >&2
    exit 1
}

# Check if the variable is a positive integer
function check_integer () {
    # Function takes the variable name unexpanded (without dollar sign) for
    # better error reporting - can give the variable name with error
    local VAR_NAME="$1"
    local VAR="${!VAR_NAME-}"

    if [[ ! "$VAR" =~ ^[0-9]+$ ]]; then
        echo "LFRIC LAUNCHER ERROR: $VAR_NAME must be a non-negative integer, got: '$VAR'" >&2
        exit 1
    fi
}

# Set default values for required variables
function utils::set_default_env () {
    # Set default values for UNSET OR EMPTY env vars (overrides "")
    : "${TOTAL_RANKS:=1}"
    : "${OMP_NUM_THREADS:=1}"
    : "${CORES_PER_NODE_OVERRIDE:=0}"
    : "${RANKS_DEPTH_PAD:=1}"
    : "${XIOS_SERVER_MODE:="False"}"
    : "${MEMORY_PROFILE:="False"}"
    : "${EXTRA_LAUNCH_OPTS:=""}"

    # Set default values for UNSET env vars - preserves empty
    : "${XIOS_SERVER_EXEC_PATH="xios_server.exe"}"
    : "${PAT_EXE_EXTEN=""}"

    # Depth value - default to OMP threads unless RANKS_DEPTH_PAD is larger
    : "${DEPTH_VAL:=${OMP_NUM_THREADS}}"
    if (( "${RANKS_DEPTH_PAD}" > "${OMP_NUM_THREADS}" )); then
        export DEPTH_VAL=${RANKS_DEPTH_PAD}
    fi

    # Set up Executable information
    : "${EXEC_PATH:="${BIN_DIR}/${EXEC_NAME}${PAT_EXE_EXTEN}"}"
    : "${NAMELIST_FILE="${EXTRA_NAMELIST_FILES:=configuration.nml}"}"

}

# Validate environment variables
function utils::validate_env () {

    # Values must match a set of pre-determined options
    exists_in_list "${TARGET_PLATFORM}" "${VALID_PLATFORMS}"
    exists_in_list "${RUN_METHOD}" "${VALID_RUN_METHODS}"
    exists_in_list "${XIOS_SERVER_MODE}" "${BOOL_VALUES}"

    # Values must be an integer but we don't care about their size
    check_integer CORES_PER_NODE
    check_integer OMP_NUM_THREADS
    check_integer XIOS_SERVER_RANKS
    check_integer xios_nodes
    check_integer mpi_parts_xios
    check_integer RANKS_DEPTH_PAD
    check_integer TOTAL_RANKS
    check_integer CORES_PER_NODE_OVERRIDE

}

# Determine cores per-node requirements
# This is used when the CORES_PER_NODE override is set.
function utils::get_cores_per_node () {
    # Set MPI ranks / PEs per Node from provided value
    # ================================================
    # Output -> `SET_CORES_PER_NODE` (used for ppn)
    # =============================================
    # Set a default value based on cores per node and OpenMP
    local CORES_PER_NODE_OVER_DEPTH
    CORES_PER_NODE_OVER_DEPTH="$(( CORES_PER_NODE / DEPTH_VAL ))"

    # Safety check: if CORES_PER_NODE does not cleanly divide by the supplied
    # depth then we can't use the result as a ppn value
    check_integer CORES_PER_NODE_OVER_DEPTH

    # Override provided
    # -----------------
    if (( "${CORES_PER_NODE_OVERRIDE}" > "0" )) ; then
        # Multiply by depth value for ppn
        local REQUESTED_PPN="$(( CORES_PER_NODE_OVERRIDE * DEPTH_VAL ))"

        # Check if this value is less than available cores per node
        if (( "${REQUESTED_PPN}" <= "${CORES_PER_NODE}" )) ; then
            export SET_CORES_PER_NODE="${CORES_PER_NODE_OVERRIDE}"
            echo "LFRIC LAUNCHER INFO: Ranks per node has been set to the supplied override value ("${SET_CORES_PER_NODE}")" >&2
        else
            # If nothing appropriate can be set, fall back on the default
            export SET_CORES_PER_NODE="${CORES_PER_NODE_OVER_DEPTH}"
            echo "LFRIC LAUNCHER INFO: Ranks per node has defaulted to standard value of maximum ranks per node ("${SET_CORES_PER_NODE}")" >&2
        fi

    # No override provided
    # --------------------
    else
        export SET_CORES_PER_NODE="${CORES_PER_NODE_OVER_DEPTH}"
    fi

    # Required ranks for job smaller than requested
    # ---------------------------------------------
    # Set the cores per node to the same as the total ranks provided if too few
    if (( "$((TOTAL_RANKS * DEPTH_VAL))" < "${SET_CORES_PER_NODE}" )) ; then
        export SET_CORES_PER_NODE="${TOTAL_RANKS}"
        echo "LFRIC LAUNCHER INFO: Ranks per node has been set equal to the required ranks for this job ("${SET_CORES_PER_NODE}")" >&2
    fi

}

# Launcher options - model
function utils::set_model_options () {
    # Define the model launcher options depending on the run method
    if [ "${RUN_METHOD}" = "executable" ]; then
        export MODEL_LAUNCHER_OPTS=""
        return 0

    elif [ "${RUN_METHOD}" = "mpiexec" ]; then
        if [ "${TARGET_PLATFORM}" = "meto-azspice" ]; then
            export MODEL_LAUNCHER_OPTS="-n ${TOTAL_RANKS}"
            return 0
        elif [ "${TARGET_PLATFORM}" = "meto-ex1a" ]; then
            export MODEL_LAUNCHER_OPTS="--cpu-bind=depth --np ${TOTAL_RANKS} --depth ${DEPTH_VAL} --ppn ${SET_CORES_PER_NODE}"
            return 0
        fi
    fi
}

function get_xios_depth () {
    # ensure ceiling rounding for division
    local XIOS_FRACTIONAL_DEPTH
    XIOS_FRACTIONAL_DEPTH="$(( CORES_PER_NODE / XIOS_CORES_PER_NODE ))"

    # Milan = 128, Genoa = 192
    if (( "${XIOS_FRACTIONAL_DEPTH}" >= "$(( CORES_PER_NODE ))" )) ; then
        export XIOS_DEPTH=" --depth $(( CORES_PER_NODE ))"

    # Milan = 64, Genoa = 96
    elif (( "${XIOS_FRACTIONAL_DEPTH}" >= "$(( CORES_PER_NODE / 2 ))" )) ; then
        export XIOS_DEPTH=" --depth $(( CORES_PER_NODE / 2 ))"

    # Milan = 32, Genoa = 48
    elif (( "${XIOS_FRACTIONAL_DEPTH}" >= "$(( CORES_PER_NODE / 4 ))" )) ; then
        export XIOS_DEPTH=" --depth $(( CORES_PER_NODE / 4 ))"

    # Milan = 16, Genoa = 24
    elif (( "${XIOS_FRACTIONAL_DEPTH}" >= "$(( CORES_PER_NODE / 8 ))" )) ; then
        export XIOS_DEPTH=" --depth $(( CORES_PER_NODE / 8 ))"

    # Milan = 8, Genoa = 12
    elif (( "${XIOS_FRACTIONAL_DEPTH}" >= "$(( CORES_PER_NODE / 16 ))" )) ; then
        export XIOS_DEPTH=" --depth $(( CORES_PER_NODE / 16 ))"

    # Milan = 8, Genoa = 8 - To better fit onto Genoa hardware this is needed
    # Duplication here for Milan but no negative effects
    elif (( "${XIOS_FRACTIONAL_DEPTH}" >= 8 )) ; then
        export XIOS_DEPTH=" --depth 8"
    elif (( "${XIOS_FRACTIONAL_DEPTH}" >= 4 )) ; then
        export XIOS_DEPTH=" --depth 4"
    elif (( "${XIOS_FRACTIONAL_DEPTH}" >= 2 )) ; then
        export XIOS_DEPTH=" --depth 2"
    else
        export XIOS_DEPTH=""
    fi
}

# Launcher options - XIOS
function utils::set_xios_options () {
    # Define the XIOS launcher options depending on run method and XIOS mode
    local XIOS_CORES_PER_NODE
    XIOS_CORES_PER_NODE="$(( (mpi_parts_xios / xios_nodes) + (mpi_parts_xios % xios_nodes > 0 ) ))"

    if [ "${RUN_METHOD}" = "executable" ]; then
        export XIOS_LAUNCHER_OPTS=""
        export XIOS_SERVER_EXEC_PATH=""
        return 0
    elif [ "${RUN_METHOD}" = "mpiexec" ]; then
        if [ "${TARGET_PLATFORM}" = "meto-ex1a" ]; then
            if [ "${XIOS_SERVER_MODE}" = "False" ]; then
                export XIOS_LAUNCHER_OPTS=""
                export XIOS_SERVER_EXEC_PATH=""
                return 0
            else
                if [ "${xios_nodes}" = "0" ]; then
                    export XIOS_LAUNCHER_OPTS=" : --cpu-bind=depth --np ${XIOS_SERVER_RANKS} "
                    return 0
                else
                    # Calculate the depth - sets XIOS_DEPTH variable
                    get_xios_depth

                    export XIOS_LAUNCHER_OPTS=" : --cpu-bind=depth --np ${XIOS_SERVER_RANKS}${XIOS_DEPTH} --ppn ${XIOS_CORES_PER_NODE} "
                    return 0
                fi
            fi
        elif [ "${TARGET_PLATFORM}" = "meto-azspice" ]; then
            export XIOS_LAUNCHER_OPTS=" : -n ${XIOS_SERVER_RANKS} "
            return 0
        fi
    fi
    return 1
}

function utils::build_exec_command () {
    # Combines all the component parts for model, xios, and memory profiling to
    # create a single executable command
    if [ "${MEMORY_PROFILE}" = "False" ] ; then
        export EXEC_COMMAND="${LAUNCHER} ${LAUNCHER_OPTS} ${EXEC_PATH} ${NAMELIST_FILE} ${XIOS_OPTS}${XIOS_SERVER_EXEC_PATH}"
        return 0
    else
        local XIOS_MEM_TIME
        XIOS_MEM_TIME=""
        if [ "$XIOS_OPTS" != "" ] ; then
            XIOS_MEM_TIME="/usr/bin/time -f max_mem_xios_server_%Mkb "
        fi
        if [ "${TARGET_PLATFORM}" = "meto-ex1a" ] ; then
            export EXEC_COMMAND="${LAUNCHER} --line-buffer --label ${LAUNCHER_OPTS} /usr/bin/time -f max_mem_${EXEC_NAME}_%Mkb ${EXEC_PATH} ${NAMELIST_FILE} ${XIOS_OPTS}${XIOS_MEM_TIME}${XIOS_SERVER_EXEC_PATH}"
        else
            export EXEC_COMMAND="${LAUNCHER} -l ${LAUNCHER_OPTS} /usr/bin/time -f max_mem_${EXEC_NAME}_%Mkb ${EXEC_PATH} ${NAMELIST_FILE} ${XIOS_OPTS}${XIOS_MEM_TIME}${XIOS_SERVER_EXEC_PATH}"
        fi
    fi

}

function utils::run_memory_profiler () {

    if [[ "${MEMORY_PROFILE}" != "True" || "${TARGET_PLATFORM}" != "meto-ex1a" ]] ; then
        return
    fi

    # report on Memory usage in job.out
    local number_of_nodes
    number_of_nodes=$(sort -u "$PBS_NODEFILE" | wc -l)
    cat > Scr_Memory.job << 'EOF'
#!/bin/bash
echo "$(uname -n) $PALS_NODEID $(( $(cat /sys/fs/cgroup/memory/pbspro.service/jobid/$PBS_JOBID/memory.max_usage_in_bytes) / 1048576 ))"
EOF
    chmod u+x Scr_Memory.job

    echo "Maximum memory usage per node (MiB)"

    mpiexec -n "$number_of_nodes" --ppn 1 ./Scr_Memory.job |
        sort -r -n -k3 |
        column -t -o" " |
        sed 's/^/  /'
}

