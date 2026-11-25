
# Helper functions
# Check if item exists in a space-separated list
function exists_in_list () {
    local ITEM=$1
    local LIST=$2

    for x in $LIST; do
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
    : "${XIOS_SERVER_RANKS:=0}"
    : "${xios_nodes:=0}"
    : "${mpi_parts_xios:=0}"
    : "${MEMORY_PROFILE:="False"}"
    : "${EXTRA_LAUNCHER_OPTS:=""}"
    : "${TEST_LAUNCH_EXE_EXEC:=""}"
    : "${TEST_MODE:=""}"
    : "${NAMELIST_FILE="${EXTRA_NAMELIST_FILES:=configuration.nml}"}"
    : "${INTERLEAVE_EXECS:="False"}"

    # Set default values for UNSET env vars - preserves empty
    : "${XIOS_SERVER_EXEC_PATH="xios_server.exe"}"
    : "${PAT_EXE_EXTEN=""}"
    : "${LAUNCHER=""}"
    : "${BIN_DIR="."}"  # default to current directory

    # Depth value - default to OMP threads unless RANKS_DEPTH_PAD is larger
    : "${DEPTH_VAL:=${OMP_NUM_THREADS}}"
    if (( "${RANKS_DEPTH_PAD}" > "${OMP_NUM_THREADS}" )); then
        export DEPTH_VAL=${RANKS_DEPTH_PAD}
    fi

    # Test Mode
    # if TEST_LAUNCH_EXE_EXEC is set, then it's value shall be the expected
    # value to test the behaviour of `launch-exe` and the script will error
    # if the value does not match.
    if [ "${TEST_LAUNCH_EXE_EXEC}" != "" ]; then
        export TEST_MODE=0
    fi
}

# Validate environment variables
function utils::validate_env () {

    # Values must match a set of pre-determined options
    exists_in_list "${TARGET_PLATFORM}" "${VALID_PLATFORMS}"
    exists_in_list "${RUN_METHOD}" "${VALID_RUN_METHODS}"

    # Values must be an integer but we don't care about their size
    check_integer CORES_PER_NODE
    check_integer CORES_PER_NODE_OVERRIDE

    if [ "${INTERLEAVE_EXECS}" = "False" ]; then
        # None of these are required by interleaving
        exists_in_list "${XIOS_SERVER_MODE}" "${BOOL_VALUES}"
        check_integer XIOS_SERVER_RANKS
        check_integer xios_nodes
        check_integer mpi_parts_xios
        check_integer TOTAL_RANKS
        check_integer RANKS_DEPTH_PAD
        check_integer OMP_NUM_THREADS
        # This line will cause a fault if EXEC_NAME is not set.
        : "${EXEC_NAME:?LFRIC LAUNCHER ERROR: EXEC_NAME must name an executable to run.}"
    else
        utils::validate_interleaving_vars
    fi
}

function utils::validate_interleaving_vars () {
    # Interleaving requires the following variables:
    # EXEC_LABELS               A space-separated list of executable names in MPMD order. Names are case-insensitive.

    # For each name listed in EXEC_LABELS, the following are required (all uppercase):
    # <EXEC>_RANKS_PER_NODE     The number of MPI ranks to give the executable.
    # <EXEC>_RANK_DEPTH         The number of CPUs to allocate to each MPI rank.
    # <EXEC>_OMP_NUM_THREADS    The number of OpenMP threads the executable will use.
    # <EXEC>_PATH               The path to use for the executable.
    # <EXEC>_TOTAL_RANKS        The total number of ranks the executable will use.

    # 0. Validate that we are on a valid platform for interleaving
    # ------------------------------------------------------------
    if [ "${TARGET_PLATFORM}" != "meto-ex1a" ]; then
        echo "LFRIC LAUNCHER ERROR: Interleaving is only valid on the EX platform" >&2
        exit 1
    fi

    # 1. Validate that we are using a valid method
    # --------------------------------------------
    # Set default interleaving method to "uniform"
    : "${INTERLEAVING_METHOD:="uniform"}"

    # Validate that the method is allowed
    exists_in_list "${INTERLEAVING_METHOD}" "${VALID_INTERLEAVING_METHODS}"

    if [ "${INTERLEAVING_METHOD}" != "uniform" ]; then
        echo "LFRIC LAUNCHER ERROR: Only uniform interleaving currently supported" >&2
        exit 1
    fi

    # 2. Validate the required variables exist with values
    # ----------------------------------------------------
    # Assert that we have executables to interleave
    : "${EXEC_LABELS:=?LFRIC LAUNCHER ERROR: No executables defined in EXEC_LABELS but INTERLEAVE_EXECS is True}"

    # Assert that the required env vars for interleaving are set
    for name in ${EXEC_LABELS}; do
        check_integer "${name^^}_RANKS_PER_NODE"
        check_integer "${name^^}_RANK_DEPTH"
        check_integer "${name^^}_OMP_NUM_THREADS"
        check_integer "${name^^}_TOTAL_RANKS"
        # Check if path is empty or unset
        local PATH_VAR="${name^^}_PATH"
        : "${!PATH_VAR:=?LFRIC LAUNCHER ERROR: No path set for ${name} executable (use "${"${name^^}_PATH"}")}"
        # Assure that exec's args are not unset - can be empty
        local ARGS_VAR="${name^^}_ARGS"
        : "${!ARGS_VAR=""}"
    done

    # Assert that an output directory exists
    # (default to CYLC_TASK_WORK_DIR, which can default to BIN_DIR)
    : "${INTERLEAVER_FILES:="${CYLC_TASK_WORK_DIR:="${BIN_DIR}"}"}"

    # And the files that interleaving depends on (defaults for if empty or unset)
    : "${HOSTFILE:="${INTERLEAVER_FILES}/hostfile.txt"}"
    : "${BINDINGS_FILE:="${INTERLEAVER_FILES}/cpu_bind.txt"}"
}

# This function exists to set variables that are derived from required
# settings but are not themselves required to be set by the user.
function utils::set_derived_vars () {

    # Set launcher (default is "")
    if [ "${RUN_METHOD}" = "mpiexec" ] ; then
        export LAUNCHER=${RUN_METHOD}
    fi

    # Set executable path (if not already set or if empty)
    : "${EXEC_PATH:="${BIN_DIR}/${EXEC_NAME}${PAT_EXE_EXTEN}"}"
}

# Determine cores per-node requirements
function utils::get_cores_per_node () {
    # Set MPI ranks / PEs per Node from provided value
    # ================================================
    # Output -> `MODEL_PPN` (used for ppn)
    # =============================================
    # Maximum number of MPI ranks that can fit on one node
    local MAX_PPN=$(( CORES_PER_NODE / DEPTH_VAL ))

    # Keep PPN even, rounding down to avoid exceeding the node
    if (( MAX_PPN % 2 != 0 )); then
        (( MAX_PPN-- ))
    fi

    # Default PPN
    if (( TOTAL_RANKS <= MAX_PPN )); then
        # Entire job fits on one node
        export MODEL_PPN="${TOTAL_RANKS}"

        echo "LFRIC LAUNCHER INFO: Ranks per node has been set equal to the required ranks for this job (${MODEL_PPN})" >&2
    else
        # Job requires multiple nodes
        export MODEL_PPN="${MAX_PPN}"
    fi

    # Override PPN if supplied and valid
    if (( CORES_PER_NODE_OVERRIDE > 0 )); then
        if (( CORES_PER_NODE_OVERRIDE * DEPTH_VAL <= CORES_PER_NODE )); then
            export MODEL_PPN="${CORES_PER_NODE_OVERRIDE}"

            echo "LFRIC LAUNCHER INFO: Ranks per node set to supplied override value (${MODEL_PPN})" >&2
        else
            echo "LFRIC LAUNCHER INFO: Ranks per node defaulted to max possible ranks per node (${MODEL_PPN})" >&2
        fi
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
            export MODEL_LAUNCHER_OPTS="--cpu-bind=depth --np ${TOTAL_RANKS} --depth ${DEPTH_VAL}"

            # Legacy support:
            # Job spans multiple nodes, or PPN was explicitly overridden.
            if (( TOTAL_RANKS * DEPTH_VAL > MODEL_PPN || CORES_PER_NODE_OVERRIDE > 0 )); then
                # XIOS with no dedicated XIOS nodes is the exception.
                if [[ "${XIOS_SERVER_MODE}" != "True" || "${xios_nodes}" != "0" ]]; then
                    export MODEL_LAUNCHER_OPTS+=" --ppn ${MODEL_PPN}"
                fi
            fi
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

    if [[ "${XIOS_SERVER_MODE}" = "False" || "${RUN_METHOD}" = "executable" ]]; then
        export XIOS_LAUNCHER_OPTS=""
        export XIOS_SERVER_EXEC_PATH=""
        return
    fi

    if [ "${RUN_METHOD}" = "mpiexec" ]; then

        if [ "${TARGET_PLATFORM}" = "meto-ex1a" ]; then
            if [ "${xios_nodes}" = "0" ]; then
                export XIOS_LAUNCHER_OPTS=" : --cpu-bind=depth --np ${XIOS_SERVER_RANKS} "
                return
            else
                # Define the XIOS launcher options depending on run method and XIOS mode
                local XIOS_CORES_PER_NODE
                XIOS_CORES_PER_NODE="$(( (mpi_parts_xios / xios_nodes) + (mpi_parts_xios % xios_nodes > 0 ) ))"

                # Calculate the depth - sets XIOS_DEPTH variable
                get_xios_depth

                export XIOS_LAUNCHER_OPTS=" : --cpu-bind=depth --np ${XIOS_SERVER_RANKS}${XIOS_DEPTH} --ppn ${XIOS_CORES_PER_NODE} "
                return
            fi

        elif [ "${TARGET_PLATFORM}" = "meto-azspice" ]; then
            export XIOS_LAUNCHER_OPTS= ": -n ${XIOS_SERVER_RANKS} "
            return
        fi
    fi
}

function utils::build_exec_command () {
    # Combines all the component parts for model, xios, and memory profiling to
    # create a single executable command
    if [ "${MEMORY_PROFILE}" = "False" ] ; then
        if [ "${EXTRA_LAUNCHER_OPTS}" != "" ]; then
            MODEL_LAUNCHER_OPTS="${MODEL_LAUNCHER_OPTS} ${EXTRA_LAUNCHER_OPTS}"
        fi
        export EXEC_COMMAND="${LAUNCHER} ${MODEL_LAUNCHER_OPTS} ${EXEC_PATH} ${NAMELIST_FILE} ${XIOS_LAUNCHER_OPTS}${XIOS_SERVER_EXEC_PATH}"
        return 0
    else
        # Set the XIOS memory profiling command
        local XIOS_MEM_TIME=""
        if [ "$XIOS_LAUNCHER_OPTS" != "" ] ; then
            XIOS_MEM_TIME="/usr/bin/time -f max_mem_xios_server_%Mkb "
        fi

        # Warn users that their provided extras are being ignored
        if [ "${EXTRA_LAUNCHER_OPTS}" != "" ]; then
            echo "LFRIC LAUNCHER INFO: Overriding EXTRA_LAUNCHER_OPTS in favour of Memory Profiling (MEMORY_PROFILE = "${MEMORY_PROFILE}")" >&2
        fi

        # Build the command
        if [ "${TARGET_PLATFORM}" = "meto-ex1a" ] ; then
            export EXEC_COMMAND="${LAUNCHER} --line-buffer --label ${MODEL_LAUNCHER_OPTS} /usr/bin/time -f max_mem_${EXEC_NAME}_%Mkb ${EXEC_PATH} ${NAMELIST_FILE} ${XIOS_LAUNCHER_OPTS}${XIOS_MEM_TIME}${XIOS_SERVER_EXEC_PATH}"
        else
            export EXEC_COMMAND="${LAUNCHER} -l ${MODEL_LAUNCHER_OPTS} /usr/bin/time -f max_mem_${EXEC_NAME}_%Mkb ${EXEC_PATH} ${NAMELIST_FILE}  ${XIOS_LAUNCHER_OPTS}${XIOS_MEM_TIME}${XIOS_SERVER_EXEC_PATH}"
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

