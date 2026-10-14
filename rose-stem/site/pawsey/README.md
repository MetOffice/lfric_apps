# Setonix

This document describes the Setonix-specific changes, dependencies, workarounds, and instructions required to build and run LFRic using Rose-stem.

## Dependencies

For **LFRic version 3.2**, update the `lfric_core` entry in `dependencies.yaml` to use the following revision. This revision includes changes required for compatibility with **CCE 18.0.1** and **Cray MPICH 9.0.1**.

```yaml
lfric_core:
  source: git@github.com:ddeeptimahanti/lfric_core.git
  ref: 666770aaf375c8888858d37095c08d7714ae9904
```

## Cylc configuration on Setonix

Before running Rose-stem, configure Cylc by creating the following files in your home directory.

### 1. Rose configuration: `~/.cylc/rose.conf`

This file selects the Setonix Rose-stem site configuration and automatically sets the `SITE` option.

```ini
site=setonix[rose-stem]
automatic-options=SITE=setonix
```

### 2. Cylc global configuration: `~/.cylc/flow/8/global.cylc`

This file configures the Setonix Slurm platform and specifies the directory used for workflow installation and execution.

```ini
[install]
    [[symlink dirs]]
        [[[localhost]]]
            run = /scratch/pawsey0835/ddeeptimahanti

[platforms]
    [[setonix,setonix-slurm]]
        hosts = localhost
        job runner = slurm
        install target = localhost
```

**Important:** Replace `/scratch/pawsey0835/ddeeptimahanti` with your own scratch directory. The configuration above uses a user-specific path and should not be copied unchanged by other users.

Create the configuration directory and files if they do not already exist:

```bash
mkdir -p ~/.cylc/flow/8
```

## Running Rose-stem with Cylc

**Cylc 8.5.4** is available as a module on Setonix. The LFRic Apps Rose-stem suite requires Cylc 8.

From the `lfric_apps` working copy, load the appropriate Cylc module and run the following command to launch the workflow:

```bash
module use /software/projects/pawsey0835/ddeeptimahanti/modules/
module load cylc/8.5.4
cylc vip \
    -S SITE="'pawsey'" \
    -z group=lfric_atm_gnu_prod,lfric_atm_cce_prod,lfric_atm_gnu_full,lfric_atm_cce_full,lfric_atm_gnu_fast,lfric_atm_cce_fast \
    -S USE_MIRRORS=false \
    -n my_rose_stem_suite \
    --set 'CYLC_WORKFLOW_SRC_DIR="'"$PWD/rose-stem"'"' \
    ./rose-stem \
    --verbose
```

The command uses the following options:

- `-S SITE="'pawsey'"`: Selects the Pawsey site configuration.
- `-z group=...`: Selects the LFRic atmosphere test groups for GNU and CCE compilers across production, full-debug, and fast-debug configurations.
- `-S USE_MIRRORS=false`: Disables the use of mirrors.
- `-n my_rose_stem_suite`: Sets the workflow name. Replace this with a suitable name for your run.
- `--set CYLC_WORKFLOW_SRC_DIR=...`: Sets the workflow source directory to the `rose-stem` directory in the current working copy.
- `--verbose`: Enables verbose output.

Ensure that the required dependencies and site configuration are available before launching the workflow.

## CCE 18.0.1 linking error

When compiling with **CCE 18.0.1**, the `ftn` compiler may report the following error during the linking stage:

```text
ftn-2116 ftn: INTERNAL
"driver" was terminated due to receipt of signal 013: Segmentation fault.
```

This error can cause Cylc to terminate the workflow.

### Workaround

A conditional workaround has been added to:

```text
rose-stem/templates/runtime/generate_runtime_build.cylc
```

This allows the workflow to continue when the specific linking error is encountered.

**Important:** This workaround allows the build task to report success despite the compiler error. Verify that the executable was generated correctly and passes the relevant tests before using it for production runs.

## Variables used before definition with CCE 18.0.1

CCE 18.0.1 reports errors in `pc2_bl_forced_cu.F90` when the following variables are used before being defined:

- `dcfl`
- `cf_forced`
- `cf_base`
- `dqcl`
- `qcl_tol`
- `qcl_forced`

The affected source file is:

```text
lfric_apps/science/physics_schemes/source/large_scale_cloud/pc2_bl_forced_cu.F90
```

### Workaround

Initialise the affected variables to `0.0` before they are used:

```fortran
dcfl       = 0.0
cf_forced  = 0.0
cf_base    = 0.0
dqcl       = 0.0
qcl_tol    = 0.0
qcl_forced = 0.0
```

These initialisations address the reported use-before-definition errors in CCE 18.0.1.
