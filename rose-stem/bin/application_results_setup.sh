#!/bin/bash

# *****************************COPYRIGHT*******************************
# (C) Crown copyright 2024 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

# Prepare results location
mkdir -p $TASK_OUTPUT_DIR/results/
if [ $LUSTRE_FILESYSTEM ]; then
  # Set Lustre striping to maximum for results (performance)
  lfs setstripe -c -1 $TASK_OUTPUT_DIR/results/
fi

# Symbolic link for each potential output file,
# from `work` to `results`
# avoiding cp copy commands, as these are very
# storage & wall clock intensive
ln -sf $TASK_OUTPUT_DIR/results/lfric_diagnostics.nc $CYLC_TASK_WORK_DIR/lfric_diagnostics.nc
ln -sf $TASK_OUTPUT_DIR/results/lfric_diag.nc $CYLC_TASK_WORK_DIR/lfric_diag.nc
ln -sf $TASK_OUTPUT_DIR/results/lfric_ver.nc $CYLC_TASK_WORK_DIR/lfric_ver.nc
ln -sf $TASK_OUTPUT_DIR/results/lfric_ver_tp0.nc $CYLC_TASK_WORK_DIR/lfric_ver_tp0.nc
ln -sf $TASK_OUTPUT_DIR/results/lfric_initial.nc $CYLC_TASK_WORK_DIR/lfric_initial.nc
ln -sf $TASK_OUTPUT_DIR/results/lfric_averages.nc $CYLC_TASK_WORK_DIR/lfric_averages.nc
