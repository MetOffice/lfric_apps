#!/bin/bash

# *****************************COPYRIGHT*******************************
# (C) Crown copyright 2024 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

# Prepare results location
mkdir -p $TASK_OUTPUT_DIR/results/
# Symbolic link to results directory, so all files can write to /work/results/
ln -sf $TASK_OUTPUT_DIR/results $CYLC_TASK_WORK_DIR/results
# Specific results directory symbolic link exceptions to support
# `lfric_diag.nc`, `lfric_averages.nc` & `lfric_initial` across many tests
# In general output files should target `name="results/fname"` in XIOS xml
ln -sf $TASK_OUTPUT_DIR/results/lfric_diag.nc $CYLC_TASK_WORK_DIR/lfric_diag.nc
ln -sf $TASK_OUTPUT_DIR/results/lfric_averages.nc $CYLC_TASK_WORK_DIR/lfric_averages.nc
ln -sf $TASK_OUTPUT_DIR/results/lfric_initial.nc $CYLC_TASK_WORK_DIR/lfric_initial.nc
if [ $LUSTRE_FILESYSTEM ]; then
  # Set Lustre striping to maximum for results (performance)
  lfs setstripe -c 32 -S 4m -p flash $TASK_OUTPUT_DIR/results/
  lfs setstripe -c 32 -S 4m -p flash $CYLC_SUITE_SHARE_DIR/data/
fi
