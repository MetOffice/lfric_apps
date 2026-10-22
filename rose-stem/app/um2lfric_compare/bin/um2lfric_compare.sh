#!/usr/bin/env bash

set -eou pipefail

# Exclude Mesh2d variable as this is a randomised integer created by XIOS
nccmp -Fdf --var-diff-count=1 --exclude=Mesh2d $KGO_DIR/checkpoint_um2lfric_000001.nc $SUITE_DIR/checkpoint_um2lfric_000001.nc
