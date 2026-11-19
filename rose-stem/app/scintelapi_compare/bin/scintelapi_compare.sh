#!/usr/bin/env bash

set -eou pipefail

# Exclude Mesh2d variable as this is a randomised integer created by XIOS
nccmp -Fdmf --var-diff-count=1 --exclude=Mesh2d $KGO_DIR/target_dump.nc $SUITE_DIR/target_dump.nc
