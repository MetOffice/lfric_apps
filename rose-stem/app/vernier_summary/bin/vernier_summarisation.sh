#!/bin/bash
set -euo pipefail

for vout in $(find $TASK_LOG_DIR/ -name vernier-output-*-collated)
do
  gzip -c $vout > $CYLC_TASK_LOG_DIR/$(basename $vout).gz
  summarise-vernier $vout > $CYLC_TASK_LOG_DIR/$(basename $vout).summary
done
