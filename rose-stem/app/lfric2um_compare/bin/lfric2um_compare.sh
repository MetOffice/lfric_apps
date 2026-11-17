#!/usr/bin/env bash

set -eou pipefail

mule-cumf --fail-if-differ $KGO_DIR/output.dump $SUITE_DIR/output.dump
