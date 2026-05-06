#!/usr/bin/env bash

set -eou pipefail

python3 rose_stem_extract_source.py
cp $SOURCE_ROOT/lfric_apps/dependencies.yaml $CYLC_WORKFLOW_RUN_DIR
cp $SOURCE_ROOT/SimSys_Scripts/github_scripts/suite_report_git.py $CYLC_WORKFLOW_RUN_DIR/bin
cp $SOURCE_ROOT/SimSys_Scripts/github_scripts/suite_data.py $CYLC_WORKFLOW_RUN_DIR/bin
cp $SOURCE_ROOT/SimSys_Scripts/github_scripts/git_bdiff.py $CYLC_WORKFLOW_RUN_DIR/bin
cp $SOURCE_DIRECTORY/SimSys_Scripts/github_scripts/get_git_sources.py $CYLC_WORKFLOW_RUN_DIR/bin
