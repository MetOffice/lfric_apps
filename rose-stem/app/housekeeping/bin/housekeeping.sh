# *****************************COPYRIGHT*******************************
# (C) Crown copyright 2024 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

# Remove everything in work and /source and /data in
# share. If work or share is a symlink then read the link
# and delete from there to ensure cylc clean continues
# working

if [ -n ${CYLC_WORKFLOW_WORK_DIR} ]; then
    if [ -L ${CYLC_WORKFLOW_WORK_DIR} ]; then
        rm -rfv "$(readlink -f $CYLC_WORKFLOW_WORK_DIR)"/*
    else
        rm -rfv $CYLC_WORKFLOW_WORK_DIR/*
    fi
fi
if [ -n ${CYLC_WORKFLOW_SHARE_DIR} ]; then
    if [ -L ${CYLC_WORKFLOW_SHARE_DIR} ]; then
        rm -rfv "$(readlink -f $CYLC_WORKFLOW_SHARE_DIR)"/source
        rm -rfv "$(readlink -f $CYLC_WORKFLOW_SHARE_DIR)"/data
    else
        rm -rfv $CYLC_WORKFLOW_SHARE_DIR/source
        rm -rfv $CYLC_WORKFLOW_SHARE_DIR/data
    fi
fi
