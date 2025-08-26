#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

import os
import subprocess

def run_command(command):
    """
    Run a subprocess command and return the result object
    Inputs:
        - command, str with command to run
    """
    result = subprocess.run(
        command.split(),
        capture_output=True,
        text=True,
        timeout=600,
        shell=False,
        check=False,
    )
    print(result.stdout)
    if result.returncode:
        raise RuntimeError(
            f"The command '{command}' failed with error:\n\n{result.stderr}"
        )


print(f"\n[INFO] LFRic Apps Source: {os.environ["SOURCE_LFRIC_APPS"]}")
print(f"[INFO] LFRic Core Source: {os.environ["SOURCE_LFRIC_CORE"]}\n")

# Copy suite_report and fcm_bdiff to bin dir
bin_dir = os.path.join(os.environ["CYLC_WORKFLOW_RUN_DIR"], "bin")
scripts_dir = os.path.join(os.environ["SOURCE_ROOT"], "SimSys_Scripts")
run_command(f"cp {os.path.join(scripts_dir, "suite_report.py")} {bin_dir}")
run_command(f"cp {os.path.join(scripts_dir, "fcm_bdiff.py")} {bin_dir}")

# Make ctldata directory
ctldata_dir = os.path.join(os.environ["SOURCE_ROOT"], "ctldata")
run_command(f'mkdir -p {ctldata_dir}')

# Read through the file defining ctldata sources
# For each source, export into ctldata directory
dests = []
ctldata_file = "ctldata_list.txt"
with open(ctldata_file) as f:
    for line in f:
        line = line.strip()
        line = line.strip("\n")
        if line.startswith("#") or len(line) == 0:
            continue
        source, dest = line.split()
        if dest in dests:
            raise RuntimeError(
                "Two sources in the ctldata_list.txt file share the same "
                f"destination: {dest}. Please deduplicate these."
            )
        dests.append(dest)
        print(f"[INFO] Source: {source}   Dest: {dest}")
        run_command(
            f"fcm export --force {source} {os.path.join(ctldata_dir, dest)}"
        )

