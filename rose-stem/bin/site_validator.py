#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Run 'cylc --validate' on the rose-stem suite for different sites and groups
Ensure the suite validates for all WORKING_CONFIGS
"""

import sys
import argparse
import subprocess

WORKING_CONFIGS = {
    "meto": [
        "all",
        "build",
        "gungho_model_baroclinic-C24_MG_azspice_gnu_fast-debug-64bit",
        "build_gungho_model_azspice_gnu_fast-debug-64bit",
    ],
    "nci": [
        "gungho_model_gadi_developer",
        "lfric_atm_gadi_developer",
        "lfricinputs_developer",
    ],
    "esnz" : ["all"],
    "ncas" : ["archer2_atm"],
    "uoe" : ["all"],
    "azngarch": ["azngarch_atm"],
}


def run_command(command):
    """
    Launch a subprocess command and return the output
    """
    result = subprocess.run(command.split(), capture_output=True, text=True)
    return result


def generate_validate_command(source, site, group):
    """
    Generate the cylc validate command for this site and group
    """

    fncas = ""
    if site == 'ncas':
        fncas = ("-S HPC_ACCOUNT='hpc_account' "
                 "-S HPC_USERNAME='hpc_username' ")

    command = (
        "cylc validate --debug --check-circular "
        f"-z g={group} "
        f"-S SITE='{site}' "
        f"{fncas}"
        f"{source}/lfric_apps/rose-stem"
    )

    return command


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate rose-stem suites for different sites"
    )
    parser.add_argument(
        "-s",
        "--source",
        help="The code Source",
        required=True,
    )
    args = parser.parse_args()

    failures = False
    for site in WORKING_CONFIGS:
        for group in WORKING_CONFIGS[site]:
            print(f"[INFO] Validating {site} with {group}")
            command = generate_validate_command(args.source, site, group)
            result = run_command(command)
            if result.returncode:
                print(f"[FAIL] {site} with {group} failed to validate")
                print(result.stdout)
                print(result.stderr, file=sys.stderr)
                failures = True
            else:
                print(f"[Pass] {site} with {group} validated successfully")

    if failures:
        sys.exit(1)
