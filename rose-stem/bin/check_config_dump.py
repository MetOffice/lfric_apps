#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Run a rose config-dump command on all relevant lfric_apps directories.
Fail if any files changed as result of command
'''

import sys
from os import path
import subprocess
import argparse

APPLICATIONS = [
    "applications/jedi_lfric_tests",
    "applications/gravity_wave",
    "applications/gungho_model",
    "applications/lfric_atm",
    "applications/linear_model",
    "applications/name_transport",
    "applications/ngarch",
    "applications/shallow_water",
    "applications/solver",
    "applications/transport",
    "science/linear",
    "rose-stem",
]


def check_config_dump(path):
    """
    Run rose config-dump for a given application in a subprocess
    Check output and return Fail/Pass
    """

    command = f"rose config-dump -C {path}".split()
    result = subprocess.run(command, capture_output=True, text=True)
    if "[INFO] M" in result.stdout:
        return False
    else:
        return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check through applications in lfric_apps and check "
                    "that rose config-dump has been run on them."
    )
    parser.add_argument(
        "-s", "--source", help="Source directory for lfric_apps"
    )
    args = parser.parse_args()

    failed_applications = []

    for application in APPLICATIONS:
        print(f"Checking Application: {application}")
        app_dir = path.join(args.source, application)
        passed = check_config_dump(app_dir)
        if passed:
            print(f"{application} passed")
        else:
            print(f"{application} failed")
            failed_applications.append(application)

    if failed_applications:
        error_string = (
            "rose config-dump needs to be run on the following "
            "directories:\n" + "\n".join(failed_applications)
        )
        sys.exit(error_string)
