#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Run cylc lint on the rose-stem suite and fail if any issues detected.
Ignore codes S007 (family all uppercase - due to not detecting jinja2
variables), S012 (line length) and S013 (should be 4 tabs)
"""

import sys
import os
import re
import argparse
import subprocess


def run_command(command):
    """
    Launch a subprocess command and return the output
    """
    return subprocess.run(command.split(), capture_output=True, text=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check the cylc lint for a rose-stem suite"
    )
    parser.add_argument(
        "-s",
        "--source",
        help="The path to the rose-stem directory",
        required=True,
    )
    args = parser.parse_args()
    source = os.path.abspath(args.source)
    command = f"cylc lint {source} -n S007 -n S012 -n S013 -n S014"
    result = run_command(command)
    try:
        nissues = int(re.search(r"found (\d+) issue", result.stdout).group(1))
        print(result.stdout, file=sys.stderr)
        sys.exit(
            f"{nissues} Errors were detected - please fix them and run again"
        )
    except AttributeError:
        print(result.stdout)
