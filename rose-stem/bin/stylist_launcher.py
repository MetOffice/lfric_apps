#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Launch stylist on list of directories. Run on all and print outputs. Fail if any
style changes required.
"""

import sys
import os
import subprocess
import argparse
from stylist import __main__

def launch_stylist(app_path, config_path):
    """
    Launch stylist as a subprocess command and check the output
    """

    command = f"stylist -verbose -configuration {config_path} {app_path}"

    result = subprocess.run(command.split(), capture_output=True, text=True)

    print(result.stdout)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run stylist on all applications. If "
        "application/stylist.py exists use that file, otherwise "
        "use one in rose-stem/app/check_style/file. "
        "Print output, raise error if any style changes required."
    )
    parser.add_argument(
        "-s",
        "--source",
        help="The top level of lfric_apps directory.",
        required=True,
    )
    args = parser.parse_args()

    failed_apps = {}
    for top_dir in ["applications", "science"]:
        applications = os.listdir(os.path.join(args.source, top_dir))
        for app in applications:
            print(f"Running on {app}\n")
            app_path = os.path.join(args.source, top_dir, app)
            config_path = os.path.join(app_path, "stylist.py")
            if not os.path.exists(os.path.join(config_path)):
                config_path = os.path.join(
                    args.source,
                    "rose-stem",
                    "app",
                    "check_style",
                    "file",
                    "stylist.py",
                )
            result = launch_stylist(app_path, config_path)
            if result.returncode:
                failed_apps[app] = result.stderr

    if failed_apps:
        error_message = ""
        for failed in failed_apps:
            error_message += f"Style Errors in {failed}:\n\n"
            error_message += f"{failed_apps[failed]}\n\n\n"
        sys.exit(error_message)
