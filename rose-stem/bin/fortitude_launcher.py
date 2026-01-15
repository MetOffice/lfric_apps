#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

"""
Launch fortitude on list of directories. Run on all and print outputs.
Fail if any style changes required.
"""

import sys
import os
import subprocess
import argparse


def launch_fortitude(config_path, app_path):
    """
    Launch fortitude as a subprocess command and check the output
    """

    command = f"fortitude --config-file {config_path} check {app_path}"
    result = subprocess.run(command.split(), capture_output=True, text=True)

    print(result.stdout)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run fortitude on all applications. If "
        "application/fortitude.toml exists use that file, otherwise "
        "use one in rose-stem/app/check_fortitude_linter/file. "
        "Print output, raise error if any changes required."
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
            config_path = os.path.join(app_path, "fortitude.toml")
            if not os.path.exists(os.path.join(config_path)):
                print('''Using universal config (toml) file.
                      (Some apps use their own config file.)''')
                config_path = os.path.join(
                    args.source,
                    "rose-stem",
                    "app",
                    "check_fortitude_linter",
                    "file",
                    "fortitude.toml",
                )
            result = launch_fortitude(config_path, app_path)
            if result.returncode:
                # prints the app run on if there are errors of any kind
                print(f"Checking: {app} \n", file=sys.stderr)
                if not result.stderr:
                    # prints if no other/config errors are found
                    print("Found lint errors:", file=sys.stderr)
                    # prints the lint errors
                    print(result.stdout, file=sys.stderr)
                if result.stderr:
                    # prints if there are other/config errors
                    print("Found non-lint errors: \n", file=sys.stderr)
                    # prints the other/config errors
                    print(result.stderr, "\n\n\n", file=sys.stderr)
                failed_apps[app] = result.stderr

    if failed_apps:
        error_message = ""
        print('''\n\n\nSummary: Fortitude found errors in
              the following repositories:\n''', file=sys.stderr)
        for failed in failed_apps:
            error_message += f"{failed}\n"
        sys.exit(error_message)
