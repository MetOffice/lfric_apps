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
import subprocess
import argparse
from pathlib import Path


def launch_fortitude(config_path: Path, app_path: Path
                     ) -> subprocess.CompletedProcess[str]:
    """
    Launch fortitude as a subprocess command and check the output
    """

    command: list[str] = ["fortitude", "--config-file", str(config_path),
                          "check", str(app_path)]
    result = subprocess.run(command, capture_output=True, text=True)

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
        "source",
        help="The top level of lfric_apps directory."
    )
    args = parser.parse_args()

    source_path: Path = Path(args.source)

    failed_apps: dict[str, str] = {}

    for top_dir_path in source_path.iterdir():
        # e.g. applications,science,interfaces
        if not top_dir_path.is_dir():  # don't try to loop over files
            continue
        for app_path in top_dir_path.iterdir():
            # e.g. adjoint_tests, adjoint, coupled_interface
            if not app_path.is_dir():
                continue
            app_name: str = app_path.name
            print(f"Running on {app_name}\n")
            config_path: Path = app_path/"fortitude.toml"
            if not config_path.exists():
                print("Using universal config (toml) file."
                      " (Some apps use their own config file.)")
                config_path: Path = (
                     source_path / "rose-stem" / "app"
                     / "check_fortitude_linter" / "file" / "fortitude.toml")
            result: subprocess.CompletedProcess[str] = launch_fortitude(
                    config_path, app_path)
            if result.returncode:
                # prints the app run on if there are errors of any kind
                print(f"Checking: {app_name} \n", file=sys.stderr)
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
                failed_apps[app_name] = result.stderr

    if failed_apps:
        error_message: str = ""
        print("\n\n\nSummary: Fortitude found errors in"
              " the following repositories:\n", file=sys.stderr)
        for failed in failed_apps:
            error_message += f"{failed}\n"
        sys.exit(error_message)
