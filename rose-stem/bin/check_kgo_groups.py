#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Launched if a checksum fails.
Fails if the required groups weren't run.
Currently targets all for meto only
"""

import sys
import argparse

SITE_GROUPS = {
    "meto": "all",
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Check required groups have been run in case of kgo failing."
    )
    parser.add_argument(
        "-s", "--site", required=True, help="The site being run at."
    )
    parser.add_argument(
        "-g", "--groups", required=True, help="The groups being run."
    )
    args = parser.parse_args()

    args.groups = args.groups.split(",")

    if args.site not in SITE_GROUPS:
        message = (
            f"{args.site} is not currently checked by the "
            "check_kgo_groups script."
        )
        print(message)
        print(message, file=sys.stderr)
        sys.exit()

    print("Groups Run With:", args.groups)

    required_group = SITE_GROUPS[args.site]
    if required_group not in args.groups:
        message = (
            f"The required group '{required_group}' was not run during "
            "this run of the test suite.\nAs there are kgo failures, "
            "before this ticket can be submitted for review it must be "
            f"run with the '{required_group}' group."
        )
        sys.exit(message)

    print("The testing groups were appropriate.")
