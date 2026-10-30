#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
##############################################################################

'''
This file contains a simple script that reads a dependencies.yaml file,
and extracts the required LFRic core version into the specified directory.
'''

import logging
from pathlib import Path
import sys

from fab.api import DependencyInfo
from fab.tools.versioning import Git


logger = logging.getLogger("fab")

# Simple main program that checkouts the core repository, given
# a dependency.yaml file (to indicate which version)
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"{sys.argv[0]} path_to_dependency_yaml_file dest_dir")
        print()
        print("This script checks out the specified LFRic core version")
        print("based on the details specified in the dependencies.yaml file.")
        print()
        raise RuntimeError(f"{sys.argv[0]} path_to_dependency_yaml_file "
                           f"dest_dir")

    # Make sure we see logging information from Fab
    logger = logging.getLogger('fab')
    logger.setLevel(logging.DEBUG)

    dep_yaml = Path(sys.argv[1])
    dest_dir = Path(sys.argv[2])
    dep_info = DependencyInfo(dep_yaml, only_repos=["lfric_core"])

    git = Git()
    dest_dir.mkdir(parents=True)
    git.init(dest_dir)
    for repo, repo_infos in dep_info.items():
        for repo_info in repo_infos:
            print(repo, repo_info)
            git.checkout(src=repo_info.source,
                         dst=str(dest_dir),
                         revision=repo_info.ref)
