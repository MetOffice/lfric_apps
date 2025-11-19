##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

import yaml
import tempfile
import os
from subprocess import run
from shutil import rmtree

def get_dependencies_file(wc_loc):
    """
    Copy the dependencies file to a temporary directory on the local machine that can be
    read.
    """

    tempdir = tempfile.mkdtemp()

    try:
        host, path = wc_loc.split(":", 1)
        path = os.path.join(path, "dependencies.yaml")
        copy_command = f"scp -o StrictHostKeyChecking=no {host}:"
    except ValueError:
        path = os.path.join(wc_loc, "dependencies.yaml")
        copy_command = "cp "
    copy_command += f"{path} {tempdir}"

    result = run(
        copy_command.split(), capture_output=True, text=True, timeout=120
    )

    # Raise an error if the returncode is positive
    if result.returncode:
        raise RuntimeError(
            f"An error occured while running the command '{copy_command}' "
            "in order to read the dependencies file. The error message is:\n\n"
            f"'{result.stderr}'"
        )

    return tempdir

def read_sources(clone_source, repo, use_heads):
    """
    Load the dependencies.yaml file as a dictionary
    """

    dependencies_file = get_dependencies_file(clone_source)

    with open(os.path.join(dependencies_file, "dependencies.yaml")) as stream:
        dependencies = yaml.safe_load(stream)

    if not dependencies[repo]["source"]:
        dependencies[repo]["source"] = clone_source

    # Populate parent, assume MetOffice is owner if not set
    for dependency, values in dependencies.items():
        if "parent" not in values:
            dependencies[dependency]["parent"] = f"MetOffice/{dependency}.git"
        if use_heads:
            dependencies[dependency]["ref"] = ""

    rmtree(dependencies_file)

    return dependencies
