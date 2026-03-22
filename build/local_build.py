#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

"""
Wrapper script for makefiles when doing a local build
Will export a copy of lfric_core using a defined source and rsync it to a
working dir so that incremental builds can occur.
It then runs the makefile for the project being made.
"""

import os
import sys
import subprocess
import argparse
import yaml
import logging
from pathlib import Path
from extract.get_git_sources import clone_and_merge


def subprocess_run(command):
    """
    Run a subprocess command with live output and check the return code
    """

    process = subprocess.Popen(
        command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    # Realtime output of stdout and stderr
    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break
        print(line.decode(), end="", flush=True)

    retcode = process.returncode
    if retcode != 0:
        sys.exit(f"\nError while running subprocess command:\n{command}\n")


def get_root_path():
    """
    Get the root path of the current working copy
    """

    return Path(__file__).absolute().parent.parent


def determine_core_source(root_dir):
    """
    Determine the core source code location from the dependencies file
    Returns an fcm url or a path
    """

    # Read through the dependencies file and populate revision and source
    # variables for requested repo
    with open(root_dir / "dependencies.yaml", "r") as stream:
        dependencies = yaml.safe_load(stream)
    return dependencies["lfric_core"]


def determine_project_path(project, root_dir):
    """
    Determine the path to the makefile for the lfric_apps project being
    built. Defaults to the makefile in the top level if none provided.
    Returns a relative path from this file to the makefile directory
    """

    # Find the project in either science/ interfaces/ or applications/
    for drc in ["science/", "interfaces/", "applications/"]:
        path = root_dir / drc
        for item in os.listdir(path):
            item_path = path / item
            if item_path and item == project:
                return item_path

    sys.exit(
        f"The project {project} could not be found in either the "
        "science/ or applications/ directories in this working copy."
    )


def build_makefile(
    root_dir,
    project_path,
    project,
    working_dir,
    ncores,
    target,
    optlevel,
    psyclone,
    verbose,
):
    """
    Call the make command to build lfric_apps project
    """

    if target == "clean":
        working_path = working_dir
    else:
        working_path = working_dir / f"{target}_{project}"

    print(f"Calling make command for makefile at {project_path}")
    make_command = (
        f"make {target} -C {project_path} -j {ncores} "
        f"WORKING_DIR={working_path} "
        f"CORE_ROOT_DIR={working_dir / 'scratch' / 'lfric_core'} "
        f"APPS_ROOT_DIR={root_dir} "
    )
    if optlevel:
        make_command += f"PROFILE={optlevel} "
    if psyclone:
        make_command += f"PSYCLONE_TRANSFORMATION={psyclone} "
    if verbose:
        make_command += "VERBOSE=1 "

    subprocess_run(make_command)


def main():
    """
    Main function
    """

    parser = argparse.ArgumentParser(
        description="Wrapper for build makefiles for lfric_apps."
    )
    parser.add_argument(
        "project",
        help="project to build. Will search in both science and projects dirs.",
    )
    parser.add_argument(
        "-c",
        "--core_source",
        default=None,
        help="Source for lfric_core. Defaults to looking in dependencies file.",
    )
    parser.add_argument(
        "-m",
        "--mirrors",
        action="store_true",
        help="If true, attempts to use local git mirrors",
    )
    parser.add_argument(
        "--mirror_loc",
        default="/data/users/gitassist/git_mirrors",
        help="Location of github mirrors",
    )
    parser.add_argument(
        "-w",
        "--working_dir",
        default=None,
        type=Path,
        help="Working directory where builds occur. Defaults to the project "
        "directory in the working copy.",
    )
    parser.add_argument(
        "-j",
        "--ncores",
        default=4,
        help="The number of cores for the build task.",
    )
    parser.add_argument(
        "-t",
        "--target",
        default="build",
        help="The makefile target, eg. unit-tests, clean, etc. Default of build.",
    )
    parser.add_argument(
        "-o",
        "--optlevel",
        default=None,
        help="The optimisation to build with, eg. fast-debug. Default of the "
        "the makefile default, usually fast-debug",
    )
    parser.add_argument(
        "-p",
        "--psyclone",
        default=None,
        help="Value passed to PSYCLONE_TRANSFORMATION variable in makefile. "
        "Defaults to the makefile default",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Request verbose output from the makefile ",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # If using mirrors, set environment variable for science extract step
    if args.mirrors:
        os.environ["USE_MIRRORS"] = args.mirror_loc

    # Find the root directory of the working copy
    root_dir = get_root_path()

    # Work out path for the makefile that we are building
    project_path = determine_project_path(args.project, root_dir)

    # Set the working dir default of the project directory
    if not args.working_dir:
        args.working_dir = Path(project_path) / "working"
    else:
        # If the working dir doesn't end in working, set that here
        if not args.working_dir.name == "working":
            args.working_dir = Path(args.working_dir) / "working"
    # Ensure that working_dir is an absolute path and make the directory
    args.working_dir = args.working_dir.resolve()
    args.working_dir.mkdir(parents=True, exist_ok=True)

    # Determine the core source if not provided
    if args.core_source is None:
        core_source = determine_core_source(root_dir)
    else:
        core_source = {"source": args.core_source, "ref": ""}

    if not isinstance(core_source, list):
        core_source = [core_source]

    core_loc = args.working_dir / "scratch" / "lfric_core"
    clone_and_merge("lfric_core", core_source, core_loc, args.mirrors, args.mirror_loc)

    # Build the makefile
    build_makefile(
        root_dir,
        project_path,
        args.project,
        args.working_dir,
        args.ncores,
        args.target,
        args.optlevel,
        args.psyclone,
        args.verbose,
    )


if __name__ == "__main__":
    main()
