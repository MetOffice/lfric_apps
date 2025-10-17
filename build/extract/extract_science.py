import argparse
import subprocess
import os
import re
import tempfile
import yaml
from shutil import rmtree
from pathlib import Path
from typing import Dict, List


def run_command(command):
    """
    Run a subprocess command and return the result object
    Inputs:
        - command, str with command to run
    """
    command = command.split()
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=120,
        shell=False,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(
            f"The command '{command}' failed with error:\n\n{result.stderr}"
        )


def load_yaml(fpath: Path) -> Dict:
    """
    Read in the dependencies.yaml file
    """

    with open(fpath) as stream:
        sources = yaml.safe_load(stream)

    return sources


def clone_dependency(values, temp_dep):
    """
    Clone the physics dependencies into a temporary directory
    """

    source = values["source"]
    ref = values["ref"]

    # Check if it's a hash
    if re.match(r"^\s*([0-9a-f]{40})\s*$", ref):
        commands = (
            f"git clone --depth 1 {source} {temp_dep}",
            f"git -C {temp_dep} fetch --depth 1 origin {ref}",
            f"git -C {temp_dep} checkout {ref}"
        )
        for command in commands:
            run_command(command)
    else: # This is a branch/tag
        command = f"git clone --branch {ref} --depth 1 {source} {temp_dep}"
        run_command(command)


def extract_files(dependency: str, values: Dict, files: List[str], working: Path):
    """
    Clone the dependency to a temporary location
    Then copy the desired files to the working directory
    Then delete the temporary directory
    """

    tempdir = Path(tempfile.mkdtemp())
    if (
        "PHYSICS_ROOT" not in os.environ
        or not Path(os.environ["PHYSICS_ROOT"]).exists()
    ):
        temp_dep = tempdir / dependency
        clone_dependency(values, temp_dep)
    else:
        temp_dep = Path(os.environ["PHYSICS_ROOT"]) / dependency

    working_dep = working / dependency

    # make the working directory location
    working_dep.mkdir(parents=True)

    for extract_file in files:
        source_file = temp_dep / extract_file
        dest_file = working_dep / extract_file
        run_command(f"mkdir -p {dest_file.parents[0]}")
        copy_command = f"cp -r {source_file} {dest_file}"
        run_command(copy_command)

    rmtree(tempdir)


def parse_args() -> argparse.Namespace:
    """
    Read command line args
    """

    parser = argparse.ArgumentParser("Extract physics code for LFRic Builds.")
    parser.add_argument(
        "-d",
        "--dependencies",
        default="./dependencies.yaml",
        help="The dependencies file for the apps working copy.",
    )
    parser.add_argument(
        "-w", "--working", default=".", help="Location to perform extract steps in."
    )
    parser.add_argument(
        "-e",
        "--extract",
        default="./extract.yaml",
        help="Path to file containing extract lists",
    )
    return parser.parse_args()


def main():
    args: argparse.Namespace = parse_args()

    extract_lists: Dict = load_yaml(args.extract)
    dependencies: Dict = load_yaml(args.dependencies)

    for dependency in dependencies:
        if dependency in extract_lists:
            extract_files(
                dependency,
                dependencies[dependency],
                extract_lists[dependency],
                Path(args.working),
            )


if __name__ == "__main__":
    main()
