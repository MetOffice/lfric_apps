import argparse
import subprocess
import os
import tempfile
import yaml
from shutil import rmtree
from pathlib import Path
from typing import Dict, List


def run_command(command):
    """
    Run a subprocess command and check output
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


def clone_dependency(values: Dict, temp_dep: Path) -> None:
    """
    Clone the physics dependencies into a temporary directory
    """

    source = values["source"]
    ref = values["ref"]

    commands = (
        f"git -C {temp_dep} init",
        f"git -C {temp_dep} remote add origin {source}",
        f"git -C {temp_dep} fetch origin {ref}",
        f"git -C {temp_dep} checkout FETCH_HEAD"
    )
    for command in commands:
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
        temp_dep.mkdir(parents=True)
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
