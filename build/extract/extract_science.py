import argparse
import subprocess
import os
import yaml
from pathlib import Path
from typing import Dict, List
from get_git_sources import get_source
import shlex


def run_command(command):
    """
    Run a subprocess command and check output
    Inputs:
        - command, str with command to run
    """
    command = shlex.split(command)
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


def extract_files(dependency: str, values: Dict, files: List[str], working: Path):
    """
    Clone the dependency to a temporary location
    Then copy the desired files to the working directory
    Then delete the temporary directory
    """

    mirror_loc: Path = os.getenv("MIRROR_LOC", "")

    if (
        "PHYSICS_ROOT" not in os.environ
        or not Path(os.environ["PHYSICS_ROOT"]).exists()
    ):
        clone_loc = working.parent / "scratch" / dependency
        get_source(
            values["source"],
            values["ref"],
            clone_loc,
            dependency,
            bool(mirror_loc),
            mirror_loc
        )
    else:
        clone_loc = Path(os.environ["PHYSICS_ROOT"]) / dependency


    # make the working directory location
    working_dir = working / dependency
    working_dir.mkdir(parents=True, exist_ok=True)

    # rsync extract files from clone loc to the working directory
    copy_command = "rsync --include='**/' "
    for extract_file in files:
        if not extract_file:
            continue
        if Path(clone_loc / extract_file).is_dir():
            extract_file = extract_file.rstrip("/")
            extract_file += "/**"
        copy_command += f"--include='{extract_file}' "
    copy_command += f"--exclude='*' -avmq {clone_loc}/ {working_dir}"
    run_command(copy_command)


def parse_args() -> argparse.Namespace:
    """
    Read command line args
    """

    parser = argparse.ArgumentParser("Extract physics code for LFRic Builds.")
    parser.add_argument(
        "-d",
        "--dependencies",
        default="./dependencies.yaml",
        help="The dependencies file for the apps working copy",
    )
    parser.add_argument(
        "-w", "--working", default=".", help="Build location"
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
