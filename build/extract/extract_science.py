import argparse
import os
import yaml
from pathlib import Path
from get_git_sources import clone_and_merge, run_command
import logging


def load_yaml(fpath: Path) -> dict:
    """
    Read in the dependencies.yaml file
    """

    with open(fpath) as stream:
        sources = yaml.safe_load(stream)

    return sources


def extract_files(dependencies: dict, extract_lists: dict, working: Path) -> None:
    """
    Clone the dependency to a temporary location
    Then copy the desired files to the working directory
    Then delete the temporary directory
    """

    mirror_loc = os.getenv("MIRROR_LOC", "")
    use_mirrors = bool(mirror_loc)
    mirror_loc = Path(mirror_loc)

    for dependency, sources in dependencies.items():
        if dependency not in extract_lists:
            continue
        files = extract_lists[dependency]

        # If the PHYSICS_ROOT environment variable is provided, then use sources there
        if "PHYSICS_ROOT" in os.environ and Path(os.environ["PHYSICS_ROOT"]).exists():
            clone_loc = Path(os.environ["PHYSICS_ROOT"]) / dependency
        else:
            clone_loc = working.parent / "scratch" / dependency
            clone_and_merge(dependency, sources, clone_loc, use_mirrors, mirror_loc)

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
    parser.add_argument("-w", "--working", default=".", help="Build location")
    parser.add_argument(
        "-e",
        "--extract",
        default="./extract.yaml",
        help="Path to file containing extract lists",
    )

    args = parser.parse_args()
    args.working = Path(args.working)
    return args


def main():
    args: argparse.Namespace = parse_args()

    logging.basicConfig(level=logging.INFO)

    extract_lists: dict = load_yaml(args.extract)
    dependencies: dict = load_yaml(args.dependencies)

    extract_files(dependencies, extract_lists, args.working)


if __name__ == "__main__":
    main()
