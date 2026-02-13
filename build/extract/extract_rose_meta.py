import argparse
import logging
import os
import tempfile
import yaml
from pathlib import Path
from shutil import copytree, rmtree
from get_git_sources import clone_and_merge


def load_yaml(fpath: Path) -> dict:
    """
    Read in the dependencies.yaml file
    """

    with open(fpath) as stream:
        sources = yaml.safe_load(stream)

    return sources


def parse_args() -> argparse.Namespace:
    """
    Read command line args
    """

    parser = argparse.ArgumentParser("Extract physics code for LFRic Builds.")
    parser.add_argument(
        "-r",
        "--repo",
        default="jules",
        help="The repository to get the rose-meta directories from",
    )
    parser.add_argument(
        "-d",
        "--dependencies",
        default="./dependencies.yaml",
        help="The dependencies file for the apps working copy",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=".",
        type=Path,
        help="Where to output the rose-meta directories to",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO)

    dependencies = load_yaml(args.dependencies)
    sources = dependencies[args.repo]

    mirror_loc = os.getenv("MIRROR_LOC", "")
    use_mirrors = bool(mirror_loc)
    mirror_loc = Path(mirror_loc)
    tidy_clone_loc = False

    # If the PHYSICS_ROOT environment variable is provided, then use sources there
    # Otherwise clone and merge based on dependencies file
    if "PHYSICS_ROOT" in os.environ and Path(os.environ["PHYSICS_ROOT"]).exists():
        clone_loc = Path(os.environ["PHYSICS_ROOT"]) / args.repo
        tidy_clone_loc = True
    else:
        clone_loc = Path(tempfile.mkdtemp()) / args.repo
        clone_and_merge(args.repo, sources, clone_loc, use_mirrors, mirror_loc)

    # Copy the rose-meta contents to the output directory
    rose_meta_orig = clone_loc / "rose-meta"
    for directory in rose_meta_orig.iterdir():
        copytree(directory, args.output / directory.name, dirs_exist_ok=True)

    if tidy_clone_loc:
        rmtree(clone_loc)


if __name__ == "__main__":
    main()
