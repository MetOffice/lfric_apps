# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
"""
Helper functions for cloning git sources in command line builds
"""

import re
import subprocess
from datetime import datetime
from typing import Optional, Union
from pathlib import Path
from shutil import rmtree
import shlex
import logging

logger = logging.getLogger(__name__)


def run_command(
    command: str, check: bool = True, capture: bool = True, timeout: int = 600
) -> Optional[subprocess.CompletedProcess]:
    """
    Run a subprocess command and return the result object
    Inputs:
        - command, str with command to run
    Outputs:
        - result object from subprocess.run
    """

    args = shlex.split(command)

    try:
        result = subprocess.run(
            args,
            capture_output=capture,
            text=capture,
            timeout=timeout,
            shell=False,
            check=False,
        )
        if check and result.returncode != 0:
            err_msg = (result.stderr or "").strip()
            logger.error(f"[FAIL] Command failed: {command}\nError: {err_msg}")
            raise subprocess.CalledProcessError(
                result.returncode, args, output=result.stdout, stderr=result.stderr
            )
        return result

    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.error(f"[FAIL] Execution error for '{args[0]}': {e}")
        raise


def validate_dependencies(dependencies: dict) -> None:
    """
    Check that the dependencies file dictionary matches format expectations.
    Each dictionary value should be a list of dictionaries (or a single dictionary)
    Those dictionaries should have a "source" and a "ref" key
    """
    for item, values in dependencies.items():
        failed = False
        if isinstance(values, dict):
            values = [values]
        if not isinstance(values, list):
            failed = True
        else:
            for entry in values:
                if not isinstance(entry, dict) or (
                    "source" not in entry or "ref" not in entry
                ):
                    failed = True
        if failed:
            raise ValueError(
                f"The dependency {item} does not contain a list of dictionaries (or a "
                "single dictionary) with keys of 'source' and 'ref'.\nPlease edit your "
                "dependencies.yaml file to satisfy this."
            )


def datetime_str() -> str:
    """
    Create and return a datetime string at the current time
    """
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def clone_and_merge(
    dependency: str, opts: Union[list, dict], loc: Path, use_mirrors: bool, mirror_loc: Path
) -> None:
    """
    Wrapper script for calling get_source and merge_source for a single dependency

    dependency: name of the dependency
    opts: dict or list of dicts for a dependency in the dependencies file
    loc: path to location to clone to
    use_mirrors: bool, use local git mirrors if true
    mirror_loc: path to local git mirrors
    """

    if not isinstance(opts, list):
        opts = [opts]

    for i, values in enumerate(opts):
        if values["ref"] is None:
            values["ref"] = ""

        # Clone the first provided source
        if i == 0:
            get_source(
                values["source"],
                values["ref"],
                loc,
                dependency,
                use_mirrors,
                mirror_loc,
            )
        # For all other sources, attempt to merge into the first
        else:
            merge_source(
                values["source"],
                values["ref"],
                loc,
                dependency,
                use_mirrors,
                mirror_loc,
            )


def get_source(
    source: str,
    ref: str,
    dest: Path,
    repo: str,
    use_mirrors: bool = False,
    mirror_loc: Path = Path(""),
) -> None:
    """
    Call functions to clone or rsync git source
    """

    if ".git" in source:
        if use_mirrors:
            logger.info(
                f"[{datetime_str()}] Cloning {repo} from {mirror_loc} at ref {ref}"
            )
            mirror_loc = Path(mirror_loc) / "MetOffice" / repo
            clone_repo_mirror(source, ref, mirror_loc, dest)
        else:
            logger.info(f"[{datetime_str()}] Cloning {repo} from {source} at ref {ref}")
            clone_repo(source, ref, dest)
    else:
        logger.info(f"[{datetime_str()}] Syncing {repo} at ref {ref}")
        sync_repo(source, ref, dest)


def merge_source(
    source: str,
    ref: str,
    dest: Path,
    repo: str,
    use_mirrors: bool = False,
    mirror_loc: Path = Path(""),
) -> None:
    """
    Merge git source into a local git clone. Assumes dest is a git clone that this
    source can be merged into.
    """

    logger.info(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Merging "
        f"{source} at ref {ref} into {repo}"
    )

    if use_mirrors:
        remote_path = Path(mirror_loc) / "MetOffice" / repo
    else:
        remote_path = source
    run_command(f"git -C {dest} remote add local {remote_path}")

    if use_mirrors:
        fetch = determine_mirror_fetch(source, ref)
    else:
        fetch = ref

    run_command(f"git -C {dest} fetch local {fetch}")
    command = f"git -C {dest} merge --no-gpg-sign FETCH_HEAD"
    result = run_command(command, check=False)
    if result.returncode:
        unmerged_files = get_unmerged(dest)
        if unmerged_files:
            handle_merge_conflicts(source, ref, dest, repo)
        else:
            raise subprocess.CalledProcessError(
                result.returncode, command, result.stdout, result.stderr
            )

    # Remove the added remote
    run_command(f"git -C {dest} remote remove local")


def handle_merge_conflicts(source: str, ref: str, loc: Path, dependency: str) -> None:
    """
    If merge conflicts are in `rose-stem/` or `dependencies.yaml` then accept the
    current changes and mark as resolved.
    If others remain then raise an error
    """

    # For suites, merge conflicts in these files/directories are unimportant so accept
    # the current changes
    for filepath in ("dependencies.yaml", "rose-stem"):
        logger.warning(f"Ignoring merge conflicts in {filepath}")
        run_command(f"git -C {loc} checkout --ours -- {filepath}")
        run_command(f"git -C {loc} add {filepath}")

    # Check if there are any remaining merge conflicts
    unmerged = get_unmerged(loc)
    if unmerged:
        files = "\n".join(f for f in unmerged)
        raise RuntimeError(
            "\nA merge conflict has been identified while merging the following branch "
            f"into the {dependency} source:\n\nsource: {source}\nref: {ref}\n\n"
            f"with conflicting files:{files}"
            "\n\nThese will need changing in the source branches to be useable together"
        )


def get_unmerged(loc: Path) -> list[str]:
    """
    Return list of unmerged files in a git clone
    """

    files = run_command(f"git -C {loc} --no-pager diff --name-only --diff-filter=U")
    return files.stdout.split()


def clone_repo_mirror(
    repo_source: str,
    repo_ref: str,
    mirror_loc: Path,
    loc: Path,
) -> None:
    """
    Clone a repo source using a local git mirror.
    Assume the mirror is set up as per the Met Office
    - repo_source: ssh url of the source repository
    - repo_ref: git ref for the source. An empty string will get the default branch
    - mirror_loc: path to the local git mirrors
    - loc: path to clone the repository to
    """

    # If the repository exists and isn't a git repo, exit now as we don't want to
    # overwrite it
    if loc.exists():
        if not Path(loc / ".git").exists():
            raise RuntimeError(
                f"The destination for the clone of {repo_source} already exists but "
                "isn't a git directory. Exiting so as to not overwrite it."
            )

    # Clone if the repo doesn't exist
    else:
        command = f"git clone {mirror_loc} {loc}"
        run_command(command)

    # If not provided a ref, pull the latest repository and return
    if not repo_ref:
        run_command(f"git -C {loc} pull")
        return

    fetch = determine_mirror_fetch(repo_source, repo_ref)
    commands = (
        f"git -C {loc} fetch origin {fetch}",
        f"git -C {loc} checkout FETCH_HEAD",
    )
    for command in commands:
        run_command(command)


def determine_mirror_fetch(repo_source: str, repo_ref: str) -> str:
    """
    Determine the fetch ref for the git mirrors
    """

    repo_source = repo_source.removeprefix("git@github.com:")
    user = repo_source.split("/")[0]
    # Check that the user is different to the Upstream User
    if "MetOffice" in user:
        user = None

    # If the ref is a hash then we don't need the fork user as part of the fetch.
    # Equally, if the user is the Upstream User, it's not needed
    if not user or re.match(r"^\s*([0-9a-f]{40})\s*$", repo_ref):
        fetch = repo_ref
    else:
        fetch = f"{user}/{repo_ref}"

    return fetch


def clone_repo(repo_source: str, repo_ref: str, loc: Path) -> None:
    """
    Clone the repo and checkout the provided ref
    Only if a remote source
    - repo_source: ssh url of the source repository
    - repo_ref: git ref for the source. An empty string will get the default branch
    - loc: path to clone the repository to
    """

    if not loc.exists():
        # Create a clean clone location
        loc.mkdir(parents=True)

        # This process is equivalent to doing a git clone
        # It saves a small amount of space by not fetching all refs
        commands = (
            f"git -C {loc} init",
            f"git -C {loc} remote add origin {repo_source}",
            f"git -C {loc} fetch origin {repo_ref}",
            f"git -C {loc} checkout FETCH_HEAD",
            f"git -C {loc} fetch origin main:main",
        )
        for command in commands:
            run_command(command)
    else:
        commands = (
            f"git -C {loc} fetch origin {repo_ref}",
            f"git -C {loc} checkout FETCH_HEAD",
        )
        for command in commands:
            run_command(command)


def sync_repo(repo_source: str, repo_ref: str, loc: Path) -> None:
    """
    Rsync a local git clone and checkout the provided ref
    """

    # Remove if this clone already exists
    if loc.exists():
        rmtree(loc)

    # Create a clean clone location
    loc.mkdir(parents=True)

    exclude_dirs = []
    try:
        host, path = repo_source.split(":", 1)
        result = run_command(f"ssh {host} git -C {path} status --ignored -s")
    except ValueError:
        # In case the path does not contain `host:` - see if it can be accessed locally
        result = run_command(f"git -C {repo_source} status --ignored -s")
    for ignore_file in result.stdout.split("\n"):
        ignore_file = ignore_file.strip()
        if not ignore_file.startswith("!!"):
            continue
        ignore_file = ignore_file.removeprefix("!!").strip()
        exclude_dirs.append(ignore_file)

    # Trailing slash required for rsync
    command = f"rsync -av {repo_source}/ {loc}"
    for item in exclude_dirs:
        command = f"{command} --exclude '{item}'"
    run_command(command)

    # Fetch the main branch from origin
    # Ignore errors - these are likely because the main branch already exists
    # Instead write them as warnings
    command = f"git -C {loc} fetch origin main:main"
    result = run_command(command, check=False)
    if result and result.returncode:
        logger.warning(
            "Fetching main from origin resulted in an error."
            "This is likely due to the main branch already existing"
            f"\nError message:\n\n{result.stderr}"
        )

    if repo_ref:
        command = f"git -C {loc} checkout {repo_ref}"
        run_command(command)


def set_https(dependencies: dict) -> dict:
    """
    Change sources in a dependencies dictionary to use https instead of ssh
    """

    logger.info("Modifying Dependencies to use https")
    for dependency, opts in dependencies.items():
        if not isinstance(opts, list):
            opts = [opts]
        for i, values in enumerate(opts):
            if values["source"].startswith("git@github.com:"):
                values["source"] = values["source"].replace(
                    "git@github.com:", "https://github.com/"
                )
                opts[i] = values
        dependencies[dependency] = opts

    return dependencies
