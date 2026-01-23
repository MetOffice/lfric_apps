# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""
Helper functions for cloning git sources in command line builds
"""

import re
import subprocess
from typing import Optional
from pathlib import Path
from shutil import rmtree
import shlex


def get_source(
    source: str,
    ref: str,
    dest: Path,
    repo: str,
    use_mirrors: bool = False,
    mirror_loc: Path = "",
) -> None:

    if ".git" in source:
        if use_mirrors:
            mirror_loc = Path(mirror_loc) / "MetOffice" / repo
            print(f"Cloning/Updating {repo} from mirror {mirror_loc} at ref {ref}")
            clone_repo_mirror(source, ref, repo, mirror_loc, dest)
        else:
            print(f"Cloning/Updating {repo} from {source} at ref {ref}")
            clone_repo(source, ref, dest)
    else:
        print(f"Syncing {repo} at ref {ref}")
        sync_repo(source, ref, dest)


def run_command(
    command: str, rval: bool = False, check: bool = True
) -> Optional[subprocess.CompletedProcess]:
    """
    Run a subprocess command and return the result object
    Inputs:
        - command, str with command to run
    Outputs:
        - result object from subprocess.run
    """
    command = shlex.split(command)
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=300,
        shell=False,
        check=False,
    )
    if check and result.returncode:
        print(result.stdout, end="\n\n\n")
        raise RuntimeError(
            f"[FAIL] Issue found running command {command}\n\n{result.stderr}"
        )
    if rval:
        return result


def clone_repo_mirror(
    repo_source: str, repo_ref: str, parent: str, mirror_loc: Path, loc: Path
) -> None:
    """
    Clone a repo source using a local git mirror.
    Assume the mirror is set up as per the Met Office
    - repo_source: ssh url of the source repository
    - repo_ref: git ref for the source. An empty string will get the default branch
    - parent: Owner of the github repository being cloned (required to construct the
              mirror path)
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

    repo_source = repo_source.removeprefix("git@github.com:")
    user = repo_source.split("/")[0]
    # Check that the user is different to the Upstream User
    if user in parent.split("/")[0]:
        user = None

    # If the ref is a hash then we don't need the fork user as part of the fetch.
    # Equally, if the user is the Upstream User, it's not needed
    if not user or re.match(r"^\s*([0-9a-f]{40})\s*$", repo_ref):
        fetch = repo_ref
    else:
        fetch = f"{user}/{repo_ref}"
    commands = (
        f"git -C {loc} fetch origin {fetch}",
        f"git -C {loc} checkout FETCH_HEAD",
    )
    for command in commands:
        run_command(command)


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

    exclude_dirs = (
        "applications/*/working",
        "applications/*/test",
        "applications/*/bin",
        "science/*/working",
        "science/*/test",
        "science/*/bin",
        "interfaces/*/working",
        "interfaces/*/test",
        "interfaces/*/bin",
        "components/*/working",
        "components/*/test",
        "components/*/bin",
        "infrastructure/*/working",
        "infrastructure/*/test",
        "infrastructure/*/bin",
        "mesh_tools/*/working",
        "mesh_tools/*/test",
        "mesh_tools/*/bin",
    )

    # Trailing slash required for rsync
    command = f"rsync -av {repo_source}/ {loc}"
    for item in exclude_dirs:
        command = f"{command} --exclude '{item}'"
    run_command(command)

    # Fetch the main branch from origin
    # Ignore errors - these are likely because the main branch already exists
    # Instead write them as warnings
    command = f"git -C {loc} fetch origin main:main"
    result = run_command(command, check=False, rval=True)
    if result.returncode:
        print("Warning - fetching main from origin resulted in an error")
        print("This is likely due to the main branch already existing")
        print(f"Error message:\n\n{result.stderr}")

    if repo_ref:
        command = f"git -C {loc} checkout {repo_ref}"
        run_command(command)
