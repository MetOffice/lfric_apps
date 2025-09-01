#!/usr/bin/env python3
"""
Clone sources for a rose-stem run for use with git bdiff module in scripts
"""

import os
import subprocess

REPOS = (
    "lfric_apps",
    "lfric_core",
    "simsys_scripts"
)


def run_command(command, shell=False, rval=False):
    """
    Run a subprocess command and return the result object
    Inputs:
        - command, str with command to run
    Outputs:
        - result object from subprocess.run
    """
    if not shell:
        command = command.split()
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=300,
        shell=shell,
        check=False,
    )
    if result.returncode:
        print(result.stdout, end="\n\n\n")
        raise RuntimeError(
            f"[FAIL] Issue found running command {command}\n\n{result.stderr}"
        )
    if rval:
        return result


def get_repo_info(repo):
    """
    Parse the <REPO>_SOURCE_CODE environment variable provided by the suite to get the
    repo location and ref required. The variable is of rose git source file format
    """

    repo_str = os.environ[f"{repo.upper()}_SOURCE_CODE"].removeprefix("git:")
    repo_str = repo_str.split("::")
    source = repo_str[0]
    try:
        ref = repo_str[2]
    except IndexError:
        ref = None
    return source, ref


def clone_repo(repo_source, repo_ref, loc):
    """
    Clone the repo and checkout the provided ref
    """

    command = f"git clone {repo_source} {loc}"
    run_command(command)
    command = f"git -C {loc} rev-parse --abbrev-ref HEAD"
    result = run_command(command, rval = True)
    branch_name = result.stdout
    parent = "trunk"
    if "simsys_scripts" in repo_source.lower():
        parent = "main"
    commands = [
        f"git -C {loc} checkout {parent}",
        f"git -C {loc} checkout {branch_name}"
    ]
    for command in commands:
        run_command(command)
    if repo_ref:
        command = f"git -C {loc} checkout {repo_ref}"
        run_command(command)


def copy_files(clone_loc):
    """
    The rose-stem suite expects some files to be copied to other locations
    """

    apps = os.path.join(clone_loc, "apps")
    scripts = os.path.join(clone_loc, "simsys_scripts")
    bin_dir = os.path.join(os.environ["CYLC_WORKFLOW_RUN_DIR"], "bin")

    commands = (
        f"cp {os.path.join(scripts, "suite_report_git", "suite_report_git.py")} {bin_dir}",
        f"cp {os.path.join(scripts, "suite_report_git", "suite_data.py")} {bin_dir}",
        f"cp {os.path.join(scripts, "bdiff", "git_bdiff.py")} {bin_dir}",
        f"cp {os.path.join(apps, "dependencies.yaml")} {os.environ["CYLC_WORKFLOW_RUN_DIR"]}"
    )
    for command in commands:
        run_command(command)


def get_ctldata():
    """
    Read in ctldata sources and export them
    #TODO - currently running from fcm for now, need to update this
    """

    # Make ctldata directory
    ctldata_dir = os.path.join(os.environ["SOURCE_ROOT"], "ctldata")
    run_command(f'mkdir -p {ctldata_dir}')

    # Read through the file defining ctldata sources
    # For each source, export into ctldata directory
    dests = []
    ctldata_file = "ctldata_list.txt"
    with open(ctldata_file) as f:
        for line in f:
            line = line.strip()
            line = line.strip("\n")
            if line.startswith("#") or len(line) == 0:
                continue
            source, dest = line.split()
            if dest in dests:
                raise RuntimeError(
                    "Two sources in the ctldata_list.txt file share the same "
                    f"destination: {dest}. Please deduplicate these."
                )
            dests.append(dest)
            print(f"[INFO] Source: {source}   Dest: {dest}")
            run_command(
                f"fcm export --force {source} {os.path.join(ctldata_dir, dest)}"
            )


def main():

    clone_loc = os.environ["SOURCE_ROOT"]

    for repo in REPOS:

        loc = os.path.join(clone_loc, repo.lower().removeprefix("lfric_"))

        repo, ref = get_repo_info(repo)

        clone_repo(repo, ref, loc)

    copy_files(clone_loc)
    get_ctldata()


if __name__ == "__main__":
    main()
