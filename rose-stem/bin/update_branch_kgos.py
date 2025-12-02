#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Copy failed checksum kgos from a rose-stem suite to a working copy
"""

import sys
import os
import subprocess
import argparse

PLATFORMS = {
    "meto": {
        "azspice": {"copy_command": "cp "},
        "ex1a": {},
        },
    "azngarch": {
        "sandbox": {"copy_command": "cp "},
        },
}

EX1A_PLATFORMS = {
    "exab": {"copy_command": "scp login.exab.sc:"},
    "excd": {"copy_command": "scp login.excd.sc:"},
}


def run_command(command):
    """
    Run a subprocess command and return the result object
    """
    return subprocess.run(
        command, shell=True, capture_output=True, text=True, timeout=60
    )


def get_ex_platform():
    """
    If KGOs are required for ex1a ask the user which platform the tests were run on and
    adapt the PLATFORM dictionary
    """
    zone = input("Which ex1a zone were the tests run on? 1 for EXAB or 2 for EXCD: ")

    if zone == "1" or zone.lower() == "exab":
        PLATFORMS["meto"]["ex1a"] = EX1A_PLATFORMS["exab"]
    elif zone == "2" or zone.lower() == "excd":
        PLATFORMS["meto"]["ex1a"] = EX1A_PLATFORMS["excd"]
    else:
        sys.exit("EX Platform choice not recognised")


def parse_status_file(status_path, job):
    """
    Parse through a job status file for the exit code
    """
    with open(status_path) as status:
        for line in status:
            if "CYLC_JOB_EXIT" in line:
                code = line.split("=")[-1]
                break
        else:
            code = ""
    if "ERR" in code:
        return 1
    if "SUCCEEDED" in line:
        return 0
    sys.exit(f"[FAIL]: Got an unknown exit code '{code}' for job '{job}'")


def find_failed_tasks(log_file_path):
    """
    Look in the cylc-run log/job/ directories for checksum tasks. If this task
    exits with an ERR status record the job as failed.
    """

    failed_checksums = set()
    for job_name in os.listdir(log_file_path):
        task_path = os.path.join(log_file_path, job_name)
        if (
            job_name.startswith("check")
            and "-v-" not in job_name
            and os.path.isdir(task_path)
        ):
            status_file = os.path.join(task_path, "NN", "job.status")
            if parse_status_file(status_file, job_name):
                failed_checksums.add(job_name)
    return failed_checksums


def get_kgo_dirs(job, flow_file):
    """
    Parse through the flow-processed.cylc file and find this kgo jobs stored
    and generated kgo directories.
    """

    new, current = None, None
    in_job_section = False
    with open(flow_file) as workflow:
        for line in workflow:
            line = line.strip()
            if in_job_section:
                if line.startswith("CURRENT_KGO"):
                    current = line.split("=")[-1].strip()
                    current = current.removeprefix("$SOURCE_ROOT/lfric_apps")
                elif line.startswith("NEW_KGO"):
                    new = line.split("=")[-1].strip()
                    new = new.removeprefix("$OUTPUT_ROOT")
                elif line.startswith("[[") and not line.startswith("[[["):
                    sys.exit(
                        "Couldn't identify KGO Directories in the suites "
                        f"flow-processed.cylc file for job '{job}'."
                    )
            if current and new:
                break
            if f"[[{job}]]" in line:
                in_job_section = True
    return current, new


def copy_checksums(job, stored_kgo, new_kgo, suite_name, working_copy, site):
    """
    Generate the command to copy the checksum file from
    """

    for platform in PLATFORMS[site]:
        if platform in job:
            break
    else:
        sys.exit(
            f"[FAIL]: Couldn't find a valid platform for job {job} at " f"site {site}"
        )

    # Join the kgo paths with the path to their base location (wc / cylc suite)
    new_kgo_path = os.path.join(
        "~", "cylc-run", suite_name, "share", "output", new_kgo.strip("/")
    )
    stored_kgo_path = os.path.join(working_copy, stored_kgo.lstrip("/"))

    # Run mkdir -p on the stored kgo path to ensure it exists
    stored_kgo_dir = stored_kgo_path.removesuffix(os.path.basename(stored_kgo_path))
    mkdir_command = f"mkdir -p {stored_kgo_dir}"
    result = run_command(mkdir_command)
    if result.returncode:
        sys.exit(
            f"[FAIL]: Failed while making the directory {stored_kgo_dir} "
            f"with error:\n\n{result.stderr}"
        )

    command = PLATFORMS[site][platform]["copy_command"]
    command += f"{new_kgo_path} {stored_kgo_path}"
    result = run_command(command)
    if result.returncode:
        sys.exit(
            f"[FAIL]: Failed to copy kgo for job {job} with error:\n\n"
            f"{result.stderr}"
        )


def parse_cl_args():
    """
    Parse command line options
    """

    parser = argparse.ArgumentParser(
        "Update failed checksum kgos from a rose-stem run. This script is "
        "expected to be run from the spice/vdi filesystem."
    )
    parser.add_argument(
        "-s",
        "--suite",
        required=True,
        help="The name of the suite being run. Will look in ~/cylc-run for log "
        "files. If the numbered run dir isn't given or is given as runN, will "
        "read the runN symlink and use the real number to work on remote "
        "platforms.",
    )
    parser.add_argument(
        "-w",
        "--working_copy",
        required=True,
        help="The working copy in which to store the kgo. This is expected to "
        "be on the spice/vdi filesystem.",
    )
    parser.add_argument(
        "-p",
        "--site",
        default="meto",
        help="The rose-stem site being run used. Needed for platform settings.",
    )
    args = parser.parse_args()
    args.working_copy = os.path.expanduser(args.working_copy)
    args.suite = args.suite.strip("/")

    # Don't allow use of runN as this symlink doesn't exist on remote platforms
    # If runN or no run provided, read the runN symlink to work out actual num
    path_base = args.suite.split("/")
    if len(path_base) > 1:
        path_base = path_base[-1]
    else:
        path_base = ""
    if path_base == "runN":
        args.suite = args.suite.removesuffix("runN")
        path_base = ""
    if "run" not in path_base:
        spath = os.path.expanduser(os.path.join("~", "cylc-run", args.suite, "runN"))
        sym_path = run_command(f"readlink {spath}")
        args.suite = os.path.join(args.suite, sym_path.stdout.strip("\n"))
    return args


if __name__ == "__main__":
    args = parse_cl_args()

    suite_path = os.path.expanduser(os.path.join("~", "cylc-run", args.suite))
    log_file = os.path.join(suite_path, "log", "job", "1")
    flow_file = os.path.join(suite_path, "log", "config", "flow-processed.cylc")

    failed_jobs = find_failed_tasks(log_file)

    if args.site == "meto" and str(failed_jobs).find("ex1a"):
        get_ex_platform()

    for failed_job in failed_jobs:
        print(f"[INFO]: Copying kgo for {failed_job}")
        stored_kgo_dir, suite_kgo_dir = get_kgo_dirs(failed_job, flow_file)
        copy_checksums(
            failed_job,
            stored_kgo_dir,
            suite_kgo_dir,
            args.suite,
            args.working_copy,
            args.site,
        )

    print(
        "[PASS]: Successfully copied all kgo for all failed tasks. Commit "
        "changes to your branch before review."
    )
