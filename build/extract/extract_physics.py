import argparse
import subprocess
import os
import tempfile
import re
import yaml
from shutil import rmtree
from collections import defaultdict


def run_command(command, shell=False):
    """
    Run a subprocess command and return the result object
    Inputs:
        - command, str with command to run
    """
    if not shell:
        command = command.split()
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=120,
        shell=shell,
        check=False,
    )
    if result.returncode:
        print_error(command, result)


def print_error(command, result):
    """
    Print command and error message if a subprocess fails
    """

    raise RuntimeError(
        f"The command '{command}' failed with error:\n\n{result.stderr}"
    )

def parse_extract_file(fpath):
    """
    Parse the extract file, getting a list of files to extract for each source defined.
    Input:
        - fpath, path to the extract file
    Output:
        - Dictionary of sources. Values is a set of extract files/directories
    """

    with open(fpath, "r") as f:
        extract_file = f.readlines()

    extract_lists = defaultdict(list)

    source = None
    for line in extract_file:
        line = line.strip()
        if not line:
            continue

        if line.startswith("extract.path-incl"):
            source = re.match(r"extract.path-incl\[(\w+)\]", line).group(1)
            # Check for filepath on same line
            line = line.split("=")[1].replace("\\", "").strip()
            if line:
                extract_lists[source].append(line)
        elif source:
            line = line.replace("\\", "").strip()
            extract_lists[source].append(line)

    return extract_lists


def read_dependencies(dependencies_file):
    """
    Read through the dependencies.yaml file, reading in each source. Return a dict of
    source: source_string where each source_string is of format:
    - repo_location::./::ref for git repos, where repo_location is the repo url or path
    to clone, ./ indicates the entire repo, and ref is either the branch name or commit
    hash.
    - repo_location@revision for fcm repos.
    """

    with open(dependencies_file) as stream:
        sources = yaml.safe_load(stream)

    parsed_sources = {}

    for source, values in sources.items():
        source_str = values["source"]
        if not source_str:
            continue
        ref = values["ref"].strip()
        source_str = f"{source_str}::./::{ref}"
        parsed_sources[source] = source_str

    return parsed_sources


def extract_files(dependency, source, files, working):
    """
    Clone the dependency to a temporary location
    Then copy the desired files to the working directory
    Then delete the temporary directory
    """

    source, _, ref = source.split("::")
    tempdir = tempfile.mkdtemp()
    temp_dep = os.path.join(tempdir, dependency)
    working_dep = os.path.join(working, dependency)

    checkout_commands = (
        f"git clone {source} {temp_dep}",
        f"git -C {temp_dep} checkout {ref}"
    )
    for command in checkout_commands:
        run_command(command)

    # make the working directory location
    run_command(f"mkdir -p {working_dep}")

    for extract_file in files:
        source_file = os.path.join(temp_dep, extract_file)
        dest_file = os.path.join(working_dep, extract_file)
        run_command(f"mkdir -p {os.path.dirname(dest_file)}")
        copy_command = f"cp -r {source_file} {dest_file}"
        run_command(copy_command)

    rmtree(tempdir)


def parse_args():
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
        "-w",
        "--working",
        default=".",
        help="Location to perform extract steps in."
    )
    parser.add_argument(
        "-e",
        "--extract",
        default="./extract.cfg",
        help="Path to file containing extract lists"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    extract_lists = parse_extract_file(args.extract)
    dependencies = read_dependencies(args.dependencies)

    for dependency in dependencies:
        if extract_lists[dependency]:
            extract_files(
                dependency,
                dependencies[dependency],
                extract_lists[dependency],
                args.working
            )


if __name__ == "__main__":
    main()
