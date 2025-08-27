import yaml
import tempfile
import os
from subprocess import run
from shutil import rmtree

def get_dependencies_file(wc_loc):
    """
    Copy the dependencies file to a temporary directory on the local machine that can be
    read.
    """

    print(wc_loc)

    tempdir = tempfile.mkdtemp()

    try:
        host, path = wc_loc.split(":")
        path = os.path.join(path, "dependencies.yaml")
        copy_command = f"scp -o StrictHostKeyChecking=no {host}:"
    except ValueError:
        path = os.path.join(wc_loc, "dependencies.yaml")
        copy_command = f"cp "
    copy_command += f"{path} {tempdir}"

    result = run(
        copy_command.split(), capture_output=True, text=True, timeout=120
    )

    # Raise an error if the returncode is positive
    if result.returncode:
        raise RuntimeError(
            f"An error occured while running the command '{copy_command}' "
            "in order to read the dependencies file. The error message is:\n\n"
            f"'{result.stderr}'"
        )

    return tempdir


def read_sources(clone_source, use_heads):
    """
    Read through the dependencies.yaml file, reading in each source. Return a dict of
    source: source_string where each source_string is of format:
    - repo_location::./::ref for git repos, where repo_location is the repo url or path
    to clone, ./ indicates the entire repo, and ref is either the branch name or commit
    hash.
    - repo_location@revision for fcm repos.
    """

    dependencies_file = get_dependencies_file(clone_source)

    with open(os.path.join(dependencies_file, "dependencies.yaml")) as stream:
        sources = yaml.safe_load(stream)

    parsed_sources = {}

    for source, values in sources.items():
        source_str = values["source"]
        if not source_str:
            continue
        if use_heads:
            ref = "trunk"
        else:
            ref = values["ref"].strip()
        source_str = f"git:{source_str}::./::{ref}"
        parsed_sources[source] = source_str

    rmtree(dependencies_file)

    return parsed_sources
