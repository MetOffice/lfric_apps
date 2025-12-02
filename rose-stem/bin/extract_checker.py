#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Script to look through build output for lfric_atm (intel fast-debug) and
check files being extracted from feeder repos.
Two modes:
 * check - '-m check' Will check that no files are extracted and unused. Will
                      raise an error if this is the case.
 * generate - '-m generate' Will generate a list of files that need to be
                            extracted for each repo. Can copy these into the
                            extract-list file manually.
"""
import os
import sys
import re
import argparse
from collections import defaultdict

# List of repos which files are extracted from
REPOS_LIST = ["um", "ukca", "jules", "socrates", "casim"]

# Files which need to be treated manually for each repo - eg. some header files
MANUAL_FILES = {
    "um": [
        "um/src/atmosphere/atmosphere_service/include/qsat_mod_qsat.h",
        "um/src/atmosphere/atmosphere_service/include/qsat_mod_qsat_mix.h",
        "um/src/atmosphere/large_scale_precipitation/include/lsp_moments.h",
        "um/src/atmosphere/large_scale_precipitation/include/lsp_subgrid_lsp_qclear.h",
    ],
}


def file_loader(fpath):
    """
    Return list of lines from file
    """
    file_contents = []
    with open(fpath, "r") as in_file:
        for line in in_file:
            file_contents.append(line.strip())
    return file_contents


def pattern_search(inlist, pattern, mode="whole"):
    """
    Search through a given input for the given pattern.
    Returns a dictionary where the key is each found item that matches and the
    value is a tally of how many times that appears.
    """
    output = defaultdict(int)
    for item in inlist:
        if (
            pattern.lower() in item.lower()
            and mode == "whole"
            or item.lower().startswith(pattern)
            and mode == "start"
        ):
            output[item] += 1
    return output


def extract_paths(inlist):
    """
    Search each string in the list for non-whitespace with /'s as a crude
    path definition. Check if this path has been found before and flag if it
    has before storing in a set based on whether there is a file extension
    at the end or not.

    Returns two sets for file paths and directory paths
    """
    file_paths = set()
    directory_paths = set()
    for item in inlist:
        search_output = re.search(r"(\S+/\S*)", item)
        if search_output is not None:
            path = search_output.group(1)
            # If we find a duplicate it suggests the regex has gone wrong
            if (path in file_paths) or (path in directory_paths):
                sys.exit(f"Unexpected duplicate path found: {path}")
            if re.search(r"\.[A-z0-9]+", path):
                file_paths.add(path)
            else:
                directory_paths.add(path)

    return file_paths, directory_paths


def things_in_A_but_not_B(setA, setB):
    """
    Return a set of elements that are in set A, but not in set B
    """
    return setA.difference(setB)


def things_only_in_A_AND_B(setA, setB):
    """
    Return a set of elements that are in both set A and set B, not just one of
    them
    """
    return setA.intersection(setB)


def things_only_in_A_OR_B(setA, setB):
    """
    Return a set of elements that are only in one of sets A or B, but not both.
    """
    return setA.symmetric_difference(setB)


def subsetting_by_pattern(Set_In, pattern):
    """
    Return a set of elements from the input set which pattern match with
    supplied "pattern".
    """
    pattern = re.compile(pattern)
    subset = set([])
    for thing in Set_In:
        if pattern.match(thing):
            subset.add(thing)
    return subset


def print_extract_list(source_repos, extract_list):
    """
    For a given repo print the list of files required to be extracted
    """
    print("\n\nCreating an extract list of " f"{source_repos} files....")

    contents = f"extract.path-incl[{source_repos}] ="
    print(f"\n\n{contents:<88}\\")

    extract_list_sorted = sorted(extract_list[source_repos])
    for elem in extract_list_sorted[:-1]:
        print(f"    {elem.split('/', 1)[-1]:<84}\\")
    if len(extract_list_sorted) > 0:
        print(f"    {extract_list_sorted[-1].split('/', 1)[-1]}")


def generate_extract_lists(set_of_copied_files, set_of_compiled_files):
    """
    Generate a set of files required to be extracted for each repo
    """
    compiled_files_to_extract = things_only_in_A_AND_B(
        set_of_copied_files, set_of_compiled_files
    )

    extract_list = defaultdict(str)

    print(
        "\n\nGenerating and printing extract lists which can be copied into "
        "lfric_atm/fcm-make/extract.yaml."
    )

    print(
        f"\n\nThere are {len(compiled_files_to_extract)} files compiled "
        "that need to be extracted\nOf which..."
    )

    for source_repos in REPOS_LIST:
        match_set = subsetting_by_pattern(
            compiled_files_to_extract, source_repos
        )
        print(f"    {len(match_set)} were in {source_repos}")
        if source_repos in MANUAL_FILES:
            for item in MANUAL_FILES[source_repos]:
                match_set.add(item)
        extract_list[source_repos] = match_set

    for source_repos in REPOS_LIST:
        print_extract_list(source_repos, extract_list)


def check_for_uncompiled_files(set_of_copied_files, set_of_compiled_files):
    """
    Check for files that are extracted but not required. If any exist,
    print them out and raise an error
    """
    # Create a set of files copied but not compiled
    copied_not_compiled = things_in_A_but_not_B(
        set_of_copied_files, set_of_compiled_files
    )

    # For each repo generate set of uncompiled files and count total found
    # Only check repos included in REPOS_LIST defined above
    uncompiled_files = {}
    num_uncompiled = 0
    for source_repos in REPOS_LIST:
        match_set = subsetting_by_pattern(copied_not_compiled, source_repos)
        uncompiled_files[source_repos] = match_set
        num_uncompiled += len(match_set)

    # If we've found uncompiled files print out errors and file paths
    if num_uncompiled > 0:
        print(
            f"There are {num_uncompiled} files copied but not compiled\n"
            "Please remove any uncompiled files\n"
            "See output for the full list",
            file=sys.stderr,
        )

        for source_repos in REPOS_LIST:
            print(
                f"    {len(uncompiled_files[source_repos])} in "
                f"{source_repos}",
                file=sys.stderr,
            )

        for source_repos in uncompiled_files:
            if len(uncompiled_files[source_repos]) > 0:
                print(
                    f"\n\nFiles in {source_repos} that were copied but not "
                    "compiled:"
                )
                for elem in sorted(uncompiled_files[source_repos]):
                    print(f"    {elem}")

        sys.exit(1)
    else:
        print("No files extracted but not compiled")


def main(filepath, mode):
    """
    Main function
    """

    # read build output
    file_contents = file_loader(filepath)

    # Make the lists of files we believe were extracted (copied) and compiled
    copied_list = pattern_search(file_contents, "src/", mode="start")
    compiled_list = pattern_search(file_contents, "ompile*")

    # Need to extract the file paths from the list of compiled file output
    set_of_compiled_files, _ = extract_paths(
        [item.split("/", 2)[-1] for item in compiled_list.keys()]
    )

    # Using path extractor to split into files and dirs, plus return sets.
    set_of_copied_files, set_of_copied_dirs = extract_paths(
        [item.split("/", 1)[-1] for item in copied_list.keys()]
    )

    if mode == "check":
        check_for_uncompiled_files(set_of_copied_files, set_of_compiled_files)
    elif mode == "generate":
        generate_extract_lists(set_of_copied_files, set_of_compiled_files)
    else:
        print(f"Unrecognised mode requested: {mode}", file=sys.stderr)
        sys.exit(1)


def debugging_prints(set_of_copied_files, set_of_compiled_files):
    """
    Some print statements useful for debugging
    """
    print("\n\nTaking a peek in the list of copied files....")
    count = 1
    for elem in iter(set_of_copied_files):
        print(f"    {elem} is an example of a copied_list entry")
        count += 1
        if count == 10:
            break

    print("\n\nTaking a peek in the list of compiled files....")
    count = 0
    for elem in iter(set_of_compiled_files):
        print(f"    {elem} is an example of a compiled_list entry")
        count += 1
        if count == 10:
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check through output from lfric_atm build task for "
        "files extracted and compiled. Two modes - one to "
        "generate list of all required extractions, one to "
        "test no unnecessary extractions are occuring.",
    )
    parser.add_argument(
        "-u",
        "--user",
        default=None,
        help="Username suite was run under (default = lfric)",
    )
    parser.add_argument(
        "-s",
        "--suite_name",
        default=None,
        help="Name of suite run "
        + "(default = lfric_atm-heads-meto-azspice-nightly)",
    )
    parser.add_argument(
        "-d", "--suite_dir", default=None, help="Directory of rose suite"
    )
    parser.add_argument(
        "-m",
        "--mode",
        default="check",
        help="Mode to run - either generate or check",
        choices=["check", "generate"],
    )
    parser.add_argument(
        "-p",
        "--platform",
        default="azspice",
        help="Which platform the build job was run from on",
    )
    args = parser.parse_args()

    suite_default = f"lfric_atm-heads-{args.platform}-nightly"

    if args.suite_name and not args.user:
        user = "~"
        suite_name = args.suite_name
    elif not args.suite_name and not args.user:
        user = "~lfric"
        suite_name = suite_default
    elif args.user and not args.suite_name:
        parser.error(
            "\nCan't handle user (-u) without suite_name (-s)\n"
            "Please provide suite_name"
        )
    else:
        user = args.user
        suite_name = args.suite_name
        if not user.startswith("~"):
            user = "~" + user

    # Create path to build output file.
    log_file_path = (
        "log/job/1/build_lfric_atm_spice_intel_fast-debug-64bit/NN/job.out"
    )
    if args.suite_dir:
        filepath = args.suite_dir
    else:
        filepath = os.path.join(user, "cylc-run", suite_name)
        print(f"Using the Suite: {filepath}")

    filepath = os.path.join(os.path.expanduser(filepath), log_file_path)

    main(filepath, args.mode)
