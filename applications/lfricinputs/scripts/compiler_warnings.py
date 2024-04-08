#!/usr/bin/env python3
#
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""
Creates a list of all compiler warnings output from rose-stem and stats
on type. Uses a branch diff to search for warnings that are specific to the
files that have changed if on a branch.

For details on inputs and how to use please reference
https://code.metoffice.gov.uk/trac/um/wiki/CompilerWarnings
"""

import re
import os
import subprocess
import argparse

# Configure the different searches needed to extract warnings from the
# log files.Each warning pattern has an associated number of lines the warning
# is printed over.
WARNING_PATTERNS = [
    (re.compile(r'warning #(?P<code>\d+): (?P<desc>.*\.)', re.IGNORECASE), 3),
    (re.compile(r'remark #(?P<code>\d+): (?P<desc>.*\.)', re.IGNORECASE), 3),
    (re.compile(r'warning:(?P<desc>.*) \[enabled by default\]'), 3),
    (re.compile(r'warning:(?P<desc>.*) \[(?P<code>.*)\]', re.IGNORECASE), 5),
    (re.compile(r'Warning:(?P<desc>.*)'), 5),
    (re.compile(r'warning:(?P<desc>.*)', re.IGNORECASE), 1),
]

FILE_PATTERNS_ALL = [
    (re.compile(r'[>>&2].*/src/.*'), 0),
    (re.compile(r'[>>&2].*/include/.*'), 0)
]

FILE_PATTERNS_LI = [
    (re.compile(r'[>>&2].*lfricinputs/src/.*'), 0)
]


class CompilerWarningError(Exception):
    pass


def setup_specific_file_search(branch,
                               print_all_trunk,
                               print_all_branch,
                               file_search=''):
    """
    Configure the different searches needed to extract file names from the log
    files. If a file search has been provided then use this pattern. If this is
    a branch (not the trunk) and if print_all_branch is False,then this is a
    list of all source files containing changed. Otherwise if either
    print_all_trunk or print_all_branch is True then search for any
    source file.
    :param branch: str
    :param print_all_trunk: bool
    :param print_all_branch: bool
    :param file_search: str
    :return: file_search_pattern: list of tuples containing a regex search and
             a flag
    """

    branch_info = get_branch_info(branch)

    print(f'branch = {branch}')

    file_search_pattern = []
    flag = 0  # 0 shows this number is a placeholder, not relevant information

    if file_search:
        file_search_pattern = [(re.compile(r'\[>>&2\].*' + file_search), 0)]
    elif 'trunk' in branch_info['url']:
        if print_all_trunk:
            file_search_pattern = FILE_PATTERNS_ALL
    else:
        if print_all_branch:
            file_search_pattern = FILE_PATTERNS_ALL
        else:
            print('Fetching branch diff')
            files = get_changed_source_files(branch, branch_info['parent'])
            if files:
                print('Outputting warnings related to the following files:\n' +
                      '\n'.join(files))

                # create custom search term for each file found in the diff.
                # Interest is only in the non-"info" lines as these contain
                # warnings and errors.

                # search pattern should be a list of tuples to match the
                # warning list
                file_search_pattern = [(re.compile(r'\[>>&2\].*' + item), flag)
                                       for item in files]

            else:
                print('No source changes found')

    return file_search_pattern


def get_branch_diff(branch):
    """
    Use the bdiff command to extract a list of files that have changed
    on the user's branch

    :param branch: str
    :return: a list of files changed on the branch
    """

    bdiff = subprocess.run('fcm bdiff ' + branch + ' --summarize',
                           shell=True,
                           stdout=subprocess.PIPE,
                           universal_newlines=True)

    bdiff.check_returncode()

    bdiff_files = bdiff.stdout.strip().split('\n')

    return bdiff_files


def get_branch_info(branch):
    """
    Function uses the fcm binfo command to extract the relative
    root name of the repository and the full branch name.

    :param branch: str
    :return: dictionary of required branch info
    """

    search = [(re.compile('^Branch Parent: (?P<parent>.*?)@'), 'parent'),
              (re.compile('^URL: (?P<url>.*)'), 'url')]

    binfo = subprocess.run('fcm binfo ' + branch,
                           stdout=subprocess.PIPE,
                           shell=True,
                           universal_newlines=True)

    binfo.check_returncode()

    branch_info = {}

    binfo_out = binfo.stdout.strip().split('\n')
    for line in binfo_out:
        for pattern, group_name in search:
            match = pattern.search(line)
            if match:
                branch_info[group_name] = match.group(group_name)

    if not branch_info:
        raise CompilerWarningError('fcm binfo returned no information')

    return branch_info


def get_changed_source_files(branch, branch_parent):
    """
    Get the branch diff and the directory path, then extract a list of
    files that have changed.

    :param branch_parent: str
    :param branch: str
    :return: a list of file names
    """

    bdiff_files = get_branch_diff(branch)

    if len(bdiff_files) == 1 and bdiff_files[0] == '':
        print('"fcm bdiff" command returned nothing. This may be an '
              'error, or this may be because you have no changes on '
              'your branch.')
        return []
    else:
        bdiff_files = [bfile[8:] for bfile in bdiff_files
                       if bfile[0] == 'M' or bfile[0] == 'A']

        bdiff_files = [os.path.relpath(bfile,
                                   branch_parent)
                   for bfile in bdiff_files]

        # Remove any folder changes (check for '.' rather than file existence
        # as this is agnostic to the directory structure.)
        bdiff_files = [file for file in bdiff_files if
                       ('.' in os.path.basename(file))]

        # Remove any non-source changes
        bdiff_files = [file for file in bdiff_files if 'src' in file]

        return bdiff_files


def print_warning(line):
    """
    Checks the line provided is not an info line and then prints it
    :param line: str
    :return:
    """
    if line.startswith('[>>&2]'):
        print(f'{line.strip()}')
    else:
        raise CompilerWarningError('Tried to print non-warning ' + line.strip())


def print_title(title):
    """
    Prints the supplied string as a title
    :param title: str
    :return:
    """
    print('\n\n' + '=' * len(title))
    print(title)
    print('=' * len(title))


def print_table(warnings_found):
    """
    Create a summary table of the number of each type of warning found
    :param warnings_found: dictionary of warning codes/descriptions
    :return:
    """

    warning_count = sum([count for _, (_, count) in warnings_found.items()])

    if warning_count:
        # Create a summary table of the number of each type of warning found
        print('\n{} warning types found over {} warnings:'
              .format(len(warnings_found), warning_count))

        print('Code                     : No     : Example Description')
        for warning, (desc, count) in sorted(warnings_found.items()):
            print(
                f'{warning:<24} : {count:<6} : {desc}')

    else:
        print('No warnings found')

    return warning_count


def search_line(line_list, search_pattern):
    """
    Searches each file line supplied for a match with each of the search
    patterns supplied.
    :param line_list: list of lines taken from the logfile
    :param search_pattern: list of tuples containing a search pattern and an
    int. If the int is 0 then this is a flag to say this number is not relevant.
    :return: match object and a count (either the integer supplied with the
    search pattern, or the index of which line in the list the match was found
    in).
    """
    for index, line in enumerate(line_list):
        for pattern, line_count in search_pattern:
            match = pattern.search(line)
            if match:
                # if line_count is non-zero assume this is the return value
                # otherwise return index.
                if line_count:
                    count = line_count
                else:
                    count = index
                return match, count

    return None, 0


def create_logfile_list(cylc_run_dir):
    """
    Search the cylc-run output for compiler output files and create a list of
    these.
    :param cylc_run_dir: str directory name of top level cylc-run folder
    :return: list of logfile paths
    """
    # Create a list of logfiles from compiler output
    logfiles = []
    logfile_path = os.path.join(cylc_run_dir, 'share')
    for directory in os.listdir(logfile_path):
        if directory.startswith('fcm_make') and 'rigorous' in directory:
            logfile = os.path.join(logfile_path, directory, 'fcm-make.log')
            if os.path.isfile(logfile):
                logfiles.append(logfile)
            logfile = os.path.join(logfile_path, directory, 'fcm-make2.log')
            if os.path.isfile(logfile):
                logfiles.append(logfile)

    logfiles.sort()

    if not logfiles:
        print('No rigorous compiler output found in '
                                   + logfile_path + '\nExiting script.')
        # No point continuing as no logfiles found. Exit without error.
        exit(0)

    return logfiles


def extract_warning_description(match, warnings_found, type_to_match):
    """
    Check if we've seen this warning before and add to dictionary formatted
    {code: description}. If no specific code found then group under "Other"
    instead.

    :param match: match object containing code and description information
    :param warnings_found: dictionary of all warning types found so far.
    :param type_to_match: string containing a warning code to check for. False
           if all warnings are relevant.
    :return: warnings_found is modified and used again outside this function.
             True if type_to_match is found in the warning code or type_to_match
             is False.
    """

    # Not all warning strings contain a code to identify the warning by. Group
    # these together with a code of "Other"
    if 'code' not in match.re.groupindex:
        warning_code = 'Other'
    else:
        warning_code = match.group('code')

    # If this is a new warning type then add to the dictionary
    if warning_code not in warnings_found:
        warnings_found[warning_code] = [match.group('desc'), 0]

    # Keep a count of how many times this warning type has occurred
    warnings_found[warning_code][1] += 1

    # If type_to_match is False then all warnings are relevant so return True
    if not type_to_match:
        return True

    # if the type_to_match string is found in the warning code then flag a match
    if type_to_match in warning_code:
        return True

    return False


def main():
    parser = argparse.ArgumentParser(usage='Search each compiler output from '
                        'the rose stem for compiler warnings. In each case '
                        'create a summary table of the types of warning found '
                        'and print to stdout the full details of each warning '
                        'as required.')
    parser.add_argument('--branch', help='working copy location',
                        default=os.path.normpath(
                                os.path.join(os.getcwd(), '..')))
    parser.add_argument('--cylc_run_dir', help='top level cylc-run directory',
                        default=None)
    parser.add_argument('--print_all_trunk', action='store_true',
                        help='Print all available warnings if on the trunk, '
                             'otherwise print no warnings.')
    parser.add_argument('--print_all_branch', action='store_true',
                        help='Print all available warnings if on a branch, '
                             'otherwise print warnings specific to files that'
                             ' have been changed.')
    parser.add_argument('--file_search', help='A string containing part of a '
                        'file path. Outputs warnings which relate to files that'
                        ' match this search criteria. Overrides print_all_trunk'
                        ' and print_all_branch.')
    parser.add_argument('--warning_code', help='A string containing the warning'
                        ' type which should be on the output.', default=False)
    args = parser.parse_args()

    # Get branch path from mirror if possible
    if "SOURCE_UM_MIRROR" in os.environ:
        args.branch = os.getenv("SOURCE_UM_MIRROR")

    # If cylc_run_dir is not specified, then assume it is in the users home
    # directory and the folder name matches the branch name
    if args.cylc_run_dir is None:
        args.cylc_run_dir = os.path.join(
            os.getenv('HOME'), 'cylc-run', os.path.basename(args.branch))

    logfiles = create_logfile_list(args.cylc_run_dir)
    file_pattern_spec = setup_specific_file_search(args.branch,
                                                   args.print_all_trunk,
                                                   args.print_all_branch,
                                                   args.file_search)

    lines_to_store = max([line_count for _, line_count in WARNING_PATTERNS])

    # Open and process each available logfile in turn
    for logfile in logfiles:
        with open(logfile, errors='replace') as log:
            print_title('Processing ' + logfile)

            # Reinitialise variables for each logfile
            warnings_found = {}
            lines_to_print = 0
            lines_in_store = 0
            line_store = []
            processing_warning = 0
            code_match = False
            skipped_warnings = 0

            # Process logfile line by line
            for line in log:
                # Each warning should be contained on consecutive lines, so
                # a non-warning triggers a clean start.
                if line.startswith('[info]'):
                    lines_to_print = 0
                    lines_in_store = 0
                    line_store = []
                    processing_warning = 0

                else:
                    # Keep a buffer of enough lines to capture a full
                    # warning for each file reference.
                    line_store.insert(0, line)
                    if len(line_store) > lines_to_store:
                        line_store.pop()

                    # A warning string has previously been found, the
                    # current line is part of this warning.
                    if processing_warning:
                        lines_in_store += 1

                    # Check current line for warning string
                    match, warning_lines = search_line([line], WARNING_PATTERNS)

                    # If match has found a "missing terminating character"
                    # warning then skip. Keep track of how many are skipped.
                    if match:
                        if "missing terminating" in line:
                            match = None
                            skipped_warnings = skipped_warnings + 1

                    # If warning string found check current store for any
                    # previous source file reference.
                    # If one is found then lines_in_store will contain
                    # details of its location and the warning will be
                    # assumed to start from here, otherwise it will be 0 to
                    # indicate the first line is the one containing the
                    # warning string.
                    # Setup lines_to_print based on which warning string was
                    # identified.
                    if match:
                        processing_warning = 1

                        _, lines_in_store = \
                            search_line(line_store, FILE_PATTERNS_ALL)

                        lines_to_print = warning_lines

                    # Once the whole warning should be in the store
                    # search for either a generic or specific file
                    # reference. If one is found then print the whole
                    # warning to stdout. NB lines_in_store is zero indexed,
                    # lines_to_print is a simple count.
                    if lines_in_store == (lines_to_print - 1):
                        lfricinp, _ = search_line(line_store, FILE_PATTERNS_LI)

                        if lfricinp:
                            code_match = extract_warning_description(match,
                                                            warnings_found,
                                                            args.warning_code)
                        if file_pattern_spec:
                            file_match, _ = search_line(line_store,
                                                        file_pattern_spec)
                        else:
                            file_match = None

                        if file_match and code_match:
                            while lines_to_print > 0:
                                print_warning(line_store[lines_to_print - 1])
                                lines_to_print -= 1

                        lines_to_print = 0
                        lines_in_store = 0
                        line_store = []
                        processing_warning = 0

            # Final summary table
            print_table(warnings_found)

            if skipped_warnings:
                print(f'\nNB: {skipped_warnings} missing terminating '
                      f'character warnings not included.')


if __name__ == "__main__":
    main()
