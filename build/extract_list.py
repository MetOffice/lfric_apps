#!/usr/bin/env python3

'''
This module contains a class that reads in fcm extract specifications.
'''

from collections import defaultdict
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Union

import yaml

from fab.steps.find_source_files import Include, Exclude

logger = logging.getLogger(__name__)


class ExtractList():
    '''
    A simple class that reads in an fcm extract.cfg file and stores
    the information about excluded and included file to be used in FAB.
    It can then produce a list of Include/Exclude directive to be used
    by Fab when finding source files.

    There is one difference between the Fcm default behaviour and
    Fab's, tracked in https://github.com/MetOffice/fab/issues/471
    FCM seems to pick the 'most specific' match, while Fab picks
    the last match. E.g.:

        extract.path-excl[shumlib] = common/src/shumlib_version.c
        extract.path-incl[shumlib] = common/src

    FCM would exclude the .c file, while Fab would include it.
    Since #471 might not get fixed (since FCM is outdated), it is
    recommended to just change the order in the extract.cfg files, so
    that the most-specific matches are at the end.

    :param filename: the name of the fcm extract file to read.
    '''

    def __init__(self,
                 filename: Path) -> None:
        self._sections: Dict[str,  List[str]]
        self._sections = defaultdict(list)
        super().__init__()

        with open(filename, "r", encoding="utf8") as stream:
            dependencies = yaml.safe_load(stream)

        for repo in dependencies.keys():
            self._sections[repo] = dependencies[repo]

    def get_all_sections(self) -> Iterable[str]:
        """
        :returns: the list of all sections defined in the FCM Extract file.
        """
        return list(self._sections.keys())

    def get_include_exclude_list(
            self,
            section: str,
            root_path: Path,
            ) -> list[Union[Exclude, Include]]:
        '''
        Converts the information from the read fcm file into a list of
        Include/Exclude directives.

        Note that Fab Include/Exclude classes only use a sub-string tests.
        For example, a line like:

            extract.path-excl[casim] = / # everything

        would actually ignore any file containing 'casim'. Therefore, it is
        required to add a `root_path`, which is the path where the suite is
        checked out in. This `root_path` will be added when specifying the
        matching pattern, e.g. if `root_path="science/casim/src"` the above
        line becomes `science/casim/src`  (and if specific files will be
        ignored, these also use the `root_path`) to avoid mismatches.

        :param section: the name of the section to convert into an
            Include/Exclude list.
        :param root_path: the path under which the suite is checked out.

        :returns: a list with the corresponding include/exclude instances.
        '''

        path_filters: list[Union[Exclude, Include]] = []
        source_file_info = self._sections[section]
        path_filters.append(Exclude(root_path))
        for path in source_file_info:
            if path == "/":
                # Appending Path("something") and  "/" using Path results
                # in just "/", so instead add an empty string
                path = ""
            path_filters.append(Include(root_path / path))

        return path_filters
