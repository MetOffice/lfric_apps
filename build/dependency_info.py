#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
##############################################################################

'''
This module contains a class that manages the dependencies specified in
a dependencies.yaml file.
'''

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union
import yaml


@dataclass
class RepoInfo:
    """
    A simple class to maintain source and ref attributes of a repository.
    It maintains a list of sources and references, the idea being that the
    first one is checked out, then the following are merged into the checked
    out version.

    This data is maintained in source_ref, a list of SourceRef instances.

    Each RepoInfo also maintains one path, the location in which the source
    code is checked out to.
    """

    @dataclass
    class SourceRef:
        """
        A simple data class that stores a source and ref, and allows
        to access and update them individually.
        """
        source: str
        ref: str

    # The URL and reference of all sources
    source_ref: list[SourceRef] = field(default_factory=list)

    # The location on the file system where the checked out files are stored.
    # Provide a default, this will be filled in later by the build script
    # with the absolute path of the checked out version.
    path: Path = Path(".")

    def append(self, source: str, ref: str) -> None:
        """
        Provide a simple append function, to add a source and reference.

        :param source: the source to add.
        :param ref: the reference to use in this source
        """
        self.source_ref.append(RepoInfo.SourceRef(source, ref))

    def __iter__(self) -> SourceRef:
        """
        This function allows to iterate over all sources/references
        of this dependency.
        """
        for item in self.source_ref:
            yield item


# ============================================================================
class DependencyInfo(dict):
    '''
    A simple dictionary-like class that stores the version information
    from a yaml file:

        casim:
            source: git@github.com:MetOffice/casim.git
            ref: 2025.12.1
        ...

    The information can be accessed as a dictionary, e.g.:
        gr = DependencyInfo("$LFRIC_APPS_SRC/dependencies.yaml")
        gr["casim"] --> {"source": "git@.../casim.git",
                         "ref": "2025.12.1"}

    A filter can be specified to restrict the repositories that
    are being handled.

    The constructor will check that each dependency has indeed
    source and ref defined (note that for lfric_apps these are
    defined, but empty, indicating to use the current directory).

    If the requested section does not exist, a key error is raised.

    :param filename: The path to the dependencies.yaml file.
    '''

    def __init__(self, filename: Optional[Union[str, Path]],
                 only_repos: Optional[list[str]]) -> None:
        super().__init__()

        # If there are no dependencies, just return (this object will
        # then represents no dependencies).
        if not filename:
            return

        with open(filename, "r", encoding="utf8") as stream:
            dependencies = yaml.safe_load(stream)

        for repo, all_deps in dependencies.items():
            if only_repos and repo not in only_repos:
                continue
            # A repo can either have a single definition, or a list
            # Support both:
            if not isinstance(all_deps, list):
                all_deps = [all_deps]

            self[repo] = RepoInfo()
            for dep in all_deps:
                if "source" not in dep:
                    raise RuntimeError(f"'{filename} does not contain a "
                                       f"'source' definition for repo "
                                       f"'{repo}'.")
                if "ref" not in dep:
                    raise RuntimeError(f"'{filename} does not contain a "
                                       f"'ref' definition for repo '{repo}'.")
                self[repo].append(dep["source"], dep["ref"])

    def get_repo_names(self) -> list[str]:
        """
        :returns: the list of all repositories stored in this object.
        """
        return list(self.keys())

    def get_repo_info(self, repo: str) -> RepoInfo:
        """
        :returns: the list of repository infos for a given dependency.

        :raises:KeyError if the repository is not defined.
        """
        return self[repo]
