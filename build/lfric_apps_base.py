#!/usr/bin/env python3

"""
This module contains the LFRicAppsBase class, i.e. the base class
for all LFRic applications in the lfric_apps repository.
"""
import argparse
import logging
import os
from pathlib import Path
import sys
from typing import Iterable, Optional, Tuple, Union

from fab.api import Exclude, git_checkout, Include

from dependency_info import DependencyInfo

logger = logging.getLogger(__name__)

# We need to be able to import LFRicBase from the core
# repository. Scan the command line arguments for --core
# (defaulting to the $LFRIC_CORE environment variable)
# and add this to the Python search path:
arg_parser = argparse.ArgumentParser(add_help=False)
arg_parser.add_argument("--core", "-c", type=str,
                        help="Root of the LFRic core repository",
                        default="$LFRIC_CORE")
args = arg_parser.parse_known_args()[0]   # Ignore element [1]=unknown args

if args.core == "$LFRIC_CORE":
    # Default to root - it will then fail during the import
    core_path = Path(os.environ.get("LFRIC_CORE", "/"))
else:
    core_path = Path(args.core)

sys.path.insert(0, str(core_path / "lfric_build"))

try:
    from lfric_base import LFRicBase
except ModuleNotFoundError as mnfe:
    msg = (f"Cannot import `lfric_base` from '{core_path}'. Specify a valid "
           f"location of the LFRic core repository using the --core command "
           f"line option.")
    logger.error(msg)
    raise RuntimeError(msg) from mnfe


class LFRicAppsBase(LFRicBase):
    '''
    This is the base class for all LFRic FAB scripts.

    :param name: the name to be used for the workspace. Note that
        the name of the compiler will be added to it.
    :param app_dir: the base directory of the application.
    :param root_symbol: the symbol (or list of symbols) of the main
        programs. Defaults to the parameter `name` if not specified.

    '''
    # pylint: disable=too-many-instance-attributes
    def __init__(self, name: str,
                 app_dir: Path,
                 root_symbol: Optional[Union[list[str], str]] = None
                 ):
        self._lfric_apps_root = Path(__file__).parents[1]
        super().__init__(name=name, app_dir=app_dir,
                         root_symbol=root_symbol)

        # Some transmute function will import helper functions from
        # this path.
        # TODO: until we upgrade lfric_core with a new version which
        #       will allow us to use:
        # self.add_python_path(self.lfric_apps_root / "interfaces" /
        self._add_python_paths.append(
            self.lfric_apps_root / "interfaces" / "physics_schemes_interface" /
            "build" / "transmute_psytrans")

        self._dependency_info = DependencyInfo(*self.get_dependencies_info())

    @property
    def lfric_apps_root(self) -> Path:
        """
        :returns: the path to the LFRic_apps root directory
        """
        return self._lfric_apps_root

    @property
    def dependency_info(self) -> DependencyInfo:
        """
        :returns: the dependency info used for this application.
        """
        return self._dependency_info

    def define_command_line_options(
            self,
            parser: Optional[argparse.ArgumentParser] = None
            ) -> argparse.ArgumentParser:
        '''
        This adds LFRic Apps specific command line options. ATM this
        adds the flag --core (which is handled when importing this file,
        but it needs to be defined, otherwise argparse would throw
        an error).

        :param parser: optional a pre-defined argument parser.

        :returns: the argument parser with the LFRic specific options added.
        '''

        parser = super().define_command_line_options(parser)

        parser.add_argument("--core", "-c", type=str,
                            help="Root of the LFRic core repository",
                            default="$LFRIC_CORE")
        return parser

    def get_dependencies_info(self) -> Tuple[Optional[Path], list[str]]:
        """
        This function returns a tuple, consisting of the the path to the
        dependency.yaml file to use, and a list of repository names (from the
        dependency.yaml file). If no dependency.yaml file is required for the
        apps, the first component can be None. If the list of repository names
        is not empty, only entries in the dependencies.yaml file that are
        contained in the specified list will be read. This is frequently used
        by applications that only need some of the repositories listed in
        dependencies.yaml (e.g. typically jules) to avoid checking out all
        other lfric_atm dependencies listed in the yaml file).

        Most applications in lfric_apps now need jules. So, add the
        dependencies.yaml file, but restrict it to the jules repository.

        :returns: the path to the dependencies.yaml file to use and a list
            of repositories to extract.

        """
        return self.lfric_apps_root / "dependencies.yaml", ["jules"]

    def grab_files_step(self) -> None:
        '''
        This method overwrites the base class grab_files_step. It includes all
        the LFRic core directories that are commonly required for building
        LFRic applications. It also grabs optimisation scripts, and the psydata
        directory for profiling, if required.
        '''
        super().grab_files_step()

        # Checkout repositories that are required:
        # ----------------------------------------
        # Call site-specific updates, which allows usage of mirrors.
        self.site_config.update_repos(self.dependency_info)

        for repo in self.dependency_info.get_repo_names():
            if repo in ["lfric_apps", "lfric_core", "SimSys_Scripts"]:
                # For now don't support checking out the apps or core repo
                # (they must be already checked out), and ignore the
                # SimSys_sripts repository, which is not needed.
                logger.info(f"Ignoring repository '{repo}'.")
                continue

            repo_infos = self.dependency_info.get_repo_info(repo)

            for repo_info in repo_infos:
                logger.info(f"Extracting '{repo}' from '{repo_info.source}' "
                            f" to 'science/{repo}', "
                            f"revisions {repo_info.ref}")
                try:
                    git_checkout(self.config,
                                 repo_info.source,
                                 dst_label=f'science/{repo}',
                                 revision=repo_info.ref)
                except RuntimeError as error:
                    logger.error(f"Cannot checkout '{repo}' from "
                                 f"'{repo_info.source}' revision "
                                 f"'{repo_info.ref}': {error}. ")
                    sys.exit(-1)

    def find_source_files_step(
            self,
            path_filters: Optional[Iterable[Union[Exclude, Include]]] = None
            ) -> None:
        """
        Most applications don't need jules, so by default exclude it (any
        # applications that need it must add it to their path_filters

        :param path_filters: optional list of path filters to be passed to
            Fab find_source_files, default is None.
        """
        if path_filters:
            path_filters = list(path_filters)
        else:
            path_filters = []

        if "jules" in self.dependency_info.get_repo_names():
            # In general, jules should not be compiled for most applications
            # (it needs very specific exclude/include specifications), since
            # most applications don't need it. So, we ignore science/jules.
            # The exclude is inserted at the beginning, to allow e.g.
            # lfric_atm to add its own Include to overwrite the exclude
            # (last entry takes precedence in case of same length)
            path_filters.insert(0, Exclude("science/jules"))
        super().find_source_files_step(path_filters)

    def configurator_step(self,
                          include_paths: Optional[list[Path]] = None) -> None:
        """
        Overwrite the configurator step of the base class to
        provide an additional include directory.

        :param include_paths: optional additional include paths
        """
        if include_paths:
            # Make a copy to avoid modifying the caller's list
            include_paths = include_paths[:]
        else:
            include_paths = []

        # Most apps now need jules when running the configurator.
        if "jules" in self.dependency_info.get_repo_names():
            include_paths.append(self.config.source_root / "science" / "jules")
        include_paths.append(self.lfric_apps_root)
        super().configurator_step(include_paths=include_paths)
