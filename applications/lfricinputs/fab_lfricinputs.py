#!/usr/bin/env python3
# ##############################################################################
#  (c) Crown copyright Met Office. All rights reserved.
#  For further details please refer to the file COPYRIGHT
#  which you should have received as part of this distribution
# ##############################################################################

'''A FAB build script for lfricinputs. It relies on the LFRicBase class
contained in the infrastructure directory.
'''

import argparse
import logging
import os
from pathlib import Path
import sys
from typing import Iterable, List, Optional, Union


from fab.steps.grab.folder import grab_folder
from fab.build_config import AddFlags

# We need to import the Apps base class:
sys.path.insert(0, str(Path(__file__).parents[2] / "build"))
from lfric_apps_base import LFRicAppsBase  # noqa: E402

logger = logging.getLogger('fab')

class FabLFRicInputs(LFRicAppsBase):
    '''
    This class builds LFRic inputs. Since LFRic inputs builds
    different binaries in the same tree, it explicitly adds
    the list of target symbols to build.

    :param name: The name of the project directory.
    :param root_symbol: The symbol of the main program(s) to be
        created.
    '''

    def __init__(self, name: str):
        this_file = Path(__file__).resolve()
        # Store the root of this apps for later
        self._this_root = this_file.parent

        # The list of all binaries to compile here, must be
        # set before calling super.__init__ (since it's used in the
        # command line options).
        self.all_binaries = ['lfric2um', 'scintelapi', 'um2lfric']

        super().__init__(name, app_dir=self._this_root)

    def define_command_line_options(
            self,
            parser: Optional[argparse.ArgumentParser] = None
            ) -> argparse.ArgumentParser:
        '''
        This adds LFRic-inputs specific command line options to the base
        class define_command_line_option. Currently, precision-related
        options are added.

        :param parser: optional a pre-defined argument parser.

        :returns: the argument parser with the LFRic specific options added.
        '''
        parser = super().define_command_line_options(parser)

        group = parser.add_argument_group(
            title="LFRic-Inputs",
            description="Arguments to select the binaries to build. If no "
                        "binary is specified, all three will be built." )

        for binary in self.all_binaries:
            group.add_argument(
                f'--{binary}', action="store_true", default=False,
                help=f"Compile '{binary}'.")
        return parser

    def handle_command_line_options(self,
                                    parser: argparse.ArgumentParser) -> None:
        '''
        Analyses the parameter for specifying the binaries to compile.
        Set the selected binaries as root symbols for Fab.

        :param argparse.ArgumentParser parser: the argument parser.
        '''
        super().handle_command_line_options(parser)

        binaries = []
        for binary in self.all_binaries:
            if getattr(self.args, binary):
                binaries.append(binary)
        if not binaries:
            binaries = self.all_binaries
        self.set_root_symbols(binaries)
        logger.info(f"Compiling '{', '.join(binaries)}'.")

    def get_linker_flags(self) -> List[str]:
        '''
        This method adds shumlib to the lfric_base class get_linker_flags
        return.

        :returns: list of flags for the linker.
        '''
        return super().get_linker_flags() + ["shumlib"]

    def grab_files_step(self):
        """
        This method overwrites the base class grab_files_step. It includes
        all source files required for LFRicInputs
        """
        super().grab_files_step()
        dirs = ['applications/lfricinputs/source/',
                'science/gungho/source',
                'science/shared/source/',
                ]

        # pylint: disable=redefined-builtin
        lfric_apps_root = self._this_root.parents[1]
        for dir in dirs:
            grab_folder(self.config, src=lfric_apps_root / dir, dst_label='')

    def get_rose_meta(self) -> Path:
        """
        :returns: The path to the rose meta data config file.
        """
        lfric_apps_root = self._this_root.parents[1]
        return (lfric_apps_root / 'science' / 'gungho' / 'rose-meta' /
                'lfric-gungho' / 'HEAD' / 'rose-meta.conf')

    def analyse_step(self,
                     ignore_dependencies: Optional[Iterable[str]] = None,
                     find_programs: bool = False
                     ) -> None:
        '''
        The method adds lfric_inputs specific list of dependencies to
        shumlib to ignore.

        :param ignore_dependencies: Third party Fortran module names in
            USE statements, 'DEPENDS ON' files and modules to be ignored.
        :param find_programs: if the analyse step should try to automatically
            find all program units to build.

        '''
        inputs_ignore_dependencies = [
            'c_shum_byteswap.o', 'f_shum_ff_status_mod', 'f_shum_field_mod',
            'f_shum_fieldsfile_mod', 'f_shum_file_mod',
            'f_shum_fixed_length_header_indices_mod',
            'f_shum_lookup_indices_mod', 'f_shum_stashmaster_mod'
            ]
        if ignore_dependencies:
            inputs_ignore_dependencies.extend(ignore_dependencies)
        super().analyse_step(ignore_dependencies=inputs_ignore_dependencies,
                             find_programs=find_programs)


# -----------------------------------------------------------------------------
if __name__ == '__main__':

    logger = logging.getLogger('fab')
    logger.setLevel(logging.DEBUG)

    fab_lfric_inputs = FabLFRicInputs(name="lfric_inputs")
    fab_lfric_inputs.build()

    # Rename binaries
    executable_folder_path = fab_lfric_inputs.config.project_workspace
    for binary in fab_lfric_inputs.all_binaries:
        bin_path = executable_folder_path / binary
        if bin_path.exists():
            logger.info(f"Renaming '{binary}' to '{binary}.exe'.")
            bin_path.rename(executable_folder_path / f"{binary}.exe")
