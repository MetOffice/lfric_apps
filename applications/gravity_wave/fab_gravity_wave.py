#!/usr/bin/env python3
# ##############################################################################
#  (c) Crown copyright Met Office. All rights reserved.
#  For further details please refer to the file COPYRIGHT
#  which you should have received as part of this distribution
# ##############################################################################

'''A FAB build script for applications/gravity_wave. It relies on the LFRicBase
class contained in the infrastructure directory.
'''

import logging
from pathlib import Path
import sys
from typing import Optional, Union

from fab.api import grab_folder

# We need to import the Apps base class:
sys.path.insert(0, str(Path(__file__).parents[2] / "build"))

from lfric_apps_base import LFRicAppsBase  # noqa: E402
from pfunit_mixin import PfUnitMixin  # noqa: E402


class FabGravityWave(PfUnitMixin, LFRicAppsBase):
    """
    A Fab-based build script for GravityWave. It relies on the LFRicAppsBase
    class to implement the actual functionality, and only provides the
    required source files.

    :param name: The name of the application.
    :param root_symbol: the symbol (or list of symbols) of the main
        programs. Defaults to the parameter `name` if not specified.
    """

    def __init__(self,
                 name: str,
                 root_symbol: Optional[Union[list[str], str]] = None) -> None:

        print("fabgravitywave")
        super().__init__(name=name,
                         apps_dir=Path(__file__).parent,
                         root_symbol=root_symbol)
        # Store the root of this apps for later
        this_file = Path(__file__).resolve()
        self._this_root = this_file.parent

    def grab_files_step(self) -> None:
        """
        Grabs the required source files and optimisation scripts.
        """
        super().grab_files_step()
        dirs = ['applications/gravity_wave/source/',
                'science/gungho/source',
                'science/shared/source/',
                ]

        # pylint: disable=redefined-builtin
        lfric_apps_root = self._this_root.parents[1]
        for dir in dirs:
            grab_folder(self.config, src=lfric_apps_root / dir, dst_label='')

        # Copy the optimisation scripts into a separate directory
        grab_folder(self.config, src=self._this_root / "optimisation",
                    dst_label='optimisation')

    def get_rose_meta(self) -> Path:
        """
        :returns: the rose-meta.conf path.
        """
        return (self._this_root / 'rose-meta' / 'lfric-gravity_wave' / 'HEAD'
                / 'rose-meta.conf')


# -----------------------------------------------------------------------------
if __name__ == '__main__':

    logger = logging.getLogger('fab')
    logger.setLevel(logging.DEBUG)
    fab_gravity_wave = FabGravityWave(name="gravity_wave")
    fab_gravity_wave.build()
