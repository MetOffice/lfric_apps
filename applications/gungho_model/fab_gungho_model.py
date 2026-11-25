#!/usr/bin/env python3
# ##############################################################################
#  (c) Crown copyright Met Office. All rights reserved.
#  For further details please refer to the file COPYRIGHT
#  which you should have received as part of this distribution
# ##############################################################################

'''A FAB build script for gungho_model. It relies on the LFRicBase class
contained in the infrastructure directory.
'''

import logging
from pathlib import Path
import sys
from typing import Optional, Union

from fab.steps.grab.folder import grab_folder

# We need to import the Apps base class:
sys.path.insert(0, str(Path(__file__).parents[2] / "build"))
from lfric_apps_base import LFRicAppsBase  # noqa: E402


class FabGungho(LFRicAppsBase):
    """
    A Fab-based build script for Gungho. It relies on the LFRicAppsBase class
    to implement the actual functionality, and only provides the required
    source files.

    :param name: The name of the application.
    :param root_symbol: the symbol (or list of symbols) of the main
        programs. Defaults to the parameter `name` if not specified.
    """

    def __init__(self,
                 name: str,
                 root_symbol: Optional[Union[list[str], str]] = None) -> None:
        this_file = Path(__file__).resolve()
        # Store the root of this apps for later
        self._this_root = this_file.parent
        super().__init__(name=name,
                         app_dir=self._this_root,
                         root_symbol=root_symbol)

    def grab_files_step(self) -> None:
        """
        Grabs the required source files and optimisation scripts.
        """
        super().grab_files_step()
        dirs = ['applications/gungho_model/source/',
                'science/gungho/source',
                'science/shared/source/',
                ]
        # pylint: disable=redefined-builtin
        lfric_apps_root = self._this_root.parents[1]
        for dir in dirs:
            grab_folder(self.config, src=lfric_apps_root / dir, dst_label='')

    def get_rose_meta(self) -> Path:
        """
        :returns: the rose-meta.conf path.
        """
        return (self._this_root / 'rose-meta' / 'lfric-gungho_model' / 'HEAD'
                / 'rose-meta.conf')


# -----------------------------------------------------------------------------
if __name__ == '__main__':

    logger = logging.getLogger('fab')
    logger.setLevel(logging.DEBUG)
    fab_gungo = FabGungho(name="gungho_model")
    fab_gungo.build()
