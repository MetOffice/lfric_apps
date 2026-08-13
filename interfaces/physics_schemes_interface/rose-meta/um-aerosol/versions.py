import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version31_32 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn32_t698(MacroUpgrade):
    """
    Upgrade macro for ticket #698 by Alan J Hewitt.

    Users can now select UKCA GLOMAP setting via namelist.
    Note that dust_and_clim is treated by UKCA as if setting (6)
    but is treated by RADAER as if setting (8).

    Users can now set dust ageing via the namelist.

    Users can now set radaer SUstrat via the namelist.

    """

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t698"

    def upgrade(self, config, meta_config=None):
        # Add settings

        glomap_mode = self.get_setting_value(config,
                                             ["namelist:aerosol",
                                              "glomap_mode"])

        if glomap_mode == "'dust_and_clim'":
            # Existing suites with dust_and_clim need mode_setup==6
            mode_setup = 'DUonly_2mode'
        elif glomap_mode == "'ukca'":
            # Existing suites with ukca need mode_setup==8
            mode_setup = 'SUBCOCSSDU_7mode'
        elif glomap_mode == "'radaer_test'":
            # Existing suites with ukca need mode_setup==8
            mode_setup = 'SUBCOCSSDU_7mode'
        elif glomap_mode == "'climatology'":
            # This is trigger ignored
            mode_setup = ''
        elif glomap_mode == "'off'":
            # This is trigger ignored
            mode_setup = ''
        else:
            # This is trigger ignored
            mode_setup = ''

        # Add new settings with the specified option
        self.add_setting( config,
                          ["namelist:aerosol","mode_setup"],
                          mode_setup )

        # Default to false since this is the setting in all existing tests
        self.add_setting( config,
                          ["namelist:aerosol", "l_dust_mp_ageing"],
                          ".false." )

        # Default to false since this is the setting in all existing tests
        self.add_setting( config,
                          ["namelist:aerosol", "l_ukca_radaer_sustrat"],
                          ".false." )

        return config, self.reports
