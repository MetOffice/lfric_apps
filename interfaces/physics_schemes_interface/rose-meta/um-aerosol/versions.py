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
    Upgrade macro for ticket #856 by Alan J Hewitt.

    Scaffolding for simpler but faster RADAER.

    Currently not for production.
    """

    BEFORE_TAG = "vn3.3"
    AFTER_TAG = "vn3.3_t856"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-aerosol

        # Default to false
        self.add_setting(
            config, ["namelist:aerosol", "l_radaer_f_not_for_production"],
            ".false."
        )

        # Default to one
        self.add_setting(
            config, ["namelist:aerosol", "n_radaer_faster_seg_size"], "1"
        )

        return config, self.reports
