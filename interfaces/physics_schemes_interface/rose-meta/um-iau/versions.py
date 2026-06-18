import re
import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version30_31 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn31_t496(MacroUpgrade):
    """Upgrade macro for ticket #496 by Samantha Pullen."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t496"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-iau
        # Add new setting to iau namelist
        self.add_setting(config, ["namelist:iau", "iau_outerloop"], ".false.")

        return config, self.reports
