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

class vn31_t401(MacroUpgrade):
    # Upgrade macro for #401 by Dan Copsey

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t401"

    def upgrade(self, config, meta_config=None):
        # Add settings
        self.add_setting(config, ["namelist:jules_hydrology", "l_inland"], ".false.")
        return config, self.reports

