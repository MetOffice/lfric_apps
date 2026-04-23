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


class vn31_t226(MacroUpgrade):
    """Upgrade macro for ticket #226 by Chris Smith."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t226"

    def upgrade(self, config, meta_config=None):
        self.add_setting(config, ["namelist:initial_temperature", "profile_variable"], "'potential'")
        self.add_setting(config, ["namelist:initial_vapour", "profile_variable"], "'mr'")
        return config, self.reports
