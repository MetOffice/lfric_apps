import re
import sys

from metomi.rose.upgrade import MacroUpgrade


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn20_t429(MacroUpgrade):
    """Upgrade macro for ticket #429 by Denis Sergeev."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t429"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        # Add the Burrows & Sharp (1999) chemistry scheme option
        # (enabled if chem_scheme='flexchem')
        nml = "namelist:chemistry"
        self.add_setting(config, [nml, "flexchem_opt"], "'bs1999'")
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t429"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        # Blank Upgrade Macro
        return config, self.reports
