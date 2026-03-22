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


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""


class vn31_t322(MacroUpgrade):
    """Upgrade macro for ticket #322 by Terence Vockerodt."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t322"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-adjoint
        # Adds new namelist entry alphabetically
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        if "namelist:adjoint" not in source:
            # Insert adjoint to configuration
            for line in source.split("\n"):
                namelist = line.strip("()")
                namelist = namelist.strip()
                if "namelist:adjoint" < namelist:
                    source = re.sub(
                        line,
                        rf" namelist:adjoint\n{line}",
                        source,
                    )
                    break
            self.change_setting_value(
                config, ["file:configuration.nml", "source"], source
            )
        # Default value
        self.add_setting(
            config, ["namelist:adjoint", "l_compute_annexed_dofs"], ".true."
        )

        return config, self.reports
