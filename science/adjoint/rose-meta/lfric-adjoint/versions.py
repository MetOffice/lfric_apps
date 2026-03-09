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
        exec_name = self.get_setting_value(
            config, ["env", "EXEC_NAME"]
        )
        # To prevent macro upgrade errors, we edit the adjoint_tests config manually
        do_not_upgrade = ["adjoint_tests"]
        if exec_name not in do_not_upgrade:
          # Adds new namelist entry alphabetically
          source = self.get_setting_value(
              config, ["file:configuration.nml", "source"]
          )
          # Insert adjoint above aerosol except for these exceptions
          exception_exec_names = ["jedi_forecast", "jedi_forecast_pseudo"]
          if exec_name in exception_exec_names :
            source = re.sub(
                r"namelist:base_mesh",
                r"namelist:adjoint" + "\n" + "  namelist:base_mesh",
                source,
            )
          else:
            source = re.sub(
                r".namelist:aerosol.",
                r" namelist:adjoint" + "\n" + " (namelist:aerosol)",
                source,
            )
          self.change_setting_value(
              config, ["file:configuration.nml", "source"], source
          )

          # Default value
          self.add_setting(
              config, ["namelist:adjoint", "l_compute_annexed_dofs"], ".true."
          )

        return config, self.reports
