import sys

from metomi.rose.upgrade import MacroUpgrade

from .version21_22 import *


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


class vn22_t885(MacroUpgrade):
    """Upgrade macro for ticket #885 by Samantha Pullen."""

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t885"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Add iau_sst_path to files namelist"""
        self.add_setting(config, ["namelist:files", "iau_sst_path"], "''")
        """Add iau_sst to section_choice namelist"""
        self.add_setting(
            config, ["namelist:section_choice", "iau_sst"], ".false."
        )
        return config, self.reports


class vn22_t4661(MacroUpgrade):
    """Upgrade macro for ticket #4661 by Denis Sergeev."""

    BEFORE_TAG = "vn2.2_t885"
    AFTER_TAG = "vn2.2_t4661"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(config, ["namelist:extrusion", "eta_values"], "''")
        return config, self.reports


class vn22_t771(MacroUpgrade):
    """Upgrade macro for ticket #771 by josephwallwork."""

    BEFORE_TAG = "vn2.2_t4661"
    AFTER_TAG = "vn2.2_t771"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        """Add new namelist options for chemistry timestep halving"""
        self.add_setting(
            config,
            ["namelist:chemistry", "i_chem_timestep_halvings"],
            value="0",
        )
        return config, self.reports


class vn22_t621(MacroUpgrade):
    """Upgrade macro for ticket #621 by Mike Hobson."""

    BEFORE_TAG = "vn2.2_t771"
    AFTER_TAG = "vn2.2_t621"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-lfric2lfric
        """Add orography, and boundary namelists"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:files",
            r"namelist:files" + "\n" + " (namelist:orography)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )

        return config, self.reports


class vn22_t887(MacroUpgrade):
    """Upgrade macro for ticket #887 by Mike Whitall."""

    BEFORE_TAG = "vn2.2_t621"
    AFTER_TAG = "vn2.2_t887"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-cloud
        nml = "namelist:cloud"
        self.add_setting(config, [nml, "dbsdtbs_turb_0"], "1.5E-4")
        self.add_setting(config, [nml, "i_pc2_erosion_numerics"], "'implicit'")

        return config, self.reports
