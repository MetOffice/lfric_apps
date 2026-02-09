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


class vn22_t771(MacroUpgrade):
    """Upgrade macro for ticket #771 by josephwallwork."""

    BEFORE_TAG = "vn2.2"
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


class vn22_t797(MacroUpgrade):
    """Upgrade macro for ticket #797 by Charlotte Norris."""

    BEFORE_TAG = "vn2.2_t771"
    AFTER_TAG = "vn2.2_t797"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        """
        Add fastjx_numwavel, fjx_solcyc_type, fjx_scat_file, fjx_solar_file,
        fastjx_dir, fjx_spec_file to namelist chemistry
        """
        self.add_setting(
            config, ["namelist:chemistry", "fastjx_numwavel"], "18"
        )
        self.add_setting(config, ["namelist:chemistry", "fjx_solcyc_type"], "0")
        self.add_setting(
            config, ["namelist:chemistry", "fjx_scat_file"], "'FJX_scat.dat'"
        )
        self.add_setting(
            config,
            ["namelist:chemistry", "fjx_solar_file"],
            "'FJX_solcyc_May17.dat'",
        )
        self.add_setting(
            config,
            ["namelist:chemistry", "fjx_spec_file"],
            "'FJX_spec_Nov11.dat'",
        )
        self.add_setting(
            config,
            ["namelist:chemistry", "fastjx_dir"],
            "'$UMDIR/vn13.9/ctldata/UKCA/fastj'",
        )
        return config, self.reports


class vn22_t618(MacroUpgrade):
    """Upgrade macro for ticket #618 by Mohit Dalvi."""

    BEFORE_TAG = "vn2.2_t797"
    AFTER_TAG = "vn2.2_t618"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        """
        Add Photolysis scheme and FastJX namelist items to chemistry namelist.
        Value of photol_scheme is based on chem_scheme choice.
        Change the value of fastjx_dir to replace hardwired path with
        link from central ctldata extract
        """
        chem_scheme = self.get_setting_value(
            config, ["namelist:chemistry", "chem_scheme"]
        )
        if chem_scheme == "'strattrop'":
            phot_scheme = "'fastjx'"
        else:
            phot_scheme = "'off'"  # Default
        self.add_setting(
            config, ["namelist:chemistry", "photol_scheme"], phot_scheme
        )
        self.add_setting(
            config, ["namelist:chemistry", "chem_timestep"], "3600"
        )
        self.add_setting(config, ["namelist:chemistry", "fastjx_mode"], "1")
        self.add_setting(
            config, ["namelist:chemistry", "fastjx_prescutoff"], "20.0"
        )
        self.add_setting(
            config, ["namelist:chemistry", "fjx_solcyc_months"], "0"
        )
        self.add_setting(config, ["namelist:chemistry", "fjx_solcyc_type"], "0")
        self.change_setting_value(
            config, ["namelist:chemistry", "fastjx_dir"], "'fastj'"
        )
        return config, self.reports


class vn22_t34(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.2_t618"
    AFTER_TAG = "vn3.0"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/um-chemistry
        # Blank Upgrade Macro
        return config, self.reports
