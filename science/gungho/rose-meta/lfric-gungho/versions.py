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


class vn32_t634(MacroUpgrade):
    """Upgrade macro for ticket #634 by Ian Boutle."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t634"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        nml = "namelist:boundaries"
        self.add_setting(config, [nml, "lbc_bal_meth"], "'keep_rho'")
        self.add_setting(config, [nml, "lbc_sort_theta"], ".true.")
        nml = "namelist:initialization"
        eos_height = self.get_setting_value(config, [nml, "model_eos_height"])
        self.remove_setting(config, [nml, "model_eos_height"])
        self.add_setting(config, [nml, "init_eos_height"], eos_height)
        self.add_setting(config, [nml, "init_exner_method"], "'hydrostatic'")
        self.add_setting(config, [nml, "init_sort_theta"], ".true.")
        return config, self.reports


class vn32_t479(MacroUpgrade):
    """Upgrade macro for ticket #479 by Shusuke Nishimoto."""

    BEFORE_TAG = "vn3.2_t634"
    AFTER_TAG = "vn3.2_t479"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        self.add_setting(config, ["namelist:mixing", "fullstress"], ".false.")

        return config, self.reports


class vn32_t46(MacroUpgrade):
    """Upgrade macro for ticket #46 by Shusuke Nishimoto."""

    BEFORE_TAG = "vn3.2_t479"
    AFTER_TAG = "vn3.2_t46"

    def upgrade(self, config, meta_config=None):
        self.remove_setting(config, ["namelist:mixing", "method"])
        self.add_setting(config, ["namelist:mixing", "leonard_tke"], ".true.")
        conf_hash = {}
        conf_list = ["field_names", "enforce_min_value", "min_value"]
        # read list type namelists
        for conf in conf_list:
            conf_val = self.get_setting_value(
                config, ["namelist:transport", conf]
            ).split(",")
            conf_hash[conf] = []
            for value in conf_val:
                if "*" in value:
                    num = int(value.split("*")[0])
                    val = value.split("*")[1]
                    for i in range(num):
                        conf_hash[conf].append(val)
                else:
                    conf_hash[conf].append(value)
        # modify namelist value according to condition
        if (
            "'con_tracer'" in conf_hash["field_names"]
            and "'adv_tracer'" in conf_hash["field_names"]
        ):
            # If both con and adv tracer exist, the former is replaced with
            # pos tracer, and the latter is replaced with gen tracer.
            i = conf_hash["field_names"].index("'con_tracer'")
            j = conf_hash["field_names"].index("'adv_tracer'")
            conf_hash["field_names"][i] = "'pos_tracer'"
            conf_hash["field_names"][j] = "'gen_tracer'"
            conf_hash["enforce_min_value"][i] = ".true."
            conf_hash["enforce_min_value"][j] = ".false."
            conf_hash["min_value"][i] = "0.0"
            conf_hash["min_value"][j] = "-99999999.0"
        elif "'con_tracer'" in conf_hash["field_names"]:
            # If only con tracer exists, it is replaced with pos tracer.
            i = conf_hash["field_names"].index("'con_tracer'")
            conf_hash["field_names"][i] = "'pos_tracer'"
            conf_hash["enforce_min_value"][i] = ".true."
            conf_hash["min_value"][i] = "0.0"
        elif "'adv_tracer'" in conf_hash["field_names"]:
            # If only adv tracer exists, it is replaced with pos tracer.
            i = conf_hash["field_names"].index("'adv_tracer'")
            conf_hash["field_names"][i] = "'pos_tracer'"
            conf_hash["enforce_min_value"][i] = ".true."
            conf_hash["min_value"][i] = "0.0"
        # change namelist value
        for conf in conf_list:
            self.change_setting_value(
                config, ["namelist:transport", conf], ",".join(conf_hash[conf])
            )

        return config, self.reports
