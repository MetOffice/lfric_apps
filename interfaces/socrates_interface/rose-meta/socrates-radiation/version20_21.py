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


class vn20_t541(MacroUpgrade):
    """Upgrade macro for ticket #541 by James Manners."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t541"

    def upgrade(self, config, meta_config=None):
        # Commands From: interfaces/socrates_interface/rose-meta/socrates-radiation
        self.add_setting(
            config, ["namelist:radiation", "cloud_entrapment"], "'zero'"
        )
        return config, self.reports


class vn20_t429(MacroUpgrade):
    """Upgrade macro for ticket #429 by Denis Sergeev."""

    BEFORE_TAG = "vn2.0_t541"
    AFTER_TAG = "vn2.0_t429"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/socrates-radiation
        # Add more radiatively active gases: Cs, K, Li, Na, Rb, TiO, VO
        nml = "namelist:radiative_gases"
        self.add_setting(config, [nml, "cs_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "cs_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "cs_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "cs_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "cs_mix_ratio"], "0")
        self.add_setting(config, [nml, "cs_rad_opt"], "'off'")
        self.add_setting(config, [nml, "k_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "k_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "k_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "k_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "k_mix_ratio"], "0")
        self.add_setting(config, [nml, "k_rad_opt"], "'off'")
        self.add_setting(config, [nml, "li_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "li_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "li_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "li_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "li_mix_ratio"], "0")
        self.add_setting(config, [nml, "li_rad_opt"], "'off'")
        self.add_setting(config, [nml, "na_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "na_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "na_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "na_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "na_mix_ratio"], "0")
        self.add_setting(config, [nml, "na_rad_opt"], "'off'")
        self.add_setting(config, [nml, "rb_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "rb_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "rb_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "rb_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "rb_mix_ratio"], "0")
        self.add_setting(config, [nml, "rb_rad_opt"], "'off'")
        self.add_setting(config, [nml, "tio_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "tio_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "tio_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "tio_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "tio_mix_ratio"], "0")
        self.add_setting(config, [nml, "tio_rad_opt"], "'off'")
        self.add_setting(config, [nml, "vo_clim_fcg_levls"], "0")
        self.add_setting(config, [nml, "vo_clim_fcg_nyears"], "0")
        self.add_setting(config, [nml, "vo_clim_fcg_rates"], "0")
        self.add_setting(config, [nml, "vo_clim_fcg_years"], "0")
        self.add_setting(config, [nml, "vo_mix_ratio"], "0")
        self.add_setting(config, [nml, "vo_rad_opt"], "'off'")
        return config, self.reports


class vn20_t649(MacroUpgrade):
    """Upgrade macro for ticket None by None."""

    BEFORE_TAG = "vn2.0_t429"
    AFTER_TAG = "vn2.1"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/socrates-radiation
        # Blank Upgrade Macro
        return config, self.reports
