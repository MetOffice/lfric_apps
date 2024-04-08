#import rose.upgrade
import re

#from .versionsXX_XX import *


'''
Do not edit any lines below this point - this is the template.
'''


class vnX.X_txxx(rose.upgrade.MacroUpgrade):
    """Upgrade macro for <PROJECT> by <Author>"""
    BEFORE_TAG = "vnX.X_txxx"
    AFTER_TAG = "vnX.X_txxx"
    def upgrade(self, config, meta_config=None):
        """Upgrade a runtime app configuration."""
        # Add settings
        return config, self.reports
