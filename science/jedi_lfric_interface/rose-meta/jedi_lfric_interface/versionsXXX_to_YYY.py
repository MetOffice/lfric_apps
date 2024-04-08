from __future__ import absolute_import
#import rose.upgrade
import re
import sys

class UpgradeError(Exception):

    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__
