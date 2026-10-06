#! /usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This module contains the app-specific configuration class for NCI Gadi.
'''


from site_specific.nci_gadi.config import Config as NciSiteConfig
from app_specific.default.config import Config as DefaultAppsConfig


class Config(DefaultAppsConfig, NciSiteConfig):
    '''
    This class is the application-specific NCI Configuration for Fab builds.
    ATM it is empty, and will rely on the default application configuration
    or the site-specific NCI-Gadi settings
    '''
