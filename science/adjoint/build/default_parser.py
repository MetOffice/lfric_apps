#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Contains shared argument parser in the PSyAD build scripts.
'''

import argparse


def default_parser() -> argparse.ArgumentParser:
    '''
    Sets up and returns the default parser.
    '''
    formatter = argparse.RawDescriptionHelpFormatter
    return argparse.ArgumentParser(description=__doc__,
                                   formatter_class=formatter)
