# -----------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
# -----------------------------------------------------------------------------
'''
Lift options list and similar from each individual script up into this file.
Aim is to allow the creation of a simple global.py which adds OMP over all
loops. The option which matches the file being worked on can be pulled in
and referenced. This reduces the number of files needed
'''

# Needs to be lifted and likely set by the build system longer term.
# The filename passed to PSyclone, this is the pre-processed FTN source.
FILE_EXTEN = ".xu90"

# Basic initialisation, will be used by global script
SCRIPT_OPTIONS_DICT = {}
