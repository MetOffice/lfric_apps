##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# There are two import files for UM source. Therefore we will have two versions of the same
# variables below so reduce redundant checks.
# However once all of FCM source is forked into LFRic Apps, we can condense these
# two variables into one under a single version of each, these being:
# PSYCLONE_PHYSICS_FILES, PSYCLONE_DIRECTORIES, PSYCLONE_PHYSICS_EXCEPTION
# File lists provided will use the transmute PSyclone method.
# https://code.metoffice.gov.uk/trac/lfric_apps/ticket/724


##### TRANSMUTE_INCLUDE_METHOD specify_include #####
# For CPU OMP, we want to choose which files get run through PSyclone,
# and preserve existing hand coded optimisations.

# Choose which files to Pre-proccess and PSyclone from physics / FCM UM source.

export PSYCLONE_PHYSICS_FILES_IMPORT =

export PSYCLONE_PHYSICS_FILES_FCM =

##### TRANSMUTE_INCLUDE_METHOD specify_include #####

##### TRANSMUTE_INCLUDE_METHOD specify_exclude #####
# For GPU, we may want to use more generic local.py transformation scripts and psyclone by directory.
# Advise which directories to pass to PSyclone.
# All files in these directories will be run through PSyclone using the transmute method.
# Also provide an optional exception list.
# These files will be filtered, and will NOT be run through PSyclone.

# Directories to psyclone
export PSYCLONE_DIRECTORIES_IMPORT =

export PSYCLONE_DIRECTORIES_FCM =

# A general file exception list
export PSYCLONE_PHYSICS_EXCEPTION_IMPORT =

export PSYCLONE_PHYSICS_EXCEPTION_FCM =

##### TRANSMUTE_INCLUDE_METHOD specify_exclude #####
