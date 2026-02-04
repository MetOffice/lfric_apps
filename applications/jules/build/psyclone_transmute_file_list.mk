##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# File lists provided will use the transmute PSyclone method.
# https://code.metoffice.gov.uk/trac/lfric_apps/ticket/724


##### TRANSMUTE_INCLUDE_METHOD specify_include #####
# For CPU OMP, we want to choose which files get run through PSyclone,
# and preserve existing hand coded optimisations.

# Choose which files to Pre-proccess and PSyclone from physics_schemes / other physics source (e.g. UKCA)

export PSYCLONE_PHYSICS_FILES =

##### TRANSMUTE_INCLUDE_METHOD specify_include #####

# List to use PSyclone explicitly without any opt script
# This will remove hand written (OMP) directives in the source
# Used by both methods, specify_include and specify_exclude
export PSYCLONE_PASS_NO_SCRIPT =

##### TRANSMUTE_INCLUDE_METHOD specify_exclude #####
# For GPU, we may want to use more generic local.py transformation scripts and psyclone by directory.
# Advise which directories to pass to PSyclone.
# All files in these directories will be run through PSyclone using the transmute method.
# Also provide an optional exception list.
# These files will be filtered, and will NOT be run through PSyclone.

# Directories to psyclone
export PSYCLONE_DIRECTORIES =

# A general file exception list
export PSYCLONE_PHYSICS_EXCEPTION =

##### TRANSMUTE_INCLUDE_METHOD specify_exclude #####
