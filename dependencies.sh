##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# The variables in this file act as overrides to the contents
# of the UM's fcm-make configuration files in order to extract
# and preprocess the required code for use in LFRic.

# Revisions and sources for dependent repositories which
# are extracted as part of the build.  A blank *_sources
# variable will result in extraction from the project trunk.

# Note that on the LFRic trunk the *_sources variables should
# always be blank (but can be used on branches during development
# to enable testing of changes)

# On commit to trunk the code reviewer should ensure all *_sources
# variables below are empty, and the revisions of any projects with
# dependent changes should be updated to the revision at which those
# changes were committed to the project's trunk

export lfric_core_rev=54094
export lfric_core_sources=
export casim_rev=apps2.2
export casim_sources=
export jules_rev=31238
export jules_sources=
export socrates_rev=1869
export socrates_sources=
export ukca_rev=7497
export ukca_sources=
