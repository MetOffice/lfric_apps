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

export lfric_core_rev=51572
export lfric_core_sources=

export casim_rev=um13.7
export casim_sources=
export jules_rev=um13.7
export jules_sources=
export shumlib_rev=um13.7
export shumlib_sources=
export socrates_rev=um13.7
export socrates_sources=
export ukca_rev=um13.7
export ukca_sources=
export um_rev=vn13.7
export um_sources=

#### Do not edit the definitions below this line without
#### consulting the owners of this build system

# The revision and location of the build configuration files
export config_revision=@$um_rev
export config_root_path=fcm:um.xm_tr
export config_type=atmos

# Set to extract then preprocess the source only (this is
# disabling the build of the atmos/recon executables)
export extract=extract
export compile_recon=preprocess-recon

export compile_atmos=preprocess-atmos

# Define any preprocessor macros that are required
# (in addition to the ones normally used by the UM)
export keys_atmos_extra='LFRIC=lfric USSPPREC_32B=true LSPREC_32B=true'
