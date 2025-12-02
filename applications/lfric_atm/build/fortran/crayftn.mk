##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the Cray Fortran compiler.
##############################################################################

$(info Project specials for Cray compiler)

export FFLAGS_UM_PHYSICS = -s real64

include $(PROJECT_DIR)/build/fortran/crayftn/$(PROFILE).mk