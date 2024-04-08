##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

include $(PROJECT_DIR)/build/fortran.mk
$(info UM physics specific compile options for $(FORTRAN_COMPILER) compiler)

include $(PROJECT_DIR)/build/fortran/$(FORTRAN_COMPILER).mk

science/%.o science/%.mod: private FFLAGS_EXTRA += $(FFLAGS_UM_PHYSICS)
