#############################################################################
#(c) Crown copyright 2024 Met Office. All rights reserved.
#The file LICENCE, distributed with this code, contains details of the terms
#under which the code may be used.
#############################################################################
# Various things specific to the NVIDIA Fortran compiler.
#############################################################################

$(info Project specials for NVIDIA compiler)

export FFLAGS_UM_PHYSICS = -r8

# The lfric_atm app defines an extra set of debug flags for
# fast-debug. For this compiler use the same as the full-debug
# settings
FFLAGS_FASTD_INIT         = $(FFLAGS_INIT) 
FFLAGS_FASTD_RUNTIME      = $(FFLAGS_RUNTIME)
