##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# These environment variables are set up:
#
# FC: Fortran compiler command
# FORTRAN_COMPILER: The compiler determined to be in use
# FFLAGS: Fortran compiler flags including "-I" module file locations
# LD: Linker command, probably the same as FC
# MPILD: MPI linker command
# LDFLAGS: Linker flags including "-L" library locations
#
###############################################################################

# The default value of FC is almost always "f77" which is of no use to us.
# An empty FC is also of no use.
ifneq "$(or $(filter default, $(origin FC)), $(filter x, x$(FC)))" ""
  $(error The FC environment variable must be set to a Fortran compiler command)
endif

# Make compiler macros available...
#
# Sometimes FC holds a full path which needs to be stripped off. It may also
# include a version number which also needs to go.
#
FORTRAN_COMPILER := $(firstword $(subst -, ,$(notdir $(FC))))

ifdef CRAY_ENVIRONMENT
  ifeq '$(PE_ENV)' 'CRAY'
    FORTRAN_COMPILER = crayftn
  else ifeq '$(PE_ENV)' 'INTEL'
    FORTRAN_COMPILER = ifort
  else ifeq '$(PE_ENV)' 'GNU'
    FORTRAN_COMPILER = gfortran
  else ifeq '$(PE_ENV)' 'PGI'
    FORTRAN_COMPILER = pgfortran
  else
    $(error Unrecognised Cray programming environment)
  endif
endif

include $(PROJECT_DIR)/build/fortran/$(FORTRAN_COMPILER).mk
export FORTRAN_COMPILER F_MOD_DESTINATION_ARG

FFLAGS += $(FFLAGS_COMPILER)
FFLAGS  := $(FFLAGS)
LDFLAGS := $(LDFLAGS)
export FFLAGS LDFLAGS
