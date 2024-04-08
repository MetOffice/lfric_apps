##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# NOTE: Import of compile options from LFRic infrastructure is temporarily
# suspended here as a workaround for #2340 in which application of the
# -qoverride-limits option was preventing compilation of a UKCA module.
# include $(LFRIC_BUILD)/compile_options.mk

$(info UM physics specific compile options)

include $(PROJECT_DIR)/build/fortran/$(FORTRAN_COMPILER).mk

science/%.o science/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
%/limited_area_constants_mod.o: private FFLAGS_EXTRA = $(FFLAGS_INTEL_FIX_ARG)

$(info Disable warnings-turned-error caused by undeclared external functions - see ifort.mk)
%mpi_mod.o %mpi_mod.mod: private FFLAGS_EXTRA = $(FFLAGS_INTEL_EXTERNALS)
science/src/socrates/src/radiance_core/%.o: private FFLAGS_EXTRA = $(FFLAGS_INTEL_EXTERNALS)
science/src/socrates/src/interface_core/%.o: private FFLAGS_EXTRA = $(FFLAGS_INTEL_EXTERNALS)
