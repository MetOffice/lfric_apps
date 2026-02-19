##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# NOTE: Import of compile options from LFRic infrastructure is temporarily
# suspended here as a workaround for #2340 in which application of the
# -qoverride-limits option was preventing compilation of a UKCA module.
# include $(LFRIC_BUILD)/compile_options.mk

$(info UM physics specific compile options)

include $(PROJECT_DIR)/build/fortran/$(FORTRAN_COMPILER).mk

casim/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
ukca/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
jules/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
socrates/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
legacy/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
AC_assimilation/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
aerosols/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
atmosphere_service/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
boundary_layer/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
carbon/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
convection/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
diffusion_and_filtering/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
dynamics/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
dynamics_advection/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
electric/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
free_tracers/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
gravity_wave_drag/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
idealised/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
large_scale_cloud/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
large_scale_precipitation/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
PWS_diagnostics/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
radiation_control/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
stochastic_physics/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
tracer_advection/%.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
%/limited_area_constants_mod.o: private FFLAGS_EXTRA = $(FFLAGS_INTEL_FIX_ARG)

$(info Disable warnings-turned-error caused by undeclared external functions - see ifort.mk)
%mpi_mod.o: private FFLAGS_EXTRA = $(FFLAGS_INTEL_EXTERNALS)
socrates/src/radiance_core/%.o: private FFLAGS_EXTRA = $(FFLAGS_INTEL_EXTERNALS)
socrates/src/interface_core/%.o: private FFLAGS_EXTRA = $(FFLAGS_INTEL_EXTERNALS)
