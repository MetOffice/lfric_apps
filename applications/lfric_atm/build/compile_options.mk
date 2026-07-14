##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

include $(PROJECT_DIR)/build/fortran.mk
$(info UM physics specific compile options for $(FORTRAN_COMPILER) compiler)

include $(PROJECT_DIR)/build/fortran/$(FORTRAN_COMPILER).mk

casim/%.o science/%.mod: private FFLAGS_EXTRA += $(FFLAGS_UM_PHYSICS)
ukca/%.o science/%.mod: private FFLAGS_EXTRA += $(FFLAGS_UM_PHYSICS)
jules/%.o jules/%.mod: private FFLAGS_EXTRA += $(FFLAGS_UM_PHYSICS)
socrates/%.o socrates/%.mod: private FFLAGS_EXTRA += $(FFLAGS_UM_PHYSICS)
legacy/%.o legacy/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
AC_assimilation/%.o AC_assimilation/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
aerosols/%.o aerosols/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
atmosphere_service/%.o atmosphere_service/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
boundary_layer/%.o boundary_layer/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
carbon/%.o carbon/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
convection/%.o convection/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
diffusion_and_filtering/%.o diffusion_and_filtering/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
dynamics/%.o dynamics/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
dynamics_advection/%.o dynamics_advection/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
electric/%.o electric/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
free_tracers/%.o free_tracers/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
gravity_wave_drag/%.o gravity_wave_drag/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
idealised/%.o idealised/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
large_scale_cloud/%.o large_scale_cloud/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
large_scale_precipitation/%.o large_scale_precipitation/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
PWS_diagnostics/%.o PWS_diagnostics/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
radiation_control/%.o radiation_control/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
stochastic_physics/%.o stochastic_physics/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
tracer_advection/%.o tracer_advection/%.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS)
