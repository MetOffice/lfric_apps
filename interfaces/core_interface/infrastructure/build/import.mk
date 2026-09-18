##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export IGNORE_DEPENDENCIES += netcdf mpi mpi_f08 yaxt mod_oasis

# Sorry about the reverse logic: "if not defined no_mpi" really means: "if mpi"
ifndef NO_MPI
 export EXTERNAL_DYNAMIC_LIBRARIES += yaxt yaxt_c
endif
export EXTERNAL_DYNAMIC_LIBRARIES += netcdff netcdf

ifdef USE_VERNIER
  export IGNORE_DEPENDENCIES += vernier_mod
  export EXTERNAL_DYNAMIC_LIBRARIES += vernier_f vernier_c vernier
  $(info Linking Vernier)
endif

TEMPLATE_TOOL = $(LFRIC_BUILD)/tools/Templaterator
TYPE_TABLE_real32 = real
TYPE_TABLE_real64 = real
TYPE_TABLE_int32 = integer

.PHONY: import-infrastructure
import-infrastructure: $(WORKING_DIR)/field/field_real32_mod.f90 \
                       $(WORKING_DIR)/field/field_real64_mod.f90 \
                       $(WORKING_DIR)/field/field_int32_mod.f90 \
                       $(WORKING_DIR)/operator/operator_real32_mod.f90 \
                       $(WORKING_DIR)/operator/operator_real64_mod.f90 \
                       $(WORKING_DIR)/scalar/scalar_real32_mod.f90 \
                       $(WORKING_DIR)/scalar/scalar_real64_mod.f90 \
                       $(WORKING_DIR)/scalar/scalar_int32_mod.f90
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source
	$Q$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/interfaces/build/psyclone_psykal.mk \
	          SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)

$(WORKING_DIR)/field/field_%_mod.f90: $(LFRIC_INFRASTRUCTURE)/source/field/field_mod.t90 \
                                      | $(WORKING_DIR)/field
	$(call MESSAGE,Templating, $<)
	$Q$(TEMPLATE_TOOL) $< -o $@ -s type=$(TYPE_TABLE_$*) -s kind=$*

$(WORKING_DIR)/field:
	$Qmkdir -p $@

$(WORKING_DIR)/operator/operator_%_mod.f90: $(LFRIC_INFRASTRUCTURE)/source/operator/operator_mod.t90 \
                                            | $(WORKING_DIR)/operator
	$(call MESSAGE,Templating, $<)
	$Q$(TEMPLATE_TOOL) $< -o $@ -s kind=$*

$(WORKING_DIR)/operator:
	$Qmkdir -p $@

$(WORKING_DIR)/scalar/scalar_%_mod.f90: $(LFRIC_INFRASTRUCTURE)/source/scalar/scalar_mod.t90 \
                                      | $(WORKING_DIR)/scalar
	$(call MESSAGE,Templating, $<)
	$Q$(TEMPLATE_TOOL) $< -o $@ -s type=$(TYPE_TABLE_$*) -s kind=$*

$(WORKING_DIR)/scalar:
	$Qmkdir -p $@
