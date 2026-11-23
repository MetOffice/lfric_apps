##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Orchestrates the building of adjoint kernels and test algorithms.
#
# Expects the following macros to be defined:
# BUILD_ADJ_TESTS   - Set if the test algorithms are needed
# PATCH_DIR         - Directory containing pre and post PSyAd patches
# PSYAD_CONFIG_FILE - Configuration file for PSyAd
# PSYAD_WDIR        - Extra scratch space for this process
# WORKING_DIR       - Build's scratch space
##############################################################################
#
# Kernels which need adjointing...
#
KERNEL_LIST =

# Variables per kernel...
#
# Nada

do_psyad: $(if $(BUILD_ADJ_TESTS),$(WORKING_DIR)/driver/gen_adj_kernel_tests_mod.f90) \
          $(foreach kernel,$(KERNEL_LIST),$(WORKING_DIR)/$(kernel).adjoint)

$(WORKING_DIR)/driver/gen_adj_kernel_tests_mod.f90: | $(WORKING_DIR)/driver
	$(call MESSAGE,PSyAd test driver,$@)
	python $(ADJOINT_BUILD)/psyad_driver $(ADJOINT_BUILD)/gen_adj_kernel_tests_mod.f90.jinja \
                                         $@ $(notdir $(KERNEL_LIST))

%.f90.adjoint: %.f90 | $(PSYAD_WDIR)
	$(call MESSAGE,PSyAd,$<)
	$(ADJOINT_BUILD)/psyad_wrapper $(if $(BUILD_ADJ_TESTS),--algorithm-dir=$(PSYAD_WDIR)) \
                                   $(PSYAD_CONFIG_FILE) \
                                   $< \
                                   $(PATCH_DIR) \
                                   $(WORKING_DIR) \
                                   $(PSYAD_WDIR) \
                                   $(ACTIVE_$(patsubst %.f90.adjoint,%,$(notdir $*)))

$(PSYAD_WDIR) $(WORKING_DIR)/adjoint $(WORKING_DIR)/driver:
	$(call MESSAGE,Create,$@)
	mkdir -p $@
