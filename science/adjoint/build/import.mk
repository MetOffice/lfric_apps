##############################################################################
# (c) Crown copyright 2023-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
PROJECT_SOURCE     = $(APPS_ROOT_DIR)/science/adjoint/source
PSYAD_CONFIG_FILE ?= $(CORE_ROOT_DIR)/etc/psyclone.cfg

.PHONY: import-adjoint
import-adjoint: export ADJOINT_BUILD   := $(APPS_ROOT_DIR)/science/adjoint/build
import-adjoint: export PSYAD_WDIR      := $(abspath $(WORKING_DIR)/../psyad)
import-adjoint: export PATCH_DIR       := $(APPS_ROOT_DIR)/science/adjoint/patches
import-adjoint:
	# Building PSyAD kernels and driver.
	#
	$Q$(MAKE) $(QUIET_ARG) -f $(ADJOINT_BUILD)/build_psyad.mk

ifeq "$(BUILD_ADJ_TESTS)" "TRUE"
	# Need this step to PSyclone the generated adjoint tests from PSyAD
	#
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(PSYAD_WDIR)/outgoing \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH) \
	          PSYCLONE_CONFIG_FILE=$(PSYAD_CONFIG_FILE)
endif

	# Standard import commands
	#
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(PROJECT_SOURCE)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(PROJECT_SOURCE) \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
