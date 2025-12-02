##############################################################################
# (c) Crown copyright 2023-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(APPS_ROOT_DIR)/science/linear/source

.PHONY: import-linear
import-linear:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(PROJECT_SOURCE)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(PROJECT_SOURCE) \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
