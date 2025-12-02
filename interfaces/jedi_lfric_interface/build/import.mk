##############################################################################
# (c) Crown copyright 2022-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(APPS_ROOT_DIR)/interfaces/jedi_lfric_interface/source

.PHONY: import-jedi_lfric_interface
import-jedi_lfric_interface:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(PROJECT_SOURCE)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
            SOURCE_DIR=$(PROJECT_SOURCE) \
            OPTIMISATION_PATH=$(OPTIMISATION_PATH)