##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
PROJECT_SOURCE = $(APPS_ROOT_DIR)/interfaces/socrates_interface/source

.PHONY: import-socrates_interface
import-socrates_interface:
    # Get a copy of the source code from the SOCRATES repository
	python $(APPS_ROOT_DIR)/build/extract/extract_science.py \
        -d $(APPS_ROOT_DIR)/dependencies.yaml \
        -w $(WORKING_DIR)/../foo \
        -e $(APPS_ROOT_DIR)/interfaces/socrates_interface/build/extract.yaml
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
        SOURCE_DIR=$(WORKING_DIR)/../foo/socrates/src \
        WORKING_DIR=$(WORKING_DIR)/socrates

    # Extract the interface code
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(PROJECT_SOURCE)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(PROJECT_SOURCE) \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
