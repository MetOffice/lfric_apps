##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(APPS_ROOT_DIR)/interfaces/socrates_interface/source

.PHONY: import-socrates_interface
import-socrates_interface:
    # Get a copy of the source code from the SCORATES repository
	python $(APPS_ROOT_DIR)/build/extract/extract_science.py -d $(APPS_ROOT_DIR)/dependencies.yaml -w $(WORKING_DIR) -e $(APPS_ROOT_DIR)/interfaces/socrates_interface/build/extract.yaml

    # Extract the interface code
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
			  SOURCE_DIR=$(PROJECT_SOURCE)
