##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PRE_PROCESS_MACROS += LFRIC UM_JULES

.PHONY: import-jules_interface
import-jules_interface: project_dir = $(APPS_ROOT_DIR)/interfaces/jules_interface
import-jules_interface:
	# Get a copy of the source code from the JULES repository
	$Qpython $(APPS_ROOT_DIR)/build/extract/extract_science.py \
	         -d $(APPS_ROOT_DIR)/dependencies.yaml \
	         -w $(WORKING_DIR)/../physbits \
	         -e $(project_dir)/build/extract.yaml
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(WORKING_DIR)/../physbits/jules/src \
	          WORKING_DIR=$(WORKING_DIR)/jules

	# Extract the interface code
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(project_dir)/source
