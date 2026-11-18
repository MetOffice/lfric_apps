##############################################################################
# (c) Crown copyright, Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PRE_PROCESS_MACROS += LFRIC UM_PHYSICS USSPPREC_32B LSPREC_32B

.PHONY: import-physics_schemes
import-physics_schemes: project_dir = $(APPS_ROOT_DIR)/science/physics_schemes
import-physics_schemes:
	# Extract external phsysics schemes (UKCA and CASIM)
	python $(APPS_ROOT_DIR)/build/extract/extract_science.py \
	       -d $(APPS_ROOT_DIR)/dependencies.yaml \
	       -w $(WORKING_DIR)/../physbits \
	       -e $(project_dir)/build/extract.yaml
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(WORKING_DIR)/../physbits/casim/src \
	          WORKING_DIR=$(WORKING_DIR)/casim
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(WORKING_DIR)/../physbits/ukca/src \
	          WORKING_DIR=$(WORKING_DIR)/ukca

	# Extract the physics schemes
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(project_dir)/source
