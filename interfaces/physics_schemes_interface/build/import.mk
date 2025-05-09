##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/source

.PHONY: import-physics_schemes_interface

import-physics_schemes_interface:
    # Extract and Psyclone the interface code
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
		SOURCE_DIR=$(PROJECT_SOURCE)

	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
		SOURCE_DIR=$(PROJECT_SOURCE) \
		OPTIMISATION_PATH=$(OPTIMISATION_PATH)
    # Extract and Psyclone the physics schemes
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
		SOURCE_DIR=$(APPS_ROOT_DIR)/science/physics_schemes/source

	# Pre-process and then PSyclone the same physics source
	# Note depending on site target requirements, either all source or specific source will be affected
	# See application overrides and files exported there
	#
	echo "Import.mk pre-processing and psyclone"

	$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/build/pre_process_phys.mk \
		SOURCE_DIR=$(WORKING_DIR) \
		PSYCLONE_PHYSICS_FILES="$(PSYCLONE_PHYSICS_FILES_IMPORT)" \
		PSYCLONE_DIRECTORIES="$(PSYCLONE_DIRECTORIES_IMPORT)" \
		PSYCLONE_PHYSICS_EXCEPTION="$(PSYCLONE_PHYSICS_EXCEPTION_IMPORT)"

	$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/build/psyclone_transmute.mk \
		SOURCE_DIR=$(WORKING_DIR) \
		OPTIMISATION_PATH_PSY=$(APPS_ROOT_DIR)/applications/$(PROJECT)/$(OPTIMISATION_PATH) \
		PSYCLONE_PHYSICS_FILES="$(PSYCLONE_PHYSICS_FILES_IMPORT)" \
		PSYCLONE_DIRECTORIES="$(PSYCLONE_DIRECTORIES_IMPORT)" \
		PSYCLONE_PHYSICS_EXCEPTION="$(PSYCLONE_PHYSICS_EXCEPTION_IMPORT)"
