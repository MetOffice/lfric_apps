##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/source

.PHONY: import-physics_schemes_interface

import-physics_schemes_interface:
    # Extract the interface code
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
		SOURCE_DIR=$(PROJECT_SOURCE)

    # Extract the physics schemes
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
		SOURCE_DIR=$(APPS_ROOT_DIR)/science/physics_schemes/source
