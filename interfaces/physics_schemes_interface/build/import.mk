##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
.PHONY: import-physics_schemes_interface
import-physics_schemes_interface: project_dir = $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface
import-physics_schemes_interface:
	# Extract the interface code
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(project_dir)/source
