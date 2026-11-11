##############################################################################
# (c) Crown copyright, Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
physics_schemes_source = $(APPS_ROOT_DIR)/science/physics_schemes/source

.PHONY: import-physics_schemes
import-physics_schemes:
	# Extract the physics schemes
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
		SOURCE_DIR=$(physics_schemes_source)
