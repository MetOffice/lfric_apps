##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# Default for file selection method for transformation files is for CPU OMP currently
# This should be overwritten in:
# rose-stem/site/<site>/common/suite_config_<target>.cylc
#
TRANSMUTE_INCLUDE_METHOD ?= specify_include

# Set the DSL Method in use to collect the correct transformation files.
DSL := transmute

# Find the specific files we wish to pre-processed and PSyclone from physics source
# Set our target dependency to the version of the file we are to generate after
# the psycloning step.
#
# For CPU OMP method, we want specific files.
SOURCE_F_FILES_PASS := $(foreach THE_FILE, $(PSYCLONE_PASS_NO_SCRIPT), $(patsubst $(SOURCE_DIR)/%.xu90, $(SOURCE_DIR)/%.F90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).xu90' -print)))

# Default make target for file
#
.PHONY: psyclone_pass

# Call this target, expect these files to be done first
#
psyclone_pass: $(SOURCE_F_FILES_PASS)

# PSyclone files back into F90 files.
# Technically as they are pre-processed and PSyclone-d, they should be f90 files,
# however the analysis step broke and so needs some work. For now they will be F90.

# Where no optimisation script exists, don't use it.
#
$(SOURCE_DIR)/%.F90: $(SOURCE_DIR)/%.xu90
	echo PSyclone pass with no optimisation applied, OMP and Clauses removed on $<
	PYTHONPATH=$(LFRIC_BUILD)/psyclone:$(abspath ../../interfaces/physics_schemes_interface/build):$$PYTHONPATH psyclone \
			-l all \
			-o $(SOURCE_DIR)/$*.F90 \
			$<
