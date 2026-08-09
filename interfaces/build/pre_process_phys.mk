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

# Build a set of "-D" argument for any pre-processor macros from core
#
MACRO_ARGS := $(addprefix -D,$(PRE_PROCESS_MACROS))

# Find the specific files we wish to transmute psyclone from physics source
# Set our target dependency to the version of the file we are to generate after
# The preprocessing step. #
# .xu90 files are to represent preprocessed source, bound for transmute psyclone,
# but are not psykal files, denoted by .x90
#
ifeq ("$(TRANSMUTE_INCLUDE_METHOD)", "specify_include")
# For CPU OMP method, we want specific files
	SOURCE_xu_FILES := $(foreach THE_FILE, $(PSYCLONE_PHYSICS_FILES), $(patsubst $(SOURCE_DIR)/%.f90, $(SOURCE_DIR)/%.xu90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).f90' -print)))
	SOURCE_xu_FILES += $(foreach THE_FILE, $(PSYCLONE_PASS_NO_SCRIPT), $(patsubst $(SOURCE_DIR)/%.f90, $(SOURCE_DIR)/%.xu90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).f90' -print)))
else ifeq ("$(TRANSMUTE_INCLUDE_METHOD)", "specify_exclude")
# For the offload method, we want to filter out specific files, and psyclone the rest
# We don't want to wildcard the whole working directory, this will cause problems. 
# We want to specifically choose directories we want to pass to the psyclone transmute method
# Therefore if nothing is present in the PSYCLONE_DIRECTORIES variable, then nothing will be 
# pre-processed.
	ifneq ($(strip $(PSYCLONE_DIRECTORIES)),)
		EXTEND_DIR_FULL_PATH := $(foreach THE_DIRECTORY, $(PSYCLONE_DIRECTORIES), $(shell find $(SOURCE_DIR) -name $(THE_DIRECTORY) -print))
		SOURCE_xu_FILES_FULL := $(strip $(foreach THE_PSY_DIR, $(EXTEND_DIR_FULL_PATH), $(patsubst $(SOURCE_DIR)/%.f90, $(SOURCE_DIR)/%.xu90, $(shell find $(THE_PSY_DIR) -name '*.f90' -print))))
		SOURCE_xu_EXCEPTION := $(strip $(foreach THE_FILE, $(PSYCLONE_PHYSICS_EXCEPTION), $(patsubst $(SOURCE_DIR)/%.f90, $(SOURCE_DIR)/%.xu90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).f90' -print))))
		SOURCE_xu_FILES := $(filter-out $(SOURCE_xu_EXCEPTION), $(SOURCE_xu_FILES_FULL))
	endif
endif


# Default make target for file
# 
.PHONY: pre_process

include $(LFRIC_BUILD)/fortran.mk

# Call this target, expect these files to be done first
#
pre_process: $(SOURCE_xu_FILES)

# Move the pre-processed source to an .xu90 file
#
$(SOURCE_DIR)/%.xu90: $(SOURCE_DIR)/%.f90
	echo Moving $< to xu90 for Transmute
	-mv $(SOURCE_DIR)/$*.f90 $@
