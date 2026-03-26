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

# Find the specific files we wish to pre process and psyclone from physics source
# Set our target dependency to the version of the file we are to generate after
# The preprocessing step. #
# .xu90 files are to represent preprocessed source, bound for psyclone,
# but are not psykal files, denoted by .x90
#
ifeq ("$(TRANSMUTE_INCLUDE_METHOD)", "specify_include")
# For CPU OMP method, we want specific files
	SOURCE_xu_FILES := $(foreach THE_FILE, $(PSYCLONE_PHYSICS_FILES), $(patsubst $(SOURCE_DIR)/%.F90, $(SOURCE_DIR)/%.xu90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).F90' -print)))
	SOURCE_xu_FILES += $(foreach THE_FILE, $(PSYCLONE_PASS_NO_SCRIPT), $(patsubst $(SOURCE_DIR)/%.F90, $(SOURCE_DIR)/%.xu90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).F90' -print)))
else ifeq ("$(TRANSMUTE_INCLUDE_METHOD)", "specify_exclude")
# For the offload method, we want to filter out specific files, and psyclone the rest
# We don't want to wildcard the whole working directory, this will cause problems. 
# We want to specifically choose directories we want to pass to the psyclone transmute method
# Therefore if nothing is present in the PSYCLONE_DIRECTORIES variable, then nothing will be 
# pre-processed.
	ifneq ($(strip $(PSYCLONE_DIRECTORIES)),)
		EXTEND_DIR_FULL_PATH := $(foreach THE_DIRECTORY, $(PSYCLONE_DIRECTORIES), $(shell find $(SOURCE_DIR) -name $(THE_DIRECTORY) -print))
		SOURCE_xu_FILES_FULL := $(strip $(foreach THE_PSY_DIR, $(EXTEND_DIR_FULL_PATH), $(patsubst $(SOURCE_DIR)/%.F90, $(SOURCE_DIR)/%.xu90, $(shell find $(THE_PSY_DIR) -name '*.F90' -print))))
		SOURCE_xu_EXCEPTION := $(strip $(foreach THE_FILE, $(PSYCLONE_PHYSICS_EXCEPTION), $(patsubst $(SOURCE_DIR)/%.F90, $(SOURCE_DIR)/%.xu90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).F90' -print))))
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

# Make a copy of target file,
# Preprocess target file,
# Remove original F90, see psyclone step
#
# For the nvidia compiler, they only output into f90,
# we need to move any f90 files to xu90 files for psyclone.
# It also seems to place them at the root of the working dir.
# See ticket Apps#624 for further context
#
ifeq ("$(FORTRAN_COMPILER)", "nvfortran")
$(SOURCE_DIR)/%.xu90: $(SOURCE_DIR)/%.F90
	echo Pre processing $<
	$(FPP) $(FPPFLAGS) $(MACRO_ARGS) -o $@ $<
	-mv $(SOURCE_DIR)/$*.f90 $@
	-mv $(shell basename $*.f90) $@
	rm $<
else
$(SOURCE_DIR)/%.xu90: $(SOURCE_DIR)/%.F90
	echo Pre processing $<
	$(FPP) $(FPPFLAGS) $(MACRO_ARGS) $< >$@
	rm $<
endif
