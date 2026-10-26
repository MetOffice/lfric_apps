##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Some of the content of this file has been produced with the assistance of
# Met Office Github Copilot Enterprise."

# Default for file selection method for transformation files is for CPU OMP currently
# This should be overwritten in:
# rose-stem/site/<site>/common/suite_config_<target>.cylc
#
TRANSMUTE_INCLUDE_METHOD ?= specify_include

# Set the DSL Method in use to collect the correct transformation files.
DSL := transmute

# Set default PSyclone transmute command additional options
PSYCLONE_TRANSMUTE_EXTRAS ?= -l all
#
# The command used to invoke PSyclone. By default this is the persistent
# server client (psyclone_client.py) which keeps a single PSyclone instance
# resident and dispatches jobs to a pool of pre-forked workers, avoiding the
# repeated cost of loading the Python libraries from disk. It is a drop-in
# replacement for the "psyclone" binary and falls back to it automatically if
# the server is unavailable. Override PSYCLONE=psyclone to bypass the server.
#
# The server's lifetime is pinned to the owning make process: lfric.mk exports
# PSYCLONE_OWNER_PID (the pid of the top-level make) and the server exits as
# soon as that process does, so no server survives the build that started it.
# Should that variable be unset the client determines the owner itself, by
# finding the outermost make process in its own ancestry.
PSYCLONE_MODE ?= standard
ifeq ($(PSYCLONE_MODE),server)
    PSYCLONE = $(LFRIC_BUILD)/psyclone/psyclone_client.py
else
	PSYCLONE = psyclone
endif

# Find the specific files we wish to pre-processed and PSyclone from physics source
# Set our target dependency to the version of the file we are to generate after
# the psycloning step.
#
ifeq ("$(TRANSMUTE_INCLUDE_METHOD)", "specify_include")
# For CPU OMP method, we want specific files.
	SOURCE_F_FILES_ALL := $(foreach THE_FILE, $(PSYCLONE_PHYSICS_FILES), $(patsubst $(SOURCE_DIR)/%.xu90, $(SOURCE_DIR)/%.f90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).xu90' -print)))
	SOURCE_F_FILES_PASS := $(foreach THE_FILE, $(PSYCLONE_PASS_NO_SCRIPT), $(patsubst $(SOURCE_DIR)/%.xu90, $(SOURCE_DIR)/%.f90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).xu90' -print)))
	SOURCE_F_FILES := $(filter-out $(SOURCE_F_FILES_PASS), $(SOURCE_F_FILES_ALL))
else ifeq ("$(TRANSMUTE_INCLUDE_METHOD)", "specify_exclude")
# For the offload method, we want to filter out specific files, and do the rest.
# We don't want to wildcard the whole working directory, this will cause problems.
# We want to specifically choose directories we want to pass to the PSyclone transmute method.
# Therefore if nothing is present in the PSYCLONE_DIRECTORIES variable, then nothing will be
# passed to PSyclone.
	ifneq ($(strip $(PSYCLONE_DIRECTORIES)),)
		EXTEND_DIR_FULL_PATH := $(foreach THE_DIRECTORY, $(PSYCLONE_DIRECTORIES), $(shell find $(SOURCE_DIR) -name $(THE_DIRECTORY) -print))
		SOURCE_F_FILES_FULL := $(strip $(foreach THE_PSY_DIR, $(EXTEND_DIR_FULL_PATH), $(patsubst $(SOURCE_DIR)/%.xu90, $(SOURCE_DIR)/%.f90, $(shell find $(THE_PSY_DIR) -name '*.xu90' -print))))
		SOURCE_EXCEPTION := $(strip $(foreach THE_FILE, $(PSYCLONE_PHYSICS_EXCEPTION), $(patsubst $(SOURCE_DIR)/%.xu90, $(SOURCE_DIR)/%.f90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).xu90' -print))))
		SOURCE_F_FILES_ALL := $(filter-out $(SOURCE_EXCEPTION), $(SOURCE_F_FILES_FULL))
		SOURCE_F_FILES_PASS := $(foreach THE_FILE, $(PSYCLONE_PASS_NO_SCRIPT), $(patsubst $(SOURCE_DIR)/%.xu90, $(SOURCE_DIR)/%.f90, $(shell find $(SOURCE_DIR) -name '$(THE_FILE).xu90' -print)))
		SOURCE_F_FILES := $(filter-out $(SOURCE_F_FILES_PASS), $(SOURCE_F_FILES_ALL))
	endif
endif


# Default make target for file
#
.PHONY: psyclone

# Call this target, expect these files to be done first
#
psyclone: $(SOURCE_F_FILES)

# PSyclone files back into f90 files.

# Where an optimisation script exists for a specific file, use it.
#
$(SOURCE_DIR)/%.f90: $(SOURCE_DIR)/%.xu90 $(OPTIMISATION_PATH)/$(DSL)/%.py
	echo PSyclone with file override script $(OPTIMISATION_PATH_PSY)/$(DSL)/$*.py on $<
	PYTHONPATH=$(LFRIC_BUILD)/psyclone:$(abspath $(OPTIMISATION_PATH)/$(DSL)):$(abspath ../../interfaces/build):$$PYTHONPATH $(PSYCLONE) \
			-s $(OPTIMISATION_PATH_PSY)/$(DSL)/$*.py \
			-o $(SOURCE_DIR)/$*.f90 \
			$(PSYCLONE_TRANSMUTE_EXTRAS) \
			$<

# Where a local optimisation script exists, use it.
#
.SECONDEXPANSION:
$(SOURCE_DIR)/%.f90: $(SOURCE_DIR)/%.xu90 $$(dir $$(OPTIMISATION_PATH_PSY)/$$(DSL)/$$*)local.py
	echo PSyclone with local script $(dir $(OPTIMISATION_PATH_PSY)/$(DSL)/$*)local.py on $<
	PYTHONPATH=$(LFRIC_BUILD)/psyclone:$(abspath $(OPTIMISATION_PATH)/$(DSL)):$(abspath ../../interfaces/build):$$PYTHONPATH $(PSYCLONE) \
			-s $(dir $(OPTIMISATION_PATH_PSY)/$(DSL)/$*)local.py \
			-o $(SOURCE_DIR)/$*.f90 \
			$(PSYCLONE_TRANSMUTE_EXTRAS) \
			$<

# Where a global optimisation script exists, use it.
#
$(SOURCE_DIR)/%.f90: $(SOURCE_DIR)/%.xu90 $(OPTIMISATION_PATH)/$(DSL)/global.py
	echo PSyclone with global script $(OPTIMISATION_PATH_PSY)/$(DSL)/global.py on $<
	PYTHONPATH=$(LFRIC_BUILD)/psyclone:$(abspath $(OPTIMISATION_PATH)/$(DSL)):$(abspath ../../interfaces/build):$$PYTHONPATH $(PSYCLONE) \
			-s $(OPTIMISATION_PATH_PSY)/$(DSL)/global.py \
			-o $(SOURCE_DIR)/$*.f90 \
			$(PSYCLONE_TRANSMUTE_EXTRAS) \
			$<
