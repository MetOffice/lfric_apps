##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Builds the PSyAD kernels using three stages.
# 1) Pre-patch: copies (and patches) tangent linear kernels from their base dir to PSYAD_WDIR.
# 2) PSyAD: generates adjoint kernel and adjoint test algorithm from the pre-patch stage tangent linear kernels.
# 3) Post-patch: copies (and patches) adjoint kernels and adjoint test algorithms from PSYAD_WDIR to WORKING_DIR.

# List of paths to files we wish to compile using PSyAD.
PSYAD_FILES_CORE := $(shell cat $(ADJOINT_BUILD)/psyad_files_list_core.txt)
PSYAD_FILES_CORE := $(addprefix $(CORE_ROOT_DIR)/,$(PSYAD_FILES_CORE))
PSYAD_FILES_APPS := $(shell cat $(ADJOINT_BUILD)/psyad_files_list_apps.txt)
PSYAD_FILES_APPS := $(addprefix $(APPS_ROOT_DIR)/,$(PSYAD_FILES_APPS))
all: export PSYAD_FILES := $(PSYAD_FILES_CORE) $(PSYAD_FILES_APPS)

# List of all kernel related subdirectories in the base directories where PSyAD files are located.
all: export KERNEL_PATHS := $(shell python $(ADJOINT_BUILD)/print_paths.py -f $(PSYAD_FILES) -opt "kernel")
all: export ALG_PATHS   := $(subst kernel,algorithm,$(KERNEL_PATHS))

# Runs the pre-patch, PSyAD compilation and post-patch scripts.
all:
	$Q$(MAKE) $(QUIET_ARG) -f $(ADJOINT_BUILD)/pre_patch.mk
	$Q$(MAKE) $(QUIET_ARG) -f $(ADJOINT_BUILD)/psyad.mk
	$Q$(MAKE) $(QUIET_ARG) -f $(ADJOINT_BUILD)/post_patch.mk
