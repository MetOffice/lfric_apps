##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Script to copy and pre-patch the kernels we wish
# to compile with psyad. The source files are
# $(PSYAD_FILES). The targets have the same
# structure as $(PSYAD_FILES), except the base
# directories are replaced with the PSyAD working directory
# $(PSYAD_WDIR). This is to simplify the rules in later scripts
# and allow developers to generate patches by editing the copies
# in $(PSYAD_WDIR) instead of the source files.

# The pre-patching consists of two phases. First,
# the source files are copied to their equivalent
# target. Then, if a patch exists, the target file
# is patched.

# Get the list of base directories where psyad kernels are located.
BASE_DIRS := $(shell python $(ADJOINT_BUILD)/print_paths.py -f $(PSYAD_FILES) -opt "base")

# Determine the target paths for the psyad kernels.
TARGETS := $(PSYAD_FILES)
$(foreach base_dir, \
          $(BASE_DIRS), \
          $(eval TARGETS := $(subst $(base_dir),$(PSYAD_WDIR),$(TARGETS)) \
          ) \
)
DIRECTORIES := $(addprefix $(PSYAD_WDIR)/,$(KERNEL_PATHS))

all: $(TARGETS)

# Use define block here to account for the subdirectories in base/source/kernel.
# $1 is kernel path, $2 is base directory.
define PRE_PATCH

#-----------------------------------------------------------------------------
# For F90s
#-----------------------------------------------------------------------------
# If a patch exists, copy source to target and patch the target.
$(PSYAD_WDIR)/$1/%.F90: $2/$1/%.F90 \
	$(PATCH_DIR)/kernel/%.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy source to target.
$(PSYAD_WDIR)/$1/%.F90: $2/$1/%.F90 \
	| $(DIRECTORIES)
	cp $$< $$@

#-----------------------------------------------------------------------------
# For f90s
#-----------------------------------------------------------------------------
# If a patch exists, copy source to target and patch the target.
$(PSYAD_WDIR)/$1/%.f90: $2/$1/%.f90 \
	$(PATCH_DIR)/kernel/%.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy source to target.
$(PSYAD_WDIR)/$1/%.f90: $2/$1/%.f90 \
	| $(DIRECTORIES)
	cp $$< $$@

endef

# Evaluating the PRE_PATCH definition
# for each base directory and kernel path
# within base directories. After this, all
# possible rules are now generated.
$(foreach base_dir, \
          $(BASE_DIRS), \
          $(foreach kernel_path, \
                    $(KERNEL_PATHS), \
                    $(eval $(call PRE_PATCH,$(kernel_path),$(base_dir))) \
          ) \
)

$(DIRECTORIES):
	mkdir -p $@

