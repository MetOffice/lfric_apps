##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# PRE-PATCH STAGE RULES
##############################################################################
# Use define block here to account for the subdirectories in base/source/kernel.
# $1 is base directory, $2 is kernel path.
define PRE_PATCH

#-----------------------------------------------------------------------------
# For F90s
#-----------------------------------------------------------------------------
# If a patch exists, copy source to target and patch the target.
$(PSYAD_WDIR)/$2/%.F90: $1/$2/%.F90 \
	$(PATCH_DIR)/kernel/%.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy source to target.
$(PSYAD_WDIR)/$2/%.F90: $1/$2/%.F90 \
	| $(DIRECTORIES)
	cp $$< $$@

#-----------------------------------------------------------------------------
# For f90s
#-----------------------------------------------------------------------------
# If a patch exists, copy source to target and patch the target.
$(PSYAD_WDIR)/$2/%.f90: $1/$2/%.f90 \
	$(PATCH_DIR)/kernel/%.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy source to target.
$(PSYAD_WDIR)/$2/%.f90: $1/$2/%.f90 \
	| $(DIRECTORIES)
	cp $$< $$@

endef # PRE_PATCH

# Evaluating the PRE_PATCH definition for each base directory 
# and kernel path within base directories.
# After this, all possible rules are now generated.
$(foreach base_dir, \
          $(BASE_DIRS), \
          $(foreach kernel_path, \
                    $(KERNEL_PATHS), \
                    $(eval $(call PRE_PATCH,$(base_dir),$(kernel_path))) \
          ) \
)
