##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# PSyAD STAGE RULES
##############################################################################

# Active variables for PSyAD compilation.
include $(ADJOINT_BUILD)/psyad_vars.mk

define PSYAD
# Use define block here to account for subdirectories.
# $2 is kernel path, $3 is algorithm path.
#-----------------------------------------------------------------------------
# For F90s
#-----------------------------------------------------------------------------
$(PSYAD_WDIR)/$2/atlt_%_alg_mod.X90: $(PSYAD_WDIR)/$1/atl_%_kernel_mod.F90 \
	$(PSYAD_WDIR)/$1/tl_%_kernel_mod.F90 | $(DIRECTORIES)

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$(word 2,$$^)))))

	# Running PSyAD
	echo "*Generating* $$(basename $$(notdir $$@)) for $$(basename $$(notdir $$(word 1,$$^)))"
	psyad -api lfric -t -otest $$@ -oad /dev/null -a $$(ACTIVE_VARS) \
		-c $$(PSYAD_CONFIG_FILE) -- $$(word 2,$$^)

$(PSYAD_WDIR)/$2/adjt_%_alg_mod.X90: $(PSYAD_WDIR)/$1/adj_%_kernel_mod.F90 \
	$(PSYAD_WDIR)/$1/%_kernel_mod.F90 | $(DIRECTORIES)

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$(word 2,$$^)))))

	# Running PSyAD
	echo "*Generating* $$(basename $$(notdir $$@)) for $$(basename $$(notdir $$(word 1,$$^)))"
	psyad -api lfric -t -otest $$@ -oad /dev/null -a $$(ACTIVE_VARS) \
		-c $$(PSYAD_CONFIG_FILE) -- $$(word 2,$$^)

# Runs PSyAD:
# For tl kernels.
$(PSYAD_WDIR)/$1/atl_%_kernel_mod.F90: $(PSYAD_WDIR)/$1/tl_%_kernel_mod.F90 \
	| $(DIRECTORIES)

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api lfric -oad $$@ -a $$(ACTIVE_VARS) \
		-c $$(PSYAD_CONFIG_FILE) -- $$<

# For regular kernels.
$(PSYAD_WDIR)/$1/adj_%_kernel_mod.F90: $(PSYAD_WDIR)/$1/%_kernel_mod.F90 \
	| $(DIRECTORIES)

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api lfric -oad $$@ -a $$(ACTIVE_VARS) \
		-c $$(PSYAD_CONFIG_FILE) -- $$<

#-----------------------------------------------------------------------------
# For f90s
#-----------------------------------------------------------------------------
$(PSYAD_WDIR)/$2/atlt_%_alg_mod.x90: $(PSYAD_WDIR)/$1/atl_%_kernel_mod.f90 \
	$(PSYAD_WDIR)/$1/tl_%_kernel_mod.f90 | $(DIRECTORIES)

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$(word 2,$$^)))))

	# Running PSyAD
	echo "*Generating* $$(basename $$(notdir $$@)) for $$(basename $$(notdir $$(word 1,$$^)))"
	psyad -api lfric -t -otest $$@ -oad /dev/null -a $$(ACTIVE_VARS) \
		-c $$(PSYAD_CONFIG_FILE) -- $$(word 2,$$^)

$(PSYAD_WDIR)/$2/adjt_%_alg_mod.x90: $(PSYAD_WDIR)/$1/adj_%_kernel_mod.f90 \
	$(PSYAD_WDIR)/$1/tl_%_kernel_mod.f90 | $(DIRECTORIES)

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$(word 2,$$^)))))

	# Running PSyAD
	echo "*Generating* $$(basename $$(notdir $$@)) for $$(basename $$(notdir $$(word 1,$$^)))"
	psyad -api lfric -t -otest $$@ -oad /dev/null -a $$(ACTIVE_VARS) \
		-c $$(PSYAD_CONFIG_FILE) -- $$(word 2,$$^)

# Runs PSyAD:
# For tl kernels.
$(PSYAD_WDIR)/$1/atl_%_kernel_mod.f90: $(PSYAD_WDIR)/$1/tl_%_kernel_mod.f90 \
	| $(DIRECTORIES)

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api lfric -oad $$@ -a $$(ACTIVE_VARS) \
		-c $$(PSYAD_CONFIG_FILE) -- $$<

# For regular kernels.
$(PSYAD_WDIR)/$1/adj_%_kernel_mod.f90: $(PSYAD_WDIR)/$1/%_kernel_mod.f90 \
	| $(DIRECTORIES)

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api lfric -oad $$@ -a $$(ACTIVE_VARS) \
		-c $$(PSYAD_CONFIG_FILE) -- $$<

endef # PSYAD

# Evaluating the PSYAD definition
# for each kernel path, with the algorithm path
# being generated from the kernel path via substitution.
# After this, all possible rules are now generated.
$(foreach kernel_path,$(KERNEL_PATHS),$(eval $(call PSYAD,$(kernel_path),$(subst kernel,algorithm,$(kernel_path)))))
