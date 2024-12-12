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
# Adjoint test actually generates with the adjoint, but having a two-target
# pattern rule causes race conditions when run in parallel. This faux rule
# exists so Make knows the adjoint test has been created after the adjoint
# is made. This implementation avoids such race conditions.
$(PSYAD_WDIR)/$2/atlt_%_alg_mod.X90: $(PSYAD_WDIR)/$1/atl_%_kernel_mod.F90 \
	| $(DIRECTORIES)
	echo $$(notdir $$@) generated with $$(notdir $$<)

$(PSYAD_WDIR)/$2/adjt_%_alg_mod.X90: $(PSYAD_WDIR)/$1/adj_%_kernel_mod.F90 \
	| $(DIRECTORIES)
	echo $$(notdir $$@) generated with $$(notdir $$<)

# Runs PSyAD:
# For tl kernels.
$(PSYAD_WDIR)/$1/atl_%_kernel_mod.F90: $(PSYAD_WDIR)/$1/tl_%_kernel_mod.F90 \
	| $(DIRECTORIES)

	# Generating ATLT_TARGET
	$$(eval ATLT_TARGET := $$(join $$(dir $$<), $$(subst tl_,atlt_,$$(notdir $$<))))
	$$(eval ATLT_TARGET := $$(subst kernel_mod,alg_mod,$$(ATLT_TARGET)))
	$$(eval ATLT_TARGET := $$(subst .F90,.X90,$$(ATLT_TARGET)))
	$$(eval ATLT_TARGET := $$(subst $$(PSYAD_WDIR)/kernel,$$(PSYAD_WDIR)/algorithm,$$(ATLT_TARGET)))

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api lfric -t -otest $$(ATLT_TARGET) -oad $$@ -a $$(ACTIVE_VARS) -- $$<

# For regular kernels.
$(PSYAD_WDIR)/$1/adj_%_kernel_mod.F90: $(PSYAD_WDIR)/$1/%_kernel_mod.F90 \
	| $(DIRECTORIES)

	# Generating ADJT_TARGET
	$$(eval ADJT_TARGET := $$(join $$(dir $$<), $$(addprefix adjt_, $$(notdir $$<))))
	$$(eval ADJT_TARGET := $$(subst kernel_mod,alg_mod,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst .F90,.X90,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst $$(PSYAD_WDIR)/kernel,$$(PSYAD_WDIR)/algorithm,$$(ADJT_TARGET)))

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api lfric -t -otest $$(ADJT_TARGET) -oad $$@ -a $$(ACTIVE_VARS) -- $$<

#-----------------------------------------------------------------------------
# For f90s
#-----------------------------------------------------------------------------
# Adjoint test actually generates with the adjoint, but having a two-target
# pattern rule causes race conditions when run in parallel. This faux rule
# exists so Make knows the adjoint test has been created after the adjoint
# is made. This implementation avoids such race conditions.
$(PSYAD_WDIR)/$2/atlt_%_alg_mod.x90: $(PSYAD_WDIR)/$1/atl_%_kernel_mod.f90 \
	| $(DIRECTORIES)
	echo $$(notdir $$@) generated with $$(notdir $$<)

$(PSYAD_WDIR)/$2/adjt_%_alg_mod.x90: $(PSYAD_WDIR)/$1/adj_%_kernel_mod.f90 \
	| $(DIRECTORIES)
	echo $$(notdir $$@) generated with $$(notdir $$<)

# Runs PSyAD:
# For tl kernels.
$(PSYAD_WDIR)/$1/atl_%_kernel_mod.f90: $(PSYAD_WDIR)/$1/tl_%_kernel_mod.f90 \
	| $(DIRECTORIES)

	# Generating ATLT_TARGET
	$$(eval ATLT_TARGET := $$(join $$(dir $$<), $$(subst tl_,atlt_,$$(notdir $$<))))
	$$(eval ATLT_TARGET := $$(subst kernel_mod,alg_mod,$$(ATLT_TARGET)))
	$$(eval ATLT_TARGET := $$(subst .f90,.x90,$$(ATLT_TARGET)))
	$$(eval ATLT_TARGET := $$(subst $$(PSYAD_WDIR)/kernel,$$(PSYAD_WDIR)/algorithm,$$(ATLT_TARGET)))

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api lfric -t -otest $$(ATLT_TARGET) -oad $$@ -a $$(ACTIVE_VARS) -- $$<

# For regular kernels.
$(PSYAD_WDIR)/$1/adj_%_kernel_mod.f90: $(PSYAD_WDIR)/$1/%_kernel_mod.f90 \
	| $(DIRECTORIES)

	# Generating ADJT_TARGET
	$$(eval ADJT_TARGET := $$(join $$(dir $$<), $$(addprefix adjt_, $$(notdir $$<))))
	$$(eval ADJT_TARGET := $$(subst kernel_mod,alg_mod,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst .f90,.x90,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst $$(PSYAD_WDIR)/kernel,$$(PSYAD_WDIR)/algorithm,$$(ADJT_TARGET)))

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api lfric -t -otest $$(ADJT_TARGET) -oad $$@ -a $$(ACTIVE_VARS) -- $$<

endef # PSYAD

# Evaluating the PSYAD definition
# for each kernel path, with the algorithm path
# being generated from the kernel path via substitution.
# After this, all possible rules are now generated.
$(foreach kernel_path,$(KERNEL_PATHS),$(eval $(call PSYAD,$(kernel_path),$(subst kernel,algorithm,$(kernel_path)))))
