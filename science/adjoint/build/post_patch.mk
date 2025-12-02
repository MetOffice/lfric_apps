##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# POST-PATCH STAGE RULES
##############################################################################

define POST_PATCH

# Use define block here to account for the subdirectories in base/source/kernel.
# $1 is kernel path, $2 is algorithm path.
#-----------------------------------------------------------------------------
# For F90s and X90s
#-----------------------------------------------------------------------------
# For atl kernels and atlt algorithms.
# If a patch exists, copy and patch the generated adjoint into working.
$(WORKING_DIR)/$1/atl_%_mod.F90: \
$(PSYAD_WDIR)/$1/atl_%_mod.F90 \
$(PATCH_DIR)/kernel/atl_%_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If a patch exists, copy and patch the generated adjoint test into working.
$(WORKING_DIR)/$2/atlt_%_alg_mod.X90: \
$(PSYAD_WDIR)/$2/atlt_%_alg_mod.X90 \
$(PATCH_DIR)/algorithm/atlt_%_alg_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy the generated adjoint into working.
$(WORKING_DIR)/$1/atl_%_mod.F90: \
$(PSYAD_WDIR)/$1/atl_%_mod.F90 | $(DIRECTORIES)
	cp $$< $$@

# If no patch exists, just copy the generated adjoint test into working.
$(WORKING_DIR)/$2/atlt_%_alg_mod.X90: \
$(PSYAD_WDIR)/$2/atlt_%_alg_mod.X90 | $(DIRECTORIES)
	cp $$< $$@

# For regular kernels and algorithms.
# If a patch exists, copy and patch the generated adjoint into working.
$(WORKING_DIR)/$1/adj_%_mod.F90: \
$(PSYAD_WDIR)/$1/adj_%_mod.F90 \
$(PATCH_DIR)/kernel/adj_%_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If a patch exists, copy and patch the generated adjoint test into working.
$(WORKING_DIR)/$2/adjt_%_alg_mod.X90: \
$(PSYAD_WDIR)/$2/adjt_%_alg_mod.X90 \
$(PATCH_DIR)/algorithm/adjt_%_alg_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy the generated adjoint into working.
$(WORKING_DIR)/$1/adj_%_mod.F90: \
$(PSYAD_WDIR)/$1/adj_%_mod.F90 | $(DIRECTORIES)
	cp $$< $$@

# If no patch exists, just copy the generated adjoint test into working.
$(WORKING_DIR)/$2/adjt_%_alg_mod.X90: \
$(PSYAD_WDIR)/$2/adjt_%_alg_mod.X90 | $(DIRECTORIES)
	cp $$< $$@

#-----------------------------------------------------------------------------
# For f90s and x90s
#-----------------------------------------------------------------------------
# For atl kernels and atlt algorithms.
# If a patch exists, copy and patch the generated adjoint into working.
$(WORKING_DIR)/$1/atl_%_mod.f90: \
$(PSYAD_WDIR)/$1/atl_%_mod.f90 \
$(PATCH_DIR)/kernel/atl_%_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If a patch exists, copy and patch the generated adjoint test into working.
$(WORKING_DIR)/$2/atlt_%_alg_mod.x90: \
$(PSYAD_WDIR)/$2/atlt_%_alg_mod.x90 \
$(PATCH_DIR)/algorithm/atlt_%_alg_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy the generated adjoint into working.
$(WORKING_DIR)/$1/atl_%_mod.f90: \
$(PSYAD_WDIR)/$1/atl_%_mod.f90 | $(DIRECTORIES)
	cp $$< $$@

# If no patch exists, just copy the generated adjoint test into working.
$(WORKING_DIR)/$2/atlt_%_alg_mod.x90: \
$(PSYAD_WDIR)/$2/atlt_%_alg_mod.x90 | $(DIRECTORIES)
	cp $$< $$@

# For regular kernels and algorithms.
# If a patch exists, copy and patch the generated adjoint into working.
$(WORKING_DIR)/$1/adj_%_mod.f90: \
$(PSYAD_WDIR)/$1/adj_%_mod.f90 \
$(PATCH_DIR)/kernel/adj_%_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If a patch exists, copy and patch the generated adjoint test into working.
$(WORKING_DIR)/$2/adjt_%_alg_mod.x90: \
$(PSYAD_WDIR)/$2/adjt_%_alg_mod.x90 \
$(PATCH_DIR)/algorithm/adjt_%_alg_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy the generated adjoint into working.
$(WORKING_DIR)/$1/adj_%_mod.f90: \
$(PSYAD_WDIR)/$1/adj_%_mod.f90 | $(DIRECTORIES)
	cp $$< $$@

# If no patch exists, just copy the generated adjoint test into working.
$(WORKING_DIR)/$2/adjt_%_alg_mod.x90: \
$(PSYAD_WDIR)/$2/adjt_%_alg_mod.x90 | $(DIRECTORIES)
	cp $$< $$@

endef # POST_PATCH

# Evaluating the POST_PATCH definition
# for each kernel path, with the algorithm path
# being generated from the kernel path via substitution.
# After this, all possible rules are now generated.
$(foreach kernel_path,$(KERNEL_PATHS),$(eval $(call POST_PATCH,$(kernel_path),$(subst kernel,algorithm,$(kernel_path)))))

