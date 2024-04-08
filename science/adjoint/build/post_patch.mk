##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# This script handles the patching of the raw outputs from PSyAD.
# The sources are the adjoint kernels and adjoint test algorithms in $(PSYAD_WDIR).
# The targets are similar to the sources except the $(PSYAD_WDIR) part of the path
# is replaced with the project $(WORKING_DIR).

# The post-patching consists of two phases. First,
# the source files are copied to their equivalent
# target. Then, if a patch exists, the target file
# is patched.

ADJ_SOURCES := $(shell find $(PSYAD_WDIR)/kernel -name 'adj_*_kernel_mod.[Ff]90' -print)
ADJT_SOURCES := $(shell find $(PSYAD_WDIR)/algorithm -name 'adjt_*_alg_mod.[Xx]90' -print)
ADJ_TARGETS := $(subst $(PSYAD_WDIR)/kernel,$(WORKING_DIR)/kernel,$(ADJ_SOURCES))
ADJT_TARGETS := $(subst $(PSYAD_WDIR)/algorithm,$(WORKING_DIR)/algorithm,$(ADJT_SOURCES))

DIRECTORIES += $(addprefix $(WORKING_DIR)/,$(KERNEL_PATHS))
DIRECTORIES += $(addprefix $(WORKING_DIR)/,$(ALG_PATHS))

all: $(ADJ_TARGETS) $(ADJT_TARGETS)

# Use define block here to account for the subdirectories in base/source/kernel.
# $1 is kernel path, $2 is algorithm path.
define POST_PATCH
#-----------------------------------------------------------------------------
# For F90s and X90s
#-----------------------------------------------------------------------------
# If a patch exists, copy and patch the generated adjoint into working.
$(WORKING_DIR)/$1/adj_%_mod.F90: \
$(PSYAD_WDIR)/$1/adj_%_mod.F90 \
$(PATCH_DIR)/kernel/adj_%_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If a patch exists, copy and patch the generated adjoint test into working.
$(WORKING_DIR)/$2/adjt_%_alg_mod.X90: $(PSYAD_WDIR)/$2/adjt_%_alg_mod.X90 \
$(PATCH_DIR)/algorithm/adjt_%_alg_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy the generated adjoint into working.
$(WORKING_DIR)/$1/adj_%_mod.F90: $(PSYAD_WDIR)/$1/adj_%_mod.F90 | $(DIRECTORIES)
	cp $$< $$@

# If no patch exists, just copy the generated adjoint test into working.
$(WORKING_DIR)/$2/adjt_%_alg_mod.X90: $(PSYAD_WDIR)/$2/adjt_%_alg_mod.X90 | $(DIRECTORIES)
	cp $$< $$@

#-----------------------------------------------------------------------------
# For f90s and x90s
#-----------------------------------------------------------------------------
# If a patch exists, copy and patch the generated adjoint into working.
$(WORKING_DIR)/$1/adj_%_mod.f90: \
$(PSYAD_WDIR)/$1/adj_%_mod.f90 \
$(PATCH_DIR)/kernel/adj_%_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If a patch exists, copy and patch the generated adjoint test into working.
$(WORKING_DIR)/$2/adjt_%_alg_mod.x90: $(PSYAD_WDIR)/$2/adjt_%_alg_mod.x90 \
$(PATCH_DIR)/algorithm/adjt_%_alg_mod.patch | $(DIRECTORIES)
	cp $$< $$@
	patch $$@ $$(word 2,$$^)

# If no patch exists, just copy the generated adjoint into working.
$(WORKING_DIR)/$1/adj_%_mod.f90: $(PSYAD_WDIR)/$1/adj_%_mod.f90 | $(DIRECTORIES)
	cp $$< $$@

# If no patch exists, just copy the generated adjoint test into working.
$(WORKING_DIR)/$2/adjt_%_alg_mod.x90: $(PSYAD_WDIR)/$2/adjt_%_alg_mod.x90 | $(DIRECTORIES)
	cp $$< $$@

endef

# Evaluating the POST_PATCH definition
# for each kernel path, with the algorithm path
# being generated from the kernel path via substitution.
# After this, all possible rules are now generated.
$(foreach kernel_path,$(KERNEL_PATHS),$(eval $(call POST_PATCH,$(kernel_path),$(subst kernel,algorithm,$(kernel_path)))))

$(DIRECTORIES):
	mkdir -p $@

