##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# This script runs PSyAD on all non-adjoint files in $(PSYAD_WDIR).

# There are tl kernel files and regular kernel files to process.
# PSYAD_WDIR may have some adjoint kernels on incremental builds,
# so we only want kernel files without the adj prefix.

include $(ADJOINT_BUILD)/psyad_vars.mk

KERNEL_SOURCES := $(shell find $(PSYAD_WDIR)/kernel -name '*_kernel_mod.[Ff]90' ! -name 'adj_*_kernel_mod.[Ff]90' -print)

# tl-prefixed kernel files.
TLK_SOURCES    := $(shell find $(PSYAD_WDIR)/kernel -name 'tl_*_kernel_mod.[Ff]90' -print)

# Non-prefixed (regular) kernel files.
REGK_SOURCES   := $(filter-out $(TLK_SOURCES), $(KERNEL_SOURCES))

# tl kernel files
TLK_ADJ_TARGETS := $(subst tl,adj,$(TLK_SOURCES))

# And their corresponding adjoint test names
TLK_ADJT_TARGETS := $(subst tl,adjt,$(TLK_SOURCES))
TLK_ADJT_TARGETS := $(subst kernel_mod,alg_mod,$(TLK_ADJT_TARGETS))
TLK_ADJT_TARGETS := $(subst .F90,.X90,$(TLK_ADJT_TARGETS))
TLK_ADJT_TARGETS := $(subst .f90,.x90,$(TLK_ADJT_TARGETS))
TLK_ADJT_TARGETS := $(subst $(PSYAD_WDIR)/kernel,$(PSYAD_WDIR)/algorithm,$(TLK_ADJT_TARGETS))

# Regular kernel files
REGK_ADJ_TARGETS := $(join $(dir $(REGK_SOURCES)), $(addprefix adj_, $(notdir $(REGK_SOURCES))))

# And their corresponding adjoint tests
REGK_ADJT_TARGETS := $(join $(dir $(REGK_SOURCES)), $(addprefix adjt_, $(notdir $(REGK_SOURCES))))
REGK_ADJT_TARGETS := $(subst kernel_mod,alg_mod,$(REGK_ADJT_TARGETS))
REGK_ADJT_TARGETS := $(subst .F90,.X90,$(REGK_ADJT_TARGETS))
REGK_ADJT_TARGETS := $(subst .f90,.x90,$(REGK_ADJT_TARGETS))
REGK_ADJT_TARGETS := $(subst $(PSYAD_WDIR)/kernel,$(PSYAD_WDIR)/algorithm,$(REGK_ADJT_TARGETS))

# Total targets
ADJ_TARGETS  := $(TLK_ADJ_TARGETS) $(REGK_ADJ_TARGETS)
ADJT_TARGETS := $(TLK_ADJT_TARGETS) $(REGK_ADJT_TARGETS)

DIRECTORIES := $(addprefix $(PSYAD_WDIR)/,$(KERNEL_PATHS))
DIRECTORIES += $(addprefix $(PSYAD_WDIR)/,$(ALG_PATHS))

all: $(ADJ_TARGETS) $(ADJT_TARGETS)

# Use define block here to account for subdirectories.
# $1 is kernel path, $2 is algorithm path.
# Runs PSyAD
define GEN_PSYAD

##############################################################################
# F90
##############################################################################
#-----------------------------------------------------------------------------
# For tl kernels
#-----------------------------------------------------------------------------
$(PSYAD_WDIR)/$1/adj_%_kernel_mod.F90 $(PSYAD_WDIR)/$2/adjt_%_alg_mod.X90: \
	$(PSYAD_WDIR)/$1/tl_%_kernel_mod.F90 | $(DIRECTORIES)

	# Generating ADJT_TARGET since Make does not
	# allow one to get it from the target list.
	$$(eval ADJT_TARGET := $$(subst tl,adjt,$$<))
	$$(eval ADJT_TARGET := $$(subst kernel_mod,alg_mod,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst .F90,.X90,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst $$(PSYAD_WDIR)/kernel,$$(PSYAD_WDIR)/algorithm,$$(ADJT_TARGET)))

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api dynamo0.3 -t -otest $$(ADJT_TARGET) -oad $$@ -a $$(ACTIVE_VARS) -- $$<

#-----------------------------------------------------------------------------
# For regular kernels
#-----------------------------------------------------------------------------
$(PSYAD_WDIR)/$1/adj_%_kernel_mod.F90 $(PSYAD_WDIR)/$2/adjt_%_alg_mod.X90: \
	$(PSYAD_WDIR)/$1/%_kernel_mod.F90 | $(DIRECTORIES)

	# Generating ADJT_TARGET since Make does not
	# allow one to get it from the target list.
	$$(eval ADJT_TARGET := $$(join $$(dir $$<), $$(addprefix adjt_, $$(notdir $$<))))
	$$(eval ADJT_TARGET := $$(subst kernel_mod,alg_mod,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst .F90,.X90,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst $$(PSYAD_WDIR)/kernel,$$(PSYAD_WDIR)/algorithm,$$(ADJT_TARGET)))

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api dynamo0.3 -t -otest $$(ADJT_TARGET) -oad $$@ -a $$(ACTIVE_VARS) -- $$<

##############################################################################
# f90
##############################################################################
#-----------------------------------------------------------------------------
# For tl kernels
#-----------------------------------------------------------------------------
$(PSYAD_WDIR)/$1/adj_%_kernel_mod.f90 $(PSYAD_WDIR)/$2/adjt_%_alg_mod.x90: \
	$(PSYAD_WDIR)/$1/tl_%_kernel_mod.f90 | $(DIRECTORIES)

	# Generating ADJT_TARGET since Make does not
	# allow one to get it from the target list.
	$$(eval ADJT_TARGET := $$(subst tl,adjt,$$<))
	$$(eval ADJT_TARGET := $$(subst kernel_mod,alg_mod,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst .f90,.x90,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst $$(PSYAD_WDIR)/kernel,$$(PSYAD_WDIR)/algorithm,$$(ADJT_TARGET)))

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api dynamo0.3 -t -otest $$(ADJT_TARGET) -oad $$@ -a $$(ACTIVE_VARS) -- $$<

#-----------------------------------------------------------------------------
# For regular kernels
#-----------------------------------------------------------------------------
$(PSYAD_WDIR)/$1/adj_%_kernel_mod.f90 $(PSYAD_WDIR)/$2/adjt_%_alg_mod.x90: \
	$(PSYAD_WDIR)/$1/%_kernel_mod.f90 | $(DIRECTORIES)

	# Generating ADJT_TARGET since Make does not
	# allow one to get it from the target list.
	$$(eval ADJT_TARGET := $$(join $$(dir $$<), $$(addprefix adjt_, $$(notdir $$<))))
	$$(eval ADJT_TARGET := $$(subst kernel_mod,alg_mod,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst .f90,.x90,$$(ADJT_TARGET)))
	$$(eval ADJT_TARGET := $$(subst $$(PSYAD_WDIR)/kernel,$$(PSYAD_WDIR)/algorithm,$$(ADJT_TARGET)))

	# Evaluating active variables located in $(ADJOINT_BUILD)/psyad_vars.mk
	$$(eval ACTIVE_VARS := $$(ACTIVE_$$(basename $$(notdir $$<))))

	# Running PSyAD
	echo "*Running* PSyAD on $$(basename $$(notdir $$<))"
	psyad -api dynamo0.3 -t -otest $$(ADJT_TARGET) -oad $$@ -a $$(ACTIVE_VARS) -- $$<

endef

# Evaluating the GEN_PSYAD definition
# for each kernel path, with the algorithm path
# being generated from the kernel path via substitution.
# After this, all possible rules are now generated.
$(foreach kernel_path,$(KERNEL_PATHS),$(eval $(call GEN_PSYAD,$(kernel_path),$(subst kernel,algorithm,$(kernel_path)))))

$(DIRECTORIES):
	mkdir -p $@

