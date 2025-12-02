##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# This script contains determination of all the targets for the three stages of PSyAD building.

# List of paths to files we wish to compile using PSyAD.
PSYAD_FILES_CORE := $(shell cat $(ADJOINT_BUILD)/psyad_files_list_core.txt)
PSYAD_FILES_CORE := $(addprefix $(CORE_ROOT_DIR)/,$(PSYAD_FILES_CORE))
PSYAD_FILES_APPS := $(shell cat $(ADJOINT_BUILD)/psyad_files_list_apps.txt)
PSYAD_FILES_APPS := $(addprefix $(APPS_ROOT_DIR)/,$(PSYAD_FILES_APPS))
PSYAD_FILES := $(PSYAD_FILES_CORE) $(PSYAD_FILES_APPS)

# Get the list of base directories where psyad kernels are located.
BASE_DIRS := $(shell python $(ADJOINT_BUILD)/print_paths.py -f $(PSYAD_FILES) -opt "base")

# List of all kernel related subdirectories in the base directories where PSyAD files are located.
KERNEL_PATHS := $(shell python $(ADJOINT_BUILD)/print_paths.py -f $(PSYAD_FILES) -opt "kernel")
ALG_PATHS    := $(subst kernel,algorithm,$(KERNEL_PATHS))

##############################################################################
# PRE-PATCH TARGETS
##############################################################################
# Determine the pre-patch target paths for the psyad kernels.
PRE_PATCH_TARGETS := $(PSYAD_FILES)
$(foreach base_dir, \
          $(BASE_DIRS), \
          $(eval PRE_PATCH_TARGETS := $(subst $(base_dir),$(PSYAD_WDIR),$(PRE_PATCH_TARGETS)) \
          ) \
)

##############################################################################
# PSyAD TARGETS
##############################################################################
# We need to filter based on whether or not the kernel working files are prefixed by tl_ or not.
# Generating the tl prefixed kernels.
KERNEL_DIR_LIST := $(addprefix $(PSYAD_WDIR)/,$(KERNEL_PATHS))
$(foreach dir_, \
          $(KERNEL_DIR_LIST), \
          $(eval WORKING_TL += $(filter $(dir_)/tl_%,$(PRE_PATCH_TARGETS)) \
          ) \
)
# Regular files have no tl_ prefix.
WORKING_REG := $(filter-out $(WORKING_TL), $(PRE_PATCH_TARGETS))

# Now we can generate the PSyAD targets based on the different prefixes.
# tl kernel files
ATL_TARGETS := $(join $(dir $(WORKING_TL)), $(subst tl_,atl_,$(notdir $(WORKING_TL))))

# And their corresponding adjoint tests.
ATLT_TARGETS := $(join $(dir $(WORKING_TL)), $(subst tl_,atlt_,$(notdir $(WORKING_TL))))
ATLT_TARGETS := $(subst kernel_mod,alg_mod,$(ATLT_TARGETS))
ATLT_TARGETS := $(subst .F90,.X90,$(ATLT_TARGETS))
ATLT_TARGETS := $(subst .f90,.x90,$(ATLT_TARGETS))
ATLT_TARGETS := $(subst $(PSYAD_WDIR)/kernel,$(PSYAD_WDIR)/algorithm,$(ATLT_TARGETS))

# Regular kernel files
ADJ_TARGETS := $(join $(dir $(WORKING_REG)), $(addprefix adj_, $(notdir $(WORKING_REG))))

# And their corresponding adjoint tests.
ADJT_TARGETS := $(join $(dir $(WORKING_REG)), $(addprefix adjt_, $(notdir $(WORKING_REG))))
ADJT_TARGETS := $(subst kernel_mod,alg_mod,$(ADJT_TARGETS))
ADJT_TARGETS := $(subst .F90,.X90,$(ADJT_TARGETS))
ADJT_TARGETS := $(subst .f90,.x90,$(ADJT_TARGETS))
ADJT_TARGETS := $(subst $(PSYAD_WDIR)/kernel,$(PSYAD_WDIR)/algorithm,$(ADJT_TARGETS))

# Total targets for PSyAD stage.
PSYADJ_TARGETS  := $(ATL_TARGETS) $(ADJ_TARGETS)
PSYADJT_TARGETS := $(ATLT_TARGETS) $(ADJT_TARGETS)
ifeq "$(BUILD_ADJ_TESTS)" "TRUE"
	PSYAD_TARGETS := $(PSYADJ_TARGETS) $(PSYADJT_TARGETS)
else
	PSYAD_TARGETS := $(PSYADJ_TARGETS)
endif

##############################################################################
# POST-PATCH TARGETS
##############################################################################
POST_ADJ_TARGETS := $(subst $(PSYAD_WDIR)/kernel,$(WORKING_DIR)/kernel,$(PSYADJ_TARGETS))
POST_ADJT_TARGETS := $(subst $(PSYAD_WDIR)/algorithm,$(WORKING_DIR)/algorithm,$(PSYADJT_TARGETS))
ifeq "$(BUILD_ADJ_TESTS)" "TRUE"
	POST_PATCH_TARGETS := $(POST_ADJ_TARGETS) $(POST_ADJT_TARGETS)
else
	POST_PATCH_TARGETS := $(POST_ADJ_TARGETS)
endif

# Necessary directories that need to be created.
# Pre-patch.
DIRECTORIES := $(addprefix $(PSYAD_WDIR)/,$(KERNEL_PATHS))
# PSyAD.
DIRECTORIES += $(addprefix $(PSYAD_WDIR)/,$(ALG_PATHS))
# Post-patch.
DIRECTORIES += $(addprefix $(WORKING_DIR)/,$(KERNEL_PATHS))
DIRECTORIES += $(addprefix $(WORKING_DIR)/,$(ALG_PATHS))
