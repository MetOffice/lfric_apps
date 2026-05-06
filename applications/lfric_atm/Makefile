##############################################################################
# (c) Crown copyright 2017-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

PROJECT_NAME = lfric_atm

export PROFILE ?= fast-debug

# Select option for PSyclone transformation. Select "none" for no transformation.
export PSYCLONE_TRANSFORMATION ?= minimum

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))


export EXTERNAL_DYNAMIC_LIBRARIES = shum

# List of functions from shumlib that may be used by developers during debugging
export IGNORE_DEPENDENCIES = c_shum_byteswap \
                             f_shum_field_mod \
                             f_shum_file_mod \
                             f_shum_is_nan_mod \
                             f_shum_is_inf_mod \
                             f_shum_is_denormal_mod


# Note - socrates_interface must come before jules_interface as there are dependencies on each other.
export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/driver \
                               $(CORE_ROOT_DIR)/components/science \
                               $(CORE_ROOT_DIR)/components/inventory \
                               $(CORE_ROOT_DIR)/components/lfric-xios \
                               $(APPS_ROOT_DIR)/science/shared \
                               $(APPS_ROOT_DIR)/interfaces/socrates_interface \
                               $(APPS_ROOT_DIR)/interfaces/jules_interface \
                               $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface

META_VN             ?= HEAD
export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

# Check PSyclone transformation option is valid
ifneq ("$(PSYCLONE_TRANSFORMATION)","none")
  export OPTIMISATION_PATH ?= $(wildcard optimisation/$(PSYCLONE_TRANSFORMATION))
  ifeq ("$(wildcard $(OPTIMISATION_PATH)/psykal/global.py)","")
    $(error PSyclone transformation dir "$(PSYCLONE_TRANSFORMATION)" does not have a global.py script)
  endif
endif

.PHONY: default
default: build
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)
include build/project.mk

# Include transmute list only if not using minimum or no transformations
ifeq ($(filter "$(PSYCLONE_TRANSFORMATION)", "none" "minimum"),)
  include build/psyclone_transmute_file_list.mk
endif

##############################################################################
# Build
#
.PHONY: build
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export CXX_LINK     = YES
build: export PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \;)))
build: export PROJECT      = lfric_atm
build: export SCRATCH_DIR := $(WORKING_DIR)/..
build: export WORKING_DIR := $(WORKING_DIR)

build: export LDFLAGS_GROUPS = OPENMP
ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION
else ifeq "$(PROFILE)" "fast-debug"
# Use fast-debug environment vars from build/fortran/*.mk configs
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS FASTD_INIT FASTD_RUNTIME SAFE_OPTIMISATION
else ifeq "$(PROFILE)" "production"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS RISKY_OPTIMISATION
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build: ALWAYS
	$(call MESSAGE,========================================)
	$(call MESSAGE,Importing internal dependencies...)
	$(call MESSAGE,========================================)
	$(Q)for SUBPROJECT in $(INTERNAL_DEPENDENCIES) ; do \
	    $(MAKE) $(QUIET_ARG) -f $$SUBPROJECT/build/import.mk ; done
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting coupled_interface component)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/interfaces/coupled_interface/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/gungho/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting UM physics )
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/build/extract/extract_physics.mk
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Generating $(PROJECT) namelist loaders)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/gungho/source \
	          META_FILE_DIR=$(META_FILE_DIR)
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning interface components)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/interfaces/coupled_interface/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/interfaces/jules_interface/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/interfaces/socrates_interface/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/gungho/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Preprocess and PSyclone Transmute)
	$(call MESSAGE,========================================)
# Note depending on site target requirements, either all source or specific source will be affected
# See application overrides and files exported there
	$Q$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/build/pre_process_phys.mk \
              SOURCE_DIR=$(WORKING_DIR) \
              PSYCLONE_PHYSICS_FILES="$(PSYCLONE_PHYSICS_FILES)" \
	          PSYCLONE_DIRECTORIES="$(PSYCLONE_DIRECTORIES)" \
	          PSYCLONE_PHYSICS_EXCEPTION="$(PSYCLONE_PHYSICS_EXCEPTION)"
	$Q$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/build/psyclone_transmute_pass.mk \
	          SOURCE_DIR=$(WORKING_DIR) \
	          OPTIMISATION_PATH_PSY=$(APPS_ROOT_DIR)/applications/$(PROJECT)/$(OPTIMISATION_PATH) \
	          PSYCLONE_PASS_NO_SCRIPT="$(PSYCLONE_PASS_NO_SCRIPT)"
	$Q$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/build/psyclone_transmute.mk \
	          SOURCE_DIR=$(WORKING_DIR) \
	          OPTIMISATION_PATH_PSY=$(APPS_ROOT_DIR)/applications/$(PROJECT)/$(OPTIMISATION_PATH) \
	          PSYCLONE_PHYSICS_FILES="$(PSYCLONE_PHYSICS_FILES)" \
	          PSYCLONE_PASS_NO_SCRIPT="$(PSYCLONE_PASS_NO_SCRIPT)" \
	          PSYCLONE_DIRECTORIES="$(PSYCLONE_DIRECTORIES)" \
	          PSYCLONE_PHYSICS_EXCEPTION="$(PSYCLONE_PHYSICS_EXCEPTION)"
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Analysing $(PROJECT) build dependencies)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk


##############################################################################
# Pat_build step - This runs after the build step as relies on exe

# Pass $(PROJECT_NAME) which is exe name lfric_atm which pat_report needs
pat_build: build
	$(call MESSAGE,=========================================================)
	$(call MESSAGE, Begin Pat Build of $(PROJECT_NAME) executable)
	$(call MESSAGE,=========================================================)
	bash $(PAT_BUILD_PATH) $(PROJECT_NAME)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE, Pat Build of $(PROJECT_NAME) complete as $(PROJECT_NAME)+pat)
	$(call MESSAGE,=========================================================)

##############################################################################
# Unit tests
#
unit-tests:
	$(call MESSAGE,Testing,"There are no unit tests.")


##############################################################################
# Integration tests
#
integration-tests:
	$(call MESSAGE,Testing,"There are no integration tests.")


##############################################################################
# Clean
#
.PHONY: clean
clean: partial-clean
	$(call MESSAGE,Removing,"UM work space")
	$(Q)-rm -r $(WORKING_DIR)

.PHONY: partial-clean pclean
partial-clean pclean: ALWAYS
	$(call MESSAGE,Removing,"LFRic atmosphere  workspace")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"LFRic atmosphere binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi
