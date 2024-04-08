##############################################################################
# (c) Crown copyright 2017-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

PROJECT_NAME = lfric_atm

export PROFILE ?= fast-debug

# Option defined in fcm_make/target-map.txt which in turn references options
# from fcm:um.xm_tr/fcm-make
export UM_FCM_TARGET_PLATFORM ?= generic-x86
# Select option for PSyclone transformation. Select "none" for no transformation.
export PSYCLONE_TRANSFORMATION ?= minimum

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/driver \
                               $(CORE_ROOT_DIR)/components/science \
                               $(CORE_ROOT_DIR)/components/inventory \
                               $(CORE_ROOT_DIR)/components/lfric-xios

export UMPHYSICS_TARGETS = extract-um-physics
export UMPHYSICS_CLEAN   = clean-um-physics
export UM_FCM=$(UM_PHYS_DIR)/

META_VN             ?= HEAD
export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

# Check PSyclone transformation option is valid
ifneq ("$(PSYCLONE_TRANSFORMATION)","none")
  export OPTIMISATION_PATH ?= $(wildcard optimisation/$(PSYCLONE_TRANSFORMATION))
  ifeq ("$(wildcard $(OPTIMISATION_PATH)/global.py)","")
    $(error PSyclone transformation dir "$(PSYCLONE_TRANSFORMATION)" does not have a global.py script)
  endif
endif

.PHONY: default
default: build
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)
include build/project.mk

##############################################################################
# Build
#
.PHONY: build
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export CXX_LINK     = YES
build: export PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \;)))
build: export PROJECT      = lfric_atm
build: export SCRATCH_DIR := $(WORKING_DIR)/../um-scratch
build: export WORKING_DIR := $(WORKING_DIR)
ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION
else ifeq "$(PROFILE)" "fast-debug"
# Use fast-debug environment vars from build/fortran/*.mk configs
build: export FFLAG_GROUPS = DEBUG WARNINGS FASTD_INIT FASTD_RUNTIME SAFE_OPTIMISATION
else ifeq "$(PROFILE)" "production"
build: export FFLAG_GROUPS = DEBUG WARNINGS RISKY_OPTIMISATION
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
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/coupled_interface/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/gungho/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting UM physics )
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/build/extract/extract_physics.mk
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/um_physics_interface/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/socrates_interface/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/jules_interface/source
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
	$(call MESSAGE,PSycloning coupled_interface component)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/coupled_interface/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/gungho/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning UM Physics)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/um_physics_interface/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/socrates_interface/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(APPS_ROOT_DIR)/science/jules_interface/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=source
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Analysing $(PROJECT) build dependencies)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk

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
