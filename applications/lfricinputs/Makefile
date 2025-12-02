##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:
export PROJECT_NAME = lfricinputs
export PROJECT_DIR_NAME = lfricinputs
export PROFILE ?= fast-debug

# Select option for PSyclone transformation. Select "none" for no transformation.
export PSYCLONE_TRANSFORMATION ?= minimum

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/science \
                               $(CORE_ROOT_DIR)/components/inventory \
                               $(CORE_ROOT_DIR)/components/driver \
                               $(CORE_ROOT_DIR)/components/lfric-xios \
                               $(APPS_ROOT_DIR)/science/gungho

export EXTERNAL_DYNAMIC_LIBRARIES = shum

export IGNORE_DEPENDENCIES = c_shum_byteswap \
                             f_shum_ff_status_mod \
                             f_shum_field_mod \
                             f_shum_fieldsfile_mod \
                             f_shum_file_mod \
                             f_shum_fixed_length_header_indices_mod \
                             f_shum_lookup_indices_mod \
                             f_shum_stashmaster_mod

META_VN ?= HEAD
export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)/$(META_VN)

# Check PSyclone transformation option is valid
ifneq ("$(PSYCLONE_TRANSFORMATION)","none")
  export OPTIMISATION_PATH ?= $(wildcard optimisation/$(PSYCLONE_TRANSFORMATION))
  ifeq ("$(wildcard $(OPTIMISATION_PATH)/global.py)","")
    $(error PSyclone transformation dir "$(PSYCLONE_TRANSFORMATION)" does not have a global.py script)
  endif
endif

.PHONY: default
default: build

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)
-include $(PROJECT_DIR)/build/project.mk

##############################################################################
# Build
#
.PHONY: build
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export CXX_LINK     = TRUE
build: export PROGRAMS    ?= $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -print)))
build: export PROJECT      = $(PROJECT_NAME)
build: export SOURCE_DIR   = source
build: export WORKING_DIR := $(WORKING_DIR)/inputs

build: export LDFLAGS_GROUPS = OPENMP
ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS SAFE_OPTIMISATION FORTRAN_STANDARD
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
	$(call MESSAGE,Extracting $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Generating $(PROJECT) namelist loaders)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
	          META_FILE_DIR=$(META_FILE_DIR)
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Analysing $(PROJECT) build dependencies)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
	          BIN_DIR=$(BIN_DIR) PROGRAMS="$(PROGRAMS)"

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
				$(call MESSAGE,Removing,"$(PROJECT_NAME) work space")
				$(Q)-if [ $(WORKING_DIR) != *[\*]* ] && [ -d $(WORKING_DIR) ] ; then rm -r $(WORKING_DIR) ; fi
				$(call MESSAGE,Removing,"$(PROJECT_NAME) binaries")
				$(Q)-if [ -d bin ] ; then rm -r bin ; fi
