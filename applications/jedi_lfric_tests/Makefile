##############################################################################
# (c) Crown copyright 2022-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_NAME = jedi_lfric_tests
export PROJECT_DIR_NAME = jedi_lfric_tests
export PROFILE ?= fast-debug

# Select option for PSyclone transformation. Select "none" for no transformation.
export PSYCLONE_TRANSFORMATION ?= minimum

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/driver \
                               $(CORE_ROOT_DIR)/components/science \
                               $(CORE_ROOT_DIR)/components/inventory \
                               $(APPS_ROOT_DIR)/interfaces/jedi_lfric_interface \
                               $(CORE_ROOT_DIR)/components/lfric-xios \
                               $(APPS_ROOT_DIR)/science/gungho \
                               $(APPS_ROOT_DIR)/science/linear \
                               $(APPS_ROOT_DIR)/science/adjoint

META_VN             ?= HEAD
export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)/$(META_VN)

# Check PSyclone transformation option is valid
ifneq ("$(PSYCLONE_TRANSFORMATION)","none")
  export OPTIMISATION_PATH ?= $(wildcard optimisation/$(PSYCLONE_TRANSFORMATION))
  ifeq ("$(wildcard $(OPTIMISATION_PATH)/psykal/global.py)","")
    $(error PSyclone transformation dir "$(PSYCLONE_TRANSFORMATION)" does not have a global.py script)
  endif
endif

.PHONY: default
default: build integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-latex document-api
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)
-include $(PROJECT_DIR)/build/project.mk

##############################################################################
# Documentation
#
.PHONY: document-uml
document-uml: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/uml
document-uml: SOURCE_DIR    = documentation/uml
document-uml: WORKING_DIR  := $(WORKING_DIR)/uml
document-uml: uml-documentation
	$(Q)echo > /dev/null

.PHONY: document-latex
document-latex: DOCUMENT_DIR  ?= $(PROJECT_DIR)/documents
document-latex: SOURCE_DIR     = documentation
document-latex: DOCUMENTS      = $(shell find $(SOURCE_DIR) -name '*.latex' -print)
document-latex: WORKING_DIR   := $(WORKING_DIR)/latex
document-latex: TEX_STUFF      = documentation/tex
document-latex: COMMON_FIGURES = documentation/common-figures
document-latex: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/latex.mk       \
                             SOURCE_DIR=$(SOURCE_DIR)         \
	                         WORKING_DIR=$(WORKING_DIR)       \
	                         DOCUMENT_DIR=$(DOCUMENT_DIR)     \
	                         TEX_STUFF=$(TEX_STUFF)           \
	                         COMMON_FIGURES=$(COMMON_FIGURES) \
	                         DOCUMENTS="$(DOCUMENTS)"

.PHONY: document-api
document-api: PROJECT       = $(PROJECT_NAME)
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null


##############################################################################
# Build
#
.PHONY: build
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export LIB_DIR     ?= $(PROJECT_DIR)/lib
build: export MOD_DIR     ?= $(PROJECT_DIR)/mod
build: export CXX_LINK     = TRUE
build: export PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -print)))
build: export PROJECT      = $(PROJECT_NAME)
build: export WORKING_DIR := $(WORKING_DIR)/jedi_lfric_tests

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
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Generating $(PROJECT) namelist loaders)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
	          PROJECT=$(PROJECT) \
	          SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR) \
	          META_FILE_DIR=$(META_FILE_DIR)
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=source \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Analysing $(PROJECT) build dependencies)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Archiving $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) libraries $(QUIET_ARG) -C $(WORKING_DIR) \
	          -f $(LFRIC_BUILD)/compile.mk \
	          LIB_DIR=$(LIB_DIR) MOD_DIR=$(MOD_DIR) \
	          PROGRAMS="$(PROGRAMS)"
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) applications $(QUIET_ARG) -C $(WORKING_DIR) \
	          -f $(LFRIC_BUILD)/compile.mk \
	          BIN_DIR=$(BIN_DIR) PROGRAMS="$(PROGRAMS)"


##############################################################################
# Unit tests
#
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export FFLAGS += -DUNIT_TEST
unit-tests/%: export IMPORT_PARTS  = $(CORE_ROOT_DIR)/infrastructure \
                                     $(CORE_ROOT_DIR)/components/driver \
                                     $(CORE_ROOT_DIR)/components/science \
                                     $(CORE_ROOT_DIR)/components/lfric-xios \
                                     $(CORE_ROOT_DIR)/components/inventory \
                                     $(APPS_ROOT_DIR)/interfaces/jedi_lfric_interface \
                                     $(APPS_ROOT_DIR)/science/gungho \
                                     $(APPS_ROOT_DIR)/science/linear
unit-tests/%: export META_FILE_DIR = rose-meta/jedi_lfric_tests/HEAD
unit-tests/%: export PROGRAMS = $(PROJECT_NAME)_unit_tests
unit-tests/%: export PROJECT = $(PROJECT_NAME)
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-test
unit-tests: unit-tests/run


##############################################################################
# Integration tests
#
integration-tests/%: export BIN_DIR      ?= $(PROJECT_DIR)/test
integration-tests/%: export FFLAGS       += -DINT_TEST
integration-tests/%: export IMPORT_PARTS  = $(CORE_ROOT_DIR)/infrastructure \
                                            $(CORE_ROOT_DIR)/components/science \
                                            $(CORE_ROOT_DIR)/components/driver \
                                            $(CORE_ROOT_DIR)/components/inventory \
                                            $(CORE_ROOT_DIR)/components/lfric-xios \
                                            $(APPS_ROOT_DIR)/interfaces/jedi_lfric_interface \
                                            $(APPS_ROOT_DIR)/science/gungho \
                                            $(APPS_ROOT_DIR)/science/linear \
                                            $(APPS_ROOT_DIR)/science/adjoint
integration-tests/%: export META_FILE_DIR = rose-meta/jedi_lfric_tests/HEAD
integration-tests/%: export PROJECT       = $(PROJECT_NAME)
integration-tests/%: export SOURCE_DIR    = source
integration-tests/%: export TEST_DIR      = integration-test
integration-tests/%: export WORKING_DIR  := $(WORKING_DIR)/integration-tests
integration-tests: integration-tests/run


##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"$(PROJECT_NAME) work space")
	$(Q)-if [ $(WORKING_DIR) != *[\*]* ] && [ -d $(WORKING_DIR) ] ; then rm -r $(WORKING_DIR) ; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) documents")
	$(Q)-if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi
	$(Q)-if [ -d test ] ; then rm -r test ; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) libraries")
	$(Q)-if [ -d lib ] ; then rm -r lib ; fi
	$(Q)-if [ -d mod ] ; then rm -r mod ; fi
