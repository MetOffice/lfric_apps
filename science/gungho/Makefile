###############################################################################
# (c) Crown copyright 2021-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
#
# Make file for science model component
# Targets provided our detailed below.
#
# all: (default) Complete build and test the science model component.

# test: Run test battery including unit tests and others.
# clean: Delete all final products and working files.
#
# The following variables may be specified to modify the build process:
#
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
# PROFILE: Set to a string representing a package of compiler options.
#          Potential profiles are 'full-debug', 'fast-debug' and 'production'.
#          Default: 'fast-debug'
# LINK_TYPE: Either 'static' or 'dynamic'.
#            Default: 'dynamic'
#
###############################################################################

PROJECT_NAME = gungho

export PROFILE ?= fast-debug

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/driver \
                               $(CORE_ROOT_DIR)/components/science \
                               $(CORE_ROOT_DIR)/components/inventory \
                               $(CORE_ROOT_DIR)/components/lfric-xios

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: unit-tests integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-api
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)

###############################################################################
# Documentation
#
.PHONY: document-api
document-api: PROJECT       = gungho
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation

.PHONY: document-uml
document-uml: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/uml
document-uml: SOURCE_DIR    = documentation/uml
document-uml: WORKING_DIR  := $(WORKING_DIR)/uml
document-uml: uml-documentation

.PHONY: document-latex
document-latex: DOCUMENT_DIR  ?= $(PROJECT_DIR)/documents
document-latex: SOURCE_DIR     = documentation
document-latex: DOCUMENTS      = $(shell find $(SOURCE_DIR) -name '*.latex' -print)
document-latex: WORKING_DIR   := $(WORKING_DIR)/latex
document-latex: TEX_STUFF      = documentation/tex
document-latex: COMMON_FIGURES = documentation/common-figures
document-latex: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/latex.mk \
	                         SOURCE_DIR=$(SOURCE_DIR) \
	                         WORKING_DIR=$(WORKING_DIR) \
	                         DOCUMENT_DIR=$(DOCUMENT_DIR) \
	                         TEX_STUFF=$(TEX_STUFF) \
	                         COMMON_FIGURES=$(COMMON_FIGURES) \
	                         DOCUMENTS="$(DOCUMENTS)"

###############################################################################
# Unit tests
#
unit-tests/%: export ADDITIONAL_EXTRACTION = \
                          $(CORE_ROOT_DIR)/components/science/unit-test/support
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export CXX_LINK = TRUE
unit-tests/%: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests/%: export IMPORT_PARTS  = $(CORE_ROOT_DIR)/infrastructure \
                                     $(CORE_ROOT_DIR)/components/science \
                                     $(CORE_ROOT_DIR)/components/lfric-xios \
                                     $(CORE_ROOT_DIR)/components/inventory
unit-tests/%: export META_FILE_DIR = rose-meta/lfric-gungho/HEAD
unit-tests/%: export PRE_PROCESS_MACROS = 'UNIT_TEST=TRUE'
unit-tests/%: export PROGRAMS = gungho_unit_tests
unit-tests/%: export PROJECT = gungho
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: unit-tests/run

###############################################################################
# Integration tests
#
integration-tests/%: export ADDITIONAL_EXTRACTION = \
                                integration-test/support/tiny_world.nc
integration-tests/%: export BIN_DIR      ?= $(PROJECT_DIR)/test
integration-tests/%: export IMPORT_PARTS  = $(INTERNAL_DEPENDENCIES)
integration-tests/%: export META_FILE_DIR = rose-meta/lfric-gungho/HEAD
integration-tests/%: export PROJECT       = gungho
integration-tests/%: export SOURCE_DIR    = source
integration-tests/%: export TEST_DIR      = integration-test
integration-tests/%: export WORKING_DIR  := $(WORKING_DIR)/integration-tests
integration-tests: integration-tests/run

###############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"Gung Ho work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"Gung Ho documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"Gung Ho binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
