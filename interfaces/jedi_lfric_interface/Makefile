##############################################################################
# (c) Crown copyright 2020-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Make file for jedi_lfric_interface model component
# Targets provided our detailed below...
#
# all: (default) Complete build and test the jedi_lfric_interface model component.

# test: Run test battery including unit tests and others
# clean: Delete all final products and working files
#
# The following variables may be specified to modify the build process...
#
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
# PROFILE: Set to a string representing a package of compiler options.
#          Potential profiles are 'full-debug', 'fast-debug' and 'production'.
#          Default: fast-debug
# LINK_TYPE: Either 'static' or 'dynamic'.
#            Default: dynamic
#
##############################################################################

PROJECT_NAME = jedi_lfric_interface

export PROFILE ?= fast-debug

export EXTERNAL_DYNAMIC_LIBRARIES = netcdff netcdf yaxt yaxt_c

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(APPS_ROOT_DIR)/science/gungho \
                               $(CORE_ROOT_DIR)/components/lfric-xios \
                               $(APPS_ROOT_DIR)/science/adjoint

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: unit-tests integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-api
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(PROJECT_DIR)/build/import.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)

##############################################################################
# Documentation
#
.PHONY: document-uml
document-uml: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/uml
document-uml: SOURCE_DIR    = documentation/uml
document-uml: WORKING_DIR  := $(WORKING_DIR)/uml
document-uml: uml-documentation
	$(Q)echo > /dev/null

.PHONY: document-api
document-api: PROJECT       = jedi_lfric_interface
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null


##############################################################################
# Unit tests
#
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export CXX_LINK = YES
unit-tests/%: export FFLAGS += -D UNIT_TEST
unit-tests/%: export IMPORT_PARTS = $(CORE_ROOT_DIR)/infrastructure \
                                    $(CORE_ROOT_DIR)/components/science \
                                    $(CORE_ROOT_DIR)/components/lfric-xios \
                                    $(CORE_ROOT_DIR)/components/driver \
                                    $(APPS_ROOT_DIR)/science/gungho
unit-tests/%: export META_FILE_DIR = rose-meta/jedi_lfric_interface/HEAD
unit-tests/%: export PROGRAMS = jedi_lfric_unit_tests
unit-tests/%: export PROJECT = jedi_lfric
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: unit-tests/run


##############################################################################
# Integration tests
#
integration-tests/%: export BIN_DIR      ?= $(PROJECT_DIR)/test
integration-tests/%: export FFLAGS       += -DINT_TEST
integration-tests/%: export IMPORT_PARTS  = $(CORE_ROOT_DIR)/infrastructure \
                                            $(CORE_ROOT_DIR)/components/science \
                                            $(APPS_ROOT_DIR)/science/gungho
integration-tests/%: export META_FILE_DIR = rose-meta/jedi_lfric_interface/HEAD
integration-tests/%: export PROGRAMS = jedi_lfric_integration_tests
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
	$(call MESSAGE,Removing,"jedi_lfric_interface work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"jedi_lfric_interface documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"jedi_lfric_interface binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
