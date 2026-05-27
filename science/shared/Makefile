###############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
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

PROJECT_NAME = shared

export PROFILE ?= fast-debug

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/science \
                               $(CORE_ROOT_DIR)/components/inventory \
                               $(CORE_ROOT_DIR)/components/lfric-xios \
                               $(APPS_ROOT_DIR)/science/gungho

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: unit-tests
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)

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
                                     $(CORE_ROOT_DIR)/components/inventory \
                                     $(APPS_ROOT_DIR)/science/gungho
unit-tests/%: export META_FILE_DIR = rose-meta/lfric-shared/HEAD
unit-tests/%: export PRE_PROCESS_MACROS = 'UNIT_TEST=TRUE'
unit-tests/%: export PROGRAMS = shared_unit_tests
unit-tests/%: export PROJECT = shared
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: unit-tests/run

###############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"shared work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"shared documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"shared binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
