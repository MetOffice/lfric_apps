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

PROJECT_NAME = physics_schemes

export PROFILE ?= fast-debug

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

COMPONENT_PREREQUISITES = $(CORE_ROOT_DIR)/infrastructure \
                          $(CORE_ROOT_DIR)/components/science \
                          $(APPS_ROOT_DIR)/interfaces/jules_interface \
                          $(APPS_ROOT_DIR)/science/shared

.PHONY: default
default: unit-tests

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(COMPONENT_PREREQUISITES:=/build/import.mk)
include $(PROJECT_DIR)/build/import.mk

###############################################################################
# Unit tests
#
unit-tests/%: export ADDITIONAL_EXTRACTION = 
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export IMPORT_PARTS  = $(COMPONENT_PREREQUISITES)
unit-tests/%: export PROGRAMS = physics_schemes_unit_tests
unit-tests/%: export PROJECT = $(PROJECT_NAME)
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: unit-tests/run

###############################################################################
# Integration tests - commented out until we have some
#
#integration-tests/%: export ADDITIONAL_EXTRACTION =
#integration-tests/%: export BIN_DIR      ?= $(PROJECT_DIR)/test
#integration-tests/%: export IMPORT_PARTS  = 
#integration-tests/%: export PROJECT       = physics_schemes_integration_tests
#integration-tests/%: export SOURCE_DIR    = source
#integration-tests/%: export TEST_DIR      = integration-test
#integration-tests/%: export WORKING_DIR  := $(WORKING_DIR)/integration-tests
#integration-tests: integration-tests/run

###############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"work space")
	$Q-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"binaires")
	$Q-rm -r test
