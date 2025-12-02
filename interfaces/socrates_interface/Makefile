##############################################################################
# (c) Crown copyright 2021-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Make file for socrates_interface
# Targets provided our detailed below.
#
# all: (default) Complete build and test the socrates_interface.

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
##############################################################################

PROJECT_NAME = socrates_interface

export PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = yaxt yaxt_c netcdff netcdf hdf5

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)-radiation/$(META_VN)

.PHONY: default
default: unit-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-api
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk

##############################################################################
# Documentation
#
.PHONY: document-api
document-api: PROJECT       = socrates_interface
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null

##############################################################################
# Unit tests
#
unit-tests/%: export ADDITIONAL_EXTRACTION = $(CORE_ROOT_DIR)/components/science/unit-test/support
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export FFLAGS += -D UNIT_TEST
unit-tests/%: export IMPORT_PARTS = $(CORE_ROOT_DIR)/infrastructure \
                                    $(CORE_ROOT_DIR)/components/science
unit-tests/%: export META_FILE_DIR = rose-meta/socrates-radiation/HEAD
unit-tests/%: export PROGRAMS = socrates_interface_unit_tests
unit-tests/%: export PROJECT = socrates_interface
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export SCRATCH_DIR = $(WORKING_DIR)/../scratch
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: unit-tests/run

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"socrates_interface work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"socrates_interface documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"socrates_interface binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
