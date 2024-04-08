##############################################################################
# (c) Crown copyright 2021-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Master make file for LFRic um_physics_interface project.
# Targets provided are detailed below.
# TODO: Currently all targets are interactive only, #2757 will add the testing
#		to a rose-stem test-suite.
#
# all: (default) Complete build and test the um_physics_interface.
#      Intended for interactive use. Currently only builds and runs unit-tests.
# unit-tests: Run unit tests interactively.
# clean: Delete all final products and working files.
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
#
##############################################################################

PROJECT_NAME = um_physics_interface

PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = MPI pfunit_mod yaxt mod_oasis
export EXTERNAL_DYNAMIC_LIBRARIES = yaxt yaxt_c
export EXTERNAL_STATIC_LIBRARIES =

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

.PHONY: all
all: default

.PHONY: default
default: unit-tests
	$(Q)echo > /dev/null

META_VN       ?= HEAD

.PHONY: documentation doc docs
documentation doc docs: document-api
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk

##############################################################################
# Documentation
#
.PHONY: document-api
document-api: PROJECT       = um_physics_interface
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null

##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests/%: export ADDITIONAL_EXTRACTION = $(CORE_ROOT_DIR)/infrastructure/source \
                                             $(CORE_ROOT_DIR)/components/science/source \
                                             $(CORE_ROOT_DIR)/components/diagnostics_infrastructure/source \
                                             $(APPS_ROOT_DIR)/science/jules_interface/source/kernel \
                                             $(CORE_ROOT_DIR)/components/science/unit-test/support \
                                             $(APPS_ROOT_DIR)/science/gungho/source/diagnostics
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests/%: export FFLAG_GROUPS = -D UNIT_TEST
unit-tests/%: export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)/$(META_VN)
unit-tests/%: export PROGRAMS = um_physics_interface_unit_tests
unit-tests/%: export PROJECT = um_physics_interface
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: unit-tests/run
unit-tests/run: unit-tests/build
unit-tests/build: create_rose-meta

##############################################################################
# Create a consolidated rose-meta file for um_physics_interface, this must be performed
# before any build commands/rose-picker to create a single config file as
# required by rose-picker.
#
create_rose-meta: rose-meta/$(PROJECT_NAME)/HEAD/rose-meta.conf
rose-meta/$(PROJECT_NAME)/HEAD/rose-meta.conf:
	$(call MESSAGE, Creating, $@)
	$(Q)mkdir -p $(dir $@)
	$(Q)for meta_file in $(addprefix science/$(PROJECT_NAME)/,$(addsuffix $(META_VN), \
	$(dir $(wildcard rose-meta/*/.)))); \
	do echo import\=$${meta_file} >> $@; done

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"um_physics_interface work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"$(PROJECT_NAME) documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"um_physics_interface generated rose-meta file")
	$(Q)if [ -d rose-meta/$(PROJECT_NAME) ] ; then rm -r rose-meta/$(PROJECT_NAME); fi
	$(call MESSAGE,Removing,"um_physics_interface binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
