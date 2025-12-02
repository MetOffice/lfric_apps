##############################################################################
# (c) Crown copyright 2021-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Master make file for LFRic physics_schemes_interface project.
# Targets provided are detailed below.
# TODO: Currently all targets are interactive only, #2757 will add the testing
#		to a rose-stem test-suite.
#
# all: (default) Complete build and test the physics_schemes_interface.
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

PROJECT_NAME = physics_schemes_interface

export PROFILE ?= fast-debug

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
document-api: PROJECT       = physics_schemes_interface
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
unit-tests/%: export ADDITIONAL_EXTRACTION = $(CORE_ROOT_DIR)/components/science/unit-test/support
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests/%: export FFLAG_GROUPS = -D UNIT_TEST
unit-tests/%: export PRE_PROCESS_INCLUDE_DIRS = \
        $(WORKING_DIR)/atmosphere_service/include \
        $(WORKING_DIR)/boundary_layer/include \
		$(WORKING_DIR)/free_tracers/include \
        $(WORKING_DIR)/large_scale_precipitation/include

unit-tests/%: export PRE_PROCESS_MACROS += UM_PHYSICS LFRIC USSPPREC_32B LSPREC_32B UM_JULES
unit-tests/%: export IMPORT_PARTS = $(CORE_ROOT_DIR)/infrastructure \
                                    $(CORE_ROOT_DIR)/components/science \
                                    $(APPS_ROOT_DIR)/interfaces/jules_interface \
                                    $(APPS_ROOT_DIR)/science/gungho \
                                    $(APPS_ROOT_DIR)/interfaces/socrates_interface \
                                    $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface
unit-tests/%: export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)/$(META_VN)
unit-tests/%: export PROGRAMS = physics_schemes_interface_unit_tests
unit-tests/%: export PROJECT = physics_schemes_interface
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export SCRATCH_DIR = $(WORKING_DIR)/../scratch
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: unit-tests/run
unit-tests/run: unit-tests/build
unit-tests/build: create_rose-meta

##############################################################################
# Create a consolidated rose-meta file for physics_schemes_interface, this must be performed
# before any build commands/rose-picker to create a single config file as
# required by rose-picker.
#
create_rose-meta: rose-meta/$(PROJECT_NAME)/HEAD/rose-meta.conf
rose-meta/$(PROJECT_NAME)/HEAD/rose-meta.conf:
	$(call MESSAGE, Creating, $@)
	$(Q)mkdir -p $(dir $@)
	${Q}for meta_file in ./rose-meta/* ; do \
		  if [[ $${meta_file} != *$(PROJECT_NAME)* ]]; then \
	        echo import\=$${meta_file##*/}/${META_VN} >> $@ ; \
		  fi \
	    done


##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"physics_schemes_interface work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"$(PROJECT_NAME) documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"physics_schemes_interface generated rose-meta file")
	$(Q)if [ -d rose-meta/$(PROJECT_NAME) ] ; then rm -r rose-meta/$(PROJECT_NAME); fi
	$(call MESSAGE,Removing,"physics_schemes_interface binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
