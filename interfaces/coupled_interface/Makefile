##############################################################################
# (c) Crown copyright 2017-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_NAME = coupling
export PROJECT_DIR_NAME = interfaces/coupled_interface
export PROFILE ?= fast-debug

# Select option for PSyclone transformation. Select "none" for no transformation
export PSYCLONE_TRANSFORMATION ?= none

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/lfric-xios

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)/$(META_VN)

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(PROJECT_DIR)/build/import.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)

##############################################################################
# Check PSyclone transformation option is valid
ifneq ("$(PSYCLONE_TRANSFORMATION)","none")
  export OPTIMISATION_PATH ?= $(wildcard optimisation/$(PSYCLONE_TRANSFORMATION))
  ifeq ("$(wildcard $(OPTIMISATION_PATH)/global.py)","")
    $(error PSyclone transformation dir "$(PSYCLONE_TRANSFORMATION)" does not have a global.py script)
  endif
endif

.PHONY: default
default: unit-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-latex document-api
	$(Q)echo > /dev/null

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
# Unit tests
#
unit-tests/%: export ADDITIONAL_EXTRACTION = \
                     $(CORE_ROOT_DIR)/components/science/unit-test/support
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
#
# It should not be necessary to explicitly link the pfunit library as it is
# done for you by the build system. This just adds a second "-lpfunit"
# argument straight after that one. However it doesn't link without it and
# the current ticket needs to go on sooner rather than later.
#
unit-tests/%: export EXTERNAL_STATIC_LIBRARIES = pfunit
unit-tests/%: export IMPORT_PARTS = $(CORE_ROOT_DIR)/infrastructure \
                                    $(CORE_ROOT_DIR)/components/science
unit-tests/%: export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)/$(META_VN)
unit-tests/%: export PROGRAMS = $(PROJECT_NAME)_unit_tests
unit-tests/%: export PROJECT = $(PROJECT_NAME)
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-test
unit-tests: unit-tests/run


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
