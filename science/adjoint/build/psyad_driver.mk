##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Makes driver routine based on psyad files.
PSYAD_FILES_CORE := $(shell cat $(ADJOINT_BUILD)/psyad_files_list_core.txt)
PSYAD_FILES_CORE := $(addprefix $(CORE_ROOT_DIR)/,$(PSYAD_FILES_CORE))
PSYAD_FILES_APPS := $(shell cat $(ADJOINT_BUILD)/psyad_files_list_apps.txt)
PSYAD_FILES_APPS := $(addprefix $(APPS_ROOT_DIR)/,$(PSYAD_FILES_APPS))
PSYAD_FILES := $(PSYAD_FILES_CORE) $(PSYAD_FILES_APPS)

DRIVER_TEMPLATE_PATH := $(ADJOINT_BUILD)/gen_adj_kernel_tests_mod.txt
DRIVER_TARGET := $(WORKING_DIR)/driver/gen_adj_kernel_tests_mod.f90
DIRECTORIES := $(WORKING_DIR)/driver

all: $(DRIVER_TARGET)

$(DRIVER_TARGET): | $(DIRECTORIES)
	cp $(DRIVER_TEMPLATE_PATH) $@
	python $(ADJOINT_BUILD)/psyad_driver.py -f $(PSYAD_FILES) -d $@

$(DIRECTORIES):
	mkdir -p $@

