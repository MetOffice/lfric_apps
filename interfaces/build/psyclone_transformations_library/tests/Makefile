.PHONY: test
test: no-optimisation/confirm \
      global-optimisation/confirm \
      local-optimisation/confirm \
      both-optimisation/confirm \
      override/confirm \
      override-global-opt/confirm \
      override-local-opt/confirm \
      override-both-opt/confirm

# Tests operation without transformation scripts. Both preprocessed and
# non-preprocessed algorithms are present.
#
no-optimisation/%: export SOURCE_DIR = source/automatic
no-optimisation/%: export WORKING_DIR = working/no-optimisation

no-optimisation/confirm: no-optimisation/invoke
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/preprocessed_alg.f90
	test $(shell grep -c 'one_proxy%data(df) = 1.0' $(WORKING_DIR)/preprocessed_alg_psy.f90) = 1
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/sub/unprocessed_alg.f90
	test $(shell grep -c 'one_proxy%data(df) = 2.0' $(WORKING_DIR)/sub/unprocessed_alg_psy.f90) = 1

no-optimisation/invoke:
	$(info *** No optimisation)
	$(MAKE) -f ../psyclone_psykal.mk

# Tests operation with a global transformation script. Both preprocessed and
# non-preprocessed algorithms are present.
#
global-optimisation/%: export SOURCE_DIR = source/automatic
global-optimisation/%: export WORKING_DIR = working/global-optimisation
global-optimisation/%: export OPTIMISATION_PATH = optimisation/global

global-optimisation/confirm: global-optimisation/invoke
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/preprocessed_alg.f90
	test $(shell grep -c 'one_proxy%data(df) = 1.0' $(WORKING_DIR)/preprocessed_alg_psy.f90) = 2
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/sub/unprocessed_alg.f90
	test $(shell grep -c 'one_proxy%data(df) = 2.0' $(WORKING_DIR)/sub/unprocessed_alg_psy.f90) = 2

global-optimisation/invoke:
	$(info *** Global optimisation)
	$(MAKE) -f ../psyclone_psykal.mk

# Tests operation with a transformation script for only one of the algorithms.
# Both preprocessed and non-preprocessed algorithms are present.
#
local-optimisation/%: export SOURCE_DIR = source/automatic
local-optimisation/%: export WORKING_DIR = working/local-optimisation
local-optimisation/%: export OPTIMISATION_PATH = optimisation/local

local-optimisation/confirm: local-optimisation/invoke
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/preprocessed_alg.f90
	test $(shell grep -c 'one_proxy%data(df) = 1.0' $(WORKING_DIR)/preprocessed_alg_psy.f90) = 2
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/sub/unprocessed_alg.f90
	test $(shell grep -c 'one_proxy%data(df) = 2.0' $(WORKING_DIR)/sub/unprocessed_alg_psy.f90) = 1

local-optimisation/invoke:
	$(info *** Local optimisation)
	$(MAKE) -f ../psyclone_psykal.mk

# Tests with both global and local transformation scripts. Both preprocessed
# and non-preprocessed algorithms are present.
#
both-optimisation/%: export SOURCE_DIR = source/automatic
both-optimisation/%: export WORKING_DIR = working/both-optimisation
both-optimisation/%: export OPTIMISATION_PATH = optimisation/both

both-optimisation/confirm: both-optimisation/invoke
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/preprocessed_alg.f90
	test $(shell grep -c 'one_proxy%data(df) = 1.0' $(WORKING_DIR)/preprocessed_alg_psy.f90) = 3
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/sub/unprocessed_alg.f90
	test $(shell grep -c 'one_proxy%data(df) = 2.0' $(WORKING_DIR)/sub/unprocessed_alg_psy.f90) = 2

both-optimisation/invoke:
	$(info *** Both optimisations)
	$(MAKE) -f ../psyclone_psykal.mk

# Tests operation with a handwriten PSy source.
#
override/%: export SOURCE_DIR = source/handwritten
override/%: export WORKING_DIR = working/override

override/confirm: override/invoke
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/unprocessed_alg.f90
	test ! -e $(WORKING_DIR)/preprocessed_alg_psy.f90
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/sub/preprocessed_alg.f90
	test ! -e $(WORKING_DIR)/sub/unprocessed_alg_psy.f90

override/invoke:
	$(info *** Override PSy layer)
	$(MAKE) -f ../psyclone_psykal.mk

# Tests operation with a handwriten PSy source and global optimisation.
#
override-global-opt/%: export SOURCE_DIR = source/handwritten
override-global-opt/%: export WORKING_DIR = working/override
override-global-opt/%: export OPTIMISATION_PATH = optimisation/global

override-global-opt/confirm: override-global-opt/invoke
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/unprocessed_alg.f90
	test ! -e $(WORKING_DIR)/preprocessed_alg_psy.f90
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/sub/preprocessed_alg.f90
	test ! -e $(WORKING_DIR)/sub/unprocessed_alg_psy.f90

override-global-opt/invoke:
	$(info *** Override PSy layer)
	$(MAKE) -f ../psyclone_psykal.mk

# Tests operation with a handwriten PSy source and local optimisation.
#
override-local-opt/%: export SOURCE_DIR = source/handwritten
override-local-opt/%: export WORKING_DIR = working/override
override-local-opt/%: export OPTIMISATION_PATH = optimisation/local

override-local-opt/confirm: override-local-opt/invoke
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/unprocessed_alg.f90
	test ! -e $(WORKING_DIR)/preprocessed_alg_psy.f90
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/sub/preprocessed_alg.f90
	test ! -e $(WORKING_DIR)/sub/unprocessed_alg_psy.f90

override-local-opt/invoke:
	$(info *** Override PSy layer)
	$(MAKE) -f ../psyclone_psykal.mk

# Tests operation with a handwriten PSy source and both global and local 
# optimisation.
#
override-both-opt/%: export SOURCE_DIR = source/handwritten
override-both-opt/%: export WORKING_DIR = working/override
override-both-opt/%: export OPTIMISATION_PATH = optimisation/both

override-both-opt/confirm: override-both-opt/invoke
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/unprocessed_alg.f90
	test ! -e $(WORKING_DIR)/preprocessed_alg_psy.f90
	grep 'CALL invoke_0(one)' $(WORKING_DIR)/sub/preprocessed_alg.f90
	test ! -e $(WORKING_DIR)/sub/unprocessed_alg_psy.f90

override-both-opt/invoke:
	$(info *** Override PSy layer)
	$(MAKE) -f ../psyclone_psykal.mk

# Clean up after ourselves.
#
.PHONY: clean
clean:
	-rm -r working

# Perform style checking on test files
#
.PHONY: style
style:
	pycodestyle .
	pylint .

# Bring in the basic bits and pieces.
#
include ../../lfric.mk
