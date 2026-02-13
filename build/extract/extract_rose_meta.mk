##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Run this file to extract source code from the UM repository.
#
# The following environment variables are used for input:
#   REPO : Repo to extract rose-meta from.
#   WORKING_DIR : Directory to hold working copies of source.
#
###############################################################################

.PHONY: extract_meta

extract_meta:
	# Extract rose-meta from REPO
	python $(APPS_ROOT_DIR)/build/extract/extract_rose_meta.py -r $(REPO) -d $(APPS_ROOT_DIR)/dependencies.yaml -o $(WORKING_DIR)/../rose-meta
