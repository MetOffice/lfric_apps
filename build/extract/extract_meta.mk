##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Run this file to extract rose-metadata from external repositories
#
# The following environment variables are used for input:
#   WORKING_DIR : Directory to hold working copies of source.
#   META_DIR: Directory to save extracted external metadata
#   APPS_ROOT_DIR: Location of apps clone
#   EXTRA_ROSE_META: List of repos to extract metadata from
#
###############################################################################

.PHONY: extract_meta

extract_meta:
	$(Q)for REPO in $(EXTRA_ROSE_META) ; do \
		python $(APPS_ROOT_DIR)/build/extract/extract_science.py \
			-r $$REPO \
			-d $(APPS_ROOT_DIR)/dependencies.yaml \
			-w $(WORKING_DIR) \
			-m $(META_DIR) \
		; done
