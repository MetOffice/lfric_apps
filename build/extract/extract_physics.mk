##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Run this file to extract source code from the UM repository.
#
# The following environment variables are used for input:
#   UM_FCM_TARGET_PLATFORM : Target identifier used to get the
#                            correct UM build configs.
#   PROFILE : Build profile used to determine optimisation level.
#   PROJECT_DIR : Full path to the current project's root directory.
#   SCRATCH_DIR : Temporary space for extracted source.
#   WORKING_DIR : Directory to hold working copies of source.
#
###############################################################################

.PHONY: extract

extract:
	# Retrieve and preprocess the UM, Jules and Socrates code
	# The UM_ENV file contains the appropriate locations and UM side
	# environment variables
	$Q. $(APPS_ROOT_DIR)/dependencies.sh \
	   && fcm make -C $(SCRATCH_DIR) -f $(APPS_ROOT_DIR)/build/extract/extract.cfg
	# Note that if wanting to modify UM source this should be done via the
	# UM repository either through a working copy or branch
	$Qrsync -acvz $(SCRATCH_DIR)/extract/ $(WORKING_DIR)/science/

	# preprocess and psyclone files that are required.
	# For forked code this occurs in the um_physics_interface
	# For FCM checked out code, it occurs here.
	# Note depending on site target requirements, either all source or specific source will be affected
	# See application overrides and files exported there
	#
	echo "Extract physics pre-processing and psyclone"

	$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/build/pre_process_phys.mk \
		SOURCE_DIR=$(WORKING_DIR) \
		PSYCLONE_PHYSICS_FILES="$(PSYCLONE_PHYSICS_FILES_FCM)" \
		PSYCLONE_DIRECTORIES="$(PSYCLONE_DIRECTORIES_FCM)" \
		PSYCLONE_PHYSICS_EXCEPTION="$(PSYCLONE_PHYSICS_EXCEPTION_FCM)"

	$(MAKE) $(QUIET_ARG) -f $(APPS_ROOT_DIR)/interfaces/physics_schemes_interface/build/psyclone_transmute.mk \
		SOURCE_DIR=$(WORKING_DIR) \
		OPTIMISATION_PATH_PSY=$(APPS_ROOT_DIR)/applications/$(PROJECT)/$(OPTIMISATION_PATH) \
		PSYCLONE_PHYSICS_FILES="$(PSYCLONE_PHYSICS_FILES_FCM)" \
		PSYCLONE_DIRECTORIES="$(PSYCLONE_DIRECTORIES_FCM)" \
		PSYCLONE_PHYSICS_EXCEPTION="$(PSYCLONE_PHYSICS_EXCEPTION_FCM)"

