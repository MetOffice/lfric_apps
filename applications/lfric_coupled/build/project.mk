##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

$(info Coupling project configuration)

export PRE_PROCESS_INCLUDE_DIRS = \
        $(WORKING_DIR)/atmosphere_service/include \
        $(WORKING_DIR)/boundary_layer/include \
        $(WORKING_DIR)/large_scale_precipitation/include \
        $(WORKING_DIR)/free_tracers/include

export PRE_PROCESS_MACROS += UM_PHYSICS LFRIC USSPPREC_32B LSPREC_32B UM_JULES

export PRE_PROCESS_MACROS += COUPLED
