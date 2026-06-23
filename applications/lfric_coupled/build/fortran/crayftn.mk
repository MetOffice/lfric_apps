##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the Cray Fortran compiler.
##############################################################################

$(info Project specials for Cray compiler)

export FFLAGS_UM_PHYSICS = -s real64

#to try and ease compile time for CCE on EXZ
ifeq ($(shell expr ${CRAYFTN_VERSION} \>= 015000000), 1)
    %bl_imp_alg_mod_psy.o %bl_imp_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %bl_imp_alg_mod_psy.o %bl_imp_alg_mod_psy.mod: private FFLAGS_DEBUG = -G0
    %ukca_emiss_mode_mod.o %ukca_emiss_mode_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %ukca_emiss_mode_mod.o %ukca_emiss_mode_mod.mod: private FFLAGS_DEBUG = -G0
    %ukca_step_control_mod.o %ukca_step_control_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %ukca_step_control_mod.o %ukca_step_control_mod.mod: private FFLAGS_DEBUG = -G0
    %bl_exp_alg_mod_psy.o %bl_exp_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %bl_exp_alg_mod_psy.o %bl_exp_alg_mod_psy.mod: private FFLAGS_DEBUG = -G0
    %aerosol_ukca_alg_mod_psy.o %aerosol_ukca_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %aerosol_ukca_alg_mod_psy.o %aerosol_ukca_alg_mod_psy.mod: private FFLAGS_DEBUG = -G0
    %init_aerosol_fields_alg_mod_psy.o init_aerosol_fields_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %init_aerosol_fields_alg_mod_psy.o init_aerosol_fields_alg_mod_psy.mod: private FFLAGS_DEBUG = -G0
    %jules_extra_kernel_mod.o %jules_extra_kernel_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %jules_extra_kernel_mod.o %jules_extra_kernel_mod.mod: private FFLAGS_DEBUG = -G0
    %conv_gr_alg_mod_psy.o %conv_gr_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %conv_gr_alg_mod_psy.o %conv_gr_alg_mod_psy.mod: private FFLAGS_DEBUG = -G0
    %gungho_model_mod.o %gungho_model_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %gungho_model_mod.o %gungho_model_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_collections_mod.o %driver_collections_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_collections_mod.o %driver_collections_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_comm_mod.o %driver_comm_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_comm_mod.o %driver_comm_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_config_mod.o %driver_config_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_config_mod.o %driver_config_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_coordinates_mod.o %driver_coordinates_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_coordinates_mod.o %driver_coordinates_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_counter_mod.o %driver_counter_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_counter_mod.o %driver_counter_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_fem_mod.o %driver_fem_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_fem_mod.o %driver_fem_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_io_mod.o %driver_io_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_io_mod.o %driver_io_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_log_mod.o %driver_log_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_log_mod.o %driver_log_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_mesh_mod.o %driver_mesh_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_mesh_mod.o %driver_mesh_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_model_data_mod.o %driver_model_data_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_model_data_mod.o %driver_model_data_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_modeldb_mod.o %driver_modeldb_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_modeldb_mod.o %driver_modeldb_mod.mod: private FFLAGS_DEBUG = -G0
    %driver_time_mod.o %driver_time_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %driver_time_mod.o %driver_time_mod.mod: private FFLAGS_DEBUG = -G0
    %io_context_collection_mod.o %io_context_collection_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %io_context_collection_mod.o %io_context_collection_mod.mod: private FFLAGS_DEBUG = -G0
    %variable_fields_mod.o %variable_fields_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %variable_fields_mod.o %variable_fields_mod.mod: private FFLAGS_DEBUG = -G0
    %conv_comorph_kernel_mod.o %conv_comorph_kernel_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %conv_comorph_kernel_mod.o %conv_comorph_kernel_mod.mod: private FFLAGS_DEBUG = -G0
    %conv_comorph_alg_mod_psy.o %conv_comorph_alg_mod_psy.mod: private FFLAGS_SAFE_OPTIMISATION = -O0
    %conv_comorph_alg_mod_psy.o %conv_comorph_alg_mod_psy.mod: private FFLAGS_DEBUG = -G0
    %parcel_ascent_5a.o %parcel_ascent_5a.mod: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS) -h vector0
    # fast-debug options
    %qsat_mod.o %qsat_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O2 -hflex_mp=strict
    %glue_conv-6a.o %glue_conv-6a.mod: private FFLAGS_SAFE_OPTIMISATION = -O2 -hflex_mp=strict
    %sci_iterative_solver_mod.o %sci_iterative_solver_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O2 -hflex_mp=strict
    %semi_implicit_timestep_alg_mod.o % %semi_implicit_timestep_alg_mod.mod: private FFLAGS_SAFE_OPTIMISATION = -O2 -hflex_mp=strict
endif
