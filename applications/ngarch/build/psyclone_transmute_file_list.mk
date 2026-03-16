##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# File lists provided will use the transmute PSyclone method.
# https://code.metoffice.gov.uk/trac/lfric_apps/ticket/724


##### TRANSMUTE_INCLUDE_METHOD specify_include #####
# For CPU OMP, we want to choose which files get run through PSyclone,
# and preserve existing hand coded optimisations.

# Choose which files to Pre-proccess and PSyclone from physics / FCM UM source.

export PSYCLONE_PHYSICS_FILES = \
                                bl_lsp \
                                bl_diags_mod \
                                bm_tau_kernel_mod \
                                btq_int \
                                conv_gr_kernel_mod \
                                ex_flux_tq \
                                ex_flux_uv \
                                fm_drag \
                                gw_ussp_mod \
                                imp_mix \
                                jules_exp_kernel_mod \
				                jules_extra_kernel_mod \
                                jules_imp_kernel_mod \
                                kmkh \
                                kmkhz_9c_wtrac \
                                lw_kernel_mod \
                                mphys_kernel_mod \
                                pc2_initiation_kernel_mod \
                                pc2_bl_forced_cu \
                                pc2_bm_initiate \
                                pc2_initiation_ctl \
                                pc2_conv_coupling_kernel_mod \
                                sw_kernel_mod \
                                sw_rad_tile_kernel_mod \
                                tr_mix \
	                            ukca_aero_ctl \
                                ukca_chemistry_ctl_full_mod \
                                ukca_main1-ukca_main1


##### TRANSMUTE_INCLUDE_METHOD specify_include #####

# List to use PSyclone explicitly without any opt script
# This will remove hand written (OMP) directives in the source
# Used by both methods, specify_include and specify_exclude
export PSYCLONE_PASS_NO_SCRIPT = ukca_abdulrazzak_ghan

##### TRANSMUTE_INCLUDE_METHOD specify_exclude #####
# For GPU, we may want to use more generic local.py transformation scripts and psyclone by directory.
# Advise which directories to pass to PSyclone.
# All files in these directories will be run through PSyclone using the transmute method.
# Also provide an optional exception list.
# These files will be filtered, and will NOT be run through PSyclone.

# Directories to psyclone
export PSYCLONE_DIRECTORIES =

# A general file exception list
export PSYCLONE_PHYSICS_EXCEPTION =

##### TRANSMUTE_INCLUDE_METHOD specify_exclude #####
