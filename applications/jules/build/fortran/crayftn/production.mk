##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
##############################################################################
# Various things specific to production lfric_atm when using the
# Cray Fortran compiler.
##############################################################################

# ==========================================================================
# DIRECTORIES
# ==========================================================================
gravity_wave_drag/%.o: private FFLAGS_RISKY_OPTIMISATION = -O2 -hflex_mp=strict

# ==========================================================================
# MODULES
# ==========================================================================
# UKCA
%ukca_emiss_mode_mod.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%ukca_step_control_mod.o: private FFLAGS_RISKY_OPTIMISATION = -O0

# UM
%parcel_ascent_5a.o: private FFLAGS_EXTRA = $(FFLAGS_UM_PHYSICS) -hvector0

# LFRic Apps
%aerosol_ukca_alg_mod_psy.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%bl_exp_alg_mod_psy.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%bl_imp_alg_mod_psy.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%conv_comorph_alg_mod_psy.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%conv_comorph_kernel_mod.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%conv_gr_alg_mod_psy.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%gungho_model_mod.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%init_aerosol_fields_alg_mod_psy.o: private FFLAGS_RISKY_OPTIMISATION = -O0
%jules_extra_kernel_mod.o: private FFLAGS_RISKY_OPTIMISATION = -O0
