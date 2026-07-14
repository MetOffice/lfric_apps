##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the Intel Fortran compiler.
##############################################################################

$(info Project specials for Intel compiler)

export FFLAGS_UM_PHYSICS = -r8
# Remove -check-all option as it causes very slow runs due to a lot of array
# temporary warnings caused by UM code 
FFLAGS_RUNTIME            = -fpe0
FFLAGS_FORTRAN_STANDARD   = 

# NOTE: The -qoverride-limits option contained in $(FFLAGS_INTEL_FIX_ARG) is
# not currently applied here. This is a temporary workaround for #3465
# which it was found to be inadvertently preventing compilation
# openMP has also been removed from this routine via the optimisation script
%bl_imp_alg_mod_psy.o %bl_imp_alg_mod_psy.mod:   private FFLAGS_EXTRA =
%aerosol_ukca_alg_mod_psy.o %aerosol_ukca_alg_mod_psy.mod:   private FFLAGS_EXTRA =
%conv_comorph_alg_mod_psy.o %conv_comorph_alg_mod_psy.mod:   private FFLAGS_EXTRA =

$(info LFRic compile options required for files with OpenMP - see Ticket 1490)
%psy.o %psy.mod:   private export FFLAGS_EXTRA += $(FFLAGS_INTEL_FIX_ARG)
# NOTE: The -qoverride-limits option contained in $(FFLAGS_INTEL_FIX_ARG) is
# not currently applied here. This is a temporary workaround for #3205
# which it was found to be inadvertently preventing compilation
# psy/%.o psy/%.mod: private export FFLAGS_EXTRA += $(FFLAGS_INTEL_FIX_ARG)
