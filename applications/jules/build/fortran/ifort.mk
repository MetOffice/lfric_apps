##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the Intel Fortran compiler.
##############################################################################

$(info Project specials for Intel compiler)

export FFLAGS_UM_PHYSICS = -r8

# Set fast-debug and full-debug options for lfric_atm Intel that
# differ from the setup in infrastructure used by other LFRic apps.
# First, the -check all option is removed from all configurations as
# it results in excessive numbers of warnings from UM code running
# within LFRic
FFLAGS_RUNTIME            = -fpe0
# Second, certain fast-debug options cause XIOS failures on the Cray
# in jobs that write diagnostic so they are removed. Note: the
# full-debug test can still use these options as it avoids such XIOS
# use.
ifdef CRAY_ENVIRONMENT
# On the Cray these options are switched off for fast-debug
FFLAGS_FASTD_INIT               = 
FFLAGS_FASTD_RUNTIME            =
else
# Otherwise, use the same as the default full-debug settings
FFLAGS_FASTD_INIT         = $(FFLAGS_INIT) 
FFLAGS_FASTD_RUNTIME      = $(FFLAGS_RUNTIME)
endif

# NOTE: The -qoverride-limits option contained in $(FFLAGS_INTEL_FIX_ARG) is
# not currently applied here. This is a temporary workaround for #3465
# which it was found to be inadvertently preventing compilation
# openMP has also been removed from this routine via the optimisation script
%bl_imp_alg_mod_psy.o %bl_imp_alg_mod_psy.mod:   private FFLAGS_EXTRA =
%aerosol_ukca_alg_mod_psy.o %aerosol_ukca_alg_mod_psy.mod:   private FFLAGS_EXTRA =
%conv_comorph_alg_mod_psy.o %conv_comorph_alg_mod_psy.mod:   private FFLAGS_EXTRA =

$(info LFRic compile options required for files with OpenMP - see Ticket 1490)
%psy.o %psy.mod:   private FFLAGS_EXTRA = $(FFLAGS_INTEL_FIX_ARG)
# NOTE: The -qoverride-limits option contained in $(FFLAGS_INTEL_FIX_ARG) is
# not currently applied here. This is a temporary workaround for #3205
# which it was found to be inadvertently preventing compilation
# psy/%.o psy/%.mod: private FFLAGS_EXTRA = $(FFLAGS_INTEL_FIX_ARG)


# -warn noexternals applied to code that imports lfric_mpi_mod to avoid
# a warning-turned-error about missing interfaces for MPI calls in
# mpi.mod, such as MPI_Allreduce - switching to mpi_f08.mod resolves
# this via polymorphic interface declarations. Some SOCRATES functions
# do not currently declare interfaces either. Flag was introduced in
# Intel Fortran v19.1.0 according to Intel release notes.
ifeq ($(shell test "$(IFORT_VERSION)" -ge 0190100; echo $$?), 0)
  $(info ** Activating externals warning override for selected source files)
  export FFLAGS_INTEL_EXTERNALS = -warn noexternals
else
  export FFLAGS_INTEL_EXTERNALS =
endif

ifeq ($(shell test "$(IFORT_VERSION)" -ge 0190103; echo $$?), 0)
  $(info ** Disabling OpenMP due to a compiler bug for intel-compiler newer than 2020.3.304 - see Ticket 3853)
  %ls_ppnc.o %ls_ppnc.mod: FFLAGS:=${FFLAGS} -qno-openmp
endif
