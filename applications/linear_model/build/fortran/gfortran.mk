##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the GNU Fortran compiler.
##############################################################################

$(info Project specials for GNU compiler)

export FFLAGS_UM_PHYSICS = -fdefault-real-8 -fdefault-double-8
# Most lfric_atm dependencies contain code with implicit lossy conversions and
# unused variables.
# We reset the FFLAGS_WARNINGS variable here in order to prevent
# -Werror induced build failures.
FFLAGS_WARNINGS          = -Wall -Werror=character-truncation -Werror=unused-value \
                           -Werror=tabs
# But, we can apply some more lfric infrastructure checking to socrates
FFLAGS_SOCRATES_WARNINGS = -Werror=unused-variable

science/src/socrates/%.o science/src/socrates/%.mod: private FFLAGS_EXTRA = $(FFLAGS_SOCRATES_WARNINGS)

# We remove bounds checking (applied by -fcheck=all) and underflow checking. The
# latter is due to regular permitting of exponents going to zero for small numbers
# to imply total extinction of radiation passing through a medium
FFLAGS_RUNTIME           = -fcheck=all,no-bounds -ffpe-trap=invalid,zero,overflow

# The lfric_atm app defines an extra set of debug flags for
# fast-debug which are the same as the full-debug settings
# except for some platforms
ifdef CRAY_ENVIRONMENT
# On the EXZ these options are switched off for fast-debug
# due to an unexpected FPE in the NetCDF library
FFLAGS_FASTD_INIT               =
FFLAGS_FASTD_RUNTIME            =
else
# Otherwise, use the same as the default full-debug settings
FFLAGS_FASTD_INIT         = $(FFLAGS_INIT)
FFLAGS_FASTD_RUNTIME      = $(FFLAGS_RUNTIME)
endif
