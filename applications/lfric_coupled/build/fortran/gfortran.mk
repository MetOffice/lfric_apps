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

# But, we can apply full lfric infrastructure checking to socrates
FFLAGS_SOCRATES_WARNINGS = -Werror=conversion -Werror=unused-variable

science/src/socrates/%.o science/src/socrates/%.mod: private export FFLAGS_EXTRA += $(FFLAGS_SOCRATES_WARNINGS)

# We remove bounds checking (applied by -fcheck=all) and underflow checking. The 
# latter is due to regular permitting of exponents going to zero for small numbers
# to imply total extinction of radiation passing through a medium
FFLAGS_RUNTIME           = -fcheck=all,no-bounds -ffpe-trap=invalid,zero,overflow
