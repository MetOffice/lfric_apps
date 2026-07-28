! *********************************COPYRIGHT************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *********************************COPYRIGHT************************************
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_vera
!
! Part of the Vera scheme for diagnosing visual range.
!
! This module contains the Vera aerosol hydration scheme, vera_hydrate.
!
! This hydration scheme uses a set of Kohler curves to define the
! equilibrium state of an ensemble population of aerosol particles.
!
! The set of coupled equations defining the Vera hydration scheme,
! F(x)=0, are solved using the MINPACK solver. This scheme is designed
! so that just the function F(x)=0 needs be defined, in which case
! the MINPACK subroutine MinpackHybrDriver forms a numerical
! approximation to the Jacobian J(x) of the system F(x)=0.
!
! Or, the analytical form of the  J(x) can be used by the subroutine
! MinpackHybrDriver in MINPACK.
!
! The definitions of F(x)=0 and J(x) are provided in the module
! vera_function_mod.F90
!
! For details of the MINPACK solver routines, look at
!
!   "User Guide for MINPACK-1"
!    Jorge J. More, Burton S. Garbow, Kenneth E. Hillstrom, 1980.
!

module vera_hydration_mod

use vera_kind_mod,     only: wp => vera_real, wi => vera_integer ,             &
                             mr => minpack_r, mi => minpack_i

use vera_global_mod,   only: vera_phantom ,                                    &
                             vera_minpack

use vera_minpack_mod,  only: MinpackUserFunction ,                             &
                             MinpackDriver       ,                             &
                             minpack_parameters

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

! Description:
!   This module contains the Vera aerosol hydration scheme, and relies on the
!   third party software MINPACK to solve a set of nonlinear equations that
!   describe the equilibrium conditions for hydrated aerosol particles.
!
! Method:
!   This module comprises one subroutine:
!
!     vera_hydrate
!       Hydrates the aerosol particles using the MINPACK
!       solver using either a numerical approximation
!       to compute the Jacobian J(x) of the system of
!       defining equations F(x)=0, or an analytical expression
!       to compute the Jacobian J(x) of the system of
!       defining equations F(x)=0.
!
!       Which Jacobian computation method is used is determined by a
!       switch in the Vera glogal data module,
!       vera_minpack % solver_flag, as follows:
!
!       vera_minpack % solver_flag ==
!
!         vera_minpack % solver_flag_numeric
!         - use numerical approximation
!
!         vera_minpack % solver_flag_analytic
!         - use the explicit analytical expression. This is the default.
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_HYDRATION_MOD'

private

! make the subroutine in this module Public,
! it's called by vera_scheme_mod.F90
public :: vera_hydrate

contains

  !=============================================================================
  !
  ! vera_hydrate - Hydrates the aerosol particles using the MINPACK
  !                solver MinpackHybrDriver.
  !
  !=============================================================================

subroutine vera_hydrate( x_values, VeraFunction )

implicit none

!===========================================================================
! arguments for vera_hydrate
!
! NOTE: the type declarations for these variables are
!       either integer (kind = 4) or REAl (kind = 8) to
!       match with the interface to the MINPACK solver.
!       These Kinds are defined in the module vera_kind_mod as
!       mi and mr
!===========================================================================

! function to use to define the hydration scheme
procedure(MinpackUserFunction):: VeraFunction

! on input x_values are the initial guess as to the solution,
! and on output, x_values contains the computed solution
real (wp), intent(in out)     :: x_values(:)

!==========================================================================
! local variables for vera_hydrate
!
!                   all these are used to interface with the
!                   MINPACK solver, hence they have specific
!                   type declarations, either mr or mi
!==========================================================================

! number of variables in the system to solve, i.e. number of aerosol
! particle species
integer (mi)                   :: n

! the variable values to use for the MINPACK solver
real    (mr)                   :: x_solver(size(x_values))

! used by the MINPACK solver to report back the final value
! of func(x_values)
real    (mr)                   :: fvec(size(x_values))

! Jacobian of function F
real    (mr)                   :: j_x(1:size(x_values), 1:size(x_values))

! leading dimension of Jacobian, set this to n_ensemble
integer (mi)                   :: ld_j_x

! reported back by the MINPACK solver - should be 0 for no problems
integer (mi)                   :: info

! DrHook variables
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_HYDRATE'

!===========================================================================
! start of the executable code for vera_hydrate
!===========================================================================

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! number of variables in the system to solve, ie number of particle species
n = size(x_values, kind = mi)

! leading dimension of Jacobian, i.e. n
ld_j_x = n

! pop the input variable values into the MINPACK solver's set
x_solver = x_values

! call the MINPACK solver
call MinpackDriver ( VeraFunction, n, ld_j_x, x_solver, info, fvec, j_x )

! pop MINPACK solver's set of variables back into this subroutine's set
if ( info == minpack_parameters % info_solution_converged ) then
  x_values = x_solver
else
  x_values = vera_minpack % no_growth
end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_hydrate

end module vera_hydration_mod
