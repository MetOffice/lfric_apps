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
! This module contains the functions defining the hydration equilibrium
! conditions for Vera, utilising Kohler curves.
!

module vera_function_mod

use vera_kind_mod,    only: mr => minpack_r, mi => minpack_i

use vera_global_mod,  only: vera_constants          ,                          &
                            vera_kohler             ,                          &
                            vera_physics            ,                          &
                            vera_thread_minpack     ,                          &
                            vera_aerosol_population

use vera_minpack_mod, only: minpack_parameters

implicit none

! Description:
!   This module contains the functions defining the hydration equilibrium
!   conditions for Vera, utilising Kohler curves. These functions form the
!   basis of the hydration scheme in Vera, i.e. hydrating the dry aerosol
!   population to be in equilibrium with the environmental humidity.
!
! Method:
!   This module consists of a subroutine, defining the function
!   used to hydrate aerosol particles:
!
!     vera_kohler_function
!       Computes the equlibrium equations F(g)=0 that comprise the Vera
!       hydration scheme, and the Jacobian J(g) of F(g).
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character(len=*), parameter, private :: ModuleName='VERA_FUNCTION_MOD'

private

! make the subroutine in this module Public, as it's referenced by
! vera_hydration_mod.F90 and vera_scheme_mod.F90 and
! called by vera_minpack_mod.F90
public :: vera_kohler_function

contains

  !=============================================================================
  !
  ! vera_kohler_function - Computes the equilibrium conditions
  !                        that define the Vera hydration scheme,
  !                        F(g)=0, and the Jacobian J(g) of F(g).
  !
  !=============================================================================

subroutine vera_kohler_function( n_ensemble, ld_j_x, x_in, iflag, f_x, j_x )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_kohler_function
!
! NOTE: the type declarations for these variables match with those
!       used by the MINPACK solver in vera_minpack_mod.F90.
!       These Kinds are defined in the module vera_kind_mod as
!       minpack_i and minpack_r
!===========================================================================

! number of aerosol species
integer (mi), intent(in)       :: n_ensemble

! growth factors at which to evaluate the function F(g) = 0
real    (mr), intent(in)       :: x_in(1:n_ensemble)

! function F(g) = 0 that describes the Vera model
real    (mr), intent(out)      :: f_x(1:n_ensemble)

! Jacobian of function F
real    (mr), intent(out)      :: j_x(1:n_ensemble, 1:n_ensemble)

! leading dimension of Jacobian, set this to n_ensemble
integer (mi), intent(in)       :: ld_j_x

! integer output variable
integer (mi), intent(in out)   :: iflag

!===========================================================================
! local variables for vera_kohler_function
!===========================================================================

! loop index
integer (mi)                   :: ii

! growth factors to use to evaluate the function F(g) = 0
real    (mr)                   :: x(1:n_ensemble)

! logarithm of the equilibrium Rh of each aerosol species computed using
! the Kohler curves
real    (mr)                   :: ln_rh(1:n_ensemble)

! normalised surface tension constant
real    (mr)                   :: a_star(1:n_ensemble)

! total liquid water
real    (mr)                   :: ql

! partial derivative d q_l / d g_i
real    (mr)                   :: partial_dqcl_dg(1:n_ensemble)

! partial derivative d H_i / d g_i
real    (mr)                   :: partial_dh_dg(1:n_ensemble)

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_KOHLER_FUNCTION'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for vera_kohler_function
!===========================================================================

! form a working copy of the initial guess growth factors
x = x_in

! if the input growth factor is less than 1, then use a growth factor of 1
do ii = vera_constants%one_i, n_ensemble
  if ( x(ii) < vera_constants%one ) x(ii) = vera_constants%one
end do

! compute the normalised surface tension constant, a_star
a_star = vera_kohler%a0 / vera_aerosol_population%rd

! compute the logarithm of the equilibrium Rh using the Kohler curve,
! for each aerosol species
ln_rh = ( a_star / x )                                      -                  &
        ( vera_aerosol_population%b0                        /                  &
          ( (x**vera_constants%three) - vera_constants%one )  )

!==========================================================================
! evaluate either the function F(g), or its Jacobian J(g),
! as directed by an integer flag:
!
!   iflag == minpack_parameters % iflag_function compute Vera function F(g)
!
!   iflag == minpack_parameters % iflag_jacobian compute Jacobian J(g)
!
!==========================================================================

if ( iflag == minpack_parameters % iflag_function ) then
  !========================================================================
  ! compute the set of n equations F(g)=0 for the equation solver
  ! where n is the number of distinct aerosol speocies, n__ensemble
  !========================================================================

  ! compute the total liquid water, in kg/kg
  ql = sum( vera_constants%four_thirds_pi                        *             &
            ( vera_aerosol_population%rd**vera_constants%three ) *             &
            ( (x**vera_constants%three) - vera_constants%one   ) *             &
            vera_physics%rho_water                               *             &
            vera_aerosol_population%nc                             )

  ! the conservation of water term
  f_x(1) = vera_thread_minpack%q_total_use -                                   &
           ( exp(ln_rh(1))*vera_thread_minpack%q_sat ) - ql

  ! compute the differences in ln(Rh) between the aerosol species
  if (n_ensemble > 1) then

    f_x(2:n_ensemble) = ln_rh( 2: n_ensemble                        ) -        &
                        ln_rh( 1: n_ensemble - vera_constants%one_i )

  end if

else if ( iflag == minpack_parameters % iflag_jacobian ) then
  !========================================================================
  ! compute the Jacobian J(g) of F(g)=0 for the equation solver
  !
  ! this is a sparse n x n matrix where n is the number of aerosol
  ! species, n_ensemble, with entries across the first row
  ! and then just two entries per row i.e. entries are:
  !
  ! J(1, i)               for i = 1,n
  !
  ! J(i, i-1) and J(i, i) for i > 1
  !
  !========================================================================

  ! initialise the Jacobian so that all entries are 0
  ! Remember: this is a sparse matrix, and as the number of
  !           aerosol species increases, then the sparser the Jacobian
  j_x = vera_constants%zero

  ! compute the partial derivative d q_l / d g_i
  partial_dqcl_dg = vera_constants%four_pi                            *        &
               ( vera_aerosol_population%rd ** vera_constants%three ) *        &
                    vera_physics%rho_water                            *        &
                    vera_aerosol_population%nc                        *        &
                    ( x ** vera_constants%two )

  ! compute the partial derivative d H_i / d g_i
  partial_dh_dg = ( ( vera_constants%three                            *        &
                      vera_aerosol_population%b0                      *        &
                      (x ** vera_constants%two) )                     /        &
                    ( ( ( x ** vera_constants%three )                 -        &
                         vera_constants%one )**vera_constants%two ) ) -        &
                         ( a_star / (x **  vera_constants%two) )

  ! compute the entry J(1,1) of the Jacobian J(g) of F(g)
  ! NOTE: for monodisperse aerosol, this is the only entry in the Jacobian
  j_x(1,1) = -(exp( ln_rh(1) ) * partial_dh_dg(1) *                            &
              vera_thread_minpack%q_sat)          -                            &
              partial_dqcl_dg(1)

  ! if there are more than one aerosol species, i.e. not monodisperse,
  ! then compute the rest of the Jacobian
  if (n_ensemble > vera_constants%one_i) then

    ! loop over the aerosol species counter
    do ii = vera_constants%two_i, n_ensemble

      ! compute the rest of the 1st row entries J(1,i) with i > 1
      j_x(1, ii)                       = - partial_dqcl_dg(ii)

      ! compute the entries J(i, i-1) with i > 1
      j_x(ii, ii-vera_constants%one_i) =                                       &
                               - partial_dh_dg(ii-vera_constants%one_i)

      ! compute the entries J(i, i) with i > 1
      j_x(ii, ii)                      =   partial_dh_dg(ii)

    end do

  end if

  ! end of checking iflag, i.e. compute F(x)=0 or J(x)
end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_kohler_function

end module vera_function_mod
