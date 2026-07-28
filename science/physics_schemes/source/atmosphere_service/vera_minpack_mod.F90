! *********************************COPYRIGHT***********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *********************************COPYRIGHT***********************************
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_vera
!
! Part of the Vera scheme for diagnosing visual range.
!
! This module contains the Vera solver, based on the MINPACK routines.
!
! For details of the MINPACK solver routines, look at
!
!   "User Guide for MINPACK-1"
!    Jorge J. More, Burton S. Garbow, Kenneth E. Hillstrom, 1980.
!

module vera_minpack_mod

! grab the data types to use - real, integer and logical
use vera_kind_mod, only: mr       => minpack_r  ,                              &
                         mi       => minpack_i  ,                              &
                         mlogical => minpack_l

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

! Description:
!
!   This module contains the Vera version of the MINPACK solver.
!
!   Reference: for more detail about MINPACK, please see
!
!     Jorge More, Burton Garbow, Kenneth Hillstrom,
!     User Guide for MINPACK-1,
!     Technical Report ANL-80-74,
!     Argonne National Laboratory, 1980.
!
!   Author:
!
!     Original Fortran 77 version by Jorge More, Burton Garbow,
!     Kenneth Hillstrom, and a Fortran 90 version by John Burkardt.
!
!     This Fortran 2003 version was written by a Met Office employee.
!
!  Method:
!
!    This module defines interfaces to the user defined functions that
!    describe the nonlinear set of coupled equations, F(x) = 0, that are to
!    be solved by MINPACK. These interfaces are:
!
!      MinpackUserFunction
!         used by the solver MinpackSolver
!
!    This module also comprises the subroutines:
!
!      MinpackDriver
!        Driver for the solver MinpackSolver.
!
!        Requires subroutine MinpackSolver.
!
!      MinpackSolver
!        Solver which can use either a numerical approximation
!        to compute the Jacobian J(x) of the system of
!        defining equations F(x)=0, or an explicit expression
!        to compute the Jacobian J(x) of the system of
!        defining equations F(x)=0.
!
!        Requires subroutines MinpackEnorm
!                             MinpackQrfac
!                             MinpackQform
!                             MinpackDogleg
!                             MinpackR1updt
!                             MinpackAQ
!                             MinpackFdJacNbyN
!
!      MinpackEnorm
!        Computes the Euclidian norm of a vector.
!
!      MinpackDogleg
!        Finds the minimizing combination of Gauss-Newton and
!        gradient steps.
!
!        Requires subroutine MinpackEnorm.
!
!      MinpackQform
!        Computes the explicit QR factorization of a matrix.
!
!      MinpackQrfac
!        Computes the QR factorization of a matrix A(m,n) using
!        Householder transformations.
!
!        Requires subroutine MinpackEnorm.
!
!      MinpackAQ
!        Computes A*Q, where Q is the product of Householder
!        transformations.
!
!      MinpackR1updt
!        Re-triangularizes a matrix after a rank one update.
!
!      MinpackFdJacNbyN
!        Estimates an n by n Jacobian matrix using forward differences.
!
!  For more detail, please refer to the Vera user guide.
!
!  Code description:
!
!    Language: Fortran 2003
!    This code is written to UMDP3 standards.
!
! MINPACK is free to use, provided the following copyright notice is displayed:
!
! Minpack Copyright Notice (1999) University of Chicago. All rights reserved
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,!
! this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! 3. The end-user documentation included with the redistribution, if any, must
! include the following acknowledgment:
!
!   "This product includes software developed by the
!   University of Chicago, as Operator of Argonne National
!   Laboratory.
!
! Alternately, this acknowledgment may appear in the software itself, if and
! wherever such third-party acknowledgments normally appear.
!
! 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
! WITHOUT WARRANTY OF any kind. THE COPYRIGHT HOLDER, THE UNITED STATES,
! THE UNITED STATES DEPARTMENT OF ENERGY, and THEIR EMPLOYEES: (1) DISCLAIM
! any WARRANTIES, EXPRESS or IMPLIED, INCLUDING BUT not LIMITED TO any
! IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
! TITLE or NON-INFRINGEMENT, (2) do not ASSUME any LEGAL LIABILITY or
! RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, or USEFULNESS OF THE
! SOFTWARE, (3) do not REPRESENT THAT use OF THE SOFTWARE WOULD not INFRINGE
! PRIVATELY OWNED RIGHTS, (4) do not WARRANT THAT THE SOFTWARE WILL function
! UNINTERRUPTED, THAT IT IS ERROR-FREE or THAT any ERRORS WILL BE CORRECTED.
!
! 5. LIMITATION OF LIABILITY. in NO EVENT WILL THE COPYRIGHT HOLDER, THE
! UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, or THEIR EMPLOYEES:
! BE LIABLE FOR any INDIRECT, INCIDENTAL, CONSEQUENTIAL, SPECIAL or PUNITIVE
! DAMAGES OF any kind or NATURE, INCLUDING BUT not LIMITED TO LOSS OF
! PROFITS or LOSS OF data, FOR any REASON WHATSOEVER, WHETHER SUCH LIABILITY
! IS ASSERTED ON THE BASIS OF CONTRACT, TORT (INCLUDING NEGLIGENCE or STRICT
! LIABILITY), or OTHERWISE, EVEN if any OF SAID PARTIES HAS BEEN WARNED OF
! THE POSSIBILITY OF SUCH LOSS or DAMAGES.
!

! DrHook parameters
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit

! name of this module
character(len=*), parameter, private :: ModuleName='VERA_MINPACK_MOD'

private

public :: MinpackUserFunction ,                                                &
          MinpackDriver       ,                                                &
          minpack_parameters

! define the interface for the minpack user defined procedure
ABSTRACT interface

  ! specify the argument list for the subroutine MinpackUserFunction
  ! - this function defines both the system for minpack to solve
  !   and its analytic Jacobian, used by procedure MinpackDriver
  subroutine MinpackUserFunction( n, ld_j_x, x, iflag, fvec, j_x )

  ! precision to use - import the minpack real and integer kinds
  import :: mr, mi

  implicit none

  ! define the types and kinds of the procedure arguments
  integer ( kind = mi ), intent (in)     :: n
  integer ( kind = mi ), intent (in)     :: ld_j_x
  real    ( kind = mr ), intent (in)     :: x(n)
  integer ( kind = mi ), intent (in out) :: iflag
  real    ( kind = mr ), intent (out)    :: fvec(n)
  real    ( kind = mr ), intent (out)    :: j_x(n, n)

  end subroutine MinpackUserFunction

end interface

! define a flag to control the behaviour of the subroutine MinpackSolver,
! whether to compute the function UserFunction at the
! required values, or compute the Jacobian at the required values.
! The default is to compute the function.
integer (mi), parameter :: iflag_function        = 1_mi
integer (mi), parameter :: iflag_jacobian        = 2_mi
integer (mi), parameter :: iflag_default         = iflag_function

! default maximum number of calls to UserFunction during an iteration
! is set to maxfev_default ( n + 1 ) where n is the order of the set
! of equations to solve
integer (mi), parameter :: maxfev_default        = 100_mi

! scaling option for procedure MinpackSolver -
!    1  variables will be scaled internally
!    2  scaling is specified by the input DIAG vector
integer (mi), parameter :: scaling_mode_internal = 1_mi
integer (mi), parameter :: scaling_mode_diag     = 2_mi
integer (mi), parameter :: scaling_mode_default  = scaling_mode_diag

! Jacobian evaluation option for procedure MinpackDriver -
!    2  use the explicit Jacobian in UserFunction
!    1  compute a numerical approximation to the Jacobian
integer (mi), parameter :: jacobian_mode_explicit  = 2_mi
integer (mi), parameter :: jacobian_mode_numerical = 1_mi
integer (mi), parameter :: jacobian_mode_default   = jacobian_mode_explicit

! Euclidian norm function to use -
!   1 - use short function, i.e. norm = sqrt( sum( x**2) )
!   2 - use more complex function that utilises scaling
integer (mi), parameter :: norm_mode_short   = 1_mi
integer (mi), parameter :: norm_mode_long    = 2_mi
integer (mi), parameter :: norm_mode_default = norm_mode_short

! maximum number iterations to cycle through during
! slow convergence of the solution
integer (mi), parameter :: max_iterations_slow           = 10_mi

! maximum number iterations to cycle through during
! very slow convergence of the solution
integer (mi), parameter :: max_iterations_very_slow      = 20_mi

! info code returned from this solver - duff inputs
integer (mi), parameter :: info_duff_inputs              = 0_mi

! info code returned from this solver - solution converged
integer (mi), parameter :: info_solution_converged       = 1_mi

! info code returned from this solver - max function iterations
integer (mi), parameter :: info_max_fn_iterations        = 2_mi

! info code returned from this solver - tolerance too small
integer (mi), parameter :: info_tolerance_too_small      = 3_mi

! info code returned from this solver - slow convergence
integer (mi), parameter :: info_slow_convergence         = 4_mi

! info code returned from this solver - very slow convergence
integer (mi), parameter :: info_very_slow_convergence    = 5_mi

! parameter used to determine if the Jacobian should be recalculated
integer (mi), parameter :: max_ncfail                    = 2_mi

! determines the initial step bound. This bound is set to the product of
! factor and the euclidean norm of diag*x if nonzero, or else to
! factor itself.  In most cases, FACTOR should lie in the interval
! (0.1, 100) with 100 the recommended value.
real (mr), parameter    :: factor_default                = 100.0_mr

! if the ratio of predicted to actual reduction in solution error
! is greater then this theshold, then the current solution iteration is
! deemed successful
real (mr), parameter    :: ratio_good_default            = 0.0001_mr

! threshold for determining slow convergence of the solution
real (mr), parameter    :: slow_convergence_default      = 0.1_mr

! threshold for determining very slow convergence of the solution
real (mr), parameter    :: very_slow_convergence_default = 0.0000001_mr

! default tolerance to use when computing solutions
real (mr), parameter    :: tolerance_default             = 0.000001_mr

! parameter for determining if the required tolerance is too
! small to allow for successful convergence of the solution
real (mr), parameter    :: tolerance_test_default        = 0.1_mr

! used in determining a suitable step length for the forward-difference
! approximation. This approximation assumes that the relative errors in the
! functions are of the order of epsilon_min. If epsilon_min is less than
! the machine precision, it is assumed that the relative errors in the
! functions are of the order of the machine precision.
real (mr), parameter    :: epsilon_min_default           = 0.0_mr

! gather the MINPACK solver parameters into a handy structure
type :: MinpackParameters
  integer (mi) :: iflag_function             = iflag_function
  integer (mi) :: iflag_jacobian             = iflag_jacobian
  integer (mi) :: iflag_default              = iflag_default
  integer (mi) :: maxfev_use                 = maxfev_default
  integer (mi) :: scaling_mode_internal      = scaling_mode_internal
  integer (mi) :: scaling_mode_diag          = scaling_mode_diag
  integer (mi) :: scaling_mode_use           = scaling_mode_default
  integer (mi) :: jacobian_mode_explicit     = jacobian_mode_explicit
  integer (mi) :: jacobian_mode_numerical    = jacobian_mode_numerical
  integer (mi) :: jacobian_mode              = jacobian_mode_default
  integer (mi) :: norm_mode_short            = norm_mode_short
  integer (mi) :: norm_mode_long             = norm_mode_long
  integer (mi) :: norm_mode                  = norm_mode_default
  integer (mi) :: max_iterations_slow        = max_iterations_slow
  integer (mi) :: max_iterations_very_slow   = max_iterations_very_slow
  integer (mi) :: info_duff_inputs           = info_duff_inputs
  integer (mi) :: info_solution_converged    = info_solution_converged
  integer (mi) :: info_max_fn_iterations     = info_max_fn_iterations
  integer (mi) :: info_tolerance_too_small   = info_tolerance_too_small
  integer (mi) :: info_slow_convergence      = info_slow_convergence
  integer (mi) :: info_very_slow_convergence = info_very_slow_convergence
  integer (mi) :: max_ncfail                 = max_ncfail
  real    (mr) :: factor_use                 = factor_default
  real    (mr) :: ratio_good                 = ratio_good_default
  real    (mr) :: slow_convergence           = slow_convergence_default
  real    (mr) :: very_slow_convergence      = very_slow_convergence_default
  real    (mr) :: tolerance                  = tolerance_default
  real    (mr) :: tolerance_test             = tolerance_test_default
  real    (mr) :: epsilon_min                = epsilon_min_default
end type MinpackParameters

type(MinpackParameters), save :: minpack_parameters

! value of zero as a real number
real    (mr), parameter :: minpack_zero_real = 0.0_mr

! value of 0.05 as a real number
real    (mr), parameter :: minpack_p05_real  = 0.05_mr

! value of 0.25 as a real number
real    (mr), parameter :: minpack_p25_real  = 0.25_mr

! value of 0.1 as a real number
real    (mr), parameter :: minpack_p1_real   = 0.1_mr

! value of 0.5 as a real number
real    (mr), parameter :: minpack_p5_real   = 0.5_mr

! value of one as a real number
real    (mr), parameter :: minpack_one_real  = 1.0_mr

! value of zero as an integer
integer (mi), parameter :: minpack_zero_int  = 0_mi

! value of one as an integer
integer (mi), parameter :: minpack_one_int   = 1_mi

! value of two as an integer
integer (mi), parameter :: minpack_two_int   = 2_mi

! largest possible real number
real    (mr), parameter :: giant = huge( giant )

! square root of smallest possible real number
real    (mr), parameter :: square_root_tiny  = sqrt(tiny( square_root_tiny  ))

! square root largest possible real number
real    (mr), parameter :: square_root_giant = sqrt(huge( square_root_giant ))

! least positive number, e, such that 1 + e > e
real    (mr), parameter :: epsmch = epsilon( epsmch )

contains

subroutine MinpackDriver ( UserFunction, n, ldfjac, x , info, fvec, fjac )

  !===========================================================================
  !
  ! MinpackDriver is a driver for the solver MinpackSolver.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackDriver
!===========================================================================

! define the procedure UserFunction using an abstract interface
procedure(MinpackUserFunction):: UserFunction

! number of variables in the system to solve
integer (mi), intent(in)      :: n

! leading dimension of fjac, must be at least n
integer (mi), intent(in)      :: ldfjac

! variable values to use for the MINPACK solver
real    (mr), intent(in out)  :: x(n)

! error code reported back by the MINPACK solver.
! If the user has terminated execution, info is set to the (negative)
! value of iflag - see description of UserFunction. Otherwise, info
! is set as follows:
!    0  improper input parameters.
!    1  algorithm estimates that the relative error between X and the
!       solution is at most tolerance.
!    2  number of calls to UserFunction with flag = 1 has reached
!       minpack_parameters % maxfev_use.
!    3  tolerance is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, iteration is not making good progress, as measured by the
!       improvement from the last five Jacobian evaluations.
!    5, iteration is not making good progress, as measured by the
!       improvement from the last ten iterations.
integer (mi), intent(out)     :: info

! used by the MINPACK solver to report back the final value
! of UserFunction(x)
real    (mr), intent(out)     :: fvec(n)

! orthogonal matrix produced by the QR factorization of the approximate
! Jacobian
real    (mr), intent(out)     :: fjac(ldfjac, n)

!===========================================================================
! local variables for MinpackDriver
!===========================================================================

! termination occurs when the number of calls to UserFunction
! is at least maxfev by the end of an iteration
integer (mi)                  :: maxfev

! tolerance to pass to solver MinpackSolver
real    (mr)                  :: xtol

! scaling option -
!    1  variables will be scaled internally
!    2  scaling is specified by the input DIAG vector
integer (mi)                  :: scaling_mode

! size of the R array, which must be no less than (n*(n+1))/2
integer (mi)                  :: lr

! the number of calls to UserFunction with iflag = 1, i.e.
! to compute the function
integer (mi)                  :: nfev

! the number of calls to UserFunction with iflag = 2, i.e.
! to compute the Jacobian
integer (mi)                  :: njev

! dummy argument for calling the solver MinpackSolver
integer (mi)                  :: ml

! dummy argument for calling the solver MinpackSolver
integer (mi)                  :: mu

! If scaling_mode = 1, then diag is set internally.
! If scaling_mode = 2, then diag must contain positive entries that
! serve as multiplicative scale factors for the variables - these
! entries are set to 1.
real    (mr)                  :: diag(n)

! determines the initial step bound. This bound is set to the product of
! factor and the euclidean norm of diag*x if nonzero, or else to
! factor itself.  In most cases, factor should lie in the interval
! (0.1, 100) with 100 the recommended value.
real    (mr)                  :: factor

! upper triangular matrix produced by the QR factorization of the final
! approximate Jacobian, stored row-wise
real    (mr)                  :: r( ( n * (n + minpack_one_int) ) /            &
                                    minpack_two_int )

! vector Q'*FVEC
real    (mr)                  :: qtf(n)

! dummy argument for calling the solver MinpackSolver
real    (mr)                  :: epsilon_min

! DrHook variables
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'MINPACKDRIVER'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackDriver
!===========================================================================

! initialise the return code info to zero
info = minpack_zero_int

! check the inputs to this procedure
info = minpack_zero_int
if ( ( n         <= minpack_zero_int ) .or.                                    &
     ( ldfjac    < n                 )      ) then
  if (lhook) call dr_hook(                                                     &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  return
end if

! set parameters to pass to the solver MinpackSolver
maxfev        = minpack_parameters % maxfev_use * ( n + minpack_one_int )
xtol          = minpack_parameters % tolerance
ml            = n - minpack_one_int
mu            = n - minpack_one_int
epsilon_min   = minpack_parameters % epsilon_min
scaling_mode  = minpack_parameters % scaling_mode_use
factor        = minpack_parameters % factor_use
diag(:)       = minpack_one_real

! set the size of R, the upper triangular matrix produced
! by the QR factorization of the final approximate Jacobian
lr = ( n * ( n + minpack_one_int ) ) / minpack_two_int

! initialise arrays returned by the solver MinpackSolver
nfev          = minpack_zero_int
njev          = minpack_zero_int
fjac(:,:)     = minpack_zero_real
r(:)          = minpack_zero_real
qtf(:)        = minpack_zero_real

! run the solver MinpackSolver
call MinpackSolver ( UserFunction, n, ldfjac, ml, mu, maxfev,                  &
                     scaling_mode, lr, xtol, factor, epsilon_min,              &
                     x, diag, info, nfev, njev, fvec, fjac, r, qtf )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine MinpackDriver

subroutine MinpackSolver ( UserFunction, n, ldfjac, ml, mu, maxfev,            &
                           scaling_mode, lr, xtol, factor, epsilon_min, x,     &
                           diag,info, nfev, njev, fvec, fjac, r, qtf       )

  !===========================================================================
  !
  ! MinpackSolver seeks a zero of N nonlinear equations in N variables.
  !
  ! Discussion:
  !
  !    MinpackSolver finds a zero of a system of N nonlinear functions
  !    in N variables by a modification of the Powell hybrid method.
  !    The user must provide a subroutine, UserFunction ,which calculates
  !    the functions. The Jacobian is can either be approximated numerically
  !    from the function set, or computed using an explicit expression.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackSolver
!===========================================================================

! define the procedure UserFunction using an abstract interface
procedure(MinpackUserFunction):: UserFunction

! number of variables in the system to solve
integer (mi), intent(in)     :: n

! leading dimension of fjac, the approximate Jacobian
integer (mi), intent(in)     :: ldfjac

! number of subdiagonals within the band of the Jacobian matrix. If the
! matrix is not banded, then set ml to at least n - 1.
integer (mi)                 :: ml

! number of superdiagonals within the band of the Jacobian matrix. If the
! matrix is not banded, then set mu to at least n - 1.
integer (mi)                 :: mu

! termination occurs when the number of calls to UserFunction
! is at least maxfev by the end of an iteration
integer (mi), intent(in)     :: maxfev

! scaling option -
!    1  variables will be scaled internally
!    2  scaling is specified by the input DIAG vector
integer (mi), intent(in)     :: scaling_mode

! size of the R array, which must be no less than (n*(n+1))/2
integer (mi), intent(in)     :: lr

! numerical tolerance of the solution computed by the MINPACK solver
real    (mr), intent(in)     :: xtol

! determines the initial step bound. This bound is set to the product of
! factor and the euclidean norm of diag*x if nonzero, or else to
! factor itself.  In most cases, factor should lie in the interval
! (0.1, 100) with 100 the recommended value.
real    (mr), intent(in)     :: factor

! used in determining a suitable step length for the forward-difference
! approximation. This approximation assumes that the relative errors in the
! functions are of the order of epsilon_min. If epsilon_min is less than
! the machine precision, it is assumed that the relative errors in the
! functions are of the order of the machine precision.
real    (mr), intent(in)     :: epsilon_min

! on input, the initial guess to the solution, on output the solution
! to the set of N nonlinear functions
real    (mr), intent(in out) :: x(n)

! If scaling_mode = 1, then diag is set internally.
! If scaling_mode = 2, then diag must contain positive entries that
! serve as multiplicative scale factors for the variables - these
! entries are set to 1.
real    (mr), intent(in out) :: diag(n)

! error code reported back by the MINPACK solver.
! If the user has terminated execution, info is set to the (negative)
! value of iflag - see description of UserFunction. Otherwise, info
! is set as follows:
!    0  improper input parameters.
!    1  algorithm estimates that the relative error between X and the
!       solution is at most tolerance.
!    2  number of calls to UserFunction with flag = 1 has reached
!       minpack_parameters % maxfev_use.
!    3  tolerance is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, iteration is not making good progress, as measured by the
!       improvement from the last five Jacobian evaluations.
!    5, iteration is not making good progress, as measured by the
!       improvement from the last ten iterations.
integer (mi), intent(out)    :: info

! number of calls to UserFunction with iflag = 1, i.e.
! to evaluate the function
integer (mi), intent(out)    :: nfev

! number of evaluations of the Jacobian
integer (mi), intent(out)    :: njev

! used by the MINPACK solver to report back the final value
! of UserFunction(x)
real    (mr), intent(out)    :: fvec(n)

! orthogonal matrix produced by the QR factorization of the approximate
! Jacobian
real    (mr), intent(out)    :: fjac(ldfjac, n)

! upper triangular matrix produced by the QR factorization of the final
! approximate Jacobian, stored rowwise
real    (mr), intent(out)    :: r(lr)

! vector Q'*FVEC
real    (mr), intent(out)    :: qtf(n)

!===========================================================================
! local variables for MinpackSolver
!===========================================================================

! loop counters
integer (mi)                 :: i
integer (mi)                 :: j

! array index
integer (mi)                 :: l

! flag that controls calls to UserFunction
integer (mi)                 :: iflag

! iteration counter
integer (mi)                 :: iter

! dummy argument used to call MinpackQrfac - value is ignored
integer (mi)                 :: iwa(minpack_one_int)

! dummy argument used to call MinpackQrfac
integer (mi)                 :: l_one = minpack_one_int

! counts failed iterations
integer (mi)                 :: ncfail

! counts slow iterations
integer (mi)                 :: nslow1

! counts very slow iterations, involving recomputing the Jacobian
integer (mi)                 :: nslow2

! counts sucessfull step increments
integer (mi)                 :: ncsuc

! controls whether pivoting is used to factorise the Jacobian matrix,
! set to false so no pivoting
logical (mlogical),parameter :: pivot = .false.

! indicates if the Jacobian has been recomputed
logical (mlogical)           :: jeval

! dummy argument used to call MinpackR1updt - value is ignored
logical (mlogical)           :: sing

! scaled actual error reduction in the iterated solution
real    (mr)                 :: actred

! step bound
real    (mr)                 :: delta

! Euclidian norm of the scaled x
real    (mr)                 :: xnorm

! Euclidian norm of a trial solution at x
real    (mr)                 :: fnorm

! Euclidian norm of a trial solution at x + p
real    (mr)                 :: fnorm1

! Euclidian norm of p
real    (mr)                 :: pnorm

! scaled, prediction reduction in the solution error
real    (mr)                 :: prered

! ratio of actual to predicted solution error
real    (mr)                 :: ratio

! partial sum used to compute the scaled, predicted reduction in
! the solution error
real    (mr)                 :: sum2

! temporary store used to compute the scaled, predicted reduction in
! the solution error
real    (mr)                 :: temp

! working arrays
real    (mr)                 :: wa1(n)
real    (mr)                 :: wa2(n)
real    (mr)                 :: wa3(n)
real    (mr)                 :: wa4(n)

! DrHook variables
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'MINPACKSOLVER'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackSolver
!===========================================================================

! initialise flags and counters to zero
info  = minpack_parameters % info_duff_inputs
iflag = minpack_parameters % iflag_default
nfev  = minpack_zero_int
njev  = minpack_zero_int

! initialise the outputs from this procedure
fvec(:)   = minpack_zero_real
fjac(:,:) = minpack_zero_real
r(:)      = minpack_zero_real
qtf(:)    = minpack_zero_real

! check the input parameters for errors
if ( n      <= minpack_zero_int  .or.                                          &
     ldfjac <  n                 .or.                                          &
     ml     <  minpack_zero_int  .or.                                          &
     mu     <  minpack_zero_int  .or.                                          &
     xtol   <  minpack_zero_real .or.                                          &
     maxfev <= minpack_zero_int  .or.                                          &
     factor <= minpack_zero_real .or.                                          &
     lr     < ( n * ( n + minpack_one_int ) ) / minpack_two_int ) then

  ! at least one input is duff, so return from this solver
  info = minpack_parameters % info_duff_inputs
  if (lhook) call dr_hook(                                                     &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  return

end if

! check for negative input scaling factors
if ( scaling_mode == minpack_parameters % scaling_mode_diag ) then
  do j = minpack_one_int, n
    if ( diag(j) <= minpack_zero_real ) then

      ! found a negative scaling factor, so return from this solver
      info = minpack_parameters % info_duff_inputs
      if (lhook) call dr_hook(                                                 &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      return

    end if
  end do
end if

! evaluate the function at the starting point and calculate its norm
iflag = minpack_parameters % iflag_function
call UserFunction ( n, ldfjac, x, iflag, fvec, fjac )
nfev = minpack_one_int

if ( iflag < minpack_zero_int ) then
  ! problem calling UserFunction, so return
  info  = iflag
  if (lhook) call dr_hook(                                                     &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  return
end if

fnorm = MinpackEnorm ( n, fvec )

! initialize iteration counter and monitors
iter   = minpack_one_int
ncsuc  = minpack_zero_int
ncfail = minpack_zero_int
nslow1 = minpack_zero_int
nslow2 = minpack_zero_int

! beginning of the outer loop
do

  jeval = .true.
  ! calculate the Jacobian matrix
  iflag = minpack_parameters % iflag_jacobian

  if ( minpack_parameters % jacobian_mode ==                                   &
       minpack_parameters % jacobian_mode_explicit ) then

    ! use the explicit Jacobian defined in UserFunction
    call UserFunction ( n, ldfjac, x, iflag, fvec, fjac )

  else

    ! use a numerical approximation to the Jacobian
    call MinpackFdJacNbyN ( UserFunction, n, ldfjac, ml, mu,                   &
                            x, fvec, epsilon_min, iflag, fjac )

  end if

  ! Jacobian has been evaluated, so increment the counter
  njev = njev + minpack_one_int

  if ( iflag < minpack_zero_int ) then
    ! problem calling UserFunction, so return
    info = iflag
    if (lhook) call dr_hook(                                                   &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    return
  end if

  ! compute the QR factorization of the Jacobian
  call MinpackQrfac ( n, n, ldfjac, l_one, pivot, fjac, iwa, wa1, wa2 )

  ! on the first iteration, if scaling_mode is internal, scale according
  ! to the norms of the columns of the initial Jacobian
  if ( iter == minpack_one_int ) then

    if ( scaling_mode == minpack_parameters % scaling_mode_internal ) then
      diag(minpack_one_int:n) = wa2(minpack_one_int:n)
      do j = minpack_one_int, n
        if ( wa2(j) == minpack_zero_real ) then
          diag(j) = minpack_one_real
        end if
      end do
    end if

    ! on the first iteration, calculate the norm of the scaled x
    ! and initialize the step bound delta
    wa3(minpack_one_int:n) = diag(minpack_one_int:n) * x(minpack_one_int:n)
    xnorm = MinpackEnorm ( n, wa3 )
    delta = factor * xnorm
    if ( delta == minpack_zero_real ) then
      delta = factor
    end if

  end if

  ! form Q'*FVEC and store in qtf
  qtf(minpack_one_int:n) = fvec(minpack_one_int:n)
  do j = minpack_one_int, n
    if ( fjac(j,j) /= minpack_zero_real ) then
      temp     = - dot_product( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
      qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
    end if
  end do

  ! copy the triangular factor of the QR factorization into r
  sing = .false.
  do j = minpack_one_int, n
    l = j
    do i = minpack_one_int, j - minpack_one_int
      r(l) = fjac(i,j)
      l    = l + n - i
    end do
    r(l) = wa1(j)
    if ( wa1(j) == minpack_zero_real ) then
      sing = .true.
    end if
  end do

  ! accumulate the orthogonal factor in fjac
  call MinpackQform ( n, n, ldfjac, fjac )

  ! rescale if necessary
  if ( scaling_mode /= minpack_parameters % scaling_mode_diag ) then
    do j = minpack_one_int, n
      diag(j) = max( diag(j), wa2(j) )
    end do
  end if

  ! beginning of the inner loop
  do

    ! determine the direction p
    call MinpackDogleg ( n, r, lr, diag, qtf, delta, wa1 )

    ! store the direction p and x + p
    wa1(minpack_one_int:n) = - wa1(minpack_one_int:n)
    wa2(minpack_one_int:n) =     x(minpack_one_int:n) +                        &
                               wa1(minpack_one_int:n)
    wa3(minpack_one_int:n) =  diag(minpack_one_int:n) *                        &
                               wa1(minpack_one_int:n)

    ! calculate the norm of p
    pnorm = MinpackEnorm ( n, wa3 )

    ! on the first iteration, adjust the initial step bound
    if ( iter == minpack_one_int ) then
      delta = min( delta, pnorm )
    end if

    ! evaluate the function at x + p and calculate its norm
    iflag = minpack_parameters % iflag_function
    call UserFunction ( n, ldfjac, wa2, iflag, wa4, fjac )
    nfev = nfev + minpack_one_int

    if ( iflag < minpack_zero_int ) then
      ! problem calling UserFunction, so return
      info = iflag
      if (lhook) call dr_hook(                                                 &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      return
    end if

    ! compute the norm at x + p
    fnorm1 = MinpackEnorm ( n, wa4 )

    ! compute the scaled actual reduction
    actred = - minpack_one_real
    if ( fnorm1 < fnorm ) then
      actred = minpack_one_real - ( fnorm1 / fnorm ) ** minpack_two_int
    end if

    ! compute the scaled predicted reduction
    l = minpack_one_int
    do i = minpack_one_int, n
      sum2 = minpack_zero_real
      do j = i, n
        sum2 = sum2 + r(l) * wa1(j)
        l    = l + minpack_one_int
      end do
      wa3(i) = qtf(i) + sum2
    end do

    temp = MinpackEnorm ( n, wa3 )
    prered = minpack_zero_real
    if ( temp < fnorm ) then
      prered = minpack_one_real - ( temp / fnorm ) ** minpack_two_int
    end if

    ! compute the ratio of the actual to the predicted reduction
    if ( prered > minpack_zero_real ) then
      ratio = actred / prered
    else
      ratio = minpack_zero_real
    end if

    ! update the step bound
    if ( ratio < minpack_p1_real ) then

      ncsuc  = minpack_zero_int
      ncfail = ncfail + minpack_one_int
      delta  = minpack_p5_real * delta

    else

      ncfail = minpack_zero_int
      ncsuc  = ncsuc + minpack_one_int

      if ( (ratio >= minpack_p5_real) .or. (ncsuc > minpack_one_int) ) then
        delta = max( delta, pnorm / minpack_p5_real )
      end if

      if ( abs( ratio - minpack_one_real ) <= minpack_p1_real ) then
        delta = pnorm / minpack_p5_real
      end if

    end if

    ! test for successful iteration
    if ( ratio >= minpack_parameters % ratio_good) then

      ! all good, so update x, wa2 and fvec
      x(   minpack_one_int:n) =  wa2(minpack_one_int:n)
      wa2( minpack_one_int:n) = diag(minpack_one_int:n) * x(minpack_one_int:n)
      fvec(minpack_one_int:n) =  wa4(minpack_one_int:n)

      ! update norms
      xnorm = MinpackEnorm ( n, wa2 )
      fnorm = fnorm1

      ! update the iteration counter
      iter  = iter + minpack_one_int

    end if

    ! determine the progress of the iteration
    nslow1 = nslow1 + minpack_one_int
    if ( actred >= minpack_parameters % very_slow_convergence ) then
      nslow1 = minpack_zero_int
    end if

    if ( jeval ) then
      nslow2 = nslow2 + minpack_one_int
    end if

    if ( actred >= minpack_parameters % slow_convergence ) then
      nslow2 = minpack_zero_int
    end if

    ! test for convergence
    if ( (delta <= xtol * xnorm) .or. (fnorm == minpack_zero_real) ) then
      ! solution has converged
      info = minpack_parameters % info_solution_converged
    end if

    if ( info == minpack_parameters % info_solution_converged ) then
      ! solution has converged, so return from this solver
      ! having successfully found a solution
      if (lhook) call dr_hook(                                                 &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      return
    end if

    ! tests for termination and stringent tolerances
    if ( nfev >= maxfev ) then
      ! have used the maximum number of calls to UserFunction
      info = minpack_parameters % info_max_fn_iterations
    end if

    if (      minpack_parameters % tolerance_test                 *            &
         max( minpack_parameters % tolerance_test * delta, pnorm )             &
         <= epsmch * xnorm ) then
      ! tolerance is too small to allow successful convergence
      info = minpack_parameters % info_tolerance_too_small
    end if

    if ( nslow2 == minpack_parameters % max_iterations_slow ) then
      ! convergence is slow
      info = minpack_parameters % info_slow_convergence
    end if

    if ( nslow1 == minpack_parameters % max_iterations_very_slow ) then
      ! convergence is very slow
      info = minpack_parameters % info_very_slow_convergence
    end if

    if ( info /= minpack_zero_int ) then
      ! return from this solver - solution has not converged, but
      ! one of the above termination tests has been triggered
      if (lhook) call dr_hook(                                                 &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      return
    end if

    ! criterion for recalculating the Jacobian is a maximum number of
    ! failures to adequately reduce the solution error despite halving
    ! the step bound with each try
    if ( ncfail == minpack_parameters % max_ncfail ) then
      ! exit this inner loop, back into the outer loop where the
      ! Jacobian recalculation is done
      exit
    end if

    ! calculate the rank one modification to the Jacobian
    ! and update qtf if necessary
    do j = minpack_one_int, n
      sum2 = dot_product( wa4(minpack_one_int:n),fjac(minpack_one_int:n,j) )
      wa2(j) = ( sum2 - wa3(j) ) / pnorm
      wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
      if ( ratio >= minpack_parameters % ratio_good) then
        qtf(j) = sum2
      end if
    end do

    ! compute the QR factorization of the updated Jacobian
    call MinpackR1updt ( n, n, lr, wa1, r, wa2, wa3, sing )
    call MinpackAQ     ( n    , n, ldfjac, wa2, wa3, fjac )
    call MinpackAQ     ( l_one, n,  l_one, wa2, wa3, qtf  )

    ! end of the inner loop
    jeval = .false.
  end do

  ! end of the outer loop
end do

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine MinpackSolver

function MinpackEnorm ( n, x ) result ( enorm_computed )

  !===========================================================================
  !
  !  MinpackEnorm computes the Euclidean norm of a vector.
  !
  !  Discussion:
  !
  !    The Euclidean norm is computed by accumulating the sum of
  !    squares in three different sums.  The sums of squares for the
  !    small and large components are scaled so that no overflows
  !    occur.  Non-destructive underflows are permitted.  Underflows
  !    and overflows do not occur in the computation of the unscaled
  !    sum of squares for the intermediate components.
  !
  !    The definitions of small, intermediate and large components
  !    depend on two constants, square_root_tiny and square_root_giant. So,
  !    for vector component x(i):
  !
  !                       | x(i) | > square_root_giant      - large
  !    square_root_tiny < | x(i) | < square_root_giant      - intermediate
  !    square_root_tiny > | x(i) |                          - small
  !
  !    The main restrictions on these constants are that square_root_tiny^2
  !    does not underflow and square_root_giant^2 does not overflow.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackEnorm
!===========================================================================

! number of components in the input vector, x
integer (mi), intent(in)     :: n

! vector for which Euclidean norm is to be computed
real    (mr), intent(in)     :: x(n)

!===========================================================================
! local variables for MinpackEnorm
!===========================================================================

! computed Euclidean norm of vector x - this is the computed result of
! this function
real    (mr)                 :: enorm_computed

! normalised largest possible real number, i.e. square_root_giant / n
real    (mr)                 :: ngiant

! sum of the squares of the large components in vector x
real    (mr)                 :: sum_large

! sum of the squares of the intermediate components in vector x
real    (mr)                 :: sum_inter

! sum of the squares of the small components in vector x
real    (mr)                 :: sum_small

! absolute value of a component in the vector x
real    (mr)                 :: xabs

! largest of the "large components" of vector x
real    (mr)                 :: large_max

! largest of the "small components" of vector x
real    (mr)                 :: small_max

! counter to loop over the components of the vector x
integer (mi)                 :: ii

! DrHook variables
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'MINPACKENORM'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackEnorm
!===========================================================================

if (minpack_parameters%norm_mode /= minpack_parameters%norm_mode_long) then

   ! use the short form of the Euclidian norm function
  enorm_computed = sqrt ( sum ( x(minpack_one_int:n) ** minpack_two_int ) )

else

  ! use the long form of the Euclidian norm function

  ! initialise the sums of the squares of the components of vector x
  sum_large = minpack_zero_real
  sum_inter = minpack_zero_real
  sum_small = minpack_zero_real

  ! initialise the largest "large" and "small" components of vector x
  large_max = minpack_zero_real
  small_max = minpack_zero_real

  ! define the "large component" value
  ngiant = square_root_giant / n

  ! loop over the components of the input vector, x, checking for
  ! large and small components, and compute the scaled sums of squares
  ! sum_large, sum_inter and sum_small
  do ii = minpack_one_int, n

    xabs = abs( x(ii) )

    if ( (xabs > square_root_tiny) .and. (xabs < ngiant) ) then

      ! update sum of the squares of intermediate components
      sum_inter = sum_inter + xabs ** minpack_two_int

    else if ( xabs <= square_root_tiny ) then

      ! update sum of the squares of small components
      if ( xabs > small_max) then
        sum_small = minpack_one_real +                                         &
                    sum_small * ( small_max / xabs ) ** minpack_two_int
        small_max = xabs
      else if ( xabs /= minpack_zero_real ) then
        sum_small = sum_small + ( xabs / small_max ) ** minpack_two_int
      end if

    else

      ! update sum of the squares of large components
      if ( xabs > large_max ) then
        sum_large = minpack_one_real +                                         &
                    sum_large * ( large_max / xabs ) ** minpack_two_int
        large_max = xabs
      else
        sum_large = sum_large + ( xabs / large_max ) ** minpack_two_int
      end if

    end if

    ! end of looping over the components of the input vector, x, forming the
    ! scaled sums
  end do

  ! calculation of the Euclidian norm
  if ( sum_large /= minpack_zero_real ) then

    ! at least one large component in vector x, so ignore small components
    enorm_computed = large_max *                                               &
                     sqrt( sum_large + ( sum_inter / large_max ) / large_max )

  else if ( sum_inter /= minpack_zero_real ) then

    ! there are intermediate components in vector x
    if ( small_max <= sum_inter ) then
      enorm_computed = sqrt(sum_inter *                                        &
                         ( minpack_one_real +                                  &
                           (small_max / sum_inter) * (small_max * sum_small) ) )
    else
      enorm_computed = sqrt( small_max *                                       &
                 ( ( sum_inter / small_max ) + ( small_max * sum_small ) ) )
    end if

  else

    ! all components in vector x are small
    enorm_computed = small_max * sqrt( sum_small )

  end if

end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end function MinpackEnorm

subroutine MinpackDogleg ( n, r, lr, diag, qtb, delta, x )

  !===========================================================================
  ! MinpackDogleg finds the minimizing combination of Gauss-Newton and
  ! gradient steps.
  !
  !  Discussion:
  !
  !    Given an M by N matrix A, an N by N nonsingular diagonal
  !    matrix D, an M-vector B, and a positive number DELTA, the
  !    problem is to determine the convex combination X of the
  !    Gauss-Newton and scaled gradient directions that minimizes
  !    (A*X - B) in the least squares sense, subject to the
  !    restriction that the euclidean norm of D*X be at most DELTA.
  !
  !    This procedure completes the solution of the problem
  !    if it is provided with the necessary information from the
  !    QR factorization of A.  That is, if A = Q*R, where Q has
  !    orthogonal columns and R is an upper triangular matrix,
  !    then MinpackDogleg expects the full upper triangle of R and
  !    the first N components of Q'*B.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackDogleg
!===========================================================================

! order of the input matrix r
integer (mi), intent(in)      :: n

! size of the matrix r, must be at least (n*(n+1))/2
integer (mi), intent(in)      :: lr

! matrix which must contain the upper triangular matrix r stored by rows
real    (mr), intent(in)      :: r(lr)

! diagonal elements of the matrix d
real    (mr), intent(in)      :: diag(n)

! first n elements of the vector transpose(q) * b
real    (mr), intent(in)      :: qtb(n)

! upper bound on the Euclidian norm of d * x
real    (mr), intent(in)      :: delta

! convex combination of the Gauss-Newton direction and the scaled gradient
! direction
real    (mr), intent(out)     :: x(n)

!===========================================================================
! local variables for MinpackDogleg
!===========================================================================

! index to the input array r
integer (mi)                  :: index_r

! loop counters
integer (mi)                  :: i
integer (mi)                  :: j
integer (mi)                  :: k
integer (mi)                  :: l

! ratio of delta to the scaled gradient
real    (mr)                  :: alpha

! Euclidian norms
real    (mr)                  :: bnorm
real    (mr)                  :: gnorm
real    (mr)                  :: qnorm
real    (mr)                  :: sgnorm

! temporary sums
real    (mr)                  :: sum2
real    (mr)                  :: temp

! working arrays
real    (mr)                  :: wa1(n)
real    (mr)                  :: wa2(n)

! DrHook variables
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'MINPACKDOGLEG'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackDogleg
!===========================================================================

! calculate the Gauss-Newton direction
index_r = ( ( n * ( n + minpack_one_int ) ) / minpack_two_int ) +              &
          minpack_one_int

do k = minpack_one_int, n

  j       = n - k + minpack_one_int
  index_r = index_r - k
  l       = index_r + minpack_one_int
  sum2    = minpack_zero_real

  do i = j + minpack_one_int, n
    sum2 = sum2 + r(l) * x(i)
    l    = l + minpack_one_int
  end do

  temp = r( index_r )

  if ( temp == minpack_zero_real ) then

    l = j
    do i = minpack_one_int, j
      temp = max ( temp, abs ( r(l)) )
      l    = l + n - i
    end do

    if ( temp == minpack_zero_real ) then
      temp = epsmch
    else
      temp = epsmch * temp
    end if

  end if

  x(j) = ( qtb(j) - sum2 ) / temp

end do

! test whether the Gauss-Newton direction is acceptable
wa1(minpack_one_int:n) = minpack_zero_real
wa2(minpack_one_int:n) = diag(minpack_one_int:n) * x(minpack_one_int:n)
qnorm                  = MinpackEnorm ( n, wa2 )

if ( qnorm <= delta ) then
  if (lhook) call dr_hook(                                                     &
                  ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  return
end if

! the Gauss-Newton direction is not acceptable, so
! calculate the scaled gradient direction
l = minpack_one_int
do j = minpack_one_int, n
  temp = qtb(j)
  do i = j, n
    wa1(i) = wa1(i) + r(l) * temp
    l      = l + minpack_one_int
  end do
  wa1(j) = wa1(j) / diag(j)
end do

! calculate the norm of the scaled gradient and
! test for the special case in which the scaled gradient is zero
gnorm  = MinpackEnorm ( n, wa1 )
sgnorm = minpack_zero_real
alpha  = delta / qnorm

if ( gnorm /= minpack_zero_real ) then

  ! calculate the point along the scaled gradient which minimizes
  ! the quadratic
  wa1(minpack_one_int:n) = ( wa1(minpack_one_int:n) / gnorm ) /                &
                           diag(minpack_one_int:n)

  l = minpack_one_int
  do j = minpack_one_int, n
    sum2 = minpack_zero_real
    do i = j, n
      sum2 = sum2 + r(l) * wa1(i)
      l    = l + minpack_one_int
    end do
    wa2(j) = sum2
  end do

  temp   = MinpackEnorm ( n, wa2 )
  sgnorm = ( gnorm / temp ) / temp


  ! test whether the scaled gradient direction is acceptable
  alpha = minpack_zero_real
  if ( sgnorm < delta ) then

    ! the scaled gradient direction is not acceptable, so
    ! calculate the point along the dogleg at which the quadratic
    ! is minimized
    bnorm = MinpackEnorm ( n, qtb )
    temp = (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta)
    temp = temp - (delta / qnorm) * ( (sgnorm / delta)                  **     &
                                       minpack_two_int )                 +     &
           sqrt( ( temp             - ( ( delta  / qnorm )              **     &
                   minpack_two_int ) )                                   +     &
                 ( minpack_one_real - ( ( delta  / qnorm )              **     &
                   minpack_two_int ) )                                   *     &
                 ( minpack_one_real - ( ( sgnorm / delta )              **     &
                   minpack_two_int ) ) )

    alpha = ( ( delta / qnorm )                                          *     &
              (minpack_one_real - (sgnorm / delta) ** minpack_two_int) ) /     &
            temp
  end if

end if

! form appropriate convex combination of the Gauss-Newton
! direction and the scaled gradient direction
temp   = ( minpack_one_real - alpha ) * min( sgnorm, delta )
x(1:n) = (temp * wa1(minpack_one_int:n)) + (alpha * x(minpack_one_int:n))

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine MinpackDogleg

subroutine MinpackQform ( m, n, ldq, q )

  !===========================================================================
  !
  ! MinpackQform computes the explicit QR factorization of a matrix.
  !
  ! Discussion:
  !
  !    The QR factorization of a matrix is usually accumulated in implicit
  !    form, that is, as a series of orthogonal transformations of the
  !    original matrix.  This routine carries out those transformations,
  !    to explicitly exhibit the factorization, Q,  constructed on an input
  !    matrix A by MinpackQrfac.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackQform
!===========================================================================

! number of rows of A and the order of Q
integer (mi), intent(in)      :: m

! number of columns of A
integer (mi), intent(in)      :: n

! leading dimension of matrix Q, cannot be less than m
integer (mi), intent(in)      :: ldq

! on input the full lower trapezoid in the first min(M,N) columns of Q
! contains the factored form. On output, Q has been accumulated into
! a square matrix.
real    (mr), intent(in out)  :: q(ldq, m)

!===========================================================================
! local variables for MinpackQform
!===========================================================================

! loop counters
integer (mi)                  :: j
integer (mi)                  :: k
integer (mi)                  :: l

! minmum of m and n
integer (mi)                  :: minmn

! temporary dot ptoduct
real    (mr)                  :: temp

! working array
real    (mr)                  :: wa(m)

! DrHook variables
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'MINPACKQFORM'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackQform
!===========================================================================

! zero out upper triangle of Q in the first min(m,n) columns
minmn = min( m, n )
do j = minpack_two_int, minmn
  q(minpack_one_int : j - minpack_one_int, j) = minpack_zero_real
end do

! initialize remaining columns to those of the identity matrix
q(minpack_one_int : m, n + minpack_one_int : m) = minpack_zero_real

do j = n + minpack_one_int, m
  q(j, j) = minpack_one_real
end do

! accumulate Q from its factored form
do l = minpack_one_int, minmn

  k        = minmn - l + minpack_one_int

  wa(k:m)  = q(k:m,k)

  q( k:m, k) = minpack_zero_real
  q( k  , k) = minpack_one_real

  if ( wa(k) /= minpack_zero_real ) then

    do j = k, m
      temp      = dot_product( wa(k:m), q(k:m, j) ) / wa(k)
      q(k:m, j) = q(k:m, j) - ( temp * wa(k:m) )
    end do

  end if

end do

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine MinpackQform

subroutine MinpackQrfac ( m, n, lda, lipvt, pivot, a, ipvt, rdiag, acnorm )

  !===========================================================================
  !
  ! MinpackQrfac computes the QR factorization of a matrix A(m,n) using
  ! Householder transformations.
  !
  ! Discussion:
  !
  !    This procedure uses Householder transformations with optional column
  !    pivoting to compute a QR factorization of the M by N matrix A.
  !    That is, MinpackQrfac determines an orthogonal matrix Q, a
  !    permutation matrix P, and an upper trapezoidal matrix R with
  !    diagonal elements of nonincreasing magnitude, such that A*P = Q*R.
  !
  !
  !    The Householder transformation for column K, K = 1,2,...,min(M,N),
  !    is of the form
  !
  !      I - ( 1 / U(K) ) * U * U'
  !
  !    where U has zeros in the first K-1 positions.
  !
  !    The form of this transformation and the method of pivoting first
  !    appeared in the corresponding LINPACK routine.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackQrfac
!===========================================================================

! number of rows of A
integer (mi), intent(in)       :: m

! number of columns of A
integer (mi), intent(in)       :: n

! leading dimension of matrix A, cannot be less than m
integer (mi), intent(in)       :: lda

! dimension of the permutation matrix P. If pivot is set to false
! then lipvt may be as small as 1, otherwise lipvt must be at least n
integer (mi), intent(in)       :: lipvt

! if pivot is set true then colomn pivoting is used, otherwise if
! set to false then no pivoting is used
logical (mlogical), intent(in) :: pivot

! on input, A contains the matrix for which the QR factorization is to
! be computed.  On output, the strict upper trapezoidal part of A contains
! the strict upper trapezoidal part of R, and the lower trapezoidal
! part of A contains a factored form of Q, the non-trivial elements of
! the U vectors described above.
real    (mr), intent(in out)   :: a(lda, n)

! defines the permutation matrix P such that A*P = Q*R.
! Column j of P is column ipvt(j) of the identity matrix.
! If pivot is false, ipvt is not referenced.
integer (mi), intent(out)      :: ipvt(lipvt)

! diagonal elements of R
real    (mr), intent(out)      ::rdiag(n)

! norms of the corresponding columns of the input matrix A.
! If this information is not needed, then acnorm can coincide with rdiag.
real    (mr), intent(out)      ::acnorm(n)

!===========================================================================
! local variables for MinpackQrfac
!===========================================================================

! loop counters
integer (mi)                  :: j
integer (mi)                  :: k

! temporary store for a pivot location
integer (mi)                  :: i4_temp

! array index
integer (mi)                  :: kmax

! minimum of m and n
integer (mi)                  :: minmn

! Euclidean norms
real    (mr)                  :: ajnorm

! temporary stores
real    (mr)                  :: temp
real    (mr)                  :: r8_temp(m)
real    (mr)                  :: wa(n)

! DrHook variables
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'MINPACKQRFAC'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackQrfac
!===========================================================================

! compute the initial column norms and initialize several arrays
do j = 1, n
  acnorm(j) = MinpackEnorm ( m, a(minpack_one_int:m,j) )
end do

rdiag(minpack_one_int:n) = acnorm(minpack_one_int:n)
wa(   minpack_one_int:n) = acnorm(minpack_one_int:n)

if ( pivot ) then
  do j = minpack_one_int, n
    ipvt(j) = j
  end do
end if

! reduce A to R with Householder transformations
minmn = min( m, n )

do j = minpack_one_int, minmn

  ! bring the column of largest norm into the pivot position
  if ( pivot ) then

    kmax = j

    do k = j, n
      if ( rdiag(kmax) < rdiag(k) ) then
        kmax = k
      end if
    end do

    if ( kmax /= j ) then

      r8_temp(minpack_one_int:m) = a(minpack_one_int:m,j)
      a(minpack_one_int:m, j   ) = a(minpack_one_int:m,kmax)
      a(minpack_one_int:m, kmax) = r8_temp(minpack_one_int:m)

      rdiag(kmax)  = rdiag(j)
      wa(kmax)     = wa(j)

      i4_temp      = ipvt(j)
      ipvt(j)      = ipvt(kmax)
      ipvt(kmax)   = i4_temp

    end if

  end if

  ! compute the Householder transformation to reduce the
  ! j-th column of A to a multiple of the j-th unit vector
  ajnorm = MinpackEnorm ( m - j + minpack_one_int, a(j,j) )

  if ( ajnorm /= minpack_zero_real ) then

    if ( a(j,j) < minpack_zero_real ) then
      ajnorm = -ajnorm
    end if

    a(j:m,j) = a(j:m,j) / ajnorm
    a(j  ,j) = a(j  ,j) + minpack_one_real

    ! apply the transformation to the remaining columns and update the norms
    do k = j + minpack_one_int, n

      temp = dot_product( a(j:m,j), a(j:m,k) ) / a(j,j)

      a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

      if ( pivot .and. rdiag(k) /= minpack_zero_real ) then

        temp     = a(j,k) / rdiag(k)
        rdiag(k) = rdiag(k) *                                                  &
                   sqrt( max( minpack_zero_real,                               &
                            (minpack_one_real - temp) ** minpack_two_int ) )

        if ( minpack_p05_real *                                                &
             ( ( rdiag(k) / wa(k) ) ** minpack_two_int ) <= epsmch ) then
          rdiag(k) = MinpackEnorm ( m-j, a(j + minpack_one_int,k) )
          wa(k)    = rdiag(k)
        end if

      end if

    end do

  end if

  rdiag(j) = - ajnorm

end do

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine MinpackQrfac

subroutine MinpackAQ ( m, n, lda, v, w, a )

  !===========================================================================
  !
  ! MinpackAQ computes A*Q, where Q is the product of Householder
  ! transformations.
  !
  ! Discussion:
  !
  !    Given an M by N matrix A, this function computes A*Q where
  !    Q is the product of 2*(N - 1) transformations
  !
  !      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
  !
  !    and GV(I), GW(I) are Givens rotations in the (I,N) plane which
  !    eliminate elements in the I-th and N-th planes, respectively.
  !    Q itself is not given, rather the information to recover the
  !    GV, GW rotations is supplied.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackAQ
!===========================================================================

! number of rows of A
integer (mi), intent(in)      :: m

! number of columns of A
integer (mi), intent(in)      :: n

! leading dimension of matrix A, cannot be less than m
integer (mi), intent(in)      :: lda

! information necessary to recover the Givens rotations GV
real    (mr), intent(in)      :: v(n)

! information necessary to recover the Givens rotations GW
real    (mr), intent(in)      :: w(n)

! on input, the matrix A to be postmultiplied by the orthogonal matrix Q.
! On output, the value of A*Q.
real    (mr), intent(in out)  :: a(lda, n)

!===========================================================================
! local variables for MinpackAQ
!===========================================================================

! loop counters
integer (mi)                  :: i
integer (mi)                  :: j

! cosine
real    (mr)                  :: c

! sine
real    (mr)                  :: s

! temporary store
real    (mr)                  :: temp

! DrHook variables
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'MINPACKAQ'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackAQ
!===========================================================================

! apply the first set of Givens rotations to A
do j = n - minpack_one_int, minpack_one_int, - minpack_one_int

  if ( abs( v(j) ) > minpack_one_real ) then
    c = minpack_one_real / v(j)
    s = sqrt( minpack_one_real - c ** minpack_two_int )
  else
    s = v(j)
    c = sqrt( minpack_one_real - s ** minpack_two_int )
  end if

  do i = minpack_one_int, m
    temp   = c * a(i,j) - s * a(i,n)
    a(i,n) = s * a(i,j) + c * a(i,n)
    a(i,j) = temp
  end do

end do

! apply the second set of Givens rotations to A
do j = minpack_one_int, n - minpack_one_int

  if ( abs( w(j) ) > minpack_one_real ) then
    c = minpack_one_real / w(j)
    s = sqrt( minpack_one_real - c ** minpack_two_int )
  else
    s = w(j)
    c = sqrt( minpack_one_real - s ** minpack_two_int )
  end if

  do i = minpack_one_int, m
    temp   =   c * a(i,j) + s * a(i,n)
    a(i,n) = - s * a(i,j) + c * a(i,n)
    a(i,j) =   temp
  end do

end do

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine MinpackAQ

subroutine MinpackR1updt ( m, n, ls, u, s, v, w, sing )

  !===========================================================================
  !
  ! MinpackR1updt re-triangularizes a matrix after a rank one update.
  !
  ! Discussion:
  !
  !    Given an M by N lower trapezoidal matrix S, an M-vector U, and an
  !    N-vector V, the problem is to determine an orthogonal matrix Q such
  !    that
  !
  !      (S + U * V' ) * Q
  !
  !    is again lower trapezoidal.
  !
  !    This function determines Q as the product of 2 * (N - 1)
  !    transformations
  !
  !      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
  !
  !    where GV(I), GW(I) are Givens rotations in the (I,N) plane
  !    which eliminate elements in the I-th and N-th planes,
  !    respectively.  Q itself is not accumulated, rather the
  !    information to recover the GV and GW rotations is returned.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackR1updt
!===========================================================================

! number of rows of S
integer (mi), intent(in)        :: m

! number of columns of S
integer (mi), intent(in)        :: n

! length of the S array, ls must be at least (n*(2*m-n+1))/2
integer (mi), intent(in)        :: ls

! m-vector U
real    (mr), intent(in)        :: u(m)

! on input, the lower trapezoidal matrix S stored by columns.
! On output S contains the lower trapezoidal matrix produced as
! described above.
real    (mr), intent(in out)    :: s(ls)

! on input, V must contain the n-vector V.
! On output V contains the information necessary to recover the
! Givens rotations GV described above.
real    (mr), intent(in out)    :: v(n)

! information necessary to recover the Givens rotations GW described above
real    (mr), intent(out)       :: w(m)

! set to true if any of the diagonal elements of the output S are zero.
! Otherwise SING is set false.
logical (mlogical), intent(out) :: sing

!===========================================================================
! local variables for MinpackR1updt
!===========================================================================

! loop counters
integer (mi)                    :: i
integer (mi)                    :: j
integer (mi)                    :: jj
integer (mi)                    :: l

! cosine
real    (mr)                    :: cs

! sine
real    (mr)                    :: sn

! tangent
real    (mr)                    :: tn

! cotangent
real    (mr)                    :: ct

! 1 / cosine
real    (mr)                    :: tau

real    (mr)                    :: temp

! DrHook variables
real    (kind=jprb)             :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter  :: RoutineName = 'MINPACKR1UPDT'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackR1updt
!===========================================================================

! initialise the logical that flags zero diagonal elements in the output S
sing = .false.

! initialize the diagonal element pointer
jj = ( n * ( minpack_two_int * m - n + minpack_one_int ) ) /                   &
     minpack_two_int - ( m - n )

! move the nontrivial part of the last column of S into W
l = jj
do i = n, m
  w(i) = s(l)
  l    = l + minpack_one_int
end do

! rotate the vector V into a multiple of the N-th unit vector
! in such a way that a spike is introduced into W
do j = n - minpack_one_int, minpack_one_int, - minpack_one_int

  jj   = jj - ( m - j + minpack_one_int )
  w(j) = minpack_zero_real

  if ( v(j) /= minpack_zero_real ) then

    ! determine a Givens rotation which eliminates the
    ! J-th element of V.
    if ( abs( v(n) ) < abs ( v(j) ) ) then
      ct    = v(n) / v(j)
      sn    = minpack_p5_real /                                                &
              sqrt(minpack_p25_real + minpack_p25_real * ct ** minpack_two_int)
      cs    = sn * ct
      tau   = minpack_one_real
      if ( abs( cs ) * giant > minpack_one_real ) then
        tau = minpack_one_real / cs
      end if
    else
      tn    = v(j) / v(n)
      cs    = minpack_p5_real /                                                &
              sqrt(minpack_p25_real + minpack_p25_real * tn ** minpack_two_int)
      sn    = cs * tn
      tau   = sn
    end if

    ! apply the transformation to V and store the information
    ! necessary to recover the Givens rotation GV
    v(n) = sn * v(j) + cs * v(n)
    v(j) = tau

    ! apply the transformation to S and extend the spike in W
    l = jj
    do i = j, m
      temp = cs * s(l) - sn * w(i)
      w(i) = sn * s(l) + cs * w(i)
      s(l) = temp
      l = l + minpack_one_int
    end do

  end if

end do

! add the spike from the rank 1 update to W
w(minpack_one_int:m) = w(minpack_one_int:m) + v(n) * u(minpack_one_int:m)

! eliminate the spike
do j = minpack_one_int, n - minpack_one_int

  if ( w(j) /= minpack_zero_real ) then

    ! determine a Givens rotation which eliminates the
    ! J-th element of the spike
    if ( abs ( s(jj) ) < abs ( w(j) ) ) then

      ct    = s(jj) / w(j)
      sn    = minpack_p5_real /                                                &
              sqrt(minpack_p25_real + minpack_p25_real * ct ** minpack_two_int)
      cs    = sn * ct

      if ( abs( cs ) * giant > minpack_one_real) then
        tau = minpack_one_real / cs
      else
        tau = minpack_one_real
      end if

    else

      tn  = w(j) / s(jj)
      cs  = minpack_p5_real /                                                  &
            sqrt( minpack_p25_real + minpack_p25_real * tn ** minpack_two_int )
      sn  = cs * tn
      tau = sn

    end if

    ! apply the transformation to S and reduce the spike in W
    l = jj
    do i = j, m
      temp =   cs * s(l) + sn * w(i)
      w(i) = - sn * s(l) + cs * w(i)
      s(l) =   temp
      l = l + minpack_one_int
    end do

    ! store information necessary to recover the Givens rotation GW
    w(j) = tau

  end if

  ! test for zero diagonal elements in the output S
  if ( s(jj) == minpack_zero_real ) then
    sing = .true.
  end if

  jj = jj + ( m - j + minpack_one_int )

end do

! move W back into the last column of the output S
l = jj
do i = n, m
  s(l) = w(i)
  l    = l + minpack_one_int
end do

if ( s(jj) == minpack_zero_real ) then
  sing = .true.
end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine MinpackR1updt

subroutine MinpackFdJacNbyN ( UserFunction, n, ldfjac, ml, mu,                 &
                              x, fn_x, epsilon_min, iflag, fjac )

  !===========================================================================
  !
  ! MinpackFdJacNbyN estimates an n by n Jacobian matrix using forward
  !                  differences.
  !
  ! Discussion:
  !
  !    This subroutine computes a forward-difference approximation
  !    to the N by N Jacobian matrix associated with a specified
  !    problem of N functions in N variables. If the Jacobian has
  !    a banded form, then function evaluations are saved by only
  !    approximating the nonzero terms.
  !
  !===========================================================================

implicit none

!===========================================================================
! arguments for MinpackFdJacNbyN
!===========================================================================

! define the procedure UserFunction using an abstract interface
procedure(MinpackUserFunction):: UserFunction

 ! number of columns of functions and variables in the nonlinear problem
integer (mi), intent(in)      :: n

! leading diemsion of array fjac, must not be less than n
integer (mi), intent(in)      :: ldfjac

! number of subdiagonals within the band of the Jacobian matrix. If the
! matrix is not banded, then set ml to at least n - 1.
integer (mi), intent(in)      :: ml

! number of superdiagonals within the band of the Jacobian matrix. If the
! matrix is not banded, then set mu to at least n - 1.
integer (mi), intent(in)      :: mu

! the point where the Jacobian is to be evaluated
real    (mr), intent(in)      :: x(n)

! the n problem functions evaluated at x
real    (mr), intent(in out)  :: fn_x(n)

! used in determining a suitable step length for the
! forward-difference approximation. This approximation assumes that the
! relative errors in the functions are of the order of epsilon_min.
! If epsilon_min is less than the machine precision, it is assumed that
! the relative errors in the functions are of the order of the machine
! precision.
real    (mr), intent(in)      :: epsilon_min

! flag that will terminate execution of MinpackFdJacNbyN if set to
! a negative integer
integer (mi), intent(out)     :: iflag

! approximation to the n by n Jacobian matrix
real    (mr), intent(out)     :: fjac(ldfjac,n)

!===========================================================================
! local variables for MinpackFdJacNbyN
!===========================================================================

! loop counters
integer (mi)                  :: i
integer (mi)                  :: j
integer (mi)                  :: k

integer (mi)                  :: msum

! dummy argument used to call UserFunction
integer (mi)                  :: dummy_ldfjac

! step length scaling to use
real    (mr)                  :: eps

! step length
real    (mr)                  :: h

! an element of the vector x
real    (mr)                  :: element_x

! problem function evluated at x + h
real    (mr)                  :: fn_x_h(n)

! copy of vector x
real    (mr)                  :: copy_x(n)

! dummy argument used to call UserFunction
real    (mr)                  :: dummy_fjac(ldfjac, n)

! DrHook variables
real    (kind=jprb)           :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*) , parameter :: RoutineName = 'MINPACKFDJACNBYN'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for MinpackFdJacNbyN
!===========================================================================

! set the dummy arguments for UserFunction
dummy_ldfjac    = n
dummy_fjac(:,:) = minpack_zero_real

! compute the step length scaling to use
eps  = sqrt( max( epsilon_min, epsmch ) )

! how many bands in the Jacobian
msum = ml + mu + minpack_one_int

! check to determine if the Jacobian is dense or banded
if ( msum >= n ) then

  ! make a copy of the input vector x
  copy_x = x

  ! computation of dense approximate Jacobian
  do j = minpack_one_int, n

    element_x = copy_x(j)

    ! compute step length h
    h = eps * abs( element_x )
    if ( h == minpack_zero_real ) then
      h = eps
    end if

    ! compute problem functions at x + h
    iflag = minpack_parameters % iflag_function
    copy_x(j)  = element_x + h
    call UserFunction ( n, dummy_ldfjac, copy_x, iflag, fn_x_h, dummy_fjac )

    if ( iflag < minpack_zero_int ) then
      ! problem calling UserFunction, so exit from this procedure
      exit
    end if

    ! compute Jacobian at x
    fjac(minpack_one_int:n,j) = ( fn_x_h( minpack_one_int:n) -                 &
                                  fn_x(   minpack_one_int:n)   ) / h
    copy_x(j)                 = element_x

  end do

else

  ! computation of banded approximate Jacobian,
  ! start by making a copy of the input vector x
  copy_x = x

  do k = minpack_one_int, msum

    do j = k, n, msum

      ! compute step length h
      h = eps * abs( x(j) )
      if ( h == minpack_zero_real ) then
        h = eps
      end if

      copy_x(j) = x(j) + h
    end do

    ! compute the function at x + h
    iflag = minpack_parameters % iflag_function
    call UserFunction ( n, dummy_ldfjac, copy_x, iflag, fn_x_h, dummy_fjac )

    if ( iflag < minpack_zero_int ) then
      exit
    end if

    do j = k, n, msum

      copy_x(j) = x(j)

      ! compute step length h
      h = eps * abs( x(j) )
      if ( h == minpack_zero_real ) then
        h = eps
      end if

      fjac(minpack_one_int:n,j) = minpack_zero_real

      do i = minpack_one_int, n
        if ( (i >= j - mu) .and. (i <= j + ml) ) then
          fjac(i,j) = ( fn_x_h(i) - fn_x(i) ) / h
        end if
      end do

    end do

    ! end of computing the banded approximate Jacobian
  end do

end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine MinpackFdJacNbyN

end module vera_minpack_mod
