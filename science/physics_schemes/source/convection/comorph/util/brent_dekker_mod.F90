! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module brent_dekker_mod

implicit none

! This module contains code to iteratively find the root f(x)=0
! where x is some input and f is the result of a subroutine which
! is passed in as an argument.  It implements the well-known
! Brent-Dekker root-finding algorithm.
! The calculation within the passed-in subroutine may have
! additional input arguments other than x;
! these can be passed in through brent_dekker_solve by
! "packaging them up" in a super-array to pass into
! f(x) itself.
! The subroutine passed in (which calculates f(x) for a given x)
! must have its argument list defined as in the interface
! f_of_x_interface below.
! Subroutine brent_dekker_solve below takes as input the
! subroutine to calculate f(x) and 2 initial guesses for x,
! which must bracket the root to be found.
! After each new guess x has been set, the input subroutine is
! called (passing in the super-array containing its arguments)
! to compute the new value of f(x).
! The iteration repeats until either f(x) is near enough zero to
! meet the convergence criteria, or the specified maximum number
! of iterations is reached.


!----------------------------------------------------------------
! Interface to the subroutine which is passed into
! brent_dekker_solve, to calculate the error f(x) for a given x
!----------------------------------------------------------------
abstract interface
  subroutine f_of_x_interface( nc, n_points,                                   &
                               n_real_sca, n_real_arr,                         &
                               args_real_sca, args_real_arr,                   &
                               x, f )
  use comorph_constants_mod, only: real_cvprec
  implicit none

  ! Number of points to operate on
  integer, intent(in) :: nc

  ! Size of arrays
  ! (maybe larger than the actual number of points, as the
  !  same arrays are re-used even when some points have converged
  !  and so have been removed from the list of points).
  integer, intent(in) :: n_points

  ! Number of each sort of argument
  integer, intent(in) :: n_real_sca   ! in real scalar arguments
  integer, intent(in) :: n_real_arr   ! in real array arguments

  ! List of scalar inputs (constants)
  real(kind=real_cvprec), intent(in) :: args_real_sca                          &
                                                  ( n_real_sca )
  ! Super-array containing array arguments
  real(kind=real_cvprec), intent(in) :: args_real_arr                          &
                                        ( n_points, n_real_arr )

  ! Input value of x to calc f
  real(kind=real_cvprec), intent(in) :: x(n_points)
  ! Output value of f(x)
  real(kind=real_cvprec), intent(out) :: f(n_points)

  end subroutine f_of_x_interface
end interface


contains


!----------------------------------------------------------------
! Main routine to find the value of x yielding f(x) = 0
!----------------------------------------------------------------
! takes as input a subroutine f_of_x, which calculates f(x).
subroutine brent_dekker_solve( n_points, n_real_sca, n_real_arr,               &
                               args_real_sca, args_real_arr,                   &
                               f_of_x, n_iter, tolerance, i_check_converge,    &
                               x_a, f_a, x_b, f_b, x_out )

use comorph_constants_mod, only: real_cvprec, newline,                         &
                                 zero, quarter, half, one, three,              &
                                 min_delta, sqrt_min_float,                    &
                                 i_check_bad_none, i_check_bad_warn,           &
                                 i_check_bad_fatal
use raise_error_mod, only: raise_fatal, raise_warning

implicit none

! Number of points
integer, intent(in) :: n_points

! Arguments to pass through to the input subroutine f_of_x:

! Number of each sort of argument
integer, intent(in) :: n_real_sca   ! in real scalar arguments
integer, intent(in) :: n_real_arr   ! in real array arguments

! List of scalar inputs (constants)
real(kind=real_cvprec), intent(in) ::    args_real_sca                         &
                                                   ( n_real_sca )
! Super-array containing array arguments
real(kind=real_cvprec), intent(in out) :: args_real_arr                        &
                                         ( n_points, n_real_arr )
! Note: the array arguments are intent(inout) even though
! this routine does not modify the values of these arguments;
! this is because it needs to rearrange the data within the
! arrays during the iteration, to make the non-converged points
! contiguous in memory.
! Any data that needs to be used after the call to this routine
! should not be passed in here, as it will get scrambled!


! Subroutine passed in to compute f(x).
procedure (f_of_x_interface) :: f_of_x
! Must have arguments arranged as in the interface block above.
! Arguments in addition to x, f and the array sizes must be
! packaged up into args_real_sca (scalar inputs) and
! args_real_arr (array inputs), declared above.


! Max number of iterations allowed
integer, intent(in) :: n_iter

! Tolerance for convergence test
real(kind=real_cvprec), intent(in) :: tolerance

! Switch indicating what to do if any points have not met the
! convergence criteria after n_iter iterations
! (0 do nothing, 1 print a warning, 2 raise a fatal error).
integer, intent(in) :: i_check_converge

! Initial 2 guesses for x, and their corresponding values of f.
! These must be set to 2 different values, such that
! f_a and f_b have opposite signs (ie the root is bracketed).
! Where these initial guesses do not bracket the root,
! an error will be raised.
real(kind=real_cvprec), intent(in out) :: x_a(n_points)
real(kind=real_cvprec), intent(in out) :: f_a(n_points)
real(kind=real_cvprec), intent(in out) :: x_b(n_points)
real(kind=real_cvprec), intent(in out) :: f_b(n_points)
! Note that these arrays get rearranged during compression
! onto non-converged points, so must not store data that's
! needed afterwards!

! Output final converged values of x
real(kind=real_cvprec), intent(out) :: x_out(n_points)


! Work arrays...

! Additional guesses for x and their errors
real(kind=real_cvprec) :: x_n(n_points)
real(kind=real_cvprec) :: f_n(n_points)
real(kind=real_cvprec) :: x_c(n_points)
real(kind=real_cvprec) :: f_c(n_points)
real(kind=real_cvprec) :: x_d(n_points)

! temporary storage for swapping values
real(kind=real_cvprec) :: tmp

! Flag for whether the algorithm used bisection on the
! previous iteration
logical :: l_bisected_last_step(n_points)

! Flag for bisection at the current step
logical :: l_bisect

! Updated number of points to keep iterating
integer :: nc
! Indices of points where we keep iterating
integer :: index_ic(n_points)

! Iteration count
integer :: iter

! Strings for storing values printed to error message when it
! all goes horribly wrong.
character(len=30+14*n_real_sca) :: args_real_sca_string
character(len=30+14*n_real_arr) :: args_real_arr_string
character(len=30+14*4) :: ab_string
character(len=13) :: write_fmt

! Small number for checking that 2 numbers are
! noticably different
real(kind=real_cvprec), parameter :: small_diff = 10.0_real_cvprec * min_delta

! Threshold ratio of f(x_a) / f(x_b) beyond-which to stop the
! scheme bisecting twice in a row
! (modification of the Brent-Dekker algorithm to prevent it
!  getting stuck in a rut repeatedly bisecting when x_b is much
!  closer to the root than x_a).
real(kind=real_cvprec), parameter :: double_bisect_thresh = 4.0_real_cvprec

! Loop counters
integer :: i_field, ic

character(len=*), parameter :: routinename = "BRENT_DECKER_SOLVE"


! Initialise number of points
nc = n_points
! Initialise indices of non-converged points to be all points
do ic = 1, nc
  index_ic(ic) = ic
end do

! Check that the initial guesses a and b bracket the root
do ic = 1, nc
  if ( f_a(ic) * f_b(ic) > zero ) then
    ! If the guess's errors are both positive or both negative,
    ! the root is not bracketed; raise a fatal error

    ! Print the arguments at the first non-bracketed point
    if ( n_real_sca > 0 ) then
      write(write_fmt,"(A3,I3,A7)") "(A,", n_real_sca, "ES14.6)"
      write(args_real_sca_string,write_fmt) "args_real_sca:",                  &
           (args_real_sca(i_field), i_field=1,n_real_sca)
    else
      write(args_real_sca_string,"(A)") "args_real_sca: (none)"
    end if

    if ( n_real_arr > 0 ) then
      write(write_fmt,"(A3,I3,A7)") "(A,", n_real_arr, "ES14.6)"
      write(args_real_arr_string,write_fmt) "args_real_arr:",                  &
           (args_real_arr(1,i_field), i_field=1,n_real_arr)
    else
      write(args_real_arr_string,"(A)") "args_real_arr: (none)"
    end if

    write(ab_string,"(A,4ES14.6)") "x_a,f_a,x_b,f_b:",                         &
          x_a(ic), f_a(ic), x_b(ic), f_b(ic)

    ! Fatal error, printing the values at the problem point:
    call raise_fatal( routinename,                                             &
           "The input initial guesses don't bracket the root, "              //&
           "(errors of the 2 initial guesses have same sign), "     //newline//&
           "so no solution can be found."                           //newline//&
           trim(adjustl(args_real_sca_string))                      //newline//&
           trim(adjustl(args_real_arr_string))                      //newline//&
           trim(adjustl(ab_string)) )

  end if  ! ( f_a(ic) * f_b(ic) > zero )
end do  ! ic = 1, nc

! Check that guesses a and b are the right way round;
! b should have the smaller error
do ic = 1, nc
  if ( abs(f_a(ic)) < abs(f_b(ic)) ) then
    ! Swap a and b
    tmp     = x_a(ic)
    x_a(ic) = x_b(ic)
    x_b(ic) = tmp
    tmp     = f_a(ic)
    f_a(ic) = f_b(ic)
    f_b(ic) = tmp
  end if
end do

! Initialise guess c, and guess d just to be safe
do ic = 1, nc
  x_c(ic) = x_a(ic)
  f_c(ic) = f_a(ic)
  x_d(ic) = x_a(ic)
end do

! Initialise bisection flag to true for first iteration
do ic = 1, nc
  l_bisected_last_step(ic) = .true.
end do

! Check whether the initial guess meets the convergence
! criteria; if it does, remove such points from the list
! so that we don't do any iterations at such points
call convergence_test( n_points, n_real_arr,                                   &
                       tolerance, nc, index_ic,                                &
                       x_a, f_a, x_b, f_b, x_c, f_c, x_d,                      &
                       l_bisected_last_step,                                   &
                       args_real_arr, x_out )


! Begin iterating
iter = 0
do while ( nc > 0 .and. iter < n_iter )
  iter = iter + 1


  if ( iter == 1 ) then
    ! On first iteration, set next guess using
    ! linear interpolation, as guess c = guess a
    do ic = 1, nc
      x_n(ic) = x_b(ic) - f_b(ic) * ( ( x_b(ic) - x_a(ic) )                    &
                                    / ( f_b(ic) - f_a(ic) ) )
    end do
  else  ! ( iter > 1 )
    ! On later iterations, try using inverse quadratic
    ! interpolation...
    do ic = 1, nc
      if ( abs(f_c(ic) - f_a(ic)) > sqrt_min_float .and.                       &
           abs(f_c(ic) - f_b(ic)) > sqrt_min_float ) then
        ! Inverse quadratic formula is only safe if f(x_c) is not equal
        ! to f(x_a) or f(x_b) (and not within a small numerical
        ! tolerance of them)
        x_n(ic) = x_a(ic) * f_b(ic) * f_c(ic)                                  &
             / ( ( f_a(ic) - f_b(ic) ) * ( f_a(ic) - f_c(ic) ) )               &
              + x_b(ic) * f_c(ic) * f_a(ic)                                    &
             / ( ( f_b(ic) - f_c(ic) ) * ( f_b(ic) - f_a(ic) ) )               &
              + x_c(ic) * f_a(ic) * f_b(ic)                                    &
             / ( ( f_c(ic) - f_a(ic) ) * ( f_c(ic) - f_b(ic) ) )
      else
        ! Where f(x_c) too close to one of the others, use
        ! linear interpolation instead
        x_n(ic) = x_b(ic) - f_b(ic) * ( ( x_b(ic) - x_a(ic) )                  &
                                      / ( f_b(ic) - f_a(ic) ) )
      end if
    end do
  end if  ! ( iter > 1 )

  ! Find points where the interpolated value is not optimal and
  ! needs to be replaced by simple bisection
  do ic = 1, nc
    tmp = ( three*x_a(ic) + x_b(ic) )*quarter
    l_bisect = .false.
    ! If new x not between x_b and the limit (stored in tmp)...
    if ( .not. ( x_n(ic) > min( tmp, x_b(ic) ) .and.                           &
                 x_n(ic) < max( tmp, x_b(ic) )       ) ) then
      if ( abs(x_n(ic)-x_b(ic)) < min_delta * abs(x_n(ic)) ) then
        ! Special case where x_n = x_b, just due to rounding error in the
        ! interpolation, which causes the algorithm to get stuck in a rut.
        ! Just nudge x_n slightly towards x_a when this happens
        ! (this is an extension to the published Brent-Dekker algorithm):
        if ( x_a(ic) > x_b(ic) ) then
          x_n(ic) = x_n(ic) * (one + min_delta)
        else
          x_n(ic) = x_n(ic) * (one - min_delta)
        end if
      else
        l_bisect = .true.
      end if
    else  ! New x between x_b and limit
      if ( l_bisected_last_step(ic) ) then
        ! or we did bisection last step, and...
        if (     abs( x_n(ic) - x_b(ic) )                                      &
              >= abs( x_b(ic) - x_c(ic) )*half                                 &
            .or. abs( x_b(ic) - x_c(ic) )                                      &
              <  small_diff * max( abs(x_b(ic)), abs(x_c(ic)) )                &
           ) then
          ! Extra condition; don't bisect twice in a row if
          ! a is an outlier (avoids getting stuck in a rut
          ! when b, c are close together but a is miles off)
          if ( abs(f_a(ic)) < double_bisect_thresh * abs(f_b(ic)) ) then
            l_bisect = .true.
          end if
          ! Note: this check to sometimes avoid bisecting twice
          ! in a row is a departure from the published Brent-Dekker
          ! algorithm.  It can greatly reduce the number of
          ! iterations required, provided that f(x) is reasonably
          ! well-behaved.
        end if
      else  ! ( .not. l_bisected_last_step(ic) )
        ! or we didn't do bisection last step, and...
        if (     abs( x_n(ic) - x_b(ic) )                                      &
              >= abs( x_c(ic) - x_d(ic) )*half                                 &
            .or. abs( x_c(ic) - x_d(ic) )                                      &
              <  small_diff * max( abs(x_c(ic)), abs(x_d(ic)) )                &
           ) then
          l_bisect = .true.
        end if
      end if  ! ( .not. l_bisected_last_step(ic) )
    end if  ! New x between x_b and limit
    ! If bisection criteria are satisfied, overwrite x_n
    ! with bisection value and set flag
    if ( l_bisect ) then
      x_n(ic) = half * ( x_a(ic) + x_b(ic) )
      l_bisected_last_step(ic) = .true.
    else
      l_bisected_last_step(ic) = .false.
    end if
  end do  ! ic = 1, nc


  ! Calculate error for the new guess using the input subroutine
  call f_of_x( nc, n_points, n_real_sca, n_real_arr,                           &
               args_real_sca, args_real_arr, x_n, f_n )


  ! Copy guesses b and c backwards, overwriting d
  do ic = 1, nc
    x_d(ic) = x_c(ic)
    x_c(ic) = x_b(ic)
    f_c(ic) = f_b(ic)
  end do
  ! Copy guess n into either a or b, so-as to retain bracketing
  ! of the root
  do ic = 1, nc
    if ( f_a(ic) * f_n(ic) < zero ) then
      x_b(ic) = x_n(ic)
      f_b(ic) = f_n(ic)
    else
      x_a(ic) = x_n(ic)
      f_a(ic) = f_n(ic)
    end if
  end do

  ! Ensure the new 'a' and 'b' are the right way round
  ! ('b' should have the smaller residual error f)
  do ic = 1, nc
    if ( abs(f_a(ic)) < abs(f_b(ic)) ) then
      ! Swap a and b
      tmp     = x_a(ic)
      x_a(ic) = x_b(ic)
      x_b(ic) = tmp
      tmp     = f_a(ic)
      f_a(ic) = f_b(ic)
      f_b(ic) = tmp
    end if
  end do

  ! Test convergence criteria:
  call convergence_test( n_points, n_real_arr,                                 &
                         tolerance, nc, index_ic,                              &
                         x_a, f_a, x_b, f_b, x_c, f_c, x_d,                    &
                         l_bisected_last_step,                                 &
                         args_real_arr, x_out )


end do  ! while ( nc > 0 .and. iter < n_iter )


! If any points haven't converged, scatter final values back
! to output...
if ( nc > 0 ) then

  do ic = 1, nc
    x_out(index_ic(ic)) = x_b(ic)
  end do

  if ( i_check_converge > i_check_bad_none ) then

    ! Print the arguments at the first non-converged point
    if ( n_real_sca > 0 ) then
      write(write_fmt,"(A3,I3,A7)") "(A,", n_real_sca, "ES14.6)"
      write(args_real_sca_string,write_fmt) "args_real_sca:",                  &
           (args_real_sca(i_field), i_field=1,n_real_sca)
    else
      write(args_real_sca_string,"(A)") "args_real_sca: (none)"
    end if

    if ( n_real_arr > 0 ) then
      write(write_fmt,"(A3,I3,A7)") "(A,", n_real_arr, "ES14.6)"
      write(args_real_arr_string,write_fmt) "args_real_arr:",                  &
           (args_real_arr(1,i_field), i_field=1,n_real_arr)
    else
      write(args_real_arr_string,"(A)") "args_real_arr: (none)"
    end if

    write(ab_string,"(A,4ES14.6)") "x_a,f_a,x_b,f_b:",                         &
          x_a(1), f_a(1), x_b(1), f_b(1)

    select case ( i_check_converge )
    case ( i_check_bad_warn )
      call raise_warning( routinename,                                         &
         "Brent-Dekker has failed to converge on a solution " //               &
         "at some points."                           //newline//               &
         trim(adjustl(args_real_sca_string))         //newline//               &
         trim(adjustl(args_real_arr_string))         //newline//               &
         trim(adjustl(ab_string)) )
    case ( i_check_bad_fatal )
      call raise_fatal( routinename,                                           &
         "Brent-Dekker has failed to converge on a solution " //               &
         "at some points."                           //newline//               &
         trim(adjustl(args_real_sca_string))         //newline//               &
         trim(adjustl(args_real_arr_string))         //newline//               &
         trim(adjustl(ab_string)) )
    end select

  end if  ! ( i_check_converge > i_check_bad_none )

end if  ! ( nc > 0 )


return
end subroutine brent_dekker_solve


!----------------------------------------------------------------
! Routine to test whether convergence has been achieved
!----------------------------------------------------------------
! Also recompresses the work arrays to keep the non-converged
! points contiguous in memory, and scatters output data from
! converged points back into the output full array x_out
subroutine convergence_test( n_points, n_real_arr,                             &
                             tolerance, nc, index_ic,                          &
                             x_a, f_a, x_b, f_b, x_c, f_c, x_d,                &
                             l_bisected_last_step,                             &
                             args_real_arr, x_out )

use comorph_constants_mod, only: real_cvprec, min_delta, sqrt_min_float

implicit none

! Number of points in the full arrays
integer, intent(in) :: n_points

! Number of input real array arguments for f_of_x
integer, intent(in) :: n_real_arr

! Fractional tolerance for convergence criterion
real(kind=real_cvprec), intent(in) :: tolerance

! Number of points where further iterations are needed
! (updated by this routine based on how many points in l_continue
!  are still set to true)
integer, intent(in out) :: nc

! Indices of points in the list
integer, intent(in out) :: index_ic(n_points)

! Guess values of x and their errors
real(kind=real_cvprec), intent(in out) :: x_a(n_points)
real(kind=real_cvprec), intent(in out) :: f_a(n_points)
real(kind=real_cvprec), intent(in out) :: x_b(n_points)
real(kind=real_cvprec), intent(in out) :: f_b(n_points)
real(kind=real_cvprec), intent(in out) :: x_c(n_points)
real(kind=real_cvprec), intent(in out) :: f_c(n_points)
real(kind=real_cvprec), intent(in out) :: x_d(n_points)

! Flags for where bisection was done last step
logical, intent(in out) :: l_bisected_last_step(n_points)

! Array arguments input to f_of_x
real(kind=real_cvprec), intent(in out) :: args_real_arr                        &
                                         ( n_points, n_real_arr )

! Full array of final output value of x
real(kind=real_cvprec), intent(in out) :: x_out(n_points)

! Array of flags indicating where iterations should continue
logical :: l_continue(nc)

! Index of first terminated point in the list
integer :: ic_first

! Store for indices of points being copied:
! Points where convergence reached
integer :: n_converged
integer :: index_converged(nc)
! Points not yet converged, where iterations must continue
integer :: n_continue
integer :: index_continue(nc)

! Small number for checking that 2 numbers are
! noticably different
real(kind=real_cvprec), parameter :: small_diff = 10.0_real_cvprec * min_delta

! Loop counters
integer :: ic, ic2, i_field


! Test for convergence; l_continue is set to false
! at points where any of the convergence criteria are met.
do ic = 1, nc
  l_continue(ic) =                                                             &
    ! Difference of errors too small for safe division
    abs( f_b(ic) - f_a(ic) ) > sqrt_min_float                                  &
      .and.                                                                    &
    ! Two guesses too close together
    abs( x_b(ic) - x_a(ic) ) > small_diff * abs(x_b(ic))                       &
      .and.                                                                    &
    ! Error of current best guess within tolerance
    abs(f_b(ic)) > tolerance * abs(x_b(ic))
end do

! See if any points have terminated
ic_first = 0
over_nc: do ic = 1, nc
  if ( .not. l_continue(ic) ) then
    ! Store index of first terminated point and exit the do-loop
    ic_first = ic
    exit over_nc
  end if
end do over_nc

! Only anything else to do if at least one point has terminated
if ( ic_first > 0 ) then

  ! Initialise lists of points where l_continue is true / false
  n_continue = ic_first - 1
  n_converged = 1
  index_converged(1) = ic_first

  ! Store indices of any converged / non-converged points occuring
  ! after ic_first
  if ( ic_first < nc ) then
    do ic = ic_first+1, nc
      if ( l_continue(ic) ) then
        n_continue = n_continue + 1
        index_continue(n_continue) = ic
      else
        n_converged = n_converged + 1
        index_converged(n_converged) = ic
      end if
    end do
  end if

  ! Copy compressed output value of x back to full array
  do ic2 = 1, n_converged
    ic = index_converged(ic2)
    x_out(index_ic(ic)) = x_b(ic)
  end do

  ! If any non-converged points were found after ic_first
  if ( n_continue >= ic_first ) then

    ! In this case, the data for non-converged points later
    ! in the list needs to be shuffled along to make
    ! all the non-converged points contiguous in memory.
    ! Copy all input and work array values forwards
    ! as appropriate
    do ic2 = ic_first, n_continue
      ic = index_continue(ic2)
      index_ic(ic2) = index_ic(ic)
      x_a(ic2) = x_a(ic)
      f_a(ic2) = f_a(ic)
      x_b(ic2) = x_b(ic)
      f_b(ic2) = f_b(ic)
      x_c(ic2) = x_c(ic)
      f_c(ic2) = f_c(ic)
      x_d(ic2) = x_d(ic)
      l_bisected_last_step(ic2) = l_bisected_last_step(ic)
    end do
    do i_field = 1, n_real_arr
      do ic2 = ic_first, n_continue
        ic = index_continue(ic2)
        args_real_arr(ic2,i_field) = args_real_arr(ic,i_field)
      end do
    end do

  end if  ! ( n_continue >= ic_first )

  ! Set final new number of non-converged points
  nc = n_continue

end if  ! ( ic_first > 0 )

return
end subroutine convergence_test


end module brent_dekker_mod
