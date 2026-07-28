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
! This module contains the work-horse routines used for
! casting phantom aerosol populations.
!
! Routines in this module are:
!
! (1) vera_lin_spacing           computes linearly spaced points
!
! (2) vera_lin_intervals         computes linearly spaced intervals
!
! (3) vera_log_spacing           computes a logarithmically spaced points
!
! (4) vera_log_normal            computes the normal curve on logarithmically
!                                spaced values
!
! (5) vera_log_normal_distn      computes a logarithmically spaced
!                                gaussian distribution
!
! (6) vera_triangle              computes triangular pdf
!
! (7) vera_linear_flat_distn     computes a flat linear distribution
!

module vera_phantom_tools_mod

use vera_kind_mod,   only: wp => vera_real, wi => vera_integer

use vera_global_mod, only: vera_constants

implicit none

! Description:
!   This module contains the work-horse subroutines used for
!   casting phantom aerosol populations. These routines are
!   invoked by the subroutines in the module vera_phantom_cast_mod
!   that cast phantom aerosol populations from the MURK aerosol field.
!
! Method:
!   The seven suboutines in this module are:
!
!     vera_lin_spacing
!       Computes linearly spaced points.
!
!     vera_lin_intervals
!       Computes linearly spaced intervals.
!
!     vera_log_spacing
!       Computes a logarithmically spaced points.
!
!     vera_log_normal
!       Computes the normal curve on logarithmically spaced values.
!
!     vera_log_normal_distn
!       Computes a logarithmically spaced gaussian distribution.
!
!     vera_triangle
!       Computes triangular pdf.
!
!     vera_linear_flat_distn
!       Computes a flat linear distribution.
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_PHANTOM_TOOLS_MOD'

private

! make all the subroutines in this module Public,
! they are all called by vera_phantom_cast_mod.F90
public :: vera_lin_spacing        ,                                            &
          vera_lin_intervals      ,                                            &
          vera_log_spacing        ,                                            &
          vera_log_normal         ,                                            &
          vera_log_normal_distn   ,                                            &
          vera_triangle           ,                                            &
          vera_linear_flat_distn


contains

  !=============================================================================
  !
  ! vera_lin_intervals -  Computes linearly spaced intervals, so that
  !
  ! there will be n_x points, linearly spaced across the range [x_min, x_max].
  !
  ! These points are the mid-points of the n intervals that the line is split
  ! into, i.e. the points are
  !
  ! ( x_min + (x_width/2) ) + ( m * x_width ) with m = 0, 1, 2, ...,  n-2, n-1
  !
  ! and x_width = (x_max - m_min)/n
  !
  !=============================================================================

subroutine vera_lin_intervals( n_x, x_max, x_min, x )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_lin_intervals
!
!===========================================================================

! number of points x
integer (wi), intent(in)  :: n_x

! maxmimum value of x
real    (wp), intent(in)  :: x_max

! minmimum value of x
real    (wp), intent(in)  :: x_min

! computed linearly spaced x values
real    (wp), intent(out) :: x(1:n_x)

!===========================================================================
! local variables for vera_lin_intervals
!===========================================================================

! counter to loop over the distribution
integer (wi)              :: ii

! width of each interval
real    (wp)              :: x_width

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_LIN_INTERVALS'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for vera_lin_intervals
!===========================================================================

! compute the x values - these are linearly spaced
!                        across the range [x_min, x_max]
!
! start by initialising the x values to integer values, counting up from 0
do ii = vera_constants%one_i, n_x

  x(ii) = ii - vera_constants%one_i

end do

! compute the interval width
x_width = (x_max - x_min) / size( x, kind = wp )

! scale the x to be in the range [x_min, x_max]
x = ( x * x_width ) + x_min + ( x_width * vera_constants%half )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_lin_intervals

!=============================================================================
!
! vera_lin_spacing -  Computes linearly spaced points, so that
!
! there will be n_x points, linearly spaced across the range [x_min, x_max].
!
!=============================================================================

subroutine vera_lin_spacing( n_x, x_max, x_min, x )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_lin_spacing
!
!===========================================================================

! number of points x
integer (wi), intent(in)  :: n_x

! maxmimum value of x
real    (wp), intent(in)  :: x_max

! minmimum value of x
real    (wp), intent(in)  :: x_min

! computed linearly spaced x values
real    (wp), intent(out) :: x(1:n_x)

!===========================================================================
! local variables for vera_lin_spacing
!===========================================================================

! counter to loop over the distribution
integer (wi)              :: ii

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*) , parameter :: RoutineName = 'VERA_LIN_SPACING'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for vera_lin_spacing
!===========================================================================

! compute the x values - these are linearly spaced
!                        across the range [x_min, x_max]
!
! start by initialising the x values to integer values, counting up from 0
do ii = vera_constants%one_i, n_x

  x(ii) = ii - vera_constants%one_i

end do

! scale the x values to be in the range [0,1]
x = x / (n_x - vera_constants%one_i)

! scale the x to be in the range [x_min, x_max]
x = ( x * (x_max - x_min) ) + x_min

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_lin_spacing

!=============================================================================
!
! vera_log_spacing -  Computes a logarithmically spaced points
!
! x that are logarithmically spaced across the range [x_min, x_max].
!
!=============================================================================

subroutine vera_log_spacing( n_x, x_max, x_min, x )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_log_spacing
!
!===========================================================================

! number of points x
integer (wi), intent(in)  :: n_x

! maxmimum value of x
real    (wp), intent(in)  :: x_max

! minmimum value of x
real    (wp), intent(in)  :: x_min

! computed logarithmically spaced x values
real    (wp), intent(out) :: x(1:n_x)

!===========================================================================
! local variables for vera_log_spacing
!===========================================================================

! counter to loop over the distribution
integer (wi)              :: ii

! scaling between successive x values, i.e. x(i+1) / x(i)
real    (wp)              :: x_scaling

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_LOG_SPACING'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for vera_log_spacing
!===========================================================================

! compute the x values - these are logarithmically spaced
!                        across the range [x_min, x_max]
!
! start by initialising the x values to integer values, counting up from 0
do ii = vera_constants%one_i, n_x

  x(ii) = ii - vera_constants%one_i

end do

! compute the scaling between successive x values, i.e. x(i+1) / x(i)
x_scaling = ( x_max / x_min )**( vera_constants%one            /               &
                                 ( n_x - vera_constants%one_i )  )

! compute the x values
x = ( x_scaling**x ) * x_min

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_log_spacing

!=============================================================================
!
! vera_log_normal        -  Computes the normal curve on
!                           logarithmically spaced values.
!
! The x are logarithmically spaced across the range [x_min, x_max] with
! the y values computed using a gaussian.
!
! This distribution is used to cast phantom aerosol size distributions.
!
!=============================================================================

subroutine vera_log_normal( x, mode, sigma, x_max, x_min, boundless, y )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_log_normal
!
!===========================================================================

! x values in the distribution
real (wp), intent(in)  :: x(:)

! mode of the normal distribution
real (wp), intent(in)  :: mode

! geometric standard deviation of the normal distribution,
! i.e. the width of the distribution
real (wp), intent(in)  :: sigma

! maxmimum value of x in the distribution
real (wp), intent(in)  :: x_max

! minmimum value of x in the distribution
real (wp), intent(in)  :: x_min

! value to compute outside of the range [x_min, x_max]
real (wp), intent(in)  :: boundless

! number concentrations in the distribution
! NOTE: this is scaled so that the total number concentration
!       is one particle per unit volume
real (wp), intent(out) :: y(1:size(x))


!===========================================================================
! local variables for vera_log_normal
!===========================================================================

! counter for a do ... if construction to replicate where
integer (wi)           :: ii

! normalisation factor for the y pdf values
real    (wp)           :: normal_sum

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_LOG_NORMAL'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for vera_log_normal
!===========================================================================

! initialise the number concentrations to 0
y = vera_constants%zero

! compute the number concentrations assuming a log-normal distribution
y = exp( -( ( log(x/mode) / log(sigma) )**vera_constants%two ) /               &
          ( vera_constants%two ) )

! compute normalisation factor for pdf values in the range x = [x_min, x_max]
normal_sum = vera_constants%zero
do ii = vera_constants%one_i, size( x )
  if ( (x(ii) >= x_min) .and. (x(ii) <= x_max) ) then
    normal_sum = normal_sum + y(ii)
  end if
end do

! normalise the computed pdf values
y = y / normal_sum

! deal with the pdf values outside of the the range x = [x_min, x_max]
do ii = vera_constants%one_i, size( x )
  if ( (x(ii) < x_min) .or. (x(ii) > x_max) ) then
    ! x is outside of the range x =[x_min, x_max], so set the computed
    ! pdf values to the constant specified by boundless
    y(ii) = boundless
  end if
end do

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_log_normal

!=============================================================================
!
! vera_log_normal_distn -  Computes a logarithmically
!                          spaced gaussian distribution.
!
! The x are logarithmically spaced across the range [x_min, x_max] with
! the y values computed using a gaussian.
!
! This distribution is used to cast phantom aerosol size distributions.
!
!=============================================================================

subroutine vera_log_normal_distn( n_x, mode, sigma, x_max, x_min, x, y )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_log_normal_distn
!
!===========================================================================

! number of particles x values in the distribution
integer (wi), intent(in)  :: n_x

! mode of the normal distribution
real    (wp), intent(in)  :: mode

! geometric standard deviation of the normal distribution,
! i.e. the width of the distribution
real    (wp), intent(in)  :: sigma

! maxmimum value of x in the distribution
real    (wp), intent(in)  :: x_max

! minmimum value of x in the distribution
real    (wp), intent(in)  :: x_min

! x values in the distribution
real    (wp), intent(out) :: x(1:n_x)

! number concentrations in the distribution
! NOTE: this is scaled so that the total number concentration
!       is one particle per unit volume
real    (wp), intent(out) :: y(1:n_x)


!===========================================================================
! local variables for vera_log_normal_distn
!===========================================================================
! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_LOG_NORMAL_DISTN'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for vera_log_normal_distn
!===========================================================================

! compute the x values - these are logarithmically spaced
!                        across the range [x_min, x_max]
!
call vera_log_spacing( n_x, x_max, x_min, x )

! compute the normal function
!
! adding a small value to x_max and taking the same small value from x_min
! ensures that the normal function is computed over the full range
! [x_min, x_max]
!
! the out of bounds value is set to vera_constants%zero
call vera_log_normal( x, mode, sigma, x_max+1.0_wp, x_min-1.0_wp,              &
                      vera_constants%zero, y )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_log_normal_distn

!=============================================================================
!
! vera_triangle        -  Computes triangular pdf.
!
! The x are spaced across the range [x_min, x_max] with
! the y values computed using a triagular function.
!
! This distribution is used to cast phantom hygroscopy distributions.
!
!=============================================================================

subroutine vera_triangle( x, x_peak, x_max, x_min, boundless, y )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_triangle
!
!===========================================================================

! x values in the distribution
real (wp), intent(in)  :: x(:)

! x value corresponding to the peak of the triangle
! NOTE: x_peak must be in the range ]x_min, x_max[
real (wp), intent(in)  :: x_peak

! maxmimum value of x in the distribution
real (wp), intent(in)  :: x_max

! minmimum value of x in the distribution
real (wp), intent(in)  :: x_min

! value to compute outside of the range [x_min, x_max]
real (wp), intent(in)  :: boundless

! pdf values for the distribution
! NOTE: pdf values are normalised so that the sum of the pdf values is 1
real (wp), intent(out) :: y(1:size(x))


!===========================================================================
! local variables for vera_triangle
!===========================================================================

! counter for a do ... if construction to replicate where
integer (wi)           :: ii

! gradient of left leg of triangle, i.e. for x = [x_min, x_peak]
real    (wp)           :: gradient_left

! gradient of right leg of triangle, i.e. for x = [x_peak, x_max]
real    (wp)           :: gradient_right

! normalisation factor for the y pdf values
real    (wp)           :: normal_sum

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_TRIANGLE'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for vera_triangle
!===========================================================================

! initialise the pdf to 0
y = vera_constants%zero

! compute the gradient of the left leg of the triangle
gradient_left  = vera_constants%two / ((x_max - x_min) * (x_peak - x_min))

! compute the gradient of the right leg of the triangle
gradient_right = vera_constants%two / ((x_max - x_min) * (x_peak - x_max))

! compute the pdf values for the x values in the range x = [x_min, x_peak]
! i.e. the left leg of the triangle
do ii = vera_constants%one_i, size( x )
  if ( (x(ii) >= x_min) .and. (x(ii) <= x_peak) ) then
    ! compute the pdf values
    y(ii) = gradient_left * ( x(ii) - x_min )
  end if
end do

! compute the pdf values for the x values in the range x = [x_peak, x_max]
! i.e. the right leg of the triangle
do ii = vera_constants%one_i, size( x )
  if ( (x(ii) > x_peak) .and. (x(ii) <= x_max) ) then
    ! compute the pdf values
    y(ii) = gradient_right * ( x(ii) - x_max )
  end if
end do

! compute normalisation factor for pdf values in the range x = [x_min, x_max]
normal_sum = vera_constants%zero
do ii = vera_constants%one_i, size( x )
  if ( (x(ii) >= x_min) .and. (x(ii) <= x_max) ) then
    normal_sum = normal_sum + y(ii)
  end if
end do

! normalise the computed pdf values
y = y / normal_sum

! deal with the pdf values outside of the the range x = [x_min, x_max]
do ii = vera_constants%one_i, size( x )
  if ( (x(ii) < x_min) .or. (x(ii) > x_max) ) then
    ! x is outside of the range x =[x_min, x_max], so set the computed
    ! pdf values to the constant specified by boundless
    y(ii) = boundless
  end if
end do

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_triangle

!=============================================================================
!
! vera_linear_flat_distn
!
! The flat distribution of x are linearly spaced across the
! range [x_min, x_max].
!
! This distribution can be used to cast a phantom aerosol
! hygroscopy distribution.
!
!=============================================================================

subroutine vera_linear_flat_distn( n_x, x_max, x_min, x, y )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_linear_flat_distn
!
!===========================================================================

! number of particles x values in the distribution
integer (wi), intent(in)  :: n_x

! maxmimum value of x in the distribution
real    (wp), intent(in)  :: x_max

! minmimum value of x in the distribution
real    (wp), intent(in)  :: x_min

! x values in the distribution
real    (wp), intent(out) :: x(1:n_x)

! number concentrations in the distribution
! NOTE: this is scaled so that the total number concentration
!       is one particle per unit volume
real    (wp), intent(out) :: y(1:n_x)

!===========================================================================
! local variables for vera_linear_flat_distn
!===========================================================================
! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_LINEAR_FLAT_DISTN'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start executable code for vera_linear_flat_distn
!===========================================================================

! compute the x values - these are linearly spaced
!                        across the range [x_min, x_max]
!
call vera_lin_spacing( n_x, x_max, x_min, x )

! compute the number concentrations
y = vera_constants%one / n_x

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_linear_flat_distn

end module vera_phantom_tools_mod
