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
! This module contains subroutines to produce normally distributed, correlated
! random series. The application for this module is to produce synthetic
! noise to model the error in meteorological fields, either UM or observed.
!

module vera_noise_mod

use vera_kind_mod,   only: wp => vera_real, wi => vera_integer

! grab data from the Vera global module
use vera_global_mod, only: minus_one_i   => vera_minus_one_i   ,               &
                           zero_i        => vera_zero_i        ,               &
                           one_i         => vera_one_i         ,               &
                           two_i         => vera_two_i         ,               &
                           zero          => vera_zero          ,               &
                           one           => vera_one           ,               &
                           two           => vera_two           ,               &
                           one_hundred   => vera_hundred       ,               &
                           one_hundred_i => vera_one_hundred_i ,               &
                           minus_two     => vera_minus_two     ,               &
                           two_pi        => vera_two_pi        ,               &
                           vera_noise_control                  ,               &
                           vera_ns

implicit none

! Description:
!   This module implements a scheme within Vera to model the noise
!   on the meteorological inputs of (P, T, q, qcl, am). These noisy
!   inputs are then pushed through the Vera scheme, and the resulting
!   population of computed visual ranges is used to estimate the
!   error in the visibility.
!
! Method:
!   This module comprises eight subroutines:
!
!   vera_seed
!     Seeds the instrinsic random number generator random_number.
!
!   vera_sort
!     Sorts an input list of numbers into order, from lowest to highest.
!
!   vera_range
!     Computes estimates of the probabilities of visibile ranges
!     from a list of noisy visibilities.
!
!   vera_centile
!     Computes centiles from a list of noisy visibilities.
!
!   vera_random_normal_list
!     Produces a list of normally distributed random numbers,
!     with mean = 0 and standard deviation = 1.
!
!     This routine uses the Box-Muller transform to convert uniformly
!     distributed random numbers in the range [0,1] to a
!     normal distribution.
!
!   vera_normal_bounded
!     Produces a list of normally distributed random numbers, with
!     mean = 0 and standard deviation = 1, but bounded by upper
!     and lower bounds, i.e. the random numbers are in the range
!     [lower_bound, upper_bound] with lower_bound < upper_bound.
!
!   vera_cholesky_pivot
!     Computes the Cholesky decomposition of a square, positive-definite
!     matrix. This is used to factor the covariance matrix of the noise
!     on the input set (P, T, q, qcl, am).
!
!   vera_generate_noise
!     Generates noise to perturb an input set (P, T, q, qcl, am).
!     The noise is gaussian with standard deviation 1, and bounded
!     by settings in the vera global module:
!
!        vera_ns%upper_bound_normal - default is +1
!        vera_ns%lower_bound_normal - default is -1
!
!     If required, the generated perturbations can be correlated,
!     specified by an input covariance matrix.
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_NOISE_MOD'

private

! make public the optional output types and the subroutines
public :: vera_noise_mod_report              ,                                 &
          vera_random_normal_list            ,                                 &
          vera_normal_bounded                ,                                 &
          vera_cholesky_pivot                ,                                 &
          vera_generate_noise                ,                                 &
          vera_sort                          ,                                 &
          vera_seed                          ,                                 &
          vera_range                         ,                                 &
          vera_centile

! derived type for optional output from vera_normal_bounded
type :: vera_noise_mod_report
  real    (wp) :: lower_bound
  real    (wp) :: upper_bound
  integer (wi) :: n_wanted
  integer (wi) :: n_found
  integer (wi) :: loops
end type vera_noise_mod_report

contains

!=============================================================================
!
! vera_seed - sets the seed for the intrinsic random number generator
!             routine random_number
!
!             There are three possible methods to set the seed, and selection
!             is controlled by setting vera_noise_control%seed_method to:
!
!             vera_noise_control%seed_method_model
!               The model run date/time data and ensemble member number is used
!               to generate the seed. If the model run is not part of an
!               ensemble, and the Environment Variable ENS_MEMBER is not set,
!               then an ensemble member number of zero is used.
!
!               This method is reproducible between reruns of the same
!               model run.
!
!               This is the default seeding method.
!
!             vera_noise_control%seed_method_computer
!               The computer runtime date/time data is used to generate the
!               seed. Note that this method is not reproducible between reruns
!               of the same model run.
!
!             vera_noise_control%seed_method_constant
!               The constant value in vera_noise_control%seed_default is used
!               to generate the seed. This method is reproducible between
!               reruns of the same model run.
!
!=============================================================================

subroutine vera_seed( )

! pick up the model run date/time data
use model_time_mod,  only: i_second, i_minute, i_hour, i_day, i_month, i_year

! use the DrHook stuff
use parkind1,        only: jpim, jprb
use yomhook,         only: lhook, dr_hook

implicit none

!==========================================================================
! local variables for vera_seed
!
!==========================================================================

! date/time data, the eight values are:
! dt(1) - year
! dt(2) - month
! dt(3) - day
! dt(4) - shift from UTC
! dt(5) - hour
! dt(6) - minute
! dt(7) - seconds
! dt(8) - miiliseconds
integer (wi)              :: dt(8)

! ensemble member number
integer (wi)              :: ensemble_member

! number of integers used to seed the random number generator
integer (wi)              :: tam

! first guess at a single integer seed
integer (wi)              :: iarg

! maximum allowed value of the seed to prevent numerical overflow
integer (wi)   :: ihuge = huge(iarg)
integer (wi)   :: max_iarg

! first attempt seed
integer (wi), allocatable :: seed_first(:)

! second attempt seed
integer (wi), allocatable :: seed_second(:)

! random numbers for calculating second attempt seed
real    (wp), allocatable :: rnum(:)

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_SEED'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!==========================================================================
! start of executable code for vera_seed
!
!==========================================================================

! enquire the size of the seed, and reserve some memory for computing
! first and second attempt seeds
call random_seed( size = tam )
if ( .not. allocated( seed_first )  ) allocate( seed_first(tam)  )
if ( .not. allocated( seed_second ) ) allocate( seed_second(tam) )
if ( .not. allocated( rnum )        ) allocate( rnum(tam)        )

! check if random_number is to be seeded using date/time data, either from the
! model run or the runtime computer clock.
if ( (vera_noise_control%seed_method ==                                        &
      vera_noise_control%seed_method_model) .or.                               &
     (vera_noise_control%seed_method ==                                        &
      vera_noise_control%seed_method_computer) ) then

  ! are runtime computer clock date/time data required?
  if ( (vera_noise_control%seed_method ==                                      &
        vera_noise_control%seed_method_computer) ) then
    call date_and_time( values = dt )
  end if

  ! are model run date/time data required?
  if ( (vera_noise_control%seed_method ==                                      &
        vera_noise_control%seed_method_model) ) then

    ensemble_member = vera_noise_control%ensemble_member

    ! put the model run date/time data into the dt array
    dt(1) = i_year
    dt(2) = i_month
    dt(3) = i_day
    dt(4) = zero_i                          ! shift from UTC
    dt(5) = i_hour + one_i                  ! set hour range [1,24]
    dt(6) = i_minute
    dt(7) = ensemble_member + i_second      ! ensemble member into seconds
    dt(8) = ensemble_member + one_hundred_i ! ensemble member into milliseconds

  end if

  ! This computation for generating an integer first guess seed using
  ! date/time data has been pinched from subroutine stph_seed_gen_from_date in
  ! stph_seed.F90, part of the UM stochastic physics scheme.
  !
  ! Use full date in randomising seed to reduce chance of
  !  recycling from one run to the next.
  ! The formula calculates the days since ~2000AD and adds
  !  in time suitably inflated to fully change the seed.
  ! Only use last two digits of year to prevent numerical
  !  overflow at some date in the future.
  ! A random number generated from this seed is used to
  !  multiply the seed again
  dt(1) = dt(1) - 100_wi*int(0.01_wp*dt(1))
  iarg = int((dt(3) - 32075_wi +                                               &
          1461_wi*(dt(1) + 4800_wi + (dt(2) - 14_wi)/12_wi)/4_wi +             &
          367_wi*(dt(2) - 2_wi - (dt(2)-14_wi)/12_wi*12_wi)/12_wi -            &
          3_wi*((dt(1)+4900_wi+(dt(2)-14_wi)/12_wi)/100_wi)/4_wi)*1000_wi +    &
          dt(8)**2.86_wp + dt(7)**3.79_wp + dt(5)**5.12_wp +                   &
          dt(6)**3.24_wp)

  ! end of generating iarg using date/time data
else

  ! this is the fall-through default action, use the constant value defined by
  ! vera_noise_control%seed_default to seed random_number
  iarg = vera_noise_control%seed_default

end if

! constrain iarg in range 0 < iarg < huge(int)/2 to prevent numerical overflow
max_iarg = floor( one * ihuge / two )
iarg = mod( iarg, max_iarg )
iarg = max( iarg, zero_i )

! set the first attempt at a seed
seed_first(:) = iarg
call random_seed( put = seed_first(1:tam) )

! compute second attempt at a seed, as semi-random fractions of the
! first attempt seed
call random_number( rnum )
seed_second(:) = int(iarg * rnum(:))

! set final seed
call random_seed( put = seed_second(1:tam) )

! tidy up
deallocate( rnum        )
deallocate( seed_second )
deallocate( seed_first  )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_seed

!=============================================================================
!
! vera_sort - sorts a list of numbers from smallest to largest
!
!=============================================================================

subroutine vera_sort( x )

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!==========================================================================
! arguments for vera_sort
!
!==========================================================================

! the list of numbers to sort
real (wp), intent(in out) :: x(:)

!==========================================================================
! local variables for vera_sort
!
!==========================================================================

! loop counters
integer (wi) :: ii_loop
integer (wi) :: jj_loop

! temporary store to enable two elements in the array to be swapped
real    (wp) :: store

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_SORT'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!==========================================================================
! start of executable code for vera_sort
!
!==========================================================================

! sort the input array, putting the smallest element in position 1
do ii_loop = size(x) - one_i, one_i, minus_one_i

  do jj_loop = one_i, ii_loop

    if ( x( jj_loop ) > x( jj_loop + one_i ) ) then

      ! swap the two elements
      store                = x( jj_loop         )
      x( jj_loop         ) = x( jj_loop + one_i )
      x( jj_loop + one_i ) = store

    end if

    ! end of the jj_loop
  end do

  ! end of the ii_loop
end do

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_sort

!=============================================================================
!
! vera_range - computes estimates of the probabilities of visibile ranges
!              from a list of computed visibilities.
!
!=============================================================================

subroutine vera_range( visibilities, range_probability, range_in, sort )

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!==========================================================================
! arguments for vera_range
!
!==========================================================================
! the input array - a list of visibilities
real    (wp), intent(in)           :: visibilities(:)

! computed probabilities
real    (wp), intent(out)          :: range_probability(:)

! switch to indicate whether to sort the input list of visibilities,
! check this against vera_noise_control % sort_on
integer (wi), intent(in), optional :: sort

! visual range probabilities to compute
real    (wp), intent(in), optional :: range_in(:)

!==========================================================================
! local variables for vera_range
!
!==========================================================================

! sorted copy of the input array visibilities
real    (wp), allocatable :: vis_sorted(:)

! visibile ranges to use [m], e.g. [ 50.0, 200.0, 600.0, 1000.0, 4000.0 ]
real    (wp), allocatable :: ranges_use(:)

! fractional part of the array subscript
real    (wp)              :: frac_sub

! number of elements in the list of visibilities
real    (wp)              :: p_real

! mask for searching the list of visibilites for the required ranges
integer (wi), allocatable :: visibilities_mask(:)

! number of elements in the list of visibilities
integer (wi)              :: p_int

! loop counter
integer (wi)              :: ii_loop

! array subscript to use
integer (wi)              :: ii

! index counter for a do ... if construction replicating where
integer (wi)              :: where_loop

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_RANGE'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!==========================================================================
! start of executable code for vera_range
!
!==========================================================================

! compute how many visibilites are in the input list
p_int  = size(visibilities)
p_real = size(visibilities, kind = wp)

! check to see if a list of ranges to use has been passed to this routine
if (present(range_in)) then

  ! use the input list of ranges
  allocate ( ranges_use(size(range_in)) )
  ranges_use = range_in

else

  ! use the default list of visible ranges in the global module
  allocate ( ranges_use(vera_noise_control%n_ranges) )
  ranges_use = vera_noise_control%ranges(one_i:                                &
                                         vera_noise_control%n_ranges)

end if

! form the sorted list of input visibilities
allocate( vis_sorted(p_int) )
vis_sorted = visibilities

! check to see if the input list of visibilities needs sorting
if (present(sort)) then

  if ( sort == vera_noise_control % sort_on ) then
    ! sort switch is set to on, so sort the list
    call vera_sort( vis_sorted )
  end if

else

  ! switch not used on input, so sort the list
  call vera_sort( vis_sorted )

end if

! create the mask used to search the visibility population for vis ranges
allocate( visibilities_mask(p_int) )

! loop over the visual ranges for which to compute probabilities
do ii_loop = one_i, size(ranges_use)

  ! reset the mask
  visibilities_mask = zero_i

  ! flag which of the computed visibilities are less than the range
  ! indicated by ranges_use(ii_loop)
  do where_loop = one_i, p_int
    if ( vis_sorted(where_loop) <=  ranges_use(ii_loop) ) then
      visibilities_mask(where_loop) = one_i
    end if
  end do

  ! find the array element in the sorted list where the computed
  ! visibility is just less than the range required
  ii = sum(visibilities_mask)

  ! interpolate between the visibilities to estimate the probability
  ! of the visibility not exceeding the required threshold
  if ( ii < one_i ) then
    range_probability(ii_loop) = zero
  end if

  if ( ii >= p_int ) then
    range_probability(ii_loop) = one
  end if

  if ( ( ii >= one_i ) .and. ( ii < p_int ) ) then
    frac_sub = log( ranges_use(ii_loop   ) / vis_sorted(ii) ) /                &
               log( vis_sorted(ii + one_i) / vis_sorted(ii) )
    range_probability(ii_loop) = ( ii + frac_sub ) / p_real
  end if

  ! end of looping over the ranges
end do

! tidy up
if ( allocated( visibilities_mask ) ) deallocate( visibilities_mask )
if ( allocated( ranges_use        ) ) deallocate( ranges_use        )
if ( allocated( vis_sorted        ) ) deallocate( vis_sorted        )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_range

!=============================================================================
!
! vera_centile - computes centiles from a list of numbers.
!
!=============================================================================

subroutine vera_centile( x, centile_value, centile, sort )

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!==========================================================================
! arguments for vera_centile
!
!==========================================================================
! the input array - a list of numbers
real (wp), intent(in)           :: x(:)

! computed centile
real (wp), intent(out)          :: centile_value(:)

! switch to indicate whether to sort the input list of visibilities,
! check this against vera_noise_control % sort_on
integer (wi), intent(in), optional :: sort

! centiles to compute
real (wp), intent(in), optional :: centile(:)

!==========================================================================
! local variables for vera_centile
!
!==========================================================================

! sorted copy of the input array x
real    (wp), allocatable :: x_sorted(:)

! centiles to use, e.g. [ 0.0, 25.0, 50.0, 75.0, 100.0 ]
real    (wp), allocatable :: centiles_use(:)

! fractional array subscript to use
real    (wp)              :: ii_fractional

! fractional part of the array subscript
real    (wp)              :: frac_sub

! number of elements in the input array x
real    (wp)              :: p_real

! number of elements in the input array x
integer (wi)              :: p_int

! loop counter
integer (wi)              :: ii_loop

! array subscript to use
integer (wi)              :: ii

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_CENTILE'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!==========================================================================
! start of executable code for vera_centile
!
!==========================================================================

! how many elements in the input array?
p_int  = size(x)
p_real = size(x, kind = wp)

! check to see if a list of centiles to use has been input
if (present(centile)) then

  ! use the input list of centiles
  allocate ( centiles_use(size(centile)) )
  centiles_use = centile

else

  ! use the default list of centiles in the global module
  allocate ( centiles_use(vera_noise_control%n_centiles) )
  centiles_use = vera_noise_control%centiles(one_i:                            &
                                      vera_noise_control%n_centiles)

end if

! form the sorted list of input values
allocate( x_sorted(p_int) )
x_sorted = x

! check to see if the input list needs sorting
if (present(sort)) then

  if ( sort == vera_noise_control % sort_on ) then
    ! sort switch is set to on, so sort the list
    call vera_sort( x_sorted )
  end if

else

  ! switch not used on input, so sort the list
  call vera_sort( x_sorted )

end if

! convert the centiles required into ratios, i.e. 75.0% to 0.75
centiles_use = centiles_use / one_hundred

! loop over the centiles to compute
do ii_loop = one_i, size(centiles_use)

  ! compute the subscript to use
  ii_fractional = ( ( centiles_use(ii_loop) *                                  &
                    ( p_real - one ) )  )   + one

  ii = int( ii_fractional )

  ! compute the interpolated value of x
  if ( ii_fractional <= one ) then
    centile_value(ii_loop) = x_sorted(one_i)
  end if

  if ( ii_fractional >= p_real ) then
    centile_value(ii_loop) = x_sorted(p_int)
  end if

  if ( ( ii_fractional > one ) .and. ( ii_fractional < p_real ) ) then
    frac_sub = ii_fractional - ii
    centile_value(ii_loop) =                                                   &
                     exp( ( log(x_sorted(ii      )) * (one-frac_sub) ) +       &
                          ( log(x_sorted(ii+one_i)) * (frac_sub    ) ) )
  end if

  ! end of looping over the centiles
end do

! tidy up
if ( allocated( x_sorted     ) ) deallocate( x_sorted     )
if ( allocated( centiles_use ) ) deallocate( centiles_use )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_centile

!=============================================================================
!
! vera_random_normal_list - produces a list of normally distributed random
!                           numbers, with mean = 0 and standard deviation = 1.
!
!                           This routine uses the Box-Muller transform to
!                           convert uniformly distributed random numbers in
!                           the range [0,1] to a normal distribution.
!
!=============================================================================

subroutine vera_random_normal_list( numbers )

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!==========================================================================
! arguments for vera_random_normal_list
!
!==========================================================================

! the output series of normally distributed random numbers
real    (wp), intent(in out) :: numbers(:)

!==========================================================================
! local variables for vera_random_normal_list
!
!==========================================================================

! how many random numbers to produce in the series?
integer (wi)                 :: n_numbers

! series of uniformly distributed random numbers in the range [0,1]
real    (wp)                 :: u_1(1: size(numbers) / two_i)

! another series of uniformly distributed random numbers in the range [0,1]
real    (wp)                 :: u_2(1: size(numbers) / two_i)

! series of normally distributed random numbers
real    (wp)                 :: z_1(1: size(numbers) / two_i)

! another series of normally distributed random numbers
real    (wp)                 :: z_2(1: size(numbers) / two_i)

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_RANDOM_NORMAL_LIST'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!==========================================================================
! start of executable code for vera_random_normal_list
!
!==========================================================================

! how many random numbers are required?
n_numbers = size(numbers)

! if not already done so, set the seed for the random number generator
if (vera_noise_control%seed_flag==vera_noise_control%seed_flag_reset) then
  ! construct and plant the seed
  call vera_seed( )

  ! update the seed planting flag
  vera_noise_control%seed_flag = vera_noise_control%seed_flag_leave
end if

! grab a series of uniformly distributed random numbers in the range [0,1]
call random_number( u_1 )

! grab a series of uniformly distributed random numbers in the range [0,1]
call random_number( u_2 )

! form the pairs of values of normally distributed random numbers
z_1 = sqrt( minus_two * log(u_1) ) * cos( two_pi * u_2 )
z_2 = sqrt( minus_two * log(u_1) ) * sin( two_pi * u_2 )

! pop the normally distubuted pairs into the output series
numbers( one_i                        : n_numbers / two_i ) = z_1
numbers( ( n_numbers / two_i ) + one_i: n_numbers         ) = z_2

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_random_normal_list

!=============================================================================
!
! vera_normal_bounded - produces a list of normally distributed random
!                       numbers, with mean = 0 and standard deviation = 1,
!                       but bounded by upper and lower bounds, i.e. the random
!                       numbers are in the range
!
!                       [lower_bound, upper_bound]
!
!                       with lower_bound < upper_bound.
!
!
!=============================================================================

subroutine vera_normal_bounded( numbers, upper_bound, lower_bound, report )

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!==========================================================================
! arguments for vera_normal_bounded
!
!==========================================================================

! the output series of normally distributed random numbers
real (wp), intent(in out)        :: numbers(:)

! upper bound of the random numbers
real (wp), intent(in) , optional :: upper_bound

! lower bound of the random numbers
real (wp), intent(in) , optional :: lower_bound

! optional reporting from this routine
type(vera_noise_mod_report), intent(out), optional :: report

!==========================================================================
! local variables  for vera_normal_bounded
!
!==========================================================================

! upper bound of the random numbers to use, i.e. default or input value?
real    (wp) :: use_upper_bound

! lower bound of the random numbers to use, i.e. default or input value?
real    (wp) :: use_lower_bound

! how many random numbers to produce in the series?
integer (wi) :: n_numbers

! how many loops have been done to try and fill the randon number list?
integer (wi) :: n_loops = zero_i

! how many random numbers in the trial list are in the required range?
integer (wi) :: n_in_range

! index for the list of computed pseudo-random numbers
integer (wi) :: ii_upper

! index for the list of computed pseudo-random numbers
integer (wi) :: ii_lower

! index counter for a do ... if construction replicating where
integer (wi) :: where_loop

! numberof pseudo-random numbers to grab from the
! normal random number generator
integer (wi) :: n_grab

! trial series of normally distributed random numbers
real    (wp), allocatable :: trial_numbers(:)

! counter for the trial numbers
integer (wi), allocatable :: trial_counter(:)

! mask for the trial numbers
logical     , allocatable :: trial_mask(:)

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_NORMAL_BOUNDED'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!==========================================================================
! start of executable code for vera_normal_bounded
!
!==========================================================================

! check to see if a upper bound has been passed to this routine
use_upper_bound = vera_ns%upper_bound_normal
if ( present( upper_bound ) ) use_upper_bound = upper_bound

! check to see if a lower bound has been passed to this routine
use_lower_bound = vera_ns%lower_bound_normal
if ( present( lower_bound ) ) use_lower_bound = lower_bound

! how many random numbers are required?
n_numbers = size(numbers)

! set the lower index of the number list to the first entry
ii_lower  = one_i

! set the loop counter to zero
n_loops   = zero_i

! keep filling up the random number list until it's full,
! or have tried max_loops number of loops
do while ( ( ii_lower <= n_numbers )             .and.                         &
           ( n_loops < vera_noise_control%max_loops ) )

  ! update the loop counter
  n_loops = n_loops + one_i

  ! grab some random numbers
  n_grab = n_numbers - ii_lower + one_i

  ! ensure n_grab is even, rounding up if necessary
  if ( mod( n_grab, two_i )  /= zero_i ) n_grab = n_grab + one_i

  ! allocate some memory for the trial random numbers
  allocate( trial_numbers( n_grab) )
  allocate( trial_mask(    n_grab) )
  allocate( trial_counter( n_grab) )

  ! grab some random numbers
  call vera_random_normal_list( trial_numbers )

  ! construct the mask to find trial numbers in the required range
  trial_mask = ( ( trial_numbers >= use_lower_bound ) .and.                    &
                 ( trial_numbers <= use_upper_bound )       )

  ! set the trial counter to zero
  trial_counter = zero_i

  ! find the entries in the trial list that fit the required range,
  ! i.e. trial_numbers = [lower_bound, upper_bound]
  do where_loop = one_i, n_grab
    if ( trial_mask(where_loop) .eqv. .true. ) then
      trial_counter(where_loop) = one_i
    end if
  end do

  ! compute how many trial numbers are in the required range
  n_in_range = sum( trial_counter )

  ! check if there are suitable candidates in the trial list to pop
  ! into the random number list
  if ( n_in_range > zero_i ) then

    ! update the upper index of the numbers list
    ii_upper = ii_lower + n_in_range - one_i

    ! cap the upper index of the numbers list
    ! (the upper index and number of .true. elements in
    !  the trial mask cannot exceed n_numbers so that when
    !  passed to pack below the result will fit into the
    !  numbers list without going out of bounds)
    if ( ii_upper > n_numbers ) then
      ! move backwards through the list from the end
      wloop: do where_loop = n_grab, one_i, minus_one_i
        ! eliminate the first element found in the mask
        if ( trial_mask(where_loop) .eqv. .true. ) then
          trial_mask(where_loop) = .false.
          ! ... until only n_numbers elements remain
          ii_upper = ii_upper - one_i
          if ( ii_upper == n_numbers ) exit wloop
        end if
      end do wloop
    end if

    ! pop random numbers into the numbers list
    numbers( ii_lower : ii_upper ) = pack(trial_numbers,  mask = trial_mask)

    ! update the lower index of the numbers list
    ii_lower = ii_upper + one_i

    ! end of checking for suitable trial numbers
  end if

  if ( allocated( trial_counter ) ) deallocate( trial_counter )
  if ( allocated( trial_mask    ) ) deallocate( trial_mask    )
  if ( allocated( trial_numbers ) ) deallocate( trial_numbers )

  ! end of filling up the random numbers list
end do

! deallocate the memory used by the trial numbers list
if ( allocated( trial_counter ) ) deallocate( trial_counter )
if ( allocated( trial_mask    ) ) deallocate( trial_mask    )
if ( allocated( trial_numbers ) ) deallocate( trial_numbers )

! if required, fill in the optional output, out
if ( present( report ) ) report =                                              &
  vera_noise_mod_report( use_lower_bound  ,                                    &
                         use_upper_bound  ,                                    &
                         n_numbers        ,                                    &
                         ii_lower - one_i ,                                    &
                         n_loops            )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_normal_bounded

!=============================================================================
!
! vera_cholesky_pivot - computes the Cholesky factorisation of the covariance
!                       matrix for the meteological inputs to Vera, i.e.
!                       (P, T, q, qcl, am).
!
!                       Complete pivoting is implemented so that the range of
!                       values for any off-diagnonal entry in the covariance
!                       matrix is Aij = [-1,1] for i /= j.
!
!                       The diagonal entries in the covariance matrix must all
!                       be set to 1, i.e. Aii = 1
!
!                       The input covariance matrix must be upper triangular
!                       in form.
!
!                       The computed factorisation matrixm RL is lower
!                       triangular in form. This means that this matrix can be
!                       used to compute correlated noise, noise_cor as
!
!                       noise_cor = matmul( rl, noise )
!
!                       where noise is a (5, m) array of uncorrelated Guassian
!                       noise.
!
!=============================================================================

subroutine vera_cholesky_pivot( au, rl )

use vera_kind_mod,   only: wp => vera_real, wi => vera_integer

use vera_global_mod, only: one_i => vera_one_i

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!==========================================================================
! arguments for vera_cholesky_pivot
!
!==========================================================================

! input upper triangular covariance matrix
real     (wp), intent(in)      :: au(:, :)

! output lower triangular Cholesky factorisation of the covariance matrix
real     (wp), intent(in out)  :: rl(:, :)

!==========================================================================
! local variables  for vera_cholesky_pivot
!
!==========================================================================

! permutation matrix used to compute Cholesky factorisation
real     (wp), allocatable     :: permutation(:, :)

! index for generating the permutation matrix
integer  (wi), allocatable     :: p(:)

! vector used to swap columns and rows
real     (wp), allocatable     :: swap(:)

! scalar used to swap entries in the index for the permutation matrix
integer  (wi)                  :: swap_p

! rank of the input upper triangular covariance matrix
integer  (wi)                  :: n

! index of the next pivot point
integer  (wi)                  :: pivot

! loop counters
integer  (wi)                  :: ii, jj, kk

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName='VERA_CHOLESKY_PIVOT'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!==========================================================================
! start of executable code for vera_cholesky_pivot
!
!==========================================================================

! compute the rank of the input matrix AU
n = size(au, dim = 1)

! copy the input covariance matrix into the output factorisation matrix
rl = au

! create the allocatable local arrays
allocate( permutation(n, n) )
allocate(           p(n)    )
allocate(        swap(n)    )

! form the index for generating the permutation matrix
do ii = 1_wi, n
  p(ii) = ii
end do

! loop over the pivots in the matrix to compute the factorisation matrix
! NOTE: the computed factorisation from this loop is upper triangular
do kk = one_i, n

  ! find the largest pivot point
  pivot = kk
  do ii = kk + one_i, n
    if ( rl(ii, ii) > rl(ii - one_i, ii - one_i) ) pivot = ii
  end do

  ! effect the pivot
  if ( pivot > kk ) then
    ! swap the columns
    swap                   = rl( one_i : n, pivot )
    rl( one_i : n, pivot ) = rl( one_i : n, kk    )
    rl( one_i : n, kk    ) = swap
    ! swap the rows
    swap                   = rl( pivot, one_i : n )
    rl( pivot, one_i : n ) = rl( kk   , one_i : n )
    rl( kk   , one_i : n ) = swap
    ! swap the permutation index entries
    swap_p   = p(pivot)
    p(pivot) = p(kk)
    p(kk)    = swap_p
  end if

  ! compute the pivot point value
  rl(kk, kk) = sqrt( rl(kk, kk) )

  ! compute the submatrix entries
  do jj = kk + one_i, n
    do ii = kk + one_i, jj
      rl(ii, jj) = rl(ii, jj) - ( rl(kk, ii) * rl(kk, jj) )
    end do
  end do

  ! end of looping over the pivots
end do

! compute the permutation matrix
permutation = 0.0_wp
do ii = one_i, n
  permutation(p(ii), ii) = 1.0_wp
end do

! use the permutation matrix to unwind the pivoting on the factorisation
! matrix Rl, i.e. compute P.RL.TrP
rl = matmul( matmul( permutation, rl), transpose(permutation) )

! compute the matrix Tr(P.RL.TrP) - this converts the factorisation
! matrix from an upper to a lower triangular form
rl = transpose(rl)

! tidy up
deallocate( swap        )
deallocate( p           )
deallocate( permutation )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_cholesky_pivot

!=============================================================================
!
! vera_generate_noise - generates noise to perturb an input set
!                       (P, T, q, qcl, am).
!
!                       The noise is gaussian with standard deviation 1,
!                       and bounded by settings in the global module
!                       vera_global_mod
!
!                          upper_bound_normal - default is +1
!                          lower_bound_normal - default is -1
!
!                       If required, the generated perturbations can be
!                       correlated, specified by an input covariance matrix.
!
!=============================================================================

subroutine vera_generate_noise( noise, covariance )

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!==========================================================================
! arguments for vera_generate_noise
!
!==========================================================================

! the noise perturbations
real    (wp), intent(in out)       :: noise(:, :)

! covariance matrix for the noise
real    (wp), intent(in), optional :: covariance(:, :)

!==========================================================================
! local variables  for vera_generate_noise
!
!==========================================================================


! how many sets of noise to produce
integer (wi)                       :: n_noise

! how many variables in a set
integer (wi)                       :: n_variables

! Cholesky factorisation of the covariance matrix for the noise
real    (wp), allocatable          :: cholesky_factor(:, :)

! uncorrelated, random, gaussian perturbations
real    (wp), allocatable          :: gauss(:)

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_GENERATE_NOISE'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!==========================================================================
! start of executable code for vera_generate_noise
!
!==========================================================================

! compute how many sets of noise are required
n_noise     = size( noise, dim = 2 )

! compute how many variables in a set, e.g 5 for the set (P, T, q, qcl, am)
n_variables = size( noise, dim = 1 )

! generate random, uncorrelated, gaussian perturbations
allocate( gauss( n_noise * n_variables ) )
call vera_normal_bounded( gauss )
noise = reshape( gauss, ( [n_variables, n_noise] ) )
! tidy up
if ( allocated( gauss ) ) deallocate( gauss )

! check to see if a covariance matrix has been passed to this routine,
! if so then use this to form correlated gaussian perturbations
if ( present( covariance ) ) then
  ! use the input covariance matrix

  ! compute the Cholesky decomposition of the covariance matrix
  allocate( cholesky_factor( n_variables, n_variables) )
  call vera_cholesky_pivot( covariance, cholesky_factor )

  ! compute the noise series as the matrix multiplication operation
  ! noise = cholesky_factor # gaussian noise
  noise = matmul( cholesky_factor, noise )

  ! tidy up
  if ( allocated( cholesky_factor ) ) deallocate( cholesky_factor )

end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_generate_noise


end module vera_noise_mod
