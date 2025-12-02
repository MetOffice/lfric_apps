!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interpolate a time-varying profile to the current time
!> @details This is intended for idealised modelling scenarios in which
!>          forcing/relaxation is applied to prognostic fields. The
!>          forcing is specified through a namelist as a sequence of
!>          vertical profiles. The validity times of these profiles are
!>          also specified in the namelist. Linear interpolation is used
!>          to evaluate the forcing profile at the current model time.
module prof_temporal_interp_mod

  use constants_mod, only: i_def, r_def

  implicit none

  private

  public :: prof_temporal_interp

contains

  !> @brief Subroutine to interpolate vertical profile data in time
  !> @param[in] time_now Current model time (seconds)
  !> @param[in] n_heights Number of points in a vertical profile
  !> @param[in] n_times Number of vertical profiles
  !> @param[in] times Validity times of the profiles (seconds)
  !> @param[in] prof_vary Sequence of vertcal profiles
  !> @param[out] profile_now Profile data interpolated to current time
  subroutine prof_temporal_interp( time_now, n_heights, n_times, &
                                   times, prof_vary, profile_now )

  implicit none

  real( kind=r_def ),    intent(in)  :: time_now
  integer( kind=i_def ), intent(in)  :: n_heights, n_times
  real( kind=r_def ),    intent(in)  :: times( n_times )
  real( kind=r_def ),    intent(in)  :: prof_vary( n_heights * n_times )
  real( kind=r_def ),    intent(out) :: profile_now( n_heights )

  integer( kind=i_def )              :: time_index
  integer( kind=i_def )              :: k
  real( kind=r_def )                 :: time_effective
  real( kind=r_def )                 :: tau

  ! If only one vertical profile is given then temporal interpolation not needed
  if ( n_times == 1 ) then
    profile_now( 1 : n_heights ) = prof_vary( 1 : n_heights )
    return
  end if

  ! Treat all times earlier than the first validity time as though they
  ! were the first validity time. Similarly for times later than the last
  ! validity time. This gives "constant extrapolation in time".
  time_effective = min( max( time_now, times(1) ), times( n_times ) )

  ! Find the index of the validity time which is at or just before time_effective
  ! Note: it is required that times(1) <= time_effective <= times(n_times)
  do time_index = 1, n_times - 1
    if ( times( time_index + 1 ) > time_effective ) exit
  end do

  ! time_index adjusted if time_effective = times(n_times)
  time_index = min( n_times - 1, time_index )

  ! Linear weight for interpolation in time
  tau = ( time_effective - times( time_index ) ) / &
        ( times( time_index + 1 ) - times( time_index ) )

  ! Interpolate profile to current time
  do k = 1, n_heights
    profile_now( k ) = ( 1.0_r_def - tau )                               &
                       * prof_vary( (time_index - 1) * n_heights + k ) + &
                       tau * prof_vary( time_index * n_heights + k )
  end do

  end subroutine prof_temporal_interp
end module prof_temporal_interp_mod
