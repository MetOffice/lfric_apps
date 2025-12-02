! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to update the prognostic precipitation fraction.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Before calling this routine,
! we store the increments to qrain (and optionally including
! qgraup as well) attributable to each of the sub-grid regions
! containing precipitation;
! Overlap of the following with the precip fraction at start-of-step:
!  - liquid-only cloud
!  - mixed-phase cloud
!  - ice-only cloud
!  - clear-sky,
! A 5th category includes the area of precip newly created during the
! current step.
! The precip mass increments from each of these 5 regions are stored
! in the dqprec arrays, which are passed in.
! We now calculate the new in-region precip mixing-ratio
! in each region.  Based on these, we compute a new effective precip
! fraction such that the in-region mean value of qp is the
! precip-mass-weighted mean over the regions, <qp*qp> / <qp>
! (where < > denotes grid-box mean).
!
! Therefore  <qp> / sigma = <qp*qp> / <qp>
!
! => sigma = <qp>^2 / <qp^2>

module lsp_update_precfrac_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_UPDATE_PRECFRAC_MOD'

contains

subroutine lsp_update_precfrac( points,                                        &
                                rain_liq, rain_mix, rain_ice,                  &
                                rain_clear, rain_new,                          &
                                dqprec_liq, dqprec_mix, dqprec_ice,            &
                                dqprec_clear, dqprec_new,                      &
                                qrain, qgraup, precfrac_k )

use um_types,             only: real_lsprec
use lsprec_mod,           only: zero, one, small_number
use mphys_inputs_mod,     only: l_subgrid_graupel_frac

! Dr Hook modules
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

implicit none

! Number of points
integer, intent(in) :: points

! Area fractions of each of the 5 precip regions
real (kind=real_lsprec), intent(in) :: rain_liq(points)
real (kind=real_lsprec), intent(in) :: rain_mix(points)
real (kind=real_lsprec), intent(in) :: rain_ice(points)
real (kind=real_lsprec), intent(in) :: rain_clear(points)
real (kind=real_lsprec), intent(in) :: rain_new(points)

! Contribution of each region to the grid-mean precip mass increment
real (kind=real_lsprec), intent(in) :: dqprec_liq(points)
real (kind=real_lsprec), intent(in) :: dqprec_mix(points)
real (kind=real_lsprec), intent(in) :: dqprec_ice(points)
real (kind=real_lsprec), intent(in) :: dqprec_clear(points)
real (kind=real_lsprec), intent(in) :: dqprec_new(points)

! Grid-mean rain and graupel mass after the above increments have been added on
real (kind=real_lsprec), intent(in) :: qrain(points)
real (kind=real_lsprec), intent(in) :: qgraup(points)

! Effective sub-grid fraction of precip to be computed
real (kind=real_lsprec), intent(in out) :: precfrac_k(points)


! Grid-mean precip mass
real (kind=real_lsprec) :: qprec(points)

! In-region precip mass at start-of-step
real (kind=real_lsprec) :: qprec_0

! New in-region precip masses in each of the 5 sub-regions
real (kind=real_lsprec) :: qprec_liq
real (kind=real_lsprec) :: qprec_mix
real (kind=real_lsprec) :: qprec_ice
real (kind=real_lsprec) :: qprec_clear
real (kind=real_lsprec) :: qprec_new

! Grid-mean of qprec and qprec^2
real (kind=real_lsprec) :: qprec_bar
real (kind=real_lsprec) :: qprec_sq_bar

! Loop counters
integer :: i

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_UPDATE_PRECFRAC'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Compute grid-mean total precip mass
if ( l_subgrid_graupel_frac ) then
  ! "Precip" mass includes graupel as well as rain
  do i = 1, points
    qprec(i) = qrain(i) + qgraup(i)
  end do
else
  ! "Precip" only includes rain
  do i = 1, points
    qprec(i) = qrain(i)
  end do
end if

! Find points where precip-mass non-zero
do i = 1, points
  if ( qprec(i) > zero ) then

    ! Calculate in-region rain-mass before the stored increments this
    ! step were added on.
    ! This is latest grid-mean qprec, minus all the increments,
    ! divided by the original rain-fraction (not including the "new" bit)
    qprec_0 = ( qprec(i) - ( dqprec_liq(i) + dqprec_mix(i) + dqprec_ice(i)     &
                           + dqprec_clear(i) + dqprec_new(i) ) )               &
            / max( rain_liq(i) + rain_mix(i) + rain_ice(i) + rain_clear(i),    &
                   small_number )

    ! Compute latest precip-mass inside each region
    qprec_liq = qprec_0 + dqprec_liq(i) / max( rain_liq(i), small_number )
    qprec_mix = qprec_0 + dqprec_mix(i) / max( rain_mix(i), small_number )
    qprec_ice = qprec_0 + dqprec_ice(i) / max( rain_ice(i), small_number )
    qprec_clear = qprec_0 + dqprec_clear(i) / max( rain_clear(i), small_number )
    ! "New" region doesn't have a contribution from the pre-existing precip
    qprec_new = dqprec_new(i) / max( rain_new(i), small_number )

    ! Check to remove negative values
    qprec_liq = max( qprec_liq, zero )
    qprec_mix = max( qprec_mix, zero )
    qprec_ice = max( qprec_ice, zero )
    qprec_clear = max( qprec_clear, zero )
    qprec_new = max( qprec_new, zero )

    ! Compute grid-mean of qprec over all 5 regions
    ! (not always equal to qprec, due to the check to remove negatives)
    qprec_bar = rain_liq(i) * qprec_liq                                        &
              + rain_mix(i) * qprec_mix                                        &
              + rain_ice(i) * qprec_ice                                        &
              + rain_clear(i) * qprec_clear                                    &
              + rain_new(i) * qprec_new

    ! Occasionally, when qprec is very small and a much larger amount
    ! has just been removed by the dqprec increments, qprec_bar
    ! can come out too small for the following division to be numerically
    ! safe.  Check here that qprec_bar > sqrt(tiny).  If not, just
    ! leave precfrac_k unchanged.
    if ( qprec_bar > small_number ) then

      ! Compute grid-mean of qprec^2
      qprec_sq_bar = rain_liq(i) * qprec_liq * qprec_liq                       &
                   + rain_mix(i) * qprec_mix * qprec_mix                       &
                   + rain_ice(i) * qprec_ice * qprec_ice                       &
                   + rain_clear(i) * qprec_clear * qprec_clear                 &
                   + rain_new(i) * qprec_new * qprec_new

      ! Compute effective fraction = <qprec>^2 / <qprec^2>
      precfrac_k(i) = qprec_bar * qprec_bar / qprec_sq_bar

      ! Rounding errors can occasionally yield values very slightly > 1;
      ! reset these to 1
      precfrac_k(i) = min( precfrac_k(i), one )

    end if
  end if

end do


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_update_precfrac
end module lsp_update_precfrac_mod
