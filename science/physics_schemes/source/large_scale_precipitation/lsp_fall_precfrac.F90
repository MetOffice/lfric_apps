! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Sedimentation of the prognostic sub-grid fraction of precipitation
! (rain and optionally graupel as well).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

module lsp_fall_precfrac_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_FALL_PRECFRAC_MOD'

contains

subroutine lsp_fall_precfrac( points, qrain, qgraup, rainrate, grauprate,      &
                              dhi, rhor, precfrac_k, precfrac_fall )

use mphys_inputs_mod, only: l_subgrid_graupel_frac
use um_types,         only: real_lsprec

use lsp_combine_precfrac_mod, only: lsp_combine_precfrac

! Dr Hook Modules
use yomhook,          only: lhook, dr_hook
use parkind1,         only: jprb, jpim

implicit none

! Number of points
integer, intent(in) :: points

! Pre-existing mass of rain and graupel on level k
real(kind=real_lsprec), intent(in) :: qrain(points)
real(kind=real_lsprec), intent(in) :: qgraup(points)

! Fall fluxes of rain and graupel from above
real(kind=real_lsprec), intent(in) :: rainrate(points)
real(kind=real_lsprec), intent(in) :: grauprate(points)

! Timestep / thickness of model layer / s m-1
real(kind=real_lsprec), intent(in) :: dhi(points)

! 1 / Air density / m3 kg-1
real(kind=real_lsprec), intent(in) :: rhor(points)

! Value of prognostic precipitation field at level k;
! updated in this routine with increment due to fall-in from above
real(kind=real_lsprec), intent(in out) :: precfrac_k(points)

! Fractional area of falling precip:
! in:  area of precip falling in from above
! out: area of precip falling down to next level below
real(kind=real_lsprec), intent(in out) :: precfrac_fall(points)

! Mass of precipitation pre-existing at level k
real(kind=real_lsprec) :: qprec_k(points)
! Mass of precipitation falling in from above
real(kind=real_lsprec) :: qprec_fall(points)

! Loop counter
integer :: i


! Stuff for DrHook
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_FALL_PRECFRAC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


if ( l_subgrid_graupel_frac ) then
  ! If including graupel in the precipitation fraction...
  do i = 1, points
    ! Set mass of precip pre-existing at level k
    qprec_k(i) = qrain(i) + qgraup(i)
    ! Compute mass of precipitation falling into level k from above
    qprec_fall(i) = ( rainrate(i) + grauprate(i) ) * dhi(i) * rhor(i)
  end do
else
  ! If the precipitation fraction is only for rain, not graupel...
  do i = 1, points
    ! Set mass of precip pre-existing at level k
    qprec_k(i) = qrain(i)
    ! Compute mass of precipitation falling into level k from above
    qprec_fall(i) = rainrate(i) * dhi(i) * rhor(i)
  end do
end if

! Routine computes combined effective fraction assuming maximal overlap
! between the falling and pre-existing precip fractions
! (precfrac_k is updated)
call lsp_combine_precfrac( points,                                             &
                           qprec_k, qprec_fall, precfrac_k, precfrac_fall )

do i = 1, points
  ! Set falling precip fraction to be used at the next level below
  precfrac_fall(i) = precfrac_k(i)
end do


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_fall_precfrac


end module lsp_fall_precfrac_mod
