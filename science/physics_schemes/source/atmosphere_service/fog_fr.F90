! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculates fraction of grid-box covered with fog.

! Description:
!   Calculates the fraction of a gridsquare with visibility
!   less than threshold, Vis_thresh, given the total water
!   mixing ratio (qt), temperature (T), pressure (p), and the
!   (Murk) aerosol mass mixing ratio (m), assuming a triangular
!   distribution of states about the median, characterised by
!   a critical relative humdity value, RHcrit.
!   NB:  Throughout, levels are counted from the bottom up,
!   i.e. the lowest level under consideration is level 1, the
!   next lowest level 2, and so on.

!   Suitable for single-column use.

! Documentation:
!    Wright, B. J., 1997: Improvements to the Nimrod Visibility
!       Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!    Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!       for Nimrod. Met. Office FR Tech Rep., No. 222.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_fog

! Code description:
!   This code is written to UMDP3 standards.

! Subroutine Interface:
module fog_fr_mod

implicit none

character(len=*), parameter, private :: ModuleName = 'FOG_FR_MOD'
contains

subroutine fog_fr(                                                             &
 p_layer,rhcrit,levels,pfield,                                                 &
 t,aerosol,l_murk,q,qcl,vis,ff,nvis)

! Modules
use qsat_mod, only: qsat_wat

use vistoqt_mod, only: vistoqt

use um_types, only: real_umphys

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) ::                                                         &
 levels,                                                                       &
                         ! in No. of levels being processed.
 pfield,                                                                       &
                         ! in No. of points in field (at one level).
 nvis                ! in No. of visibility thresholds
real(kind=real_umphys), intent(in) ::                                          &
 p_layer(pfield,levels),                                                       &
                            ! in pressure (Pa) at processed levels.
 rhcrit(levels),                                                               &
                         ! in Critical relative humidity.  See the
!                          !    the paragraph incorporating eqs P292.11
!                          !    to P292.14; the values need to be tuned
!                          !    for the given set of levels.
   q(pfield,levels),                                                           &
                           ! in Specific Humidity
!                          !    (kg per kg air).
   qcl(pfield,levels),                                                         &
                           ! Cloud liquid water content at
!                          !     processed levels (kg per kg air).
   t(pfield,levels),                                                           &
                           ! in Temperature (K).
   aerosol(pfield,levels),                                                     &
                              ! in Aerosol mixing ratio(ug/kg)
   vis(pfield,levels,nvis)  ! Visibility thresholds
logical, intent(in) ::                                                         &
   l_murk               ! in : Aerosol present

real(kind=real_umphys), intent(out) ::                                         &
 ff(pfield,levels,nvis)   ! out Vis prob at processed levels
!                          !     (decimal fraction).

! --------------------------------------------------------------------
!    Workspace usage----------------------------------------------------
real(kind=real_umphys) ::                                                      &
                         ! "Automatic" arrays on Cray.
 p(pfield),                                                                    &
 qt(pfield),                                                                   &
                         ! total of cloud water and vapour
 qs(pfield),                                                                   &
                         ! Saturated spec humidity for temp T
 qt_thresh(pfield),                                                            &
                         ! modified qt
 bs
!  Local, including save'd, storage------------------------------------
integer :: k,i,j     ! Loop counters: K - vertical level index.
!                       !                I - horizontal field index.
                        !                J - Vis threshold index.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='FOG_FR'

!-----------------------------------------------------------------------
!  Subroutine structure :
!  Loop round levels to be processed.
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
do k = 1, levels

  !-----------------------------------------------------------------------
  !  1. Calculate Pressure and initialise temporary arrays
  !-----------------------------------------------------------------------

  do i = 1, pfield
    p(i)=p_layer(i,k)
    qt(i)=q(i,k)+qcl(i,k)
  end do ! Loop over points

  !-----------------------------------------------------------------------
  !  2.  Calculate total water threshold corresponding to visibility
  !      Since Qs is needed more than once, pre-calculate and pass it
  !-----------------------------------------------------------------------

  call qsat_wat(qs,t(:,k),p,pfield)

  do j = 1, nvis

    call vistoqt( vis(1,k,j), qs, aerosol(1,k), l_murk,                        &
                pfield, qt_thresh )

    !-----------------------------------------------------------------------
    !  3.  Calculate the width of the distribution in total water space, bs:
    !
    !            bs = ( 1 - RHcrit ) * qs(T)
    !
    !-----------------------------------------------------------------------

    do i = 1, pfield

      bs = (1.0-rhcrit(k)) * qs(i)

      !=======================================================================
      !  4.  Calculate the fraction of states in a triangular
      !      distribution which exceed the total water threshold.
      !=======================================================================

      !-----------------------------------------------------------------------
      !  4.1 If total water threshold value is less than the total water value
      !      minus the width of the distribution, then all of the states have
      !      a total water value exceeding the threshold, so set the
      !      visibility fraction to 1.0
      !-----------------------------------------------------------------------

      if ( qt_thresh(i)  <=  qt(i)-bs ) then

        ff(i,k,j) = 1.0

        !-----------------------------------------------------------------------
        !  4.2 If total water threshold value is greater than the total water
        !      value minus the width of the distribution, but less than the
        !      total water value then the visibility fraction, VF, is given by:
        !
        !                                                     2
        !                              ( qt       - qt + bs  )
        !             VF = 1.0 - 0.5 * (    thresh           )
        !                              ( ------------------- )
        !                              (          bs         )
        !
        !-----------------------------------------------------------------------

      else if ( qt_thresh(i)  >   qt(i)-bs .and.                               &
                qt_thresh(i)  <=  qt(i) ) then

        ff(i,k,j) = 1.0 - 0.5 *                                                &
             (( qt_thresh(i) - qt(i) + bs )/ bs)**2

        !-----------------------------------------------------------------------
        !  4.3 If total water threshold value is greater than the total water
        !      value, but less than the total water value plus the width of the
        !      distribution, then the visibility fraction, VF, is given by:
        !
        !                                               2
        !                        ( qt + bs - qt        )
        !             VF = 0.5 * (             thresh  )
        !                        ( ------------------- )
        !                        (          bs         )
        !
        !-----------------------------------------------------------------------

      else if ( qt_thresh(i)  >   qt(i) .and.                                  &
                qt_thresh(i)  <=  qt(i)+bs    ) then

        ff(i,k,j)= 0.5 * (( qt(i) + bs - qt_thresh(i))/bs)**2

        !-----------------------------------------------------------------------
        !  4.4 If total water threshold value is greater than the total water
        !      value plus the width of the distribution, then non of the states
        !      have a total water value exceeding the threshold, so set the
        !      visibility fraction to 0.0
        !-----------------------------------------------------------------------

      else

        ff(i,k,j) = 0.0

      end if

    end do ! Loop over PFIELD I

  end do ! Loop over VIS J

end do ! Loop over levels

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine fog_fr
end module fog_fr_mod
