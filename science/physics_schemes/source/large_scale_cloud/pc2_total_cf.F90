! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  PC2 Cloud Scheme: Change in total cloud fraction

module pc2_total_cf_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_TOTAL_CF_MOD'
contains

subroutine pc2_total_cf(                                                       &
!      Number of points
 points,                                                                       &
!      Input fields
 cfl, cff, deltacl, deltacf,                                                   &
!      Input Output fields
 cf)

use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

implicit none

! Purpose:
!   Update the total cloud fraction due to changes in
!   liquid and ice cloud fractions.
!
! Method:
!   Assumes a random overlap of the forced cloud with already existing
!   conditions in the gridbox. See Annex D of the PC2 cloud scheme
!   documentation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Subroutine Arguments:------------------------------------------------
integer ::                                                                     &
                      !, intent(in)
 points
!       No. of points being processed.

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 cfl(points),                                                                  &
!       Liquid cloud fraction
   cff(points),                                                                &
!       Ice cloud fraction
   deltacl(points),                                                            &
!       Change in liquid cloud fraction
   deltacf(points)
!       Change in ice cloud fraction

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 cf(points)
!       Total cloud fraction

!  External functions:

!  Local scalars--------------------------------------------------------

!  (a)  Scalars effectively expanded to workspace by the Cray (using
!       vector registers).
real(kind=real_umphys) ::                                                      &
 deltac_1, & ! Change in total cloud fraction due to liquid cloud
 deltac_2, & ! Change in total cloud fraction due to ice cloud
 cf_1,     & ! Input cf + deltac_1 (Used to avoid rounding error)
 cf_2        ! Input cf + deltac_2 (Used to avoid rounding error)

!  (b)  Others
integer :: i ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_TOTAL_CF'

!- End of Header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------
! Points_do1:
do i = 1, points

  ! ----------------------------------------------------------------------
  ! 1. Update total cloud fraction.
  ! ----------------------------------------------------------------------

  ! Calculate change in total cloud fraction due to a change in liquid
  ! cloud fraction. This depends upon the sign of the change of liquid
  ! cloud fraction.

  if (deltacl(i) > 0.0 .and. cfl(i) < 1.0) then

    !     Random overlap
    !     deltac_1 = deltacl(i) *(1.0 - cf(i))/(1.0 - cfl(i))

    !     Minimum overlap
    deltac_1 = min(         deltacl(i), (1.0-cf(i)) )
    cf_1     = min( cf(i) + deltacl(i),  1.0        )

  else if (deltacl(i) < 0.0 .and. cfl(i) > 0.0) then

    !     Random overlap
    !     deltac_1 = deltacl(i) * (cf(i)-cff(i)) / cfl(i)

    !     Minimum overlap
    deltac_1 = max(         deltacl(i), (cff(i)-cf(i)) )
    cf_1     = max( cf(i) + deltacl(i),  cff(i)        )

  else
    deltac_1 = 0.0
    cf_1     = cf(i)
  end if

  ! Calculate change in total cloud fraction due to a change in ice
  ! cloud fraction. This depends upon the sign of the change of ice
  ! cloud fraction.

  if (deltacf(i) > 0.0 .and. cff(i) < 1.0) then

    !     Random overlap
    !     deltac_2 = deltacf(i) *(1.0 - cf(i))/(1.0 - cff(i))

    !     Minimum overlap
    deltac_2 = min(         deltacf(i), (1.0-cf(i)) )
    cf_2     = min( cf(i) + deltacf(i),  1.0        )

  else if (deltacf(i)  <   0.0 .and. cff(i)  >   0.0) then

    !     Random overlap
    !     deltac_2 = deltacf(i) * (cf(i)-cfl(i)) / cff(i)

    !     Minimum overlap
    deltac_2 = max(         deltacf(i), (cfl(i)-cf(i)) )
    cf_2     = max( cf(i) + deltacf(i),  cfl(i)        )

  else
    deltac_2 = 0.0
    cf_2     = cf(i)
  end if

  ! Sum the two changes

  cf(i) = cf(i) + deltac_1 + deltac_2

  ! Avoid rounding error in limiting cases

  if (deltac_1 == 0.0) cf(i) = cf_2
  if (deltac_2 == 0.0) cf(i) = cf_1

  ! For minimum overlap we need to check that the total cloud
  ! fraction is constrained within 0 and 1

  cf(i) = max(  min( cf(i),1.0 )  , 0.0)

  ! Points_do1:
end do

! End of the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_total_cf
end module pc2_total_cf_mod
