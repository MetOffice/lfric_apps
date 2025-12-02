! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate whether downdraught is able to continue
!
module termdd_mod

use um_types, only: real_umphys

implicit none

! Description: Calculate whether downdraught is able to continue
!              Calculate buoyancy
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.

character(len=*), parameter, private :: ModuleName = 'TERMDD_MOD'

contains

subroutine termdd (npnts, k, bdd_start                                         &
                   , thdd_k, qdd_k, the_k, qe_k, ppn_mix_dd                    &
                   , b_dd_end, bdd_on)

use cv_run_mod, only:                                                          &
    dd_opt

use planet_constants_mod, only: c_virtual

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

integer, intent(in) ::                                                         &
  npnts                &  ! Vector length
 ,k                       ! Present model layer

logical, intent(in) ::                                                         &
  bdd_start(npnts)        ! Mask for those points where downdraught may occur
                          ! in layer k-1

real(kind=real_umphys), intent(in) ::                                          &
  thdd_k(npnts)        & ! Potential temperature of downdraught in layer k (K)

 ,qdd_k(npnts)         & ! Mixing ratio of downdraught in layer k  (kg/kg)

 ,the_k(npnts)         & ! Potential temperature of environment in layer k (K)

 ,qe_k(npnts)          & ! Mixing ratio of environment in layer k   (kg/kg)

 ,ppn_mix_dd(npnts)      ! precipitation mixing ratio (kg/kg)

logical, intent(out) ::                                                        &
  b_dd_end(npnts)       & ! Mask for those points where downdraught is
                          ! terminating
 ,bdd_on(npnts)           ! Mask for those points where downdraught continues
                          ! layer k-1 (as bdd_start here)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i             ! Loop counters

real(kind=real_umphys) ::                                                      &
  buoy1      &  ! Buoyancy of parcel

 ,thdd_v     &  ! Used in calculation of buoyancy

 ,the_v         ! Used in calculation of buoyancy


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='TERMDD'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (dd_opt == 1) then
  !-----------------------------------------------------------------------
  ! Check if parcel still negatively buoyant such that downdraught
  ! can continue to next layer
  !-----------------------------------------------------------------------

  do i=1,npnts
    thdd_v = thdd_k(i)*(1.0+c_virtual*qdd_k(i)) /(1.0+ppn_mix_dd(i))
    the_v  = the_k(i) *(1.0+c_virtual*qe_k(i))
    buoy1  = thdd_v - the_v

    !-----------------------------------------------------------------------
    ! Calculate state of downdraught
    !-----------------------------------------------------------------------

    if (bdd_start(i) .and. buoy1 >  0.5) then
      bdd_on(i) = .false.
    else if (buoy1 >  0.5 .or. k == 2) then
      b_dd_end(i) = .true.
    end if
  end do

else
  !-----------------------------------------------------------------------
  ! Check if parcel still negatively buoyant such that downdraught
  ! can continue to next layer
  !-----------------------------------------------------------------------

  do i=1,npnts
    thdd_v = thdd_k(i)*(1.0+c_virtual*qdd_k(i))
    the_v  = the_k(i)* (1.0+c_virtual*qe_k(i))
    buoy1  = thdd_v - the_v

    !-----------------------------------------------------------------------
    ! Calculate state of downdraught
    !-----------------------------------------------------------------------

    if (bdd_start(i) .and. buoy1 >  0.5) then
      bdd_on(i) = .false.
    else if (buoy1 >  0.5 .or. k == 2) then
      b_dd_end(i) = .true.
    end if
  end do

end if         ! dd_opt test

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine termdd

end module termdd_mod
