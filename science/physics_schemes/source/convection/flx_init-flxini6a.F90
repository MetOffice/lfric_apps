! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate initial downdraught mass flux
!
module flx_init_6a_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'FLX_INIT_6A_MOD'
contains

subroutine flx_init_6a(npnts,kct,iccb,icct,flx,flx_dd_k,bddi                   &
                   ,flx_strt)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

!-----------------------------------------------------------------------------
! Description :
!    Calculate initial downdraught mass flux
!
!   Method    : See Unified Model documentation paper 27.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.3.
!-----------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts                & ! Vector length
 ,kct                    ! Convective cloud top

integer, intent(in) ::                                                         &
  iccb(npnts)          & ! Convective cloud base level
 ,icct(npnts)            ! Convective cloud top level


logical, intent(in) ::                                                         &
  bddi(npnts)            ! Mask for points where downdraught may initiate

real(kind=real_umphys), intent(in) ::                                          &
  flx(npnts,kct)       ! updraught mass flux (Pa/s)

real(kind=real_umphys), intent(out) ::                                         &
  flx_dd_k(npnts)       & ! Downdraught mass flux of layer k (Pa/s)
 ,flx_strt(npnts)         ! Updraught mass flux at level downdraught starts
                          ! (Pa/s)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer ::                                                                     &
  i               & ! loop counter
 ,kddref            ! reference level for downdraught massflux

real(kind=real_umphys) :: flxscale
                    ! The scaling factor for the initial downdraught massflux


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='FLX_INIT_6A'

!----------------------------------------------------------------------
! Calculate downdraught massflux based of a reference level which is
! 3/4 cloud depth
!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


flxscale=0.1

do i=1,npnts
  if (bddi(i)) then
    kddref = int(iccb(i) + 0.75*(icct(i) - iccb(i)))
    if (kddref  >=  icct(i)-1) then
      kddref=icct(i)-1
    end if
    flx_strt(i) = flx(i,kddref)
    flx_dd_k(i) = flx_strt(i) * flxscale
  end if
end do




if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine flx_init_6a
end module flx_init_6a_mod
