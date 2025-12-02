! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates the amount of water loading for the parcel and the environment

module water_loading_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Calculates the amount of water loading for the parcel and the environment
!
! Method:
!   Calculates the amount of water loading for the parcel and the environment
!   depending upon the water loading option selected.
!   1) Ignore water loading
!   2) Include all the condensate in the water loading
!   3) Not yet coded: Include only the convective condensate estimated to
!      remain after precipitation
!
!   The arguments are named in terms of layer k but the routine can also
!   be applied to layer k+1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


character(len=*), parameter, private :: ModuleName='WATER_LOADING_MOD'

contains

subroutine water_loading (npnts, qclek, qcfek, qclpk, qcfpk,                   &
                          watldek,watldpk,idx,ni)

use cv_run_mod, only: cnv_wat_load_opt
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

!-----------------------------------------------------------------------
! Subroutine arguments:
!-----------------------------------------------------------------------

! Arguments with intent in:

integer, intent(in) :: npnts ! Vector length

real(kind=real_umphys),intent(in) :: qclek(npnts) ! env. liquid condensate in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfek(npnts) ! env. frozen condensate in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qclpk(npnts) ! parcel liquid condensate in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfpk(npnts) ! parcel frozen condensate in layer k (kg/kg)

! Arguments with intent out:

real(kind=real_umphys), intent(out) :: watldek(npnts)
                                    ! Environment water loading in layer k+1
real(kind=real_umphys), intent(out) :: watldpk(npnts)
                                    ! Parcel water loading in layer k+1

! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i,m               ! loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WATER_LOADING'


!---------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

select case (cnv_wat_load_opt)
case (0)
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    watldek(i) = 0.0
    watldpk(i) = 0.0
  end do
case (1)
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    watldek(i) = qclek(i) + qcfek(i)
    watldpk(i) = qclpk(i) + qcfpk(i)
  end do
end select


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine water_loading
end module water_loading_mod
