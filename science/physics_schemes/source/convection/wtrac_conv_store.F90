! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_conv_store_mod

use um_types, only: real_umphys

implicit none

! Description:
!  Water tracer (working array) structures used to store water fields in the
!  convection scheme plus related routines
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! Structure to store copies of normal water fields prior to phase changes
type :: conv_old_wtrac_type

  ! Used for phase changes in lift_par and cloud_w
  real(kind=real_umphys), allocatable :: qpkp1(:)
  real(kind=real_umphys), allocatable :: qclpkp1(:)
  real(kind=real_umphys), allocatable :: qcfpkp1(:)

  ! Used for phase changes in devap, pevp_bcb, crs_frzl and chg_phse
  real(kind=real_umphys), allocatable :: q_km1(:)
  real(kind=real_umphys), allocatable :: rain(:)
  real(kind=real_umphys), allocatable :: snow(:)
  real(kind=real_umphys), allocatable :: precip_frez(:)

end type

character(len=*), parameter, private :: ModuleName =  'WTRAC_CONV_STORE_MOD'

contains

! Subroutine Interface:
subroutine wtrac_alloc_conv_store1(npnts, ni, idx, qpkp1, qclpkp1, qcfpkp1,    &
                                   wtrac_conv_old)

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
!
! Description:
! Allocate space and store the values of water prior to the phase change in
! the convection routines lift_par and cloud_w, for use in water tracer
! calculations.
!
! Subroutine arguments
integer, intent(in) :: npnts    ! Vector length
integer, intent(in) :: ni       ! No. of working points

integer,intent(in) :: idx(ni)        ! Working points index

real(kind=real_umphys), intent(in) :: qpkp1(npnts)   ! Parcel q at k+1
real(kind=real_umphys), intent(in) :: qclpkp1(npnts) ! Parcel qcl at k+1
real(kind=real_umphys), intent(in) :: qcfpkp1(npnts) ! Parcel qcf at k+1

! Water tracer structure to store water values prior to phase change
type(conv_old_wtrac_type), intent(in out) :: wtrac_conv_old

! Local variables

integer :: i, m       ! Counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_CONV_STORE1'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate arrays in structure to store old values of 'normal' water
! Note, these are compressed fields (convecting points only)
allocate(wtrac_conv_old%qpkp1(ni))
allocate(wtrac_conv_old%qclpkp1(ni))
allocate(wtrac_conv_old%qcfpkp1(ni))

! Set values

do m = 1, ni
  i = idx(m)
  wtrac_conv_old%qpkp1(m)   = qpkp1(i)
  wtrac_conv_old%qclpkp1(m) = qclpkp1(i)
  wtrac_conv_old%qcfpkp1(m) = qcfpkp1(i)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_alloc_conv_store1

! -----------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_conv_store1(wtrac_conv_old)

!
! Description:
! Deallocate arrays used by water tracers in wtrac_conv_old which were
! allocated in wtrac_alloc_conv_store1
!
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
type(conv_old_wtrac_type), intent(in out) :: wtrac_conv_old

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_CONV_STORE1'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

deallocate(wtrac_conv_old%qcfpkp1)
deallocate(wtrac_conv_old%qclpkp1)
deallocate(wtrac_conv_old%qpkp1)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_dealloc_conv_store1

! ------------------------------------------------------------------------
! Subroutine Interface:
subroutine wtrac_alloc_conv_store2(npnts, q_km1, rain, snow, wtrac_conv_old)

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
!
! Description:
! Allocate space and store the values of water prior to the phase change in the
! convection routines devap, pevp_bcb, crs_frzl and chg_phse, for use in water
! tracer calculations.
!
! Subroutine arguments
integer, intent(in) :: npnts    ! Vector length

real(kind=real_umphys), intent(in) :: q_km1(npnts) ! Mixing ratio
                                                   ! in layer k-1 (kg/kg)
real(kind=real_umphys), intent(in) :: rain(npnts) ! Amount of rain descending
                                                  ! from k-1 to k-2 (kg/m**2/s)
real(kind=real_umphys), intent(in) :: snow(npnts) ! Amount of snow descending
                                                  ! from k-1 to k-2 (kg/m**2/s)

! Water tracer structure to store water values prior to phase change
type(conv_old_wtrac_type), intent(in out) :: wtrac_conv_old

! Local variables

integer :: i    ! Counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_CONV_STORE2'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate arrays in structure to store old values of 'normal' water

allocate(wtrac_conv_old%q_km1(npnts))
allocate(wtrac_conv_old%rain(npnts))
allocate(wtrac_conv_old%snow(npnts))
allocate(wtrac_conv_old%precip_frez(npnts))

! Set values
do i = 1, npnts
  wtrac_conv_old%q_km1(i)       = q_km1(i)
  wtrac_conv_old%rain(i)        = rain(i)
  wtrac_conv_old%snow(i)        = snow(i)
  wtrac_conv_old%precip_frez(i) = 0.0
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_alloc_conv_store2

! -----------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_conv_store2(wtrac_conv_old)

!
! Description:
! Deallocate arrays used by water tracers in wtrac_conv_old which were
! allocated in wtrac_alloc_conv_store2
!
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
type(conv_old_wtrac_type), intent(in out) :: wtrac_conv_old

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_CONV_STORE2'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

deallocate(wtrac_conv_old%q_km1)
deallocate(wtrac_conv_old%rain)
deallocate(wtrac_conv_old%snow)
deallocate(wtrac_conv_old%precip_frez)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_dealloc_conv_store2

end module wtrac_conv_store_mod
