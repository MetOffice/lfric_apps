! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_pc2_mod

use um_types, only: real_umphys

implicit none

! Description:
!  Water tracer (working array) structure used in the cloud scheme (PC2) plus
!  related routines
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

type :: wtrac_pc2_type

  ! This contains WATER fields
  real(kind=real_umphys), allocatable :: q_cond(:,:,:)
                                         ! Amount of water condensing
  real(kind=real_umphys), allocatable :: q_old(:,:,:)
                                         ! Water vapour prior to phase change
  real(kind=real_umphys), allocatable :: qcl_old(:,:,:)
                                         ! Water liquid prior to phase change

end type

type(wtrac_pc2_type) :: wtrac_pc2

character(len=*), parameter, private :: ModuleName = 'WTRAC_PC2_MOD'
contains

! Subroutine Interface:
subroutine wtrac_pc2_store(q, qcl)

use atm_fields_bounds_mod,  only: tdims
use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim

implicit none

! Description:
!  Store the q and qcl fields which are input to pc2_homog_plus_turb and
!  pc2_initiate for use in the water tracer calculations

! q before PC2 phase change
real(kind=real_umphys), intent(in) :: q(tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end)

! qcl before PC2 phase change
real(kind=real_umphys), intent(in) :: qcl(tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end,           &
                                          1:tdims%k_end)

integer :: i, j, k      ! Counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_PC2_STORE'

! End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate and set arrays
allocate(wtrac_pc2%q_old(tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                         1:tdims%k_end))

allocate(wtrac_pc2%qcl_old(tdims%i_start:tdims%i_end,                          &
                           tdims%j_start:tdims%j_end,                          &
                           1:tdims%k_end))

allocate(wtrac_pc2%q_cond(tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end,                           &
                          1:tdims%k_end))

!$OMP  PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                              &
!$OMP  SHARED(tdims, wtrac_pc2, q, qcl) private(i,j,k)
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      wtrac_pc2%q_old(i,j,k)   = q(i,j,k)
      wtrac_pc2%qcl_old(i,j,k) = qcl(i,j,k)
      wtrac_pc2%q_cond(i,j,k)  = 0.0
    end do
  end do
end do
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_pc2_store

! -------------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_pc2_store_dealloc()

use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim

implicit none

! Description:
!  Deallocate PC2 water tracer working arrays

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_PC2_STORE_DEALLOC'

! End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

deallocate(wtrac_pc2%q_cond)
deallocate(wtrac_pc2%qcl_old)
deallocate(wtrac_pc2%q_old)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_pc2_store_dealloc

end module wtrac_pc2_mod
