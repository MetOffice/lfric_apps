! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_bl_mod

use um_types, only: r_bl

implicit none

! Description:
!  Water tracer (working array) structure used in the boundary layer plus
!  related routines
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName='WTRAC_BL_MOD'

type :: bl_wtrac_type

  ! WATER TRACER fields used in both the explicit and implicit sections:

  real(r_bl), allocatable :: qcf_use(:,:,:) ! qcf used in bl scheme
                                            ! (depends on l_noice_in_turb)
  real(r_bl), allocatable :: qw(:,:,:)      ! Total water content (q+qcl+qcf)
  real(r_bl), allocatable :: fqw(:,:,:)     ! Water tracer flux between layers
                                            ! (kg per square metre per sec).
                                            !  FQW(,1) is total water flux
                                            !  from surface, 'E'.
  real(r_bl), allocatable :: fqw_surft(:,:)       ! Surface fqw for land tiles
                                            ! (Not used in model at present)
  real(r_bl), allocatable :: fqw_sicat(:,:,:)     ! Surface fqw over sea ice
                                            ! (Not used in model at present)

  ! WATER TRACER fields used in explicit section only:

  real(r_bl), allocatable :: grad_q_adj(:,:)! Humidity gradient adjustment
                                            ! for non-local mixing in unstable
                                            ! turbulent boundary layer.
  real(r_bl), allocatable :: fq_nt(:,:,:)   ! Non-turubulent moisture flux
  real(r_bl), allocatable :: fq_nt_dscb(:,:)! Non-turbulent moisture flux at
                                            ! base of the DSC layer
  real(r_bl), allocatable :: totqf_zh(:,:)  ! Total moisture fluxes at
  real(r_bl), allocatable :: totqf_zhsc(:,:)! inversions

  ! WATER TRACER fields used in implicit section only:

  real(r_bl), allocatable :: dqw(:,:,:)     ! BL increment to q field
  real(r_bl), allocatable :: dqw_nt(:,:,:)  ! Non-turbulent increment to qw
  real(r_bl), allocatable :: dqw1_1(:,:)    ! Coefficient neede for implicit
                                            ! coupling at level k_blend_tq

  real(r_bl), allocatable :: qcf_latest_use(:,:,:)
  real(r_bl), allocatable :: q_earliest(:,:,:)    ! Fields used for storage in
  real(r_bl), allocatable :: qcl_earliest(:,:,:)  ! ni_imp_ctl
  real(r_bl), allocatable :: qcf_earliest(:,:,:)

end type

! Water tracer structure containing boundary layer fields
type(bl_wtrac_type), allocatable, save :: wtrac_bl(:)

! Logical that controls the water tracer calculations in the implicit part
! of the boundary layer scheme
logical :: l_wtrac_bl = .false.

contains

! Subroutine Interface:
subroutine wtrac_alloc_bl(land_pts, ntiles, bl_levels)

!
! Description:
! Allocate water tracer structure wtrac_bl which contains fields
! used in the boundary layer.  Allocate fields that are used in both the
! explicit and implicit parts and initialise water tracer total water field.
!

use atm_fields_bounds_mod,   only: tdims, pdims, tdims_l
use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use jules_sea_seaice_mod,    only: nice_use

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: land_pts    ! Number of land points being processed,
                                   ! can be 0  !!
integer, intent(in) :: ntiles      ! Number of land tiles

integer, intent(in) :: bl_levels   ! Max no. of "boundary" levels

integer :: i_wt    ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_BL'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

allocate(wtrac_bl(n_wtrac))

if (l_wtrac) then
  do i_wt = 1, n_wtrac

    allocate(wtrac_bl(i_wt)%qcf_use(tdims_l%i_start:tdims_l%i_end,             &
             tdims_l%j_start:tdims_l%j_end, tdims_l%k_start:bl_levels))

    allocate(wtrac_bl(i_wt)%qw(tdims%i_start:tdims%i_end,                      &
             tdims%j_start:tdims%j_end, bl_levels))

    allocate(wtrac_bl(i_wt)%fqw(pdims%i_start:pdims%i_end,                     &
             pdims%j_start:pdims%j_end, bl_levels))

    allocate(wtrac_bl(i_wt)%fqw_surft(land_pts,ntiles))

    allocate(wtrac_bl(i_wt)%fqw_sicat(pdims%i_start:pdims%i_end,               &
             pdims%j_start:pdims%j_end, nice_use))

  end do
else
  do i_wt = 1, n_wtrac
    allocate(wtrac_bl(i_wt)%qcf_use(1,1,1))
    allocate(wtrac_bl(i_wt)%qw(1,1,1))
    allocate(wtrac_bl(i_wt)%fqw(1,1,1))
    allocate(wtrac_bl(i_wt)%fqw_surft(1,1))
    allocate(wtrac_bl(i_wt)%fqw_sicat(1,1,1))
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

return
end subroutine wtrac_alloc_bl

!--------------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_bl()

!
! Description:
! Deallocate water tracer arrays in structure wtrac_bl which contains fields
! used in the boundary layer
!

use free_tracers_inputs_mod, only: n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer :: i_wt    ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_BL'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac

  deallocate(wtrac_bl(i_wt)%fqw_sicat)
  deallocate(wtrac_bl(i_wt)%fqw_surft)
  deallocate(wtrac_bl(i_wt)%fqw)
  deallocate(wtrac_bl(i_wt)%qw)
  deallocate(wtrac_bl(i_wt)%qcf_use)

end do

deallocate(wtrac_bl)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

return
end subroutine wtrac_dealloc_bl

! --------------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_alloc_bl_exp(bl_levels)

!
! Description:
! Allocate water tracer arrays in structure wtrac_bl that are used only in
! the explicit part of the boundary layer
!

use atm_fields_bounds_mod,   only: pdims
use free_tracers_inputs_mod, only: l_wtrac, n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: bl_levels     ! Max no. of "boundary" levels

integer :: i_wt       ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_BL_EXP'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_wtrac) then
  do i_wt = 1, n_wtrac
    allocate(wtrac_bl(i_wt)%grad_q_adj(pdims%i_start:pdims%i_end,              &
                                       pdims%j_start:pdims%j_end))
    allocate(wtrac_bl(i_wt)%fq_nt(pdims%i_start:pdims%i_end,                   &
                                  pdims%j_start:pdims%j_end,bl_levels+1))
    allocate(wtrac_bl(i_wt)%fq_nt_dscb(pdims%i_start:pdims%i_end,              &
                                       pdims%j_start:pdims%j_end))
    allocate(wtrac_bl(i_wt)%totqf_zh(pdims%i_start:pdims%i_end,                &
                                     pdims%j_start:pdims%j_end))
    allocate(wtrac_bl(i_wt)%totqf_zhsc(pdims%i_start:pdims%i_end,              &
                                       pdims%j_start:pdims%j_end))
  end do

else

  do i_wt = 1, n_wtrac
    allocate(wtrac_bl(i_wt)%grad_q_adj(1,1))
    allocate(wtrac_bl(i_wt)%fq_nt(1,1,1))
    allocate(wtrac_bl(i_wt)%fq_nt_dscb(1,1))
    allocate(wtrac_bl(i_wt)%totqf_zh(1,1))
    allocate(wtrac_bl(i_wt)%totqf_zhsc(1,1))
  end do

end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

return
end subroutine wtrac_alloc_bl_exp

!--------------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_bl_exp()

!
! Description:
! Deallocate water tracer arrays in structure wtrac_bl which are only
! used in the explicit part of the boundary layer
!

use free_tracers_inputs_mod, only: n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer :: i_wt    ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_BL_EXP'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  deallocate(wtrac_bl(i_wt)%totqf_zhsc)
  deallocate(wtrac_bl(i_wt)%totqf_zh)
  deallocate(wtrac_bl(i_wt)%fq_nt_dscb)
  deallocate(wtrac_bl(i_wt)%fq_nt)
  deallocate(wtrac_bl(i_wt)%grad_q_adj)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

return
end subroutine wtrac_dealloc_bl_exp

! --------------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_alloc_bl_imp(bl_levels)

!
! Description:
! Allocate water tracer arrays in structure wtrac_bl that are used only in
! the implicit part of the boundary layer
!

use atm_fields_bounds_mod,   only: tdims
use free_tracers_inputs_mod, only: n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: bl_levels     ! Max no. of "boundary" levels

integer :: i_wt    ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_BL_IMP'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_wtrac_bl) then
  do i_wt = 1, n_wtrac

    allocate(wtrac_bl(i_wt)%dqw(tdims%i_start:tdims%i_end,                     &
             tdims%j_start:tdims%j_end,bl_levels))

    allocate(wtrac_bl(i_wt)%dqw_nt(tdims%i_start:tdims%i_end,                  &
             tdims%j_start:tdims%j_end,bl_levels))

    allocate(wtrac_bl(i_wt)%dqw1_1(tdims%i_start:tdims%i_end,                  &
             tdims%j_start:tdims%j_end))

    allocate(wtrac_bl(i_wt)%qcf_latest_use(tdims%i_start:tdims%i_end,          &
             tdims%j_start:tdims%j_end,bl_levels))

    allocate(wtrac_bl(i_wt)%q_earliest(tdims%i_start:tdims%i_end,              &
             tdims%j_start:tdims%j_end,tdims%k_end))

    allocate(wtrac_bl(i_wt)%qcl_earliest(tdims%i_start:tdims%i_end,            &
             tdims%j_start:tdims%j_end,tdims%k_end))

    allocate(wtrac_bl(i_wt)%qcf_earliest(tdims%i_start:tdims%i_end,            &
             tdims%j_start:tdims%j_end,tdims%k_end))

  end do

else

  do i_wt = 1, n_wtrac
    allocate(wtrac_bl(i_wt)%dqw(1,1,1))
    allocate(wtrac_bl(i_wt)%dqw_nt(1,1,1))
    allocate(wtrac_bl(i_wt)%dqw1_1(1,1))
    allocate(wtrac_bl(i_wt)%qcf_latest_use(1,1,1))
    allocate(wtrac_bl(i_wt)%q_earliest(1,1,1))
    allocate(wtrac_bl(i_wt)%qcl_earliest(1,1,1))
    allocate(wtrac_bl(i_wt)%qcf_earliest(1,1,1))
  end do
end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

return
end subroutine wtrac_alloc_bl_imp

!--------------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_bl_imp()

!
! Description:
! Deallocate water tracer arrays in structure wtrac_bl which are only
! used in the implicit part of the boundary layer
!

use free_tracers_inputs_mod, only: n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer :: i_wt    ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_BL_IMP'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac

  deallocate(wtrac_bl(i_wt)%qcf_earliest)
  deallocate(wtrac_bl(i_wt)%qcl_earliest)
  deallocate(wtrac_bl(i_wt)%q_earliest)
  deallocate(wtrac_bl(i_wt)%qcf_latest_use)
  deallocate(wtrac_bl(i_wt)%dqw1_1)
  deallocate(wtrac_bl(i_wt)%dqw_nt)
  deallocate(wtrac_bl(i_wt)%dqw)

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

return
end subroutine wtrac_dealloc_bl_imp

end module wtrac_bl_mod
