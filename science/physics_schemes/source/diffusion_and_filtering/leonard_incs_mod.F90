! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to store increments from the Leonard terms.

module leonard_incs_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Module stores increments from the Leonard terms, so that they can
!   just be calculated on the first EG outer loop call to atmos_physics2
!   and then just re-used on subsequent calls.
!
! Method:
!   Module contains the allocatable increment arrays.  If Leonard terms
!   are active, the arrays are allocated by the contained alloc subroutine
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: diffusion_and_filtering
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName =                           &
                                        'LEONARD_INCS_MOD'

! Arrays to store Leonard term increments of:
real(kind=real_umphys), allocatable :: u_inc_leonard(:,:,:) ! Zonal wind
real(kind=real_umphys), allocatable :: v_inc_leonard(:,:,:) ! Meridional wind
real(kind=real_umphys), allocatable :: w_inc_leonard(:,:,:) ! Vertical wind
real(kind=real_umphys), allocatable :: thetal_inc_leonard(:,:,:)
                                                        ! Water potential temp
real(kind=real_umphys), allocatable :: qw_inc_leonard(:,:,:)
                                                        ! Total water content

contains

! Subroutine to allocate the above arrays; called in the first
! call to atmos_physics2, by atmo_physics2_alloc
subroutine leonard_incs_alloc()

use atm_fields_bounds_mod, only: udims, vdims, wdims, tdims
use nlsizes_namelist_mod,  only: bl_levels

use parkind1, only: jpim, jprb       !DrHook
use yomhook,  only: lhook, dr_hook   !DrHook

implicit none

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LEONARD_INCS_ALLOC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (.not. allocated (u_inc_leonard) ) then
  allocate( u_inc_leonard      ( udims%i_start:udims%i_end,                    &
                                 udims%j_start:udims%j_end, bl_levels ) )
end if

if (.not. allocated (v_inc_leonard) ) then
  allocate( v_inc_leonard      ( vdims%i_start:vdims%i_end,                    &
                                 vdims%j_start:vdims%j_end, bl_levels ) )
end if

if (.not. allocated (w_inc_leonard) ) then
  allocate( w_inc_leonard      ( wdims%i_start:wdims%i_end,                    &
                                 wdims%j_start:wdims%j_end, bl_levels ) )
end if

if (.not. allocated (thetal_inc_leonard) ) then
  allocate( thetal_inc_leonard ( tdims%i_start:tdims%i_end,                    &
                                 tdims%j_start:tdims%j_end, bl_levels ) )
end if

if (.not. allocated (qw_inc_leonard) ) then
  allocate( qw_inc_leonard     ( tdims%i_start:tdims%i_end,                    &
                                 tdims%j_start:tdims%j_end, bl_levels ) )
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine leonard_incs_alloc

end module leonard_incs_mod
