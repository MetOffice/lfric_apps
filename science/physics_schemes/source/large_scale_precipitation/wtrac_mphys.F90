! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_mphys_mod

! Water tracer microphysics: type definitions and related routines

use um_types,                only: real_umphys, real_lsprec

implicit none

! Description:
!  Water tracer (working array) structures used in the microphysics plus
!  related routines.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! Microphysics working fields (these are all WATER TRACER fields)
! These are 3D fields that are updated each microphysics sub-timestep

type :: mp_wtrac_type

  real(kind=real_umphys), allocatable :: q(:,:,:)         ! Vapour
  real(kind=real_umphys), allocatable :: qcl(:,:,:)       ! Liq condensate
  real(kind=real_umphys), allocatable :: qcf(:,:,:)       ! Ice condensate
  real(kind=real_umphys), allocatable :: qr(:,:,:)        ! Prognostic rain
  real(kind=real_umphys), allocatable :: lsrain(:,:)      ! Largescale sfc rain
  real(kind=real_umphys), allocatable :: lssnow(:,:)      ! Largescale sfc snow
  real(kind=real_umphys), allocatable :: droplet_flux(:,:)! Water drop flux
                                                          ! (kg m-2 s-1)

end type

! Microphysics compressed fields (these are all WATER TRACER fields).
!
! These are compressed 1D arrays containing only the grid points where
! microphysics needs to be calculated.  These fields are set up on individual
! OMP threads in ls_ppnc.F90.
! They are then used and updated in the microphysics calculations (in lsp_ice)
! before being 'scattered' back to the 3D working fields (wtrac_mp).
! The equivalent water fields are held in the super array 'spr'.
!
! Note, these fields use 'real_lsprec' precision, which is the
! scheme-specific precision for the microphysics and can differ from the
! precision used the other UM physics ('real_umphys').

type :: mp_cpr_wtrac_type

  real(kind=real_lsprec), allocatable :: q(:)             ! Vapour
  real(kind=real_lsprec), allocatable :: qcl(:)           ! Liq condensate
  real(kind=real_lsprec), allocatable :: qcf(:)           ! Ice condensate
  real(kind=real_lsprec), allocatable :: qrain(:)         ! Prognostic rain
  real(kind=real_lsprec), allocatable :: lsrain(:)        ! Largescale sfc rain
  real(kind=real_lsprec), allocatable :: lssnow(:)        ! Largescale sfc snow
  real(kind=real_lsprec), allocatable :: droplet_flux(:)  ! Water drop flux

  real(kind=real_lsprec), allocatable :: rainratet(:)     ! Cummulative fall
  real(kind=real_lsprec), allocatable :: snowratet(:)     ! out of hydrometeor
                                                          ! within iterations

end type

! Structure to store copies of WATER fields prior to phase changes
! and also the amount of water changing phase for various processes.
! These are compressed arrays with the same compression method as detailed
! above for the arrays in mp_cpr_wtrac_type.

type :: mp_cpr_old_wtrac_type

  real(kind=real_lsprec), allocatable :: q(:)       ! Vapour
  real(kind=real_lsprec), allocatable :: qcl(:)     ! Liq condensate
  real(kind=real_lsprec), allocatable :: qcf(:)     ! Ice condensate
  real(kind=real_lsprec), allocatable :: qrain(:)   ! Prognostic rain
  real(kind=real_lsprec), allocatable :: t(:)       ! Temperature
  real(kind=real_lsprec), allocatable :: droplet_flux(:)
                                                    ! Water drop flux
  real(kind=real_lsprec), allocatable :: qchange(:) ! Phase change amount
  real(kind=real_lsprec), allocatable :: het_q(:)   ! Amount of heterogenous
                                                    ! nucleation of vapour
  real(kind=real_lsprec), allocatable :: hom_qr(:)  ! Amount of homogeneous
                                                    ! nucleation of rain
  real(kind=real_lsprec), allocatable :: depos_l(:) ! Amount of deposition
                                                    ! (qcl -> qcf)
  real(kind=real_lsprec), allocatable :: depos_v(:) ! Amount of deposition
                                                    ! (q -> qcf)

end type

character(len=*),   parameter, private :: ModuleName='WTRAC_MPHYS_MOD'

contains

! Subroutine Interface:
subroutine mphys_init_wtrac(wtrac_n, wtrac_mp)

! Description:
!   Set up water tracer working arrays for microphysics

use atm_fields_bounds_mod,   only: tdims
use mphys_inputs_mod,        only: l_mcr_qrain
use free_tracers_inputs_mod, only: n_wtrac
use water_tracers_mod,       only: wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Water tracer structures
type(wtrac_type),    intent(in)     :: wtrac_n(n_wtrac)

type(mp_wtrac_type), intent(in out) :: wtrac_mp(n_wtrac)


! Local variables
integer :: i,j,k,i_wt            ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MPHYS_INIT_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate arrays
do i_wt = 1, n_wtrac
  allocate(wtrac_mp(i_wt)%q(tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end,                       &
                            1 : tdims%k_end))
  allocate(wtrac_mp(i_wt)%qcl(tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                              1 : tdims%k_end))
  allocate(wtrac_mp(i_wt)%qcf(tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                              1 : tdims%k_end))
  if (l_mcr_qrain) then
    allocate(wtrac_mp(i_wt)%qr(tdims%i_start : tdims%i_end,                    &
                               tdims%j_start : tdims%j_end,                    &
                               1 : tdims%k_end))
  else
    allocate(wtrac_mp(i_wt)%qr(1,1,1))
  end if
  allocate(wtrac_mp(i_wt)%lsrain(tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end))
  allocate(wtrac_mp(i_wt)%lssnow(tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end))
  allocate(wtrac_mp(i_wt)%droplet_flux(tdims%i_start : tdims%i_end,            &
                                       tdims%j_start : tdims%j_end))
end do

! Set to 'n' values
do i_wt = 1, n_wtrac
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( i_wt, tdims, wtrac_mp, wtrac_n, l_mcr_qrain) private( i, j, k )
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        wtrac_mp(i_wt)%q(i,j,k)   = wtrac_n(i_wt)%q(i,j,k)
        wtrac_mp(i_wt)%qcl(i,j,k) = wtrac_n(i_wt)%qcl(i,j,k)
        wtrac_mp(i_wt)%qcf(i,j,k) = wtrac_n(i_wt)%qcf(i,j,k)
        if (l_mcr_qrain) wtrac_mp(i_wt)%qr(i,j,k) = wtrac_n(i_wt)%qr(i,j,k)
      end do
    end do
  end do
!$OMP end PARALLEL do
end do


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mphys_init_wtrac

! -------------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_mphys_update_inc(wtrac_mp, wtrac_n, wtrac_as)

! Description:
!   Update water tracer increments fields for microphysics

use atm_fields_bounds_mod,   only: tdims
use mphys_inputs_mod,        only: l_mcr_qrain
use water_tracers_mod,       only: wtrac_type
use wtrac_atm_step_mod,      only: atm_step_wtrac_type
use free_tracers_inputs_mod, only: n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Water tracer structures
type(wtrac_type),          intent(in)     :: wtrac_n(n_wtrac)
type(mp_wtrac_type),       intent(in out) :: wtrac_mp(n_wtrac)
type(atm_step_wtrac_type), intent(in out) :: wtrac_as(n_wtrac)

! Local variables
integer :: i,j,k,i_wt            ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_MPHYS_UPDATE_INC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Update increments for water tracers

do i_wt = 1, n_wtrac
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( i_wt, tdims, wtrac_as, wtrac_mp, wtrac_n, l_mcr_qrain)           &
!$OMP private(i, j, k)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

        ! q_inc = q_work - q_n + q_inc
        wtrac_as(i_wt)%q_star(i,j,k) = wtrac_mp(i_wt)%q(i,j,k)                 &
                - wtrac_n(i_wt)%q(i,j,k) + wtrac_as(i_wt)%q_star(i,j,k)
        wtrac_as(i_wt)%qcl_star(i,j,k) = wtrac_mp(i_wt)%qcl(i,j,k)             &
                - wtrac_n(i_wt)%qcl(i,j,k) + wtrac_as(i_wt)%qcl_star(i,j,k)
        wtrac_as(i_wt)%qcf_star(i,j,k) = wtrac_mp(i_wt)%qcf(i,j,k)             &
                - wtrac_n(i_wt)%qcf(i,j,k) + wtrac_as(i_wt)%qcf_star(i,j,k)
        if (l_mcr_qrain) then
          wtrac_as(i_wt)%qr_star(i,j,k) = wtrac_mp(i_wt)%qr(i,j,k)             &
               - wtrac_n(i_wt)%qr(i,j,k) + wtrac_as(i_wt)%qr_star(i,j,k)
        end if

      end do
    end do
  end do
!$OMP end PARALLEL do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_mphys_update_inc

! -------------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_mphys_dealloc(wtrac_mp)

! Description:
!   Deallocate water tracer working arrays used for microphysics
!   (in microphys_ctl)

use free_tracers_inputs_mod, only: n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Water tracer structures
type(mp_wtrac_type), intent(in out) :: wtrac_mp(n_wtrac)

! Local variable
integer :: i_wt            ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_MPHYS_DEALLOC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Dellocate arrays
do i_wt = 1, n_wtrac
  deallocate(wtrac_mp(i_wt)%droplet_flux)
  deallocate(wtrac_mp(i_wt)%lssnow)
  deallocate(wtrac_mp(i_wt)%lsrain)
  deallocate(wtrac_mp(i_wt)%qr)
  deallocate(wtrac_mp(i_wt)%qcf)
  deallocate(wtrac_mp(i_wt)%qcl)
  deallocate(wtrac_mp(i_wt)%q)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_mphys_dealloc

! -------------------------------------------------------------------------

! Subroutine Interface:
subroutine ls_ppnc_gather_wtrac(ip, i1, ix, level, wtrac_mp, wtrac_mp_cpr)

! Description:
!  Set up microphysics compressed arrays for water tracer variables

use atm_fields_bounds_mod,   only: tdims
use lsprec_mod,              only: zero
use mphys_inputs_mod,        only: l_mcr_qrain
use free_tracers_inputs_mod, only: n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: ip       ! Number of points
integer, intent(in) :: i1       ! Starting point
integer, intent(in) :: ix (  tdims%i_len * tdims%j_len , 2 )
                                ! Gather/scatter index
integer, intent(in) :: level    ! Vertical level

! Water tracer structures
type(mp_wtrac_type),  intent(in)     :: wtrac_mp(n_wtrac)

type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)


! Local variables
integer :: i,ii,ij,i_wt            ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LS_PPNC_GATHER_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate arrays
do i_wt = 1, n_wtrac
  allocate(wtrac_mp_cpr(i_wt)%q(ip))
  allocate(wtrac_mp_cpr(i_wt)%qcl(ip))
  allocate(wtrac_mp_cpr(i_wt)%qcf(ip))
  allocate(wtrac_mp_cpr(i_wt)%qrain(ip)) ! Used as diagnostic if l_mcr_rain=F
  allocate(wtrac_mp_cpr(i_wt)%lsrain(ip))
  allocate(wtrac_mp_cpr(i_wt)%lssnow(ip))
  allocate(wtrac_mp_cpr(i_wt)%droplet_flux(ip))
  allocate(wtrac_mp_cpr(i_wt)%rainratet(ip))
  allocate(wtrac_mp_cpr(i_wt)%snowratet(ip))
end do

do i=1, ip

  ii   =  ix(i+i1-1,1)
  ij   =  ix(i+i1-1,2)

  do i_wt = 1, n_wtrac
    wtrac_mp_cpr(i_wt)%q(i)   =                                                &
                real(wtrac_mp(i_wt)%q(ii,ij,level), kind=real_lsprec)
    wtrac_mp_cpr(i_wt)%qcl(i) =                                                &
                real(wtrac_mp(i_wt)%qcl(ii,ij,level), kind=real_lsprec)
    wtrac_mp_cpr(i_wt)%qcf(i) =                                                &
                real(wtrac_mp(i_wt)%qcf(ii,ij,level), kind=real_lsprec)
    wtrac_mp_cpr(i_wt)%lsrain(i) =                                             &
                real(wtrac_mp(i_wt)%lsrain(ii,ij), kind=real_lsprec)
    wtrac_mp_cpr(i_wt)%lssnow(i) =                                             &
                real(wtrac_mp(i_wt)%lssnow(ii,ij), kind=real_lsprec)
    wtrac_mp_cpr(i_wt)%droplet_flux(i) =                                       &
                real(wtrac_mp(i_wt)%droplet_flux(ii,ij), kind=real_lsprec)
    if (l_mcr_qrain) then
      wtrac_mp_cpr(i_wt)%qrain(i) =                                            &
                real(wtrac_mp(i_wt)%qr(ii,ij,level), kind=real_lsprec)
    else
      wtrac_mp_cpr(i_wt)%qrain(1) = zero
    end if
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_ppnc_gather_wtrac

! -------------------------------------------------------------------------

! Subroutine Interface:
subroutine ls_ppnc_scatter_wtrac(ip, i1, ix, level, wtrac_mp_cpr, wtrac_mp)

! Description:
!   Scatter water tracer fields and deallocate compressed arrays

use atm_fields_bounds_mod,   only: tdims
use mphys_inputs_mod,        only: l_mcr_qrain
use free_tracers_inputs_mod, only: n_wtrac

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: ip       ! Number of points
integer, intent(in) :: i1       ! Starting point
integer, intent(in) :: ix (  tdims%i_len * tdims%j_len , 2 )
                                ! Gather/scatter index
integer, intent(in) :: level    ! Vertical level

type(mp_wtrac_type),  intent(in out) :: wtrac_mp(n_wtrac)

type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)


! Local variables
integer :: i,ii,ij,i_wt            ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LS_PPNC_SCATTER_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i=1, ip

  ii   =  ix(i+i1-1,1)
  ij   =  ix(i+i1-1,2)

  do i_wt = 1, n_wtrac
    wtrac_mp(i_wt)%q(ii,ij,level)   =                                          &
                real(wtrac_mp_cpr(i_wt)%q(i), kind=real_umphys)
    wtrac_mp(i_wt)%qcl(ii,ij,level) =                                          &
                real(wtrac_mp_cpr(i_wt)%qcl(i), kind=real_umphys)
    wtrac_mp(i_wt)%qcf(ii,ij,level) =                                          &
                real(wtrac_mp_cpr(i_wt)%qcf(i), kind=real_umphys)
    wtrac_mp(i_wt)%lsrain(ii,ij) =                                             &
                real(wtrac_mp_cpr(i_wt)%lsrain(i), kind=real_umphys)
    wtrac_mp(i_wt)%lssnow(ii,ij) =                                             &
                real(wtrac_mp_cpr(i_wt)%lssnow(i), kind=real_umphys)
    wtrac_mp(i_wt)%droplet_flux(ii,ij) =                                       &
                real(wtrac_mp_cpr(i_wt)%droplet_flux(i), kind=real_umphys)
    if (l_mcr_qrain) then
      wtrac_mp(i_wt)%qr(ii,ij,level) =                                         &
                real(wtrac_mp_cpr(i_wt)%qrain(i), kind=real_umphys)
    end if
  end do
end do

! Deallocate compressed arrays
do i_wt = 1, n_wtrac
  deallocate(wtrac_mp_cpr(i_wt)%snowratet)
  deallocate(wtrac_mp_cpr(i_wt)%rainratet)
  deallocate(wtrac_mp_cpr(i_wt)%droplet_flux)
  deallocate(wtrac_mp_cpr(i_wt)%lssnow)
  deallocate(wtrac_mp_cpr(i_wt)%lsrain)
  deallocate(wtrac_mp_cpr(i_wt)%qrain)
  deallocate(wtrac_mp_cpr(i_wt)%qcf)
  deallocate(wtrac_mp_cpr(i_wt)%qcl)
  deallocate(wtrac_mp_cpr(i_wt)%q)
end do


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_ppnc_scatter_wtrac

end module wtrac_mphys_mod
