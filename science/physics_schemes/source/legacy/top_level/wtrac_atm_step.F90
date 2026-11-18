! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_atm_step_mod

use um_types, only: real_umphys

implicit none

! Description:
! Water tracer (working array) structure used in top_level code plus
! related routines

! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

type :: atm_step_wtrac_type

  ! 'Star' fields for specific humidity (kg/kg)
  real(kind=real_umphys), allocatable :: q_star(:,:,:)    ! q (vapour)
  real(kind=real_umphys), allocatable :: qcl_star(:,:,:)  ! qcl (liq condensate)
  real(kind=real_umphys), allocatable :: qcf_star(:,:,:)  ! qcf (ice condensate)
  real(kind=real_umphys), allocatable :: qcf2_star(:,:,:) ! qcf (2nd ice)
  real(kind=real_umphys), allocatable :: qr_star(:,:,:)   ! qrain (rain)
  real(kind=real_umphys), allocatable :: qgr_star(:,:,:)  ! qgraupel (graupel)

  ! 'Star' fields for mass mixing ratio (kg/kg)
  real(kind=real_umphys), allocatable :: m_star(:,:,:)
  real(kind=real_umphys), allocatable :: mcl_star(:,:,:)
  real(kind=real_umphys), allocatable :: mcf_star(:,:,:)
  real(kind=real_umphys), allocatable :: mcf2_star(:,:,:)
  real(kind=real_umphys), allocatable :: mr_star(:,:,:)
  real(kind=real_umphys), allocatable :: mgr_star(:,:,:)

  ! Precipitation diagnostics (ls = large scale, conv = convection)
  real(kind=real_umphys), allocatable :: ls_rain(:,:)
  real(kind=real_umphys), allocatable :: ls_snow(:,:)
  real(kind=real_umphys), allocatable :: conv_rain(:,:)
  real(kind=real_umphys), allocatable :: conv_snow(:,:)

  ! Microphysics tendencies  (Used in explicit boundary layer code)
  real(kind=real_umphys), allocatable :: micro_tends(:,:,:)

end type

type(atm_step_wtrac_type), allocatable :: wtrac_as(:)

! Logical to control the water tracer calculations in the ENDGAME dynamics
! outer loop.  Water tracer calculations (i.e. advection, implicit boundary
! layer and convection) only need to be done on the final outer loop.
! This logical is set within the outer loop in atm_step_4A.

logical :: l_wtrac_final_loop = .false.

character(len=*), parameter, private :: ModuleName='WTRAC_ATM_STEP_MOD'

contains

! Subroutine Interface:
subroutine init_star_wtrac()

!
! Description:
! Allocate water tracer arrays in structure wtrac_as which contains fields
! used in atm_step
!

use atm_fields_bounds_mod,   only: tdims
use mphys_inputs_mod,        only: l_mcr_qcf2, l_mcr_qrain,                    &
                                   graupel_option, no_graupel
use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use nlsizes_namelist_mod,    only: bl_levels

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer :: i_wt    ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='INIT_STAR_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_wtrac) then
  allocate(wtrac_as(n_wtrac))
else
  allocate(wtrac_as(1))
end if

if (l_wtrac) then
  do i_wt = 1, n_wtrac

    allocate(                                                                  &
    ! q_star
          wtrac_as(i_wt)%q_star(tdims%i_start:tdims%i_end,                     &
                 tdims%j_start:tdims%j_end,                                    &
                 tdims%k_start:tdims%k_end),                                   &
          wtrac_as(i_wt)%qcl_star(tdims%i_start:tdims%i_end,                   &
                 tdims%j_start:tdims%j_end,                                    &
                 tdims%k_start:tdims%k_end),                                   &
          wtrac_as(i_wt)%qcf_star(tdims%i_start:tdims%i_end,                   &
                 tdims%j_start:tdims%j_end,                                    &
                 tdims%k_start:tdims%k_end),                                   &
    ! m_star
          wtrac_as(i_wt)%m_star(tdims%i_start:tdims%i_end,                     &
                 tdims%j_start:tdims%j_end,                                    &
                 tdims%k_start:tdims%k_end),                                   &
          wtrac_as(i_wt)%mcl_star(tdims%i_start:tdims%i_end,                   &
                 tdims%j_start:tdims%j_end,                                    &
                 tdims%k_start:tdims%k_end),                                   &
          wtrac_as(i_wt)%mcf_star(tdims%i_start:tdims%i_end,                   &
                 tdims%j_start:tdims%j_end,                                    &
                 tdims%k_start:tdims%k_end),                                   &
    ! ls_rain, ls_snow, conv_rain, conv_snow
          wtrac_as(i_wt)%ls_rain(tdims%i_start:tdims%i_end,                    &
                 tdims%j_start:tdims%j_end),                                   &
          wtrac_as(i_wt)%ls_snow(tdims%i_start:tdims%i_end,                    &
                 tdims%j_start:tdims%j_end),                                   &
          wtrac_as(i_wt)%conv_rain(tdims%i_start:tdims%i_end,                  &
                 tdims%j_start:tdims%j_end),                                   &
          wtrac_as(i_wt)%conv_snow(tdims%i_start:tdims%i_end,                  &
                 tdims%j_start:tdims%j_end),                                   &
    ! micro_tends
          wtrac_as(i_wt)%micro_tends(tdims%i_start:tdims%i_end,                &
                 tdims%j_start:tdims%j_end,bl_levels) )

  end do

else
  do i_wt = 1, n_wtrac
    allocate(wtrac_as(i_wt)%q_star(1,1,1))
    allocate(wtrac_as(i_wt)%qcl_star(1,1,1))
    allocate(wtrac_as(i_wt)%qcf_star(1,1,1))
    allocate(wtrac_as(i_wt)%m_star(1,1,1))
    allocate(wtrac_as(i_wt)%mcl_star(1,1,1))
    allocate(wtrac_as(i_wt)%mcf_star(1,1,1))
    allocate(wtrac_as(i_wt)%ls_rain(1,1))
    allocate(wtrac_as(i_wt)%ls_snow(1,1))
    allocate(wtrac_as(i_wt)%conv_rain(1,1))
    allocate(wtrac_as(i_wt)%conv_snow(1,1))
    allocate(wtrac_as(i_wt)%micro_tends(1,1,1))
  end do
end if

  ! Allocate additional microphysics variables to full size
  ! if in use, otherwise allocate minimum amount of space

if ( l_mcr_qcf2 .and. l_wtrac) then
  do i_wt = 1, n_wtrac
    allocate( wtrac_as(i_wt)%qcf2_star(tdims%i_start:tdims%i_end,              &
                       tdims%j_start:tdims%j_end,                              &
                       tdims%k_start:tdims%k_end) )
    allocate( wtrac_as(i_wt)%mcf2_star(tdims%i_start:tdims%i_end,              &
                       tdims%j_start:tdims%j_end,                              &
                       tdims%k_start:tdims%k_end) )
  end do
else
  do i_wt = 1, n_wtrac
    allocate( wtrac_as(i_wt)%qcf2_star(1,1,1) )
    allocate( wtrac_as(i_wt)%mcf2_star(1,1,1) )
  end do
end if

if ( l_mcr_qrain .and. l_wtrac ) then
  do i_wt = 1, n_wtrac
    allocate(  wtrac_as(i_wt)%qr_star(tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,                             &
                        tdims%k_start:tdims%k_end) )
    allocate(  wtrac_as(i_wt)%mr_star(tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,                             &
                        tdims%k_start:tdims%k_end) )
  end do
else
  do i_wt = 1, n_wtrac
    allocate( wtrac_as(i_wt)%qr_star(1,1,1) )
    allocate( wtrac_as(i_wt)%mr_star(1,1,1) )
  end do
end if

if ( graupel_option > no_graupel .and. l_wtrac ) then
  do i_wt = 1, n_wtrac
    allocate( wtrac_as(i_wt)%qgr_star(tdims%i_start:tdims%i_end,               &
                         tdims%j_start:tdims%j_end,                            &
                         tdims%k_start:tdims%k_end) )
    allocate( wtrac_as(i_wt)%mgr_star(tdims%i_start:tdims%i_end,               &
                         tdims%j_start:tdims%j_end,                            &
                         tdims%k_start:tdims%k_end) )
  end do
else
  do i_wt = 1, n_wtrac
    allocate( wtrac_as(i_wt)%qgr_star(1,1,1) )
    allocate( wtrac_as(i_wt)%mgr_star(1,1,1) )
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

return
end subroutine init_star_wtrac

end module wtrac_atm_step_mod

