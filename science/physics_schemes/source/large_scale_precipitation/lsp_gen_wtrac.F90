! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsp_gen_wtrac_mod

use um_types, only: real_lsprec

implicit none

! Description:
!  Generic routine to update water tracers for a single phase change in
!  the microphysics.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'LSP_GEN_WTRAC_MOD'
contains

! Subroutine Interface:
subroutine lsp_gen_wtrac(points, q1, q2, t, process, q1_old, q2_old, t_old,    &
                         q_change, wtrac_mp_cpr)

! Description:
!  Generic routine to update water tracers for a single phase change in
!  the microphysics from phase type 1 to phase type 2.
!  Phase changes are caused by:
!  1. 'Rim'      = Riming (qcl -> qcf)
!  2. 'Rim orog' = Enhanced riming by orographic water (q -> qcf)
!  3. 'Capture'  = Capture of raindrops by ice (qrain -> qcf)
!  4. 'Evapsnow' = Sublimation of melting snow (qcf -> q)
!  5. 'Melt'     = Melting of ice (qcf -> qrain)
!  6. 'Evaprain' = Evaporation of rain (qrain -> q) (which includes
!                  isotopic fractionation and exchange)
!  7. 'Accretion'= Accretion of cloud droplets by rain drops (qcl -> qrain)
!  8. 'Accr orog'= Enhanced accretion by orographic water (q -> qrain)
!  9. 'Autoc'    = Autoconversion of liquid to rain (qcl -> qrain)
!

use free_tracers_inputs_mod,    only: n_wtrac
use lsprec_mod,                 only: zero
use wtrac_all_phase_chg_mod,    only: wtrac_all_phase_chg
use wtrac_mphys_mod,            only: mp_cpr_wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: points      ! No. of points

! 'Normal' water after phase change (q1 -> q2)
real (kind=real_lsprec), intent(in) :: q1(points)   ! Phase type 1
real (kind=real_lsprec), intent(in) :: q2(points)   ! Phase type 2
real (kind=real_lsprec), intent(in) :: t(points)    ! Temperature

character(len=*), intent(in) :: process     ! Process to be modelled

! 'Normal' water before phase change (in) and after phase change (out)
real (kind=real_lsprec), intent(in out) :: q1_old(points)    ! Phase type 1
real (kind=real_lsprec), intent(in out) :: q2_old(points)    ! Phase type 2
real (kind=real_lsprec), intent(in out) :: t_old(points)     ! Temperature

! Amount of water changing phase (set to 0 at end of routine)
real (kind=real_lsprec), intent(in out) :: q_change(points)

! Structure containing water tracer fields
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

! Local variables
integer :: i, i_wt         ! Loop counters

real (kind=real_lsprec) :: q1_wtrac(points, n_wtrac)  ! Water tracer for
                                                      ! phase type 1
real (kind=real_lsprec) :: q2_wtrac(points, n_wtrac)  ! Water tracer for
                                                      ! phase type 2

character(len=3) :: q1_type, q2_type                  ! Phase type names
                                                      ! (e.g. 'liq','ice')

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_GEN_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------
!    Set up fields depending on process
!------------------------------------------------------------------------

if (process == 'Rim') then           ! Riming (qcl -> qcf)
  ! (Note, there can also be a very small phase change direction qcf -> qcl)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcl(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
    end do
  end do
  q1_type = 'liq'
  q2_type = 'ice'

else if (process == 'Rim orog') then  ! Enhanced Riming due to orog (q -> qcf)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%q(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
    end do
  end do
  q1_type = 'vap'
  q2_type = 'ice'

else if (process == 'Capture') then  ! Capture of raindrops by ice
                                     ! (qrain -> qcf)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qrain(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
    end do
  end do
  q1_type = 'rai'
  q2_type = 'ice'

else if (process == 'Evapsnow') then  ! Sublimation of melting snow (qcf -> q)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%q(i)
    end do
  end do
  q1_type = 'ice'
  q2_type = 'vap'

else if (process == 'Melt') then      ! Melting of ice (qcf -> qrain)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qrain(i)
    end do
  end do
  q1_type = 'ice'
  q2_type = 'rai'

else if (process == 'Evaprain') then  ! Evaporation of rain (qrain -> q)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qrain(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%q(i)
    end do
  end do
  q1_type = 'rai'
  q2_type = 'vap'

else if (process == 'Accretion') then ! Accretion of cloud droplets by rain
                                      ! drops (qcl -> qrain)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcl(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qrain(i)
    end do
  end do
  q1_type = 'liq'
  q2_type = 'rai'

else if (process == 'Accr orog') then ! Enhanced accretion due to orog
                                      ! (q -> qrain)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%q(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qrain(i)
    end do
  end do
  q1_type = 'vap'
  q2_type = 'rai'

else if (process ==  'Autoc') then    ! Autoconversion of liquid to rain
                                       ! (qcl -> qrain)
  do i_wt = 1, n_wtrac
    do i = 1, points
      q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcl(i)
      q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qrain(i)
    end do
  end do
  q1_type = 'liq'
  q2_type = 'rai'

end if

! Calculate phase change for water tracers

call wtrac_all_phase_chg(points, n_wtrac, q1_old, q2_old, q_change,            &
                         q1, q2, q1_type, q2_type, 'two_way',                  &
                         q1_wtrac, q2_wtrac)

! Update main water tracer structure and store current normal water values

do i_wt = 1, n_wtrac
  if (process == 'Rim') then           ! Riming (qcl -> qcf)
    do i = 1, points
      wtrac_mp_cpr(i_wt)%qcl(i) = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%qcf(i) = q2_wtrac(i,i_wt)
    end do

  else if (process == 'Rim orog') then ! Enhanced riming due to orog (q -> qcf)
    do i = 1, points
      wtrac_mp_cpr(i_wt)%q(i)   = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%qcf(i) = q2_wtrac(i,i_wt)
    end do

  else if (process == 'Capture') then  ! Capture of raindrops by ice
    do i = 1, points                                       ! (qrain -> qcf)
      wtrac_mp_cpr(i_wt)%qrain(i) = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%qcf(i)   = q2_wtrac(i,i_wt)
    end do

  else if (process == 'Evapsnow') then  ! Sublimation of melting snow (qcf -> q)
    do i = 1, points
      wtrac_mp_cpr(i_wt)%qcf(i) = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%q(i)   = q2_wtrac(i,i_wt)
    end do

  else if (process == 'Melt') then      ! Melting of ice (qcf -> qrain)
    do i = 1, points
      wtrac_mp_cpr(i_wt)%qcf(i)   = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%qrain(i) = q2_wtrac(i,i_wt)
    end do

  else if (process == 'Evaprain') then  ! Evaporation of rain (qrain -> q)
    do i = 1, points
      wtrac_mp_cpr(i_wt)%qrain(i) = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%q(i)     = q2_wtrac(i,i_wt)
    end do

  else if (process == 'Accretion') then ! Accretion of cloud droplets by rain
                                       ! drops (qcl -> qrain)
    do i = 1, points
      wtrac_mp_cpr(i_wt)%qcl(i)   = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%qrain(i) = q2_wtrac(i,i_wt)
    end do

  else if (process == 'Accr orog') then ! Enhanced accretion due to orog
                                        ! (q -> qrain)
    do i = 1, points
      wtrac_mp_cpr(i_wt)%q(i)     = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%qrain(i) = q2_wtrac(i,i_wt)
    end do

  else if (process ==  'Autoc') then    ! Autoconversion of liquid to rain
                                       ! (qcl -> qrain)
    do i = 1, points
      wtrac_mp_cpr(i_wt)%qcl(i)   = q1_wtrac(i,i_wt)
      wtrac_mp_cpr(i_wt)%qrain(i) = q2_wtrac(i,i_wt)
    end do

  end if
end do

! Store current normal water values
do i = 1, points
  q1_old(i)   = q1(i)
  q2_old(i)   = q2(i)
  t_old(i)    = t(i)
  q_change(i) = zero
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_gen_wtrac

end module lsp_gen_wtrac_mod
