! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lift_par_phase_chg_wtrac_mod

use um_types, only: real_umphys

implicit none

! Description:
!  Updates water tracers for the phase changes calculated in lift_par
!
! Method:
!  Calculates the amount of water changing phase and then the
!  corresponding amount of water tracer.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private ::                                        &
                                 ModuleName = 'LIFT_PAR_PHASE_CHG_WTRAC_MOD'

contains

! Subroutine interface:
subroutine lift_par_phase_chg_wtrac(k, npnts, n_wtrac, ni, idx, q1, q2,        &
                                    qchange, process, q1_old, q2_old, wtrac_p)

use wtrac_conv_mod,          only: conv_p_wtrac_type
use wtrac_all_phase_chg_mod, only: wtrac_all_phase_chg

use parkind1,       only: jprb, jpim
use yomhook,        only: lhook, dr_hook

implicit none

! Subroutine arguments:
integer,intent(in) :: k              ! Present model layer number
integer,intent(in) :: npnts          ! Vector length
integer,intent(in) :: n_wtrac        ! No. of water tracers
integer,intent(in) :: ni             ! No. of working points

integer,intent(in) :: idx(ni)        ! Working points index

real(kind=real_umphys),intent(in) :: q1(npnts)
                                     ! Parcel mixing ratio in layer k+1
                                     !  for phase type 1
real(kind=real_umphys),intent(in) :: q2(npnts)
                                     ! Parcel mixing ratio in layer k+1
                                     !  for phase type 2
real(kind=real_umphys),intent(in) :: qchange(npnts)
                                     ! Amount of water changing phase

character(len=*), intent(in) :: process  ! Process to be modelled

real(kind=real_umphys),intent(in out) :: q1_old(ni)
                                     ! Parcel mixing ratio in layer k+1
                                     ! for phase type 1 before phase change
                                     ! (compressed to working points only)
real(kind=real_umphys),intent(in out) :: q2_old(ni)
                                     ! Parcel mixing ratio in layer k+1
                                     ! for phase type 2 before phase change
                                     ! (compressed to working points only)

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)
                                     ! Structure containing parcel
                                     ! water tracer fields
! Local variables
integer :: i, m, i_wt       ! Loop counters

! Work arrays for water tracer (compress to working points only)
real(kind=real_umphys) :: q1_new(ni)
                                     ! Parcel mixing ratio in layer k+1
                                     ! for phase type 1 after phase change
real(kind=real_umphys) :: q2_new(ni)
                                     ! Parcel mixing ratio in layer k+1
                                     ! for phase type 2 after phase change
real(kind=real_umphys) :: qchange_c(ni)
                                     ! Amount of water changing phase
real(kind=real_umphys) :: q1_wtrac(ni,n_wtrac)
                                     ! Parcel water tracer(kg/kg) at layer k+1
                                     ! for phase type 1
real(kind=real_umphys) :: q2_wtrac(ni,n_wtrac)
                                     ! Parcel water tracer(kg/kg) at layer k+1
                                     ! for phase type 2
real(kind=real_umphys) :: qchange_wtrac(ni,n_wtrac)
                                     ! Water tracer change caused by phase
                                     ! change

character(len=3) :: q1_type, q2_type

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LIFT_PAR_PHASE_CHG_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Working arrays for k+1 values
if (process == 'Frezkp1') then
  ! Freezing (liq -> ice) or melting (ice -> liq)
  do m = 1, ni
    i = idx(m)
    q1_new(m)    = q1(i)
    q2_new(m)    = q2(i)
    qchange_c(m) = qchange(i)
    do i_wt = 1, n_wtrac
      q1_wtrac(m,i_wt) = wtrac_p(i_wt)%qcl(i,k+1)
      q2_wtrac(m,i_wt) = wtrac_p(i_wt)%qcf(i,k+1)
    end do
  end do
  q1_type = 'liq'
  q2_type = 'ice'

else if (process == 'Qlkp1') then

  ! Evaporation (liq -> vap) or Condensation (vap -> liq)
  do m = 1, ni
    i = idx(m)
    q1_new(m)    = q1(i)
    q2_new(m)    = q2(i)
    qchange_c(m) = qchange(i)
    do i_wt = 1, n_wtrac
      q1_wtrac(m,i_wt) = wtrac_p(i_wt)%q(i,k+1)
      q2_wtrac(m,i_wt) = wtrac_p(i_wt)%qcl(i,k+1)
    end do
  end do
  q1_type = 'vap'
  q2_type = 'liq'

else if (process == 'Qfkp1') then

  ! Sublimation (ice -> vap) or Deposition (vap -> ice)
  do m = 1, ni
    i = idx(m)
    q1_new(m)    = q1(i)
    q2_new(m)    = q2(i)
    qchange_c(m) = qchange(i)
    do i_wt = 1, n_wtrac
      q1_wtrac(m,i_wt) = wtrac_p(i_wt)%q(i,k+1)
      q2_wtrac(m,i_wt) = wtrac_p(i_wt)%qcf(i,k+1)
    end do
  end do
  q1_type = 'vap'
  q2_type = 'ice'
end if

! -----------------------------------------------------------------
! Calculate change in water tracer for phase change
! -----------------------------------------------------------------

call wtrac_all_phase_chg(ni, n_wtrac,  q1_old, q2_old, qchange_c,              &
                         q1_new, q2_new, q1_type, q2_type, 'two_way',          &
                         q1_wtrac, q2_wtrac, qchange_wtrac = qchange_wtrac)

! Update k+1 values in main structure

if (process == 'Frezkp1') then
  do i_wt = 1, n_wtrac
    do m = 1, ni
      i = idx(m)
      wtrac_p(i_wt)%qcl(i,k+1) = q1_wtrac(m,i_wt)
      wtrac_p(i_wt)%qcf(i,k+1) = q2_wtrac(m,i_wt)
      wtrac_p(i_wt)%Frezkp1(i) = qchange_wtrac(m,i_wt)
    end do
  end do
end if

if (process == 'Qlkp1') then
  do i_wt = 1, n_wtrac
    do m = 1, ni
      i = idx(m)
      wtrac_p(i_wt)%q(i,k+1)   = q1_wtrac(m,i_wt)
      wtrac_p(i_wt)%qcl(i,k+1) = q2_wtrac(m,i_wt)
      wtrac_p(i_wt)%Qlkp1(i)   = qchange_wtrac(m,i_wt)
    end do
  end do
end if

if (process == 'Qfkp1') then
  do i_wt = 1, n_wtrac
    do m = 1, ni
      i = idx(m)
      wtrac_p(i_wt)%q(i,k+1)   = q1_wtrac(m,i_wt)
      wtrac_p(i_wt)%qcf(i,k+1) = q2_wtrac(m,i_wt)
      wtrac_p(i_wt)%Qfkp1(i)   = qchange_wtrac(m,i_wt)
    end do
  end do
end if

if (process == 'Frezkp1') then
  ! Store current values of water ahead of next call to this routine
  do m = 1, ni
    q1_old(m) = q1_new(m)
    q2_old(m) = q2_new(m)
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lift_par_phase_chg_wtrac

end module lift_par_phase_chg_wtrac_mod
