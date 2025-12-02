! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsp_nucleation_wtrac_mod

use um_types, only: real_lsprec

implicit none

! Description:
!  Update water tracers for nucleation of ice particles
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'LSP_NUCLEATION_WTRAC_MOD'
contains

! Subroutine Interface:
subroutine lsp_nucleation_wtrac(points, q, qcl, qcf, qrain, t,                 &
                                wtrac_mp_cpr_old, wtrac_mp_cpr)

use free_tracers_inputs_mod, only: n_wtrac
use lsprec_mod,              only: zero
use wtrac_all_phase_chg_mod, only: wtrac_all_phase_chg
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type, mp_cpr_wtrac_type

use yomhook,                 only: lhook, dr_hook
use parkind1,                only: jprb, jpim

implicit none

integer, intent(in) :: points      ! No. of points

! 'Normal' water after nucleation phase changes
real (kind=real_lsprec), intent(in) :: q(points)
real (kind=real_lsprec), intent(in) :: qcl(points)
real (kind=real_lsprec), intent(in) :: qcf(points)
real (kind=real_lsprec), intent(in) :: qrain(points)
real (kind=real_lsprec), intent(in) :: t(points)

! Structure containing normal water fields before call to lsp_nucleation
type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Structure containing water tracer fields
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

! Local variables
integer :: i, i_wt            ! Loop counters

real (kind=real_lsprec) :: q1_wtrac(points, n_wtrac)
                              ! Water tracer for phase type 1 (q, qcl or qrain)
                              ! which undergoes phase change to phase type 2
real (kind=real_lsprec) :: q2_wtrac(points, n_wtrac)
                              ! Water tracer for phase type 2 (always qcf here)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_NUCLEATION_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Note, if isotopic fractionation is included here in the future, need to
! look carefully at the qcf values that are being used due to several
! phase changes causing changes.

!------------------------------------------------------------------------
!   Homogeneous or heterogeneous nucleation
!   qcl -> qcf
!-------------------------------------------------------------------------

do i_wt = 1, n_wtrac
  do i = 1, points
    q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcl(i)
    q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
  end do
end do

call wtrac_all_phase_chg(points, n_wtrac,                                      &
                               wtrac_mp_cpr_old%qcl, wtrac_mp_cpr_old%qcf,     &
                               wtrac_mp_cpr_old%qchange,                       &
                               qcl, qcf, 'liq', 'ice', 'one_way',              &
                               q1_wtrac, q2_wtrac)

do i_wt = 1, n_wtrac
  do i = 1, points
    wtrac_mp_cpr(i_wt)%qcl(i) = q1_wtrac(i,i_wt)
    wtrac_mp_cpr(i_wt)%qcf(i) = q2_wtrac(i,i_wt)
  end do
end do

!-----------------------------------------------------------------------
!   Homogeneous nucleation of rain
!   qrain -> qcf
!-----------------------------------------------------------------------

do i_wt = 1, n_wtrac
  do i = 1, points
    q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qrain(i)
    q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
  end do
end do

call wtrac_all_phase_chg(points, n_wtrac,                                      &
                               wtrac_mp_cpr_old%qrain, wtrac_mp_cpr_old%qcf,   &
                               wtrac_mp_cpr_old%hom_qr,                        &
                               qrain, qcf, 'rai', 'ice', 'one_way',            &
                               q1_wtrac, q2_wtrac)

do i_wt = 1, n_wtrac
  do i = 1, points
    wtrac_mp_cpr(i_wt)%qrain(i) = q1_wtrac(i,i_wt)
    wtrac_mp_cpr(i_wt)%qcf(i)   = q2_wtrac(i,i_wt)
  end do
end do

! ---------------------------------------------------------------------
! Heterogeneous nucleation (only uses vapour when all liquid has gone)
! q -> qcf
! ---------------------------------------------------------------------

do i_wt = 1, n_wtrac
  do i = 1, points
    q1_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%q(i)
    q2_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
  end do
end do

call wtrac_all_phase_chg(points, n_wtrac,                                      &
                               wtrac_mp_cpr_old%q, wtrac_mp_cpr_old%qcf,       &
                               wtrac_mp_cpr_old%het_q,                         &
                               q, qcf, 'vap', 'ice', 'one_way',                &
                               q1_wtrac, q2_wtrac)

do i_wt = 1, n_wtrac
  do i = 1, points
    wtrac_mp_cpr(i_wt)%q(i)   = q1_wtrac(i,i_wt)
    wtrac_mp_cpr(i_wt)%qcf(i) = q2_wtrac(i,i_wt)
  end do
end do


! Store current values of q, qcl, qcf and qrain for future use
do i = 1, points
  wtrac_mp_cpr_old%q(i)       = q(i)
  wtrac_mp_cpr_old%qcl(i)     = qcl(i)
  wtrac_mp_cpr_old%qcf(i)     = qcf(i)
  wtrac_mp_cpr_old%qrain(i)   = qrain(i)
  wtrac_mp_cpr_old%t(i)       = t(i)
  wtrac_mp_cpr_old%qchange(i) = zero
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_nucleation_wtrac

end module lsp_nucleation_wtrac_mod
