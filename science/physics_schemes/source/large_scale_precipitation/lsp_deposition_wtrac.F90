! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsp_deposition_wtrac_mod

use um_types, only: real_lsprec

implicit none

! Description:
!  Update water tracers for deposition/sublimation.
!
! Method:
!  Note, as deposition and sublimation can both occur at the same point and
!  as isotopic fractionation will occur for deposition but not
!  sublimation, these processes have to be done separately here.
!  They are done sequentially here to reduce the potential of numerical
!  error building up.  This is different to how sublimation and deposition
!  are calculated in lsp_deposition where they are done in parallel.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private ::                                        &
                           ModuleName = 'LSP_DEPOSITION_WTRAC_MOD'
contains

! Subroutine Interface:
subroutine lsp_deposition_wtrac(points, q, qcl, qcf, t,                        &
                                wtrac_mp_cpr_old, wtrac_mp_cpr)

use free_tracers_inputs_mod,    only: n_wtrac
use lsprec_mod,                 only: zero
use wtrac_all_phase_chg_mod,    only: wtrac_all_phase_chg
use wtrac_mphys_mod,            only: mp_cpr_old_wtrac_type, mp_cpr_wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: points      ! No. of points

! 'Normal' water after all lsp_deposition phase changes
real (kind=real_lsprec), intent(in) :: q(points)
real (kind=real_lsprec), intent(in) :: qcl(points)
real (kind=real_lsprec), intent(in) :: qcf(points)
real (kind=real_lsprec), intent(in) :: t(points)

! Structure containing normal water fields before call to lsp_deposition
type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Structure containing water tracer fields
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

! Local variables
integer :: i, i_wt                  ! Loop counters


real (kind=real_lsprec) :: q_wtrac(points, n_wtrac)
                                    ! Water tracer vapour
real (kind=real_lsprec) :: qcl_wtrac(points, n_wtrac)
                                    ! Water tracer liquid condensate
real (kind=real_lsprec) :: qcf_wtrac(points, n_wtrac)
                                    ! Water tracer ice condensate
real (kind=real_lsprec) :: q_new(points)
                                    ! q updated sequentially
real (kind=real_lsprec) :: qcl_new(points)
                                    ! qcl updated sequentially
real (kind=real_lsprec) :: qcf_new(points)
                                    ! qcf updated sequentially
real (kind=real_lsprec) :: q_old(points)
                                    ! q updated for sublimation
real (kind=real_lsprec) :: qcf_old(points)
                                    ! qcf updated sublimation and qcl->qcf

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_DEPOSITION_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  do i = 1, points
    q_wtrac(i,i_wt)   = wtrac_mp_cpr(i_wt)%q(i)
    qcl_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcl(i)
    qcf_wtrac(i,i_wt) = wtrac_mp_cpr(i_wt)%qcf(i)
  end do
end do

!------------------------------------------------------------------------
!    Sublimation qcf -> q
!------------------------------------------------------------------------

! Calculate qcf after just sublimation
do i = 1, points
  qcf_new(i) =  wtrac_mp_cpr_old%qcf(i) - wtrac_mp_cpr_old%qchange(i)
  q_new(i)   =  wtrac_mp_cpr_old%q(i)   + wtrac_mp_cpr_old%qchange(i)
end do

call wtrac_all_phase_chg(points, n_wtrac,                                      &
                        wtrac_mp_cpr_old%qcf, wtrac_mp_cpr_old%q,              &
                        wtrac_mp_cpr_old%qchange,                              &
                        qcf_new, q_new, 'ice', 'vap', 'one_way',               &
                        qcf_wtrac, q_wtrac)

!-------------------------------------------------------------------
!   Deposition  qcl -> qcf
!   Possibly some fractionation here for isotopes as this is the
!   Bergeron process (qcl -> q -> qcf)??
!-------------------------------------------------------------------

! Calculate qcf after sublimation and liquid deposition
do i = 1, points
  qcl_new(i) =  wtrac_mp_cpr_old%qcl(i) - wtrac_mp_cpr_old%depos_l(i)
  qcf_new(i) =  qcf_new(i)              + wtrac_mp_cpr_old%depos_l(i)
end do

! (Note any change to qcl here is just due to deposition)
call wtrac_all_phase_chg(points, n_wtrac,                                      &
                       wtrac_mp_cpr_old%qcl, wtrac_mp_cpr_old%qcf,             &
                       wtrac_mp_cpr_old%depos_l,                               &
                       qcl_new, qcf_new, 'liq', 'ice', 'one_way',              &
                       qcl_wtrac, qcf_wtrac)

!-------------------------------------------------------------------------
!   Deposition  q -> qcf
!   Fractionation will occur here for isotopes
!-------------------------------------------------------------------------

! Store water fields before deposition and calculate water quantities
! after deposition
do i = 1, points
  q_old(i)     =  q_new(i)
  qcf_old(i)   =  qcf_new(i)
  q_new(i)     =  q_old(i)    - wtrac_mp_cpr_old%depos_v(i)
  qcf_new(i)   =  qcf_old(i)  + wtrac_mp_cpr_old%depos_v(i)
end do

call wtrac_all_phase_chg(points, n_wtrac,                                      &
                         q_old, qcf_old,                                       &
                         wtrac_mp_cpr_old%depos_v,                             &
                         q_new, qcf_new, 'vap', 'ice', 'one_way',              &
                         q_wtrac, qcf_wtrac)

! Now update main fields
do i_wt = 1, n_wtrac
  do i = 1, points
    wtrac_mp_cpr(i_wt)%q(i)   = q_wtrac(i,i_wt)
    wtrac_mp_cpr(i_wt)%qcl(i) = qcl_wtrac(i,i_wt)
    wtrac_mp_cpr(i_wt)%qcf(i) = qcf_wtrac(i,i_wt)
  end do
end do


! Store current values of q, qcl, and qcf for future use
do i = 1, points
  wtrac_mp_cpr_old%q(i)       = q(i)
  wtrac_mp_cpr_old%qcl(i)     = qcl(i)
  wtrac_mp_cpr_old%qcf(i)     = qcf(i)
  wtrac_mp_cpr_old%t(i)       = t(i)
  wtrac_mp_cpr_old%qchange(i) = zero
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_deposition_wtrac

end module lsp_deposition_wtrac_mod
