! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module pc2_checks_wtrac_mod

use um_types, only: real_umphys

implicit none

! Description:
!  Remove small water tracer condensate amounts which are below a minimum
!  level using same method as in pc2_checks.
! (Note, water tracers are updated directly in pc2_checks for any phase
!  changes in that routine.)
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'PC2_CHECKS_WTRAC_MOD'
contains

! Subroutine Interface:
subroutine pc2_checks_wtrac(wtrac)

use atm_fields_bounds_mod,    only: tdims
use free_tracers_inputs_mod,  only: n_wtrac
use mphys_inputs_mod,         only: l_mcr_qcf2
use water_tracers_mod,        only: wtrac_type, wtrac_info

use yomhook,                  only: lhook, dr_hook
use parkind1,                 only: jprb, jpim

implicit none

! Water tracer structure
type(wtrac_type), intent(in out) :: wtrac(n_wtrac)

! Local parameters

integer :: i,j,k,i_wt        ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_CHECKS_WTRAC'

! End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Additional checks on water tracers to ensure they remain positive
do i_wt = 1, n_wtrac
!$OMP  PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                              &
!$OMP  SHARED(i_wt, tdims, wtrac, l_mcr_qcf2, wtrac_info) private(i, j, k)
  do k = 1, tdims%k_end
    do j = tdims%j_start,tdims%j_end                  ! Not over halos
      do i = tdims%i_start,tdims%i_end
        if (wtrac(i_wt)%qcl(i,j,k) < wtrac_info(i_wt)%qcf_limit) then
          wtrac(i_wt)%q(i,j,k) = wtrac(i_wt)%q(i,j,k) + wtrac(i_wt)%qcl(i,j,k)
          wtrac(i_wt)%qcl(i,j,k) = 0.0
        end if
        if (wtrac(i_wt)%qcf(i,j,k) < wtrac_info(i_wt)%qcf_limit) then
          wtrac(i_wt)%q(i,j,k) = wtrac(i_wt)%q(i,j,k) + wtrac(i_wt)%qcf(i,j,k)
          wtrac(i_wt)%qcf(i,j,k) = 0.0
        end if
        if (l_mcr_qcf2) then
          if (wtrac(i_wt)%qcf2(i,j,k) < wtrac_info(i_wt)%qcf_limit) then
            wtrac(i_wt)%q(i,j,k) = wtrac(i_wt)%q(i,j,k) +                      &
                                   wtrac(i_wt)%qcf2(i,j,k)
            wtrac(i_wt)%qcf2(i,j,k) = 0.0
          end if
        end if
      end do
    end do
  end do
!$OMP end PARALLEL do
end do


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_checks_wtrac

end module pc2_checks_wtrac_mod
