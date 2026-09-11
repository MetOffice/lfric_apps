! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module finalise_init_parcel_mod

implicit none

contains

! This routine does a couple of safety-checks on the
! initiating parcel properties (e.g. avoid negative q).
subroutine finalise_init_parcel( n_points, nc, index_ic,                       &
                                 q_vap_k,                                      &
                                 par_mean, par_core )

use comorph_constants_mod, only: real_cvprec, zero, one, l_par_core,           &
                                 max_qpert, par_gen_core_fac
use fields_type_mod, only: n_fields, i_q_vap

implicit none

! Total number of points in the par_gen arrays
integer, intent(in) :: n_points

! Points where any initiation mass-sources occur
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Grid-mean water-vapour mixing-ratio from level k
real(kind=real_cvprec), intent(in) :: q_vap_k(n_points)

! Parcel mean and core properties averaged over regions
real(kind=real_cvprec), intent(in out) :: par_mean                             &
                                          ( n_points, n_fields )
real(kind=real_cvprec), intent(in out) :: par_core                             &
                                          ( n_points, n_fields )

! Loop counters
integer :: ic, ic2


! Safety-check; don't allow initiating parcel q_vap
! to exceed the source-layer q_vap by more than a certain
! ratio, otherwise we will get too strongly CFL-limited
! by the mass-flux restriction to avoid creating negative q
do ic2 = 1, nc
  ic = index_ic(ic2)
  par_mean(ic,i_q_vap) = min( max( par_mean(ic,i_q_vap), zero ),               &
                              (one + max_qpert) * q_vap_k(ic) )
end do
if ( l_par_core ) then
  do ic2 = 1, nc
    ic = index_ic(ic2)
    par_core(ic,i_q_vap) = min( max( par_core(ic,i_q_vap), zero ),             &
             (one + max_qpert*par_gen_core_fac) * q_vap_k(ic) )
  end do
end if


return
end subroutine finalise_init_parcel


end module finalise_init_parcel_mod
