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
subroutine finalise_init_parcel( n_points, n_points_super, n_fields_tot,       &
                                 nc, index_ic, l_tracer,                       &
                                 fields_k,                                     &
                                 par_mean, par_core )

use comorph_constants_mod, only: real_cvprec, zero, one,                       &
                                 l_par_core, n_tracers,                        &
                                 max_qpert, par_gen_core_fac
use fields_type_mod, only: i_q_vap

use tracer_source_mod, only: tracer_gen

implicit none

! Total number of points in the par_gen arrays
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Total number of fields in the input super-arrays,
! including tracers if they are used
integer, intent(in) :: n_fields_tot

! Points where any initiation mass-sources occur
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Flag for whether tracers are in use
logical, intent(in) :: l_tracer

! Primary model-fields at level k
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                      ( n_points_super, n_fields_tot )

! Parcel mean and core properties averaged over regions
real(kind=real_cvprec), intent(in out) :: par_mean                             &
                                          ( n_points, n_fields_tot )
real(kind=real_cvprec), intent(in out) :: par_core                             &
                                          ( n_points, n_fields_tot )

! Loop counters
integer :: ic, ic2


! Safety-check; don't allow initiating parcel q_vap
! to exceed the source-layer q_vap by more than a certain
! ratio, otherwise we will get too strongly CFL-limited
! by the mass-flux restriction to avoid creating negative q
do ic2 = 1, nc
  ic = index_ic(ic2)
  par_mean(ic,i_q_vap) = min( max( par_mean(ic,i_q_vap), zero ),               &
                              (one + max_qpert) * fields_k(ic,i_q_vap) )
end do
if ( l_par_core ) then
  do ic2 = 1, nc
    ic = index_ic(ic2)
    par_core(ic,i_q_vap) = min( max( par_core(ic,i_q_vap), zero ),             &
             (one + max_qpert*par_gen_core_fac) * fields_k(ic,i_q_vap) )
  end do
end if

if ( l_tracer .and. n_tracers > 0 ) then
  ! Set parcel initial tracer concentrations
  ! (done in a separate routine to handle any species-dependence)
  call tracer_gen( n_points, n_points_super, n_fields_tot,                     &
                   fields_k, par_mean, par_core )
end if


return
end subroutine finalise_init_parcel


end module finalise_init_parcel_mod
