! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module normalise_init_parcel_mod

implicit none

contains

! Subroutine to normalise the dry-mass-flux weighted initiating
! parcel properties, after computing mass-weighted contributions
! from each sub-grid region.  Involves dividing by the
! initiating mass-flux.
! Note this normalisation is not needed for the winds or
! tracer fields, since they are not averaged over sub-grid
! regions (they are assumed equal in all regions).
! This routine also does a couple of safety-checks on the
! initiating parcel properties (e.g. avoid negative q)
subroutine normalise_init_parcel( n_points, nc, index_ic,                      &
                                  q_vap_k, par_massflux_d,                     &
                                  par_mean, par_core )

use comorph_constants_mod, only: real_cvprec, one, l_par_core,                 &
                     max_qpert, par_gen_core_fac
use fields_type_mod, only: n_fields, i_q_vap, i_temperature

implicit none

! Total number of points in the par_gen arrays
integer, intent(in) :: n_points

! Points where any initiation mass-sources occur
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Grid-mean water-vapour mixing-ratio from level k
real(kind=real_cvprec), intent(in) :: q_vap_k(n_points)
! Initiating mass-flux summed over sub-grid regions
real(kind=real_cvprec), intent(in) :: par_massflux_d(n_points)

! Parcel mean and core properties averaged over regions
real(kind=real_cvprec), intent(in out) :: par_mean                             &
                                         ( n_points, n_fields )
real(kind=real_cvprec), intent(in out) :: par_core                             &
                                         ( n_points, n_fields )

! Loop counters
integer :: ic, ic2, i_field


! Normalise the mean initiating parcel properties over
! all the regions
do i_field = i_temperature, n_fields
  do ic2 = 1, nc
    ic = index_ic(ic2)
    par_mean(ic,i_field) = par_mean(ic,i_field)                                &
                         / par_massflux_d(ic)
  end do
end do

! Same for parcel core if used
if ( l_par_core ) then
  do i_field = i_temperature, n_fields
    do ic2 = 1, nc
      ic = index_ic(ic2)
      par_core(ic,i_field) = par_core(ic,i_field)                              &
                           / par_massflux_d(ic)
    end do
  end do
end if

! Safety-check; don't allow initiating parcel q_vap
! to exceed the source-layer q_vap by more than a certain
! ratio, otherwise we will get too strongly CFL-limited
! by the mass-flux restriction to avoid creating negative q
do ic2 = 1, nc
  ic = index_ic(ic2)
  par_mean(ic,i_q_vap) = min( par_mean(ic,i_q_vap),                            &
                              (one + max_qpert) * q_vap_k(ic) )
end do
if ( l_par_core ) then
  ! Apply scaled-up version of the same limit to parcel core
  do ic2 = 1, nc
    ic = index_ic(ic2)
    par_core(ic,i_q_vap) = min( par_core(ic,i_q_vap),                          &
             (one + max_qpert*par_gen_core_fac) * q_vap_k(ic) )
  end do
end if


return
end subroutine normalise_init_parcel


end module normalise_init_parcel_mod
