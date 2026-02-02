! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module add_region_parcel_mod

implicit none

contains

! Subroutine to add initiation sources from the current
! sub-grid region to the initiating parcel properties,
! for either updraft or downdraft.
! This routine only handles parcel properties which differ
! between the sub-grid regions (T,q,qc,cf); winds and tracers
! are assumed equal in all regions and so are calculated
! earlier, in set_par_fields.
! This routine also applies the CFL limit to the initiating
! mass-flux from each level / region.
subroutine add_region_parcel( n_points, nc, index_ic,                          &
                              init_mass, fields_par,                           &
                              pert_tl, pert_qt,                                &
                              par_super, par_mean, par_core )

use comorph_constants_mod, only: real_cvprec, par_gen_core_fac, l_par_core
use fields_type_mod, only: n_fields, i_temperature, i_q_vap,                   &
                           i_qc_first, i_qc_last
use parcel_type_mod, only: n_par, i_massflux_d
use calc_virt_temp_mod, only: calc_virt_temp

implicit none

! Total number of points in the par_gen arrays
integer, intent(in) :: n_points

! Points where initiation mass-sources occur in current region
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Initiating mass-flux from current region
real(kind=real_cvprec), intent(in out) :: init_mass(nc)

! Unperturbed initiating parcel properties
real(kind=real_cvprec), intent(in) :: fields_par                               &
                            ( n_points, i_temperature:n_fields )
! Perturbations applied to initiating parcel Tl and qt
real(kind=real_cvprec), intent(in) :: pert_tl(nc)
real(kind=real_cvprec), intent(in) :: pert_qt(nc)

! Parcel super-array for initiating mass-flux summed over sub-grid regions
real(kind=real_cvprec), intent(in out) :: par_super                            &
                                          ( n_points, n_par )
! Parcel mean and core properties averaged over regions
real(kind=real_cvprec), intent(in out) :: par_mean                             &
                                          ( n_points, n_fields )
real(kind=real_cvprec), intent(in out) :: par_core                             &
                                          ( n_points, n_fields )

! Virtual temperature of unperturbed parcel fields
real(kind=real_cvprec) :: virt_temp(nc)

! Loop counters
integer :: ic, ic2, i_field


! Add contribution to total initiation mass-source
! summed over all regions
do ic2 = 1, nc
  ic = index_ic(ic2)
  par_super(ic,i_massflux_d) = par_super(ic,i_massflux_d) + init_mass(ic2)
end do

! Calculate virtual temperature of the unperturbed parcel fields
call calc_virt_temp( nc, n_points,                                             &
                     fields_par(:,i_temperature),                              &
                     fields_par(:,i_q_vap),                                    &
                     fields_par(:,i_qc_first:i_qc_last),                       &
                     virt_temp )

! Store mass-flux-weighted contribution in the
! parcel mean-fields array
do ic2 = 1, nc
  ! Perturbations applied to T,q
  ic = index_ic(ic2)
  par_mean(ic,i_temperature) = par_mean(ic,i_temperature)                      &
      + ( fields_par(ic2,i_temperature) + pert_tl(ic2) )                       &
        * init_mass(ic2)
  par_mean(ic,i_q_vap) = par_mean(ic,i_q_vap)                                  &
      + ( fields_par(ic2,i_q_vap) + pert_qt(ic2) )                             &
        * init_mass(ic2)
end do
do i_field = i_qc_first, n_fields
  ! Other fields unperturbed
  do ic2 = 1, nc
    ic = index_ic(ic2)
    par_mean(ic,i_field) = par_mean(ic,i_field)                                &
      + fields_par(ic2,i_field) * init_mass(ic2)
  end do
end do

! Add contributions to parcel core if used, with
! the perturbations scaled up by par_gen_core_fac
if ( l_par_core ) then
  do ic2 = 1, nc
    ! Perturbations applied to T,q
    ic = index_ic(ic2)
    par_core(ic,i_temperature) = par_core(ic,i_temperature)                    &
      + ( fields_par(ic2,i_temperature)                                        &
        + pert_tl(ic2) * par_gen_core_fac )                                    &
        * init_mass(ic2)
    par_core(ic,i_q_vap) = par_core(ic,i_q_vap)                                &
      + ( fields_par(ic2,i_q_vap)                                              &
        + pert_qt(ic2) * par_gen_core_fac )                                    &
        * init_mass(ic2)
  end do
  do i_field = i_qc_first, n_fields
    ! Other fields unperturbed
    do ic2 = 1, nc
      ic = index_ic(ic2)
      par_core(ic,i_field) = par_core(ic,i_field)                              &
        + fields_par(ic2,i_field) * init_mass(ic2)
    end do
  end do
end if


return
end subroutine add_region_parcel


end module add_region_parcel_mod
