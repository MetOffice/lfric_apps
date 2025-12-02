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
! earlier, in set_par_winds.
! This routine also applies the CFL limit to the initiating
! mass-flux from each level / region.
subroutine add_region_parcel( n_points, nc, index_ic,                          &
                              init_mass, fields_par,                           &
                              pert_tl, pert_qt,                                &
                              max_ent_frac, l_within_bl,                       &
                              layer_mass, frac_r_k,                            &
                              par_massflux_d,                                  &
                              par_mean, par_core )

use comorph_constants_mod, only: real_cvprec, par_gen_core_fac, l_par_core,    &
                     comorph_timestep, i_cfl_check,                            &
                     i_cfl_check_local, i_cfl_check_hybrid
use fields_type_mod, only: n_fields, i_temperature, i_q_vap,                   &
                           i_qc_first

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

! Max fraction of mass in region which can initiate (CFL limit)
real(kind=real_cvprec), intent(in) :: max_ent_frac

! Flag for whether each point is within the boundary-layer
logical, intent(in) :: l_within_bl(n_points)

! Total dry-mass on level k
real(kind=real_cvprec), intent(in) :: layer_mass(n_points)

! Fraction occupied by current region
real(kind=real_cvprec), intent(in) :: frac_r_k(n_points)

! Initiating mass-flux summed over sub-grid regions
real(kind=real_cvprec), intent(in out) :: par_massflux_d(n_points)
! Parcel mean and core properties averaged over regions
real(kind=real_cvprec), intent(in out) :: par_mean                             &
                                         ( n_points, n_fields )
real(kind=real_cvprec), intent(in out) :: par_core                             &
                                         ( n_points, n_fields )

! Loop counters
integer :: ic, ic2, i_field


! First, check that the initiating mass doesn't exceed the
! CFL limit for numerical stability...
select case(i_cfl_check)
case (i_cfl_check_local)
  ! Always impose vertically local limit
  do ic2 = 1, nc
    ic = index_ic(ic2)
    init_mass(ic2) = min( init_mass(ic2),                                      &
                          max_ent_frac * layer_mass(ic)                        &
                          * frac_r_k(ic) / comorph_timestep )
  end do
case (i_cfl_check_hybrid)
  ! Only apply vertically local limit above the BL-top
  do ic2 = 1, nc
    ic = index_ic(ic2)
    if ( .not. l_within_bl(ic) ) then
      init_mass(ic2) = min( init_mass(ic2),                                    &
                            max_ent_frac * layer_mass(ic)                      &
                            * frac_r_k(ic) / comorph_timestep )
    end if
  end do
end select


! Add contribution to total initiation mass-source
! summed over all regions
do ic2 = 1, nc
  ic = index_ic(ic2)
  par_massflux_d(ic) = par_massflux_d(ic) + init_mass(ic2)
end do

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
