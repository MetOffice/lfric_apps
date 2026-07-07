! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module cor_init_mass_liq_1_mod

implicit none

contains

! Subroutine to reduce the initiation mass-sources from liquid-cloud,
! due to the reduction in liquid-cloud fraction over the timestep
! forced by the increments from removal of initiating mass
! (i.e. this is an approximate implicit correction).

subroutine cor_init_mass_liq_1( n_points, n_points_super,                      &
                                nc_up, index_ic_up, nc_dn, index_ic_dn,        &
                                frac_r_k, layer_mass_k, q_cl_loc_k,            &
                                cp_tot, L_con, dqsatdt, dqsatdt_t,             &
                                fields_par_up, fields_par_dn,                  &
                                pert_tl_up_t, pert_qt_up_t,                    &
                                pert_tl_dn_t, pert_qt_dn_t,                    &
                                grid_km1, grid_k, grid_kp1,                    &
                                fields_km1, fields_kp1,                        &
                                init_mass_up_t, init_mass_dn_t,                &
                                imp_coef )

use comorph_constants_mod, only: real_cvprec, zero, one, two, sqrt_min_float,  &
                                 n_updraft_types, n_dndraft_types,             &
                                 comorph_timestep
use fields_type_mod, only: n_fields, i_temperature
use grid_type_mod, only: n_grid

use calc_qss_forcing_init_mod, only: calc_qss_forcing_init

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Points where updrafts and downdrafts initiate in current region
integer, intent(in) :: nc_up
integer, intent(in) :: index_ic_up(n_points)
integer, intent(in) :: nc_dn
integer, intent(in) :: index_ic_dn(n_points)

! Area-fraction of current region (liquid-only or mixed-phase cloud)
real(kind=real_cvprec), intent(in) :: frac_r_k(n_points)

! Dry-mass per unit surface area on the current model-level
real(kind=real_cvprec), intent(in) :: layer_mass_k(n_points)

! Local liquid-cloud water content within the liquid-cloud
real(kind=real_cvprec), intent(in) :: q_cl_loc_k(n_points)

! Env k total heat capacity
real(kind=real_cvprec), intent(in) :: cp_tot(n_points)
! Env k latent heat of condensation
real(kind=real_cvprec), intent(in) :: L_con(n_points)
! dqsat/dT w.r.t. liquid water at level k
real(kind=real_cvprec), intent(in) :: dqsatdt(n_points)   ! Value at Tl
real(kind=real_cvprec), intent(in) :: dqsatdt_t(n_points) ! Value at T

! Unperturbed initiating parcel properties from current region,
! for updrafts and downdrafts
! (doesn't include winds, since they are assumed equal in all
!  the sub-grid regions)
real(kind=real_cvprec), intent(in) :: fields_par_up                            &
                                      ( n_points, i_temperature:n_fields )
real(kind=real_cvprec), intent(in) :: fields_par_dn                            &
                                      ( n_points, i_temperature:n_fields )

! Tl and qt perturbations for current region, for updrafts and downdrafts
! (can be set separately for each convection type)
real(kind=real_cvprec), intent(in) :: pert_tl_up_t                             &
                                      ( n_points, n_updraft_types )
real(kind=real_cvprec), intent(in) :: pert_qt_up_t                             &
                                      ( n_points, n_updraft_types )
real(kind=real_cvprec), intent(in) :: pert_tl_dn_t                             &
                                      ( n_points, n_dndraft_types )
real(kind=real_cvprec), intent(in) :: pert_qt_dn_t                             &
                                      ( n_points, n_dndraft_types )

! Height and pressure at k and neighbouring full levels
real(kind=real_cvprec), intent(in) :: grid_km1                                 &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_k                                   &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kp1                                 &
                                      ( n_points_super, n_grid )

! Primary model-fields at level k and neighbouring full levels
real(kind=real_cvprec), intent(in) :: fields_km1                               &
                                      ( n_points_super, n_fields )
real(kind=real_cvprec), intent(in) :: fields_kp1                               &
                                      ( n_points_super, n_fields )

! Initiating updraft and downdraft mass-flux for each convection type
real(kind=real_cvprec), intent(in out) :: init_mass_up_t                       &
                                          ( n_points, n_updraft_types )
real(kind=real_cvprec), intent(in out) :: init_mass_dn_t                       &
                                          ( n_points, n_dndraft_types )

! Rescaling of mass-sources so implicit w.r.t. cloud fraction
real(kind=real_cvprec), intent(out) :: imp_coef(n_points)

! Sum of updraft and downdraft initiation mass-sources
real(kind=real_cvprec) :: mass_forc(n_points)

! Sum of initiation forcing of qss = q_vap+q_cl - qsat(Tl)
real(kind=real_cvprec) :: qss_forc(n_points)

! Loop counters
integer :: ic, ic2, i_type


! If the input liquid-cloud fraction is f_exp, and explicit estimates of
! the initiating mass-fluxes are M_up_exp, M_dn_exp, we rescale them consistent
! with an implicit solution for the cloud-fraction f_imp:
!
! M_up_imp = M_up_exp f_imp/f_exp
! M_dn_imp = M_dn_exp f_imp/f_exp
!
! When initiating from large-scale liquid cloud (regions liq and mph),
! the main feedback limiting the mass-flux is that the convective heating
! and drying will cause the cloud-scheme to remove the cloud.
! Exactly how quickly will depend on the details of the cloud-scheme being
! used.  Here we estimate the convective forcing of supersaturation
! qss = qw - qsat
! by the initiating mass-source, and use this to parameterise the expected
! rate of removal of liquid cloud by the convection...
!
! Convection also removes liquid cloud directly due to the cloud-volume
! which is extracted to form the initiating parcel.  The implicitness
! rescaling also accounts for this effect.

! Initialise sums over updraft and downdraft sources
do ic = 1, n_points
  ! Sum of updraft and downdraft initiating mass-sources
  ! from the current region:
  mass_forc(ic) = zero
  ! Sum of forcing of qss by the initiating mass-sources
  ! from the current region:
  qss_forc(ic) = zero
end do

! Sum the contributions to the forcing of mass and
! qss = qw - q_sat
! from both updrafts and downdrafts initiating in the
! current region / model-level
if ( nc_up > 0 ) then
  do i_type = 1, n_updraft_types
    do ic2 = 1, nc_up
      ic = index_ic_up(ic2)
      mass_forc(ic) = mass_forc(ic) + init_mass_up_t(ic2,i_type)
    end do
    call calc_qss_forcing_init(                                                &
           n_points, n_points_super, nc_up, index_ic_up,                       &
           dqsatdt_t, init_mass_up_t(:,i_type), fields_par_up,                 &
           pert_tl_up_t(:,i_type), pert_qt_up_t(:,i_type),                     &
           grid_k, grid_kp1, fields_kp1,                                       &
           qss_forc )
  end do
end if
if ( nc_dn > 0 ) then
  do i_type = 1, n_dndraft_types
    do ic2 = 1, nc_dn
      ic = index_ic_dn(ic2)
      mass_forc(ic) = mass_forc(ic) + init_mass_dn_t(ic2,i_type)
    end do
    call calc_qss_forcing_init(                                                &
           n_points, n_points_super, nc_dn, index_ic_dn,                       &
           dqsatdt_t, init_mass_dn_t(:,i_type), fields_par_dn,                 &
           pert_tl_dn_t(:,i_type), pert_qt_dn_t(:,i_type),                     &
           grid_k, grid_km1, fields_km1,                                       &
           qss_forc )
  end do
end if

! Apply limit to use explicit solution in case where the initiation
! mass-sources actually act to increase the liquid cloud instead of
! decreasing it (i.e. only retain negative values of dqss/dt)
do ic = 1, n_points
  qss_forc(ic) = min( qss_forc(ic), zero )
end do

! Note: technically qss_forc now needs to be normalised by dividing it by
! layer_mass.  This is done in the formula for imp_coef below...

! Calculate the implicit rescaling factor
! imp_coef = f_imp/f_exp.
! Parameterise the forcing of liquid cloud-fraction as:
!
! (df/dt)_exp = -M_init / M_layr
!             + f_exp * dqss/dt
!                     / ( 2 (1 + dqs/dT Lc/cp) qc_incloud )
!
! (1st term from direct fractional cloud-volume removal,
!  2nd term from homogeneous forcing of remaining cloud)
!
! The implicit liquid cloud fraction is given by:
!
! f_imp = f_exp + (df/dt)_exp * (f_imp/f_exp) * dt
!
! => (f_imp/f_exp) ( 1 - ( (df/dt)_exp / f_exp ) * dt ) = 1
! => (f_imp/f_exp)
!       = 1 / ( 1 - ( (df/dt)_exp / f_exp ) * dt )
!
! Substituting our above expression for df/dt, we get:
do ic = 1, n_points
  imp_coef(ic) = one / ( one + (                                               &
    ! Term from volume-removal:
      mass_forc(ic) / max( frac_r_k(ic), sqrt_min_float )                      &
    ! Term from homogeneous forcing:
    - qss_forc(ic) / ( two * (one + dqsatdt(ic)*L_con(ic)/cp_tot(ic))          &
                           * max( q_cl_loc_k(ic), sqrt_min_float ) )           &
                               )                                               &
    ! All terms have scaling of dt/M_layr
      * comorph_timestep / layer_mass_k(ic) )
end do

! Modify the initiating mass-fluxes accordingly
do ic2 = 1, nc_up
  ic = index_ic_up(ic2)
  do i_type = 1, n_updraft_types
    init_mass_up_t(ic2,i_type) = init_mass_up_t(ic2,i_type) * imp_coef(ic)
  end do
end do
do ic2 = 1, nc_dn
  ic = index_ic_dn(ic2)
  do i_type = 1, n_dndraft_types
    init_mass_dn_t(ic2,i_type) = init_mass_dn_t(ic2,i_type) * imp_coef(ic)
  end do
end do


return
end subroutine cor_init_mass_liq_1

end module cor_init_mass_liq_1_mod
