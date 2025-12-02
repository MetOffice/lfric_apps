! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_qss_forcing_init_mod

implicit none

contains

! Subroutine to calculate the forcing of supersaturation
! qss = q_vap - q_sat(T) by the initiating convective
! mass-sources.
! This is used to compute an implicit solution for the
! initiating mass-fluxes out of liquid cloud, to avoid
! noisy behaviour where the heating / drying due to the
! initiating convection causes the cloud-scheme to fully
! remove the cloud so that convection shuts off at the
! next timestep.
! This routine calculates the contribution to the qss forcing
! used in the implicit solve, for either updrafts
! or downdrafts (it must be called twice to sum the
! contributions from updrafts and downdrafts)
subroutine calc_qss_forcing_init( n_points, n_points_super,                    &
                                  nc, index_ic,                                &
                                  init_mass, fields_par,                       &
                                  pert_tl, pert_qt,                            &
                                  grid_k, grid_kpdk,                           &
                                  fields_k, fields_kpdk,                       &
                                  qss_forc )

use comorph_constants_mod, only: real_cvprec
use fields_type_mod, only: n_fields, i_temperature, i_q_vap,                   &
                           i_qc_first, i_qc_last, i_q_cl
use grid_type_mod, only: n_grid, i_pressure

use dry_adiabat_mod, only: dry_adiabat
use set_cp_tot_mod, only: set_cp_tot
use lat_heat_mod, only: lat_heat_incr, i_phase_change_evp
use set_qsat_mod, only: set_qsat_liq
use set_dqsatdt_mod, only:  set_dqsatdt_liq

implicit none

! Total number of points in the par_gen arrays
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Points where initiation mass-sources occur
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Initiation mass-sources from current region
real(kind=real_cvprec), intent(in) :: init_mass(nc)

! Unperturbed initiating parcel properties
real(kind=real_cvprec), intent(in) :: fields_par                               &
                            ( n_points, i_temperature:n_fields )

! Parcel perturbations to Tl and qt
real(kind=real_cvprec), intent(in) :: pert_tl(nc)
real(kind=real_cvprec), intent(in) :: pert_qt(nc)

! Height and pressure from level k and the next full level
real(kind=real_cvprec), intent(in) :: grid_k                                   &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kpdk                                &
                                      ( n_points_super, n_grid )

! Grid-mean fields from level k and the next full level
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                    ( n_points_super, n_fields )
real(kind=real_cvprec), intent(in) :: fields_kpdk                              &
                                    ( n_points_super, n_fields )

! Forcing of supersaturation qss = q_vap - q_sat by initiation
real(kind=real_cvprec), intent(in out) :: qss_forc(n_points)


! Compressed pressure
real(kind=real_cvprec) :: pressure_k(nc)
real(kind=real_cvprec) :: pressure_kpdk(nc)

! Compressed copy of grid-mean fields from next full level
real(kind=real_cvprec) :: fields_kpdk_cmpr                                     &
                          ( nc, i_temperature:i_qc_last )

! Array to store Tl, qw of the initiating parcel
real(kind=real_cvprec) :: fields_init                                          &
                          ( nc, i_temperature:i_q_vap )

! Saturation water-vapour mixing-ratio w.r.t. liquid at the
! grid-mean temperature, and its gradient with temperature
real(kind=real_cvprec) :: qsat_liq(nc)
real(kind=real_cvprec) :: dqsatdT_liq(nc)

! Compressed grid-mean temperature, for qsat call
real(kind=real_cvprec) :: temperature_k(nc)

! Exner ratio for pressure change
real(kind=real_cvprec) :: exner_ratio(nc)
! Total heat capacity for phase-change
real(kind=real_cvprec) :: cp_tot(nc)

! Convective forcing of qw=q_vap+q_cl, Tl, and supersaturation
real(kind=real_cvprec) :: qw_forc, tl_forc

! Loop counters
integer :: ic, ic2, i_field


! Compress pressures from current and next full levels
do ic2 = 1, nc
  ic = index_ic(ic2)
  pressure_k(ic2)    = grid_k(ic,i_pressure)
  pressure_kpdk(ic2) = grid_kpdk(ic,i_pressure)
end do

! Compress grid-mean fields from next level interface
do i_field = i_temperature, i_qc_last
  do ic2 = 1, nc
    ic = index_ic(ic2)
    fields_kpdk_cmpr(ic2,i_field) = fields_kpdk(ic,i_field)
  end do
end do

! Change temperature to what it would be at level k
call dry_adiabat( nc, nc,                                                      &
                  pressure_kpdk, pressure_k,                                   &
                  fields_kpdk_cmpr(:,i_q_vap),                                 &
                  fields_kpdk_cmpr(:,i_qc_first:i_qc_last),                    &
                  exner_ratio )
do ic2 = 1, nc
  fields_kpdk_cmpr(ic2,i_temperature)                                          &
    = fields_kpdk_cmpr(ic2,i_temperature) * exner_ratio(ic2)
end do

! Convert T,q to Tl,qw (i.e. evaporate the liquid cloud)
call set_cp_tot( nc, nc, fields_kpdk_cmpr(:,i_q_vap),                          &
                 fields_kpdk_cmpr(:,i_qc_first:i_qc_last),                     &
                 cp_tot )
! Add latent heating or cooling
call lat_heat_incr( nc, nc, i_phase_change_evp, cp_tot,                        &
                    fields_kpdk_cmpr(:,i_temperature),                         &
                    dq=fields_kpdk_cmpr(:,i_q_cl) )
! Update water-vapour
do ic2 = 1, nc
  fields_kpdk_cmpr(ic2,i_q_vap) = fields_kpdk_cmpr(ic2,i_q_vap)                &
                                + fields_kpdk_cmpr(ic2,i_q_cl)
end do
! fields_kpdk_cmpr now stores grid-mean qw and Tl from k+dk


! Copy initiating parcel T,q and convert them to Tl,qw
do ic2 = 1, nc
  fields_init(ic2,i_temperature) = fields_par(ic2,i_temperature)               &
                                 + pert_tl(ic2)
  fields_init(ic2,i_q_vap) = fields_par(ic2,i_q_vap)                           &
                           + pert_qt(ic2)
end do
call set_cp_tot( nc, n_points, fields_init(:,i_q_vap),                         &
                 fields_par(:,i_qc_first:i_qc_last), cp_tot )
call lat_heat_incr( nc, nc, i_phase_change_evp, cp_tot,                        &
                    fields_init(:,i_temperature),                              &
                    dq=fields_par(:,i_q_cl) )
do ic2 = 1, nc
  fields_init(ic2,i_q_vap) = fields_init(ic2,i_q_vap)                          &
                           + fields_par(ic2,i_q_cl)
end do
! fields_init now stores parcel qw and Tl


! Calculate qsat_liq and dqsat/dT for grid-mean fields at k
do ic2 = 1, nc
  ic = index_ic(ic2)
  temperature_k(ic2) = fields_k(ic,i_temperature)
end do
call set_qsat_liq( nc, temperature_k, pressure_k, qsat_liq )
call set_dqsatdt_liq( nc, temperature_k, qsat_liq, dqsatdT_liq )

do ic2 = 1, nc
  ic = index_ic(ic2)

  ! Estimate the convective tendency of vapour + liquid:
  ! Contribution from removal of initiating air and
  ! replacing it with air from the next model-level
  ! q_forc = q_kpdk - q_init
  qw_forc = fields_kpdk_cmpr(ic2,i_q_vap)                                      &
          - fields_init(ic2,i_q_vap)

  ! Estimate the convective tendency of Tl:
  ! (same method as for qw)
  tl_forc = fields_kpdk_cmpr(ic2,i_temperature)                                &
          - fields_init(ic2,i_temperature)

  ! Estimate contribution to forcing of supersaturation qss
  qss_forc(ic) = qss_forc(ic) + init_mass(ic2)                                 &
                     * ( qw_forc - dqsatdT_liq(ic2) * tl_forc )

end do


return
end subroutine calc_qss_forcing_init


end module calc_qss_forcing_init_mod
