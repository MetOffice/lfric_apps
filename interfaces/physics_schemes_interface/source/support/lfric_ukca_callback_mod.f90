!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------

!> @brief Module for providing callback routines to service UKCA calls to
!   parent-specific procedures.
!
!   The routine(s) listed below are provided for passing to UKCA as
!   callback routines via the ukca_setup argument list.
!   Other callback routines should be added to this module as needed.
!
!   bl_tracer_mix - calls UM TR_MIX routine to do boundary layer mixing
!                   for a tracer
! ----------------------------------------------------------------------

module lfric_ukca_callback_mod

implicit none

private

! public procedures
public bl_tracer_mix

contains

! ----------------------------------------------------------------------
subroutine bl_tracer_mix(row_length, rows, bl_levels,                          &
                         r_theta_levels, r_rho_levels,                         &
                         nlev_ent_tr_mix,                                      &
                         kent, kent_dsc, surf_em, zhnl, zhsc,                  &
                         we_lim, t_frac, zrzi,                                 &
                         we_lim_dsc, t_frac_dsc, zrzi_dsc,                     &
                         z_uv, rhokh_rdz, dtrdz, field)
! ----------------------------------------------------------------------
! Description:
!   UKCA-compatible wrapper for calling tr_mix routine to do boundary
!   layer mixing for a tracer after applying emission(s).
! ----------------------------------------------------------------------

use constants_mod, only : r_um, i_um
! UM modules
use tr_mix_mod, only: tr_mix
use bl_option_mod, only: alpha_cd
use atm_fields_bounds_mod, only: pdims

implicit none

integer(i_um), intent(in) :: row_length
integer(i_um), intent(in) :: rows
integer(i_um), intent(in) :: bl_levels
integer(i_um), intent(in) :: nlev_ent_tr_mix

! Inputs and outputs to the interface 2D i by j
real(r_um), intent(in) :: r_theta_levels(1:row_length,1:rows,0:bl_levels)
  ! height of theta levels from centre of earth
real(r_um), intent(in) :: r_rho_levels(1:row_length,1:rows,bl_levels)
  ! height of rho levels from centre of earth
integer, intent(in) :: kent(row_length, rows)
  ! grid level of surface mixed layer inversion
integer, intent(in) :: kent_dsc(row_length, rows)
  ! grid level of decoupled stratocumulus inversion
real(r_um), intent(in) :: surf_em(row_length, rows)
  ! emission flux into surface level (kg/m^2/s)
real(r_um), intent(in) :: zhnl(row_length, rows)
  ! atmosphere_boundary layer thickness (m) {stashcode:00025}
real(r_um), intent(in) :: zhsc(row_length, rows)
  ! height of top of decoupled stratocumulus layer (m) {stashcode:03073}
real(r_um), intent(in) :: we_lim(row_length, rows, nlev_ent_tr_mix)
  ! density * entrainment rate implied by placing of subsidence at surface mixed
  ! layer inversion (kg/m^2/s) {stashcode:03066}
real(r_um), intent(in) :: t_frac(row_length, rows, nlev_ent_tr_mix)
  ! fraction of timestep surface mixed layer inversion is above level
  ! {stashcode:03067}
real(r_um), intent(in) :: zrzi(row_length, rows, nlev_ent_tr_mix)
  ! level height as fraction of surface mixed layer inversion height above ml
  ! base {stashcode:03068}
real(r_um), intent(in) :: we_lim_dsc(row_length, rows, nlev_ent_tr_mix)
  ! density * entrainment rate implied by placing of subsidence at decoupled
  ! stratocumulus inversion (kg/m^2/s) {stashcode:03070}
real(r_um), intent(in) :: t_frac_dsc(row_length, rows, nlev_ent_tr_mix)
  ! fraction of timestep decoupled stratocumulus inversion is above level
  ! {stashcode:03071}
real(r_um), intent(in) :: zrzi_dsc(row_length, rows, nlev_ent_tr_mix)
  ! level height as fraction of decoupled stratocumulus inversion height above
  ! dsc ml base {stashcode:03072}
real(r_um), intent(in) :: z_uv(row_length, rows, bl_levels)
  ! height at rho levels (m)
real(r_um), intent(in) :: rhokh_rdz(row_length, rows, 2:bl_levels)
  ! mixing coefficient above surface:
  ! (scalar eddy diffusivity * density) / dz (kg/m^2/s) {stashcode:03060}
real(r_um), intent(in) :: dtrdz(row_length, rows, bl_levels)
  ! dt/(density*radius*radius*dz) for scalar flux divergence (s/kg)
  ! {stashcode:03064}
real(r_um), intent(in out) :: field(row_length, rows, bl_levels)
  ! tracer mixing ratio (kg/kg)

! Convert to 1D, i, for tr_mix call
real(r_um) :: r_theta_levels_1D(1:row_length,0:bl_levels)
  ! height of theta levels from centre of earth
real(r_um) :: r_rho_levels_1D(1:row_length,bl_levels)
  ! height of rho levels from centre of earth
integer :: kent_1D(row_length)
  ! grid level of surface mixed layer inversion
integer :: kent_dsc_1D(row_length)
  ! grid level of decoupled stratocumulus inversion
real(r_um) :: surf_em_1D(row_length)
  ! emission flux into surface level (kg/m^2/s)
real(r_um) :: zhnl_1D(row_length)
  ! atmosphere_boundary layer thickness (m) {stashcode:00025}
real(r_um) :: zhsc_1D(row_length)
  ! height of top of decoupled stratocumulus layer (m) {stashcode:03073}
real(r_um) :: we_lim_1D(row_length, nlev_ent_tr_mix)
  ! density * entrainment rate implied by placing of subsidence at surface mixed
  ! layer inversion (kg/m^2/s) {stashcode:03066}
real(r_um) :: t_frac_1D(row_length, nlev_ent_tr_mix)
  ! fraction of timestep surface mixed layer inversion is above level
  ! {stashcode:03067}
real(r_um) :: zrzi_1D(row_length, nlev_ent_tr_mix)
  ! level height as fraction of surface mixed layer inversion height above ml
  ! base {stashcode:03068}
real(r_um) :: we_lim_dsc_1D(row_length, nlev_ent_tr_mix)
  ! density * entrainment rate implied by placing of subsidence at decoupled
  ! stratocumulus inversion (kg/m^2/s) {stashcode:03070}
real(r_um) :: t_frac_dsc_1D(row_length, nlev_ent_tr_mix)
  ! fraction of timestep decoupled stratocumulus inversion is above level
  ! {stashcode:03071}
real(r_um) :: zrzi_dsc_1D(row_length, nlev_ent_tr_mix)
  ! level height as fraction of decoupled stratocumulus inversion height above
  ! dsc ml base {stashcode:03072}
real(r_um) :: z_uv_1D(row_length, bl_levels)
  ! height at rho levels (m)
real(r_um) :: rhokh_rdz_1D(row_length, 2:bl_levels)
  ! mixing coefficient above surface:
  ! (scalar eddy diffusivity * density) / dz (kg/m^2/s) {stashcode:03060}
real(r_um) :: dtrdz_1D(row_length, bl_levels)
  ! dt/(density*radius*radius*dz) for scalar flux divergence (s/kg)
  ! {stashcode:03064}
real(r_um) :: field_1D(row_length, bl_levels)
  ! tracer mixing ratio (kg/kg)

! local variables
real(r_um) :: rhokh_1(row_length)                  ! surface exchange coeff.
real(r_um) :: res_factor(row_length)               ! dry deposition coeff.
real(r_um) :: f_field(row_length, rows, bl_levels) ! tracer flux from tr_mix
real(r_um) :: f_field_1D(row_length, bl_levels)    ! tracer flux from tr_mix
real(r_um) :: surf_dep_flux(row_length, rows)      ! surf. deposition flux from
                                                   ! tr_mix
real(r_um) :: surf_dep_flux_1D(row_length)         ! surf. deposition flux from
                                                   ! tr_mix
integer(i_um) :: i, k   ! Local loop indexes

! rhokh_1(:) = 0.0_r_um
! res_factor(:) = 0.0_r_um

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP  private( i, k )
!$OMP  do SCHEDULE(STATIC)
do i = 1, row_length
  ! Local inits
  rhokh_1(i) = 0.0_r_um
  res_factor(i) = 0.0_r_um

! Copies to different array sizes
  r_theta_levels_1D(i,0) = r_theta_levels(i,1,0)
  kent_1D(i) = kent(i,1)
  kent_dsc_1D(i) = kent_dsc(i,1)
  surf_em_1D(i) = surf_em(i,1)
  zhnl_1D(i) = zhnl(i,1)
  zhsc_1D(i) = zhsc(i,1)
end do
!$OMP end do
!$OMP  do SCHEDULE(STATIC)
do k = 1, bl_levels
  do i = 1, row_length
    r_theta_levels_1D(i,k) = r_theta_levels(i,1,k)
    r_rho_levels_1D(i,k) = r_rho_levels(i,1,k)
    z_uv_1D(i,k) = z_uv(i,1,k)
    dtrdz_1D(i,k) = dtrdz(i,1,k)
    field_1D(i,k) = field(i,1,k)
  end do
end do
!$OMP end do
!$OMP  do SCHEDULE(STATIC)
do k = 2, bl_levels
  do i = 1, row_length
    rhokh_rdz_1D(i,k) = rhokh_rdz(i,1,k)
  end do
end do
!$OMP end do
!$OMP  do SCHEDULE(STATIC)
do k = 1, nlev_ent_tr_mix
  do i = 1, row_length
    we_lim_1D(i,k) = we_lim(i,1,k)
    t_frac_1D(i,k) = t_frac(i,1,k)
    zrzi_1D(i,k) = zrzi(i,1,k)
    we_lim_dsc_1D(i,k) = we_lim_dsc(i,1,k)
    t_frac_dsc_1D(i,k) = t_frac_dsc(i,1,k)
    zrzi_dsc_1D(i,k) = zrzi_dsc(i,1,k)
  end do
end do
!$OMP end do
!$OMP end PARALLEL

call tr_mix( r_theta_levels_1D, r_rho_levels_1D, pdims, bl_levels,             &
             alpha_cd, rhokh_rdz_1D, rhokh_1, dtrdz_1D, surf_em_1D,            &
             res_factor, kent_1D, we_lim_1D, t_frac_1D, zrzi_1D,               &
             kent_dsc_1D, we_lim_dsc_1D, t_frac_dsc_1D, zrzi_dsc_1D,           &
             zhnl_1D, zhsc_1D, z_uv_1D,                                        &
             ! Output fields
             field_1D, f_field_1D, surf_dep_flux_1D)

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP  private( i, k )
!$OMP  do SCHEDULE(STATIC)
do k = 1, bl_levels
  do i = 1, row_length
    field(i,1,k) = field_1D(i,k)
    ! f_field(i,1,k) = f_field_1D(i,k)
    ! surf_dep_flux(i,1) = surf_dep_flux_1D(i)
  end do
end do
!$OMP end do
!$OMP end PARALLEL

return
end subroutine bl_tracer_mix

end module lfric_ukca_callback_mod
