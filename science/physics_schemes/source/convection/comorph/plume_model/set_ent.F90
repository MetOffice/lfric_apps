! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_ent_mod

implicit none

contains

! Subroutine sets the entrained mass from the current layer,
! and sets the properties of the entrained air
subroutine set_ent( n_points, n_fields_tot, max_points,                        &
                    max_ent_frac,                                              &
                    par_conv_mean_fields, ent_fields,                          &
                    grid_prev_super, grid_next_super,                          &
                    par_conv_super,                                            &
                    l_within_bl, core_mean_ratio,                              &
                    layer_mass_step, sum_massflux,                             &
                    ent_mass_d, core_ent_ratio )

use comorph_constants_mod, only: real_cvprec, min_float, zero, one,            &
                                 ent_coef, comorph_timestep,                   &
                                 core_ent_fac, l_core_ent_cmr,                 &
                                 i_cfl_local, i_cfl_local_all,                 &
                                 i_cfl_local_nobl
use fields_type_mod, only: i_temperature, i_q_vap
use grid_type_mod, only: n_grid, i_height, i_pressure
use env_half_mod, only: n_env_half
use parcel_type_mod, only: n_par, i_massflux_d, i_radius
use calc_rho_dry_mod, only: calc_rho_dry

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of fields (including tracres if applicable)
integer, intent(in) :: n_fields_tot

! Number of points in the compressed super-arrays
! (to avoid repeatedly deallocating and reallocating these,
!  they are dimensioned with the biggest size they will need,
!  which will often be bigger than the number of points here)
integer, intent(in) :: max_points

! Maximum allowed entrained fraction of layer mass for current draft
real(kind=real_cvprec), intent(in) :: max_ent_frac

! Super-array containing the parcel mean primary fields
real(kind=real_cvprec), intent(in) :: par_conv_mean_fields                     &
                                      ( max_points, n_fields_tot )
! Entrained primary field values
real(kind=real_cvprec), intent(in) :: ent_fields                               &
                                      ( n_points, n_fields_tot )

! Height and pressure of parcel before and after the level-step
real(kind=real_cvprec), intent(in) :: grid_prev_super ( max_points, n_grid )
real(kind=real_cvprec), intent(in) :: grid_next_super ( max_points, n_grid )

! Super-arrays containing parcel mass-flux, radius, etc
real(kind=real_cvprec), intent(in) :: par_conv_super ( max_points, n_par )

! Flag for whether each point is within the boundary-layer
logical, intent(in) :: l_within_bl(n_points)

! Ratio of parcel core buoyancy over parcel mean buoyancy
real(kind=real_cvprec), intent(in) :: core_mean_ratio(n_points)

! Dry-mass on the current half-level-step, per m2 at surface.
real(kind=real_cvprec), intent(in) :: layer_mass_step(n_points)

! Sum of previous level-interface mass-fluxes over all convection types/layers
real(kind=real_cvprec), intent(in) :: sum_massflux(n_points)

! Rate of entrainment of dry-mass from current level / kg m-2 s-1
real(kind=real_cvprec), intent(out) :: ent_mass_d(n_points)

! Weight used to calculate properties of air entrained into the core
real(kind=real_cvprec), intent(out) :: core_ent_ratio(n_points)


! Dry-density of the entrained air and the parcel
real(kind=real_cvprec) :: ent_rho_dry(n_points)
real(kind=real_cvprec) :: par_rho_dry(n_points)

! Maximum allowed entrainment for numerical stability
real(kind=real_cvprec) :: max_ent(n_points)

! Loop counters
integer :: ic


!------------------------------------------------------------------------------
! 1) Calculate 1/R "mixing" entrainment rate
!------------------------------------------------------------------------------

! Set the fractional volume entrainment rate in m-1
do ic = 1, n_points
  ent_mass_d(ic) = ent_coef / par_conv_super(ic,i_radius)
end do

! Scale by env-relative distance travelled by parcel to get
! volume fraction entrained over the layer
do ic = 1, n_points
  ent_mass_d(ic) = ent_mass_d(ic)                                              &
          * abs( grid_next_super(ic,i_height) - grid_prev_super(ic,i_height) )
  ! Note: for downdrafts, next_height < prev_height, so need
  ! ABS call to ensure dz is positive.
end do

! Note: here we could account for horizontal distance travelled
! due to slant-wise ascent as well, but not implemented yet.

! Compute dry-density of the entrained air and the parcel
call calc_rho_dry( n_points,                                                   &
                   ent_fields(:,i_temperature),                                &
                   ent_fields(:,i_q_vap),                                      &
                   grid_prev_super(:,i_pressure), ent_rho_dry )
call calc_rho_dry( n_points,                                                   &
                   par_conv_mean_fields(:,i_temperature),                      &
                   par_conv_mean_fields(:,i_q_vap),                            &
                   grid_prev_super(:,i_pressure), par_rho_dry )

! Scale entrainment rate by ratio of dry-densities, to convert
! from volume fraction to dry-mass fraction, then scale by
! the dry-mass flux to convert to actual entrained dry-mass
do ic = 1, n_points
  ent_mass_d(ic) = ent_mass_d(ic)                                              &
                 * ( ent_rho_dry(ic) / par_rho_dry(ic) )                       &
                 * par_conv_super(ic,i_massflux_d)
end do

! Set fractional reduction of entrainment of environment air by the
! parcel core compared to the parcel mean...
if ( l_core_ent_cmr ) then
  ! Fractional contribution from environment depends on PDF-shape;
  ! a more skewed PDF implies a less dilute core.
  do ic = 1, n_points
    core_ent_ratio(ic) = min( core_ent_fac / core_mean_ratio(ic), one )
  end do
else
  ! Fractional contribution from environment is a constant factor
  do ic = 1, n_points
    core_ent_ratio(ic) = core_ent_fac
  end do
end if


!------------------------------------------------------------------------------
! 2) Apply CFL limit to entrainment for numerical stability
!------------------------------------------------------------------------------

! Compute maximum allowed total entrainment rate
! = max_ent_frac times available mass on layer
do ic = 1, n_points
  max_ent(ic) = max_ent_frac * layer_mass_step(ic) / comorph_timestep
end do

! There maybe multiple types / layers of convection entraining
! from the same model-level, in which case it is their sum
! which must not exceed max_ent.  Therefore, scale the
! max ent available to the current type / layer in proportion
! to its fraction of the total mass-flux
do ic = 1, n_points
  max_ent(ic) = max_ent(ic)                                                    &
    * ( par_conv_super(ic,i_massflux_d) / max( sum_massflux(ic), min_float ) )
end do

! Finally, check that entrained dry-mass doesn't exceed the
! amount of dry-mass available on level k
! (numerical stability constraint)
! Only do this if applying a vertically-local CFL limit;
! if not, we'll rescale the closure afterwards to avoid
! hitting this limit on any level
select case ( i_cfl_local )
case ( i_cfl_local_all )
  ! Apply local CFL limit on all levels
  do ic = 1, n_points
    ent_mass_d(ic) = min( ent_mass_d(ic), max_ent(ic) )
  end do
case ( i_cfl_local_nobl )
  ! Only apply local CFL limit above the BL-top
  do ic = 1, n_points
    if ( .not. l_within_bl(ic) ) then
      ent_mass_d(ic) = min( ent_mass_d(ic), max_ent(ic) )
    end if
  end do
end select


return
end subroutine set_ent

end module set_ent_mod
