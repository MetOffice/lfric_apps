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
subroutine set_ent( n_points, n_fields_tot,                                    &
                    n_points_par_super, n_points_ent_super,                    &
                    l_down, l_fallback,                                        &
                    par_fields, ent_fields,                                    &
                    pressure, prev_height, next_height,                        &
                    par_radius, massflux_d,                                    &
                    layer_mass, l_within_bl, sum_massflux,                     &
                    ent_mass_d )

use comorph_constants_mod, only: real_cvprec, min_float,                       &
                     ent_coef, comorph_timestep, i_cfl_check,                  &
                     i_cfl_check_local, i_cfl_check_hybrid,                    &
                     max_ent_frac_up, max_ent_frac_up_fb,                      &
                     max_ent_frac_dn, max_ent_frac_dn_fb
use fields_type_mod, only: i_temperature, i_q_vap
use calc_rho_dry_mod, only: calc_rho_dry

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of fields (including tracres if applicable)
integer, intent(in) :: n_fields_tot

! Number of points in the compressed parcel and entrained fields super-arrays
! (to avoid repeatedly deallocating and reallocating these,
!  they are dimensioned with the biggest size they will need,
!  which will often be bigger than the number of points here)
integer, intent(in) :: n_points_par_super
integer, intent(in) :: n_points_ent_super

! Flag for downdraft vs updraft
logical, intent(in) :: l_down
! Flag for fall-back flow calculation
logical, intent(in) :: l_fallback

! Super-array containing the parcel mean primary fields
real(kind=real_cvprec), intent(in) :: par_fields                               &
                                      ( n_points_par_super, n_fields_tot )
! Entrained primary field values
real(kind=real_cvprec), intent(in) :: ent_fields                               &
                                      ( n_points_ent_super, n_fields_tot )

! Pressure of entrained and parcel air
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Height of parcel before and after the level-step
real(kind=real_cvprec), intent(in) :: prev_height(n_points)
real(kind=real_cvprec), intent(in) :: next_height(n_points)

! Parcel radius used to set rate of mixing with the environment
real(kind=real_cvprec), intent(in) :: par_radius(n_points)

! Mass-flux present before entrainment from the current level k
real(kind=real_cvprec), intent(in) :: massflux_d(n_points)

! Dry-mass on the current model-level, per m2 at surface.
real(kind=real_cvprec), intent(in) :: layer_mass(n_points)

! Flag for whether each point is within the boundary-layer
logical, intent(in) :: l_within_bl(n_points)

! Sum of previous level-interface mass-fluxes over all
! convection types / layers
real(kind=real_cvprec), intent(in) :: sum_massflux(n_points)

! Rate of entrainment of dry-mass from current level / kg m-2 s-1
real(kind=real_cvprec), intent(out) :: ent_mass_d(n_points)


! Dry-density of the entrained air and the parcel
real(kind=real_cvprec) :: ent_rho_dry(n_points)
real(kind=real_cvprec) :: par_rho_dry(n_points)

! Max fraction of model-level-mass allowed to be entrained
real(kind=real_cvprec) :: max_ent_frac

! Maximum allowed entrainment for numerical stability
real(kind=real_cvprec) :: max_ent(n_points)

! Loop counters
integer :: ic


! Set the fractional volume entrainment rate in m-1
! Currently only including contribution from mixing
do ic = 1, n_points
  ent_mass_d(ic) = ent_coef / par_radius(ic)
end do

! Scale by env-relative distance travelled by parcel to get
! volume fraction entrained over the layer
do ic = 1, n_points
  ent_mass_d(ic) = ent_mass_d(ic)                                              &
                 * abs( next_height(ic) - prev_height(ic) )
  ! Note: for downdrafts, next_height < prev_height, so need
  ! abs call to ensure dz is positive.
end do

! Note: here we could account for horizontal distance travelled
! due to slant-wise ascent as well, but not implemented yet.

! Compute dry-density of the entrained air and the parcel
call calc_rho_dry( n_points,                                                   &
                   ent_fields(:,i_temperature),                                &
                   ent_fields(:,i_q_vap),                                      &
                   pressure, ent_rho_dry )
call calc_rho_dry( n_points,                                                   &
                   par_fields(:,i_temperature),                                &
                   par_fields(:,i_q_vap),                                      &
                   pressure, par_rho_dry )

! Scale entrainment rate by ratio of dry-densities, to convert
! from volume fraction to dry-mass fraction, then scale by
! the dry-mass flux to convert to actual entrained dry-mass
do ic = 1, n_points
  ent_mass_d(ic) = ent_mass_d(ic)                                              &
                 * ( ent_rho_dry(ic) / par_rho_dry(ic) )                       &
                 * massflux_d(ic)
end do


! Calculate max allowed entrainment for numerical stability...

! Select the correct max entrainment fraction for this draft
if ( l_down ) then
  if ( l_fallback ) then
    max_ent_frac = max_ent_frac_dn_fb
  else
    max_ent_frac = max_ent_frac_dn
  end if
else
  if ( l_fallback ) then
    max_ent_frac = max_ent_frac_up_fb
  else
    max_ent_frac = max_ent_frac_up
  end if
end if

! Compute maximum allowed total entrainment rate
do ic = 1, n_points
  max_ent(ic) = max_ent_frac * layer_mass(ic) / comorph_timestep
end do

! There maybe multiple types / layers of convection entraining
! from the same model-level, in which case it is their sum
! which must not exceed max_ent.  Therefore, scale the
! max ent available to the current type / layer in proportion
! to its fraction of the total mass-flux
do ic = 1, n_points
  max_ent(ic) = max_ent(ic)                                                    &
              * ( massflux_d(ic) / max( sum_massflux(ic), min_float ) )
end do

! Finally, check that entrained dry-mass doesn't exceed the
! amount of dry-mass available on level k
! (numerical stability constraint)
! Only do this if applying a vertically-local CFL limit;
! if not, we'll rescale the closure afterwards to avoid
! hitting this limit on any level
select case(i_cfl_check)
case (i_cfl_check_local)
  ! Vertically local CFL limiting; always use local limit:
  do ic = 1, n_points
    ent_mass_d(ic) = min( ent_mass_d(ic), max_ent(ic) )
  end do
case (i_cfl_check_hybrid)
  ! Hybrid CFL limit; only apply local limit when above
  ! the boundary-layer top
  do ic = 1, n_points
    if ( .not. l_within_bl(ic) ) then
      ent_mass_d(ic) = min( ent_mass_d(ic), max_ent(ic) )
    end if
  end do
end select


return
end subroutine set_ent

end module set_ent_mod
