
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module momentum_eqn_mod

implicit none

contains

subroutine momentum_eqn( n_points, n_fields_tot, l_res_source,                 &
                         n_points_env, n_points_next, n_points_res,            &
                         l_within_bl, layer_mass_step, sum_massflux,           &
                         massflux_d, par_radius, q_tot, delta_t,               &
                         env_k_winds, par_next_winds, res_source_fields )

use comorph_constants_mod, only: real_cvprec, one, three_over_eight,           &
                                 comorph_timestep, drag_coef_par,              &
                                 l_homog_conv_bl
use fields_type_mod, only: i_wind_u, i_wind_v, i_wind_w

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of primary fields in the resolved-scale source term array
! (may or may not include tracers)
integer, intent(in) :: n_fields_tot

! Switch for whether resolved-scale source terms need to be calculated
logical, intent(in) :: l_res_source

! Number of points in the parcel and environment super-arrays
! (maybe larger than needed here, to save having to reallocate)
integer, intent(in) :: n_points_env
integer, intent(in) :: n_points_next
integer, intent(in) :: n_points_res

! Flag for whether each point is within the homogenised mixed-layer
logical, intent(in) :: l_within_bl(n_points)

! Dry-mass on the current model-level step
real(kind=real_cvprec), intent(in) :: layer_mass_step(n_points)

! Sum of parcel mass-fluxes over types / layers
real(kind=real_cvprec), intent(in) :: sum_massflux(n_points)

! Dry-mass flux of the parcel
real(kind=real_cvprec), intent(in) :: massflux_d(n_points)

! Radius length-scale of the parcel
real(kind=real_cvprec), intent(in) :: par_radius(n_points)

! Total-water mixing ratio of the parcel
real(kind=real_cvprec), intent(in) :: q_tot(n_points)

! Time taken for the parcel to ascend across the current level-step
real(kind=real_cvprec), intent(in) :: delta_t(n_points)

! Environment winds at current full-level
real(kind=real_cvprec), intent(in) :: env_k_winds                              &
                                      ( n_points_env, i_wind_u:i_wind_w )

! Parcel winds to be updated
real(kind=real_cvprec), intent(in out) :: par_next_winds                       &
                                          ( n_points_next, i_wind_u:i_wind_w )

! Resolved-scale source terms
real(kind=real_cvprec), optional, intent(in out) :: res_source_fields          &
                                          ( n_points_res, n_fields_tot )

! Store for drag factor common to all wind components
real(kind=real_cvprec) :: drag_fac(n_points)

! Store for reaction force term  1 + M dt / (rho dz)
real(kind=real_cvprec) :: reaction_term(n_points)

! Increment to current parcel wind component
real(kind=real_cvprec) :: dwindp(n_points)

! Term in the implicit solution of pressure drag
real(kind=real_cvprec) :: alpha

! Loop counters
integer :: ic, i_field

! PRESSURE DRAG ON THE PARCEL
! The parcel u,v are adjusted towards the environment u,v by
! a horizontal drag force.
! Remember that this creates an equal and opposite reaction force
! on the environment u,v, so while the parcel u,v go towards the env,
! the env values will adjust towards the parcel at the same time.
! Therefore it is easy to numericaly overshoot when the mass-flux
! is large, making the parcel and env values "swap places" instead of
! moving closer together, causing a spurious numerical oscillation.
! To avoid this instability, we use an implicit-in-time discretisation,
! so that the drag depends on estimates of the parcel and env u,v
! after the increments have been applied...
!
! We assume a classical quadratic drag law:
!
! dup/dt = 3/8 coef 1/R abs( ue - up ) ( ue - up )
! => dup = alpha ( ue - up )
! where alpha = 3/8 coef 1/R abs( ue - up ) dz/w
!
! Now the reaction force on the env:
! due = -M dt / (rho dz) dup
!
! Combining:
! dup = alpha ( ue_n - M dt / (rho dz) dup - up_n - dup )
! => dup ( 1 + alpha ( 1 + M dt / (rho dz) ) ) = alpha ( ue_n - up_n )
! => dup = alpha ( ue_n - up_n ) / ( 1 + alpha ( 1 + M dt / (rho dz) ) )

! Precalculate and store the term  3/8 coef 1/R dz/w
! (same for all wind compenents)
do ic = 1, n_points
  drag_fac(ic) = three_over_eight * drag_coef_par * delta_t(ic)                &
               / par_radius(ic)
end do

! Also precalculate the reaction force term 1 + M dt / (rho dz)
do ic = 1, n_points
  reaction_term(ic) = one + sum_massflux(ic) * comorph_timestep                &
                          / layer_mass_step(ic)
  ! Note: using sum of mass-fluxes over all layers/types here, to ensure
  ! numerical stability when multiple parcels are interacting with the
  ! environment simultaneously.
end do

if ( l_homog_conv_bl ) then
  ! If we are homogenising the convective source terms within the surface
  ! mixed-layer, the reaction term is wrong below the BL-top,
  ! as the momentum source terms will be reset to a vertically-uniform profile
  ! that conserves momentum, so-as not to double-count the turbulent
  ! wind stresses calculated by the boundary-layer scheme.
  ! In this case, just treat the env u,v profiles explicitly;
  ! reset the reaction force term to 1.
  do ic = 1, n_points
    if ( l_within_bl(ic) ) reaction_term(ic) = one
  end do
end if

do i_field = i_wind_u, i_wind_v
  ! For each horizontal wind component...
  ! (momentum equation for the vertical winds is still under development)

  do ic = 1, n_points

    ! Calculate alpha term
    alpha = drag_fac(ic) * abs( env_k_winds(ic,i_field)                        &
                              - par_next_winds(ic,i_field) )

    ! Compute increment to parcel wind
    dwindp(ic) = alpha * ( env_k_winds(ic,i_field)                             &
                         - par_next_winds(ic,i_field) )                        &
               / ( one + alpha * reaction_term(ic) )

    ! Update parcel wind
    par_next_winds(ic,i_field) = par_next_winds(ic,i_field) + dwindp(ic)

  end do

  ! If resolved-scale source terms are needed
  if ( l_res_source .and. present(res_source_fields) ) then
    do ic = 1, n_points
      ! Add contribution to the resolved momentum source terms due to the
      ! reaction force from the drag on the parcel
      res_source_fields(ic,i_field) = res_source_fields(ic,i_field)            &
                                    - dwindp(ic) * massflux_d(ic)              &
                                                 * ( one + q_tot(ic) )
      ! Momentum tendency scales with the wet-mass flux, so scaling
      ! dry-mass flux by (1 + qt)
    end do
  end if

end do  ! i_field = i_wind_u, i_wind_v

return
end subroutine momentum_eqn

end module momentum_eqn_mod
