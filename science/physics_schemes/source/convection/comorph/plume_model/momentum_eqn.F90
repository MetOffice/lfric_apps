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

subroutine momentum_eqn( n_points, n_fields_tot, l_res_source, l_down,         &
                         n_points_env, n_points_next,                          &
                         n_points_res, n_points_sublevs,                       &
                         l_within_bl, layer_mass_step, sum_massflux,           &
                         massflux_d, par_radius, q_tot, delta_t, Nsq_dry,      &
                         sublevs, i_next, j_wex,                               &
                         env_k_winds, par_next_winds, par_w_drag,              &
                         res_source_fields )

use comorph_constants_mod, only: real_cvprec, zero, one, two,                  &
                                 three_over_eight, half,                       &
                                 comorph_timestep, drag_coef_par, wavedrag_fac,&
                                 l_homog_conv_bl
use fields_type_mod, only: i_wind_u, i_wind_v, i_wind_w
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev, j_height, j_env_w

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of primary fields in the resolved-scale source term array
! (may or may not include tracers)
integer, intent(in) :: n_fields_tot

! Switch for whether resolved-scale source terms need to be calculated
logical, intent(in) :: l_res_source

! Flag for downdraft vs updraft
logical, intent(in) :: l_down

! Number of points in the parcel and environment super-arrays
! (maybe larger than needed here, to save having to reallocate)
integer, intent(in) :: n_points_env
integer, intent(in) :: n_points_next
integer, intent(in) :: n_points_res
integer, intent(in) :: n_points_sublevs

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

! Environment dry static stability
real(kind=real_cvprec), intent(in) :: Nsq_dry(n_points)

! Super-array storing the parcel buoyancies and other properties
! at up to 4 sub-level heights within the current level-step:
! a) Start of the level-step
! b) End of the level-step
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(in out) :: sublevs                              &
                     ( n_points_sublevs, n_sublev_vars, max_sublevs )

! Address of next model-level interface in sublevs
integer, intent(in) :: i_next(n_points)

! Address of mean or core w-excess in sublevs, as appropriate
integer, intent(in) :: j_wex

! Environment winds at current full-level
real(kind=real_cvprec), intent(in) :: env_k_winds                              &
                                      ( n_points_env, i_wind_u:i_wind_w )

! Parcel winds to be updated
real(kind=real_cvprec), intent(in out) :: par_next_winds                       &
                                          ( n_points_next, i_wind_u:i_wind_w )

! Drag on the parcel vertical velocity / s-1
! Output for use in updating wind_w later, in the detrainment calculation
real(kind=real_cvprec), intent(out) :: par_w_drag(n_points)

! Resolved-scale source terms
real(kind=real_cvprec), optional, intent(in out) :: res_source_fields          &
                                          ( n_points_res, n_fields_tot )

! Store for drag factor common to all wind components
real(kind=real_cvprec) :: drag_fac(n_points)

! Store for reaction force term  1 + M dt / (rho dz)
real(kind=real_cvprec) :: reaction_term(n_points)

! Increment to each parcel wind component
real(kind=real_cvprec) :: dwindp ( n_points, i_wind_u:i_wind_w )

! Parcel wind minus environment wind (vector magnitude)
real(kind=real_cvprec) :: wind_ex(n_points)

! Term in the implicit solution of pressure drag
real(kind=real_cvprec) :: alpha

! Vertical interpolation weight
real(kind=real_cvprec) :: interp

! Loop counters
integer :: ic, i_field, i_lev


! Initialise parcel wind increments to zero
! (no explicit wind increments added yet, but will be added here)
do i_field = i_wind_u, i_wind_w
  do ic = 1, n_points
    dwindp(ic,i_field) = zero
  end do
end do


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


do ic = 1, n_points

  ! Compute magnitude of vector wind difference between parcel and environment
  wind_ex(ic) = sqrt(                                                          &
      ( par_next_winds(ic,i_wind_u) - env_k_winds(ic,i_wind_u) )**2            &
    + ( par_next_winds(ic,i_wind_v) - env_k_winds(ic,i_wind_v) )**2            &
    + ( par_next_winds(ic,i_wind_w) - env_k_winds(ic,i_wind_w) )**2 )

  ! Precalculate and store the term  3/8 coef 1/R
  ! (same for all wind compenents)
  drag_fac(ic) = three_over_eight * drag_coef_par * delta_t(ic)                &
               / par_radius(ic)

  ! Also precalculate the reaction force term 1 + M dt / (rho dz)
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

  do ic = 1, n_points

    ! Calculate alpha term
    alpha = drag_fac(ic) * abs( env_k_winds(ic,i_field)                        &
                              - par_next_winds(ic,i_field) )
    ! Note: we should really use the vector magnitude of the wind difference
    ! here, not just the difference of the current wind component;
    ! this is a bug and will make the CMT sensitive to grid orientation.

    ! Compute increment to parcel wind (accounting for any explicit increments)
    dwindp(ic,i_field) = dwindp(ic,i_field)                                    &
             - alpha * ( par_next_winds(ic,i_field) + dwindp(ic,i_field)       &
                       - env_k_winds(ic,i_field) )                             &
                     / ( one + alpha * reaction_term(ic) )

  end do

end do  ! i_field = i_wind_u, i_wind_v

! Update parcel winds with explicit increments + horizontal drag increments
do i_field = i_wind_u, i_wind_w
  do ic = 1, n_points
    par_next_winds(ic,i_field) = par_next_winds(ic,i_field)                    &
                               + dwindp(ic,i_field)
  end do
end do

! If resolved-scale source terms are needed
if ( l_res_source .and. present(res_source_fields) ) then
  do i_field = i_wind_u, i_wind_w
    do ic = 1, n_points
      ! Add contribution to the resolved momentum source terms due to the
      ! reaction force from the drag on the parcel
      res_source_fields(ic,i_field) = res_source_fields(ic,i_field)            &
                                    - dwindp(ic,i_field) * massflux_d(ic)      &
                                                         * ( one + q_tot(ic) )
      ! Momentum tendency scales with the wet-mass flux, so scaling
      ! dry-mass flux by (1 + qt)
    end do
  end do
end if

! Now for the vertical wind component of drag...
! Note we will use a different numerical method to solve the w-equation
! (vs the u,v equations), because w naturally couples to delta_t
! and will be modified by the change in buoyancy in the implicit
! solution of detrainment.
! Therefore, in this routine we only calculate the drag coefficient;
! the parcel w is updated in the detrainment routine set_det.

! We have:
! Dwp/Dt = b' + 3/8 coef 1/R abs( we - wp )
!               1/2 ( 1 + sqrt( 1 + s N^2 R^2 / (we-wp)^2 ) )
!               ( we - wp )
!
!        = b' + 3/8 coef 1/R 1/2 ( abs(we-wp) + sqrt( (we-wp)^2 + s N^2 R^2 ) )
!               ( we - wp )
!
! This is the same equation as for the horizontal momentum, except that
! we have extra terms:
! b'
! is acceleration due to buoyancy and pressure gradient forces.
! 1/2 ( 1 + sqrt( 1 + s N^2 R^2 / (we-wp)^2 ) )
! is an extra factor in the drag to account for gravity-wave drag.
!
! We abbreviate this to:
! Dwp/Dt = b' - d ( wp - we )
!
! where the drag coefficient / s-1 is
! d = 3/8 coef 1/R 1/2 ( abs(we-wp) + sqrt( (we-wp)^2 + s N^2 R^2 ) )

! Calculate the drag coefficient in s-1
do ic = 1, n_points
  par_w_drag(ic) = ( drag_fac(ic) / delta_t(ic) )                              &
        * half * ( wind_ex(ic)                                                 &
                 + sqrt( wind_ex(ic)**2 + wavedrag_fac * max(Nsq_dry(ic),zero) &
                                                       * par_radius(ic)**2 ) )
end do
! Reset drag to zero if the w-excess has the wrong sign
! (shouldn't be possible e.g. for an updraft to be descending relative to
!  its environment, but can happen if not using the w-excess for detrainment).
! Where the drag is used, it is assumed to always accelerate updrafts
! downwards and downdrafts upwards, which results in a run-away positive
! feedback if w-excess doesn't have the expected sign!
if ( l_down ) then
  do ic = 1, n_points
    if ( par_next_winds(ic,i_wind_w) > env_k_winds(ic,i_wind_w) ) then
      par_w_drag(ic) = zero
    end if
  end do
else
  do ic = 1, n_points
    if ( par_next_winds(ic,i_wind_w) < env_k_winds(ic,i_wind_w) ) then
      par_w_drag(ic) = zero
    end if
  end do
end if

! Initialise sub-level w excesses to values accounting for entrainment
! but not yet accounting for buoyancy and drag
! (the buoyancies need to account for the subsidence increment to Tv
!  and so need to be implicitly solved in the detrainment calculation).
do ic = 1, n_points
  do i_lev = i_prev+1, i_next(ic)
    ! Linearly ramp the entrainment already added onto par_next_winds:
    interp = (sublevs(ic,j_height,i_lev)      - sublevs(ic,j_height,i_prev))   &
           / (sublevs(ic,j_height,i_next(ic)) - sublevs(ic,j_height,i_prev))
    sublevs(ic,j_wex,i_lev)                                                    &
      = (one-interp)*( sublevs(ic,j_env_w,i_prev) + sublevs(ic,j_wex,i_prev) ) &
      +      interp *par_next_winds(ic,i_wind_w)                               &
      - sublevs(ic,j_env_w,i_lev)
  end do
end do


return
end subroutine momentum_eqn

end module momentum_eqn_mod
