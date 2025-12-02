! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module conv_level_step_mod

implicit none

contains

! Subroutine to perform one half-level step of the convective
! parcel ascent / descent, for a single updraft or downdraft type
! (either from a full-level to a half-level, or vice-versa)
subroutine conv_level_step(                                                    &
             n_points, max_points, n_points_res,                               &
             n_points_diag, n_diags_super, n_fields_tot,                       &
             l_down, l_tracer, l_last_level, l_to_full_level,                  &
             l_fallback, l_output_fallback, l_within_bl,                       &
             sum_massflux, index_ij, ij_first, ij_last,                        &
             cmpr, k, draft_string,                                            &
             grid_prev_super, grid_next_super,                                 &
             env_prev_super, env_next_super,                                   &
             env_k_fields, layer_mass_step, frac_level_step,                   &
             par_conv_super, par_conv_mean_fields, par_conv_core_fields,       &
             res_source_super, res_source_fields, convcloud_super, fields_2d,  &
             plume_model_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, zero, one, min_float,            &
                                 l_par_core, l_cv_rain, l_cv_snow, l_cv_graup, &
                                 l_precip_to_cloud, l_inter_impl_det,          &
                                 i_convcloud, i_convcloud_none,                &
                                 comorph_timestep, alpha_detrain,              &
                                 core_ent_fac, l_core_ent_cmr,                 &
                                 l_calc_cape, l_calc_mfw_cape,                 &
                                 min_cmr, max_cmr, name_length
use cmpr_type_mod, only: cmpr_type
use grid_type_mod, only: n_grid, i_height, i_pressure
use env_half_mod, only: n_env_half, i_virt_temp
use res_source_type_mod, only: n_res
use cloudfracs_type_mod, only: n_convcloud
use fields_2d_mod, only: n_fields_2d
use parcel_type_mod, only: n_par, i_radius, i_massflux_d, i_edge_virt_temp
use buoyancy_mod, only: i_mean_buoy, i_core_buoy, i_prev
use fields_type_mod, only: i_wind_u, i_temperature,                            &
                           i_q_vap, i_qc_first, i_qc_last,                     &
                           i_q_cl, i_q_rain, i_q_cf, i_q_snow, i_q_graup,      &
                           n_fields, fields_k_conserved_vars
use fields_diags_type_mod, only: fields_diags_copy
use plume_model_diags_type_mod, only: plume_model_diags_type

use calc_virt_temp_mod, only: calc_virt_temp
use dry_adiabat_mod, only: dry_adiabat
use set_ent_mod, only: set_ent
use set_det_mod, only: set_det
use entrain_fields_mod, only: entrain_fields
use entdet_res_source_mod, only: entdet_res_source
use parcel_dyn_mod, only: parcel_dyn, i_call_mean, i_call_core, i_call_det
use set_edge_virt_temp_mod, only: set_edge_virt_temp
use update_par_radius_mod, only: update_par_radius
use set_diag_conv_cloud_mod, only: set_diag_conv_cloud
use calc_cape_mod, only: calc_cape

implicit none

! Number of points in the compression list
integer, intent(in) :: n_points

! Number of points in the compressed super-arrays
! (to avoid repeatedly deallocating and reallocating these,
!  they are dimensioned with the biggest size they will need,
!  which will often be bigger than the number of points here)
integer, intent(in) :: max_points

! Dimensions of the resolved-scale source-term super-arrays
integer, intent(in) :: n_points_res

! Dimensions of the diagnostics super-array
integer, intent(in) :: n_points_diag
integer, intent(in) :: n_diags_super

! Total number of primary fields (including tracers if applicable)
integer, intent(in) :: n_fields_tot

! Flag for downdraft versus updraft
logical, intent(in) :: l_down

! Flag for whether to include passive tracers in the calculations
logical, intent(in) :: l_tracer

! Flag for whether the current model-level is the last allowed
! level (top model-level for updrafts, bottom model-level for
! downdrafts); all mass is forced to detrain if this is true
logical, intent(in) :: l_last_level

! Flag for half-level ascent from half-level to full-level
! (as opposed to full-level to half-level)
logical, intent(in) :: l_to_full_level

! Flag for whether this is the call for the fall-backs;
! Maybe used to specify a different detrainment calculation for
! fall-backs vs primary updrafts / downdrafts
logical, intent(in) :: l_fallback

! Flag for whether to ouput fall-back mass sources based
! on detrained mass which is sufficiently negatively buoyant
! to subsequently fall down below its detrainment level.
logical, intent(in) :: l_output_fallback

! Flag for whether each point is within the boundary-layer
logical, intent(in) :: l_within_bl(n_points)

! Sum of parcel mass-fluxes over types / layers
! at start of this half-level step
real(kind=real_cvprec), intent(in) :: sum_massflux(n_points)

! Indices used to reference a collapsed horizontal coordinate from the parcel
integer, intent(in) :: index_ij(n_points)

! Collapsed horizontal indices of the first and last point available
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Stuff used to make error messages more informative:
! Structure storing i,j indices of compression list points
type(cmpr_type), intent(in) :: cmpr
! Current model-level index
integer, intent(in) :: k
! String identifying what sort of draft (updraft, downdraft, etc)
character(len=name_length), intent(in) :: draft_string

! Super-arrays containing model grid fields (height and pressure)
! at start and end of the current half-level step
real(kind=real_cvprec), target, intent(in) :: grid_prev_super                  &
                                              ( max_points, n_grid )
real(kind=real_cvprec), target, intent(in) :: grid_next_super                  &
                                              ( max_points, n_grid )

! Super-arrays containing any environment fields required
! at start and end of the current half-level step
real(kind=real_cvprec), target, intent(in) :: env_prev_super                   &
                                              ( max_points, n_env_half )
real(kind=real_cvprec), target, intent(in) :: env_next_super                   &
                                              ( max_points, n_env_half )

! Super-array containing the environment primary fields
! at the current thermodynamic level k;
real(kind=real_cvprec), intent(in) :: env_k_fields                             &
                                      ( max_points, n_fields_tot )
! These are the fields used for the properties of entrained air.

! Mass of the current model-level-step / kg m-2
real(kind=real_cvprec), intent(in) :: layer_mass_step( n_points )
! Fraction of full-level-step performed by this call
real(kind=real_cvprec), intent(in) :: frac_level_step( n_points )

! Parcel properties...
! Mass-flux and parcel radius, in-parcel w-weighted means of all
! primary fields, and in-parcel core values of all primary fields.
! in at start of this half-level-step, out at end of it
real(kind=real_cvprec), intent(in out) :: par_conv_super                       &
                                          ( max_points, n_par )
real(kind=real_cvprec), intent(in out) :: par_conv_mean_fields                 &
                                          ( max_points, n_fields_tot )
real(kind=real_cvprec), intent(in out) :: par_conv_core_fields                 &
                                          ( max_points, n_fields_tot )
! NOTE: on input and output the primary-field arrays are not
! in conserved variable form
! (e.g. temperature really does store actual temperature in K).


! Resolved-scale source terms at the current thermodynamic level k...

! Source and sink of dry-mass
real(kind=real_cvprec), intent(in out) :: res_source_super                     &
                                          ( n_points_res, n_res )
! Source terms for primary fields
real(kind=real_cvprec), intent(in out) :: res_source_fields                    &
                                          ( n_points_res, n_fields_tot )
! Convective cloud fields
real(kind=real_cvprec), intent(in out) :: convcloud_super                      &
                                          ( n_points_res, n_convcloud )

! Super-array for 2D fields that might be used elsewhere in comorph
real(kind=real_cvprec), intent(in out) :: fields_2d                            &
                                          ( ij_first:ij_last, n_fields_2d )

! Structure storing flags and super-array addresses for
! various diagnostics
type(plume_model_diags_type), intent(in) :: plume_model_diags
! Super-array to contain output diagnostics
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                          ( n_points_diag, n_diags_super )


! Saved parcel properties at the start of the current level-step:
real(kind=real_cvprec) :: par_prev_super                                       &
                          ( n_points, n_par )
real(kind=real_cvprec) :: par_prev_mean_fields                                 &
                          ( n_points, i_temperature:n_fields )
real(kind=real_cvprec) :: par_prev_core_fields                                 &
                          ( n_points, i_temperature:n_fields )

! Parcel mean and core virtual temperatures
real(kind=real_cvprec) :: par_mean_virt_temp(n_points)
real(kind=real_cvprec) :: par_core_virt_temp(n_points)
! Environment virtual temperature at end of half-level-step, after subsidence
real(kind=real_cvprec) :: env_next_virt_temp(n_points)

! Ratio of parcel core buoyancy over parcel mean buoyancy;
! used in the assumed-PDF in the detrainment calculation.
real(kind=real_cvprec) :: core_mean_ratio(n_points)

! Exner ratio to use when adjusting entrained air temperature
! due to pressure change
real(kind=real_cvprec) :: exner_ratio(n_points)

! Estimated change in environment virtual temperature due to subsidence
real(kind=real_cvprec) :: delta_tv_k(n_points)

! Super-arrays storing mean primary field properties of
! entrained and detrained air
real(kind=real_cvprec) :: ent_fields                                           &
                          ( n_points, n_fields_tot )
real(kind=real_cvprec) :: det_fields                                           &
                          ( n_points, n_fields_tot )

! Entrained and detrained dry-mass
real(kind=real_cvprec) :: ent_mass_d(n_points)
real(kind=real_cvprec) :: det_mass_d(n_points)

! Weight used to calculate properties of air entrained into the core
real(kind=real_cvprec) :: weight(n_points)

! Super-array storing the buoyancies of the parcel mean and core
! at various sub-levels within the current level-step...

! Number of vars in the super-array:
! 1) height
! 2) mean_buoy
! 3) core_buoy
integer, parameter :: n_buoy_vars = 3
! Max number of different heights at which to store these:
! 1) Start of current half-level step
! 2) End of current half-level step
! 3) Height where parcel mean properties first hit saturation
! 4) Height where parcel core first hit saturation
! (note they are stored in height order, not the order above).
integer, parameter :: max_buoy_heights = 4
! Declare the super-array:
real(kind=real_cvprec) :: buoyancy_super                                       &
                     ( n_points, n_buoy_vars, max_buoy_heights )

! Address of end of half-level step in the buoyancy super-array
! (prev is always in address 1)
integer :: i_next(n_points)

! Addresses of saturation heights in the buoyancy super-array
! (0 if didn't hit saturation this step)
integer :: i_mean_sat(n_points)
integer :: i_core_sat(n_points)

! Flag input to conserved variable conversion routine
logical :: l_reverse
! Flag input to entdet_res_source to indicate entrainment vs detrainment
logical :: l_ent
! Flag input to fields_diags_copy indicating fields are in conserved form
logical :: l_conserved_form

! Flag passed into parcel_dyn indicating calls for core vs mean
integer :: i_call

! Pointers to grid and env fields at level k (either prev or next)
real(kind=real_cvprec), pointer :: pressure_k(:)
real(kind=real_cvprec), pointer :: virt_temp_k(:)

! Temporary stores for virtual temperature excesses
real(kind=real_cvprec) :: tmp1, tmp2

! Loop counters
integer :: ic, i_field, i_diag, i_lev


!----------------------------------------------------------------
! 1) Initialisations
!----------------------------------------------------------------

! Save copies of the parcel properties at start of the half-level-step
do i_field = 1, n_par
  do ic = 1, n_points
    par_prev_super(ic,i_field) = par_conv_super(ic,i_field)
  end do
end do
do i_field = i_temperature, n_fields
  do ic = 1, n_points
    par_prev_mean_fields(ic,i_field) = par_conv_mean_fields(ic,i_field)
  end do
end do
if ( l_par_core ) then
  do i_field = i_temperature, n_fields
    do ic = 1, n_points
      par_prev_core_fields(ic,i_field) = par_conv_core_fields(ic,i_field)
    end do
  end do
end if

! Initialise diagnostics to zero
if ( plume_model_diags % n_diags > 0 ) then
  do i_field = 1, plume_model_diags % n_diags_super
    do ic = 1, n_points
      diags_super(ic,i_field) = zero
    end do
  end do
end if

! Initialise heights for start and end buoyancies (prev, next),
! so far assuming the parcel won't hit saturation
do ic = 1, n_points
  i_next(ic) = i_prev + 1
  i_mean_sat(ic) = 0
  i_core_sat(ic) = 0
  buoyancy_super(ic,i_height,i_prev)     = grid_prev_super(ic,i_height)
  buoyancy_super(ic,i_height,i_next(ic)) = grid_next_super(ic,i_height)
end do

! Calculate parcel virtual temperature at start of half-level-step
call calc_virt_temp( n_points, n_points,                                       &
                     par_prev_mean_fields(:,i_temperature),                    &
                     par_prev_mean_fields(:,i_q_vap),                          &
                     par_prev_mean_fields(:,i_qc_first:i_qc_last),             &
                     par_mean_virt_temp )
! Store virtual temperature excess at previous level
do ic = 1, n_points
  buoyancy_super(ic,i_mean_buoy,i_prev) = par_mean_virt_temp(ic)               &
                                        - env_prev_super(ic,i_virt_temp)
end do

if ( l_par_core ) then

  ! Calculate parcel core virtual temperature if needed
  call calc_virt_temp( n_points, n_points,                                     &
                       par_prev_core_fields(:,i_temperature),                  &
                       par_prev_core_fields(:,i_q_vap),                        &
                       par_prev_core_fields(:,i_qc_first:i_qc_last),           &
                       par_core_virt_temp )
  do ic = 1, n_points
    buoyancy_super(ic,i_core_buoy,i_prev) = par_core_virt_temp(ic)             &
                                          - env_prev_super(ic,i_virt_temp)
  end do

  ! Calculate ratio of core buoyancy over mean buoyancy,
  ! relative to the previous existing parcel edge Tv
  do ic = 1, n_points
    ! Initialise to max value
    core_mean_ratio(ic) = max_cmr
  end do
  ! Safety checks depend on whether updraft or downdraft
  if ( l_down ) then
    ! Downdraft: require core and mean to be negatively buoyant
    do ic = 1, n_points
      tmp1 = par_mean_virt_temp(ic) - par_prev_super(ic,i_edge_virt_temp)
      tmp2 = par_core_virt_temp(ic) - par_prev_super(ic,i_edge_virt_temp)
      if ( tmp1 < zero .and. tmp2 < zero )  core_mean_ratio(ic) = tmp2 / tmp1
    end do
  else
    ! Updraft: require core and mean to be positively buoyant
    do ic = 1, n_points
      tmp1 = par_mean_virt_temp(ic) - par_prev_super(ic,i_edge_virt_temp)
      tmp2 = par_core_virt_temp(ic) - par_prev_super(ic,i_edge_virt_temp)
      if ( tmp1 > zero .and. tmp2 > zero )  core_mean_ratio(ic) = tmp2 / tmp1
    end do
  end if
  ! Don't let this go outside specified bounds,
  ! for safety in the detrainment calculations.
  do ic = 1, n_points
    core_mean_ratio(ic) = max( min( core_mean_ratio(ic), max_cmr ), min_cmr )
  end do

end if  ! ( l_par_core )


!----------------------------------------------------------------
! 2) Perform entrainment of environment air from level k
!    into the parcel
!----------------------------------------------------------------

! Set entrained air properties the same as the mean environment at level k
do i_field = 1, n_fields_tot
  do ic = 1, n_points
    ent_fields(ic,i_field) = env_k_fields(ic,i_field)
  end do
end do
! Could do something clever here to make convection entrain
! pre-moistened air which is moister than the env mean fields,
! but not implemented yet.


if ( l_to_full_level ) then
  ! If this is the first of the two half-level-steps
  ! (from previous half-level to level k),
  ! then the environment fields to entrain are at level k, but the
  ! parcel is at a different pressure, at the previous half-level.
  ! Adjust the environment temperature to what it would be at the
  ! start of the level-step, so that we entrain it into the parcel
  ! consistently...
  call dry_adiabat( n_points, n_points,                                        &
                    grid_next_super(:,i_pressure),                             &
                    grid_prev_super(:,i_pressure),                             &
                    ent_fields(:,i_q_vap),                                     &
                    ent_fields(:,i_qc_first:i_qc_last),                        &
                    exner_ratio )
  do ic = 1, n_points
    ent_fields(ic,i_temperature) = ent_fields(ic,i_temperature)                &
                                 * exner_ratio(ic)
  end do
end if

! Set amount of entrained dry-mass over the current half-level-step
call set_ent( n_points, n_fields_tot,                                          &
              max_points, n_points,                                            &
              l_down, l_fallback,                                              &
              par_conv_mean_fields, ent_fields,                                &
              grid_prev_super(:,i_pressure),                                   &
              grid_prev_super(:,i_height), grid_next_super(:,i_height),        &
              par_conv_super(:,i_radius), par_conv_super(:,i_massflux_d),      &
              layer_mass_step, l_within_bl, sum_massflux,                      &
              ent_mass_d )

! Add the entrained mass onto the mass-flux
do ic = 1, n_points
  par_conv_super(ic,i_massflux_d) = par_conv_super(ic,i_massflux_d)            &
                                  + ent_mass_d(ic)
end do

! Convert the entrained fields and the existing parcel
! properties to conserved variable form
l_reverse = .false.
call fields_k_conserved_vars( n_points, n_points,                              &
                              n_fields_tot, l_reverse,                         &
                              ent_fields )
call fields_k_conserved_vars( n_points, max_points,                            &
                              n_fields_tot, l_reverse,                         &
                              par_conv_mean_fields )

! If using parcel core properties, perform entrainment of
! parcel mean properties into the parcel core
if ( l_par_core ) then

  ! Convert core properties to conserved variable form
  l_reverse = .false.
  call fields_k_conserved_vars( n_points, max_points,                          &
                                n_fields_tot, l_reverse,                       &
                                par_conv_core_fields )

  ! The core uses the same fractional entrainment rate as the
  ! parcel mean, but entrains a contribution from the parcel
  ! mean fields, rather than pure environment air
  ! (which makes it far less dilute than the parcel mean).
  ! NOTE: the entrainment of parcel mean air into the core
  ! is calculated BEFORE the parcel mean is modified by
  ! entrainment.  This is important for avoiding vertical
  ! resolution sensitivity; the parcel mean fields after
  ! entrainment but before detrainment are not in balance
  ! (since those 2 processes often have large compensating
  !  effects on parcel buoyancy); the parcel mean fields after
  ! entrainment may have a temporary buoyancy reduction which
  ! scales with the model-level thickness.

  ! Calculate properties of air entrained into the core
  ! (store in detrained fields array, just to save memory).
  if ( l_core_ent_cmr ) then
    ! Fractional contribution from environment depends on PDF-shape;
    ! a more skewed PDF implies a less dilute core.
    do ic = 1, n_points
      weight(ic) = min( core_ent_fac / core_mean_ratio(ic), one )
    end do
  else
    ! Fractional contribution from environment is a constant factor
    do ic = 1, n_points
      weight(ic) = core_ent_fac
    end do
  end if
  do i_field = 1, n_fields_tot
    do ic = 1, n_points
      det_fields(ic,i_field) = weight(ic)  * ent_fields(ic,i_field)            &
                        + (one-weight(ic)) * par_conv_mean_fields(ic,i_field)
    end do
  end do

  ! Calculate effect of entrainment on parcel core properties
  call entrain_fields( n_points, max_points,                                   &
                       n_points, n_fields_tot,                                 &
                       ent_mass_d, par_conv_super(:,i_massflux_d),             &
                       det_fields, par_conv_core_fields )

  ! Convert the core fields back from conserved variable form to normal form
  l_reverse = .true.
  call fields_k_conserved_vars( n_points, max_points,                          &
                                n_fields_tot, l_reverse,                       &
                                par_conv_core_fields )

end if

! Calculate effect of entrainment of environment air on the
! parcel mean properties.
call entrain_fields( n_points, max_points,                                     &
                     n_points, n_fields_tot,                                   &
                     ent_mass_d, par_conv_super(:,i_massflux_d),               &
                     ent_fields, par_conv_mean_fields )

! Convert the mean fields back from conserved variable form to normal form
l_reverse = .true.
call fields_k_conserved_vars( n_points, max_points,                            &
                              n_fields_tot, l_reverse,                         &
                              par_conv_mean_fields )

if ( l_to_full_level ) then
  ! Set entrained temperature back to level k pressure
  do ic = 1, n_points
    ent_fields(ic,i_temperature) = ent_fields(ic,i_temperature)                &
                                 / exner_ratio(ic)
  end do
end if

! Add contribution from entrainment to the resolved-scale source terms
l_ent = .true.
call entdet_res_source( n_points, n_points,                                    &
                        n_points_res, n_fields_tot, l_ent,                     &
                        ent_mass_d, ent_fields,                                &
                        res_source_super, res_source_fields )

! Save diagnostics of entrained mass and air properties
! (in conserved variable form ready for finding mean over types)
if ( plume_model_diags % ent_mass_d % flag ) then
  i_diag = plume_model_diags % ent_mass_d % i_super
  do ic = 1, n_points
    diags_super(ic,i_diag) = ent_mass_d(ic)
  end do
end if
if ( plume_model_diags % ent_fields % n_diags > 0 ) then
  ! Need pressure and Tv at level k; for 2nd half-level-step (from k to the
  ! next model-level interface), "prev" fields are at k.  Otherwise "next" are.
  if ( l_to_full_level ) then
    pressure_k => grid_next_super(:,i_pressure)
    virt_temp_k => env_next_super(:,i_virt_temp)
  else
    pressure_k => grid_prev_super(:,i_pressure)
    virt_temp_k => env_prev_super(:,i_virt_temp)
  end if
  l_conserved_form = .true.
  call fields_diags_copy( n_points, n_points, n_points_diag,                   &
                          n_diags_super, n_fields_tot, l_conserved_form,       &
                          plume_model_diags % ent_fields, ent_fields,          &
                          virt_temp_k, pressure_k, diags_super )
end if


!----------------------------------------------------------------
! 3) Update parcel with moist thermodynamic and dynamical
!    processes occuring over this half-level-step
!----------------------------------------------------------------

if ( l_par_core ) then
  ! If using the parcel core, update the core properties first...

  ! Call parcel dynamics routine; does dry-lifting up to
  ! next level, moist processes, and momentum equations
  i_call = i_call_core
  ! Note: setting i_call_core disables calculation of
  ! resolved-scale sources terms and diagnostics in parcel_dyn.

  call parcel_dyn( n_points, n_points, max_points,                             &
                   max_points, n_points_res, n_points,                         &
                   1, 1, n_fields_tot,                                         &
                   max_buoy_heights, n_buoy_vars, l_down, l_tracer,            &
                   cmpr, k, draft_string, i_call,                              &
                   grid_prev_super, grid_next_super,                           &
                   l_within_bl, layer_mass_step, sum_massflux,                 &
                   par_conv_super(:,i_massflux_d), par_conv_super(:,i_radius), &
                   env_k_fields(:,i_wind_u:i_temperature),                     &
                   env_prev_super(:,i_virt_temp),env_next_super(:,i_virt_temp),&
                   par_prev_core_fields(:,i_temperature:i_qc_last),            &
                   par_conv_core_fields,                                       &
                   buoyancy_super, i_next, i_core_sat )

end if  ! ( l_par_core )

! Call parcel dynamics routine; does dry-lifting up to
! next level, moist processes, and momentum equations
i_call = i_call_mean
! Note: setting i_call_mean enables calculation of both
! resolved-scale sources terms and diagnostics in parcel_dyn.
call parcel_dyn( n_points, n_points, max_points,                               &
                 max_points, n_points_res, n_points,                           &
                 n_points_diag, n_diags_super, n_fields_tot,                   &
                 max_buoy_heights, n_buoy_vars, l_down, l_tracer,              &
                 cmpr, k, draft_string, i_call,                                &
                 grid_prev_super, grid_next_super,                             &
                 l_within_bl, layer_mass_step, sum_massflux,                   &
                 par_conv_super(:,i_massflux_d), par_conv_super(:,i_radius),   &
                 env_k_fields(:,i_wind_u:i_temperature),                       &
                 env_prev_super(:,i_virt_temp), env_next_super(:,i_virt_temp), &
                 par_prev_mean_fields(:,i_temperature:i_qc_last),              &
                 par_conv_mean_fields,                                         &
                 buoyancy_super, i_next, i_mean_sat,                           &
                 res_source_fields = res_source_fields,                        &
                 plume_model_diags = plume_model_diags,                        &
                 diags_super       = diags_super,                              &
                 i_core_sat        = i_core_sat,                               &
                 core_next_q_cl    = par_conv_core_fields(:,i_q_cl),           &
                 core_mean_ratio   = core_mean_ratio )


!---------------------------------------------------------------
! 4) Calculate detrainment, based on an assumed-PDF of buoyancy
!---------------------------------------------------------------

! Convert fields to conserved variable form for detrainment calculations
l_reverse = .false.
call fields_k_conserved_vars( n_points, max_points,                            &
                              n_fields_tot, l_reverse,                         &
                              par_conv_mean_fields )
if ( l_par_core ) then
  call fields_k_conserved_vars( n_points, max_points,                          &
                                n_fields_tot, l_reverse,                       &
                                par_conv_core_fields )
end if

! Estimate change in environment virtual temperature
! due to compensating subsidence.

! First, need exner ratio for subsidence of the environment air.
! Technically, we should take the env q_vap and qc at the end of
! the level-step here, but these are not available, and are
! unlikely to differ much from the values at k, so using those:
call dry_adiabat( n_points, max_points,                                        &
                  grid_next_super(:,i_pressure),                               &
                  grid_prev_super(:,i_pressure),                               &
                  env_k_fields(:,i_q_vap),                                     &
                  env_k_fields(:,i_qc_first:i_qc_last),                        &
                  exner_ratio )

! Now calculate delta_tv_k
! This is the subsided value of Tv from the end of the level-step
! minus the value at the start,
! * dt * massflux / layer-mass
! We also scale by the implicitness weight alpha_detrain.
do ic = 1, n_points
  delta_tv_k(ic) = ( env_next_super(ic,i_virt_temp) * exner_ratio(ic)          &
                   - env_prev_super(ic,i_virt_temp) )                          &
                 * comorph_timestep * alpha_detrain                            &
                 / layer_mass_step(ic)
end do

! For safety, reset delta_tv_k to zero in statically-unstable
! layers so we just use explicit value of env Tv for detrainment.
! Expected sign depends on whether this is updraft or downdraft
if ( l_down ) then
  do ic = 1, n_points
    delta_tv_k(ic) = min( delta_tv_k(ic), zero )
  end do
else
  do ic = 1, n_points
    delta_tv_k(ic) = max( delta_tv_k(ic), zero )
  end do
end if

if ( l_inter_impl_det ) then
  ! If using interdependent implicit detrainment...

  ! Scale by the ratio of the total mass-flux over the mass-flux
  ! from the current convection type / layer, to approximately
  ! account for the subsidence from the other types / layers.
  do ic = 1, n_points
    delta_tv_k(ic) = delta_tv_k(ic)                                            &
                   * sum_massflux(ic)                                          &
                   * ( par_conv_super(ic,i_massflux_d)                         &
                / max( par_prev_super(ic,i_massflux_d), min_float ) )
  end do
  ! NOTE: what we ought to be able to do here is have
  !  delta_tv_k(ic) = delta_tv_k(ic) * par_conv_massflux
  ! outside the if block, and then additionally scale by
  !  ( sum_massflux(ic) / par_prev_massflux )
  ! if l_inter_impl_det is true.
  ! However, the above rearranged order of calculation for l_inter_impl_det
  ! turned out to be needed to avoid a rare floating-point overflow error
  ! (which happened when both par_prev_massflux and par_conv_massflux
  !  were extremely small, but sum_massflux was very large).
  ! Taking the ratio of the 2 very small terms inside the brackets first
  ! avoids the problem.

else  ! ( l_inter_impl_det )
  ! For independent implicit detrainment, only scale by
  ! the mass-flux from the current type / layer
  do ic = 1, n_points
    delta_tv_k(ic) = delta_tv_k(ic)                                            &
                   * par_conv_super(ic,i_massflux_d)
  end do

end if  ! ( l_inter_impl_det )

! The estimated environment Tv after compensating subsidence
! can now be calculated as
!   virt_temp_impl = env_k_virt_temp + delta_tv_k * frac
!
! where frac is the non-detrained fraction of the mass-flux.
!
! Note that this estimate of the final Tv at level k does not
! account for the increments from:
! - entrainment and detrainment.
! - other convective drafts.
! - boundary-layer mixing...
! It is therefore not expected to yield an accurate implicit
! solution.  It is only designed to be a relatively simple
! modification of the explicit detrainment calculation, to
! introduce just enough implicitness to smooth out unphysical
! noisy / intermittent behaviour from one timestep to the next.


! Call routine to set detrained mass and detrained air properties,
! and update the mass-flux and parcel mean properties due to detrainment
call set_det( n_points, n_fields_tot, max_points, max_points,                  &
              max_buoy_heights, n_buoy_vars,                                   &
              l_down, l_last_level, l_to_full_level,                           &
              cmpr, k, draft_string,                                           &
              buoyancy_super, i_next,                                          &
              delta_tv_k, core_mean_ratio,                                     &
              par_conv_mean_fields, par_conv_core_fields,                      &
              par_prev_super(:,i_massflux_d),                                  &
              par_conv_super(:,i_massflux_d),                                  &
              det_mass_d, det_fields )

! Store the final-guess total environment virtual temperature
! increment in delta_tv_k (just needs to be scaled down by the
! final non-detrained fraction)
do ic = 1, n_points
  delta_tv_k(ic) = delta_tv_k(ic)                                              &
                 * par_conv_super(ic,i_massflux_d)                             &
            / max( par_conv_super(ic,i_massflux_d)+det_mass_d(ic), min_float )
end do

! Adjust buoyancies consistent with final guess
! end-of-timestep environment virtual temperature
do i_lev = 1, minval(i_next)
  ! Buoyancy sub-levels used at all points
  do i_field = i_mean_buoy, n_buoy_vars
    do ic = 1, n_points
      buoyancy_super(ic,i_field,i_lev) = buoyancy_super(ic,i_field,i_lev)      &
                                       - delta_tv_k(ic)
    end do
  end do
end do
do i_lev = minval(i_next)+1, maxval(i_next)
  ! Buoyancy sub-levels used at only some points (requires if test
  ! to avoid referencing uninitialised values in buoyancy_super)
  do i_field = i_mean_buoy, n_buoy_vars
    do ic = 1, n_points
      if ( i_lev <= i_next(ic) ) then
        buoyancy_super(ic,i_field,i_lev) = buoyancy_super(ic,i_field,i_lev)    &
                                         - delta_tv_k(ic)
      end if
    end do
  end do
end do

! Set estimated next environment virtual temperature after subsidence
do ic = 1, n_points
  env_next_virt_temp(ic) = env_next_super(ic,i_virt_temp) + delta_tv_k(ic)
end do

! Convert fields back from conserved variable form to normal form
l_reverse = .true.
call fields_k_conserved_vars( n_points, max_points,                            &
                              n_fields_tot, l_reverse,                         &
                              par_conv_mean_fields )
if ( l_par_core ) then
  call fields_k_conserved_vars( n_points, max_points,                          &
                                n_fields_tot, l_reverse,                       &
                                par_conv_core_fields )
end if

if ( .not. l_to_full_level ) then
  ! If this is the second of the two half-level-steps
  ! (from level k to the next half-level),
  ! the detrained air properties have now been calculated at the next
  ! half-level, but need to be detrained to the environment at level k.
  ! So we need to adjust the detrained parcel's thermodynamic fields
  ! from the end of the level-step back to level k...

  ! Convert detrained fields back from conserved variable form to normal form
  l_reverse = .true.
  call fields_k_conserved_vars( n_points, n_points,                            &
                                n_fields_tot, l_reverse,                       &
                                det_fields )

  ! Copy a few of the existing detrained fields for use as the
  ! "prev" fields passed into parcel_dyn
  ! (using the ent_fields array as work-space)
  do i_field = i_temperature, i_qc_last
    do ic = 1, n_points
      ent_fields(ic,i_field) = det_fields(ic,i_field)
    end do
  end do

  ! Call parcel dynamics routine to adjust parcel back to level k.
  i_call = i_call_det
  ! Note: setting i_call_det suppresses the updating of the buoyancy
  ! super-array inside parcel_dyn
  call parcel_dyn( n_points, n_points, n_points,                               &
                   max_points, n_points_res, n_points,                         &
                   1, 1, n_fields_tot,                                         &
                   max_buoy_heights, n_buoy_vars, l_down, l_tracer,            &
                   cmpr, k, draft_string, i_call,                              &
                   grid_next_super, grid_prev_super,                           &
                   l_within_bl, layer_mass_step, sum_massflux,                 &
                   det_mass_d, par_conv_super(:,i_radius),                     &
                   env_k_fields(:,i_wind_u:i_temperature),                     &
                   env_next_super(:,i_virt_temp),env_prev_super(:,i_virt_temp),&
                   ent_fields(:,i_temperature:i_qc_last),                      &
                   det_fields,                                                 &
                   buoyancy_super, i_next, i_mean_sat,                         &
                   res_source_fields = res_source_fields )

  ! Convert detrained fields back into conserved variable form
  l_reverse = .false.
  call fields_k_conserved_vars( n_points, n_points,                            &
                                n_fields_tot, l_reverse,                       &
                                det_fields )

end if  ! ( .not. l_to_full_level )

if ( l_precip_to_cloud ) then
  ! Under this switch, convert any detrained precip into
  ! cloud.  The large-scale microphysics will then
  ! re-decide the partitioning of condensed water into precip
  if ( l_cv_rain ) then
    do ic = 1, n_points
      det_fields(ic,i_q_cl) = det_fields(ic,i_q_cl) + det_fields(ic,i_q_rain)
      det_fields(ic,i_q_rain) = zero
    end do
  end if
  if ( l_cv_snow ) then
    do ic = 1, n_points
      det_fields(ic,i_q_cf) = det_fields(ic,i_q_cf) + det_fields(ic,i_q_snow)
      det_fields(ic,i_q_snow) = zero
    end do
  end if
  if ( l_cv_graup ) then
    do ic = 1, n_points
      det_fields(ic,i_q_cf) = det_fields(ic,i_q_cf) + det_fields(ic,i_q_graup)
      det_fields(ic,i_q_graup) = zero
    end do
  end if
end if

! Compute resolved-scale source terms due to the detrainment
! of the updated detrained air properties at level k
l_ent = .false.
call entdet_res_source( n_points, n_points,                                    &
                        n_points_res, n_fields_tot, l_ent,                     &
                        det_mass_d, det_fields,                                &
                        res_source_super, res_source_fields )

! Save diagnostics of detrained mass and air properties
! (in conserved variable form ready for finding mean over types)
if ( plume_model_diags % det_mass_d % flag ) then
  i_diag = plume_model_diags % det_mass_d % i_super
  do ic = 1, n_points
    diags_super(ic,i_diag) = det_mass_d(ic)
  end do
end if
if ( plume_model_diags % det_fields % n_diags > 0 ) then
  ! Need pressure and Tv at level k; for 2nd half-level-step (from k to the
  ! next model-level interface), "prev" fields are at k.  Otherwise "next" are.
  if ( l_to_full_level ) then
    pressure_k => grid_next_super(:,i_pressure)
    virt_temp_k => env_next_super(:,i_virt_temp)
  else
    pressure_k => grid_prev_super(:,i_pressure)
    virt_temp_k => env_prev_super(:,i_virt_temp)
  end if
  l_conserved_form = .true.
  call fields_diags_copy( n_points, n_points, n_points_diag,                   &
                          n_diags_super, n_fields_tot, l_conserved_form,       &
                          plume_model_diags % det_fields, det_fields,          &
                          virt_temp_k, pressure_k, diags_super )
end if
! Diagnostic of core_mean_ratio
if ( plume_model_diags % core_mean_ratio % flag ) then
  i_diag = plume_model_diags % core_mean_ratio % i_super
  do ic = 1, n_points
    diags_super(ic,i_diag) = core_mean_ratio(ic)
  end do
end if


!---------------------------------------------------------------
! 5) Set other final parcel properties
!---------------------------------------------------------------

! Recalculate the parcel mean buoyancy at the end of the level-step,
! accounting for modification of the parcel mean properties by detrainment.
call calc_virt_temp( n_points, max_points,                                     &
                     par_conv_mean_fields(:,i_temperature),                    &
                     par_conv_mean_fields(:,i_q_vap),                          &
                     par_conv_mean_fields(:,i_qc_first:i_qc_last),             &
                     par_mean_virt_temp )
do ic = 1, n_points
  buoyancy_super(ic,i_mean_buoy,i_next(ic)) = par_mean_virt_temp(ic)           &
                                            - env_next_virt_temp(ic)
end do

! Set new parcel edge virtual temperature at next model-level interface
call set_edge_virt_temp( n_points, n_buoy_vars, max_buoy_heights,              &
                         i_next, buoyancy_super, l_down,                       &
                         env_next_virt_temp, core_mean_ratio,                  &
                         par_conv_super(:,i_edge_virt_temp) )

! Set new parcel radius (accounting for volume change due to
! expansion / contraction and entrainment / detrainment)
call update_par_radius( n_points, det_mass_d,                                  &
  par_prev_mean_fields(:,i_temperature), par_conv_mean_fields(:,i_temperature),&
  par_prev_mean_fields(:,i_q_vap),       par_conv_mean_fields(:,i_q_vap),      &
  grid_prev_super(:,i_pressure),         grid_next_super(:,i_pressure),        &
  par_prev_super(:,i_massflux_d),        par_conv_super(:,i_massflux_d),       &
  par_prev_super(:,i_radius),            par_conv_super(:,i_radius) )

! NOTE: Processes which may modify the parcel radius but are
! not yet represented:
! - Splitting and merging with height
! - Selective detrainment of smaller thermals from the ensemble.

if ( i_convcloud > i_convcloud_none ) then
  ! Set diagnosed convective cloud-fraction
  ! and liquid & ice mixing ratios
  call set_diag_conv_cloud(                                                    &
         n_points, max_points, n_points_res,                                   &
         n_points, max_points, n_fields_tot,                                   &
         max_buoy_heights, n_buoy_vars,                                        &
         grid_prev_super, grid_next_super, frac_level_step,                    &
         env_prev_super(:,i_virt_temp),  env_next_super(:,i_virt_temp),        &
         par_prev_super(:,i_massflux_d), par_conv_super(:,i_massflux_d),       &
         par_prev_super(:,i_radius),     par_conv_super(:,i_radius),           &
         par_prev_mean_fields(:,i_temperature:n_fields), par_conv_mean_fields, &
         buoyancy_super, i_next, i_mean_sat,                                   &
         convcloud_super )
end if

if ( l_calc_cape .or. l_calc_mfw_cape ) then
  ! Calculate contribution to CAPE from current half-level-step
  call calc_cape(                                                              &
         n_points, max_buoy_heights, n_buoy_vars,                              &
         buoyancy_super, i_next, l_within_bl,                                  &
         index_ij, ij_first, ij_last,                                          &
         env_prev_super(:,i_virt_temp), env_next_super(:,i_virt_temp),         &
         par_prev_super(:,i_massflux_d), par_conv_super(:,i_massflux_d),       &
         fields_2d )
end if


return
end subroutine conv_level_step

end module conv_level_step_mod
