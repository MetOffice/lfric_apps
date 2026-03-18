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
             l_within_bl,                                                      &
             sum_massflux, max_ent_frac, index_ij, ij_first, ij_last,          &
             cmpr, k, draft_string,                                            &
             grid_prev_super, grid_next_super,                                 &
             env_prev_super, env_next_super,                                   &
             env_k_fields, layer_mass_step, frac_level_step,                   &
             delta_tv,                                                         &
             par_conv_super, par_conv_mean_fields, par_conv_core_fields,       &
             res_source_super, res_source_fields,                              &
             convcloud_super, fields_2d,                                       &
             plume_model_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, zero, one, min_float,            &
                                 l_par_core, i_impl_det, i_impl_det_indep,     &
                                 i_impl_det_inter,                             &
                                 i_convcloud, i_convcloud_none,                &
                                 i_cf_conv, i_cf_conv_coma,                    &
                                 l_calc_cape, l_calc_mfw_cape,                 &
                                 name_length
use cmpr_type_mod, only: cmpr_type
use grid_type_mod, only: n_grid, i_height, i_pressure
use env_half_mod, only: n_env_half, i_virt_temp
use res_source_type_mod, only: n_res
use cloudfracs_type_mod, only: n_convcloud
use fields_2d_mod, only: n_fields_2d
use parcel_type_mod, only: n_par, i_radius, i_massflux_d
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev,                     &
                       j_mean_buoy, j_core_buoy,                               &
                       j_env_tv, j_delta_tv, j_massflux_d
use fields_type_mod, only: i_wind_w, i_temperature,                            &
                           i_q_vap, i_qc_first, i_qc_last,                     &
                           n_fields, fields_k_conserved_vars
use fields_diags_type_mod, only: fields_diags_copy
use plume_model_diags_type_mod, only: plume_model_diags_type

use calc_virt_temp_mod, only: calc_virt_temp
use dry_adiabat_mod, only: dry_adiabat
use calc_core_mean_ratio_mod, only: calc_core_mean_ratio
use calc_env_nsq_mod, only: calc_env_nsq
use set_ent_mod, only: set_ent
use set_det_mod, only: set_det
use entrain_fields_mod, only: entrain_fields
use entdet_res_source_mod, only: entdet_res_source
use init_sublevs_mod, only: init_sublevs
use parcel_dyn_mod, only: parcel_dyn, i_call_mean, i_call_core, i_call_det
use set_par_cloudfrac_mod, only: set_par_cloudfrac
use calc_rho_dry_mod, only: calc_rho_dry
use update_edge_virt_temp_mod, only: update_edge_virt_temp
use update_par_radius_mod, only: update_par_radius
use set_diag_conv_cloud_a_mod, only: set_diag_conv_cloud_a
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

! Flag for whether each point is within the boundary-layer
logical, intent(in) :: l_within_bl(n_points)

! Sum of parcel mass-fluxes over types/layers at start of this half-level step
real(kind=real_cvprec), intent(in) :: sum_massflux(n_points)

! Maximum allowed entrained fraction of layer mass for current draft
real(kind=real_cvprec), intent(in) :: max_ent_frac

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

! Difference in virtual temperature between level k and next full-level,
! divided by layer-mass (used to pre-estimate subsidence increment to Tv)
real(kind=real_cvprec), intent(in) :: delta_tv( n_points )

! Parcel properties...
! Mass-flux and parcel radius, in-parcel w-weighted means of all
! primary fields, and in-parcel core values of all primary fields.
! IN at start of this half-level-step, OUT at end of it
real(kind=real_cvprec), intent(in out) :: par_conv_super                       &
                                          ( max_points, n_par )
real(kind=real_cvprec), intent(in out) :: par_conv_mean_fields                 &
                                          ( max_points, n_fields_tot )
real(kind=real_cvprec), intent(in out) :: par_conv_core_fields                 &
                                          ( max_points, n_fields_tot )
! NOTE: on input and output the primary-field arrays are NOT
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
                          ( n_points, n_fields )
real(kind=real_cvprec) :: par_prev_core_fields                                 &
                          ( n_points, n_fields )

! Parcel mean and core virtual temperatures
real(kind=real_cvprec) :: par_mean_virt_temp(n_points)
real(kind=real_cvprec) :: par_core_virt_temp(n_points)

! Parcel mean and core vertical velocity drag / s-1
real(kind=real_cvprec) :: par_mean_w_drag(n_points)
real(kind=real_cvprec) :: par_core_w_drag(n_points)

! Ratio of parcel core buoyancy over parcel mean buoyancy;
! used in the assumed-PDF in the detrainment calculation.
real(kind=real_cvprec) :: core_mean_ratio(n_points)

! Environment dry static stability
real(kind=real_cvprec) :: Nsq_dry(n_points)

! Exner ratio to use when adjusting entrained air temperature
! due to pressure change
real(kind=real_cvprec) :: exner_ratio(n_points)

! Super-arrays storing mean primary field properties of
! entrained and detrained air
real(kind=real_cvprec) :: ent_fields                                           &
                          ( n_points, n_fields_tot )
real(kind=real_cvprec) :: det_fields                                           &
                          ( n_points, n_fields_tot )

! Entrained and detrained dry-mass / kg m-2 s-1
real(kind=real_cvprec) :: ent_mass_d(n_points)
real(kind=real_cvprec) :: det_mass_d(n_points)

! Weight used to calculate properties of air entrained into the core
real(kind=real_cvprec) :: core_ent_ratio(n_points)

! Supersaturation the parcel would have with all liquid cloud evaporated
! (values for parcel mean and core, and previous and next model-level)
real(kind=real_cvprec) :: par_prev_mean_ss(n_points)
real(kind=real_cvprec) :: par_next_mean_ss(n_points)
real(kind=real_cvprec) :: par_prev_core_ss(n_points)
real(kind=real_cvprec) :: par_next_core_ss(n_points)

! Virtual temperature the parcel would have with all liquid cloud evaporated
real(kind=real_cvprec) :: par_prev_mean_tvl(n_points)
real(kind=real_cvprec) :: par_next_mean_tvl(n_points)
real(kind=real_cvprec) :: par_prev_core_tvl(n_points)
real(kind=real_cvprec) :: par_next_core_tvl(n_points)

! Total mass-flux at prev used in the implicit detrainment calculation
! (optionally includes contributions from other layers / types)
real(kind=real_cvprec) :: sum_massflux_det(n_points)

! Super-array storing the parcel buoyancies and other properties
! at various sub-levels within the current level-step...
real(kind=real_cvprec) :: sublevs                                              &
                          ( n_points, n_sublev_vars, max_sublevs )

! Address of end of half-level step in the sub-levels super-array
! (prev is always in address 1)
integer :: i_next(n_points)

! Address of saturation height in the sub-levels super-array
! (0 if didn't hit saturation this step)
integer :: i_sat(n_points)
integer :: i_core_sat(n_points)

! Max sub-level height address being used (max value in i_next)
integer :: i_next_max

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
do i_field = 1, n_fields
  do ic = 1, n_points
    par_prev_mean_fields(ic,i_field) = par_conv_mean_fields(ic,i_field)
  end do
end do
if ( l_par_core ) then
  do i_field = 1, n_fields
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

! Need pressure and Tv at level k; for 2nd half-level-step (from k to the
! next model-level interface), "prev" fields are at k.  Otherwise "next" are.
if ( l_to_full_level ) then
  pressure_k => grid_next_super(:,i_pressure)
  virt_temp_k => env_next_super(:,i_virt_temp)
else
  pressure_k => grid_prev_super(:,i_pressure)
  virt_temp_k => env_prev_super(:,i_virt_temp)
end if

! Calculate parcel virtual temperature at start of half-level-step
call calc_virt_temp( n_points, n_points,                                       &
                     par_prev_mean_fields(:,i_temperature),                    &
                     par_prev_mean_fields(:,i_q_vap),                          &
                     par_prev_mean_fields(:,i_qc_first:i_qc_last),             &
                     par_mean_virt_temp )

if ( l_par_core ) then

  ! Calculate parcel core virtual temperature if needed
  call calc_virt_temp( n_points, n_points,                                     &
                       par_prev_core_fields(:,i_temperature),                  &
                       par_prev_core_fields(:,i_q_vap),                        &
                       par_prev_core_fields(:,i_qc_first:i_qc_last),           &
                       par_core_virt_temp )

  ! Calculate ratio of core buoyancy over mean buoyancy,
  ! relative to the previous existing parcel edge Tv
  call calc_core_mean_ratio( n_points, n_points, l_down,                       &
                             par_mean_virt_temp, par_core_virt_temp,           &
                             par_prev_super, core_mean_ratio )

end if  ! ( l_par_core )

! Calculate environment static stability
call calc_env_nsq( n_points, max_points, l_to_full_level,                      &
                   env_prev_super(:,i_virt_temp),env_next_super(:,i_virt_temp),&
                   grid_prev_super, grid_next_super, env_k_fields,             &
                   Nsq_dry )

! Implicit detrainment options...
select case ( i_impl_det )
case ( i_impl_det_indep )
  ! Solve detrainment for each convective layer / type independently.

  ! So set effective mass-flux to just the current layer / type
  do ic = 1, n_points
    sum_massflux_det(ic) = par_prev_super(ic,i_massflux_d)
  end do

case ( i_impl_det_inter )
  ! Solve detrainment accounting for the subsidence heating from other
  ! layers / types, assuming they have the same fractional mass-flux
  ! change over the level-step as the current one.

  ! So set effective mass-flux to the sum over all layers / types
  do ic = 1, n_points
    sum_massflux_det(ic) = sum_massflux(ic)
  end do

end select  ! i_impl_det


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

if ( l_to_full_level ) then
  ! If this is the first of the two half-level-steps
  ! (from previous half-level to level k),
  ! then the environment fields to entrain are at level k, but the
  ! parcel is at a different pressure, at the previous half-level.
  ! Adjust the environment temperature to what it would be at the
  ! start of the level-step, so that we entrain it into the parcel
  ! consistently...
  do ic = 1, n_points
    exner_ratio(ic) = one
  end do
  call dry_adiabat( n_points, n_points,                                        &
                  grid_next_super(:,i_pressure), grid_prev_super(:,i_pressure),&
                    ent_fields(:,i_q_vap),                                     &
                    ent_fields(:,i_qc_first:i_qc_last),                        &
                    exner_ratio )
  do ic = 1, n_points
    ent_fields(ic,i_temperature) = ent_fields(ic,i_temperature)                &
                                 * exner_ratio(ic)
  end do
end if

! Set amount of entrained dry-mass over the current half-level-step
call set_ent( n_points, n_fields_tot, max_points,                              &
              max_ent_frac,                                                    &
              par_conv_mean_fields, ent_fields,                                &
              grid_prev_super, grid_next_super,                                &
              par_conv_super,                                                  &
              l_within_bl, core_mean_ratio,                                    &
              layer_mass_step, sum_massflux,                                   &
              ent_mass_d, core_ent_ratio )

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
  do i_field = 1, n_fields_tot
    do ic = 1, n_points
      det_fields(ic,i_field)                                                   &
        =      core_ent_ratio(ic)  * ent_fields(ic,i_field)                    &
        + (one-core_ent_ratio(ic)) * par_conv_mean_fields(ic,i_field)
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

! Core environment entrainment ratio diagnostic
if ( plume_model_diags % core_ent_ratio % flag ) then
  i_diag = plume_model_diags % core_ent_ratio % i_super
  do ic = 1, n_points
    diags_super(ic,i_diag) = core_ent_ratio(ic)
  end do
end if

! Save diagnostics of entrained mass and air properties
! (in conserved variable form ready for finding mean over types)
if ( plume_model_diags % ent_mass_d % flag ) then
  i_diag = plume_model_diags % ent_mass_d % i_super
  do ic = 1, n_points
    diags_super(ic,i_diag) = ent_mass_d(ic)
  end do
end if
if ( plume_model_diags % ent_fields % n_diags > 0 ) then
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

! Initialise a super-array which stores various fields at sub-level steps
! between the previous and next model-level interfaces
call init_sublevs ( n_points, max_points, l_down,                              &
                    grid_prev_super(:,i_height), grid_next_super(:,i_height),  &
                    env_prev_super, env_next_super,                            &
                    par_mean_virt_temp, par_core_virt_temp,                    &
                    par_prev_mean_fields(:,i_wind_w),                          &
                    par_prev_core_fields(:,i_wind_w),                          &
                    par_prev_super, par_conv_super,                            &
                    delta_tv, sum_massflux_det,                                &
                    sublevs, i_next, i_sat, i_core_sat )

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
                   l_down, l_tracer,                                           &
                   cmpr, k, draft_string, i_call,                              &
                   grid_prev_super, grid_next_super,                           &
                   l_within_bl, layer_mass_step, sum_massflux,                 &
                   par_conv_super(:,i_massflux_d), par_conv_super(:,i_radius), &
                   env_k_fields, Nsq_dry,                                      &
                   par_prev_core_fields(:,i_temperature:n_fields),             &
                   par_conv_core_fields,                                       &
                   sublevs, i_next, i_core_sat, par_core_w_drag,               &
                   par_prev_core_ss, par_next_core_ss,                         &
                   par_prev_core_tvl, par_next_core_tvl )

end if  ! ( l_par_core )

! Call parcel dynamics routine; does dry-lifting up to
! next level, moist processes, and momentum equations
i_call = i_call_mean
! Note: setting i_call_mean enables calculation of both
! resolved-scale sources terms and diagnostics in parcel_dyn.
call parcel_dyn( n_points, n_points, max_points,                               &
                 max_points, n_points_res, n_points,                           &
                 n_points_diag, n_diags_super, n_fields_tot,                   &
                 l_down, l_tracer,                                             &
                 cmpr, k, draft_string, i_call,                                &
                 grid_prev_super, grid_next_super,                             &
                 l_within_bl, layer_mass_step, sum_massflux,                   &
                 par_conv_super(:,i_massflux_d), par_conv_super(:,i_radius),   &
                 env_k_fields, Nsq_dry,                                        &
                 par_prev_mean_fields(:,i_temperature:n_fields),               &
                 par_conv_mean_fields,                                         &
                 sublevs, i_next, i_sat, par_mean_w_drag,                      &
                 par_prev_mean_ss, par_next_mean_ss,                           &
                 par_prev_mean_tvl, par_next_mean_tvl,                         &
                 res_source_fields = res_source_fields,                        &
                 plume_model_diags = plume_model_diags,                        &
                 diags_super       = diags_super,                              &
                 i_core_sat        = i_core_sat,                               &
                 next_core_fields  = par_conv_core_fields,                     &
                 core_mean_ratio   = core_mean_ratio )


!---------------------------------------------------------------
! 4) Calculate detrainment, based on an assumed-PDF of buoyancy
!---------------------------------------------------------------

! Find max sub-level address that we need to loop over in the sublevs array
i_next_max = maxval(i_next)

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

! Call routine to set detrained mass and detrained air properties,
! and update the mass-flux and parcel mean properties due to detrainment
call set_det( n_points, max_points, n_points_res, n_fields_tot,                &
              l_down, l_last_level, l_to_full_level,                           &
              cmpr, k, draft_string,                                           &
              sublevs, i_next, i_next_max,                                     &
              core_mean_ratio,                                                 &
              par_mean_w_drag, par_core_w_drag, res_source_fields,             &
              par_conv_mean_fields, par_conv_core_fields,                      &
              par_prev_super, par_conv_super,                                  &
              det_mass_d, det_fields )

! TEMPORARY CODE TO PRESERVE KGO; TO BE REMOVED SOON:
! (should be done more accurately inside set_det but changes answers).
! Store the final-guess total environment virtual temperature
! increment in delta_tv_k (just needs to be scaled down by the
! final non-detrained fraction)
do ic = 1, n_points
  ! Scale the value at prev (not altered by set_det)
  sublevs(ic,j_delta_tv,i_prev) = sublevs(ic,j_delta_tv,i_prev)                &
                 * par_conv_super(ic,i_massflux_d)                             &
            / max( par_conv_super(ic,i_massflux_d)+det_mass_d(ic), min_float )
  ! Overwrite other sub-levels with the constant-with-height value from prev
  do i_lev = i_prev+1, i_next(ic)
    sublevs(ic,j_delta_tv,i_lev) = sublevs(ic,j_delta_tv,i_prev)
  end do
end do

! Subtract the environment subsidence Tv increment from both the mean and
! the core buoyancies, to obtain final implicitly-solved buoyancies
do i_lev = i_prev, i_next_max
  do ic = 1, n_points
    sublevs(ic,j_mean_buoy,i_lev) = sublevs(ic,j_mean_buoy,i_lev)              &
                                  - sublevs(ic,j_delta_tv,i_lev)
  end do
end do
if ( l_par_core ) then
  do i_lev = i_prev, i_next_max
    do ic = 1, n_points
      sublevs(ic,j_core_buoy,i_lev) = sublevs(ic,j_core_buoy,i_lev)            &
                                    - sublevs(ic,j_delta_tv,i_lev)
    end do
  end do
end if

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

  ! Copy the existing detrained fields for use as the
  ! "prev" fields passed into parcel_dyn
  ! (using the ent_fields array as work-space)
  do i_field = 1, n_fields_tot
    do ic = 1, n_points
      ent_fields(ic,i_field) = det_fields(ic,i_field)
    end do
  end do

  ! Call parcel dynamics routine to adjust parcel back to level k.
  i_call = i_call_det
  ! Note: setting i_call_det suppresses the updating of the sub-levels
  ! super-array inside parcel_dyn.
  call parcel_dyn( n_points, n_points, n_points,                               &
                   max_points, n_points_res, n_points,                         &
                   1, 1, n_fields_tot,                                         &
                   (.not. l_down), l_tracer,                                   &
                   cmpr, k, draft_string, i_call,                              &
                   grid_next_super, grid_prev_super,                           &
                   l_within_bl, layer_mass_step, sum_massflux,                 &
                   det_mass_d, par_conv_super(:,i_radius),                     &
                   env_k_fields, Nsq_dry,                                      &
                   ent_fields(:,i_temperature:n_fields), det_fields,           &
                   sublevs, i_next, i_sat, par_mean_w_drag,                    &
                   par_prev_mean_ss, par_next_mean_ss,                         &
                   par_prev_mean_tvl, par_next_mean_tvl,                       &
                   res_source_fields = res_source_fields )

  ! Convert detrained fields back into conserved variable form
  l_reverse = .false.
  call fields_k_conserved_vars( n_points, n_points,                            &
                                n_fields_tot, l_reverse,                       &
                                det_fields )

end if  ! ( .NOT. l_to_full_level )

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

! TEMPORARY CODE TO PRESERVE KGO; TO BE REMOVED SOON:
! (should be done more accurately inside set_det and should have been done
!  at saturation height as well as end-of-step, but changes answers).
! Recalculate the parcel mean buoyancy at the end of the level-step,
! accounting for modification of the parcel mean properties by detrainment.
call calc_virt_temp( n_points, max_points,                                     &
                     par_conv_mean_fields(:,i_temperature),                    &
                     par_conv_mean_fields(:,i_q_vap),                          &
                     par_conv_mean_fields(:,i_qc_first:i_qc_last),             &
                     par_mean_virt_temp )
do ic = 1, n_points
  sublevs(ic,j_mean_buoy,i_next(ic)) = par_mean_virt_temp(ic)                  &
    - ( sublevs(ic,j_env_tv,i_next(ic)) + sublevs(ic,j_delta_tv,i_next(ic)) )
  ! Also reset massflux at next
  sublevs(ic,j_massflux_d,i_next(ic)) = par_conv_super(ic,i_massflux_d)
end do

! Update edge virtual temperature
call update_edge_virt_temp( n_points, max_points, l_down,                      &
                            core_mean_ratio,                                   &
                            sublevs, i_next, par_conv_super )

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
  if ( i_cf_conv == i_cf_conv_coma ) then
    ! Use old convective cloud amount calculation used in comorph A
    call set_diag_conv_cloud_a(                                                &
         n_points, max_points, n_points_res,                                   &
         n_points, max_points, n_fields_tot,                                   &
         grid_prev_super, grid_next_super, frac_level_step,                    &
         env_prev_super(:,i_virt_temp),  env_next_super(:,i_virt_temp),        &
         par_prev_super(:,i_massflux_d), par_conv_super(:,i_massflux_d),       &
         par_prev_super(:,i_radius),     par_conv_super(:,i_radius),           &
         par_prev_mean_fields(:,i_temperature:n_fields), par_conv_mean_fields, &
         sublevs, i_next, i_sat,                                               &
         convcloud_super )
  else
    ! Use other options...
    ! To be inserted here.
  end if
end if

if ( l_calc_cape .or. l_calc_mfw_cape ) then
  ! Calculate contribution to CAPE from current half-level-step
  call calc_cape(                                                              &
         n_points, sublevs, i_next, l_within_bl,                               &
         index_ij, ij_first, ij_last,                                          &
         fields_2d )
end if


return
end subroutine conv_level_step

end module conv_level_step_mod
