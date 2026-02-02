! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module parcel_dyn_mod

implicit none

! Indicators for which call to parcel_dyn this is
! (the input i_call will have one of these values;
!  certain behaviours of parcel_dyn need to differ between
!  the calls).

! Test ascent in the initiation mass-source calculation
integer, parameter :: i_call_genesis = 1

! Lifting of parcel mean properties
integer, parameter :: i_call_mean = 2

! Lifting of parcel core properties
integer, parameter :: i_call_core = 3

! Adjusting of detrained air to level k
integer, parameter :: i_call_det = 4

contains

!----------------------------------------------------------------
! Subroutine does moist physics and dynamical processes for a
! parcel lifting
!----------------------------------------------------------------
subroutine parcel_dyn( n_points, n_points_prev, n_points_next,                 &
                       n_points_env, n_points_res, n_points_sublevs,           &
                       n_points_diag, n_diags_super, n_fields_tot,             &
                       l_down, l_tracer,                                       &
                       cmpr, k, draft_string, i_call,                          &
                       grid_prev, grid_next,                                   &
                       l_within_bl, layer_mass_step, sum_massflux,             &
                       massflux_d, par_radius,                                 &
                       env_k_fields, Nsq_dry,                                  &
                       par_prev_fields, par_next_fields,                       &
                       sublevs, i_next, i_sat, par_w_drag,                     &
                       prev_ss, next_ss, prev_tvl, next_tvl,                   &
                       res_source_fields,                                      &
                       plume_model_diags, diags_super,                         &
                       i_core_sat, next_core_fields,                           &
                       core_mean_ratio )

use comorph_constants_mod, only: real_cvprec, newline,                         &
                                 zero, one, n_cond_species, n_tracers,         &
                                 i_check_bad_values_cmpr, i_check_bad_none,    &
                                 name_length,                                  &
                                 wind_w_fac, w_min,                            &
                                 l_par_core, l_cv_cf, l_cv_cloudfrac,          &
                                 l_tracer_scav,                                &
                                 par_vert_len_fac,                             &
                                 i_cond_cl, i_cond_cf,                         &
                                 i_mean_q_cl, i_mean_q_cl_hom
use grid_type_mod, only: n_grid, i_height, i_pressure
use fields_type_mod, only: i_temperature, i_q_vap, i_q_cl, i_q_cf,             &
                           i_qc_first, i_qc_last,                              &
                           i_cf_liq, i_cf_bulk,                                &
                           i_wind_u, i_wind_v, i_wind_w, i_tracers,            &
                           n_fields, field_names, field_positive
use cmpr_type_mod, only: cmpr_type

use linear_qs_mod, only: linear_qs_set_ref,                                    &
                         n_linear_qs_fields, i_ref_temp
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev,                     &
                       j_mean_buoy, j_core_buoy, j_env_tv,                     &
                       j_mean_wex, j_core_wex
use plume_model_diags_type_mod, only: plume_model_diags_type
use moist_proc_diags_type_mod, only: moist_proc_diags_type

use calc_rho_dry_mod, only: calc_rho_dry
use dry_adiabat_mod, only: dry_adiabat
use precip_res_source_mod, only: precip_res_source
use moist_proc_mod, only: moist_proc
use sat_adjust_mod, only: sat_adjust
use set_cp_tot_mod, only: set_cp_tot
use calc_q_tot_mod, only: calc_q_tot
use set_par_cloudfrac_mod, only: set_par_cloudfrac
use calc_qvl_supersat_mod, only: calc_qvl_supersat
use calc_mean_q_cl_with_core_mod, only: calc_mean_q_cl_with_core
use tracer_source_mod, only: tracer_source
use calc_sat_height_mod, only: calc_sat_height
use calc_virt_temp_mod, only: calc_virt_temp
use momentum_eqn_mod, only: momentum_eqn
use check_bad_values_mod, only: check_bad_values_cmpr
use raise_error_mod, only: raise_fatal

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points in the parcel and environment super-arrays
! (maybe larger than needed here, to save having to reallocate)
integer, intent(in) :: n_points_prev
integer, intent(in) :: n_points_next
integer, intent(in) :: n_points_env
integer, intent(in) :: n_points_res
integer, intent(in) :: n_points_sublevs
integer, intent(in) :: n_points_diag
integer, intent(in) :: n_diags_super

! Number of primary fields stored in par_next_fields
! (may or may not include tracers)
integer, intent(in) :: n_fields_tot

! Flag for downdrafts versus updrafts
logical, intent(in) :: l_down

! Flag for whether tracers are included in the primary fields arrays
logical, intent(in) :: l_tracer

! Stuff used to make error messages more informative:
! Structure storing compression indices of the current points
type(cmpr_type), intent(in) :: cmpr
! Current model-level index
integer, intent(in) :: k
! String identifying what sort of draft (updraft, downdraft, etc)
character(len=name_length), intent(in) :: draft_string

! Integer indicator for which call to parcel_dyn this is
integer, intent(in) :: i_call

! Height and pressure...
!   At previous level
real(kind=real_cvprec), intent(in) :: grid_prev                                &
                                      ( n_points_env, n_grid )
!   At next level which we are lifting up to
real(kind=real_cvprec), intent(in) :: grid_next                                &
                                      ( n_points_env, n_grid )

! Flag for whether each point is within the homogenised mixed-layer
logical, intent(in) :: l_within_bl(n_points)

! Dry-mass on the current model-level step
real(kind=real_cvprec), intent(in) :: layer_mass_step(n_points)

! Sum of parcel mass-fluxes over types / layers
real(kind=real_cvprec), intent(in) :: sum_massflux(n_points)

! Dry-mass flux of the parcel (for computing resolved-scale
! source terms)
real(kind=real_cvprec), intent(in) :: massflux_d(n_points)

! Radius length-scale of the parcel
real(kind=real_cvprec), intent(in) :: par_radius(n_points)

! Super-array containing environment fields at k; properties
! of source air for precip which falls into the parcel
real(kind=real_cvprec), intent(in) :: env_k_fields                             &
                                      ( n_points_env, n_fields )

! Environment dry static stability
real(kind=real_cvprec), intent(in) :: Nsq_dry(n_points)

! Super-array containing parcel fields at "previous" level
! These are valid at the previous model-level interface
! (except in the genesis test-ascent where they are at k)
real(kind=real_cvprec), intent(in) :: par_prev_fields                          &
                                      ( n_points_prev, i_temperature:n_fields )

! Super-array containing the parcel fields to be updated.
! On input they are consistent with the pressure at the start of
! the current level-step.
! On output they are valid at the end of the level-step.
real(kind=real_cvprec), intent(in out) :: par_next_fields                      &
                                ( n_points_next, n_fields_tot )

! Super-array storing the parcel buoyancies and other properties
! at up to 4 sub-level heights within the current level-step:
! a) Start of the level-step
! b) End of the level-step
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(in out) :: sublevs                              &
                     ( n_points_sublevs, n_sublev_vars, max_sublevs )

! Address of next model-level interface in sublevs
integer, intent(in out) :: i_next(n_points)
! Address of saturation height in sublevs
integer, intent(in out) :: i_sat(n_points)

! Drag on the parcel vertical velocity / s-1
! Output from momentum_eqn for use in correcting w due to changes in buoyancy
real(kind=real_cvprec), intent(out) :: par_w_drag(n_points)

! Supersaturation the parcel would have with all liquid cloud evaporated
real(kind=real_cvprec), intent(out) :: prev_ss(n_points)
real(kind=real_cvprec), intent(out) :: next_ss(n_points)

! Virtual temperature the parcel would have with all liquid cloud evaporated
real(kind=real_cvprec), intent(out) :: prev_tvl(n_points)
real(kind=real_cvprec), intent(out) :: next_tvl(n_points)

! Super-array containing resolved-scale source terms for each primary field
real(kind=real_cvprec), optional, intent(in out) :: res_source_fields          &
                                                ( n_points_res, n_fields_tot )

! Structure containing diagnostic flags and meta-data
type(plume_model_diags_type), optional, intent(in) :: plume_model_diags

! Super-array to store the diagnostics
real(kind=real_cvprec), optional, intent(in out) :: diags_super                &
                                              ( n_points_diag, n_diags_super )

! Parcel core properties which may influence the parcel mean
! (only passed in when this is the mean parcel-lifting):

! Address of separate core saturation height in buoyancy_super
! (only passed in for separate means ascents modified by core)
integer, optional, intent(in out) :: i_core_sat(n_points)

! Core liquid-cloud mixing-ratio
real(kind=real_cvprec), optional, intent(in) :: next_core_fields               &
                                                ( n_points_next, n_fields_tot )

! Ratio of core buoyancy over mean buoyancy
real(kind=real_cvprec), optional, intent(in) :: core_mean_ratio(n_points)


! Description of which call to parcel_dyn this is
character(len=name_length) :: call_string

! Air density
real(kind=real_cvprec) :: rho_dry(n_points)

! Height interval for the parcel lifting
real(kind=real_cvprec) :: delta_z(n_points)

! Time interval for the parcel lifting
real(kind=real_cvprec) :: delta_t(n_points)

! Mean vertical velocity excess over the level-step
real(kind=real_cvprec) :: wind_w_excess(n_points)

! Parcel vertical length-scale, set to twice the radius
real(kind=real_cvprec) :: vert_len(n_points)

! Factor dt / ( rho_dry lz ), used for converting between
! precip fall-fluxes and in-parcel mixing-ratio increments
real(kind=real_cvprec) :: dt_over_rhod_lz(n_points)

! Super-array containing qsat at a reference temperature,
! and dqsat/dT, for linearised qsat calculations
real(kind=real_cvprec) :: linear_qs_super                                      &
                            ( n_points, n_linear_qs_fields )

! Super-array to store precip fall fluxes
real(kind=real_cvprec) :: flux_cond( n_points, n_cond_species )

! Total heat capacity for phase-change
real(kind=real_cvprec) :: cp_tot(n_points)
! Total water mixing ratio
real(kind=real_cvprec) :: q_tot(n_points)

! Parcel virtual temperature at next level
real(kind=real_cvprec) :: par_next_virt_temp(n_points)

! Total rate of conversion of cloud-water to precip,
! for calculating tracer scavenging
real(kind=real_cvprec) :: dq_prec(n_points)

! Description of where we are in the code, for error messages
character(len=name_length) :: where_string

! Flag for whether to calculate resolved-scale source terms
! from non-dry-mass-exchanging processes (precipitation fall
! and pressure forces on momentum)
logical :: l_res_source

! Master switch for whether or not to calculate any diagnostics
! (e.g. parcel_dyn might be called within an iteration, where
!  we only want to output diags for the final call).
logical :: l_diags

! Flag passed into precip_res_source indicating whether fall-flux is
! from environment to parcel or vice-versa
logical :: l_fall_in

! Flag passed into sat_adjust to tell it not to update
! q_vap and q_cl when calculating a saturated reference T
logical, parameter :: l_update_q_false = .false.

! Flag passed into check_bad_values to indicate whether fields must be positive
logical :: l_positive

! Flag for parcel mean ascent with an accompanying core
logical :: l_mean_with_core

! Dummy diagnostics structure to pass into moist_proc in the
! case where diagnostics are not requested
type(moist_proc_diags_type) :: moist_proc_diags_dummy
! Dummy diagnostics array
real(kind=real_cvprec) :: diags_super_dummy(1,1)

! Loop counters
integer :: ic, i_cond, i_diag, i_field, j_buoy, j_wex

character(len=*), parameter :: routinename = "PARCEL_DYN"


! Set character string description for which call to parcel_dyn
! this is, for use in error-messages.
select case(i_call)
case (i_call_genesis)
  call_string = trim(adjustl(draft_string)) // ", genesis"
case (i_call_mean)
  call_string = trim(adjustl(draft_string)) // ", par mean lifting"
case (i_call_core)
  call_string = trim(adjustl(draft_string)) // ", par core lifting"
case (i_call_det)
  call_string = trim(adjustl(draft_string)) // ", detrained air"
end select

! Set flag for whether to calculate resolved-scale source terms
! due to non-mass-exchanging processes.  Only needs to be done
! for the lifting of the mean parcel and the detrained
! air; the core ascent does not contribute
l_res_source = ( i_call == i_call_mean .or. i_call == i_call_det )

! Set flag for whether to save diagnostics.  These are currently
! only ouput for the parcel mean ascent
l_diags = ( i_call == i_call_mean )

! Set flag for a parcel mean ascent with an accompanying core
! (i.e. if this is the call for the mean ascent, and required optional
!  arguments are present)
l_mean_with_core = ( i_call == i_call_mean .and. l_par_core                    &
                     .and. present(next_core_fields)                           &
                     .and. present(core_mean_ratio) )

! Check the required optional arguments for res_source and diag
! calculations are present if the flags are true.
if ( l_res_source .and. ( .not. present(res_source_fields) ) ) then
  call raise_fatal( routinename,                                               &
         "Flag l_res_source has been set to true to calculate "             // &
         "resolved-scale source terms, but the required optional " //newline// &
         "argument res_source_fields is not present.  "                      //&
         "For call: " // call_string )
end if
if ( l_diags .and. ( .not. ( present(plume_model_diags) .and.                  &
                             present(diags_super) ) ) ) then
  call raise_fatal(routinename,                                                &
         "Flag l_diags has been set to true to calculate "                  // &
         "diagnostics, but the required optional "                 //newline// &
         "arguments plume_model_diags and diags_super "                     // &
         "are not present.  For call: " // call_string )
end if


! Check inputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  ! Check parcel properties
  where_string = "Start of parcel_dyn call for "                            // &
                 trim(adjustl(call_string))   // "; "                       // &
                 "par_next_fields"
  do i_field = 1, n_fields_tot
    call check_bad_values_cmpr( cmpr, k,                                       &
                                par_next_fields(:,i_field),                    &
                                where_string,                                  &
                                field_names(i_field),                          &
                                field_positive(i_field) )
  end do
  if ( present( res_source_fields ) ) then
    ! Check resolved-scale source-terms
    where_string = "Start of parcel_dyn call for "                          // &
                   trim(adjustl(call_string))   // "; "                     // &
                   "res_source_fields"
    l_positive = .false.  ! Source terms may be positive or negative
    do i_field = 1, n_fields_tot
      call check_bad_values_cmpr( cmpr, k,                                     &
                                  res_source_fields(:,i_field),                &
                                  where_string,                                &
                                  field_names(i_field),                        &
                                  l_positive )
    end do
  end if
end if


!----------------------------------------------------------------
! 1) Set grid-related fields needed in the calculations
!----------------------------------------------------------------

! Select address for core or mean properties in the sub-levels super-array
if ( i_call == i_call_core ) then
  j_buoy = j_core_buoy
  j_wex  = j_core_wex
else
  j_buoy = j_mean_buoy
  j_wex  = j_mean_wex
end if

! Set height interval
! Note: for downdrafts, this will be negative!
do ic = 1, n_points
  delta_z(ic) = grid_next(ic,i_height) - grid_prev(ic,i_height)
end do

! Calculate dry density of the air (using value at start of current level-step)
call calc_rho_dry( n_points,                                                   &
                   par_prev_fields(:,i_temperature),                           &
                   par_prev_fields(:,i_q_vap),                                 &
                   grid_prev(:,i_pressure),                                    &
                   rho_dry )

! For now, set vertical velocity over the level-step to a
! prescribed value (sign depends on whether updraft or downdraft).
!
! Constant wind_w_fac is the draft speed in m s-1
do ic = 1, n_points
  wind_w_excess(ic) = sign( wind_w_fac, delta_z(ic) )
end do

! Set time interval dt = dz / w
do ic = 1, n_points
  delta_t(ic) = delta_z(ic) / wind_w_excess(ic)
end do

! Set vertical length-scale of parcel used for precip fall;
! this should be parcel volume over cross-sectional area = 4/3 R,
! but can be set to other values, e.g. if updrafts aren't spherical.
do ic = 1, n_points
  vert_len(ic) = par_vert_len_fac * par_radius(ic)
end do

! Calculate factor dt / (rho_dry lz), used for converting
! precip fall between flux and in-parcel mixing-ratio increment
do ic = 1, n_points
  dt_over_rhod_lz(ic) = delta_t(ic)                                            &
                      / ( rho_dry(ic) * vert_len(ic) )
end do


!----------------------------------------------------------------
! 2) Perform dry-lifting of the parcel to new pressure
!----------------------------------------------------------------

call dry_adiabat( n_points, n_points_next,                                     &
                  grid_prev(:,i_pressure), grid_next(:,i_pressure),            &
                  par_next_fields(:,i_q_vap),                                  &
                  par_next_fields(:,i_qc_first:i_qc_last),                     &
                  par_next_fields(:,i_temperature) )


!----------------------------------------------------------------
! 3) Do moist physical processes
!----------------------------------------------------------------

! Set reference temperature and compute qsat and dqsat/dT...
! Initialise ref temp equal to parcel temp before
! phase-changes.
do ic = 1, n_points
  linear_qs_super(ic,i_ref_temp)                                               &
    = par_next_fields(ic,i_temperature)
end do
! Adjust the reference temperature to saturation
call set_cp_tot( n_points, n_points_next,                                      &
                 par_next_fields(:,i_q_vap),                                   &
                 par_next_fields(:,i_qc_first:i_qc_last), cp_tot )
call sat_adjust( n_points, l_update_q_false,                                   &
                 grid_next(:,i_pressure),                                      &
                 linear_qs_super(:,i_ref_temp),                                &
                 par_next_fields(:,i_q_vap),                                   &
                 par_next_fields(:,i_q_cl), cp_tot )

! Call routine to compute qsat and dqsat/dT at the final pressure
call linear_qs_set_ref( n_points,                                              &
                        grid_next(:,i_pressure), linear_qs_super)


! FALL-IN OF PRECIP FROM THE ENVIRONMENT INTO THE PARCEL...

! The fall-in flux will be calculated here, but this is not
! yet implemented.
! For now, set fall-in flux of precip to zero
! (it is needed as an input to moist_proc below).
do i_cond = 1, n_cond_species
  do ic = 1, n_points
    flux_cond(ic,i_cond) = zero
  end do
end do

! Calculate resolved-scale source terms from fall-in
! of precip from the environment into the parcel
if ( l_res_source .and. .false. ) then
  ! Note: this block of code is disabled pending implementation
  ! of a calculation of precip fall-in above.

  ! Calculate resolved-scale source terms from fall-in
  ! of precip from the environment into the parcel
  l_fall_in = .true.
  call precip_res_source( n_points, n_points_env, n_points_res,                &
                          l_fall_in, dt_over_rhod_lz, massflux_d,              &
                          flux_cond, env_k_fields,                             &
                          res_source_fields )

end if  ! ( l_res_source )


if ( l_tracer .and. l_tracer_scav ) then
  ! Store pre-existing + falling-in cloud-condensate + vapour.
  ! This will be used later to compute the total conversion of
  ! water into precip inside the parcel.
  ! Note that flux_cond is defined positive for a downwards flux of precip,
  ! so scaling it by dt/(rhod lz) yields the positive increment to in-parcel
  ! qc that we will get from precip falling into the parcel.
  if ( l_cv_cf ) then
    do ic = 1, n_points
      dq_prec(ic) = par_next_fields(ic,i_q_vap)                                &
                  + par_next_fields(ic,i_q_cl) + par_next_fields(ic,i_q_cf)    &
                  + ( flux_cond(ic,i_cond_cl) + flux_cond(ic,i_cond_cf) )      &
                    * dt_over_rhod_lz(ic)
    end do
  else
    do ic = 1, n_points
      dq_prec(ic) = par_next_fields(ic,i_q_vap) + par_next_fields(ic,i_q_cl)   &
                  + flux_cond(ic,i_cond_cl) * dt_over_rhod_lz(ic)
    end do
  end if
end if

! Call moist process routine to update temperature and
! water species mixing ratios due to condensation,
! evaporation and other microphysical processes.
! Call with real diagnostics passed out if doing diags,
! otherwise call with dummy diagnostics
if ( l_diags ) then
  call moist_proc( n_points, n_points_next,  linear_qs_super,                  &
                   delta_z, delta_t, vert_len, wind_w_excess,                  &
                   dt_over_rhod_lz, grid_next(:,i_pressure),                   &
                   par_prev_fields(:,i_temperature),                           &
                   env_k_fields(:,i_wind_u),                                   &
                   env_k_fields(:,i_wind_v),                                   &
                   env_k_fields(:,i_wind_w),                                   &
                   env_k_fields(:,i_temperature),                              &
                   par_next_fields(:,i_wind_u),                                &
                   par_next_fields(:,i_wind_v),                                &
                   par_next_fields(:,i_wind_w),                                &
                   par_next_fields(:,i_temperature),                           &
                   par_next_fields(:,i_q_vap),                                 &
                   par_next_fields(:,i_qc_first:i_qc_last),                    &
                   flux_cond, cmpr, k, call_string, l_diags,                   &
                   plume_model_diags % moist_proc,                             &
                   n_points_diag, n_diags_super, diags_super )
else
  call moist_proc( n_points, n_points_next, linear_qs_super,                   &
                   delta_z, delta_t, vert_len, wind_w_excess,                  &
                   dt_over_rhod_lz, grid_next(:,i_pressure),                   &
                   par_prev_fields(:,i_temperature),                           &
                   env_k_fields(:,i_wind_u),                                   &
                   env_k_fields(:,i_wind_v),                                   &
                   env_k_fields(:,i_wind_w),                                   &
                   env_k_fields(:,i_temperature),                              &
                   par_next_fields(:,i_wind_u),                                &
                   par_next_fields(:,i_wind_v),                                &
                   par_next_fields(:,i_wind_w),                                &
                   par_next_fields(:,i_temperature),                           &
                   par_next_fields(:,i_q_vap),                                 &
                   par_next_fields(:,i_qc_first:i_qc_last),                    &
                   flux_cond, cmpr, k, call_string, l_diags,                   &
                   moist_proc_diags_dummy,                                     &
                   1, 1, diags_super_dummy )
end if


if ( l_tracer .and. l_tracer_scav ) then
  ! Calculate final + falling-out cloud-condensate + vapour,
  ! and subtract this from the preexisting + falling-in cloud condensate
  ! saved earlier to compute the total precip production inside
  ! the parcel.
  if ( l_cv_cf ) then
    do ic = 1, n_points
      dq_prec(ic) = dq_prec(ic)                                                &
                  - ( par_next_fields(ic,i_q_vap)                              &
                    + par_next_fields(ic,i_q_cl) + par_next_fields(ic,i_q_cf)  &
                    + ( flux_cond(ic,i_cond_cl) + flux_cond(ic,i_cond_cf) )    &
                       * dt_over_rhod_lz(ic) )
    end do
  else
    do ic = 1, n_points
      dq_prec(ic) = dq_prec(ic)                                                &
                  - ( par_next_fields(ic,i_q_vap) + par_next_fields(ic,i_q_cl) &
                    + flux_cond(ic,i_cond_cl) * dt_over_rhod_lz(ic) )
    end do
  end if
end if

if ( i_call == i_call_genesis ) then

  ! For the genesis test ascents, don't let the precip fall out
  ! of the parcel.  This is a temporary fix to avoid
  ! spurious convective triggering from columns with large
  ! rain mixing-ratios (the parcel becomes buoyant at the
  ! next level just because the rain falls out, so it loses
  ! the water loading; in reality, rain should be falling
  ! in as fast as it is falling out, but fall-in isn't
  ! implemented yet).
  do i_cond = 1, n_cond_species
    i_field = i_qc_first-1 + i_cond
    do ic = 1, n_points
      ! Add the fall-out flux back onto the parcel
      par_next_fields(ic,i_field) = par_next_fields(ic,i_field)                &
        + flux_cond(ic,i_cond) * dt_over_rhod_lz(ic)
    end do
  end do

end if  ! ( i_call == i_call_genesis )

! Update in-parcel cloud-fractions and precipitation fraction
if ( l_cv_cloudfrac ) then
  call set_par_cloudfrac( n_points, n_points_next,                             &
                          par_next_fields(:,i_q_cl),                           &
                          par_next_fields(:,i_q_cf),                           &
                          par_next_fields(:,i_cf_liq:i_cf_bulk) )
end if

! If using separate parcel mean and core properties,
! attempt to approximately account for effect of sub-plume
! moisture variability on condensation.  At the moment, this
! just entails retrospectively adding more condensation to
! the parcel mean if the core has cloud in it but the mean
! has disproportionately less (the mean is supposed to include the core).
if ( i_mean_q_cl > i_mean_q_cl_hom .and. l_mean_with_core ) then
  call calc_mean_q_cl_with_core( n_points, n_points_next, n_fields_tot,        &
                                 next_core_fields(:,i_q_cl),                   &
                                 core_mean_ratio, par_next_fields )
end if

if ( l_res_source ) then

  ! Add on resolved-scale source term contributions from
  ! precipitation fall-out
  l_fall_in = .false.
  call precip_res_source( n_points, n_points_next, n_points_res,               &
                          l_fall_in, dt_over_rhod_lz, massflux_d,              &
                          flux_cond, par_next_fields,                          &
                          res_source_fields )

end if

! Apply any in-parcel source terms for tracers
! (e.g. in-plume scavenging of aerosols by precipitation)
if ( l_tracer ) then
  call tracer_source( n_points, n_points_next,                                 &
                      massflux_d, dq_prec, par_next_fields(:,1:n_fields),      &
                      par_next_fields(:,i_tracers(1):i_tracers(n_tracers)) )
end if


!----------------------------------------------------------------
! 4) Find parcel virtual temperature excess used for detrainment
!----------------------------------------------------------------

if ( .not. i_call == i_call_det ) then
  ! Don't need to do this part when just adjusting the detrained air
  ! to level k...

  ! Calculate supersaturation and virtual temperature the parcel would have
  ! if all liquid-cloud evaporated, at both prev and next level...
  call calc_qvl_supersat( n_points, n_points_prev,                             &
                          grid_prev(:,i_pressure),                             &
                          par_prev_fields(:,i_temperature),                    &
                          par_prev_fields(:,i_q_vap),                          &
                          par_prev_fields(:,i_qc_first:i_qc_last),             &
                          prev_ss, prev_tvl )
  call calc_qvl_supersat( n_points, n_points_next,                             &
                          grid_next(:,i_pressure),                             &
                          par_next_fields(:,i_temperature),                    &
                          par_next_fields(:,i_q_vap),                          &
                          par_next_fields(:,i_qc_first:i_qc_last),             &
                          next_ss, next_tvl,                                   &
                          linear_qs=linear_qs_super )

  ! Find new parcel virtual temperature and excess
  call calc_virt_temp( n_points, n_points_next,                                &
                       par_next_fields(:,i_temperature),                       &
                       par_next_fields(:,i_q_vap),                             &
                       par_next_fields(:,i_qc_first:i_qc_last),                &
                       par_next_virt_temp )

  ! Store virtual temperature excess at the
  ! next model-level interface
  do ic = 1, n_points
    sublevs(ic,j_buoy,i_next(ic)) = par_next_virt_temp(ic)                     &
                                  - sublevs(ic,j_env_tv,i_next(ic))
  end do

  ! At the level where the parcel first hits saturation, need to
  ! find the buoyancy at the point of saturation, where there
  ! will be a kink in the parcel buoyancy profile in-between
  ! the model-levels.  This needs to be accounted for in the
  ! detrainment calculation, and also in the vertical momentum
  ! equation.
  ! In updrafts, this usually corresponds to the minimum buoyancy of the
  ! parcel as it rises through the BL-top inversion, and is
  ! therefore critical for getting the detrainment below
  ! cloud-base (and hence the mass-flux at cloud-base) right.
  call calc_sat_height(                                                        &
         n_points, n_points_sublevs, l_mean_with_core, l_down, j_buoy,         &
         prev_ss, next_ss, prev_tvl, next_tvl,                                 &
         par_prev_fields(:,i_q_cl), par_next_fields(:,i_q_cl),                 &
         i_next, i_sat, sublevs,                                               &
         i_core_sat=i_core_sat )

end if  ! ( .NOT. i_call == i_call_det )


!---------------------------------------------------------------
! 5) Do dynamical processes (momentum equations)
!---------------------------------------------------------------

! Only needed in the main parcel mean and core ascents
if ( i_call == i_call_mean .or. i_call == i_call_core ) then

  ! Calculate updated parcel total-water mixing ratio
  call calc_q_tot( n_points, n_points_next,                                    &
                   par_next_fields(:,i_q_vap),                                 &
                   par_next_fields(:,i_qc_first:i_qc_last), q_tot )

  call momentum_eqn( n_points, n_fields_tot, l_res_source, l_down,             &
                     n_points_env, n_points_next,                              &
                     n_points_res, n_points_sublevs,                           &
                     l_within_bl, layer_mass_step, sum_massflux,               &
                     massflux_d, par_radius, q_tot, delta_t, Nsq_dry,          &
                     sublevs, i_next, j_wex,                                   &
                     env_k_fields(:,i_wind_u:i_wind_w),                        &
                     par_next_fields(:,i_wind_u:i_wind_w),                     &
                     par_w_drag,                                               &
                     res_source_fields = res_source_fields )

end if

!---------------------------------------------------------------
! 6) Diagnostics and output-checking
!---------------------------------------------------------------

! Save diagnostics from the plume-model
if ( l_diags ) then

  ! Time interval to traverse the model-level
  if ( plume_model_diags % delta_t % flag ) then
    i_diag = plume_model_diags % delta_t % i_super
    do ic = 1, n_points
      diags_super(ic,i_diag) = delta_t(ic)
    end do
  end if

  ! Dry-mass-flux in which the diags are valid
  if ( plume_model_diags % massflux_d_k % flag ) then
    i_diag = plume_model_diags % massflux_d_k % i_super
    do ic = 1, n_points
      diags_super(ic,i_diag) = massflux_d(ic)
    end do
  end if

end if


! Check outputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  ! Check parcel properties
  where_string = "End of parcel_dyn call for "                              // &
                 trim(adjustl(call_string))   // "; "                       // &
                 "par_next_fields"
  do i_field = 1, n_fields_tot
    call check_bad_values_cmpr( cmpr, k,                                       &
                                par_next_fields(:,i_field),                    &
                                where_string,                                  &
                                field_names(i_field),                          &
                                field_positive(i_field) )
  end do
  if ( present( res_source_fields ) ) then
    ! Check resolved-scale source-terms
    where_string = "End of parcel_dyn call for "                            // &
                   trim(adjustl(call_string))   // "; "                     // &
                   "res_source_fields"
    l_positive = .false.  ! Source terms may be positive or negative
    do i_field = 1, n_fields_tot
      call check_bad_values_cmpr( cmpr, k,                                     &
                                  res_source_fields(:,i_field),                &
                                  where_string,                                &
                                  field_names(i_field),                        &
                                  l_positive )
    end do
  end if
end if


return
end subroutine parcel_dyn


end module parcel_dyn_mod
