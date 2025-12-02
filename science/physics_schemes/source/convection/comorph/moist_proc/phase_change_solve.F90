! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module phase_change_solve_mod

implicit none


contains

! Subroutine takes the coefficients computed for the exhange of
! heat and moisture between hydrometeors and the parcel air, and
! computes an implicit solution for the parcel temperature and
! water-vapour mixing ratio at end-of-step, accounting for
! condensation/evaporation, deposition/sublimation, and melting
! of each condensed water species.
!
! Freezing and riming processes are not handled here; they are
! computed explicitly and their increments added on before we get
! here.  Note that the explicit riming processes assume water
! freezes on contact with ice.  If the ambient temperature and
! riming rate are such that freezing shouldn't have happened on
! contact, the implicit thermodynamic solve performed by this
! routine will remelt the wrongly frozen liquid.

subroutine phase_change_solve( n_points, n_points_super,                       &
             nc, index_ic, cmpr, k, call_string, linear_qs,                    &
             delta_z, delta_t, wind_w, prev_temp,                              &
             kq_cond, kt_cond, wf_cond, dq_frz_cond,                           &
             q_loc_cond, q_cond, cp_tot, temperature, q_vap,                   &
             l_diags, moist_proc_diags, n_points_diag, n_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, zero, one, n_cond_species,       &
                                 n_cond_species_liq, n_cond_species_ice,       &
                                 indi_thresh, name_length, cond_params,        &
                                 sqrt_min_delta,                               &
                                 i_check_imp_consistent,                       &
                                 i_check_bad_values_cmpr, i_check_bad_none
use linear_qs_mod, only: n_linear_qs_fields, i_ref_temp,                       &
                         i_qsat_liq_ref, i_qsat_ice_ref,                       &
                         i_dqsatdT_liq, i_dqsatdT_ice
use phase_change_coefs_mod, only: n_coefs
use moist_proc_diags_type_mod, only: moist_proc_diags_type
use cmpr_type_mod, only: cmpr_type

use lat_heat_mod, only: set_l_con, set_l_sub, set_l_fus,                       &
                        lat_heat_incr,                                         &
                        i_phase_change_con, i_phase_change_dep,                &
                        i_phase_change_mlt
use calc_phase_change_coefs_mod, only: calc_phase_change_coefs
use solve_tq_mod, only: solve_tq
use melt_ctl_mod, only: melt_ctl
use proc_incr_mod, only: proc_incr
use check_negatives_mod, only: check_negatives
use calc_cond_temp_mod, only: calc_cond_temp
use moist_proc_consistency_mod, only: moist_proc_consistency_1,                &
                                      moist_proc_consistency_2
use check_bad_values_mod, only: check_bad_values_cmpr

implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points in the condensed water species super-array
! (this might be reused on subsequent calls with bigger size
!  than needed, to save having to re-allocate it)
integer, intent(in) :: n_points_super

! Number of points where each condensed water species is non-zero
integer, intent(in) :: nc( n_cond_species )
! Indices of those points
integer, intent(in) :: index_ic( n_points, n_cond_species )

! Stuff used to make error messages more informative:
! Structure storing compression indices of the current points
type(cmpr_type), intent(in) :: cmpr
! Current model-level index
integer, intent(in) :: k
! Description of which call to moist_proc this is
character(len=name_length), intent(in) :: call_string

! Super-array containing qsat at a reference temperature,
! and dqsat/dT, for linearised qsat calculations
real(kind=real_cvprec), intent(in) :: linear_qs                                &
                          ( n_points, n_linear_qs_fields )

! Height interval for this step; = the vertical distance
! between the current point and the previous point where
! prev_temp is defined.
! If integrating downwards (as is conventional for
! Eulerian calculations), it should be negative.
real(kind=real_cvprec), intent(in) :: delta_z(n_points)

! Time interval for converting process rates to increments.
! For Eulerian calculations, this is the model timestep length,
! but for Lagrangian ascents, it is the time taken for the
! parcel to rise over the height interval delta_z,
! so delta_t = delta_z/wind_w
real(kind=real_cvprec), intent(in) :: delta_t(n_points)

! Vertical wind velocity
! For Euelerian calculations, this is the vertical wind-speed.
! For Lagrangian ascents, it is the vertical velocity of the
! parcel relative to the environment, such that
! wind_w = delta_z/delta_t
real(kind=real_cvprec), intent(in) :: wind_w(n_points)

! Parcel temperature at the previous model-level
real(kind=real_cvprec), intent(in) :: prev_temp(n_points)

! Vapour exchange coefficient for each hydrometeor species
real(kind=real_cvprec), intent(in out) :: kq_cond                              &
                          ( n_points, n_cond_species )
! Heat exchange coefficient for each hydrometeor species
real(kind=real_cvprec), intent(in out) :: kt_cond                              &
                          ( n_points, n_cond_species )
! These have intent inout because they occasionally need to be
! limited for numerical safety reasons.  Also we scale them
! by the timestep delta_t in this routine.

! Fall-speed of each hydrometeor species
real(kind=real_cvprec), intent(in) :: wf_cond                                  &
                          ( n_points, n_cond_species )

! Total freezing increment onto each ice hydrometeor species
! (includes homogeneous and heterogeneous freezing and riming)
! Needed for the hydrometeor surface heat budget, important for
! determining the melting rate
real(kind=real_cvprec), intent(in) :: dq_frz_cond                              &
             ( n_points, n_cond_species_liq+1 : n_cond_species )

! Local mixing ratio of each condensed water species,
! implicitly accounting for fall-out from current level / parcel
real(kind=real_cvprec), intent(in) :: q_loc_cond                               &
                          ( n_points, n_cond_species )

! Total available mixing ratio of each condensed water species
! (includes amount that falls through during this step, which
!  maybe considerably larger than the amount actually present
!  at a given instant).
! These are the values updated here; fall-out is calculated
! after this routine.
real(kind=real_cvprec), intent(in out) :: q_cond                               &
                          ( n_points_super, n_cond_species )

! Total heat capacity incremented by phase-changes
real(kind=real_cvprec), intent(in out) :: cp_tot(n_points)

! Parcel air temperature and water vapour mixing ratio
real(kind=real_cvprec), intent(in out) :: temperature(n_points)
real(kind=real_cvprec), intent(in out) :: q_vap(n_points)

! Master switch for whether or not to calculate any diagnostics
logical, intent(in) :: l_diags
! Structure containing diagnostic flags etc.
type(moist_proc_diags_type), intent(in) :: moist_proc_diags
! Number of points in the diagnostics super-array
integer, intent(in) :: n_points_diag
! Total number of diagnostics in the super-array
integer, intent(in) :: n_diags
! Super-array to store all the diagnostics
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                         ( n_points_diag, n_diags )


! Latent heats of condensation, sublimation and fusion
real(kind=real_cvprec) :: L_con(n_points)
real(kind=real_cvprec) :: L_sub(n_points)
real(kind=real_cvprec) :: L_fus(n_points)

! Estimated temperature and water-vapour mixing-ratio after
! implicit phase-changes (minus reference values)
real(kind=real_cvprec) :: imp_temp(n_points)
real(kind=real_cvprec) :: imp_q_vap(n_points)

! Super-arrays containing coefficients in the implicit formula
! for temperature and q_vap after phase-changes
real(kind=real_cvprec) :: coefs_temp( n_points, n_coefs )
real(kind=real_cvprec) :: coefs_q_vap( n_points, n_coefs )

! Super-array containing the 3 condensation / deposition rate
! coefficients for each condensed water species
real(kind=real_cvprec) :: coefs_cond                                           &
                          ( n_points, n_coefs, n_cond_species )

! Coefficients for melting rate of ice species
real(kind=real_cvprec) :: coefs_melt                                           &
   ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )

! Coefficients for deposition rate if ice is melting
real(kind=real_cvprec) :: coefs_cond_m                                         &
   ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )

! Condensation / deposition mixing-ratio increment for each
! species (has a negative value for evaporation / sublimation)
real(kind=real_cvprec) :: dq_cond( n_points, n_cond_species )

! Mixing ratio increment due to melting, for each ice species
real(kind=real_cvprec) :: dq_melt                                              &
            ( n_points, n_cond_species_liq+1 : n_cond_species )

! Flag for whether each ice species is melting
logical :: l_melt                                                              &
            ( n_points, n_cond_species_liq+1 : n_cond_species )

! Flag to pass into solve_tq
logical, parameter :: l_full_do_true = .true.

! Flags for whether to do full-field do-loops instead of
! indirect indexing for each condensate species
logical :: l_full_do(n_cond_species)

! Flag input to calc_cond_temp
logical :: l_ice
! Address of dqsat/dT w.r.t. current water phase in linear_qs
integer :: i_dqsatdt

! Max allowable dimensionless condensation coefficient, for
! numerical stability / avoiding loss of precision
real(kind=real_cvprec), parameter :: max_cond_coef = one / sqrt_min_delta

! Name of a field (for error message)
character(len=name_length) :: field_name
! Flag for whether field is positive-only
logical :: l_positive
! Description of where we are in the code, for error messages
character(len=name_length) :: where_string

! Loop counters
integer :: ic, ic2, i_cond, i_liq, i_ice, i_super, i_coef


! Check inputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  where_string = "Start of phase_change_solve call for "      //               &
                 trim(adjustl(call_string))
  l_positive = .true.
  field_name = "temperature"
  call check_bad_values_cmpr( cmpr, k, temperature,                            &
         where_string, field_name, l_positive )
  field_name = "q_vap"
  call check_bad_values_cmpr( cmpr, k, q_vap,                                  &
         where_string, field_name, l_positive )
  do i_cond = 1, n_cond_species
    field_name = "q_" //                                                       &
                 trim(adjustl( cond_params(i_cond)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, q_cond(:,i_cond),                     &
           where_string, field_name, l_positive )
  end do
  do i_cond = 1, n_cond_species
    field_name = "kq_" //                                                      &
                 trim(adjustl( cond_params(i_cond)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, kq_cond(:,i_cond),                    &
           where_string, field_name, l_positive )
  end do
  do i_cond = 1, n_cond_species
    field_name = "kt_" //                                                      &
                 trim(adjustl( cond_params(i_cond)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, kt_cond(:,i_cond),                    &
           where_string, field_name, l_positive )
  end do
  do i_cond = 1, n_cond_species
    field_name = "wf_" //                                                      &
                 trim(adjustl( cond_params(i_cond)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, wf_cond(:,i_cond),                    &
           where_string, field_name, l_positive )
  end do
  do i_cond = 1, n_cond_species
    field_name = "q_loc_" //                                                   &
                 trim(adjustl( cond_params(i_cond)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, q_loc_cond(:,i_cond),                 &
           where_string, field_name, l_positive )
  end do
  do i_ice = n_cond_species_liq+1, n_cond_species
    field_name = "dq_frz_" //                                                  &
                 trim(adjustl( cond_params(i_ice)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, dq_frz_cond(:,i_ice),                 &
           where_string, field_name, l_positive )
  end do
end if


!----------------------------------------------------------------
! 1) Initialisations
!----------------------------------------------------------------

! Set flags for whether indirect indexing is required for each
! condensate species
do i_cond = 1, n_cond_species
  ! Use full-field calculations instead of indirect indexing whenever
  ! fraction of points is above threshold,
  ! i.e. the calculations for this species will not be sparse enough
  ! to justify indirect indexing.
  l_full_do(i_cond) =   real(nc(i_cond),real_cvprec)                           &
                      > indi_thresh * real(n_points,real_cvprec)
end do

! Initialise process rate increments to zero and melting flags to false.
do i_cond = 1, n_cond_species
  do ic = 1, n_points
    dq_cond(ic,i_cond) = zero
  end do
end do
do i_ice = n_cond_species_liq+1, n_cond_species
  do ic = 1, n_points
    dq_melt(ic,i_ice) = zero
    l_melt(ic,i_ice) = .false.
  end do
end do

! Compute latent heats
call set_l_con( n_points, temperature, L_con )
call set_l_sub( n_points, temperature, L_sub )
call set_l_fus( n_points, temperature, L_fus )


do i_cond = 1, n_cond_species
  if ( nc(i_cond) > 0 ) then
    if ( l_full_do(i_cond) ) then

      ! Scale the vapour and heat exchange coefficients by the
      ! time interval.  Everything else now calculated will be
      ! in terms of increments rather than tendencies
      do ic = 1, n_points
        kq_cond(ic,i_cond) = kq_cond(ic,i_cond) * delta_t(ic)
        kt_cond(ic,i_cond) = kt_cond(ic,i_cond) * delta_t(ic)
      end do

      ! Numerical safety-check; don't let kq * delta_t
      ! (which is the dimensionless coefficient for relaxation
      ! towards saturation) exceed 1/sqrt(epsilon);
      ! beyond this point, the process-rate formulae may have
      ! unacceptably poor precision due to cancellation of large
      ! opposing terms from
      ! q - qsat(T_ref) and dqdsat/dT (T - T_ref).
      do ic = 1, n_points
        kq_cond(ic,i_cond) = min( kq_cond(ic,i_cond),                          &
                                  max_cond_coef )
      end do

    else
      ! Compressed version of the above when only needed at
      ! a small fraction of points

      do ic2 = 1, nc(i_cond)
        ic = index_ic(ic2,i_cond)
        kq_cond(ic,i_cond) = kq_cond(ic,i_cond) * delta_t(ic)
        kt_cond(ic,i_cond) = kt_cond(ic,i_cond) * delta_t(ic)
        kq_cond(ic,i_cond) = min( kq_cond(ic,i_cond),                          &
                                  max_cond_coef )
      end do

    end if
  end if
end do


!----------------------------------------------------------------
! 2) Calculate condensation/evaporation and melting coefficients
!    for each condensed water species
!----------------------------------------------------------------

! This routine also adds the condensation/evaporation
! contributions for each hydrometeor species onto the main
! coefficients for the implicit solve of T,q, assuming no
! melting occurs.  Melting will be tested for later...
! NOTE: the coefficients calculated here are all scaled
! by the timestep delta_t, so that wherever the coefficients
! are used, the result is an increment, not a process rate.
call calc_phase_change_coefs( n_points, nc, index_ic, l_full_do,               &
                              linear_qs(:,i_ref_temp),                         &
                              linear_qs(:,i_qsat_liq_ref),                     &
                              linear_qs(:,i_qsat_ice_ref),                     &
                              linear_qs(:,i_dqsatdT_liq),                      &
                              linear_qs(:,i_dqsatdT_ice),                      &
                              L_con, L_sub, L_fus, cp_tot,                     &
                              delta_t, delta_z, wind_w,                        &
                              prev_temp,                                       &
                              kq_cond, kt_cond, wf_cond,                       &
                              q_loc_cond, dq_frz_cond,                         &
                              q_vap, coefs_cond,                               &
                              coefs_melt, coefs_cond_m,                        &
                              coefs_temp, coefs_q_vap )


! Check outputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  where_string = "phase_change_solve, call for "             //                &
                 trim(adjustl(call_string))                  //                &
                 "; after calc_phase_change_coefs"
  l_positive = .false.
  ! Check condensation and melting coefficients
  do i_cond = 1, n_cond_species
    do i_coef = 1, n_coefs
      write(field_name,"(I1)") i_coef
      field_name = "coefs_cond_" //                                            &
               trim(adjustl( cond_params(i_cond)%pt % cond_name )) //          &
               "_c" // trim(adjustl( field_name ))
      call check_bad_values_cmpr( cmpr, k,                                     &
             coefs_cond(:,i_coef,i_cond),                                      &
             where_string, field_name, l_positive )
    end do
  end do
  do i_ice = n_cond_species_liq+1, n_cond_species
    do i_coef = 1, n_coefs
      write(field_name,"(I1)") i_coef
      field_name = "coefs_melt_" //                                            &
               trim(adjustl( cond_params(i_ice)%pt % cond_name ))  //          &
               "_c" // trim(adjustl( field_name ))
      call check_bad_values_cmpr( cmpr, k,                                     &
             coefs_melt(:,i_coef,i_ice),                                       &
             where_string, field_name, l_positive )
    end do
  end do

end if


!----------------------------------------------------------------
! 3) Compute initial-guess implicit solution for T and q
!----------------------------------------------------------------

! Use the coefficients to estimate T and q after phase-changes,
! based on a linear system.
! At this point, solving the system assuming there is no melting,
! and negative hydrometeor mixing ratios may occur.
call solve_tq( n_points, l_full_do_true,                                       &
               linear_qs(:,i_ref_temp),                                        &
               linear_qs(:,i_qsat_liq_ref),                                    &
               coefs_temp, coefs_q_vap,                                        &
               temperature, q_vap, imp_temp, imp_q_vap,                        &
               kq_cond, kt_cond,                                               &
               coefs_cond, coefs_cond_m, coefs_melt )


!----------------------------------------------------------------
! 4) Check whether ice species need to melt, and modify the
!    implicit solve to include melting coefficients where
!    appropriate
!----------------------------------------------------------------

if ( n_cond_species_ice > 0 ) then
  ! If any points contain ice
  if ( maxval( nc(n_cond_species_liq+1:n_cond_species) )                       &
       > 0 ) then

    ! This routine also calculates the melting increments dq_melt
    call melt_ctl( n_points, nc, index_ic,                                     &
                   linear_qs, cp_tot, L_sub, L_fus,                            &
                   temperature, q_vap, kq_cond, kt_cond,                       &
                   coefs_melt, coefs_cond, coefs_cond_m,                       &
                   dq_melt, l_melt, coefs_temp, coefs_q_vap,                   &
                   imp_temp, imp_q_vap )

  end if
end if


! Check outputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  where_string = "phase_change_solve, call for "             //                &
                 trim(adjustl(call_string))                  //                &
                 "; after melt_ctl"
  l_positive = .false.
  ! Check condensation and melting coefficients
  do i_cond = 1, n_cond_species
    do i_coef = 1, n_coefs
      write(field_name,"(I1)") i_coef
      field_name = "coefs_cond_" //                                            &
               trim(adjustl( cond_params(i_cond)%pt % cond_name )) //          &
               "_c" // trim(adjustl( field_name ))
      call check_bad_values_cmpr( cmpr, k,                                     &
             coefs_cond(:,i_coef,i_cond),                                      &
             where_string, field_name, l_positive )
    end do
  end do
  do i_ice = n_cond_species_liq+1, n_cond_species
    do i_coef = 1, n_coefs
      write(field_name,"(I1)") i_coef
      field_name = "coefs_melt_" //                                            &
               trim(adjustl( cond_params(i_ice)%pt % cond_name ))  //          &
               "_c" // trim(adjustl( field_name ))
      call check_bad_values_cmpr( cmpr, k,                                     &
             coefs_melt(:,i_coef,i_ice),                                       &
             where_string, field_name, l_positive )
    end do
  end do

end if


!----------------------------------------------------------------
! 5) Calculate condensation / evaporation rates for all species
!    consistent with the current implicit solve
!----------------------------------------------------------------
do i_cond = 1, n_cond_species
  if ( nc(i_cond) > 0 ) then
    call proc_incr( n_points, nc(i_cond), index_ic(:,i_cond),                  &
                    l_full_do(i_cond),                                         &
                    imp_temp, imp_q_vap,                                       &
                    coefs_cond(:,:,i_cond), dq_cond(:,i_cond) )
  end if
end do


!----------------------------------------------------------------
! 6) Calculate total tendency for each hydrometeor species and
!    check to avoid negative hydrometeor mixing ratios
!----------------------------------------------------------------

! The solution for T, q is only implicit w.r.t. T and q;
! it is explicit w.r.t. the hydrometeor mixing ratios.  As such,
! there is nothing to stop the evaporation and melting rates
! from "overshooting", and consuming more of a hydrometeor
! species' mass than is actually present.
! In this event, the affected processes don't need to be solved
! implicitly, as the known solution is just to evaporate / melt
! exactly all of the mixing ratio remaining.
! However, to retain the correct implicit solution for any
! other condensed water species which haven't "run out", we
! need to ensure the implicit coefficients for T,q account
! for the truncation of "overshot" evaporation etc.
! This routine modifies the implicit coefficients and repeats
! the solve if needed to prevent it from producing negative
! mixing-ratios, while maintaining the consistency of the
! implicit solution...

call check_negatives( n_points, n_points_super, nc, index_ic,                  &
                      linear_qs, cp_tot, L_con, L_sub, L_fus,                  &
                      temperature, q_vap, q_cond,                              &
                      kq_cond, kt_cond,                                        &
                      coefs_melt, coefs_cond, coefs_cond_m,                    &
                      dq_cond, dq_melt, l_melt,                                &
                      coefs_temp, coefs_q_vap,                                 &
                      imp_temp, imp_q_vap )


if ( i_check_imp_consistent > i_check_bad_none ) then
  ! Check self-consistency of the final implicit solution
  call moist_proc_consistency_1( cmpr, k, call_string,                         &
                                 nc, index_ic, l_full_do,                      &
                                 linear_qs, imp_temp, imp_q_vap,               &
                                 coefs_cond, kq_cond, dq_cond )
end if


!----------------------------------------------------------------
! 7) Add on increments for all implicit phase-changes
!----------------------------------------------------------------

! For each present liquid species
do i_liq = 1, n_cond_species_liq
  if ( nc(i_liq) > 0 ) then
    if ( l_full_do(i_liq) ) then
      ! Full-field condensation / evaporation
      do ic = 1, n_points
        q_cond(ic,i_liq) = q_cond(ic,i_liq) + dq_cond(ic,i_liq)
        q_vap(ic)        = q_vap(ic)        - dq_cond(ic,i_liq)
      end do
      call lat_heat_incr( n_points, nc(i_liq),                                 &
                          i_phase_change_con,                                  &
                          cp_tot, temperature,                                 &
                          dq=dq_cond(:,i_liq) )
    else
      ! Indirect indexing condensation / evaporation
      do ic2 = 1, nc(i_liq)
        ic = index_ic(ic2,i_liq)
        q_cond(ic,i_liq) = q_cond(ic,i_liq) + dq_cond(ic,i_liq)
        q_vap(ic)        = q_vap(ic)        - dq_cond(ic,i_liq)
      end do
      call lat_heat_incr( n_points, nc(i_liq),                                 &
                          i_phase_change_con,                                  &
                          cp_tot, temperature,                                 &
                          index_ic=index_ic(:,i_liq),                          &
                          dq=dq_cond(:,i_liq) )
    end if
  end if
end do

! For each present ice species
if ( n_cond_species_ice > 0 ) then
  do i_ice = n_cond_species_liq+1, n_cond_species
    if ( nc(i_ice) > 0 ) then
      if ( l_full_do(i_ice) ) then
        ! Full-field deposition / sublimation
        do ic = 1, n_points
          q_cond(ic,i_ice) = q_cond(ic,i_ice) + dq_cond(ic,i_ice)
          q_vap(ic)        = q_vap(ic)        - dq_cond(ic,i_ice)
        end do
        call lat_heat_incr( n_points, nc(i_ice),                               &
                            i_phase_change_dep,                                &
                            cp_tot, temperature,                               &
                            dq=dq_cond(:,i_ice) )
        ! Full-field melting
        i_liq = cond_params(i_ice)%pt % i_cond_frzmlt
        do ic = 1, n_points
          q_cond(ic,i_ice) = q_cond(ic,i_ice) - dq_melt(ic,i_ice)
          q_cond(ic,i_liq) = q_cond(ic,i_liq) + dq_melt(ic,i_ice)
        end do
        call lat_heat_incr( n_points, nc(i_ice),                               &
                            i_phase_change_mlt,                                &
                            cp_tot, temperature,                               &
                            dq=dq_melt(:,i_ice) )
      else
        ! Indirect indexing deposition / sublimation
        do ic2 = 1, nc(i_ice)
          ic = index_ic(ic2,i_ice)
          q_cond(ic,i_ice) = q_cond(ic,i_ice) + dq_cond(ic,i_ice)
          q_vap(ic)        = q_vap(ic)        - dq_cond(ic,i_ice)
        end do
        call lat_heat_incr( n_points, nc(i_ice),                               &
                            i_phase_change_dep,                                &
                            cp_tot, temperature,                               &
                            index_ic=index_ic(:,i_ice),                        &
                            dq=dq_cond(:,i_ice) )
        ! Indirect indexing melting
        i_liq = cond_params(i_ice)%pt % i_cond_frzmlt
        do ic2 = 1, nc(i_ice)
          ic = index_ic(ic2,i_ice)
          q_cond(ic,i_ice) = q_cond(ic,i_ice) - dq_melt(ic,i_ice)
          q_cond(ic,i_liq) = q_cond(ic,i_liq) + dq_melt(ic,i_ice)
        end do
        call lat_heat_incr( n_points, nc(i_ice),                               &
                            i_phase_change_mlt,                                &
                            cp_tot, temperature,                               &
                            index_ic=index_ic(:,i_ice),                        &
                            dq=dq_melt(:,i_ice) )
      end if
    end if
  end do
end if

! Check to avoid negative condensed-water mixing-ratios;
! Very slight negative values can occasionally occur due to
! rounding errors in the process-rates.  Just reset to
! zero when this occurs
do i_cond = 1, n_cond_species
  if ( nc(i_cond) > 0 ) then
    do ic2 = 1, nc(i_cond)
      ic = index_ic(ic2,i_cond)
      q_cond(ic,i_cond) = max( q_cond(ic,i_cond), zero )
    end do
  end if
end do
! Note: if non-negligible negatives are removed by this check,
! then conservation of total-water will be broken.  Use the
! conservation-checking option to catch this potential problem.


if ( i_check_imp_consistent > i_check_bad_none ) then
  ! Check whether final T,q are consistent with their implicit
  ! solution values
  call moist_proc_consistency_2( cmpr, k, call_string,                         &
                                 imp_temp, imp_q_vap,                          &
                                 temperature, q_vap,                           &
                                 linear_qs(:,i_ref_temp),                      &
                                 linear_qs(:,i_qsat_liq_ref),                  &
                                 dq_cond, dq_melt, cp_tot )
end if


!----------------------------------------------------------------
! 8) Calculate diagnostics
!----------------------------------------------------------------

if ( l_diags ) then

  ! Condensation / evaporation increment
  do i_cond = 1, n_cond_species
    if ( nc(i_cond) > 0 ) then
      ! If condensation / evaporation increment diag requested
      if ( moist_proc_diags % diags_cond(i_cond)%pt                            &
           % dq_cond % flag ) then
        ! Extract super-array address
        i_super = moist_proc_diags % diags_cond(i_cond)%pt                     &
                  % dq_cond % i_super
        ! Copy diagnostic...
        if ( l_full_do(i_cond) ) then
          do ic = 1, n_points
            diags_super(ic,i_super) = dq_cond(ic,i_cond)
          end do
        else
          do ic2 = 1, nc(i_cond)
            ic = index_ic(ic2,i_cond)
            diags_super(ic,i_super) = dq_cond(ic,i_cond)
          end do
        end if
      end if
    end if
  end do

  ! Melting increment
  do i_ice = n_cond_species_liq+1, n_cond_species
    if ( nc(i_ice) > 0 ) then

      ! If melting increment diag requested for species i_ice
      if ( moist_proc_diags % diags_cond(i_ice)%pt                             &
           % dq_frzmlt % flag ) then
        ! Extract super-array address
        i_super = moist_proc_diags % diags_cond(i_ice)%pt                      &
                  % dq_frzmlt % i_super
        ! Update diagnostic...
        if ( l_full_do(i_ice) ) then
          do ic = 1, n_points
            diags_super(ic,i_super) = diags_super(ic,i_super)                  &
                                    - dq_melt(ic,i_ice)
          end do
        else
          do ic2 = 1, nc(i_ice)
            ic = index_ic(ic2,i_ice)
            diags_super(ic,i_super) = diags_super(ic,i_super)                  &
                                    - dq_melt(ic,i_ice)
          end do
        end if
      end if

      ! Find liquid species the ice converts to when it melts
      i_liq = cond_params(i_ice)%pt % i_cond_frzmlt

      ! If melting increment diag requested for species i_liq
      if ( moist_proc_diags % diags_cond(i_liq)%pt                             &
           % dq_frzmlt % flag ) then
        ! Extract super-array address
        i_super = moist_proc_diags % diags_cond(i_liq)%pt                      &
                  % dq_frzmlt % i_super
        ! Update diagnostic...
        if ( l_full_do(i_ice) ) then
          do ic = 1, n_points
            diags_super(ic,i_super) = diags_super(ic,i_super)                  &
                                    + dq_melt(ic,i_ice)
          end do
        else
          do ic2 = 1, nc(i_ice)
            ic = index_ic(ic2,i_ice)
            diags_super(ic,i_super) = diags_super(ic,i_super)                  &
                                    + dq_melt(ic,i_ice)
          end do
        end if
      end if

    end if
  end do

  ! Hydrometeor temperatures...
  ! For all present condensed water species
  do i_cond = 1, n_cond_species
    if ( nc(i_cond) > 0 ) then
      ! If temperature of current liquid species is requested
      if ( moist_proc_diags % diags_cond(i_cond)%pt                            &
           % cond_temp % flag ) then
        ! Extract super-array address
        i_super = moist_proc_diags % diags_cond(i_cond)%pt                     &
                  % cond_temp % i_super
        ! Differences for liquid / ice species calculations
        if ( i_cond > n_cond_species_liq ) then
          l_ice = .true.
          i_dqsatdt = i_dqsatdt_ice
        else
          l_ice = .false.
          i_dqsatdt = i_dqsatdt_liq
        end if
        ! Calculate hydrometeor temperature
        call calc_cond_temp( n_points, nc(i_cond),                             &
                             index_ic(:,i_cond), l_full_do(i_cond), l_ice,     &
                             linear_qs(:,i_qsat_liq_ref),                      &
                             linear_qs(:,i_qsat_ice_ref),                      &
                             linear_qs(:,i_dqsatdt),                           &
                             kq_cond(:,i_cond),                                &
                             imp_temp, imp_q_vap,                              &
                             coefs_cond(:,:,i_cond),                           &
                             diags_super(:,i_super) )
        ! calc_cond_temp outputs T_cond - T_ref;
        ! add on reference temperature to get actual temperature
        do ic2 = 1, nc(i_cond)
          ic = index_ic(ic2,i_cond)
          diags_super(ic,i_super) = diags_super(ic,i_super)                    &
                                  + linear_qs(ic,i_ref_temp)
        end do
      end if
    end if
  end do

end if  ! ( l_diags )


! Check outputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  where_string = "End of phase_change_solve call for "        //               &
                 trim(adjustl(call_string))
  l_positive = .true.

  ! Check temperature, vapour and condensed water mixing-ratios
  field_name = "temperature"
  call check_bad_values_cmpr( cmpr, k, temperature,                            &
         where_string, field_name, l_positive )
  field_name = "q_vap"
  call check_bad_values_cmpr( cmpr, k, q_vap,                                  &
         where_string, field_name, l_positive )
  do i_cond = 1, n_cond_species
    field_name = "q_" //                                                       &
                 trim(adjustl( cond_params(i_cond)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, q_cond(:,i_cond),                     &
           where_string, field_name, l_positive )
  end do

  ! Check condensation and melting increments
  do i_ice = n_cond_species_liq+1, n_cond_species
    field_name = "dq_melt_" //                                                 &
                 trim(adjustl( cond_params(i_ice)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, dq_melt(:,i_ice),                     &
           where_string, field_name, l_positive )
  end do
  l_positive = .false.
  do i_cond = 1, n_cond_species
    field_name = "dq_cond_" //                                                 &
                 trim(adjustl( cond_params(i_cond)%pt % cond_name ))
    call check_bad_values_cmpr( cmpr, k, dq_cond(:,i_cond),                    &
           where_string, field_name, l_positive )
  end do

  ! Check condensation and melting coefficeints
  do i_cond = 1, n_cond_species
    do i_coef = 1, n_coefs
      write(field_name,"(I1)") i_coef
      field_name = "coefs_cond_" //                                            &
               trim(adjustl( cond_params(i_cond)%pt % cond_name )) //          &
               "_c" // trim(adjustl( field_name ))
      call check_bad_values_cmpr( cmpr, k,                                     &
             coefs_cond(:,i_coef,i_cond),                                      &
             where_string, field_name, l_positive )
    end do
  end do
  do i_ice = n_cond_species_liq+1, n_cond_species
    do i_coef = 1, n_coefs
      write(field_name,"(I1)") i_coef
      field_name = "coefs_melt_" //                                            &
               trim(adjustl( cond_params(i_ice)%pt % cond_name ))  //          &
               "_c" // trim(adjustl( field_name ))
      call check_bad_values_cmpr( cmpr, k,                                     &
             coefs_melt(:,i_coef,i_ice),                                       &
             where_string, field_name, l_positive )
    end do
  end do

end if


return
end subroutine phase_change_solve

end module phase_change_solve_mod
