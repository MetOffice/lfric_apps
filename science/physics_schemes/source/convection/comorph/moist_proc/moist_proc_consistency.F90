! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module moist_proc_consistency_mod

implicit none

! Max length of string containing coordinates and values of
! any inconsistent values found
integer, parameter :: str_len = 200

contains


!----------------------------------------------------------------
! Subroutine to check that the calculated hydrometeor
! temperatures and condensation increments are consistent with
! eachother and with the implicit solution for T and q
!----------------------------------------------------------------
subroutine moist_proc_consistency_1( cmpr, k, call_string,                     &
                                     nc, index_ic, l_full_do,                  &
                                     linear_qs,                                &
                                     imp_temp, imp_q_vap,                      &
                                     coefs_cond, kq_cond,                      &
                                     dq_cond )

use comorph_constants_mod, only: real_cvprec, n_cond_species,                  &
                                 n_cond_species_liq,                           &
                                 newline, cond_params, i_check_imp_consistent, &
                                 name_length, min_delta, min_float
use linear_qs_mod, only: n_linear_qs_fields,                                   &
                         i_qsat_liq_ref, i_qsat_ice_ref,                       &
                         i_dqsatdT_liq, i_dqsatdT_ice
use phase_change_coefs_mod, only: n_coefs, i_q, i_t, i_0
use cmpr_type_mod, only: cmpr_type

use calc_cond_temp_mod, only: calc_cond_temp

use check_bad_values_mod, only: check_bad_values_print

implicit none


! Stuff used to make error messages more informative:
! Structure containing compression list indices
type(cmpr_type), intent(in) :: cmpr
! Current model-level
integer, intent(in) :: k
! Character string describing which call to moist_proc this is
character(len=name_length), intent(in) :: call_string

! Number of points where each condensed water species is non-zero
integer, intent(in) :: nc( n_cond_species )
! Indices of those points
integer, intent(in) :: index_ic( cmpr%n_points, n_cond_species )
! Flags for whether to do full-field do-loops instead of
! indirect indexing for each condensate species
logical :: l_full_do(n_cond_species)

! Super-array containing qsat at a reference temperature,
! and dqsat/dT, for linearised qsat calculations
real(kind=real_cvprec), intent(in) :: linear_qs                                &
                               ( cmpr%n_points, n_linear_qs_fields )

! Implicit solution for temperature and water-vapour mixing-ratio
! (minus reference values)
real(kind=real_cvprec), intent(in) :: imp_temp(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: imp_q_vap(cmpr%n_points)

! Condensation / evaporation rate coefficients
real(kind=real_cvprec), intent(in) :: coefs_cond                               &
                          ( cmpr%n_points, n_coefs, n_cond_species )

! Vapour exchange coefficient for each species
real(kind=real_cvprec), intent(in) :: kq_cond                                  &
                                   ( cmpr%n_points, n_cond_species )

! Condensation / evaporation increment consistent with the
! current implicit solution
real(kind=real_cvprec), intent(in) :: dq_cond                                  &
                                   ( cmpr%n_points, n_cond_species )

! Temperature of hydrometeor
real(kind=real_cvprec) :: cond_temp(cmpr%n_points)

! Estimate of condensation / evaporation increment based on
! hydrometeor temperature, minus actual dq_cond
real(kind=real_cvprec) :: dq_error(cmpr%n_points)

! Maximum acceptable error, given the expected loss of precision
! in the calculation
real(kind=real_cvprec) :: max_error(cmpr%n_points)

! Array size for terms below
integer, parameter :: n_terms = 7

! List of terms used to estimate max_error
real(kind=real_cvprec) :: terms(n_terms)

! Flag to pass into calc_cond_temp
logical :: l_ice

! Address of dqsat w.r.t. current water phase in linear_qs
integer :: i_qsat
integer :: i_dqsatdt

! Flag for points where inconsistency found
logical :: l_bad(cmpr%n_points)

! Scaling factor for numerical tolerance in error of dq_cond
real(kind=real_cvprec), parameter :: tolerance = 10.0_real_cvprec

! Loop counters
integer :: ic, ic2, i_cond


! Loop over all present species
do i_cond = 1, n_cond_species
  if ( nc(i_cond) > 0 ) then

    ! Differences for liquid / ice species calculations
    if ( i_cond > n_cond_species_liq ) then
      l_ice = .true.
      i_qsat = i_qsat_ice_ref
      i_dqsatdt = i_dqsatdt_ice
    else
      l_ice = .false.
      i_qsat = i_qsat_liq_ref
      i_dqsatdt = i_dqsatdt_liq
    end if

    ! Calculate hydrometeor temperature:
    call calc_cond_temp( cmpr%n_points, nc(i_cond),                            &
                         index_ic(:,i_cond), l_full_do(i_cond), l_ice,         &
                         linear_qs(:,i_qsat_liq_ref),                          &
                         linear_qs(:,i_qsat_ice_ref),                          &
                         linear_qs(:,i_dqsatdt),                               &
                         kq_cond(:,i_cond), imp_temp, imp_q_vap,               &
                         coefs_cond(:,:,i_cond), cond_temp )

    if ( l_full_do(i_cond) ) then

      ! Calculate condensation / evaporation increment based on
      ! the diagnosed hydrometeor temperature
      ! (note cond_temp stores T_cond - T_ref)
      do ic = 1, cmpr%n_points
        dq_error(ic) = kq_cond(ic,i_cond)                                      &
          * ( imp_q_vap(ic) + ( linear_qs(ic,i_qsat_liq_ref)                   &
                              - linear_qs(ic,i_qsat) )                         &
            - linear_qs(ic,i_dqsatdt) * cond_temp(ic)                          &
            )        - dq_cond(ic,i_cond)
      end do

      ! Work out maximum expected error; loss of precision should scale
      ! with epsilon times the size of the largest of the terms that
      ! have been added or subtracted to calculate dq_error above:
      do ic = 1, cmpr%n_points
        ! Terms in the above formula
        terms(1) = kq_cond(ic,i_cond) * imp_q_vap(ic)
        terms(2) = kq_cond(ic,i_cond) * ( linear_qs(ic,i_qsat_liq_ref)         &
                                        - linear_qs(ic,i_qsat) )
        terms(3) = kq_cond(ic,i_cond) * linear_qs(ic,i_dqsatdt) * cond_temp(ic)
        terms(4) = dq_cond(ic,i_cond)
        ! Terms used in calc_cond_temp to compute cond_temp
        terms(5) = coefs_cond(ic,i_q,i_cond) * imp_q_vap(ic)
        terms(6) = coefs_cond(ic,i_t,i_cond) * imp_temp(ic)
        terms(7) = coefs_cond(ic,i_0,i_cond)
        max_error(ic) = maxval(abs(terms)) * tolerance * min_delta
        ! Some compilers allow numbers less than tiny to be stored
        ! with reduced precision, which means the above threshold based
        ! on epsilon goes wrong.  Avoid this problem by not letting
        ! max_error fall below tolerance * tiny
        max_error(ic) = max( max_error(ic), tolerance * min_float )
      end do

      ! Flag points where inconsistency found
      do ic = 1, cmpr%n_points
        l_bad(ic) = abs(dq_error(ic)) > max_error(ic)
      end do

    else
      ! Indirect indexing version of the same calculation

      do ic2 = 1, nc(i_cond)
        ic = index_ic(ic2,i_cond)
        dq_error(ic) = kq_cond(ic,i_cond)                                      &
          * ( imp_q_vap(ic) + ( linear_qs(ic,i_qsat_liq_ref)                   &
                              - linear_qs(ic,i_qsat) )                         &
            - linear_qs(ic,i_dqsatdt) * cond_temp(ic)                          &
            )        - dq_cond(ic,i_cond)
        terms(1) = kq_cond(ic,i_cond) * imp_q_vap(ic)
        terms(2) = kq_cond(ic,i_cond) * ( linear_qs(ic,i_qsat_liq_ref)         &
                                        - linear_qs(ic,i_qsat) )
        terms(3) = kq_cond(ic,i_cond) * linear_qs(ic,i_dqsatdt) * cond_temp(ic)
        terms(4) = dq_cond(ic,i_cond)
        terms(5) = coefs_cond(ic,i_q,i_cond) * imp_q_vap(ic)
        terms(6) = coefs_cond(ic,i_t,i_cond) * imp_temp(ic)
        terms(7) = coefs_cond(ic,i_0,i_cond)
        max_error(ic) = maxval(abs(terms)) * tolerance * min_delta
        max_error(ic) = max( max_error(ic), tolerance * min_float )
      end do

      do ic = 1, cmpr%n_points
        l_bad(ic) = .false.
      end do
      do ic2 = 1, nc(i_cond)
        ic = index_ic(ic2,i_cond)
        l_bad(ic) = abs(dq_error(ic)) > max_error(ic)
      end do

    end if

    ! Call routine to print coordinates and values at "bad" points
    call check_bad_values_print(                                               &
           cmpr, k, l_bad, i_check_imp_consistent,                             &
           "(" // trim(adjustl(call_string)) // ")"            //newline//     &
           "Hydrometeor temperature and condensation / "                //     &
           "evaporation rate are inconsistent.  "              //newline//     &
           "Condensed water species: "                                  //     &
           trim(adjustl(cond_params(i_cond)%pt%cond_name))     //newline//     &
           "dq_cond, error, tolerance: ",                                      &
           field1=dq_cond(:,i_cond), field2=dq_error, field3=max_error )

  end if  ! ( nc(i_cond) > 0 )
end do  ! i_cond = 1, n_cond_species


return
end subroutine moist_proc_consistency_1


!----------------------------------------------------------------
! Subroutine to check that the final incremented values of T,q
! are consistent with the implicit solution.
!----------------------------------------------------------------
subroutine moist_proc_consistency_2( cmpr, k, call_string,                     &
                                     imp_temp, imp_q_vap,                      &
                                     temperature, q_vap,                       &
                                     ref_temp, qsat_liq_ref,                   &
                                     dq_cond, dq_melt, cp_tot )

use comorph_constants_mod, only: real_cvprec, zero, newline,                   &
                                 i_check_imp_consistent,                       &
                                 name_length, min_delta,                       &
                                 n_cond_species, n_cond_species_liq
use lat_heat_mod, only: lat_heat_incr, i_phase_change_evp, i_phase_change_sub, &
                        i_phase_change_frz
use cmpr_type_mod, only: cmpr_type

use check_bad_values_mod, only: check_bad_values_print

implicit none

! Stuff used to make error messages more informative:
! Structure containing compression list indices
type(cmpr_type), intent(in) :: cmpr
! Current model-level
integer, intent(in) :: k
! Character string describing which call to moist_proc this is
character(len=name_length), intent(in) :: call_string

! Implicit solution for temperature and water-vapour mixing-ratio
! (minus reference values)
real(kind=real_cvprec), intent(in) :: imp_temp(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: imp_q_vap(cmpr%n_points)

! Current updated temperature and water-vapour mixing-ratio
real(kind=real_cvprec), intent(in) :: temperature(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: q_vap(cmpr%n_points)

! Reference temperature and qsat value
real(kind=real_cvprec), intent(in) :: ref_temp(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: qsat_liq_ref(cmpr%n_points)

! Increments to each hydrometeor species due to:
! condensation / evaporation
real(kind=real_cvprec), intent(in) :: dq_cond                                  &
                                      ( cmpr%n_points, n_cond_species )
! melting
real(kind=real_cvprec), intent(in) :: dq_melt                                  &
               ( cmpr%n_points, n_cond_species_liq+1 : n_cond_species )

! Total heat capacity per unit dry-mass (used as a work array)
real(kind=real_cvprec), intent(in out) :: cp_tot(cmpr%n_points)

! Temperature and q_vap increment scales used for finding max error
real(kind=real_cvprec) :: delta(cmpr%n_points)
! Work array used to find the above
real(kind=real_cvprec) :: work(cmpr%n_points)

! Maximum acceptable error
real(kind=real_cvprec) :: max_error(cmpr%n_points)

! Actual error found in final T,q compared to implicit solution
real(kind=real_cvprec) :: final_error(cmpr%n_points)

! Logical flags set to true where inconsistency found
logical :: l_bad(cmpr%n_points)

! Tolerances for error of implicit solution
real(kind=real_cvprec), parameter :: tolerance_1                               &
                        = 0.01_real_cvprec
real(kind=real_cvprec), parameter :: tolerance_2                               &
                        = 10.0_real_cvprec * min_delta

! Loop counters
integer :: ic, i_cond


!------------------------------------------------------------------------------
! Check the final temperature is consistent with the implicit solution
!------------------------------------------------------------------------------

! Find temperature increment scale...
! Initialise to zero
do ic = 1, cmpr%n_points
  delta(ic) = zero
  work(ic) = zero
end do
! Sum up the T increments from all processes in work; we want the maximum
! amplitude cumulative increment applied, so take max after each
! process has been added on.  Note the sign of the phase-changes
! has been reversed here, as we are working backwards from the
! value of cp_tot after all the increments have already been added on.
do i_cond = 1, n_cond_species_liq
  ! Condensation for liquid species
  call lat_heat_incr( cmpr%n_points, cmpr%n_points, i_phase_change_evp,        &
                      cp_tot, work, dq=dq_cond(:,i_cond) )
  do ic = 1, cmpr%n_points
    delta(ic) = max( delta(ic), abs(work(ic)) )
  end do
end do
do i_cond = n_cond_species_liq+1, n_cond_species
  ! Deposition for ice species
  call lat_heat_incr( cmpr%n_points, cmpr%n_points, i_phase_change_sub,        &
                      cp_tot, work, dq=dq_cond(:,i_cond) )
  do ic = 1, cmpr%n_points
    delta(ic) = max( delta(ic), abs(work(ic)) )
  end do
  ! Melting for ice species
  call lat_heat_incr( cmpr%n_points, cmpr%n_points, i_phase_change_frz,        &
                      cp_tot, work, dq=dq_melt(:,i_cond) )
  do ic = 1, cmpr%n_points
    delta(ic) = max( delta(ic), abs(work(ic)) )
  end do
end do

! Find points where solution for temperature is out by more than tolerance
do ic = 1, cmpr%n_points

  ! Maximum allowable error = either 1% of the increment,
  ! or EPS * value (max of actual and reference), whichever is greater
  max_error(ic) = max( tolerance_1 * delta(ic),                                &
                       tolerance_2 * max( temperature(ic), ref_temp(ic) ) )

  ! Find error in final T;
  ! imp_temp stores T_imp - T_ref, so need to add on T_ref
  final_error(ic) = temperature(ic) - (ref_temp(ic) + imp_temp(ic))

  ! Set flag if error exceeds tolerance
  l_bad(ic) = abs( final_error(ic) ) > max_error(ic)

end do  ! ic = 1, cmpr%n_points

! Call routine to print coordinates and values at "bad" points
call check_bad_values_print(                                                   &
           cmpr, k, l_bad, i_check_imp_consistent,                             &
           "(" // trim(adjustl(call_string)) // ")"            //newline//     &
           "Final temperature from phase_change_solve "                 //     &
           "is inconsistent with the implicit solution."       //newline//     &
           "T, dT, error, tolerance: ",                                        &
           field1=temperature, field2=delta,                                   &
           field3=final_error, field4=max_error )


!------------------------------------------------------------------------------
! Check the final q_vap is consistent with the implicit solution
!------------------------------------------------------------------------------

! Find q_vap increment scale...
! Initialise to zero
do ic = 1, cmpr%n_points
  delta(ic) = zero
  work(ic) = zero
end do
! Sum up the q_vap increments from all processes in work; we want the maximum
! amplitude cumulative increment applied, so take max after each
! process has been added on.
do i_cond = 1, n_cond_species
  do ic = 1, cmpr%n_points
    work(ic) = work(ic) + dq_cond(ic,i_cond)
    delta(ic) = max( delta(ic), abs(work(ic)) )
  end do
end do

! Find points where solution for q_vap is out by more than tolerance
do ic = 1, cmpr%n_points

  ! Maximum allowable error = either 1% of the increment,
  ! or EPS * value (max of actual and reference), whichever is greater
  max_error(ic) = max( tolerance_1 * delta(ic),                                &
                       tolerance_2 * max( q_vap(ic), qsat_liq_ref(ic) ) )

  ! Find error in final q_vap;
  ! imp_q_vap stores q_imp - qsat_ref, so need to add on qsat
  final_error(ic) = q_vap(ic) - (qsat_liq_ref(ic) + imp_q_vap(ic))

  ! Set flag if error exceeds tolerance
  l_bad(ic) = abs( final_error(ic) ) > max_error(ic)

end do  ! ic = 1, cmpr%n_points

! Call routine to print coordinates and values at "bad" points
call check_bad_values_print(                                                   &
           cmpr, k, l_bad, i_check_imp_consistent,                             &
           "(" // trim(adjustl(call_string)) // ")"            //newline//     &
           "Final q_vap from phase_change_solve "                       //     &
           "is inconsistent with the implicit solution."       //newline//     &
           "q, dq, error, tolerance: ",                                        &
           field1=q_vap,       field2=delta,                                   &
           field3=final_error, field4=max_error )


return
end subroutine moist_proc_consistency_2


end module moist_proc_consistency_mod
