! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module moist_proc_conservation_mod

implicit none

! Number of conserved variables to check
integer, parameter :: n_cons_vars = 5

! Addresses in the conserved variables super-array
integer, parameter :: i_heat_tot = 1
integer, parameter :: i_q_tot = 2
integer, parameter :: i_mom_u = 3
integer, parameter :: i_mom_v = 4
integer, parameter :: i_mom_w = 5

! Name of each conserved variable, for error reporting
integer, parameter :: var_names_len = 11
character(len=var_names_len), parameter :: cons_var_names(n_cons_vars)         &
 = [ "total heat ", "total water", "u momentum ", "v momentum ", "w momentum "]

! Max length of string containing coordinates and values of
! any non-conserved values found
integer, parameter :: str_len = 200

contains


! Subroutine to check that moist_proc conserves heat,
! total water and momentum per unit dry-mass

subroutine moist_proc_conservation( cmpr, k, call_string, n_points_super,      &
             dt_over_rhod_lz,                                                  &
             q_cond, q_vap, temperature, wind_u, wind_v, wind_w,               &
             flux_cond, fall_temp,                                             &
             fall_wind_u, fall_wind_v, fall_wind_w,                            &
             cons_vars )

use comorph_constants_mod, only: real_cvprec, zero, one, n_cond_species,       &
                     n_cond_species_liq, min_delta, min_float,                 &
                     cp_dry, cp_vap, cp_liq, cp_ice,                           &
                     L_con_0, L_sub_0, newline, name_length,                   &
                     i_check_conservation
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

! Number of points in the condensed water species super-array
! (this might be reused on subsequent calls with bigger size
!  than needed, to save having to re-allocate it)
integer, intent(in) :: n_points_super

! Factor delta_t / ( rho_dry vert_len ), needed
! for converting fall-fluxes to mixing-ratio increments
real(kind=real_cvprec), intent(in) :: dt_over_rhod_lz(cmpr%n_points)

! Current condensed water mixing ratios
real(kind=real_cvprec), intent(in) :: q_cond                                   &
                                      ( n_points_super, n_cond_species )
! Current water-vapour mixing-ratio
real(kind=real_cvprec), intent(in) :: q_vap(cmpr%n_points)
! Current temperature and wind vector components
real(kind=real_cvprec), intent(in) :: temperature(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: wind_u(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: wind_v(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: wind_w(cmpr%n_points)

! Flux of condensed water species in or out
real(kind=real_cvprec), intent(in) :: flux_cond                                &
                                      ( cmpr%n_points, n_cond_species )
! Temperature and velocity of the transfered condensed water
real(kind=real_cvprec), intent(in) :: fall_temp(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: fall_wind_u(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: fall_wind_v(cmpr%n_points)
real(kind=real_cvprec), intent(in) :: fall_wind_w(cmpr%n_points)

! Super-array to store total heat, water mixing-ratio and
! momentum per unit dry-mass
real(kind=real_cvprec), allocatable, intent(in out) :: cons_vars(:,:)
! When this routine is called at the start of moist_proc,
! these will not be allocated; this routine then allocates
! them and uses them to store the conserved variables at
! the start of moist_proc.
! When this routine is called at the end of moist_proc,
! these are then compared with updated values calculated
! using the fields at the end.  If they differ by more
! than some tolerance (which depends on the precision),
! a warning is raised.

! Work array for calculating new values of conserved variables
real(kind=real_cvprec) :: cons_vars_new(cmpr%n_points,n_cons_vars)

! Mixing-ratio increment due to precip fall
real(kind=real_cvprec) :: dq_cond_fall(cmpr%n_points)

! Maximum error expected just from rounding etc
real(kind=real_cvprec) :: max_error(cmpr%n_points)

! Actual conservation error found
real(kind=real_cvprec) :: cons_error(cmpr%n_points)

! Logical flags set to true where non-conservation detected
logical :: l_bad(cmpr%n_points)

! Heat capacity and latent heat for a condensed-water species
real(kind=real_cvprec) :: cp_cond
real(kind=real_cvprec) :: l_0

! Factor used for computing numerical tolerance of conservation check
real(kind=real_cvprec), parameter :: tolerance = 10.0_real_cvprec

! Loop counters
integer :: ic, i_cond, i_cv


! Initialise conserved variables to zero
do i_cv = 1, n_cons_vars
  do ic = 1, cmpr%n_points
    cons_vars_new(ic,i_cv) = zero
  end do
end do


! Compute contributions for each condensed water species...
do i_cond = 1, n_cond_species

  ! Select heat capacity and latent heat for the current species
  if ( i_cond > n_cond_species_liq ) then
    cp_cond = cp_ice
    l_0 = L_sub_0
  else
    cp_cond = cp_liq
    l_0 = L_con_0
  end if

  ! Add on contributions from existing condensed water
  do ic = 1, cmpr%n_points

    ! Add heat-capacity and latent heat terms to existing air
    cons_vars_new(ic,i_heat_tot) = cons_vars_new(ic,i_heat_tot)                &
       + q_cond(ic,i_cond) * ( cp_cond * temperature(ic) - l_0 )

    ! Add to total water
    cons_vars_new(ic,i_q_tot) = cons_vars_new(ic,i_q_tot)                      &
                              + q_cond(ic,i_cond)

    ! Add to total momentum
    cons_vars_new(ic,i_mom_u) = cons_vars_new(ic,i_mom_u)                      &
                              + q_cond(ic,i_cond) * wind_u(ic)
    cons_vars_new(ic,i_mom_v) = cons_vars_new(ic,i_mom_v)                      &
                              + q_cond(ic,i_cond) * wind_v(ic)
    cons_vars_new(ic,i_mom_w) = cons_vars_new(ic,i_mom_w)                      &
                              + q_cond(ic,i_cond) * wind_w(ic)

  end do

  ! Convert fall-flux to mixing ratio increment
  do ic = 1, cmpr%n_points
    dq_cond_fall(ic) = flux_cond(ic,i_cond) * dt_over_rhod_lz(ic)
  end do

  ! Add contributions from the falling in/out condensed water
  do ic = 1, cmpr%n_points

    cons_vars_new(ic,i_heat_tot) = cons_vars_new(ic,i_heat_tot)                &
          + dq_cond_fall(ic) * ( cp_cond * fall_temp(ic) - l_0 )

    cons_vars_new(ic,i_q_tot) = cons_vars_new(ic,i_q_tot)                      &
                              + dq_cond_fall(ic)

    cons_vars_new(ic,i_mom_u) = cons_vars_new(ic,i_mom_u)                      &
                              + dq_cond_fall(ic) *fall_wind_u(ic)
    cons_vars_new(ic,i_mom_v) = cons_vars_new(ic,i_mom_v)                      &
                              + dq_cond_fall(ic) *fall_wind_v(ic)
    cons_vars_new(ic,i_mom_w) = cons_vars_new(ic,i_mom_w)                      &
                              + dq_cond_fall(ic) *fall_wind_w(ic)

  end do

end do  ! i_cond = 1, n_cond_species


! Add on contributions from dry air and water-vapour
! (these are usually much larger than the contributions from
!  condensed water, so are added on last to improve the
!  precision of the condensed water calculations)
do ic = 1, cmpr%n_points

  cons_vars_new(ic,i_heat_tot) = cons_vars_new(ic,i_heat_tot)                  &
               + temperature(ic) * ( cp_dry + cp_vap*q_vap(ic) )

  cons_vars_new(ic,i_q_tot) = cons_vars_new(ic,i_q_tot)                        &
                            + q_vap(ic)

  cons_vars_new(ic,i_mom_u) = cons_vars_new(ic,i_mom_u)                        &
                            + wind_u(ic) * ( one + q_vap(ic) )
  cons_vars_new(ic,i_mom_v) = cons_vars_new(ic,i_mom_v)                        &
                            + wind_v(ic) * ( one + q_vap(ic) )
  cons_vars_new(ic,i_mom_w) = cons_vars_new(ic,i_mom_w)                        &
                            + wind_w(ic) * ( one + q_vap(ic) )

end do


! If this is the first call, at the start of moist_proc
! (deciphered by checking whether the work array is allocated)
if ( .not. allocated(cons_vars) ) then

  ! Allocate conserved variable store array
  allocate( cons_vars(cmpr%n_points,n_cons_vars) )

  ! Copy current values into the output array
  do i_cv = 1, n_cons_vars
    do ic = 1, cmpr%n_points
      cons_vars(ic,i_cv) = cons_vars_new(ic,i_cv)
    end do
  end do

  ! If this is the 2nd call, at the end of moist_proc
else  ! ( allocated(cons_vars) )

  ! Compare new values of conserved variables with the old ones,
  ! and raise an error if they differ by more than the tolerance

  ! Loop over conserved variables
  do i_cv = 1, n_cons_vars

    ! Compute conservation error
    do ic = 1, cmpr%n_points
      cons_error(ic) = cons_vars_new(ic,i_cv) - cons_vars(ic,i_cv)
    end do

    ! Estimate max expected error just from rounding etc;
    ! this should scale with the magnitude of the conserved variable
    ! times epsilon, times a tolerance factor to account for the fact
    ! that moist_proc will have sequentially added multiple different
    ! increments, over-which rounding errors may have accumulated.
    do ic = 1, cmpr%n_points
      max_error(ic) = tolerance * min_delta                                    &
                      * max( abs(cons_vars(ic,i_cv)),                          &
                             abs(cons_vars_new(ic,i_cv)) )
      ! If the cons_vars get VERY small, max_error can go to zero
      ! due to floating-point underflow, resulting in wrongly
      ! diagnosing non-conservation.  Prevent this by not
      ! letting max_error fall below tolerance * tiny
      max_error(ic) = max( max_error(ic), tolerance * min_float )
    end do

    ! Look for non-conservation exceeding tolerance
    do ic = 1, cmpr%n_points
      l_bad(ic) = abs(cons_error(ic)) > max_error(ic)
    end do

    ! Call routine to print coordinates and values at "bad" points
    call check_bad_values_print(                                               &
           cmpr, k, l_bad, i_check_conservation,                               &
           "(" // trim(adjustl(call_string)) // ")"            //newline//     &
           "Non-conservation of " // trim(adjustl(cons_var_names(i_cv)))//     &
           " detected in moist_proc."                          //newline//     &
           "value, error, tolerance: ",                                        &
           field1=cons_vars(:,i_cv), field2=cons_error, field3=max_error )

  end do  ! i_cv = 1, n_cons_vars

  ! Deallocate work array
  deallocate( cons_vars )

end if  ! ( allocated(cons_vars) )


return
end subroutine moist_proc_conservation

end module moist_proc_conservation_mod
