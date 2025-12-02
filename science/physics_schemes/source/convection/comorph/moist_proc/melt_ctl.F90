! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module melt_ctl_mod

implicit none

contains


! Subroutine to work out which ice species should undergo
! melting, and which ones shouldn't, and update the phase-change
! coefficients accordingly. The melting increment of each ice
! species is also output.

! Method: first see if the solution with no melting yields
! melting of any of the species.  If it does, change
! coefficients to turn on melting of all ice species and
! recompute the implicit solution for T,q.  Then compute the
! melting rates; if any are negative, turn them off again and
! recompute.
! Doing it this way round, when we turn off a melting rate
! (due to it having gone negative), we remove spurious freezing
! and so cool the parcel, making melting even less likely.
! Therefore, there is no way that a melting rate that we already
! turned off should need to be turned back on due to subsequent
! checks on the melting rates of other ice species.
! However, cooling due to turning off negative melting for one
! ice species increases the chance of other melting species
! ending up with a negative melting rate.  Therefore, we require
! some fancy logic to recompute and recheck the melting rates of
! other species each time we switch a species' negative melting
! rate off.

subroutine melt_ctl( n_points, nc, index_ic,                                   &
                     linear_qs, cp_tot, L_sub, L_fus,                          &
                     temperature, q_vap, kq_cond, kt_cond,                     &
                     coefs_melt, coefs_cond, coefs_cond_m,                     &
                     dq_melt, l_melt, coefs_temp, coefs_q_vap,                 &
                     imp_temp, imp_q_vap,                                      &
                     nc_in, index_ic_in, l_trunc_mlt )

use comorph_constants_mod, only: real_cvprec, n_cond_species,                  &
                                 n_cond_species_liq, n_cond_species_ice,       &
                                 melt_check_temp, melt_temp, zero, indi_thresh,&
                                 cond_params
use linear_qs_mod, only: n_linear_qs_fields, i_ref_temp,                       &
                         i_qsat_liq_ref, i_qsat_ice_ref,                       &
                         i_dqsatdT_ice
use phase_change_coefs_mod, only: n_coefs

use calc_cond_temp_mod, only: calc_cond_temp
use toggle_melt_mod, only: toggle_melt
use solve_tq_mod, only: solve_tq
use proc_incr_mod, only: proc_incr

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where each condensed water species is non-zero
! (and indices of those points)
integer, intent(in) :: nc( n_cond_species )
integer, intent(in) :: index_ic( n_points, n_cond_species )

! Super-array containing qsat at a reference temperature,
! and dqsat/dT, for linearised qsat calculations
real(kind=real_cvprec), intent(in) :: linear_qs                                &
                          ( n_points, n_linear_qs_fields )

! Total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(in) :: cp_tot(n_points)

! Latent heats of sublimation and fusion
real(kind=real_cvprec), intent(in) :: L_sub(n_points)
real(kind=real_cvprec), intent(in) :: L_fus(n_points)


! Parcel air temperature and water vapour mixing ratio before
! implicit phase changes
real(kind=real_cvprec), intent(in) :: temperature(n_points)
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Vapour and heat exchange coefficients for each species
real(kind=real_cvprec), intent(in out) :: kq_cond                              &
                                  ( n_points, n_cond_species )
real(kind=real_cvprec), intent(in out) :: kt_cond                              &
                                  ( n_points, n_cond_species )
! In-Out just in case they need to be adjusted by solve_tq

! Coefficients for melting of each ice species
real(kind=real_cvprec), intent(in out) :: coefs_melt                           &
  ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )

! Coefficients for deposition / sublimation
real(kind=real_cvprec), intent(in out) :: coefs_cond                           &
                         ( n_points, n_coefs, n_cond_species )
! Coefficients for deposition rate if ice is melting
real(kind=real_cvprec), intent(in out) :: coefs_cond_m                         &
   ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )
! Note: these are intent inout because their values get swapped
! at points where melting is switched on or off, such that
! coefs_cond always stores the coefficients actually in use


! Increments to each ice species due to melting
real(kind=real_cvprec), intent(in out) :: dq_melt                              &
           ( n_points, n_cond_species_liq+1 : n_cond_species )

! Flag for whether each ice species is melting
logical, intent(in out) :: l_melt                                              &
           ( n_points, n_cond_species_liq+1 : n_cond_species )

! Super-arrays containing coefficients in the implicit formula
! for temperature and q_vap after phase-changes
real(kind=real_cvprec), intent(in out) :: coefs_temp                           &
                                         ( n_points, n_coefs )
real(kind=real_cvprec), intent(in out) :: coefs_q_vap                          &
                                         ( n_points, n_coefs )

! Estimated temperature and water vapour mixing ratio after
! implicit solution of phase-changes (minus reference values)
real(kind=real_cvprec), intent(in out) :: imp_temp(n_points)
real(kind=real_cvprec), intent(in out) :: imp_q_vap(n_points)

! Optional input list of points; if present, calculations
! in this routine are restricted to points in this list
integer, optional, intent(in) :: nc_in
integer, optional, intent(in) :: index_ic_in(n_points)

! Flag for ice species whose melting coefficients have been used to
! modify the coefficients for liquid species they melt into, so-as to
! yield exactly zero condensate mass after melt-source is accounted for.
logical, optional, intent(in out) :: l_trunc_mlt                               &
         ( n_points, n_cond_species_liq+1 : n_cond_species )

! Points where we check for each ice species
integer :: nc_check( n_cond_species_liq+1 : n_cond_species )
integer :: index_ic_check                                                      &
         ( n_points, n_cond_species_liq+1 : n_cond_species )

! Points where we find negative melting rate and have to correct
integer :: nc_cor
integer :: index_ic_cor(n_points)

! Work variable stores number of points in list
integer :: nc_tmp
integer :: index_ic_tmp(n_points)

! Flags for whether to use full do-loops for each species
logical :: l_full_do( n_cond_species_liq+1 : n_cond_species )

! Flag for points where we checked for negative melting rates
logical :: l_checked(n_points)
! Array of flags indicating where melting needs to be turned on
logical :: l_check_again(n_points)
! Also gets used to flag points where we need to check for
! negative melting rates again, due to the implicit solution
! for T,q having changed,

! Temperature of ice hydrometeors (minus the reference temperature)
real(kind=real_cvprec) :: cond_temp(n_points)

! Flag to pass into calc_cond_temp; always called for
! ice species in this routine, never liquid
logical, parameter :: l_ice = .true.

! Flag to pass into proc_incr; always called in indirect indexing mode
! as calculating non-zero melting rates at points where l_melt
! is false would be incorrect
logical, parameter :: l_full_do_false = .false.

! Flag input to toggle_melt
logical :: l_switch_on

! Loop counters
integer :: ic, ic2, i_ice, i_ice2, iter, i_liq_mlt


!----------------------------------------------------------------
! 1) Test for ice species that aren't melting but should be
!----------------------------------------------------------------

! If a list of points has been input
if ( present(index_ic_in) ) then

  ! For each species, find points in this list where:
  ! a) temperature above melt_check_temp
  ! b) ice is present
  ! c) melting is off
  nc_check(:) = 0
  do ic2 = 1, nc_in
    ic = index_ic_in(ic2)
    if ( linear_qs(ic,i_ref_temp) > melt_check_temp ) then
      do i_ice = n_cond_species_liq+1, n_cond_species
        if ( kq_cond(ic,i_ice) > zero ) then
          if ( .not. l_melt(ic,i_ice) ) then
            nc_check(i_ice) = nc_check(i_ice) + 1
            index_ic_check(nc_check(i_ice),i_ice) = ic
          end if
        end if
      end do
    end if
  end do

  ! Initialise flag for switching melting on
  do ic2 = 1, nc_in
    ic = index_ic_in(ic2)
    l_check_again(ic) = .false.
  end do

  ! In this case (melt_ctl being called from check_negatives),
  ! we expect the calculations to always be sparse, so just hardwire
  ! the flags for full-field do-loops to false
  do i_ice = n_cond_species_liq+1, n_cond_species
    l_full_do(i_ice) = .false.
  end do

  ! No list input
else

  ! Check all points where there is condensate
  ! and temperature above melt_check_temp
  ! (assuming melting is off at all points on input)
  do i_ice = n_cond_species_liq+1, n_cond_species
    nc_check(i_ice) = 0
    do ic2 = 1, nc(i_ice)
      ic = index_ic(ic2,i_ice)
      if ( linear_qs(ic,i_ref_temp) > melt_check_temp ) then
        nc_check(i_ice) = nc_check(i_ice) + 1
        index_ic_check(nc_check(i_ice),i_ice) = ic
      end if
    end do
  end do

  ! Initialise flag for switching melting on
  do ic = 1, n_points
    l_check_again(ic) = .false.
  end do

  ! Set flags for using full-field calculations
  do i_ice = n_cond_species_liq+1, n_cond_species
    l_full_do(i_ice) =   real(nc_check(i_ice),real_cvprec)                     &
                       > indi_thresh * real(n_points,real_cvprec)
  end do

end if


! Loop over ice species
do i_ice = n_cond_species_liq+1, n_cond_species
  if ( nc_check(i_ice) > 0 ) then

    ! For each ice species, at points where melting isn't
    ! currently switched on, compute the hydrometeor temperature
    ! and see if it exceeds the melting point.
    ! Set l_check_again to true at points where any of the
    ! not-yet-melting ice species have temperature > Tm and so
    ! need to have melting switched on...

    call calc_cond_temp( n_points, nc_check(i_ice),                            &
                         index_ic_check(:,i_ice), l_full_do(i_ice), l_ice,     &
                         linear_qs(:,i_qsat_liq_ref),                          &
                         linear_qs(:,i_qsat_ice_ref),                          &
                         linear_qs(:,i_dqsatdt_ice),                           &
                         kq_cond(:,i_ice), imp_temp, imp_q_vap,                &
                         coefs_cond(:,:,i_ice), cond_temp )

    ! Set l_check_again to true at points where hydrometeor
    ! temperature exceeds the melting point
    do ic2 = 1, nc_check(i_ice)
      ic = index_ic_check(ic2,i_ice)
      l_check_again(ic) = l_check_again(ic) .or.                               &
                     cond_temp(ic) + linear_qs(ic,i_ref_temp)                  &
                     > melt_temp
      ! Note: cond_temp stores T_cond-T_ref, so need to add T_ref
    end do

  end if
end do

! l_check_again is now true at points where any of the ice species
! need to have melting switched on


!----------------------------------------------------------------
! 2) Adjust the implicit solve for T,q to add on melting
!----------------------------------------------------------------

! Loop over ice species
do i_ice = n_cond_species_liq+1, n_cond_species
  if ( nc_check(i_ice) > 0 ) then

    ! Refine list of points to those where l_check_again is true
    nc_tmp = nc_check(i_ice)
    nc_check(i_ice) = 0
    do ic2 = 1, nc_tmp
      ic = index_ic_check(ic2,i_ice)
      if ( l_check_again(ic) ) then
        nc_check(i_ice) = nc_check(i_ice) + 1
        index_ic_check(nc_check(i_ice),i_ice) = ic
      end if
    end do

    ! If any points of this species which need melting
    ! switched on...
    if ( nc_check(i_ice) > 0 ) then

      if ( .not. present(index_ic_in) ) then
        ! Reset flags for full do-loops
        l_full_do(i_ice) =   real(nc_check(i_ice),real_cvprec)                 &
                           > indi_thresh * real(n_points,real_cvprec)
      end if

      ! Find liquid species that this ice melts into
      i_liq_mlt = cond_params(i_ice)%pt % i_cond_frzmlt

      ! Change the coefficients in the implicit solve for T,q,
      ! consistent with switching melting on
      l_switch_on = .true.
      call toggle_melt( n_points, nc_check(i_ice),                             &
                        index_ic_check(:,i_ice), l_switch_on,                  &
                        L_sub, L_fus, cp_tot,                                  &
                        coefs_melt(:,:,i_ice),                                 &
                        coefs_cond(:,:,i_ice),                                 &
                        coefs_cond_m(:,:,i_ice),                               &
                        coefs_temp, coefs_q_vap,                               &
                        l_melt(:,i_ice),                                       &
                        coefs_cond(:,:,i_liq_mlt), i_ice,                      &
                        l_trunc_mlt=l_trunc_mlt )

      ! Recompute the implicit solution for T,q
      call solve_tq( n_points, l_full_do(i_ice),                               &
                     linear_qs(:,i_ref_temp),                                  &
                     linear_qs(:,i_qsat_liq_ref),                              &
                     coefs_temp, coefs_q_vap,                                  &
                     temperature, q_vap, imp_temp, imp_q_vap,                  &
                     kq_cond, kt_cond,                                         &
                     coefs_cond, coefs_cond_m, coefs_melt,                     &
                     nc=nc_check(i_ice),                                       &
                     index_ic=index_ic_check(:,i_ice) )

    end if

  end if  ! ( nc_check(i_ice) > 0 )
end do  ! i_ice = n_cond_species_liq+1, n_cond_species


!----------------------------------------------------------------
! 3) Calculate the melting increments for all species
!----------------------------------------------------------------

! For the initial call (index_ic_in not input), the lists
! index_ic_check now contain all points where melting is on.
! But for the latter calls inside check_negatives, there
! maybe some points where melting was already on at input
! to this routine; these points have been excluded from the
! lists so far.  Therefore, if index_ic_in is input,
! recreate the lists to pick out points in the index_ic_in
! list where melting is now on
if ( present(index_ic_in) ) then
  do i_ice = n_cond_species_liq+1, n_cond_species
    nc_check(i_ice) = 0
    do ic2 = 1, nc_in
      ic = index_ic_in(ic2)
      if ( l_melt(ic,i_ice) ) then
        nc_check(i_ice) = nc_check(i_ice) + 1
        index_ic_check(nc_check(i_ice),i_ice) = ic
      end if
    end do
  end do
end if

! Loop over ice species where there's a melting rate to calculate
do i_ice = n_cond_species_liq+1, n_cond_species
  if ( nc_check(i_ice) > 0 ) then

    ! Calc melting increment
    call proc_incr( n_points, nc_check(i_ice),                                 &
                    index_ic_check(:,i_ice), l_full_do_false,                  &
                    imp_temp, imp_q_vap,                                       &
                    coefs_melt(:,:,i_ice), dq_melt(:,i_ice) )

  end if
end do


!----------------------------------------------------------------
! 4) Check whether any points / ice species have melting switched
! on but ought to have it switched off.
!----------------------------------------------------------------

! This is tested by seeing if the melting increment has come out
! negative.  Melting should only be left on where the implied
! melting rate is positive, otherwise we'll be imposing
! spurious freezing.

! If a negative melting increment is found, the implicit
! coefficients for T,q must be corrected consistent with turning
! the melting of that species off.  The implicit solution for T,q
! must then be recomputed with the modified coefficients.
! If this is done, the melting increments for all other species
! then need to be recalculated consistent with the updated
! solution.
! When using multiple ice species, there is a slight chance that
! the cooling of the parcel due to switching off spurious
! freezing (when a negative melting rate is found for one
! species) will cause another species to end up with a negative
! melting rate, even if we already checked it and found it to
! be positive before.
! For this reason, if we turn off melting for a latter species
! but didn't turn it off for a former one (in the order of
! checking below), then we need to repeat the negative melting
! increment check for the former species.
! This is achieved by looping over all ice species multiple
! times when doing the negative melting increment check.
! The loop must continue until no negative melting is found.
! In theory this could get quite expensive, but in practice
! it is very unlikely to have to go round the loop more than a
! couple of times.

! Repeat checks only while any points still need checking, and
! the number of repeats is less than the number of species
! (in theory n_cond_species_ice is the max number of repeats we
!  should ever need, as this allows all the species to be
!  checked in every concievable order)
iter = 0
do while ( maxval(nc_check)>0 .and. iter<n_cond_species_ice )
  iter = iter + 1


  ! Reset checking flags to false at all points where
  ! we're about to check
  do i_ice = n_cond_species_liq+1, n_cond_species
    if ( nc_check(i_ice) > 0 ) then
      do ic2 = 1, nc_check(i_ice)
        ic = index_ic_check(ic2,i_ice)
        l_checked(ic) = .false.
        l_check_again(ic) = .false.
      end do
    end if
  end do


  ! Loop over all ice species
  do i_ice = n_cond_species_liq+1, n_cond_species
    ! If any points of this species to check
    if ( nc_check(i_ice) > 0 ) then

      ! Find points where the melting increment is negative
      nc_cor = 0
      do ic2 = 1, nc_check(i_ice)
        ic = index_ic_check(ic2,i_ice)
        if ( dq_melt(ic,i_ice) < zero ) then
          nc_cor = nc_cor + 1
          index_ic_cor(nc_cor) = ic
        end if
      end do

      ! If any points need correcting
      if ( nc_cor > 0 ) then

        ! Set check again flag to true at points where we apply
        ! corrections, having already set the l_checked flag
        ! for an earlier species
        do ic2 = 1, nc_cor
          ic = index_ic_cor(ic2)
          l_check_again(ic) = l_check_again(ic)                                &
                         .or. l_checked(ic)
        end do

        ! Find liquid species that this ice melts into
        i_liq_mlt = cond_params(i_ice)%pt % i_cond_frzmlt

        ! Turn melting off at these points
        l_switch_on = .false.
        call toggle_melt( n_points, nc_cor, index_ic_cor,                      &
                          l_switch_on,                                         &
                          L_sub, L_fus, cp_tot,                                &
                          coefs_melt(:,:,i_ice),                               &
                          coefs_cond(:,:,i_ice),                               &
                          coefs_cond_m(:,:,i_ice),                             &
                          coefs_temp, coefs_q_vap,                             &
                          l_melt(:,i_ice),                                     &
                          coefs_cond(:,:,i_liq_mlt), i_ice,                    &
                          l_trunc_mlt=l_trunc_mlt )

        ! Recompute the implicit solve for T,q consistent with
        ! the modified coefficients
        call solve_tq( n_points, l_full_do_false,                              &
                       linear_qs(:,i_ref_temp),                                &
                       linear_qs(:,i_qsat_liq_ref),                            &
                       coefs_temp, coefs_q_vap,                                &
                       temperature, q_vap, imp_temp, imp_q_vap,                &
                       kq_cond, kt_cond,                                       &
                       coefs_cond, coefs_cond_m, coefs_melt,                   &
                       nc=nc_cor, index_ic=index_ic_cor )

        ! Reset melting increment for current species to zero
        do ic2 = 1, nc_cor
          ic = index_ic_cor(ic2)
          dq_melt(ic,i_ice) = zero
        end do

        ! Recalculate all the other melting rates consistent
        ! with the new implicit solution
        do i_ice2 = n_cond_species_liq+1, n_cond_species
          if ( .not. i_ice2 == i_ice ) then
            ! Find _cor points where other species is melting
            nc_tmp = 0
            do ic2 = 1, nc_cor
              ic = index_ic_cor(ic2)
              if ( l_melt(ic,i_ice2) ) then
                nc_tmp = nc_tmp + 1
                index_ic_tmp(nc_tmp) = ic
              end if
            end do
            ! If any points
            if ( nc_tmp > 0 ) then
              ! Recalculate melting increments
              call proc_incr( n_points, nc_tmp, index_ic_tmp,                  &
                              l_full_do_false,                                 &
                              imp_temp, imp_q_vap,                             &
                              coefs_melt(:,:,i_ice2),                          &
                              dq_melt(:,i_ice2) )
            end if
          end if
        end do

      end if  ! ( nc_cor > 0 )

      ! Set flag to indicate we checked something at this point
      do ic2 = 1, nc_check(i_ice)
        ic = index_ic_check(ic2,i_ice)
        l_checked(ic) = .true.
      end do

    end if  ! ( nc_check(i_ice) > 0 )
  end do  ! i_ice = n_cond_species_liq+1, n_cond_species


  ! Reset the checking lists to only points where l_check_again
  ! has been set to true.
  do i_ice = n_cond_species_liq+1, n_cond_species
    if ( nc_check(i_ice) > 0 ) then
      nc_tmp = nc_check(i_ice)
      nc_check(i_ice) = 0
      do ic2 = 1, nc_tmp
        ic = index_ic_check(ic2,i_ice)
        if ( l_check_again(ic) ) then
          nc_check(i_ice) = nc_check(i_ice) + 1
          index_ic_check(nc_check(i_ice),i_ice) = ic
        end if
      end do
    end if
  end do


end do  ! while (maxval(nc_check)>0.and.iter<n_cond_species_ice)


return
end subroutine melt_ctl

end module melt_ctl_mod
