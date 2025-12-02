! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module check_negatives_mod

implicit none

contains


! Subroutine calculates the increments for each condensation /
! evaporation process rate,
! checks to avoid ending up with a negative mixing ratio,
! and corrects the increments (and implicit solve coefficients)
! if needed to prevent processes "overshooting" to make negative
! mixing ratios.  Note that where the implicit solve is corrected
! to avoid negative condensate, any melting increments also
! need to be corrected to account for the modified T,q.

subroutine check_negatives( n_points, n_points_super,                          &
                            nc, index_ic,                                      &
                            linear_qs, cp_tot, L_con,L_sub,L_fus,              &
                            temperature, q_vap, q_cond,                        &
                            kq_cond, kt_cond,                                  &
                            coefs_melt, coefs_cond, coefs_cond_m,              &
                            dq_cond, dq_melt, l_melt,                          &
                            coefs_temp, coefs_q_vap,                           &
                            imp_temp, imp_q_vap )

use comorph_constants_mod, only: real_cvprec, cond_params, n_cond_species,     &
                                 n_cond_species_liq, n_cond_species_ice,       &
                                 zero
use linear_qs_mod, only: n_linear_qs_fields, i_ref_temp,                       &
                         i_qsat_liq_ref
use phase_change_coefs_mod, only: n_coefs

use modify_coefs_liq_mod, only: modify_coefs_liq
use modify_coefs_ice_mod, only: modify_coefs_ice
use solve_tq_mod, only: solve_tq
use proc_incr_mod, only: proc_incr
use melt_ctl_mod, only: melt_ctl

implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points in the condensed water species super-array
! (this might be reused on subsequent calls with bigger size
!  than needed, to save having to re-allocate it)
integer, intent(in) :: n_points_super

! Number of points where each condensed water species is non-zero
integer, intent(in) :: nc( n_cond_species )
integer, intent(in) :: index_ic( n_points, n_cond_species )

! Super-array containing qsat at a reference temperature,
! and dqsat/dT, for linearised qsat calculations
real(kind=real_cvprec), intent(in) :: linear_qs                                &
                          ( n_points, n_linear_qs_fields )

! Total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(in) :: cp_tot(n_points)

! Latent heats of condensation, sublimation and fusion
real(kind=real_cvprec), intent(in) :: L_con(n_points)
real(kind=real_cvprec), intent(in) :: L_sub(n_points)
real(kind=real_cvprec), intent(in) :: L_fus(n_points)

! Parcel air temperature and water vapour mixing ratio before
! implicit phase changes
real(kind=real_cvprec), intent(in) :: temperature(n_points)
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Super-array containing condensed water species mixing-ratios
real(kind=real_cvprec), intent(in) :: q_cond                                   &
                            ( n_points_super, n_cond_species )

! Vapour  and heat exchange coefficient for each species
real(kind=real_cvprec), intent(in out) :: kq_cond                              &
                                  ( n_points, n_cond_species )
real(kind=real_cvprec), intent(in out) :: kt_cond                              &
                                  ( n_points, n_cond_species )
! In-Out just in case they need to be adjusted by solve_tq

! Coefficients for melting of each ice species
real(kind=real_cvprec), intent(in out) :: coefs_melt                           &
  ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )

! Coefficients for condensation / evaporation
real(kind=real_cvprec), intent(in out) :: coefs_cond                           &
                         ( n_points, n_coefs, n_cond_species )
! Coefficients for condensation rate if ice is melting
real(kind=real_cvprec), intent(in out) :: coefs_cond_m                         &
   ( n_points, n_coefs, n_cond_species_liq+1 : n_cond_species )
! Note: these are intent inout because their values get swapped
! at points where melting is switched on or off, such that
! coefs_cond always stores the coefficients actually in use

! Increments due to condensation / evaporation and melting
real(kind=real_cvprec), intent(in out) :: dq_cond                              &
                                  ( n_points, n_cond_species )
real(kind=real_cvprec), intent(in out) :: dq_melt                              &
           ( n_points, n_cond_species_liq+1 : n_cond_species )

! Flag for points where each ice species is melting
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


! Number of points of each species where calculations are
! still needed
integer :: nc_check( n_cond_species )
integer :: index_ic_check( n_points, n_cond_species )

! Points where corrections are applied for a given species
integer :: nc_cor
integer :: index_ic_cor( n_points )

! Flag for points where we checked for negative values
logical :: l_checked(n_points)
! Flag for points where we applied corrections, having already
! checked another species, in which case we need to go back and
! check the former species again as the implicit solution for
! T,q will have changed
logical :: l_check_again(n_points)

! Flag for ice species whose melting coefficients have been used to
! modify the coefficients for liquid species they melt into, so-as to
! yield exactly zero condensate mass after melt-source is accounted for.
logical :: l_trunc_mlt ( n_points, n_cond_species_liq+1 : n_cond_species )

! Work array for summing the total source of a liquid species
! from melting of ice species
real(kind=real_cvprec) :: melt_source(n_points)

! Flag for whether any melt sources found
logical :: l_any_melt

! Flag to pass into solve_tq and proc_incr, indicating to always use
! indirect indexing instead of full do-loops where called from here.
! Since we expect the calculations here to nearly always be very sparse,
! there's no point retaining the option to use full do-loops.
logical, parameter :: l_full_do_false = .false.

! Loop counters
integer :: ic, ic2, i_cond, i_liq, i_ice, iter


! Initialise lists of points to check; copy all points where
! the species' mixing-ratios are non-zero
do i_cond = 1, n_cond_species
  nc_check(i_cond) = nc(i_cond)
  index_ic_check( 1:nc(i_cond), i_cond )                                       &
    = index_ic( 1:nc(i_cond), i_cond )
end do

! Initialise flags for ice species melting into liquid species whose
! process rates are truncated to avoid negative condensate
do i_cond = n_cond_species_liq+1, n_cond_species
  if ( nc(i_cond) > 0 ) then
    do ic2 = 1, nc(i_cond)
      l_trunc_mlt(index_ic(ic2,i_cond),i_cond) = .false.
    end do
  end if
end do

! Now, for each condensed water species, we need to  check
! whether the combined increments from condensation / evaporation
! and melting lead to any negative species mixing-ratio values.
! If a negative value is found, the implicit coefficients must
! be corrected so-as to yield exactly zero for that species
! instead.  The implicit solution for T,q must then be
! recomputed with the modified coefficients.  If this is done,
! the condensation and melting increments for all other species
! then need to be recalculated consistent with the updated
! solution.
! When using multiple hydrometeor species, there is a
! slight chance that the heating / drying of the parcel due to
! truncating the evaporation rate to avoid a negative value for
! one species (since we have removed spurious over-evaporation)
! will cause another species to end up with a negative value,
! even if we already checked it and found it to be positive
! before.
! For this reason, if we truncate evaporation or melting
! to avoid a negative value for a latter species but didn't
! do so when we checked a former one (in the order of
! checking below), then we need to repeat the negative value
! check for the former species.
! This is achieved by looping over all species multiple
! times when doing the negative value check.
! The loop must continue until no negative values are found.
! In theory this could get quite expensive, but in practice
! it is very unlikely to have to go round the loop more than a
! couple of times.

! Repeat checks only while any points still need checking, and
! the number of repeats is less than the number of species
! (in theory n_cond_species is the max number of repeats we
!  should ever need, as this allows all the species to be
!  checked in every concievable order)
iter = 0
do while ( maxval(nc_check)>0 .and. iter<n_cond_species )
  iter = iter + 1


  ! Reset checking flags to false at all points where
  ! we're about to check
  do i_cond = 1, n_cond_species
    if ( nc_check(i_cond) > 0 ) then
      do ic2 = 1, nc_check(i_cond)
        ic = index_ic_check(ic2,i_cond)
        l_checked(ic) = .false.
        l_check_again(ic) = .false.
      end do
    end if
  end do



  ! Loop over ice species first
  if ( n_cond_species_ice > 0 ) then
    do i_ice = n_cond_species_liq+1, n_cond_species
      ! If any points of this species to check
      if ( nc_check(i_ice) > 0 ) then

        ! Find points where the current deposition and melting
        ! rates will yield a negative hydrometeor mixing-ratio...
        nc_cor = 0
        do ic2 = 1, nc_check(i_ice)
          ic = index_ic_check(ic2,i_ice)
          if ( q_cond(ic,i_ice) + dq_cond(ic,i_ice)                            &
                                - dq_melt(ic,i_ice) < zero ) then
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
            l_check_again(ic) = l_check_again(ic)                              &
                           .or. l_checked(ic)
          end do

          ! Modify the coefficients for the implicit solve, to
          ! ensure negative q_cond is avoided
          call modify_coefs_ice( n_points, nc_cor, index_ic_cor,               &
                                 L_sub, L_fus, cp_tot,                         &
                                 q_cond(:,i_ice),                              &
                                 dq_cond(:,i_ice),                             &
                                 dq_melt(:,i_ice),                             &
                                 l_melt(:,i_ice),                              &
                                 coefs_cond(:,:,i_ice),                        &
                                 coefs_melt(:,:,i_ice),                        &
                                 coefs_cond_m(:,:,i_ice),                      &
                                 coefs_temp, coefs_q_vap )

          ! Recompute the implicit solve for T,q consistent with
          ! the modified coefficients
          call solve_tq( n_points, l_full_do_false,                            &
                         linear_qs(:,i_ref_temp),                              &
                         linear_qs(:,i_qsat_liq_ref),                          &
                         coefs_temp, coefs_q_vap,                              &
                         temperature, q_vap, imp_temp, imp_q_vap,              &
                         kq_cond, kt_cond,                                     &
                         coefs_cond, coefs_cond_m, coefs_melt,                 &
                         nc=nc_cor, index_ic=index_ic_cor )

          ! Recheck for melting and recalculate melting rates
          call melt_ctl( n_points, nc_check, index_ic_check,                   &
                         linear_qs, cp_tot, L_sub, L_fus,                      &
                         temperature, q_vap, kq_cond, kt_cond,                 &
                         coefs_melt, coefs_cond, coefs_cond_m,                 &
                         dq_melt, l_melt, coefs_temp, coefs_q_vap,             &
                         imp_temp, imp_q_vap,                                  &
                         nc_in=nc_cor, index_ic_in=index_ic_cor,               &
                         l_trunc_mlt=l_trunc_mlt )

          ! Recalculate all the condensation / evaporation
          ! increments consistent with the new implicit solution
          do i_cond = 1, n_cond_species
            call proc_incr( n_points, nc_cor, index_ic_cor,                    &
                            l_full_do_false,                                   &
                            imp_temp, imp_q_vap,                               &
                            coefs_cond(:,:,i_cond),                            &
                            dq_cond(:,i_cond) )
          end do

        end if  ! ( nc_cor > 0 )

        ! Set flag to indicate we checked something at this point
        do ic2 = 1, nc_check(i_ice)
          ic = index_ic_check(ic2,i_ice)
          l_checked(ic) = .true.
        end do

      end if  ! ( nc_check(i_ice) > 0 )
    end do  ! i_ice = n_cond_species_liq+1, n_cond_species
  end if  ! ( n_cond_species_ice > 0 )


  ! Loop over liquid species
  do i_liq = 1, n_cond_species_liq
    ! If any points of this species to check
    if ( nc_check(i_liq) > 0 ) then

      ! Calculate total source of this liquid species from
      ! melting of ice...

      ! Initialise source from melting to zero
      do ic2 = 1, nc_check(i_liq)
        melt_source(ic2) = zero
      end do

      ! Loop over ice species
      l_any_melt = .false.
      if ( n_cond_species_ice > 0 ) then
        do i_ice = n_cond_species_liq+1, n_cond_species
          ! If the current ice species may melt to form the
          ! current liquid species
          if ( cond_params(i_ice)%pt % i_cond_frzmlt                           &
               == i_liq ) then
            ! If any of this ice species present
            if ( nc(i_ice) > 0 ) then
              l_any_melt = .true.
              ! Add on contribution from melting of this species
              do ic2 = 1, nc_check(i_liq)
                ic = index_ic_check(ic2,i_liq)
                melt_source(ic2) = melt_source(ic2)                            &
                                 + dq_melt(ic,i_ice)
              end do
            end if
          end if
        end do
      end if

      ! Find points where the combinatination of evaporation
      ! and melting will lead to a negative mixing ratio for
      ! the current liquid species
      nc_cor = 0
      do ic2 = 1, nc_check(i_liq)
        ic = index_ic_check(ic2,i_liq)
        if ( q_cond(ic,i_liq) + dq_cond(ic,i_liq)                              &
                              + melt_source(ic2) < zero ) then
          nc_cor = nc_cor + 1
          index_ic_cor(nc_cor) = ic
        end if
      end do

      ! If negative value will occur anywhere
      if ( nc_cor > 0 ) then

        ! Set check again flag to true at points where we apply
        ! corrections, having already set the l_checked flag for
        ! an earlier species
        do ic2 = 1, nc_cor
          ic = index_ic_cor(ic2)
          l_check_again(ic) = l_check_again(ic)                                &
                         .or. l_checked(ic)
        end do

        ! Modify the coefficients for the implicit solve, to
        ! ensure negative q_cond is avoided
        call modify_coefs_liq( n_points, nc_cor, index_ic_cor,                 &
                               l_any_melt, i_liq, L_con, cp_tot,               &
                               q_cond(:,i_liq), l_melt, l_trunc_mlt,           &
                               coefs_melt, coefs_cond(:,:,i_liq),              &
                               coefs_temp, coefs_q_vap )

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

        ! Recheck for melting and recalculate melting rates
        if ( n_cond_species_ice > 0 ) then
          if ( maxval( nc_check(n_cond_species_liq+1:                          &
                                n_cond_species) ) > 0 ) then
            call melt_ctl( n_points, nc_check, index_ic_check,                 &
                   linear_qs, cp_tot, L_sub, L_fus,                            &
                   temperature, q_vap, kq_cond, kt_cond,                       &
                   coefs_melt, coefs_cond, coefs_cond_m,                       &
                   dq_melt, l_melt, coefs_temp, coefs_q_vap,                   &
                   imp_temp, imp_q_vap,                                        &
                   nc_in=nc_cor, index_ic_in=index_ic_cor,                     &
                   l_trunc_mlt=l_trunc_mlt )
          end if
        end if

        ! Recalculate all the condensation / evaporation
        ! increments consistent with the new implicit solution
        do i_cond = 1, n_cond_species
          call proc_incr( n_points, nc_cor, index_ic_cor,                      &
                          l_full_do_false,                                     &
                          imp_temp, imp_q_vap,                                 &
                          coefs_cond(:,:,i_cond),                              &
                          dq_cond(:,i_cond) )
        end do

      end if  ! ( nc_cor > 0 )

      ! Set flag to indicate we checked something at this point
      do ic2 = 1, nc_check(i_liq)
        ic = index_ic_check(ic2,i_liq)
        l_checked(ic) = .true.
      end do

    end if  ! ( nc_check(i_liq) > 0 )
  end do  ! i_liq = 1, n_cond_species_liq



  ! Reset the checking lists to only points where l_check_again
  ! has been set to true.
  do i_cond = 1, n_cond_species
    if ( nc_check(i_cond) > 0 ) then
      nc_cor = nc_check(i_cond)
      nc_check(i_cond) = 0
      do ic2 = 1, nc_cor
        ic = index_ic_check(ic2,i_cond)
        if ( l_check_again(ic) ) then
          nc_check(i_cond) = nc_check(i_cond) + 1
          index_ic_check(nc_check(i_cond),i_cond) = ic
        end if
      end do
    end if
  end do


end do  ! while ( maxval(nc_check)>0 .and. iter<n_cond_species )


return
end subroutine check_negatives

end module check_negatives_mod
