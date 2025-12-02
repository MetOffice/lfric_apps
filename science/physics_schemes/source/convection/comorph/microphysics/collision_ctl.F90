! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module collision_ctl_mod

implicit none

contains

! Subroutine to calculate process increments for collisions
! between the different types of condensed water species
subroutine collision_ctl( n_points, n_points_super, nc, index_ic,              &
                          q_loc_cond, n_cond, r_cond, wf_cond,                 &
                          delta_t, rho_dry, rho_wet,                           &
                          q_cond, cp_tot, temperature,                         &
                          dq_frz_cond, kq_cond, kt_cond, l_diags,              &
                          moist_proc_diags, n_points_diag, n_diags,            &
                          diags_super )

use comorph_constants_mod, only: real_cvprec, zero,                            &
                     n_cond_species, n_cond_species_liq,                       &
                     cond_params, l_cv_cf, l_cv_graup,                         &
                     i_cond_rain, i_cond_cf, i_cond_graup
use moist_proc_diags_type_mod, only: moist_proc_diags_type
use collision_rate_mod, only: collision_rate_cmpr
use lat_heat_mod, only: lat_heat_incr, i_phase_change_frz
use ice_rain_to_graupel_mod, only: ice_rain_to_graupel

implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points in the condensed water species super-array
! (this might be reused on subsequent calls with bigger size
!  than needed, to save having to re-allocate it)
integer, intent(in) :: n_points_super

! Number of points where each species is nonzero
integer, intent(in out) :: nc( n_cond_species )
! Index lists used for compressing onto points where the mixing
! ratio of each condensed water species is nonzero
integer, intent(in out) :: index_ic( n_points, n_cond_species )
! Note: these can be modified for graupel, in the case where
! graupel is created by rain-ice collisions at points where there
! was no graupel already

! Local mixing ratio of each condensed water species,
! accounting for fall-out
real(kind=real_cvprec), intent(in) :: q_loc_cond                               &
                           ( n_points, n_cond_species )

! Number concentration per unit dry-mass
real(kind=real_cvprec), intent(in) :: n_cond                                   &
                          ( n_points, n_cond_species )

! Particle radii of each hydrometeor species / m
real(kind=real_cvprec), intent(in) :: r_cond                                   &
                          ( n_points, n_cond_species )

! Fall-speed of each hydrometeor species
real(kind=real_cvprec), intent(in) :: wf_cond                                  &
                          ( n_points, n_cond_species )

! Time interval
real(kind=real_cvprec), intent(in) :: delta_t(n_points)
! Dry-density
real(kind=real_cvprec), intent(in) :: rho_dry(n_points)
! Wet-density
real(kind=real_cvprec), intent(in) :: rho_wet(n_points)

! Condensed water species mixng ratios to be updated with
! increments due to collisions
real(kind=real_cvprec), intent(in out) :: q_cond                               &
                        ( n_points_super, n_cond_species )

! Total heat capacity and temperature of the air; incremented
! by latent heating from riming
real(kind=real_cvprec), intent(in out) :: cp_tot(n_points)
real(kind=real_cvprec), intent(in out) :: temperature(n_points)

! Total amount of freezing onto each ice hydrometeor species
! (includes homogeneous and heterogeneous freezing and riming)
! Needed for the hydrometeor surface heat budget, important for
! determining the melting rate
real(kind=real_cvprec), intent(in out) :: dq_frz_cond                          &
             ( n_points, n_cond_species_liq+1 : n_cond_species )

! Vapour exchange coefficient for each hydrometeor species
real(kind=real_cvprec), intent(in out) :: kq_cond                              &
                          ( n_points, n_cond_species )

! Heat exchange coefficient for each hydrometeor species
real(kind=real_cvprec), intent(in out) :: kt_cond                              &
                          ( n_points, n_cond_species )

! Note: the vapour and heat exchange coefficients have intent
! inout as they need to be modified in the case where new
! graupel is formed by rain-ice collisions at points with
! zero pre-existing graupel.

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


! Compression indices for points where collisions are occuring
integer :: nc_col(n_cond_species)
integer :: index_ic_col( n_points, n_cond_species )

! Compression indices used for reset to avoid negative q_cond
integer :: nc_tmp
integer :: index_ic_tmp(n_points)

! Increments due to collection of one species by each other
! species.  The super-array address for the collected species
! stores the sum-total collection by all other species.
real(kind=real_cvprec) :: dq_col_cond                                          &
                          ( n_points, n_cond_species )

! Rescaling factor applied to increments to avoid negative q_cond
real(kind=real_cvprec) :: factor(n_points)

! Flag for whether any collison processes occur
logical :: l_collect

! Loop counters
integer :: ic, ic2, i_liq, i_cond, i_super


if ( l_diags ) then
  ! If collision process increment diagnostics requested, save
  ! mixing-ratios at start
  do i_cond = 1, n_cond_species
    if ( nc(i_cond) > 0 ) then
      if ( moist_proc_diags % diags_cond(i_cond)%pt                            &
           % dq_col % flag ) then
        ! Extract super-array address
        i_super = moist_proc_diags % diags_cond(i_cond)%pt                     &
                  % dq_col % i_super
        ! Save value before it gets modified by collisions
        do ic2 = 1, nc(i_cond)
          ic = index_ic(ic2,i_cond)
          diags_super(ic,i_super) = q_cond(ic,i_cond)
        end do
      end if
    end if
  end do
end if  ! ( l_diags )


! Loop over all species which have non-zero mixing ratio
! (currently only doing collection of liquids though)
do i_liq = 1, n_cond_species_liq
  if ( nc(i_liq) > 0 ) then

    ! Initialise total collection increment to species i_liq
    ! to zero
    do ic2 = 1, nc(i_liq)
      ic = index_ic(ic2,i_liq)
      dq_col_cond(ic,i_liq) = zero
    end do
    ! Initialise flag for whether or not species i_liq has
    ! experienced any collection.
    l_collect = .false.

    ! Loop over all species after the current one.
    ! Note loop from i_liq+1 ensures all pairs of species
    ! are compared once but not twice (except we don't
    ! currently bother with ice-ice collisions).
    do i_cond = i_liq+1, n_cond_species
      ! Initialise number of points where collection occurs
      nc_col(i_cond) = 0
      if ( nc(i_cond) > 0 ) then

        ! Find points where both species coincide
        if ( nc(i_liq) == n_points ) then
          ! Species i_liq present at all points;
          ! => coinciding points are just i_cond points
          nc_col(i_cond) = nc(i_cond)
          index_ic_col(1:nc_col(i_cond),i_cond)                                &
            = index_ic(1:nc_col(i_cond),i_cond)
        else if ( nc(i_cond) == n_points ) then
          ! Species i_cond present at all points;
          ! => coinciding points are just i_liq points
          nc_col(i_cond) = nc(i_liq)
          index_ic_col(1:nc_col(i_cond),i_cond)                                &
            = index_ic(1:nc_col(i_cond),i_liq)
        else
          ! Both species have partial coverage; find overlap
          do ic2 = 1, nc(i_liq)
            ic = index_ic(ic2,i_liq)
            if ( q_cond(ic,i_cond) > zero ) then
              nc_col(i_cond) = nc_col(i_cond) + 1
              index_ic_col(nc_col(i_cond),i_cond) = ic
            end if
          end do
        end if

        ! If any points:
        if ( nc_col(i_cond) > 0 ) then

          ! Call routine to calculate increment to species
          ! i_liq due to collection by species i_cond
          call collision_rate_cmpr( n_points,                                  &
                 nc_col(i_cond), index_ic_col(:,i_cond),                       &
                 q_loc_cond(:,i_liq), n_cond(:,i_cond),                        &
                 r_cond(:,i_liq), r_cond(:,i_cond),                            &
                 wf_cond(:,i_liq), wf_cond(:,i_cond),                          &
                 delta_t, rho_dry, rho_wet,                                    &
                 cond_params(i_liq)%pt % rho,                                  &
                 cond_params(i_cond)%pt % rho,                                 &
                 cond_params(i_liq)%pt % area_coef,                            &
                 cond_params(i_cond)%pt % area_coef,                           &
                 dq_col_cond(:,i_cond) )

          ! Add on contribution to total increment to species
          ! i_liq due to collection
          do ic2 = 1, nc_col(i_cond)
            ic = index_ic_col(ic2,i_cond)
            dq_col_cond(ic,i_liq) = dq_col_cond(ic,i_liq)                      &
                                  + dq_col_cond(ic,i_cond)
          end do

          ! Set flag indicating some collection has been done
          l_collect = .true.

        end if  ! ( nc_col(i_cond) > 0 )

      end if  ! ( nc(i_cond) > 0 )
    end do  ! i_cond = i_liq+1, n_cond_species


    ! If any collection was done...
    if ( l_collect ) then

      ! Check to avoid collecting more of species i_liq's
      ! mixing ratio than actually exists; rescale the increments
      ! down if this occurs
      nc_tmp = 0
      do ic2 = 1, nc(i_liq)
        ic = index_ic(ic2,i_liq)
        ! If the increment to q_cond(i_liq) is greater than its value
        if ( dq_col_cond(ic,i_liq) > q_cond(ic,i_liq) ) then
          ! Save indices of points
          nc_tmp = nc_tmp + 1
          index_ic_tmp(nc_tmp) = ic
        end if
      end do

      ! If any points were found to have negative q_cond(i_liq):
      if ( nc_tmp > 0 ) then
        do ic2 = 1, nc_tmp
          ic = index_ic_tmp(ic2)
          ! Calculate a scaling factor by which the increments
          ! need to be reduced to avoid creating negative q_cond
          !  = pre-existing mixing-ratio q / increment (dq)
          factor(ic) = q_cond(ic,i_liq) / dq_col_cond(ic,i_liq)
          ! Reset total collection increment to current value
          dq_col_cond(ic,i_liq) = q_cond(ic,i_liq)
        end do
        ! Reduce the collection increments accordingly at
        ! points where we removed too much of species i_liq
        do i_cond = i_liq+1, n_cond_species
          if ( nc_col(i_cond) > 0 ) then
            do ic2 = 1, nc_tmp
              ic = index_ic_tmp(ic2)
              dq_col_cond(ic,i_cond) = dq_col_cond(ic,i_cond)                  &
                                       * factor(ic)
            end do
          end if
        end do
      end if

      ! Decrement the collected species' mixing ratio
      do ic2 = 1, nc(i_liq)
        ic = index_ic(ic2,i_liq)
        q_cond(ic,i_liq) = q_cond(ic,i_liq)                                    &
                         - dq_col_cond(ic,i_liq)
        ! dq_col(i_liq) stores sum of increments from collection
        ! by all other species
      end do

      ! Add on increments for each collecting species...
      do i_cond = i_liq+1, n_cond_species
        if ( nc_col(i_cond) > 0 ) then

          ! Increment the collecting species' mixing ratio
          do ic2 = 1, nc_col(i_cond)
            ic = index_ic_col(ic2,i_cond)
            q_cond(ic,i_cond) = q_cond(ic,i_cond)                              &
                              + dq_col_cond(ic,i_cond)
          end do

          ! If this collection process entails a phase-change
          ! (i.e. riming)
          if ( cond_params(i_cond)%pt % l_ice .and.                            &
               (.not. cond_params(i_liq)%pt % l_ice) ) then
            ! Add on latent heating
            call lat_heat_incr( n_points, nc_col(i_cond),                      &
                                i_phase_change_frz,                            &
                                cp_tot, temperature,                           &
                                index_ic=index_ic_col(:,i_cond),               &
                                dq=dq_col_cond(:,i_cond) )
            ! Increment the total rate of freezing onto
            ! species i_cond
            do ic2 = 1, nc_col(i_cond)
              ic = index_ic_col(ic2,i_cond)
              dq_frz_cond(ic,i_cond) = dq_frz_cond(ic,i_cond)                  &
                                     + dq_col_cond(ic,i_cond)
            end do
          end if

        end if
      end do

      ! If i_liq is rain, and ice & graupel are both on
      if ( i_liq == i_cond_rain .and. l_cv_cf                                  &
                                .and. l_cv_graup ) then
        ! If any collisions between rain and ice
        if ( nc_col(i_cond_cf) > 0 ) then

          ! Call routine to convert collided ice-cloud and rain
          ! into graupel...
          call ice_rain_to_graupel( n_points,                                  &
                 nc_col(i_cond_cf), index_ic_col(:,i_cond_cf),                 &
                 nc(i_cond_graup), index_ic(:,i_cond_graup),                   &
                 dq_col_cond(:,i_cond_cf),                                     &
                 dq_col_cond(:,i_liq), q_cond(:,i_liq),                        &
                 kq_cond(:,i_liq), kt_cond(:,i_liq),                           &
                 n_cond(:,i_cond_cf), n_cond(:,i_liq),                         &
                 q_loc_cond(:,i_cond_cf), q_loc_cond(:,i_liq),                 &
                 q_cond(:,i_cond_cf), q_cond(:,i_cond_graup),                  &
                 dq_frz_cond(:,i_cond_cf),                                     &
                 dq_frz_cond(:,i_cond_graup),                                  &
                 kq_cond(:,i_cond_graup),                                      &
                 kt_cond(:,i_cond_graup) )

        end if  ! ( nc_col(i_cond_cf) > 0 )
      end if  ! ( i_liq == i_cond_rain .and. l_cv_cf
              !                        .and. l_cv_graup )

    end if  ! ( l_collect )

  end if  ! ( nc(i_liq) > 0 )
end do  ! i_liq = 1, n_cond_species_liq


if ( l_diags ) then
  ! If collision process increment diagnostics requested,
  ! subtract saved mixing-ratios from start to compute increments
  do i_cond = 1, n_cond_species
    if ( nc(i_cond) > 0 ) then
      if ( moist_proc_diags % diags_cond(i_cond)%pt                            &
           % dq_col % flag ) then
        ! Extract super-array address
        i_super = moist_proc_diags % diags_cond(i_cond)%pt                     &
                  % dq_col % i_super
        ! Calculate difference
        do ic2 = 1, nc(i_cond)
          ic = index_ic(ic2,i_cond)
          diags_super(ic,i_super) = q_cond(ic,i_cond)                          &
                                  - diags_super(ic,i_super)
        end do
      end if
    end if
  end do
end if  ! ( l_diags )


return
end subroutine collision_ctl


end module collision_ctl_mod
