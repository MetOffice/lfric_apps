! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module microphysics_1_mod

implicit none

contains

! This is the main subroutine for BOTEMS
! (Back-Of-The-Envelope Microphysical Scheme)
!
! Subroutine to calculate microphysical properties
! and the following explicit processes:
!  - Activation of new condensation
!  - Ice nucleation
!  - Collision processes (accretion and riming).
!
! This routine also outputs coefficients for the exchange
! of water vapour and heat between each condensed water
! species and the surrounding air.  These are passed out
! for use in calculating an implicit solution of the
! phase-change processes (condensation/evaporation,
! deposition/sublimation, and melting).
!
! Another output from this routine is the fall-speed for each
! condensed water species, used to calculate the fall-out
! flux of each species from the current level or parcel
! (note the sedimentation itself is not done in here).
!
! Note that some microphysical processes need to be calculated
! after the phase-changes (e.g. autoconversion of ice to graupel
! may require the vapour deposition rate).
! Therefore, autoconversion processes are computed in a 2nd
! BOTEMS subroutine (microphysics_2) which is called after the
! phase-change and sedimentation routines.

subroutine microphysics_1( n_points, n_points_super,  nc, index_ic,            &
                           ref_temp, qsat_liq_ref, dqsatdT_liq,                &
                           delta_t, vert_len, rho_dry, rho_wet,                &
                           cp_tot, temperature, q_vap, q_cond,                 &
                           q_loc_cond, wf_cond, kq_cond, kt_cond,              &
                           dq_frz_cond, l_diags, moist_proc_diags,             &
                           n_points_diag, n_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, cond_params, zero,               &
                     n_cond_species,                                           &
                     n_cond_species_liq, n_cond_species_ice
use moist_proc_diags_type_mod, only: moist_proc_diags_type

use activate_cond_mod, only: activate_cond
use ice_nucleation_mod, only: ice_nucleation
use calc_cond_properties_mod, only: calc_cond_properties_cmpr
use collision_ctl_mod, only: collision_ctl

implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points in the condensed water species super-array
! (this might be reused on subsequent calls with bigger size
!  than needed, to save having to re-allocate it)
integer, intent(in) :: n_points_super


! Number of points where each hydrometeor species is non-zero
integer, intent(in out) :: nc(n_cond_species)
! Indices of those points
integer, intent(in out) :: index_ic(n_points,n_cond_species)
! (these are intent inout because homogeneous freezing may entirely
!  remove liquid from some points)

! Reference temperature used for linearised qsat calculations
! and ice nucleation
real(kind=real_cvprec), intent(in) :: ref_temp(n_points)

! Saturation water vapour mixing ratio qsat and dqsat/dT
! w.r.t. liquid water
real(kind=real_cvprec), intent(in) :: qsat_liq_ref(n_points)
real(kind=real_cvprec), intent(in) :: dqsatdT_liq(n_points)

! Time interval for converting process rates to increments.
real(kind=real_cvprec), intent(in) :: delta_t(n_points)

! Vertical length-scale of the parcel.
real(kind=real_cvprec), intent(in) :: vert_len(n_points)

! Air dry-density and wet-density
real(kind=real_cvprec), intent(in) :: rho_dry(n_points)
real(kind=real_cvprec), intent(in) :: rho_wet(n_points)

! Total heat capacity of the air per unit dry-mass
real(kind=real_cvprec), intent(in out) :: cp_tot(n_points)

! Air temperature and water-vapour mixing ratio
! (updated in here by explicit increments)
real(kind=real_cvprec), intent(in out) :: temperature(n_points)
real(kind=real_cvprec), intent(in out) :: q_vap(n_points)

! Super-array containing all the condensed water sepcies
! mixing-ratios to be updated
real(kind=real_cvprec), intent(in out) :: q_cond                               &
                              ( n_points_super, n_cond_species )

! Initial guess end-of-timestep mixing ratio for each hydrometeor
! species, accounting for implicit solution of fall-out
real(kind=real_cvprec), intent(out) :: q_loc_cond                              &
                                   ( n_points, n_cond_species )

! Fall-speed of each hydrometeor species
real(kind=real_cvprec), intent(out) :: wf_cond                                 &
                                   ( n_points, n_cond_species )

! Coefficients for exchange of water vapour and heat between
! each condensed water species and the surrounding air
real(kind=real_cvprec), intent(out) :: kq_cond                                 &
                                   ( n_points, n_cond_species )
real(kind=real_cvprec), intent(out) :: kt_cond                                 &
                                   ( n_points, n_cond_species )

! Total amount of freezing onto each ice hydrometeor species
! (includes homogeneous and heterogeneous freezing and riming)
! Needed for the hydrometeor surface heat budget, important for
! determining the melting rate
real(kind=real_cvprec), intent(out) :: dq_frz_cond                             &
            ( n_points, n_cond_species_liq+1 : n_cond_species )

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


! Store for delta_t / vert_len, used in implicit fall-out
! calculations in calc_cond_properties
real(kind=real_cvprec) :: dt_over_lz(n_points)

! Number concentration per unit dry-mass
real(kind=real_cvprec) :: n_cond                                               &
                          ( n_points, n_cond_species )

! Particle radii of each hydrometeor species / m
real(kind=real_cvprec) :: r_cond                                               &
                          ( n_points, n_cond_species )

! Loop counters
integer :: ic, i_liq, i_ice, i_cond


!----------------------------------------------------------------
! 1) Activation of new condensation
!----------------------------------------------------------------

! Loop over activatable liquid condensed water species
do i_liq = 1, n_cond_species_liq
  if ( cond_params(i_liq)%pt % r_min > zero ) then

    ! Call routine to compute activation
    call activate_cond( n_points, nc(i_liq), index_ic(:,i_liq),                &
                        qsat_liq_ref, dqsatdT_liq, ref_temp,                   &
                        temperature, q_vap, q_cond(:,i_liq) )

  end if
end do


! If there are now any species with non-zero mixing-ratio...
if ( maxval(nc) > 0 ) then

  !--------------------------------------------------------------
  ! 2) Freezing by ice nucleation...
  !--------------------------------------------------------------

  ! Done at this point so that any liquid cloud (existing or
  ! added by activation above) below the homogeneous freezing
  ! threshold will be instantly frozen.  Crucially, the
  ! calculation of explicit coefficients for other processes
  ! comes after this, so will not see any liquid
  ! (and therefore can't create any more liquid at points
  !  where T is below homnuc_temp, eg by condensation)

  ! Initialise freezing rates to zero
  do i_ice = n_cond_species_liq+1, n_cond_species
    do ic = 1, n_points
      dq_frz_cond(ic,i_ice) = zero
    end do
  end do

  ! Can only freeze if at least one liquid and ice species is on
  if ( n_cond_species_ice > 0 ) then

    ! For each liquid species
    do i_liq = 1, n_cond_species_liq
      if ( nc(i_liq) > 0 ) then

        ! Select the ice species which this liquid species
        ! will freeze into
        i_ice = cond_params(i_liq)%pt % i_cond_frzmlt

        ! Call routine to do freezing from species
        ! i_liq to species i_cond_frz
        call ice_nucleation( n_points,                                         &
                             nc(i_liq), index_ic(:,i_liq),                     &
                             nc(i_ice), index_ic(:,i_ice),                     &
                             delta_t, ref_temp,                                &
                             q_cond(:,i_liq), q_cond(:,i_ice),                 &
                             temperature, cp_tot,                              &
                             dq_frz_cond(:,i_ice), l_diags,                    &
                             i_liq, i_ice, moist_proc_diags,                   &
                             n_points_diag, n_diags, diags_super )
        ! Note: currently using q_cond to calculate the
        ! heterogeneous freezing rate in here, but q_cond maybe
        ! more of a numerical rather than physical quantity at
        ! this point, as fall-in has been added to q_cond but
        ! fall-out has not.  This should be fine provided that
        ! the heter freeze rate comes out very small and is
        ! dwarfed by other ice formation processes.  But if
        ! we wish to use a heter freeze formulation where the
        ! actual rate of heter freeze is significant, we'll
        ! need to move heter freeze to after calc_cond_properties
        ! and compute the rate using q_loc_cond (which accounts
        ! for fall-out) instead.

      end if
    end do  ! i_liq = 1, n_cond_species_liq

  end if  ! ( n_cond_species_ice > 0 )


  !--------------------------------------------------------------
  ! 3) Calculate local properties of each hydrometeor species,
  !    implicitly accounting for precipitation fall-out
  !--------------------------------------------------------------

  ! Set factor delta_t / vert_len used in precip
  ! fall-out calculation
  do ic = 1, n_points
    dt_over_lz(ic) = delta_t(ic) / vert_len(ic)
  end do

  ! Loop over all condensed water species
  do i_cond = 1, n_cond_species

    ! Initialise outputs to zero
    do ic = 1, n_points
      q_loc_cond(ic,i_cond) = zero
      n_cond(ic,i_cond)     = zero
      r_cond(ic,i_cond)     = zero
      wf_cond(ic,i_cond)    = zero
      kq_cond(ic,i_cond)    = zero
      kt_cond(ic,i_cond)    = zero
    end do

    ! If any points
    if ( nc(i_cond) > 0 ) then

      ! Routine implicitly solves the fall-speed / fall-out
      ! and number concentration / particle radius relationship,
      ! and calculates moisture and heat exchange coefficients
      call calc_cond_properties_cmpr(                                          &
             n_points, nc(i_cond), index_ic(:,i_cond),                         &
             cond_params(i_cond)%pt, ref_temp, q_cond(:,i_cond),               &
             rho_dry, rho_wet, dt_over_lz,                                     &
             q_loc_cond(:,i_cond), n_cond(:,i_cond),                           &
             r_cond(:,i_cond), wf_cond(:,i_cond),                              &
             kq_cond(:,i_cond), kt_cond(:,i_cond) )

    end if
  end do


  !--------------------------------------------------------------
  ! 4) Collision processes:
  !--------------------------------------------------------------

  ! These are treated explicitly, with increments added on inside
  ! this routine:
  call collision_ctl( n_points, n_points_super, nc, index_ic,                  &
                      q_loc_cond, n_cond, r_cond, wf_cond,                     &
                      delta_t, rho_dry, rho_wet,                               &
                      q_cond, cp_tot, temperature,                             &
                      dq_frz_cond, kq_cond, kt_cond, l_diags,                  &
                      moist_proc_diags, n_points_diag, n_diags, diags_super )


end if  ! ( maxval(nc) > 0 )


return
end subroutine microphysics_1

end module microphysics_1_mod
