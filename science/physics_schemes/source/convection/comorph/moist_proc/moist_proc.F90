! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module moist_proc_mod

implicit none

contains


! Routine to update the thermodynamic variables due to
! water phase changes and microphysical processes.
subroutine moist_proc( n_points, n_points_super, linear_qs,                    &
                       delta_z, delta_t, vert_len, wind_w_excess,              &
                       dt_over_rhod_lz, pressure, prev_temp,                   &
                       fallin_wind_u, fallin_wind_v,                           &
                       fallin_wind_w, fallin_temp,                             &
                       wind_u, wind_v, wind_w,                                 &
                       temperature, q_vap, q_cond,                             &
                       flux_cond, cmpr, k, call_string, l_diags,               &
                       moist_proc_diags, n_points_diag, n_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, zero, one, cond_params,          &
                     n_cond_species, n_cond_species_liq,                       &
                     i_check_conservation, i_check_bad_values_cmpr,            &
                     i_check_bad_none, name_length, min_float
use linear_qs_mod, only: n_linear_qs_fields, i_ref_temp,                       &
                         i_qsat_liq_ref, i_dqsatdT_liq
use cmpr_type_mod, only: cmpr_type

use moist_proc_conservation_mod, only: moist_proc_conservation
use set_cp_tot_mod, only: set_cp_tot
use calc_rho_dry_mod, only: calc_rho_dry
use calc_q_tot_mod, only: calc_q_tot
use fall_in_mod, only: fall_in
use fall_out_mod, only: fall_out_flux
use phase_change_solve_mod, only: phase_change_solve
use microphysics_1_mod, only: microphysics_1
use microphysics_2_mod, only: microphysics_2
use check_bad_values_mod, only: check_bad_values_cmpr

use moist_proc_diags_type_mod, only: moist_proc_diags_type

implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points in the condensed water species super-array
! (this might be reused on subsequent calls with bigger size
!  than needed, to save having to re-allocate it)
integer, intent(in) :: n_points_super
! Number of points in the diagnostics super-array
integer, intent(in) :: n_points_diag

! Super-array containing qsat at a reference temperature,
! and dqsat/dT, for linearised qsat calculations
real(kind=real_cvprec), intent(in) :: linear_qs                                &
                          ( n_points, n_linear_qs_fields )
! The contained reference temperature is also used for
! calculating ice nucleation.

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
! so delta_t = delta_z/wind_w_excess
real(kind=real_cvprec), intent(in) :: delta_t(n_points)

! Vertical length-scale of the parcel.
! For Eulerian calculations, this is the thickness of the
! current model-level,
! but for Lagrangian ascents it is an independent vertical
! length-scale for the parcel.
! It is used to calculate the fraction of precip that falls
! out of the current model-level / parcel during the
! time-interval delta_t.
real(kind=real_cvprec), intent(in) :: vert_len(n_points)

! Vertical wind velocity
! For Euelerian calculations, this is the vertical wind-speed.
! For Lagrangian ascents, it is the vertical velocity of the
! parcel relative to the environment, such that
! wind_w_excess = delta_z/delta_t
real(kind=real_cvprec), intent(in) :: wind_w_excess(n_points)
! Note that wind_w_excess does not modify the fall velocities
! used in sedimentation.
! For Lagrangian ascents, the frame of reference moves with
! the vertical wind, and sedimentation is relative to the
! moving parcel.  The primary impact of wind_w_excess is its
! affect on the time interval delta_t.
! For Eulerian calculations, it is assumed that a separate
! advection scheme is in operation to add on vertical
! transport by the wind.  In this case, wind_w_excess only
! affects the solution via its role in the
! (wind_w_excess - fall_speed) dT/dz term in the
! diagnostic formula for the hydrometeor temperatures in
! phase_change_solve.

! Factor delta_t / ( rho_dry vert_len ), used to convert precip
! flux through the top and bottem of the parcel / level into
! condensed-water mixing-ratio increments within
! the parcel / level.
real(kind=real_cvprec), intent(in) :: dt_over_rhod_lz(n_points)

! Current ambient pressure
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Parcel temperature at the previous model-level
real(kind=real_cvprec), intent(in) :: prev_temp(n_points)

! Wind velocity and temperature of the air from which any
! falling-in precipitation originates.  Precipitation
! carries this momentum and heat with it when it falls
! into the present air.
! For Eulerian calculations, these are the wind and temperature
! of the previous layer (i.e. fallin_temp = prev_temp).
! For Lagrangian ascents, they correspond to the environment
! properties on the current layer.
! Modification of the air's temperature and winds by
! transfer of precipitation is expected to be small,
! but needs to be included to yield accurate conservation
real(kind=real_cvprec), intent(in) :: fallin_wind_u(n_points)
real(kind=real_cvprec), intent(in) :: fallin_wind_v(n_points)
real(kind=real_cvprec), intent(in) :: fallin_wind_w(n_points)
real(kind=real_cvprec), intent(in) :: fallin_temp(n_points)

! Properties of the current air:
! in:  -updated by advection and mixing, but not yet any phase
!       changes or microphysical processes.
! out: -fully updated by phase changes and microphysical
!       processes
real(kind=real_cvprec), intent(in out) :: wind_u(n_points)
real(kind=real_cvprec), intent(in out) :: wind_v(n_points)
real(kind=real_cvprec), intent(in out) :: wind_w(n_points)
real(kind=real_cvprec), intent(in out) :: temperature(n_points)
real(kind=real_cvprec), intent(in out) :: q_vap(n_points)
! Super-array containing all the condensed water species
real(kind=real_cvprec), intent(in out) :: q_cond                               &
                              ( n_points_super, n_cond_species )

! Super-array containing fall-fluxes of all condened
! water species / kg m-2 s-1
! in:  - fluxes falling into the current air
! out: - fluxes falling out of the current air
real(kind=real_cvprec), intent(in out) :: flux_cond                            &
                                   ( n_points, n_cond_species )

! Stuff used to make error messages more informative:
! Structure storing compression indices of the current points
type(cmpr_type), intent(in) :: cmpr
! Current model-level index
integer, intent(in) :: k
! Description of which call to moist_proc this is
character(len=name_length), intent(in) :: call_string

! Master switch for whether or not to calculate any diagnostics
! (moist_proc might be called within an iteration, where
!  we only want to output diags for the final call).
logical, intent(in) :: l_diags
! Structure containing diagnostic flags etc.
type(moist_proc_diags_type), intent(in) :: moist_proc_diags
! Total number of diagnostics in the super-array
integer, intent(in) :: n_diags
! Super-array to store all the diagnostics
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                         ( n_points_diag, n_diags )


! Super-arrays containing properties of each condensed water
! species, which are output by the microphysics call and required
! by the implicit phase-change solve...

! Initial guess end-of-timestep mixing ratio for each hydrometeor
! species, accounting for implicit solution of fall-out
real(kind=real_cvprec) :: q_loc_cond                                           &
                          ( n_points, n_cond_species )

! Fall-speed of each hydrometeor species
real(kind=real_cvprec) :: wf_cond                                              &
                          ( n_points, n_cond_species )

! Vapour exchange coefficient for each hydrometeor species
real(kind=real_cvprec) :: kq_cond                                              &
                          ( n_points, n_cond_species )

! Heat exchange coefficient for each hydrometeor species
real(kind=real_cvprec) :: kt_cond                                              &
                          ( n_points, n_cond_species )

! Total amount of freezing onto each ice hydrometeor species
! (includes homogeneous and heterogeneous freezing and riming)
! Needed for the hydrometeor surface heat budget, important for
! determining the melting rate
real(kind=real_cvprec) :: dq_frz_cond                                          &
             ( n_points, n_cond_species_liq+1 : n_cond_species )


! Total heat capacity (per unit dry-mass);
! this gets incremented by phase-changes
real(kind=real_cvprec) :: cp_tot(n_points)

! Dry-density
real(kind=real_cvprec) :: rho_dry(n_points)

! Wet density
! (also used to store total-water q_vap + sum q_cond)
real(kind=real_cvprec) :: rho_wet(n_points)

! Work array used for conservation checks
real(kind=real_cvprec), allocatable :: cons_vars(:,:)

! Index lists used for compressing onto points where the mixing
! ratio of each condensed water species is nonzero
integer :: index_ic( n_points, n_cond_species )
! Number of points where each species is nonzero
integer :: nc( n_cond_species )

! Name of a field (for error message)
character(len=name_length) :: field_name
! Flag for whether field is positive-only
logical :: l_positive
! Description of where we are in the code, for error messages
character(len=name_length) :: where_string

! Loop counter
integer :: ic, ic2, i_cond, i_super


! Check inputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  where_string = "Start of moist_proc call for "              //               &
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
end if


!----------------------------------------------------------------
! 1) Initialisations
!----------------------------------------------------------------

if ( i_check_conservation > i_check_bad_none ) then
  ! Calculate conserved variables on entry to this routine
  call moist_proc_conservation( cmpr, k, call_string, n_points_super,          &
         dt_over_rhod_lz,                                                      &
         q_cond, q_vap, temperature, wind_u, wind_v, wind_w,                   &
         flux_cond, fallin_temp,                                               &
         fallin_wind_u, fallin_wind_v, fallin_wind_w,                          &
         cons_vars )
end if

! Set the total heat capacity of the parcel
call set_cp_tot( n_points, n_points_super,                                     &
                 q_vap, q_cond, cp_tot )

! Set total-water mixing-ratio (stored in rho_wet)
call calc_q_tot( n_points, n_points_super,                                     &
                 q_vap, q_cond, rho_wet )


!----------------------------------------------------------------
! 2) Fall-in of each hydrometeor species from the environment...
!----------------------------------------------------------------

! Loop over all the hydrometeor species
do i_cond = 1, n_cond_species

  ! Count points where fall-in flux is nonzero, and store
  ! their indices
  nc(i_cond) = 0
  do ic = 1, n_points
    if ( flux_cond(ic,i_cond) > zero ) then
      nc(i_cond) = nc(i_cond) + 1
      index_ic(nc(i_cond),i_cond) = ic
    end if
  end do

  ! If any points
  if ( nc(i_cond) > 0 ) then

    ! Call routine to do fall-in calculations
    call fall_in( n_points, nc(i_cond), index_ic(:,i_cond),                    &
                  cond_params(i_cond)%pt % cp,                                 &
                  dt_over_rhod_lz, flux_cond(:,i_cond),                        &
                  fallin_wind_u, fallin_wind_v, fallin_wind_w,                 &
                  fallin_temp,                                                 &
                  cp_tot, rho_wet, wind_u, wind_v, wind_w,                     &
                  temperature, q_cond(:,i_cond),                               &
                  l_diags, i_cond, moist_proc_diags,                           &
                  n_points_diag, n_diags, diags_super )

  end if

end do


!----------------------------------------------------------------
! 3) Call microphysics scheme to compute explicit processes and
!    species properties needed for implicit phase-change solve
!----------------------------------------------------------------

! Find points where there is now non-zero mixing-ratio
! of each condensed water species
do i_cond = 1, n_cond_species
  nc(i_cond) = 0
  do ic = 1, n_points
    if ( q_cond(ic,i_cond) > zero ) then
      nc(i_cond) = nc(i_cond) + 1
      index_ic(nc(i_cond),i_cond) = ic
    end if
  end do
end do

! Calculate air dry-density
call calc_rho_dry( n_points, temperature, q_vap, pressure,                     &
                   rho_dry )

! Complete calculation of rho_wet (currently stores q_tot;
! dry-mass to wet-mass conversion factor = 1 + q_tot)
do ic = 1, n_points
  rho_wet(ic) = rho_dry(ic) * ( one + rho_wet(ic) )
end do

! Call microphysics routine:
! Back-Of-The-Envelope Microphysics Scheme (BOTEMS)
call microphysics_1( n_points, n_points_super, nc, index_ic,                   &
                     linear_qs(:,i_ref_temp),                                  &
                     linear_qs(:,i_qsat_liq_ref),                              &
                     linear_qs(:,i_dqsatdT_liq),                               &
                     delta_t, vert_len, rho_dry, rho_wet,                      &
                     cp_tot, temperature, q_vap, q_cond,                       &
                     q_loc_cond, wf_cond, kq_cond, kt_cond,                    &
                     dq_frz_cond, l_diags, moist_proc_diags,                   &
                     n_points_diag, n_diags, diags_super )

! Note: the above call to microphysics_1 needs to be done even
! if all the condensed water species are currently zero
! everywhere, as the test for activation of new condesation
! is done inside microphysics_1.

! If any condensed water species are non-zero anywhere
if ( maxval(nc) > 0 ) then

  !--------------------------------------------------------------
  ! 4) Implicitly solve condensation / evaporation and melting.
  !--------------------------------------------------------------

  ! This routine solves for the exchange of water vapour and heat
  ! between the parcel air and all the hydrometeor species at
  ! once, implicitly with respect to the parcel's temperature
  ! and water vapour mixing ratio (ie saturation):

  call phase_change_solve( n_points, n_points_super,                           &
             nc, index_ic, cmpr, k, call_string, linear_qs,                    &
             delta_z, delta_t, wind_w_excess, prev_temp,                       &
             kq_cond, kt_cond, wf_cond, dq_frz_cond,                           &
             q_loc_cond, q_cond, cp_tot, temperature, q_vap,                   &
             l_diags, moist_proc_diags, n_points_diag, n_diags, diags_super )


  ! do not add any more processes that alter T or q after this
  ! point; it will ruin the implicit solve!

  ! Now you maybe wondering, what happened to stage number 5)
  ! in this routine?  If you knew that, you would have a far
  ! more complete understanding of the nature of life, the
  ! universe and everything.


  !--------------------------------------------------------------
  ! 6) Microphysical processes occuring after phase-changes
  !--------------------------------------------------------------
  ! Currently only includes autoconversion of liquid cloud
  ! into rain
  call microphysics_2( n_points, n_points_super, nc, index_ic,                 &
                       delta_t, vert_len, wf_cond, q_cond,                     &
                       l_diags, moist_proc_diags,                              &
                       n_points_diag, n_diags, diags_super )


  !--------------------------------------------------------------
  ! 7) Fall-out of each hydrometeor species
  !--------------------------------------------------------------

  ! Loop over all present hydrometeor species
  do i_cond = 1, n_cond_species
    if ( nc(i_cond) > 0 ) then

      ! Compute fall-out
      call fall_out_flux( n_points, nc(i_cond),                                &
                          index_ic(:,i_cond),                                  &
                          wf_cond(:,i_cond), delta_t, vert_len,                &
                          dt_over_rhod_lz,                                     &
                          q_cond(:,i_cond), flux_cond(:,i_cond),               &
                          l_diags, i_cond, moist_proc_diags,                   &
                          n_points_diag, n_diags, diags_super )

      ! Tidy-up by removing any VERY small condensate mixing-ratios
      ! left-over.
      do ic2 = 1, nc(i_cond)
        ic = index_ic(ic2,i_cond)
        ! Some compilers allow values less than tiny to exist with
        ! reduced precision; the reduced precision can make some
        ! calculations unsafe, so reset to zero when this happens.
        ! Note: just removing water here, but shouldn't matter for
        ! conservation because the amounts removed will be less than
        ! epsilon * q_tot and therefore not representable increments.
        if ( q_cond(ic,i_cond) < min_float ) then
          q_cond(ic,i_cond) = zero
        end if
      end do

    end if
  end do


  !--------------------------------------------------------------
  ! 8) Save any diagnostics not already grabbed inside
  !    phase_change_solve or the microphysics calculations
  !--------------------------------------------------------------
  if ( l_diags ) then

    do i_cond = 1, n_cond_species
      if ( nc(i_cond) > 0 ) then
        ! Fall-speed diagnostic
        if ( moist_proc_diags % diags_cond(i_cond)%pt                          &
             % fall_speed % flag ) then
          ! Extract super-array address
          i_super = moist_proc_diags % diags_cond(i_cond)%pt                   &
                    % fall_speed % i_super
          ! Copy field into diags super-array
          do ic2 = 1, nc(i_cond)
            ic = index_ic(ic2,i_cond)
            diags_super(ic,i_super) = wf_cond(ic,i_cond)
          end do
        end if
      end if
    end do

  end if


end if  ! ( maxval(nc) > 0 )


if ( i_check_conservation > i_check_bad_none ) then
  ! Calculate conserved variables on exit from this routine
  ! and check whether everything adds up.
  call moist_proc_conservation( cmpr, k, call_string, n_points_super,          &
         dt_over_rhod_lz,                                                      &
         q_cond, q_vap, temperature, wind_u, wind_v, wind_w,                   &
         flux_cond, temperature,                                               &
         wind_u, wind_v, wind_w,                                               &
         cons_vars )
end if


! Check outputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  where_string = "End of moist_proc call for "                //               &
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
end if


return
end subroutine moist_proc


end module moist_proc_mod
