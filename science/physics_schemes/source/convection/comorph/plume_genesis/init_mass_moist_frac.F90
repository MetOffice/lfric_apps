! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module init_mass_moist_frac_mod

implicit none

contains

! Subroutine to calculate initiating dry-mass per unit surface
! area and initiating parcel properties on a given level.
! This version decomposes the grid area into 4 separate regions:
! - clear-sky         (dry)
! - liquid-only cloud (liq)
! - mixed-phase cloud (mph)
! - ice and rain only (icr)
! The differing temperature, humidity and condensed water
! contents of these 4 regions are calculated in
! calc_env_regions.
!
! With the properties of each region defined, we then do a
! 1-level moist parcel ascent/descent from each of the
! regions (dry,liq,mph,icr) to calculate a separate moist static
! stability Nsq for each region.  The mass-flux emanating
! from any region with negative Nsq is then calculated as:
! M_init = sqrt(-Nsq)
!
! The total mass-flux emanating from this level is then the
! area-weighted sum of mass-fluxes from the 4 regions.

subroutine init_mass_moist_frac( n_points, n_points_super,                     &
                                 l_tracer, n_fields_tot,                       &
                                 cmpr_init, k, l_within_bl,                    &
                                 layer_mass_k, turb_len_k, par_radius_amp,     &
                                 turb_kmh, turb_kph,                           &
                                 grid_km1, grid_kmh, grid_k,                   &
                                 grid_kph, grid_kp1,                           &
                                 fields_km1, fields_k,                         &
                                 fields_kp1, cloudfracs_k,                     &
                                 virt_temp_km1, virt_temp_k,                   &
                                 virt_temp_kp1,                                &
                                 updraft_par_gen, dndraft_par_gen,             &
                                 genesis_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, name_length,                     &
                                 zero, one, two,                               &
                                 n_updraft_types, n_dndraft_types,             &
                                 n_cond_species, i_cond_cl,                    &
                                 l_homog_conv_bl, min_delta, sqrt_min_float,   &
                                 max_ent_frac_up, max_ent_frac_dn,             &
                                 par_gen_core_fac, l_turb_par_gen,             &
                                 i_cfl_local, i_cfl_local_none,                &
                                 i_check_bad_values_cmpr, i_check_bad_none

use grid_type_mod, only: n_grid, i_pressure
use parcel_type_mod, only: parcel_type, i_massflux_d, i_radius,                &
                           parcel_check_bad_values
use fields_type_mod, only: n_fields, i_wind_u, i_temperature,                  &
                           i_q_vap, i_qc_first, i_qc_last, i_q_cl,             &
                           field_names, field_positive
use turb_type_mod, only: n_turb
use cloudfracs_type_mod, only: n_cloudfracs
use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_dealloc
use subregion_mod, only: n_regions, i_liq, i_mph, region_names
use genesis_diags_type_mod, only: genesis_diags_type
use check_bad_values_mod, only: check_bad_values_cmpr

use calc_env_regions_mod, only: calc_env_regions
use set_cp_tot_mod, only: set_cp_tot
use lat_heat_mod, only: lat_heat_incr, i_phase_change_evp,                     &
                        set_l_con
use set_qsat_mod, only: set_qsat_liq
use set_dqsatdt_mod, only: set_dqsatdt_liq
use calc_turb_parcel_mod, only: calc_turb_parcel
use set_par_fields_mod, only: set_par_fields
use region_parcel_calcs_mod, only: region_parcel_calcs
use cor_init_mass_liq_1_mod, only: cor_init_mass_liq_1
use cfl_limit_init_mass_mod, only: cfl_limit_init_mass
use add_region_parcel_mod, only: add_region_parcel
use normalise_init_parcel_mod, only: normalise_init_parcel

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Flag for whether tracers are being processed in this call
logical, intent(in) :: l_tracer

! Total number of fields in the input super-arrays,
! including tracers if they are used
integer, intent(in) :: n_fields_tot

! Structure storing indices of points where convective
! initiation might happen
type(cmpr_type), intent(in) :: cmpr_init

! Current model-level index (for error reporting)
integer, intent(in) :: k

! Flag for whether the current level is below the boundary-layer
! top at each point
logical, intent(in) :: l_within_bl(n_points)

! Dry-mass per unit surface area on the current model-level
real(kind=real_cvprec), intent(in) :: layer_mass_k(n_points)

! Turbulence length-scale at full-level k
real(kind=real_cvprec), intent(in) :: turb_len_k(n_points)

! parcel radius amplification factor
real(kind=real_cvprec), intent(in) :: par_radius_amp(n_points)

! Super-array containing turbulence fields at the
! model-level interfaces
real(kind=real_cvprec), intent(in) :: turb_kmh                                 &
                                      ( n_points_super, n_turb )
real(kind=real_cvprec), intent(in) :: turb_kph                                 &
                                      ( n_points_super, n_turb )

! Height and pressure at k, model-level interfaces, and
! neighbouring full levels
real(kind=real_cvprec), intent(in) :: grid_km1                                 &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kmh                                 &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_k                                   &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kph                                 &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kp1                                 &
                                      ( n_points_super, n_grid )

! Primary model-fields at level k and neighbouring full levels
real(kind=real_cvprec), intent(in) :: fields_km1                               &
                                ( n_points_super, n_fields_tot )
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                ( n_points_super, n_fields_tot )
real(kind=real_cvprec), intent(in) :: fields_kp1                               &
                                ( n_points_super, n_fields_tot )

! Compressed super-arrays containing the 3 cloud-fraction fields
! cf_liq, cf_ice, cf_bulk, and the rain-fraction, at level k
real(kind=real_cvprec), intent(in) :: cloudfracs_k                             &
                                ( n_points_super, n_cloudfracs )

! Environment virtual temperature at level k and
! neighbouring full levels
real(kind=real_cvprec), intent(in) :: virt_temp_km1(n_points)
real(kind=real_cvprec), intent(in) :: virt_temp_k(n_points)
real(kind=real_cvprec), intent(in) :: virt_temp_kp1(n_points)

! Structures containing initiating mass-source properties
! for updrafts and downdrafts initiating from level k
type(parcel_type), intent(in out) :: updraft_par_gen(n_updraft_types)
type(parcel_type), intent(in out) :: dndraft_par_gen(n_dndraft_types)

! Structure storing diagnostics and assocaited meta-data
type(genesis_diags_type), intent(in out) :: genesis_diags

! Super-array storing diagnostics to be output
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                          ( n_points_super, genesis_diags % n_diags )


! Arrays for properties of the 4 sub-grid regions
! Fractional area of each sub-region
real(kind=real_cvprec) :: frac_r_k ( n_points, n_regions )
! Temperature in each sub-region
real(kind=real_cvprec) :: temperature_r_k ( n_points, n_regions )
! Water vapour mixing ratio in each sub-region
real(kind=real_cvprec) :: q_vap_r_k ( n_points, n_regions )
! Local condensate mixing ratios within the regions in which
! each species is allowed to be non-zero
real(kind=real_cvprec) :: q_cond_loc_k(n_points,n_cond_species)

! Neutrally buoyant perturbations to T,q needed for a fractional increase of RH
real(kind=real_cvprec) :: delta_temp_neut(n_points)
real(kind=real_cvprec) :: delta_qvap_neut(n_points)

! Tl and qt perturbations for current region, for updrafts and downdrafts
! (can be set separately for each convection type)
real(kind=real_cvprec) :: pert_tl_up_t ( n_points, n_updraft_types )
real(kind=real_cvprec) :: pert_qt_up_t ( n_points, n_updraft_types )
real(kind=real_cvprec) :: pert_tl_dn_t ( n_points, n_dndraft_types )
real(kind=real_cvprec) :: pert_qt_dn_t ( n_points, n_dndraft_types )

! Turbulence-based perturbations to u,v,w,Tl,qt at level k
real(kind=real_cvprec) :: turb_pert_k                                          &
                          ( n_points, i_wind_u:i_q_vap )

! Initiating updraft and downdraft mass-flux for each convection type
real(kind=real_cvprec) :: init_mass_up_t ( n_points, n_updraft_types )
real(kind=real_cvprec) :: init_mass_dn_t ( n_points, n_dndraft_types )

! Unperturbed initiating parcel properties from current region,
! for updrafts and downdrafts
! (doesn't include winds, since they are assumed equal in all
!  the sub-grid regions)
real(kind=real_cvprec) :: fields_par_up                                        &
                          ( n_points, i_temperature:n_fields )
real(kind=real_cvprec) :: fields_par_dn                                        &
                          ( n_points, i_temperature:n_fields )

! Env k total heat capacity
real(kind=real_cvprec) :: cp_tot(n_points)
! Env k latent heat of condensation
real(kind=real_cvprec) :: L_con(n_points)
! Level k liquid water temperature
real(kind=real_cvprec) :: tl_k(n_points)
! qsat at level k liquid water temperature
real(kind=real_cvprec) :: qsat_tl_k(n_points)
! dqsat/dT w.r.t. liquid water at level k
real(kind=real_cvprec) :: dqsatdt(n_points)   ! Value at Tl
real(kind=real_cvprec) :: dqsatdt_t(n_points) ! Value at T

! Non-turbulent RH perturbation applied for each convection type
real(kind=real_cvprec) :: updraft_rhpert_t ( n_points, n_updraft_types )
real(kind=real_cvprec) :: dndraft_rhpert_t ( n_points, n_dndraft_types )
! Triggering area fraction of each convection type in each sub-grid region
real(kind=real_cvprec) :: updraft_frac_r_t                                     &
                          ( n_points, n_regions, n_updraft_types )
real(kind=real_cvprec) :: dndraft_frac_r_t                                     &
                          ( n_points, n_regions, n_dndraft_types )

! Rescaling of mass-sources so implicit w.r.t. cloud fraction
real(kind=real_cvprec) :: imp_coef(n_points)

! Work arrays for bad value checking
real(kind=real_cvprec) :: work1(n_points)
real(kind=real_cvprec) :: work2(n_points)

! Flag for updraft vs downdraft calls
logical :: l_down

! Flags for whether updrafts and downdrafts are active on this model-level
logical :: l_updraft
logical :: l_dndraft

! Points where current region has nonzero fraction
integer :: nc
integer :: index_ic(n_points)

! Points for recompression onto points not in the boundary-layer
integer :: ic2_first
integer :: nc2
integer :: index_ic2(n_points)

! Points where updrafts and downdrafts initiate in current region
integer :: nc_up
integer :: index_ic_up(n_points)
integer :: nc_dn
integer :: index_ic_dn(n_points)

! Character string for error messages
character(len=name_length) :: call_string
type(cmpr_type) :: cmpr_check
logical :: l_positive
character(len=name_length) :: field_name

! Loop counters
integer :: ic, ic2, i_region, i_field, i_type

!CHARACTER(LEN=*), PARAMETER :: routinename                                    &
!                               = "INIT_MASS_MOIST_FRAC"


if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  ! Check input fields for bad values
  do i_field = 1, n_fields_tot
    call_string = "Start of init_mass_moist_frac, fields_km1"
    call check_bad_values_cmpr( cmpr_init, k, fields_km1(:,i_field),           &
                                call_string, field_names(i_field),             &
                                field_positive(i_field) )
    call_string = "Start of init_mass_moist_frac, fields_k"
    call check_bad_values_cmpr( cmpr_init, k, fields_k(:,i_field),             &
                                call_string, field_names(i_field),             &
                                field_positive(i_field) )
    call_string = "Start of init_mass_moist_frac, fields_kp1"
    call check_bad_values_cmpr( cmpr_init, k, fields_kp1(:,i_field),           &
                                call_string, field_names(i_field),             &
                                field_positive(i_field) )
  end do
end if

! Set flags for whether updrafts and downdrafts are allowed from this
! model-level.  Firstly, they are only alowed if updrafts / downdrafts
! are actually switched on (i.e. n_up/dndraft_types is not set to zero).
! Then, in conv_genesis_ctl, updrafts are disabled at the model-top
! (and downdrafts disabled at the model-bottom) by setting the number
! of points in the compression list to zero; account for this as well.
l_updraft = .false.
if ( n_updraft_types > 0 ) then
  if ( updraft_par_gen(1) % cmpr % n_points > 0 )  l_updraft = .true.
end if
l_dndraft = .false.
if ( n_dndraft_types > 0 ) then
  if ( dndraft_par_gen(1) % cmpr % n_points > 0 )  l_dndraft = .true.
end if

! Calculate turbulence-based parcel perturbations
call calc_turb_parcel( n_points, n_points_super, cmpr_init, k,                 &
                       grid_km1, grid_kmh, grid_k,                             &
                       grid_kph, grid_kp1,                                     &
                       fields_km1, fields_k, fields_kp1,                       &
                       turb_kmh, turb_kph, turb_pert_k )

! Calculate properties of the sub-grid cloudy, rainy, icy and
! clear-sky regions of the grid-box at level k
call_string = "genesis_k"
call calc_env_regions(                                                         &
         n_points, n_points_super, cmpr_init, k, call_string,                  &
         grid_k(:,i_pressure), cloudfracs_k,                                   &
         fields_k(:,i_temperature), fields_k(:,i_q_vap),                       &
         fields_k(:,i_qc_first:i_qc_last),                                     &
         frac_r_k, temperature_r_k, q_vap_r_k,                                 &
         q_cond_loc_k,                                                         &
         genesis_diags, diags_super,                                           &
         delta_temp_neut, delta_qvap_neut )

! Compute dqsat/dT at grid-mean T
call set_qsat_liq( n_points, fields_k(:,i_temperature), grid_k(:,i_pressure),  &
                   work1 )
call set_dqsatdt_liq( n_points, fields_k(:,i_temperature), work1,              &
                      dqsatdt_t )

! Calculate liquid water temperature Tl from level k, and qsat(Tl(k))
do ic = 1, n_points
  tl_k(ic) = fields_k(ic,i_temperature)
end do
call set_cp_tot( n_points, n_points_super, fields_k(:,i_q_vap),                &
                 fields_k(:,i_qc_first:i_qc_last),                             &
                 cp_tot )
call lat_heat_incr( n_points, n_points, i_phase_change_evp,                    &
                    cp_tot, tl_k,                                              &
                    dq=fields_k(:,i_q_cl) )
call set_qsat_liq( n_points, tl_k, grid_k(:,i_pressure), qsat_tl_k )
call set_dqsatdt_liq( n_points, tl_k, qsat_tl_k, dqsatdt )
call set_l_con( n_points, tl_k, L_con )

! Set parcel initial radius, winds and tracers
! (these are assumed equal in all sub-grid regions, so can set them now).
if ( l_updraft ) then
  ! Updrafts
  l_down = .false.
  do i_type = 1, n_updraft_types
    call set_par_fields( n_points, n_points_super, n_fields_tot,               &
                         cmpr_init, k, l_tracer, l_down, i_type,               &
                         grid_k, fields_k, frac_r_k,                           &
                         par_radius_amp, turb_pert_k, turb_len_k,              &
                         updraft_par_gen(i_type) % par_super,                  &
                         updraft_par_gen(i_type) % mean_super,                 &
                         updraft_par_gen(i_type) % core_super,                 &
                         updraft_rhpert_t(:,i_type),                           &
                         updraft_frac_r_t(:,:,i_type) )
  end do
end if
if ( l_dndraft ) then
  ! Downdrafts
  l_down = .true.
  do i_type = 1, n_dndraft_types
    call set_par_fields( n_points, n_points_super, n_fields_tot,               &
                         cmpr_init, k, l_tracer, l_down, i_type,               &
                         grid_k, fields_k, frac_r_k,                           &
                         par_radius_amp, turb_pert_k, turb_len_k,              &
                         dndraft_par_gen(i_type) % par_super,                  &
                         dndraft_par_gen(i_type) % mean_super,                 &
                         dndraft_par_gen(i_type) % core_super,                 &
                         dndraft_rhpert_t(:,i_type),                           &
                         dndraft_frac_r_t(:,:,i_type) )
  end do
end if


! For each of the sub-grid regions i_dry, i_liq, i_mph, i_icr...
do i_region = 1, n_regions


  ! Find points where the current region has nonzero fraction
  nc = 0
  do ic = 1, n_points
    if ( frac_r_k(ic,i_region) > min_delta ) then
      nc = nc + 1
      index_ic(nc) = ic
    end if
  end do

  ! If any points in the current region...
  if ( nc > 0 ) then

    ! Initialise number of updraft and downdraft initiating
    ! points to zero
    nc_up = 0
    nc_dn = 0

    ! If updrafts active
    if ( l_updraft ) then

      ! Calculate initiating mass-flux and parcel properties
      ! for updrafts initiating in the current region.
      ! Only the parcel properties which differ between the
      ! regions (T,q,qc,cf) are calculated here.  The parcel
      ! winds and tracers were calculated early and are assumed
      ! equal in all regions.
      l_down = .false.
      call region_parcel_calcs( n_points, n_points_super, n_updraft_types,     &
                                nc, index_ic, i_region,                        &
                                l_down, cmpr_init, k,                          &
                                grid_k, grid_kph, grid_kp1, fields_k,          &
                                frac_r_k, temperature_r_k,                     &
                                q_vap_r_k, q_cond_loc_k,                       &
                                virt_temp_k, virt_temp_kp1,                    &
                                cloudfracs_k, layer_mass_k,                    &
                                updraft_par_gen(1) % par_super(:,i_radius),    &
                                updraft_frac_r_t, updraft_rhpert_t,            &
                                turb_pert_k(:,i_temperature),                  &
                                turb_pert_k(:,i_q_vap),                        &
                                delta_temp_neut, delta_qvap_neut,              &
                                nc_up, index_ic_up, init_mass_up_t,            &
                                fields_par_up, pert_tl_up_t, pert_qt_up_t,     &
                                genesis_diags, diags_super )


    end if  ! ( l_updraft )

    ! If downdrafts active
    if ( l_dndraft ) then

      if ( l_homog_conv_bl ) then
        ! If using homogenisation of convective source terms
        ! within the boundary-layer, suppress downdraft
        ! initiation mass-sources within the BL,
        ! since there is no way they can get out of the BL!
        ! This should make no difference to the run,
        ! but avoids redundant downdraft calculations.

        ! Recompress onto only points NOT within the
        ! boundary-layer...

        ic2_first = 0
        over_nc: do ic2 = 1, nc
          ! First see if any points in the list ARE within the BL
          if ( l_within_bl(index_ic(ic2)) ) then
            ic2_first = ic2
            exit over_nc
          end if
        end do over_nc
        if ( ic2_first > 0 ) then
          ! If at least one point is within the BL,
          ! store indices of those that are not
          nc2 = 0
          do ic2 = 1, nc
            if ( .not. l_within_bl(index_ic(ic2)) ) then
              nc2 = nc2 + 1
              index_ic2(nc2) = ic2
            end if
          end do
        else
          ! If all points are initiating, no need to recompress
          nc2 = nc
        end if

        ! If any points NOT within the BL...
        if ( nc2 > 0 ) then
          ! If some points in this region are within the BL...
          if ( nc2 < nc ) then
            ! Recompress the work compression indices in-situ.
            ! Note: only need to loop from the first non-included
            ! point, ic2_first.
            do ic2 = ic2_first, nc2
              index_ic(ic2) = index_ic(index_ic2(ic2))
            end do
          end if
        end if

        ! Reset number of points to the new more compressed value
        nc = nc2

      end if  ! ( l_homog_conv_bl )

      ! If still any points in the list
      if ( nc > 0 ) then

        ! Calculate initiating mass-flux and parcel properties
        ! for downdrafts initiating in the current region.
        l_down = .true.
        call region_parcel_calcs( n_points, n_points_super, n_dndraft_types,   &
                                nc, index_ic, i_region,                        &
                                l_down, cmpr_init, k,                          &
                                grid_k, grid_kmh, grid_km1, fields_k,          &
                                frac_r_k, temperature_r_k,                     &
                                q_vap_r_k, q_cond_loc_k,                       &
                                virt_temp_k, virt_temp_km1,                    &
                                cloudfracs_k, layer_mass_k,                    &
                                dndraft_par_gen(1) % par_super(:,i_radius),    &
                                dndraft_frac_r_t, dndraft_rhpert_t,            &
                                turb_pert_k(:,i_temperature),                  &
                                turb_pert_k(:,i_q_vap),                        &
                                delta_temp_neut, delta_qvap_neut,              &
                                nc_dn, index_ic_dn, init_mass_dn_t,            &
                                fields_par_dn, pert_tl_dn_t, pert_qt_dn_t,     &
                                genesis_diags, diags_super )

      end if

    end if  ! ( l_dndraft )


    ! If initiation mass-sources have occured from a
    ! liquid cloud region
    if ( ( i_region == i_liq .or. i_region == i_mph ) .and.                    &
         ( nc_up > 0 .or. nc_dn > 0 ) ) then

      ! Apply implicit correction to the initiating masses to account for
      ! reduction of liquid-cloud fraction by the increments from initiation.
      call cor_init_mass_liq_1( n_points, n_points_super,                      &
                                nc_up, index_ic_up, nc_dn, index_ic_dn,        &
                                frac_r_k(:,i_region), layer_mass_k,            &
                                q_cond_loc_k(:,i_cond_cl),                     &
                                cp_tot, L_con, dqsatdt, dqsatdt_t,             &
                                fields_par_up, fields_par_dn,                  &
                                pert_tl_up_t, pert_qt_up_t,                    &
                                pert_tl_dn_t, pert_qt_dn_t,                    &
                                grid_km1, grid_k, grid_kp1,                    &
                                fields_km1, fields_kp1,                        &
                                init_mass_up_t, init_mass_dn_t,                &
                                imp_coef )

    end if  ! ( i_region == i_liq .OR. i_region == i_mph )

    if ( i_cfl_local > i_cfl_local_none ) then
      ! Impose CFL limit on the initiating masses from the current region...
      if ( nc_up > 0 ) then
        call cfl_limit_init_mass( n_points, nc_up, index_ic_up,                &
                                  n_updraft_types,                             &
                                  max_ent_frac_up, l_within_bl,                &
                                  layer_mass_k, frac_r_k(:,i_region),          &
                                  init_mass_up_t )
      end if
      if ( nc_dn > 0 ) then
        call cfl_limit_init_mass( n_points, nc_dn, index_ic_dn,                &
                                  n_dndraft_types,                             &
                                  max_ent_frac_dn, l_within_bl,                &
                                  layer_mass_k, frac_r_k(:,i_region),          &
                                  init_mass_dn_t )
      end if
    end if

    if ( i_check_bad_values_cmpr > i_check_bad_none ) then
      ! Check region mass-flux and primary fields for bad values

      if ( nc_up > 0 ) then
        call_string = "init_mass_moist_frac, region: " //                      &
                      trim(adjustl(region_names(i_region))) // " updraft"
        call cmpr_alloc( cmpr_check, nc_up )
        cmpr_check % n_points = nc_up
        do ic2 = 1, nc_up
          ic = index_ic_up(ic2)
          cmpr_check % index_i(ic2) = cmpr_init % index_i(ic)
          cmpr_check % index_j(ic2) = cmpr_init % index_j(ic)
        end do
        do i_type = 1, n_updraft_types
          write(field_name,"(A,I3)") "massflux_d type", i_type
          call check_bad_values_cmpr( cmpr_check, k, init_mass_up_t(:,i_type), &
                                      call_string, field_name, l_positive )
        end do
        do i_field = i_temperature, n_fields
          call check_bad_values_cmpr( cmpr_check, k, fields_par_up(:,i_field), &
                                      call_string, field_names(i_field),       &
                                      field_positive(i_field) )
        end do
        do i_type = 1, n_updraft_types
          do ic2 = 1, nc_up
            work1(ic2) = fields_par_up(ic2,i_temperature)                      &
                       + par_gen_core_fac * pert_tl_up_t(ic2,i_type)
            work2(ic2) = fields_par_up(ic2,i_q_vap)                            &
                       + par_gen_core_fac * pert_qt_up_t(ic2,i_type)
          end do
          write(field_name,"(A,I3)") "temperature + tl_pert conv type ", i_type
          call check_bad_values_cmpr( cmpr_check, k, work1,                    &
                                      call_string, field_name, l_positive )
          write(field_name,"(A,I3)") "q_vap + qt_pert conv type ", i_type
          call check_bad_values_cmpr( cmpr_check, k, work2,                    &
                                      call_string, field_name, l_positive )
        end do
        call cmpr_dealloc( cmpr_check )
      end if  ! ( nc_up > 0 )

      if ( nc_dn > 0 ) then
        call_string = "init_mass_moist_frac, region: " //                      &
                      trim(adjustl(region_names(i_region))) // " dndraft"
        call cmpr_alloc( cmpr_check, nc_dn )
        cmpr_check % n_points = nc_dn
        do ic2 = 1, nc_dn
          ic = index_ic_dn(ic2)
          cmpr_check % index_i(ic2) = cmpr_init % index_i(ic)
          cmpr_check % index_j(ic2) = cmpr_init % index_j(ic)
        end do
        do i_type = 1, n_dndraft_types
          write(field_name,"(A,I3)") "massflux_d type", i_type
          call check_bad_values_cmpr( cmpr_check, k, init_mass_dn_t(:,i_type), &
                                      call_string, field_name, l_positive )
        end do
        do i_field = i_temperature, n_fields
          call check_bad_values_cmpr( cmpr_check, k, fields_par_dn(:,i_field), &
                                      call_string, field_names(i_field),       &
                                      field_positive(i_field) )
        end do
        do i_type = 1, n_dndraft_types
          do ic2 = 1, nc_dn
            work1(ic2) = fields_par_dn(ic2,i_temperature)                      &
                       + par_gen_core_fac * pert_tl_dn_t(ic2,i_type)
            work2(ic2) = fields_par_dn(ic2,i_q_vap)                            &
                       + par_gen_core_fac * pert_qt_dn_t(ic2,i_type)
          end do
          write(field_name,"(A,I3)") "temperature + tl_pert conv type ", i_type
          call check_bad_values_cmpr( cmpr_check, k, work1,                    &
                                      call_string, field_name, l_positive )
          write(field_name,"(A,I3)") "q_vap + qt_pert conv type ", i_type
          call check_bad_values_cmpr( cmpr_check, k, work2,                    &
                                      call_string, field_name, l_positive )
        end do
        call cmpr_dealloc( cmpr_check )
      end if  ! ( nc_dn > 0 )

    end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )

    ! If any updrafts have initiated in the current region
    if ( nc_up > 0 ) then
      do i_type = 1, n_updraft_types
        ! Add mass-weighted contribution to updraft initiating
        ! mass and parcel properties
        call add_region_parcel( n_points, nc_up, index_ic_up,                  &
                                init_mass_up_t(:,i_type), fields_par_up,       &
                                pert_tl_up_t(:,i_type), pert_qt_up_t(:,i_type),&
                                updraft_par_gen(i_type) % par_super,           &
                                updraft_par_gen(i_type) % mean_super,          &
                                updraft_par_gen(i_type) % core_super )
      end do
    end if

    ! If any downdrafts have initiated in the current region
    if ( nc_dn > 0 ) then
      do i_type = 1, n_dndraft_types
        ! Add mass-weighted contribution to downdraft initiating
        ! mass and parcel properties
        call add_region_parcel( n_points, nc_dn, index_ic_dn,                  &
                                init_mass_dn_t(:,i_type), fields_par_dn,       &
                                pert_tl_dn_t(:,i_type), pert_qt_dn_t(:,i_type),&
                                dndraft_par_gen(i_type) % par_super,           &
                                dndraft_par_gen(i_type) % mean_super,          &
                                dndraft_par_gen(i_type) % core_super )
      end do
    end if

  end if  ! ( nc > 0 )


end do  ! i_region = 1, n_regions
! End of loop over sub-grid regions


! Normalise the mass-weighted means, and do some miscellaneous checks...

if ( l_updraft ) then
  do i_type = 1, n_updraft_types
    ! Find points where the updraft initiating mass-flux is non-zero
    nc = 0
    do ic = 1, n_points
      if ( updraft_par_gen(i_type) % par_super(ic,i_massflux_d) > zero ) then
        nc = nc + 1
        index_ic(nc) = ic
      end if
    end do
    ! If any points have updraft initiating mass-sources,
    ! normalise the mass-weighted means
    if ( nc > 0 ) then
      l_down = .false.
      call normalise_init_parcel( n_points, nc, index_ic,                      &
                                  fields_k(:,i_q_vap),                         &
                                  updraft_par_gen(i_type) % par_super,         &
                                  updraft_par_gen(i_type) % mean_super,        &
                                  updraft_par_gen(i_type) % core_super )
    end if
  end do
end if

if ( l_dndraft ) then
  do i_type = 1, n_dndraft_types
    ! Find points where the downdraft initiating mass-flux is non-zero
    nc = 0
    do ic = 1, n_points
      if ( dndraft_par_gen(i_type) % par_super(ic,i_massflux_d) > zero ) then
        nc = nc + 1
        index_ic(nc) = ic
      end if
    end do
    ! If any points have downdraft initiating mass-sources,
    ! normalise the mass-weighted means
    if ( nc > 0 ) then
      l_down = .true.
      call normalise_init_parcel( n_points, nc, index_ic,                      &
                                  fields_k(:,i_q_vap),                         &
                                  dndraft_par_gen(i_type) % par_super,         &
                                  dndraft_par_gen(i_type) % mean_super,        &
                                  dndraft_par_gen(i_type) % core_super )
    end if
  end do
end if

if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  ! Check output initiating parcels for bad values (NaN, Inf, etc)
  if ( l_updraft ) then
    call_string = "End of init_mass_moist_frac; updraft_par_gen"
    do i_type = 1, n_updraft_types
      call parcel_check_bad_values( updraft_par_gen(i_type), n_fields_tot,     &
                                    k, call_string )
    end do
  end if
  if ( l_dndraft ) then
    call_string = "End of init_mass_moist_frac; dndraft_par_gen"
    do i_type = 1, n_dndraft_types
      call parcel_check_bad_values( dndraft_par_gen(i_type), n_fields_tot,     &
                                    k, call_string )
    end do
  end if
end if


return
end subroutine init_mass_moist_frac


end module init_mass_moist_frac_mod
