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
! calc_env_partitions.
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
                                 cmpr_init, k,                                 &
                                 l_within_bl, layer_mass_k, turb_len_k,        &
                                 par_radius_amp, turb_kmh, turb_kph,           &
                                 grid_km1, grid_kmh, grid_k,                   &
                                 grid_kph, grid_kp1,                           &
                                 fields_km1, fields_k,                         &
                                 fields_kp1, cloudfracs_k,                     &
                                 virt_temp_km1, virt_temp_k,                   &
                                 virt_temp_kp1,                                &
                                 updraft_par_gen,                              &
                                 dndraft_par_gen )

use comorph_constants_mod, only: real_cvprec, name_length,                     &
                     zero, one, two, comorph_timestep,                         &
                     n_cond_species, i_cond_cl,                                &
                     l_homog_conv_bl, min_delta, sqrt_min_float,               &
                     max_ent_frac_up, max_ent_frac_dn, par_gen_core_fac,       &
                     i_check_bad_values_cmpr, i_check_bad_none

use grid_type_mod, only: n_grid, i_pressure
use parcel_type_mod, only: parcel_type, i_massflux_d,                          &
                           parcel_check_bad_values
use fields_type_mod, only: n_fields, i_wind_u, i_temperature,                  &
                           i_q_vap, i_qc_first, i_qc_last, i_q_cl,             &
                           field_names, field_positive
use turb_type_mod, only: n_turb
use cloudfracs_type_mod, only: n_cloudfracs
use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_dealloc
use subregion_mod, only: n_regions, i_liq, i_mph, region_names
use check_bad_values_mod, only: check_bad_values_cmpr

use calc_env_partitions_mod, only: calc_env_partitions
use set_cp_tot_mod, only: set_cp_tot
use lat_heat_mod, only: lat_heat_incr, i_phase_change_evp,                     &
                        set_l_con
use set_qsat_mod, only: set_qsat_liq
use set_dqsatdt_mod, only: set_dqsatdt_liq
use calc_turb_parcel_mod, only: calc_turb_parcel
use set_par_winds_mod, only: set_par_winds
use region_parcel_calcs_mod, only: region_parcel_calcs
use calc_qss_forcing_init_mod, only: calc_qss_forcing_init
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
type(parcel_type), intent(in out) :: updraft_par_gen
type(parcel_type), intent(in out) :: dndraft_par_gen


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

! Neutrally buoyant perturbations to T,q so-as to yield
! a specified RH perturbation
real(kind=real_cvprec) :: delta_temp_neut(n_points)
real(kind=real_cvprec) :: delta_qvap_neut(n_points)

! Tl and qt perturbations for current region, for updrafts and downdrafts
real(kind=real_cvprec) :: pert_tl_up(n_points)
real(kind=real_cvprec) :: pert_qt_up(n_points)
real(kind=real_cvprec) :: pert_tl_dn(n_points)
real(kind=real_cvprec) :: pert_qt_dn(n_points)

! Turbulence-based perturbations to u,v,w,Tl,qt at level k
real(kind=real_cvprec) :: turb_pert_k                                          &
                          ( n_points, i_wind_u:i_q_vap )
! Parcel initial radius at level k
real(kind=real_cvprec) :: par_radius_k(n_points)

! Updraft and downdraft mass-flux initiating from current region
real(kind=real_cvprec) :: init_mass_up(n_points)
real(kind=real_cvprec) :: init_mass_dn(n_points)

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
! Level k q_vap + q_cl
real(kind=real_cvprec) :: qt_k(n_points)
! qsat at level k liquid water temperature
real(kind=real_cvprec) :: qsat_tl_k(n_points)
! dqsat/dT w.r.t. liquid water at level k
real(kind=real_cvprec) :: dqsatdt(n_points)

! Variables used in implicit rescaling of mass-sources
! from liquid-cloud:
! Sum of updraft and downdraft initiation mass-sources
real(kind=real_cvprec) :: mass_forc(n_points)

! Sum of initiation forcing of qss = q_vap+q_cl - qsat(Tl)
real(kind=real_cvprec) :: qss_forc(n_points)

! Rescaling of mass-sources so implicit w.r.t. cloud fraction
real(kind=real_cvprec) :: imp_coef(n_points)

! Work arrays for bad value checking
real(kind=real_cvprec) :: work1(n_points)
real(kind=real_cvprec) :: work2(n_points)

! Flag for updraft vs downdraft calls
logical :: l_down

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
integer :: ic, ic2, i_region, i_field

!character(len=*), parameter :: routinename                                    &
!                               = "INIT_MASS_MOIST_FRAC"


! Calculate properties of the sub-grid cloudy, rainy, icy and
! clear-sky regions of the grid-box at level k
call_string = "genesis_k"
call calc_env_partitions(                                                      &
         n_points, n_points_super, cmpr_init, k, call_string,                  &
         grid_k(:,i_pressure), cloudfracs_k,                                   &
         fields_k(:,i_temperature), fields_k(:,i_q_vap),                       &
         fields_k(:,i_qc_first:i_qc_last),                                     &
         frac_r_k, temperature_r_k, q_vap_r_k,                                 &
         q_cond_loc_k,                                                         &
         delta_temp_neut, delta_qvap_neut )

! Calculate liquid water temperature Tl from level k, and qsat(Tl(k))
do ic = 1, n_points
  tl_k(ic) = fields_k(ic,i_temperature)
  qt_k(ic) = fields_k(ic,i_q_vap) + fields_k(ic,i_q_cl)
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


! Calculate turbulence-based parcel perturbations and radius
call calc_turb_parcel( n_points, n_points_super, cmpr_init, k,                 &
                       grid_km1, grid_kmh, grid_k,                             &
                       grid_kph, grid_kp1,                                     &
                       fields_km1, fields_k, fields_kp1,                       &
                       turb_len_k, par_radius_amp, turb_kmh, turb_kph,         &
                       turb_pert_k, par_radius_k )


! Set parcel initial winds and tracers
! (these are assumed equal in all sub-grid regions, so can set them now).
if ( updraft_par_gen % cmpr % n_points > 0 ) then
  l_down = .false.
  ! Updrafts
  call set_par_winds( n_points, n_points_super, n_fields_tot,                  &
                      cmpr_init, k, l_tracer, l_down,                          &
                      fields_k, turb_pert_k, par_radius_k,                     &
                      updraft_par_gen % par_super,                             &
                      updraft_par_gen % mean_super,                            &
                      updraft_par_gen % core_super )
end if
if ( dndraft_par_gen % cmpr % n_points > 0 ) then
  l_down = .true.
  ! Downdrafts
  call set_par_winds( n_points, n_points_super, n_fields_tot,                  &
                      cmpr_init, k, l_tracer, l_down,                          &
                      fields_k, turb_pert_k, par_radius_k,                     &
                      dndraft_par_gen % par_super,                             &
                      dndraft_par_gen % mean_super,                            &
                      dndraft_par_gen % core_super )
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
    if ( updraft_par_gen % cmpr % n_points > 0 ) then

      ! Calculate initiating mass-flux and parcel properties
      ! for updrafts initiating in the current region.
      ! Only the parcel properties which differ between the
      ! regions (T,q,qc,cf) are calculated here.  The parcel
      ! winds and tracers were calculated early and are assumed
      ! equal in all regions.
      l_down = .false.
      call region_parcel_calcs( n_points, n_points_super,                      &
                                nc, index_ic, i_region,                        &
                                l_down, cmpr_init, k,                          &
                                grid_k, grid_kph, grid_kp1,                    &
                                fields_k,                                      &
                                frac_r_k, temperature_r_k,                     &
                                q_vap_r_k, q_cond_loc_k,                       &
                                virt_temp_k, virt_temp_kp1,                    &
                                cloudfracs_k, layer_mass_k,                    &
                                par_radius_k,                                  &
                                turb_pert_k(:,i_temperature),                  &
                                turb_pert_k(:,i_q_vap),                        &
                                delta_temp_neut, delta_qvap_neut,              &
                                nc_up, index_ic_up, init_mass_up,              &
                                fields_par_up, pert_tl_up, pert_qt_up )


    end if  ! ( updraft_par_gen % cmpr % n_points > 0 )

    ! If downdrafts active
    if ( dndraft_par_gen % cmpr % n_points > 0 ) then

      if ( l_homog_conv_bl ) then
        ! If using homogenisation of convective source terms
        ! within the boundary-layer, suppress downdraft
        ! initiation mass-sources within the BL,
        ! since there is no way they can get out of the BL!
        ! This should make no difference to the run,
        ! but avoids redundant downdraft calculations.

        ! Recompress onto only points not within the
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

        ! If any points not within the BL...
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
        call region_parcel_calcs( n_points, n_points_super,                    &
                                nc, index_ic, i_region,                        &
                                l_down, cmpr_init, k,                          &
                                grid_k, grid_kmh, grid_km1,                    &
                                fields_k,                                      &
                                frac_r_k, temperature_r_k,                     &
                                q_vap_r_k, q_cond_loc_k,                       &
                                virt_temp_k, virt_temp_km1,                    &
                                cloudfracs_k, layer_mass_k,                    &
                                par_radius_k,                                  &
                                turb_pert_k(:,i_temperature),                  &
                                turb_pert_k(:,i_q_vap),                        &
                                delta_temp_neut, delta_qvap_neut,              &
                                nc_dn, index_ic_dn, init_mass_dn,              &
                                fields_par_dn, pert_tl_dn, pert_qt_dn )

      end if

    end if  ! ( dndraft_par_gen % cmpr % n_points > 0 )


    ! If initiation mass-sources have occured from a
    ! liquid cloud region
    if ( ( i_region == i_liq .or. i_region == i_mph ) .and.                    &
         ( nc_up > 0 .or. nc_dn > 0 ) ) then

      ! Modify the initiation mass-sources consistent with
      ! an implicit time discretisation, accounting for the
      ! reduction in liquid cloud fraction over the timestep
      ! due to the convective heating / drying...

      ! M_up_imp = M_up_exp f_imp/f_exp
      ! M_dn_imp = M_dn_exp f_imp/f_exp
      !
      ! When initiating from large-scale liquid cloud
      ! (regions liq and mph), the main feedback limiting
      ! the mass-flux is that the convective heating and drying
      ! will cause the cloud-scheme to remove the cloud.
      ! Exactly how quickly will depend on the details of the
      ! cloud-scheme being used.  Here we estimate the convective
      ! forcing of supersaturation qss = qw - qsat
      ! by the initiating mass-source, and use this to
      ! parameterise the expected rate of removal of liquid cloud
      ! by the convection...
      ! Convection also removes liquid cloud directly due to
      ! the cloud-volume which is extracted to form the
      ! initiating parcel.  The implicitness rescaling also
      ! accounts for this effect.

      ! Initialise sums over updraft and downdraft sources
      do ic = 1, n_points
        ! Sum of updraft and downdraft initiating mass-sources
        ! from the current region:
        mass_forc(ic) = zero
        ! Sum of forcing of qss by the initiating mass-sources
        ! from the current region:
        qss_forc(ic) = zero
      end do

      ! Sum the contributions to the forcing of mass and
      ! qss = qw - q_sat
      ! from both updrafts and downdrafts initiating in the
      ! current region / model-level
      if ( nc_up > 0 ) then
        do ic2 = 1, nc_up
          ic = index_ic_up(ic2)
          mass_forc(ic) = mass_forc(ic) + init_mass_up(ic2)
        end do
        call calc_qss_forcing_init( n_points, n_points_super,                  &
                                    nc_up, index_ic_up,                        &
                                    init_mass_up, fields_par_up,               &
                                    pert_tl_up, pert_qt_up,                    &
                                    grid_k, grid_kp1,                          &
                                    fields_k, fields_kp1,                      &
                                    qss_forc )
      end if
      if ( nc_dn > 0 ) then
        do ic2 = 1, nc_dn
          ic = index_ic_dn(ic2)
          mass_forc(ic) = mass_forc(ic) + init_mass_dn(ic2)
        end do
        call calc_qss_forcing_init( n_points, n_points_super,                  &
                                    nc_dn, index_ic_dn,                        &
                                    init_mass_dn, fields_par_dn,               &
                                    pert_tl_dn, pert_qt_dn,                    &
                                    grid_k, grid_km1,                          &
                                    fields_k, fields_km1,                      &
                                    qss_forc )
      end if

      ! Apply limit to use explicit solution in case where
      ! the initiation mass-sources actually act to increase
      ! the liquid cloud instead of decreasing it
      ! (i.e. only retain negative values of dqss/dt)
      do ic = 1, n_points
        qss_forc(ic) = min( qss_forc(ic), zero )
      end do

      ! Note: technically qss_forc now needs to be normalised
      ! by dividing it by layer_mass.  This is done in the
      ! formula for imp_coef below...

      ! Calculate the implicit rescaling factor
      ! imp_coef = f_imp/f_exp.
      ! Parameterise the forcing of liquid cloud-fraction as:
      !
      ! (df/dt)_exp = -M_init / M_layr
      !             + f_exp * dqss/dt
      !                     / ( 2 (1 + dqs/dT Lc/cp) qc_incloud )
      !
      ! (1st term from direct fractional cloud-volume removal,
      !  2nd term from homogeneous forcing of remaining cloud)
      !
      ! The implicit liquid cloud fraction is given by:
      !
      ! f_imp = f_exp + (df/dt)_exp * (f_imp/f_exp) * dt
      !
      ! => (f_imp/f_exp) ( 1 - ( (df/dt)_exp / f_exp ) * dt ) = 1
      ! => (f_imp/f_exp)
      !       = 1 / ( 1 - ( (df/dt)_exp / f_exp ) * dt )
      !
      ! Substituting our above expression for df/dt, we get:
      do ic = 1, n_points
        imp_coef(ic) = one / ( one - (                                         &
          ! Term from volume-removal:
            -mass_forc(ic) / max(frac_r_k(ic,i_region), sqrt_min_float)        &
          +                                                                    &
          ! Term from homogeneous forcing:
            qss_forc(ic)                                                       &
            / ( two * (one + dqsatdt(ic)*L_con(ic)/cp_tot(ic))                 &
                    * max( q_cond_loc_k(ic,i_cond_cl), sqrt_min_float ) )      &
                                      )                                        &
          ! All terms have scaling of dt/M_layr
            * comorph_timestep / layer_mass_k(ic) )
      end do

      ! Modify the initiating mass-fluxes accordingly
      do ic2 = 1, nc_up
        ic = index_ic_up(ic2)
        init_mass_up(ic2) = init_mass_up(ic2) * imp_coef(ic)
      end do
      do ic2 = 1, nc_dn
        ic = index_ic_dn(ic2)
        init_mass_dn(ic2) = init_mass_dn(ic2) * imp_coef(ic)
      end do

    end if  ! ( i_region == i_liq .or. i_region == i_mph )

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
        field_name = "massflux_d"
        l_positive = .true.
        call check_bad_values_cmpr( cmpr_check, k, init_mass_up,               &
                                    call_string, field_name, l_positive )
        do i_field = i_temperature, n_fields
          call check_bad_values_cmpr( cmpr_check, k, fields_par_up(:,i_field), &
                                      call_string, field_names(i_field),       &
                                      field_positive(i_field) )
        end do
        do ic2 = 1, nc_up
          work1(ic2) = fields_par_up(ic2,i_temperature)                        &
                     + par_gen_core_fac * pert_tl_up(ic2)
          work2(ic2) = fields_par_up(ic2,i_q_vap)                              &
                     + par_gen_core_fac * pert_qt_up(ic2)
        end do
        field_name = "temperature + tl_pert"
        call check_bad_values_cmpr( cmpr_check, k, work1,                      &
                                    call_string, field_name, l_positive )
        field_name = "q_vap + qt_pert"
        call check_bad_values_cmpr( cmpr_check, k, work2,                      &
                                    call_string, field_name, l_positive )
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
        field_name = "massflux_d"
        l_positive = .true.
        call check_bad_values_cmpr( cmpr_check, k, init_mass_dn,               &
                                    call_string, field_name, l_positive )
        do i_field = i_temperature, n_fields
          call check_bad_values_cmpr( cmpr_check, k, fields_par_dn(:,i_field), &
                                      call_string, field_names(i_field),       &
                                      field_positive(i_field) )
        end do
        do ic2 = 1, nc_dn
          work1(ic2) = fields_par_dn(ic2,i_temperature)                        &
                     + par_gen_core_fac * pert_tl_dn(ic2)
          work2(ic2) = fields_par_dn(ic2,i_q_vap)                              &
                     + par_gen_core_fac * pert_qt_dn(ic2)
        end do
        field_name = "temperature + tl_pert"
        call check_bad_values_cmpr( cmpr_check, k, work1,                      &
                                    call_string, field_name, l_positive )
        field_name = "q_vap + qt_pert"
        call check_bad_values_cmpr( cmpr_check, k, work2,                      &
                                    call_string, field_name, l_positive )
        call cmpr_dealloc( cmpr_check )
      end if  ! ( nc_dn > 0 )

    end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )

    ! If any updrafts have initiated in the current region
    if ( nc_up > 0 ) then
      ! Add mass-weighted contribution to updraft initiating
      ! mass and parcel properties
      call add_region_parcel( n_points, nc_up, index_ic_up,                    &
                              init_mass_up, fields_par_up,                     &
                              pert_tl_up, pert_qt_up,                          &
                              max_ent_frac_up, l_within_bl,                    &
                              layer_mass_k, frac_r_k(:,i_region),              &
                              updraft_par_gen % par_super                      &
                                                (:,i_massflux_d),              &
                              updraft_par_gen % mean_super,                    &
                              updraft_par_gen % core_super )
    end if

    ! If any downdrafts have initiated in the current region
    if ( nc_dn > 0 ) then
      ! Add mass-weighted contribution to downdraft initiating
      ! mass and parcel properties
      call add_region_parcel( n_points, nc_dn, index_ic_dn,                    &
                              init_mass_dn, fields_par_dn,                     &
                              pert_tl_dn, pert_qt_dn,                          &
                              max_ent_frac_dn, l_within_bl,                    &
                              layer_mass_k, frac_r_k(:,i_region),              &
                              dndraft_par_gen % par_super                      &
                                                (:,i_massflux_d),              &
                              dndraft_par_gen % mean_super,                    &
                              dndraft_par_gen % core_super )
    end if

  end if  ! ( nc > 0 )


end do  ! i_region = 1, n_regions
! End of loop over sub-grid regions


if ( updraft_par_gen % cmpr % n_points > 0 ) then
  ! Find points where the total updraft initiating mass-flux
  ! is non-zero
  nc = 0
  do ic = 1, n_points
    if ( updraft_par_gen % par_super(ic,i_massflux_d)                          &
         > zero ) then
      nc = nc + 1
      index_ic(nc) = ic
    end if
  end do
  ! If any points have updraft initiating mass-sources,
  ! normalise the mass-weighted means
  if ( nc > 0 ) then
    call normalise_init_parcel( n_points, nc, index_ic,                        &
                                fields_k(:,i_q_vap),                           &
                                updraft_par_gen % par_super                    &
                                                (:,i_massflux_d),              &
                                updraft_par_gen % mean_super,                  &
                                updraft_par_gen % core_super )
  end if
end if

if ( dndraft_par_gen % cmpr % n_points > 0 ) then
  ! Find points where the total downdraft initiating mass-flux
  ! is non-zero
  nc = 0
  do ic = 1, n_points
    if ( dndraft_par_gen % par_super(ic,i_massflux_d)                          &
         > zero ) then
      nc = nc + 1
      index_ic(nc) = ic
    end if
  end do
  ! If any points have downdraft initiating mass-sources,
  ! normalise the mass-weighted means
  if ( nc > 0 ) then
    call normalise_init_parcel( n_points, nc, index_ic,                        &
                                fields_k(:,i_q_vap),                           &
                                dndraft_par_gen % par_super                    &
                                                (:,i_massflux_d),              &
                                dndraft_par_gen % mean_super,                  &
                                dndraft_par_gen % core_super )
  end if
end if

if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  ! Check output initiating parcels for bad values (NaN, Inf, etc)
  if ( updraft_par_gen % cmpr % n_points > 0 ) then
    call_string = "End of init_mass_moist_frac; updraft_par_gen"
    call parcel_check_bad_values( updraft_par_gen, n_fields_tot,               &
                                  k, call_string )
  end if
  if ( dndraft_par_gen % cmpr % n_points > 0 ) then
    call_string = "End of init_mass_moist_frac; dndraft_par_gen"
    call parcel_check_bad_values( dndraft_par_gen, n_fields_tot,               &
                                  k, call_string )
  end if
end if


return
end subroutine init_mass_moist_frac


end module init_mass_moist_frac_mod
