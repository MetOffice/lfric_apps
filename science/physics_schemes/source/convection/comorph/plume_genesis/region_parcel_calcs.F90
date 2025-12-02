! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module region_parcel_calcs_mod

implicit none

contains

! Subroutine to calculate the initiating parcel properties
! for either an updraft or downdraft initiating from
! one of the four sub-grid regions
! (dry/clear, liquid-cloud, mixed-phase cloud, ice/rain).
! Uses the parcel lifting code from the main convective ascent
! (parcel_dyn) to estimate the local moist static stability,
! which is used to parameterise the mass-flux initiating from
! the current model-level.
subroutine region_parcel_calcs( n_points, n_points_super,                      &
                                nc, index_ic, i_region,                        &
                                l_down, cmpr_init, k,                          &
                                grid_k, grid_next, grid_kpdk,                  &
                                fields_k,                                      &
                                frac_r_k, temperature_r_k,                     &
                                q_vap_r_k, q_cond_loc_k,                       &
                                virt_temp_k, virt_temp_kpdk,                   &
                                cloudfracs_k, layer_mass_k,                    &
                                par_radius,                                    &
                                turb_tl, turb_qt,                              &
                                delta_temp_neut, delta_qvap_neut,              &
                                nc2, index_ic2, init_mass,                     &
                                fields_par, pert_tl, pert_qt )

use comorph_constants_mod, only: real_cvprec, zero, one,                       &
                                 name_length, n_cond_species,                  &
                                 i_check_bad_values_cmpr, i_check_bad_none
use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_dealloc
use fields_type_mod, only: n_fields, i_wind_u, i_wind_w,                       &
                           i_temperature, i_q_vap, i_q_cl,                     &
                           i_qc_first, i_qc_last, field_names, field_positive
use grid_type_mod, only: n_grid, i_height, i_pressure
use subregion_mod, only: n_regions, region_names
use cloudfracs_type_mod, only: n_cloudfracs, i_frac_liq, i_frac_ice
use buoyancy_mod, only: i_buoy=>i_mean_buoy, i_prev

use set_region_cond_fields_mod, only: set_region_cond_fields
use set_cp_tot_mod, only: set_cp_tot
use sat_adjust_mod, only: sat_adjust
use parcel_dyn_mod, only: parcel_dyn, i_call_genesis
use calc_virt_temp_mod, only: calc_virt_temp
use calc_init_mass_mod, only: calc_init_mass
use calc_init_par_fields_mod, only: calc_init_par_fields
use check_bad_values_mod, only: check_bad_values_cmpr

implicit none

! Total number of points
integer, intent(in) :: n_points

! Size of input super-arrays (maybe larger than n_points,
! as arrays are dimensioned with the max size they could need
! on any level and reused for all levels)
integer, intent(in) :: n_points_super

! Number of points where the current region has nonzero fraction
! (and indices of those points)
integer, intent(in) :: nc
integer, intent(in) :: index_ic(nc)

! Index of the current sub-grid region (dry,liq,mph,icr)
integer, intent(in) :: i_region

! Flag for downdraft versus updrafts
logical, intent(in) :: l_down

! Structure storing i,j indices of all the points
! (for error reporting)
type(cmpr_type), intent(in) :: cmpr_init
! Current model-level index
integer, intent(in) :: k

! Height and pressure at level k, the next model-level interface,
! and the next full level
real(kind=real_cvprec), intent(in) :: grid_k                                   &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_next                                &
                                      ( n_points_super, n_grid )
real(kind=real_cvprec), intent(in) :: grid_kpdk                                &
                                      ( n_points_super, n_grid )

! Grid-mean primary fields at level k
real(kind=real_cvprec), intent(in) :: fields_k                                 &
                                    ( n_points_super, n_fields )

! Subregion properties at level k
real(kind=real_cvprec), intent(in) :: frac_r_k                                 &
                                      ( n_points, n_regions )
real(kind=real_cvprec), intent(in) :: temperature_r_k                          &
                                      ( n_points, n_regions )
real(kind=real_cvprec), intent(in) :: q_vap_r_k                                &
                                      ( n_points, n_regions )
real(kind=real_cvprec), intent(in) :: q_cond_loc_k                             &
                                      ( n_points, n_cond_species)

! Environment virtual temperatures from current and next
! full model-levels
real(kind=real_cvprec), intent(in) :: virt_temp_k(n_points)
real(kind=real_cvprec), intent(in) :: virt_temp_kpdk(n_points)

! Cloud and rain fractions from level k
real(kind=real_cvprec), intent(in) :: cloudfracs_k                             &
                                ( n_points_super, n_cloudfracs )

! Dry-mass on level k
real(kind=real_cvprec), intent(in) :: layer_mass_k(n_points)

! Parcel radius
real(kind=real_cvprec), intent(in) :: par_radius(n_points)

! Turbulence-based perturbations to Tl and qt
! at next model-level interface
real(kind=real_cvprec), intent(in) :: turb_tl(n_points)
real(kind=real_cvprec), intent(in) :: turb_qt(n_points)

! Neutrally buoyant perturbations to T,q so-as to yield
! a specified RH perturbation
real(kind=real_cvprec), intent(in) :: delta_temp_neut(n_points)
real(kind=real_cvprec), intent(in) :: delta_qvap_neut(n_points)

! Number of initiating points
integer, intent(out) :: nc2
! Compression indices for initiating points
integer, intent(out) :: index_ic2(n_points)

! Initiating mass from the current region
real(kind=real_cvprec), intent(out) :: init_mass(n_points)

! Unperturbed initiating parcel properties at level k
real(kind=real_cvprec), intent(out) :: fields_par                              &
                            ( n_points, i_temperature:n_fields )
! Perturbations applied to Tl and qt of the initiating parcel
real(kind=real_cvprec), intent(out) :: pert_tl(n_points)
real(kind=real_cvprec), intent(out) :: pert_qt(n_points)


! Primary fields from current region
! (used as work array by the parcel_dyn calls)
real(kind=real_cvprec) :: fields_par_next( nc, n_fields )

! Compressed copy of grid-mean fields from level k
real(kind=real_cvprec) :: fields_k_cmpr ( nc, i_wind_u:i_temperature )

! Grid fields compressed onto region points
real(kind=real_cvprec) :: grid_k_cmpr ( nc, n_grid )
real(kind=real_cvprec) :: grid_kpdk_cmpr ( nc, n_grid )

! Grid-mean virtual temperatures compressed onto region points
real(kind=real_cvprec) :: virt_temp_k_cmpr(nc)
real(kind=real_cvprec) :: virt_temp_next_cmpr(nc)
real(kind=real_cvprec) :: virt_temp_kpdk_cmpr(nc)

! Parcel virtual temperature before lifting
real(kind=real_cvprec) :: par_prev_virt_temp(nc)

! Parcel radius for parcel_dyn call
real(kind=real_cvprec) :: par_radius_cmpr(nc)

! Weights for interpolating onto the next model-level interface
real(kind=real_cvprec) :: interp(nc)

! Parcel total heat capacity for sat_adjust call
real(kind=real_cvprec) :: cp_tot(nc)

! Dummy for unused input to parcel_dyn
logical :: l_within_bl(nc)

! Parcel buoyancy super-array needed by parcel_dyn...
! Number of vars in the super-array:
! 1) height
! 2) buoyancy
integer, parameter :: n_buoy_vars = 2
! Max number of different heights at which to store these:
! 1) Previous level (starting level k)
! 2) Next model-level interface
! 3) Height where parcel first hit saturation
integer, parameter :: max_buoy_heights = 3
! Declare the super-array:
real(kind=real_cvprec) :: buoyancy_super                                       &
                     ( nc, n_buoy_vars, max_buoy_heights )

! Address of next level and saturation height in buoyancy_super
integer :: i_next(nc)
integer :: i_sat(nc)

! Points for further compression onto initiating points
integer :: ic2_first

! Structure storing i,j indices of points currently being worked
! on (for error reporting)
type(cmpr_type) :: cmpr
logical :: l_positive
character(len=name_length) :: field_name

! Strings used for error reporting
character(len=name_length) :: draft_string
character(len=name_length) :: call_string

! Flag passed into sat_adjust to make it update q_vap and q_cl
logical, parameter :: l_update_q_true = .true.

! Flag passed into parcel_dyn to indicate no tracer calculations
logical, parameter :: l_tracer_false = .false.

! Loop counters
integer :: ic, ic2, i_field

!character(len=*), parameter :: routinename                                    &
!                               = "REGION_PARCEL_CALCS"


! Set string indicating whether this is updraft or downdraft
if ( l_down ) then
  draft_string = "downdraft"
else
  draft_string = "updraft"
end if

! Initialise output initiating mass to zero at all points
do ic = 1, n_points
  init_mass(ic) = zero
end do

! Set vertical interpolation weight
do ic2 = 1, nc
  ic = index_ic(ic2)
  interp(ic2) = ( grid_next(ic,i_height) - grid_k(ic,i_height) )               &
              / ( grid_kpdk(ic,i_height) - grid_k(ic,i_height) )
end do


! Compress the primary field properties for the current
! region at level k into a fields super-array

! Copy the grid-mean winds
do i_field = i_wind_u, i_wind_w
  do ic2 = 1, nc
    fields_par_next(ic2,i_field) =                                             &
       fields_k(index_ic(ic2),i_field)
  end do
end do

! Copy temperature and water-vapour from the current region
do ic2 = 1, nc
  ic = index_ic(ic2)
  fields_par_next(ic2,i_temperature)=temperature_r_k(ic,i_region)
  fields_par_next(ic2,i_q_vap)      =q_vap_r_k(ic,i_region)
end do

! Set local condensed water species mixing ratios and
! cloud-fractions based on current region index:
call set_region_cond_fields( n_points, nc, index_ic, nc,                       &
                             i_region, q_cond_loc_k,                           &
                             cloudfracs_k(:,i_frac_ice),                       &
                             frac_r_k, fields_par_next )


! SETUP COMPRESSED ARRAYS FOR PARCEL_DYN call...

! Compress i,j indices used for error reporting
cmpr % n_points = nc
call cmpr_alloc( cmpr, nc )
do ic2 = 1, nc
  ic = index_ic(ic2)
  cmpr % index_i(ic2) = cmpr_init % index_i(ic)
  cmpr % index_j(ic2) = cmpr_init % index_j(ic)
end do

! Compress grid fields onto current region points
do i_field = 1, n_grid
  do ic2 = 1, nc
    ic = index_ic(ic2)
    grid_k_cmpr(ic2,i_field)    = grid_k(ic,i_field)
    grid_kpdk_cmpr(ic2,i_field) = grid_kpdk(ic,i_field)
  end do
end do

! Compress parcel radius
do ic2 = 1, nc
  par_radius_cmpr(ic2) = par_radius(index_ic(ic2))
end do

! Calculate environment winds and temperature seen by the
! rising test parcel; set them to values at k
do i_field = i_wind_u, i_temperature
  do ic2 = 1, nc
    fields_k_cmpr(ic2,i_field)                                                 &
      = fields_k(index_ic(ic2),i_field)
  end do
end do

! Compress environment virtual temperatures from level k
! and interpolate to next model-level interface as well
do ic2 = 1, nc
  ic = index_ic(ic2)
  virt_temp_k_cmpr(ic2)    = virt_temp_k(ic)
  virt_temp_kpdk_cmpr(ic2) = virt_temp_kpdk(ic)
  virt_temp_next_cmpr(ic2)                                                     &
    = (one-interp(ic2)) * virt_temp_k(ic)                                      &
    +      interp(ic2)  * virt_temp_kpdk(ic)
end do

! Copy current region properties at level k into the output
! parcel properties array.  We also use this as the previous
! level parcel properties input to parcel_dyn.
do i_field = i_temperature, n_fields
  do ic2 = 1, nc
    fields_par(ic2,i_field) = fields_par_next(ic2,i_field)
  end do
end do

! Ensure any liquid cloud in the parcel is in moist equilibrium
! (i.e. condense or evaporate q_cl to obtain liquid saturation).
! The liq and mph regions should already be liquid-saturated,
! but occasionally they can be supersaturated due to the
! input profiles being silly.
! parcel_dyn will adjust the parcel to saturation after lifting,
! so to get an accurate estimate of the moist adiabatic lapse
! rate from it, the parcel needs to be saturated beforehand too.
call set_cp_tot( nc, nc,                                                       &
                 fields_par_next(:,i_q_vap),                                   &
                 fields_par_next(:,i_qc_first:i_qc_last), cp_tot)
call sat_adjust( nc, l_update_q_true, grid_k_cmpr(:,i_pressure),               &
                 fields_par_next(:,i_temperature),                             &
                 fields_par_next(:,i_q_vap),                                   &
                 fields_par_next(:,i_q_cl), cp_tot )

! Calculate parcel virtual temperature before the lifting
call calc_virt_temp( nc, nc,                                                   &
                     fields_par_next(:,i_temperature),                         &
                     fields_par_next(:,i_q_vap),                               &
                     fields_par_next(:,i_qc_first:i_qc_last),                  &
                     par_prev_virt_temp )

! Initialise buoyancy super-array
do ic2 = 1, nc
  i_next(ic2) = 2
  i_sat(ic2) = 0
  buoyancy_super(ic2,i_height,i_prev)                                          &
    = grid_k_cmpr(ic2,i_height)
  buoyancy_super(ic2,i_height,i_next(ic2))                                     &
    = grid_kpdk_cmpr(ic2,i_height)
  ! Calculate virtual temperature excess at level k
  ! (should be zero by design, but not precisely zero
  !  due to error of linearising Tv in calc_env_partitions)
  buoyancy_super(ic2,i_buoy,i_prev)                                            &
    = par_prev_virt_temp(ic2) - virt_temp_k_cmpr(ic2)
end do


! Call parcel_dyn to calculate the parcel properties after
! lifting or subsiding from level k to the next
! model-level interface (half a level-step).
call_string = trim(adjustl(draft_string)) // "_genesis_" //                    &
              trim(adjustl(region_names(i_region)))
call parcel_dyn( nc, n_points, nc, nc, nc, nc, 1, 1, n_fields,                 &
                 max_buoy_heights, n_buoy_vars, l_down, l_tracer_false,        &
                 cmpr, k, call_string, i_call_genesis,                         &
                 grid_k_cmpr, grid_kpdk_cmpr,                                  &
                 l_within_bl, init_mass, init_mass,                            &
                 init_mass, par_radius_cmpr,                                   &
                 fields_k_cmpr,                                                &
                 virt_temp_k_cmpr, virt_temp_kpdk_cmpr,                        &
                 fields_par(:,i_temperature:i_qc_last),                        &
                 fields_par_next,                                              &
                 buoyancy_super, i_next, i_sat )


! Calculate intiating mass-source based on the static stability N^2
call calc_init_mass( n_points, nc, index_ic,                                   &
                     i_next, i_sat, n_buoy_vars, max_buoy_heights,             &
                     buoyancy_super, layer_mass_k, frac_r_k(:,i_region),       &
                     grid_next(:,i_height), virt_temp_next_cmpr,               &
                     init_mass )


! Recompress onto only points where convection is initiating...
ic2_first = 0
over_nc: do ic2 = 1, nc
  ! First see if any points in the list not initiating
  if ( .not. init_mass(ic2) > zero ) then
    ic2_first = ic2
    exit over_nc
  end if
end do over_nc
if ( ic2_first > 0 ) then
  ! If at least one point not initiating,
  ! store indices of those that are initiating.
  nc2 = 0
  do ic2 = 1, nc
    if ( init_mass(ic2) > zero ) then
      nc2 = nc2 + 1
      index_ic2(nc2) = ic2
    end if
  end do
else
  ! If all points are initiating, no need to recompress
  nc2 = nc
end if

! If any initiating points
if ( nc2 > 0 ) then

  ! If not all points in this region initiated...
  if ( nc2 < nc ) then
    ! Recompress all required work arrays in-situ.
    ! Note: only need to loop from the first non-initiating
    ! point, ic2_first.
    cmpr % n_points = nc2
    do ic2 = ic2_first, nc2
      cmpr % index_i(ic2) = cmpr % index_i(index_ic2(ic2))
      cmpr % index_j(ic2) = cmpr % index_j(index_ic2(ic2))
    end do
    do i_field = i_temperature, n_fields
      do ic2 = ic2_first, nc2
        fields_par(ic2,i_field)                                                &
          = fields_par(index_ic2(ic2),i_field)
      end do
    end do
    do i_field = 1, n_grid
      do ic2 = ic2_first, nc2
        grid_k_cmpr(ic2,i_field)                                               &
          = grid_k_cmpr(index_ic2(ic2),i_field)
      end do
    end do
    do ic2 = ic2_first, nc2
      init_mass(ic2) = init_mass(index_ic2(ic2))
    end do
    ! Convert the indices index_ic2 to reference the input
    ! fields (n_points)
    do ic2 = 1, nc2
      index_ic2(ic2) = index_ic(index_ic2(ic2))
    end do
  else  ! ( nc2 < nc )
    ! All region points still included; still need to
    ! copy the compression indices
    do ic2 = 1, nc2
      index_ic2(ic2)  = index_ic(ic2)
    end do
  end if  ! ( nc2 < nc )


  ! Calculate initiating parcel properties and perturbations
  call calc_init_par_fields( n_points, nc2, index_ic2, l_down, i_region,       &
                             grid_k_cmpr(:,i_pressure),                        &
                             cloudfracs_k(:,i_frac_liq),                       &
                             turb_tl, turb_qt,                                 &
                             delta_temp_neut, delta_qvap_neut,                 &
                             fields_par, pert_tl, pert_qt )


  if ( i_check_bad_values_cmpr > i_check_bad_none ) then
    ! Check outputs for bad values (NaN, Inf, etc).

    call_string = "On output from region_parcel_calcs, region: " //            &
                   trim(adjustl(region_names(i_region))) // " " //             &
                   trim(adjustl(draft_string))
    l_positive = .true.
    field_name = "massflux_d"
    call check_bad_values_cmpr( cmpr, k, init_mass,                            &
                                call_string, field_name, l_positive )
    do i_field = i_temperature, n_fields
      call check_bad_values_cmpr( cmpr, k, fields_par(:,i_field),              &
                                  call_string, field_names(i_field),           &
                                  field_positive(i_field) )
    end do
    l_positive = .false.
    field_name = "pert_tl"
    call check_bad_values_cmpr( cmpr, k, pert_tl,                              &
                                call_string, field_name, l_positive )
    field_name = "pert_qt"
    call check_bad_values_cmpr( cmpr, k, pert_qt,                              &
                                call_string, field_name, l_positive )

  end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )

end if  ! ( nc2 > 0 )


call cmpr_dealloc( cmpr )


return
end subroutine region_parcel_calcs


end module region_parcel_calcs_mod
