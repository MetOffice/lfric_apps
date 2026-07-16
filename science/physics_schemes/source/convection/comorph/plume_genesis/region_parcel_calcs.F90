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
subroutine region_parcel_calcs( n_points, n_points_super, n_conv_types,        &
                                nc, index_ic, i_region,                        &
                                l_down, cmpr_init, k,                          &
                                grid_k, grid_next, grid_kpdk, fields_k,        &
                                frac_r_k, temperature_r_k,                     &
                                q_vap_r_k, q_cond_loc_k,                       &
                                virt_temp_k, virt_temp_kpdk,                   &
                                cloudfracs_k, layer_mass_k,                    &
                                par_radius,                                    &
                                frac_r_t, rhpert_t,                            &
                                turb_tl, turb_qt,                              &
                                delta_temp_neut, delta_qvap_neut,              &
                                nc2, index_ic2, init_mass_t,                   &
                                fields_par, pert_tl_t, pert_qt_t,              &
                                genesis_diags, diags_super )

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
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev,                     &
                       j_height, j_mean_buoy, j_env_tv
use genesis_diags_type_mod, only: genesis_diags_type

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

! Number of convection types
integer, intent(in) :: n_conv_types

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
! full model-levels, at start-of-timestep
real(kind=real_cvprec), intent(in) :: virt_temp_k(n_points)
real(kind=real_cvprec), intent(in) :: virt_temp_kpdk(n_points)

! Cloud and rain fractions from level k
real(kind=real_cvprec), intent(in) :: cloudfracs_k                             &
                                ( n_points_super, n_cloudfracs )

! Dry-mass on level k
real(kind=real_cvprec), intent(in) :: layer_mass_k(n_points)

! Parcel radius
real(kind=real_cvprec), intent(in) :: par_radius(n_points)

! Triggering area fraction of each convection type in each sub-grid region
real(kind=real_cvprec), intent(in) :: frac_r_t                                 &
                                      ( n_points, n_regions, n_conv_types )
! Non-turbulent RH perturbations for each convection type
real(kind=real_cvprec), intent(in) :: rhpert_t                                 &
                                      ( n_points, n_conv_types )

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

! Initiating mass from the current region, for each convection type
real(kind=real_cvprec), intent(out) :: init_mass_t(n_points,n_conv_types)

! Unperturbed initiating parcel properties at level k
real(kind=real_cvprec), intent(out) :: fields_par                              &
                            ( n_points, i_temperature:n_fields )
! Perturbations applied to Tl and qt of the initiating parcel
real(kind=real_cvprec), intent(out) :: pert_tl_t ( n_points, n_conv_types )
real(kind=real_cvprec), intent(out) :: pert_qt_t ( n_points, n_conv_types )

! Structure storing diagnostics and associated meta-data
type(genesis_diags_type), intent(in out) :: genesis_diags
! Super-array storing diagnostics to be output
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                          ( n_points_super, genesis_diags % n_diags )

! Primary fields from current region
! (used as work array by the parcel_dyn calls)
real(kind=real_cvprec) :: fields_par_next( nc, n_fields )

! Compressed copy of grid-mean fields from level k
real(kind=real_cvprec) :: fields_k_cmpr ( nc, n_fields )

! Grid fields compressed onto region points
real(kind=real_cvprec) :: grid_k_cmpr ( nc, n_grid )
real(kind=real_cvprec) :: grid_kpdk_cmpr ( nc, n_grid )

! Grid-mean virtual temperature at next, compressed onto region points
real(kind=real_cvprec) :: virt_temp_next_cmpr(nc)

! Parcel virtual temperature before lifting
real(kind=real_cvprec) :: par_prev_virt_temp(nc)

! Parcel radius for parcel_dyn call
real(kind=real_cvprec) :: par_radius_cmpr(nc)

! Local static stability from test ascent, for diagnostic
real(kind=real_cvprec) :: Nsq_cmpr(nc)

! Weights for interpolating onto the next model-level interface
real(kind=real_cvprec) :: interp(nc)

! Parcel total heat capacity for sat_adjust call
real(kind=real_cvprec) :: cp_tot(nc)

! Dummy for unused arguments to parcel_dyn
logical :: l_within_bl(nc)
real(kind=real_cvprec) :: dummy_zeros(nc)
! Need separate arrays for these as not intent(in) to parcel_dyn
real(kind=real_cvprec) :: par_w_drag(nc)
real(kind=real_cvprec) :: prev_ss(nc)
real(kind=real_cvprec) :: next_ss(nc)
real(kind=real_cvprec) :: prev_tvl(nc)
real(kind=real_cvprec) :: next_tvl(nc)

! Sub-levels super-array containing buoyancies output by parcel_dyn
real(kind=real_cvprec) :: sublevs ( nc, n_sublev_vars, max_sublevs )

! Address of next level and saturation height in sublevs
integer :: i_next(nc)
integer :: i_sat(nc)

! Points for further compression onto initiating points
integer :: ic2_first

! Flag for points where convection initiating
logical :: l_init(nc)

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
integer :: ic, ic2, i_field, i_type, i_super

!CHARACTER(LEN=*), PARAMETER :: routinename                                    &
!                               = "REGION_PARCEL_CALCS"


! Set string indicating whether this is updraft or downdraft
if ( l_down ) then
  draft_string = "downdraft"
else
  draft_string = "updraft"
end if

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


! SETUP COMPRESSED ARRAYS FOR PARCEL_DYN CALL...

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
  ic = index_ic(ic2)
  par_radius_cmpr(ic2) = par_radius(ic)
end do

! Calculate environment winds and temperature seen by the
! rising test parcel; set them to values at k.
! Note: other primary fields (water mixing-ratios, cloud-fractions, etc)
! are passed into parcel_dyn but are not used in the call below,
! so no need to set them in the compressed array.
! We also pass this array in place of the compressed cloudfracs super-array,
! since that too is not used in the call below.
do i_field = i_wind_u, i_temperature
  do ic2 = 1, nc
    fields_k_cmpr(ic2,i_field)                                                 &
      = fields_k(index_ic(ic2),i_field)
  end do
end do

do ic2 = 1, nc
  ic = index_ic(ic2)
  ! Interpolate environment virtual temperature onto next model-level interface
  virt_temp_next_cmpr(ic2) = (one-interp(ic2)) * virt_temp_k(ic)               &
    +                             interp(ic2)  * virt_temp_kpdk(ic)
end do

! Copy current region properties at level k into the output
! parcel properties array.  We also use this as the previous
! level parcel properties input to parcel_dyn.
do i_field = i_temperature, n_fields
  do ic2 = 1, nc
    fields_par(ic2,i_field) = fields_par_next(ic2,i_field)
  end do
end do

! Set dummy array of zeros for unused input arguments
do ic2 = 1, nc
  dummy_zeros(ic2) = zero
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
call sat_adjust( nc, l_update_q_true,                                          &
                 grid_k_cmpr(:,i_pressure),                                    &
                 fields_par_next(:,i_temperature),                             &
                 fields_par_next(:,i_q_vap),                                   &
                 fields_par_next(:,i_q_cl), cp_tot )

! Calculate parcel virtual temperature before the lifting
call calc_virt_temp( nc, nc,                                                   &
                     fields_par_next(:,i_temperature),                         &
                     fields_par_next(:,i_q_vap),                               &
                     fields_par_next(:,i_qc_first:i_qc_last),                  &
                     par_prev_virt_temp )

! Initialise sub-levels super-array used by parcel_dyn
do i_field = 1, n_sublev_vars
  do ic2 = 1, nc
    sublevs(ic2,i_field,1) = zero
    sublevs(ic2,i_field,2) = zero
  end do
end do
do ic2 = 1, nc
  ic = index_ic(ic2)
  i_next(ic2) = 2
  i_sat(ic2) = 0
  sublevs(ic2,j_height,i_prev)      = grid_k_cmpr(ic2,i_height)
  sublevs(ic2,j_height,i_next(ic2)) = grid_kpdk_cmpr(ic2,i_height)
  sublevs(ic2,j_env_tv,i_prev)      = virt_temp_k(ic)
  sublevs(ic2,j_env_tv,i_next(ic2)) = virt_temp_kpdk(ic)
  ! Calculate virtual temperature excess at level k
  ! (not precisely zero if comparing latest with start-of-timestep env)
  sublevs(ic2,j_mean_buoy,i_prev) = par_prev_virt_temp(ic2) - virt_temp_k(ic)
end do


! Call parcel_dyn to calculate the parcel properties after
! lifting or subsiding from level k to the next
! model-level interface (half a level-step).
call_string = trim(adjustl(draft_string)) // "_genesis_" //                    &
              trim(adjustl(region_names(i_region)))
call parcel_dyn( nc, n_points, nc, nc, nc, nc, 1, 1, n_fields,                 &
                 l_down, l_tracer_false,                                       &
                 cmpr, k, call_string, i_call_genesis,                         &
                 grid_k_cmpr, grid_kpdk_cmpr,                                  &
                 l_within_bl, dummy_zeros, dummy_zeros,                        &
                 dummy_zeros, par_radius_cmpr,                                 &
                 fields_k_cmpr, dummy_zeros,                                   &
                 fields_par, fields_par_next,                                  &
                 sublevs, i_next, i_sat, par_w_drag,                           &
                 prev_ss, next_ss, prev_tvl, next_tvl )


! Calculate intiating mass-source based on the static stability N^2
call calc_init_mass( n_points, nc, index_ic, n_conv_types, i_region,           &
                     i_next, i_sat, sublevs,                                   &
                     layer_mass_k, grid_next(:,i_height), virt_temp_next_cmpr, &
                     frac_r_t, init_mass_t, Nsq_cmpr )

! Nsq_cmpr now stores the calculated moist static stability from the
! test parcel ascent or descent.
! Save for diagnostic output if requested...
if ( genesis_diags % subregion_diags(i_region) % n_diags > 0 ) then
  i_super = 0
  if ( l_down ) then
    if ( genesis_diags % subregion_diags(i_region) % Nsq_dn % flag ) then
      i_super = genesis_diags % subregion_diags(i_region) % Nsq_dn % i_super
    end if
  else
    if ( genesis_diags % subregion_diags(i_region) % Nsq_up % flag ) then
      i_super = genesis_diags % subregion_diags(i_region) % Nsq_up % i_super
    end if
  end if
  if ( i_super > 0 ) then
    do ic2 = 1, nc
      diags_super(index_ic(ic2),i_super) = Nsq_cmpr(ic2)
    end do
  end if
end if

! Find points where at least one convection type has nonzero initiating mass
do ic2 = 1, nc
  l_init(ic2) = .false.
  do i_type = 1, n_conv_types
    if ( init_mass_t(ic2,i_type) > zero )  l_init(ic2) = .true.
  end do
end do

! Recompress onto only points where convection is initiating...
ic2_first = 0
over_nc: do ic2 = 1, nc
  ! First see if any points in the list NOT initiating
  if ( .not. l_init(ic2) ) then
    ic2_first = ic2
    exit over_nc
  end if
end do over_nc
if ( ic2_first > 0 ) then
  ! If at least one point not initiating,
  ! store indices of those that are initiating.
  nc2 = 0
  do ic2 = 1, nc
    if ( l_init(ic2) ) then
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
      do i_type = 1, n_conv_types
        init_mass_t(ic2,i_type) = init_mass_t(index_ic2(ic2),i_type)
      end do
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


  do i_type = 1, n_conv_types
    ! Set initiating parcel perturbations for each convection type
    call calc_init_par_fields( n_points, nc2, index_ic2, l_down, i_region,     &
                               grid_k_cmpr(:,i_pressure),                      &
                               cloudfracs_k(:,i_frac_liq),                     &
                               turb_tl, turb_qt,                               &
                               delta_temp_neut, delta_qvap_neut,               &
                               rhpert_t(:,i_type), fields_par,                 &
                               pert_tl_t(:,i_type), pert_qt_t(:,i_type) )
  end do

  if ( i_check_bad_values_cmpr > i_check_bad_none ) then
    ! Check outputs for bad values (NaN, Inf, etc).

    call_string = "On output from region_parcel_calcs, region: " //            &
                   trim(adjustl(region_names(i_region))) // " " //             &
                   trim(adjustl(draft_string))
    l_positive = .true.
    do i_type = 1, n_conv_types
      write(field_name,"(A,I3)") "massflux_d conv type ", i_type
      call check_bad_values_cmpr( cmpr, k, init_mass_t(:,i_type),              &
                                  call_string, field_name, l_positive )
    end do
    do i_field = i_temperature, n_fields
      call check_bad_values_cmpr( cmpr, k, fields_par(:,i_field),              &
                                  call_string, field_names(i_field),           &
                                  field_positive(i_field) )
    end do
    l_positive = .false.
    do i_type = 1, n_conv_types
      write(field_name,"(A,I3)") "pert_tl conv type ", i_type
      call check_bad_values_cmpr( cmpr, k, pert_tl_t(:,i_type),                &
                                  call_string, field_name, l_positive )
      write(field_name,"(A,I3)") "pert_qt conv type ", i_type
      call check_bad_values_cmpr( cmpr, k, pert_qt_t(:,i_type),                &
                                  call_string, field_name, l_positive )
    end do

  end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )

end if  ! ( nc2 > 0 )


call cmpr_dealloc( cmpr )


return
end subroutine region_parcel_calcs


end module region_parcel_calcs_mod
