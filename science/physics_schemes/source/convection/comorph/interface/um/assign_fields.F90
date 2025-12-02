! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module assign_fields_mod

use um_types, only: real_umphys

implicit none

contains

! Routine to assign pointers (held in comorph's derived-type structures)
! to the corresponding UM fields, to pass them into comorph
subroutine assign_fields( z_theta, z_rho, p_layer_centres, p_layer_boundaries, &
                          r_theta_levels,                                      &
                          rho_dry_th, w_var_rh, ftl, fqw, fu_rh, fv_rh,        &
                          turb_len, par_radius_amp_um, zh_eff,                 &
                          cca0, ccw0, frac_bulk_conv,                          &
                          u_th_n, v_th_n, w, temperature_n,                    &
                          m_v, m_cl, m_cf,                                     &
                          m_cf2, m_r, m_gr,                                    &
                          cf_liquid_n, cf_frozen_n, bulk_cf_n,                 &
                          u_th_np1, v_th_np1, r_w, theta_star,                 &
                          q_star, qcl_star, qcf_star,                          &
                          qcf2_star, qrain_star, qgraup_star,                  &
                          cf_liquid_star, cf_frozen_star, bulk_cf_star,        &
                          precfrac_star, l_conv_inc_w,  l_temporary_snow,      &
                          l_temporary_rain, l_temporary_graup,                 &
                          w_work, q_snow_work, q_rain_work, q_graup_work,      &
                          grid, turb, cloudfracs, fields_n, fields_np1 )

use atm_fields_bounds_mod, only: tdims, tdims_s, tdims_l, pdims, wdims, wdims_s
use nlsizes_namelist_mod, only: bl_levels

use grid_type_mod, only: grid_type
use turb_type_mod, only: turb_type
use cloudfracs_type_mod, only: cloudfracs_type
use fields_type_mod, only: fields_type

use mphys_inputs_mod, only: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup

use comorph_constants_mod, only: l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,    &
                                 l_cv_cloudfrac, l_turb_par_gen,               &
                                 i_convcloud,                                  &
                                 i_convcloud_liqonly, i_convcloud_mph

use raise_error_mod, only: raise_fatal
use umPrintMgr, only: newline

implicit none

! GRID-RELATED FIELDS

! Model-level heights above surface
real(kind=real_umphys), target, intent(in) :: z_theta                          &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: z_rho                            &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              pdims%k_start:pdims%k_end )
! Model-level pressures
real(kind=real_umphys), target, intent(in) :: p_layer_centres                  &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              tdims%k_start:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: p_layer_boundaries               &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              pdims%k_start:pdims%k_end+1 )
real(kind=real_umphys), target, intent(in) :: r_theta_levels                   &
     (tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,0:tdims%k_end)
! Dry-density on theta-levels
real(kind=real_umphys), target, intent(in) :: rho_dry_th                       &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end-1 )

! TURBULENCE FIELDS

! Turbulent vertical velocity variance on rho-levels
real(kind=real_umphys), target, intent(in) :: w_var_rh                         &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              1:bl_levels )
! Turbulent fluxes of Tl, qw, u, v on rho-levels
real(kind=real_umphys), target, intent(in) :: ftl                              &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              1:bl_levels )
real(kind=real_umphys), target, intent(in) :: fqw                              &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              1:bl_levels )
real(kind=real_umphys), target, intent(in) :: fu_rh                            &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              1:bl_levels )
real(kind=real_umphys), target, intent(in) :: fv_rh                            &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end,                       &
                              1:bl_levels )
! Turbulence lengthscale on theta-levels
real(kind=real_umphys), target, intent(in) :: turb_len                         &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:bl_levels )
! Scaling factor applied to parcel initial radius
real(kind=real_umphys), target, intent(in) :: par_radius_amp_um                &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end )
! Effective boundary-layer top height
real(kind=real_umphys), target, intent(in) :: zh_eff                           &
                            ( pdims%i_start:pdims%i_end,                       &
                              pdims%j_start:pdims%j_end )

! CONVECTIVE CLOUD FIELDS

! Convective liquid cloud fraction and water content
real(kind=real_umphys), target, intent(in) :: cca0                             &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: ccw0                             &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
! Convective bulk cloud fraction
real(kind=real_umphys), target, intent(in) :: frac_bulk_conv                   &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )

! START-OF-TIMESTEP PRIMARY FIELDS

! Winds
real(kind=real_umphys), target, intent(in) :: u_th_n                           &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end-1 )
real(kind=real_umphys), target, intent(in) :: v_th_n                           &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end-1 )
real(kind=real_umphys), target, intent(in) :: w                                &
                            ( wdims_s%i_start:wdims_s%i_end,                   &
                              wdims_s%j_start:wdims_s%j_end,                   &
                              wdims_s%k_start:wdims_s%k_end )
! Temperature
real(kind=real_umphys), target, intent(in) :: temperature_n                    &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
! Water species mixing-ratios
real(kind=real_umphys), target, intent(in) :: m_v                              &
                            ( tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end )
real(kind=real_umphys), target, intent(in) :: m_cl                             &
                            ( tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end )
real(kind=real_umphys), target, intent(in) :: m_cf                             &
                            ( tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end )
real(kind=real_umphys), target, intent(in) :: m_cf2                            &
                            ( tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end )
real(kind=real_umphys), target, intent(in) :: m_r                              &
                            ( tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end )
real(kind=real_umphys), target, intent(in) :: m_gr                             &
                            ( tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end )
! Cloud fractions
real(kind=real_umphys), target, intent(in) :: cf_liquid_n                      &
                            ( tdims_l%i_start:tdims_l%i_end,                   &
                              tdims_l%j_start:tdims_l%j_end,                   &
                              tdims_l%k_start:tdims_l%k_end )
real(kind=real_umphys), target, intent(in) :: cf_frozen_n                      &
                            ( tdims_l%i_start:tdims_l%i_end,                   &
                              tdims_l%j_start:tdims_l%j_end,                   &
                              tdims_l%k_start:tdims_l%k_end )
real(kind=real_umphys), target, intent(in) :: bulk_cf_n                        &
                            ( tdims_l%i_start:tdims_l%i_end,                   &
                              tdims_l%j_start:tdims_l%j_end,                   &
                              tdims_l%k_start:tdims_l%k_end )

! LATEST PRIMARY FIELDS

! Winds
real(kind=real_umphys), target, intent(in) :: u_th_np1                         &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end-1 )
real(kind=real_umphys), target, intent(in) :: v_th_np1                         &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end-1 )
real(kind=real_umphys), target, intent(in) :: r_w                              &
                            ( wdims%i_start:wdims%i_end,                       &
                              wdims%j_start:wdims%j_end,                       &
                              1:wdims%k_end )
! Temperature
real(kind=real_umphys), target, intent(in) :: theta_star                       &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
! Water species mixing-ratios
real(kind=real_umphys), target, intent(in) :: q_star                           &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: qcl_star                         &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: qcf_star                         &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: qcf2_star                        &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: qrain_star                       &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: qgraup_star                      &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
! Cloud fractions
real(kind=real_umphys), target, intent(in) :: cf_liquid_star                   &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: cf_frozen_star                   &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
real(kind=real_umphys), target, intent(in) :: bulk_cf_star                     &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )
! Precipitation fraction
real(kind=real_umphys), target, intent(in) :: precfrac_star                    &
                            ( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end,                       &
                              1:tdims%k_end )


! Switch for whether comorph should increment w
logical, intent(in) :: l_conv_inc_w
! Flags for temporary work arrays for optional water species
logical, intent(in) :: l_temporary_snow
logical, intent(in) :: l_temporary_rain
logical, intent(in) :: l_temporary_graup

! Optional temporary work arrays for certain fields
real(kind=real_umphys), allocatable, target, intent(in) :: w_work(:,:,:)
real(kind=real_umphys), allocatable, target, intent(in) :: q_snow_work(:,:,:)
real(kind=real_umphys), allocatable, target, intent(in) :: q_rain_work(:,:,:)
real(kind=real_umphys), allocatable, target, intent(in) :: q_graup_work(:,:,:)


! Comorph derived-type structures containing pointers, to be assigned
! to the relevant UM fields by this routine...

! Structure storing pointers to the model-grid fields and dry-rho
type(grid_type), intent(in out) :: grid

! Structure containing pointers to turbulence fields passed
! in from the boundary-layer scheme
type(turb_type), intent(in out) :: turb

! Structure containing diagnostic area fractions for
! cloud and rain
type(cloudfracs_type), intent(in out) :: cloudfracs

! Structures storing pointers to the primary fields,
! at start-of-timestep and latest fields
type(fields_type), intent(in out) :: fields_n
type(fields_type), intent(in out) :: fields_np1

! Array bounds, for bounds-specification when assigning pointer
integer :: lb(3)

! Flag for error due to required field not present
logical :: l_error

character(len=*), parameter :: routinename = "ASSIGN_FIELDS"


! assign POINTERS FOR GRID-RELATED FIELDS

! The derived-type structure "grid" contains pointers to the
! required fields.  Point them at the appropriate data...

! Radius of surface from Earth centre
lb = lbound(r_theta_levels)
grid % r_surf( lb(1):, lb(2): ) => r_theta_levels(:,:,0)
! Note: this pointer assignment is technically referencing
! an array section, rather than a full array.  In fortran,
! this results in the bounds specs being lost and
! resetting the array lower-bounds to 1 in the pointer.
! Therefore, the above use of bounds specification is needed
! in order to retain the knowledge of the halos in the
! array pointer.

! "Full" levels correspond to theta-levels
grid % height_full => z_theta
grid % pressure_full => p_layer_centres

! "Half" levels correspond to rho-levels, except at the bottom
! of the model where the bottom half-level is the surface
grid % height_half => z_rho
grid % pressure_half => p_layer_boundaries

! Dry density on full-levels (theta-levels)
grid % rho_dry => rho_dry_th


! assign POINTERS FOR TURBULENCE FIELDS

! If turbulence fields needed
if ( l_turb_par_gen ) then
  ! Assign pointers in turb structure to the required fields
  turb % w_var    => w_var_rh
  turb % f_templ  => ftl
  turb % f_q_tot  => fqw
  turb % f_wind_u => fu_rh
  turb % f_wind_v => fv_rh
  turb % lengthscale => turb_len
  turb % par_radius_amp => par_radius_amp_um(:,:)
end if

! Boundary-layer top height needed even if not using
! turbulence-based perturbations
turb % z_bl_top => zh_eff


! assign POINTERS FOR CLOUD and PRECIP FRACTIONS

! If cloud fractions aren't included as primary fields,
! pass them in through the separate cloudfracs structure
! instead:
if ( .not. l_cv_cloudfrac ) then
  cloudfracs % frac_liq => cf_liquid_star
  cloudfracs % frac_ice => cf_frozen_star
  cloudfracs % frac_bulk => bulk_cf_star
end if

! Precipitation fraction; this is needed if using rain and
! initiation mass-sources based on moist instability
cloudfracs % frac_precip => precfrac_star


! assign POINTERS FOR CONVECTIVE CLOUD FIELDS

! What diagnosed convective cloud options are set in CoMorph?
select case( i_convcloud )
case ( i_convcloud_liqonly )

  ! CoMorph will output only convective liquid cloud
  ! (no diagnostic contribution will be added for ice)
  cloudfracs % frac_liq_conv  => cca0
  cloudfracs % frac_bulk_conv => frac_bulk_conv
  cloudfracs % q_cl_conv      => ccw0
  ! Note: the UM's prognostics cca and ccw correspond to only
  ! the convective liquid cloud, not the ice.
  ! The convective bulk cloud fraction is also output
  ! for the purpose of calculating correct convective
  ! cloud base and top height diagnostics which include the ice,
  ! but this is not yet passed to the rest of the model.

case ( i_convcloud_mph )

  ! CoMorph wants to output separate liquid and ice
  ! convective cloud fields.
  ! But the UM doesn't yet have these fields, so throw a wobbly
  call raise_fatal( routinename,                                               &
         "Calling the CoMorph convection scheme with "        //               &
         "separate convective liquid and ice cloud " //newline//               &
         "fields, but the UM currently only supports a "      //               &
         "single convective cloud fraction field." )

end select


! assign POINTERS FOR PRIMARY FIELDS

! The derived-type structure "fields_n" contains pointers to the
! required fields at start-of-timestep, and "fields_np1" contains
! pointers to the latest fields to be updated.
! Point them at the appropriate data...

! Start-of-timestep winds
fields_n % wind_u => u_th_n
fields_n % wind_v => v_th_n
fields_n % wind_w => w

! Latest winds
fields_np1 % wind_u => u_th_np1
fields_np1 % wind_v => v_th_np1
if ( l_conv_inc_w ) then
  ! Only point to the actual array for w if w increments are switched on
  fields_np1 % wind_w => r_w
else
  ! Otherwise point to a temporary copy
  fields_np1 % wind_w => w_work
end if

! Temperature (note theta_star now holds latest T not theta)
fields_n % temperature => temperature_n
fields_np1 % temperature => theta_star


! Water species...

! Vapour and liquid are always on in both the UM and CoMorph
fields_n % q_vap => m_v
fields_n % q_cl => m_cl
fields_np1 % q_vap => q_star
fields_np1 % q_cl => qcl_star

! Species that are optional in comorph...

! Ice-cloud (always in use in the UM)
! Which UM field maps onto CoMorph's q_cf depends on UM's 2nd ice category
if ( l_cv_cf ) then
  if ( l_mcr_qcf2 ) then
    ! The UM has 2nd ice category switched on
    ! In this case, CoMorph's q_cf corresponds to the UM's qcf2 (crystals)
    fields_n % q_cf => m_cf2
    fields_np1 % q_cf => qcf2_star
  else
    ! The UM has 2nd ice category switched off
    ! CoMorph's q_cf maps onto the UM's qcf
    fields_n % q_cf => m_cf
    fields_np1 % q_cf => qcf_star
  end if
end if

! Species that are optional in the UM; if in use in comorph, point to
! the UM's prognostic array if the species is in use in the UM too,
! or point to a temporary work array if not....
l_error = .false.

! Snow
if ( l_cv_snow ) then
  if ( l_temporary_snow ) then
    fields_n % q_snow => q_snow_work
    fields_np1 % q_snow => q_snow_work
  else if ( l_mcr_qcf2 ) then
    ! The UM has 2nd ice category switched on;
    ! In this case, comorph's q_snow corresponds to the UM's qcf (aggregates)
    fields_n % q_snow => m_cf
    fields_np1 % q_snow => qcf_star
  else
    l_error = .true.
  end if
end if

! Rain
if ( l_cv_rain ) then
  if ( l_temporary_rain ) then
    fields_n % q_rain => q_rain_work
    fields_np1 % q_rain => q_rain_work
  else if ( l_mcr_qrain ) then
    fields_n % q_rain => m_r
    fields_np1 % q_rain => qrain_star
  else
    l_error = .true.
  end if
end if

! Graupel
if ( l_cv_graup ) then
  if ( l_temporary_graup ) then
    fields_n % q_graup => q_graup_work
    fields_np1 % q_graup => q_graup_work
  else if ( l_mcr_qgraup ) then
    fields_n % q_graup => m_gr
    fields_np1 % q_graup => qgraup_star
  else
    l_error = .true.
  end if
end if

if ( l_error ) then
  call raise_fatal( routinename,                                               &
         "At least one water species field expected by comorph is "         // &
         "not available in the UM (snow, rain or graupel)." )
end if


! Cloud-fractions
if ( l_cv_cloudfrac ) then
  fields_n % cf_liq => cf_liquid_n
  fields_n % cf_ice => cf_frozen_n
  fields_n % cf_bulk => bulk_cf_n
  fields_np1 % cf_liq => cf_liquid_star
  fields_np1 % cf_ice => cf_frozen_star
  fields_np1 % cf_bulk => bulk_cf_star
end if


return
end subroutine assign_fields

end module assign_fields_mod
