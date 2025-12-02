! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     SUBROUTINES LS_PPN and LS_PPNC------------------------------------
!    Purpose:
!            LS_PPN and LS_PPNC:
!             Calculate large-scale (dynamical) precipitation.
!             LS_PPNC is the gather/scatter routine which then
!             calls LSP_ICE.
!    Note: in all cases, level counters (incl subscripts) run from 1
!          (lowest model layer) to LEVELS (topmost  model layer)
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Documentation: UM Documentation Paper 26.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation
module ls_ppnc_mod

  ! General atmosphere modules
use planet_constants_mod, only: g

use murk_inputs_mod, only: l_murk_advect

! Microphysics modules
use mphys_diags_mod,       only: l_aggfr_diag,                                 &
                                 l_psdep_diag, l_psaut_diag,                   &
                                 l_psacw_diag, l_psacr_diag,                   &
                                 l_psaci_diag, l_psmlt_diag,                   &
                                 l_psmltevp_diag, l_praut_diag,                &
                                 l_pracw_diag, l_prevp_diag,                   &
                                 l_pgaut_diag, l_pgacw_diag,                   &
                                 l_pgacs_diag, l_pgmlt_diag,                   &
                                 l_pifrw_diag, l_piprm_diag,                   &
                                 l_piprr_diag,                                 &
                                 l_pidep_diag, l_piacw_diag,                   &
                                 l_piacr_diag, l_pimlt_diag,                   &
                                 l_pimltevp_diag, l_pifall_diag,               &
                                 l_psfall_diag, l_prfall_diag,                 &
                                 l_pgfall_diag, l_plset_diag,                  &
                                 l_plevpset_diag, l_pifrr_diag,                &
                                 l_sfwater_diag, l_sfrain_diag,                &
                                 l_sfsnow_diag,                                &
                                 l_vtbranch_diag, l_vm_cry_diag,               &
                                 l_vm_agg_diag,                                &
                                 l_vm_used_diag,                               &
                                 !All REALs in the module have kind=real_lsprec
                                 psdep, psaut, frac_agg,                       &
                                 psacw, psacr, psaci, psmlt,                   &
                                 psmltevp, praut, pracw, prevp,                &
                                 pgaut, pgacw, pgacs, pgmlt, pifrw,            &
                                 piprm, piprr, pidep,                          &
                                 piacw, piacr, pimlt,                          &
                                 pimltevp, pifall, psfall, prfall,             &
                                 pgfall, plset, plevpset, pifrr,               &
                                 vm_cry, vm_agg, vtbranch_flag,                &
                                 vm_used, dbz_tot, dbz_g, dbz_i,               &
                                 dbz_i2, dbz_l, dbz_r,                         &
                                 sfwater, sfrain, sfsnow

use mphys_inputs_mod,      only: l_mcr_qcf2, l_mcr_qgraup,                     &
                                 l_mcr_qrain, l_mcr_precfrac, l_fsd_generator, &
                                 l_progn_tnuc

use lsprec_mod,            only: zero

use tuning_segments_mod,   only: precip_segment_size

use mphys_bypass_mod,      only: l_ref_diag

use cloud_inputs_mod,      only:  i_cld_vn
use pc2_constants_mod,     only:  i_cld_pc2

! Grid bounds module

use atm_fields_bounds_mod, only: tdims, pdims, tdims_l, pdims_s

! Dr Hook modules
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

! Large scale precipitation modules
use lsp_ice_mod,           only: lsp_ice
use lsp_scav_mod,          only: lsp_scav

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use wtrac_mphys_mod,         only: mp_wtrac_type, mp_cpr_wtrac_type,           &
                                   ls_ppnc_gather_wtrac, ls_ppnc_scatter_wtrac

! Use in KINDs
! real_umphys is the generic kind for physics variables (generally real64)
! real_lsprec is defined at compile time (generally real64 or real32) and
!             controls the precision for all the variables (and therefore
!             calculations) for science called from this module.
!             Many diagnostics are kept at real_lsprec as they are not needed
!             other than for passing into copydiag
use um_types,              only: real_lsprec, real_umphys

implicit none

integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1
character(len=*),   parameter, private :: ModuleName='LS_PPNC_MOD'

!-----------------------------------------------------------------------
!  Enumerators with the indices for the super arrays
!
!  The variables in an enumerator have sequentially increasing integer
!  values. These variables represent the indices of the various data
!  in the super array used by a segment.
!
!  How to add a variable passed from ls_ppn to ls_ppnc.
!
!  1) Add  the original field as argument into the ls_ppnc routine.
!
!  2) Add an extra element in the enumerator structure before the
!  last_index_??? element and after the element set as "= 1". The
!  enumerator that needs to be changed depends on the type of the
!  variable. For example, to add the logical variable "var", the
!  "var_i" element need to be added after "bland_i" into the logical
!  enumerator. The size of the logical super-array will increase
!  automatically since the value of last_index_??? is the total number
!  of elements in the enumerator.
!
!  3) Update the gather/scatter routines to take as argument the
!  original field and the new compressed array. Update the lsp
!  subroutines to use the new compressed sub-array. In the previous
!  example, the new var compressed array located in the logical
!  super-array will be passed to the sub-routines using the code
!  "spl(:,var_i)".
!
!  At the moment there are 2 super arrays: spr for real values and spl
!  for logical values. If a new type of compressed arrays is needed,
!  for example INTEGERS, then a new super array "spi" and a new
!  integer enumerator should be crated.
!
!-----------------------------------------------------------------------

! Logical super-array enumerator
enum, bind(c)
  enumerator ::                                                                &

    bland_i = 1,                                                               &
                ! gathered land/sea mask
    last_index_log
                ! Used to find the total number of variables in the
                ! supper array, i.e. this value minus 1.
end enum

! Real super-array enumerator
enum, bind(c)
  enumerator ::                                                                &
   cf_i = 1,                                                                   &
                        ! gathered Cloud fraction.
   q_i,                                                                        &
                         ! gathered Specific humidity (kg water/kg air).
   qcf_i,                                                                      &
                         ! gathered Cloud ice (kg per kg air).
   qcl_i,                                                                      &
                         ! gathered Cloud liquid water (kg per kg air).
   tnuc_new_i,                                                                 &
                         ! gathered nucleation temperature as function of dust
                         ! deg cel
   qcf2_i,                                                                     &
                          ! gathered cloud ice2 (kg per kg air).
   qrain_i,                                                                    &
                           ! gathered rain (kg per kg air).
   qgraup_i,                                                                   &
                            ! gathered graupel (kg per kg air).
   t_i,                                                                        &
                         ! gathered Temperature (K).
   uk_i,                                                                       &
                         ! gathered u wind on level k
   vk_i,                                                                       &
                         ! gathered v wind on level k
   ukp1_i,                                                                     &
                         ! gathered u wind on level k+1
   vkp1_i,                                                                     &
                         ! gathered v wind on level k+1
   r_theta_levels_i,                                                           &
                         ! gathered distance from centre of Earth
   r_theta_surf_i,                                                             &
                         ! ...and at surface
   g_cos_theta_latitude_i,                                                     &
                         ! gathered cos of latitude
   aero_i,                                                                     &
                         ! gathered Aerosol.
   lsrain_i,                                                                   &
                   !gathered Surface rainfall rate (kg per sq m per s).
   lssnow_i,                                                                   &
                   !gathered Surface snowfall rate (kg per sq m per s).
   lssnow2_i,                                                                  &
                    !gathered layer snowfall rate (kg per sq m per s).
   lsgraup_i,                                                                  &
                    !gathered layer graupel fall rate (kg/sq m/s)
   droplet_flux_i,                                                             &
                         ! gathered water droplet flux / kg m-2 s-1
   cttemp_i,                                                                   &
                              !gathered ice cloud top temperature.
   rainfrac_i,                                                                 &
                              !gathered rain fraction.
   rainfrac_impr_i,                                                            &
                              !gathered improved rain fraction
   precfrac_k_i,                                                               &
                              !gathered prognostic precipitation fraction at k
   precfrac_fall_i,                                                            &
                              !gathered falling precipitation fraction at k+-1/2
   frac_ice_above_i,                                                           &
                              !gathered fraction of ice in layer above
   frac_agg_i,                                                                 &
                         ! gathered aggregate fraction
   cfl_i,                                                                      &
                         ! gathered Cloud liquid fraction.
   cff_i,                                                                      &
                         ! gathered Cloud ice fraction.
   vfall_i,                                                                    &
                         ! gathered fall velocity (m per s).
   vfall2_i,                                                                   &
                          ! gathered fall velocity for qcf2 (m per s).
   vfall_rain_i,                                                               &
                          ! gathered fall velocity for qcf2 (m per s).
   vfall_graup_i,                                                              &
                          ! gathered fall velocity for qcf2 (m per s).
   rhc_i,                                                                      &
                          ! gathered RH_crit value at points.
   n_drop_tpr_i,                                                               &
   n_drop_out_i,                                                               &
                          ! gathered droplet numbers
   hmteff_i,                                                                   &
                          ! gathered effective mountain height
   zb_i,                                                                       &
                          ! gathered blocked layer depth for seederfeeder
   f_arr1_i,                                                                   &
   f_arr2_i,                                                                   &
   f_arr3_i,                                                                   &

    deltaz_i,                                                                  &
                         ! Thickness of layer (m)
     rhodz_dry_i,                                                              &
    ! Dry air density * layer thickness (kg m-2)
    rhodz_moist_i,                                                             &
                         ! Moist air density * layer thickness (kg m-2)

     p_i,                                                                      &
                    ! WORK Used for pressure at successive levels.



  ! Microphysical process rate diagnostics (compressed arrays)
    psdep_i,                                                                   &
                    ! Deposition of vapour to snow aggregates
    psaut_i,                                                                   &
                    ! Autoconversion of aggregates from crystals
    psacw_i,                                                                   &
                    ! Accretion of liq. water by snow aggregates
    psacr_i,                                                                   &
                    ! Collection of rain by snow aggregates
    psaci_i,                                                                   &
                    ! Collection of ice crystals by aggregates
    psmlt_i,                                                                   &
                    ! Melting of snow aggregates
    psmltevp_i,                                                                &
                    ! Evaporation of melting aggregates
    praut_i,                                                                   &
                    ! Autoconversion of cloud drops to rain
    pracw_i,                                                                   &
                    ! Accretion of liq. water by rain
    prevp_i,                                                                   &
                    ! Evaporation of rain
    pgaut_i,                                                                   &
                    ! Autoconversion of graupel from aggregates
    pgacw_i,                                                                   &
                    ! Accretion of liq. water by graupel
    pgacs_i,                                                                   &
                    ! Collection of snow aggregates by graupel
    pgmlt_i,                                                                   &
                    ! Melting of graupel
    pifrw_i,                                                                   &
                    ! Homogeneous freezing nucleation
    pifrr_i,                                                                   &
                    ! Homogeneous freezing of rain
    piprm_i,                                                                   &
                    ! Heterogeneous (primary) nucleation
    piprr_i,                                                                   &
                    ! Heterogeneous nucleation of rain
    pidep_i,                                                                   &
                    ! Deposition of vapour to ice crystals
    piacw_i,                                                                   &
                    ! Accretion of liq. water by ice crystals
    piacr_i,                                                                   &
                    ! Collection of rain by ice crystals
    pimlt_i,                                                                   &
                    ! Melting of ice crystals
    pimltevp_i,                                                                &
                    ! Evaporation of melting ice crystals
    pifall_i,                                                                  &
                    ! Sedimentation of ice crystals
    psfall_i,                                                                  &
                    ! Sedimentation of aggregates
    prfall_i,                                                                  &
                    ! Sedimentation of rain
    pgfall_i,                                                                  &
                    ! Sedimentation of graupel
    plset_i,                                                                   &
                    ! Droplet settling of liquid water
    plevpset_i,                                                                &
                    ! Evaporated settled droplets
    vm_agg_i,                                                                  &
                 ! Mass-weighted fallspeed with aggregate parameters
    vm_cry_i,                                                                  &
                 ! Mass-weighted fallspeed with crystal parameters
    vtbranch_flag_i,                                                           &
                 ! Flag indicating which vt-D is used
    vm_used_i,                                                                 &
                 ! Mass-weighted fallspeed used

    sfwater_i,                                                                 &
                 ! Seeder feeder orographic water mixing ratio
    sfrain_i,                                                                  &
                 ! Seeder feeder orographic rain production
    sfsnow_i,                                                                  &
                 ! Seeder feeder orographic snow production

    dbz_tot_i,                                                                 &
                 ! Total reflectivity (dBZ)
    dbz_g_i,                                                                   &
                ! Graupel reflectivity (dBZ)
    dbz_i_i,                                                                   &
                ! Ice Agg. reflectivity (dBZ)
    dbz_i2_i,                                                                  &
                ! Ice Cry. reflectivity (dBZ)
    dbz_l_i,                                                                   &
                ! Cloud liquid reflectivity (dBZ)
    dbz_r_i,                                                                   &
                ! Rain reflectivity (dBZ)

    last_index_real
                ! Used to find the total number of variables in the
                ! supper array, i.e. this value minus 1.
end enum

! A LSP segment data is composed of:
! 1) segment_size * lsp_num_real  real    values
! 2) segment_size * lsp_num_log   logical values


! The total number of variables in the LSP super array.
integer, parameter :: lsp_num_real = last_index_real - 1
integer, parameter :: lsp_num_log  = last_index_log  - 1

contains

!-----------------------------------------------------------------------
!  Main routine
!-----------------------------------------------------------------------

subroutine ls_ppnc( level, ix, n,                                              &
 lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                                  &
 cf,cfl,cff,                                                                   &
 qcf,qcl,tnuc_new,t, qcf2,qrain,qgraup,                                        &
 n_drop_tpr, n_drop_out,                                                       &
 aerosol,                                                                      &
 hmteff, zb,                                                                   &
!--------------------------------------------------------
! Layer thicknesses and variables passed from layer above
!--------------------------------------------------------
 q, p_theta_levels, r_theta_levels,                                            &
 deltaz, rhodz_dry, rhodz_moist,                                               &
 rhc_row_length, rhc_rows, bland, rhcrit, f_arr, cos_theta_latitude,           &
 vfall, vfall2, vfall_rain, vfall_graup,                                       &
 frac_ice_above, cttemp, rainfrac, rainfrac_impr, precfrac_k, precfrac_fall,   &
 niters_mp,                                                                    &
 uk, vk, ukp1, vkp1, onelevel,                                                 &
!--------------------------------------------------------
! COSP flag and water tracers
!--------------------------------------------------------
 l_cosp_lsp, wtrac_mp                                                          &
 )


implicit none

integer ::                                                                     &
 level,                                                                        &
              ! level number
 n,                                                                            &
              ! in Number of points where pptn non-zero from above
!                    or where CF>CFMIN
   niters_mp,                                                                  &
                ! Number of total iterations of microphysics
   rhc_row_length,rhc_rows,                                                    &
   ix (  tdims%i_len  *                                                        &
         tdims%j_len  , 2 )
                                ! in gather/scatter index
real(kind=real_umphys) ::                                                      &
 cf (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                    &
                            ! in Cloud fraction.
 cfl(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                    &
                            ! in Cloud liquid fraction.
 cff(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                    &
                            ! in Cloud ice fraction.
! CF, CFL and CFF are in/out if the PC2 cloud scheme is in use.

    p_theta_levels( tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end ),                               &
    r_theta_levels( tdims_l%i_start:tdims_l%i_end,                             &
                    tdims_l%j_start:tdims_l%j_end,                             &
                    tdims_l%k_start:tdims_l%k_end ),                           &
                                ! in height of theta levels from centre of Earth
    deltaz( tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end ),                                       &
                                          ! in thickness of layer (m)
    rhodz_dry(tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end ),                                     &
                                          ! Dry air density
                                          ! * layer thickness (kg m-2)
    rhodz_moist(tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end ),                                   &
                                          ! Moist air density
                                          ! * layer thickness (kg m-2)

    rhcrit(rhc_row_length,rhc_rows),                                           &
!                       in Critical humidity for cloud formation.
    f_arr(3, tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),            &
                      ! in parameters used in fractional standard deviation calc
    cos_theta_latitude(pdims_s%i_start:pdims_s%i_end,                          &
                       pdims_s%j_start:pdims_s%j_end)
                      ! in cosine of the current latitude

logical, intent(in) ::                                                         &
     bland(tdims%i_start:tdims%i_end,                                          &
           tdims%j_start:tdims%j_end )
                        !in Land/sea mask

logical, intent(in) ::  l_cosp_lsp
                        ! Switch for COSP LS diagnostics

real(kind=real_umphys), intent(in out) ::                                      &
     q(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end),                                             &
                          ! INOUT Specific humidity (kg water/kg air).
      qcf(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end),                                          &
                          ! INOUT Cloud ice (kg per kg air).
      qcl(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end),                                          &
                          ! INOUT Cloud liquid water (kg per kg air).
      qcf2(tdims%i_start:tdims%i_end,                                          &
           tdims%j_start:tdims%j_end),                                         &
                          ! INOUT Cloud ice2 (kg per kg air).
      qrain(tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end),                                        &
                          ! INOUT Rain water (kg per kg air).
      qgraup(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end),                                       &
                          ! INOUT Graupel water (kg per kg air).
      t(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end),                                            &
                          ! INOUT Temperature (K).
      aerosol(tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end),                                      &
                          ! INOUT Aerosol (K).
      lsrain(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end),                                       &
                          !INOUT Surface rainfall rate (kg m^-2 s^-1).
      lssnow(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end),                                       &
                          !INOUT Surface snowfall rate (kg m^-2 s^-1).
      lssnow2(tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end),                                      &
                          !INOUT layer snowfall rate (kg m^-2 s^-1).
      lsgraup(tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end),                                      &
                          !INOUT layer graupelfall rate (kg m^-2 s^-1)
      droplet_flux(tdims%i_start:tdims%i_end,                                  &
                   tdims%j_start:tdims%j_end),                                 &
                          !INOUT water droplet flux / kg m^-2 s^-1
      cttemp(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end),                                       &
                          ! INOUT Ice cloud top temperature (K)
      rainfrac(tdims%i_start:tdims%i_end,                                      &
               tdims%j_start:tdims%j_end),                                     &
                          ! INOUT Rain fraction.
      precfrac_k(tdims%i_start:tdims%i_end,                                    &
                 tdims%j_start:tdims%j_end),                                   &
                          ! INOUT Precipitation fraction at theta-level k
      precfrac_fall(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end),                                &
                          ! INOUT Falling precipitation fraction at rho-levels
      frac_ice_above(tdims%i_start:tdims%i_end,                                &
                     tdims%j_start:tdims%j_end),                               &
                          ! INOUT Ice fraction from layer above
 vfall(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end),                                             &
                                   ! INOUT fall velocity of ice (m per
 vfall2(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end),                                            &
                                    ! INOUT fall vel. of rain (m/s)
 vfall_rain(tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end),                                        &
                                    ! INOUT fall vel. of rain (m/s)
 vfall_graup(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end)
                                    ! INOUT fall vel. of graupel (m/s)
real (kind=real_umphys), intent(in) ::                                         &
    tnuc_new(tdims%i_start:tdims%i_end,                                        &
    tdims%j_start:tdims%j_end)
! tnuc_new as function of dust (deg cel)

real(kind=real_umphys), intent(in) ::                                          &
                    rainfrac_impr( tdims%i_start:tdims%i_end,                  &
                                   tdims%j_start:tdims%j_end )
                                    ! Improved rain fraction

real(kind=real_umphys),intent(in) ::                                           &
! U and V wind at levels k and k+1 on P grid, so defined using pdims.
   uk  (pdims%i_start : pdims%i_end,                                           &
        pdims%j_start : pdims%j_end),                                          &
   ukp1(pdims%i_start : pdims%i_end,                                           &
        pdims%j_start : pdims%j_end),                                          &
   vk  (pdims%i_start : pdims%i_end,                                           &
        pdims%j_start : pdims%j_end),                                          &
   vkp1(pdims%i_start : pdims%i_end,                                           &
        pdims%j_start : pdims%j_end)
! Level we are interested in for r_theta_level
integer, intent(in) :: onelevel

real(kind=real_umphys),intent(in) ::                                           &
     hmteff(tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end)
real(kind=real_umphys),intent(in) ::                                           &
     zb(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end)

real(kind=real_umphys) ::                                                      &
                    !, Intent(in)
      n_drop_tpr( tdims%i_start : tdims%i_end,                                 &
                  tdims%j_start : tdims%j_end )
                             ! Tapered droplet number
real(kind=real_umphys) ::                                                      &
                    ! intent(in out)
      n_drop_out( tdims%i_start : tdims%i_end,                                 &
                  tdims%j_start : tdims%j_end )
                    ! droplet number from autoconversion

! Water tracers working arrays
type(mp_wtrac_type), intent(in out) :: wtrac_mp(n_wtrac)

!    Workspace usage ---------------------------------------------------

! Super arrays containing all the data used by a segment.
real (kind=real_lsprec), allocatable :: spr(:,:)
logical,                 allocatable :: spl(:,:)

! Water tracer compressed structure
type(mp_cpr_wtrac_type), allocatable :: wtrac_mp_cpr(:)

integer :: jj,i1,i2,ip

real(kind=jprb)             :: zhook_handle
character(len=*), parameter :: RoutineName='LS_PPNC'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(SHARED) private(jj, i1, i2, ip, spr, spl, wtrac_mp_cpr)

!$OMP do SCHEDULE(DYNAMIC)
do jj = 1, n, precip_segment_size

  i1 = jj
  i2 = min(jj+precip_segment_size-1, n)
  ip = i2-i1+1

  !-----------------------------------------------------------------------
  !  Allocate super arrays
  !-----------------------------------------------------------------------

  allocate( spl (ip, lsp_num_log)    )
  allocate( spr (ip, lsp_num_real)   )
  allocate(wtrac_mp_cpr(n_wtrac))

  !-----------------------------------------------------------------------
  !  Gather variables using index
  !-----------------------------------------------------------------------

  call ls_ppnc_gather( level, ix, ip, i1,                                      &
    lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                               &
    cf,cfl,cff,                                                                &
    qcf,qcl,tnuc_new,t, qcf2,qrain,qgraup,                                     &
    n_drop_tpr, n_drop_out,                                                    &
    aerosol,                                                                   &
    hmteff, zb,                                                                &
    q, p_theta_levels, r_theta_levels,                                         &
    deltaz, rhodz_dry, rhodz_moist,                                            &
    rhc_row_length, rhc_rows, bland, rhcrit, f_arr, cos_theta_latitude,        &
    vfall, vfall2, vfall_rain, vfall_graup,                                    &
    frac_ice_above, cttemp, rainfrac, rainfrac_impr, precfrac_k, precfrac_fall,&
    uk, vk, ukp1, vkp1, onelevel, wtrac_mp,                                    &
    !Compressed super arrays
    spr,spl,wtrac_mp_cpr)

  !-----------------------------------------------------------------------
  ! Call subroutine lsp_ice to control transfers between different
  ! microphysical species
  !-----------------------------------------------------------------------

  call lsp_ice(                                                                &
    spr(:,p_i), spr(:,deltaz_i), spr(:,rhodz_dry_i),                           &
    spr(:,rhodz_moist_i), ip, spr(:,rhc_i),                                    &
    spr(:,hmteff_i), spr(:,zb_i), spr(:,qcf_i),                                &
    spr(:,qcl_i), spr(:,tnuc_new_i), spr(:,q_i), spr(:,qcf2_i), spr(:,qrain_i),&
    spr(:,qgraup_i), spr(:,n_drop_tpr_i), spr(:,n_drop_out_i),                 &
    spr(:,lsrain_i), spr(:,vfall_rain_i), spr(:,lssnow_i),                     &
    spr(:,vfall_i), spr(:,lssnow2_i), spr(:,vfall2_i),                         &
    spr(:,lsgraup_i), spr(:,vfall_graup_i), spr(:,droplet_flux_i),             &
    spr(:,frac_ice_above_i), spr(:,frac_agg_i), spr(:,cttemp_i),               &
    spr(:,rainfrac_i), spr(:,rainfrac_impr_i),                                 &
    spr(:,precfrac_k_i), spr(:,precfrac_fall_i), spr(:,t_i),                   &
    spr(:,cf_i), spr(:,cfl_i), spr(:,cff_i), spl(:,bland_i),                   &
    spr(:,psdep_i), spr(:,psaut_i), spr(:,psacw_i), spr(:,psacr_i),            &
    spr(:,psaci_i), spr(:,psmlt_i), spr(:,psmltevp_i),                         &
    spr(:,praut_i), spr(:,pracw_i), spr(:,prevp_i), spr(:,pgaut_i),            &
    spr(:,pgacw_i), spr(:,pgacs_i), spr(:,pgmlt_i), spr(:,pifrw_i),            &
    spr(:,pifrr_i), spr(:,piprm_i), spr(:,piprr_i), spr(:,pidep_i),            &
    spr(:,piacw_i), spr(:,piacr_i), spr(:,pimlt_i), spr(:,pimltevp_i),         &
    spr(:,pifall_i), spr(:,psfall_i), spr(:,prfall_i),                         &
    spr(:,pgfall_i), spr(:,plset_i), spr(:,plevpset_i),                        &
    spr(:,dbz_tot_i), spr(:,dbz_g_i), spr(:,dbz_i_i),                          &
    spr(:,dbz_i2_i), spr(:,dbz_l_i), spr(:,dbz_r_i),                           &
    spr(:,sfwater_i), spr(:,sfrain_i), spr(:,sfsnow_i), niters_mp,             &
    spr(:,uk_i), spr(:,vk_i), spr(:,ukp1_i), spr(:,vkp1_i),                    &
    spr(:,r_theta_levels_i), spr(:,g_cos_theta_latitude_i),                    &
    spr(:,r_theta_surf_i), spr(:,f_arr1_i), spr(:,f_arr2_i),                   &
    spr(:,f_arr3_i), spr(:,vm_cry_i), spr(:,vm_agg_i),                         &
    spr(:,vtbranch_flag_i), spr(:,vm_used_i), wtrac_mp_cpr                     &
    )

  !-----------------------------------------------------------------------
  ! Lose Murk aerosol by scavenging: call LSP_SCAV
  !-----------------------------------------------------------------------

  if (l_murk_advect) then
    call lsp_scav(ip,spr(:,lsrain_i), spr(:,lssnow_i), spr(:,droplet_flux_i),  &
                  spr(:,aero_i))
  end if

  !-----------------------------------------------------------------------
  ! Scatter back arrays which will have been changed.
  !
  ! Default:    variables are changed from real_lsprec to real_umphys
  ! Exceptions: diagnostics in mphys_diags_mod which are always real_lsprec
  !-----------------------------------------------------------------------

  call ls_ppnc_scatter( level, ix, ip, i1,                                     &
    lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                               &
    cf,cfl,cff,                                                                &
    qcf,qcl,t, qcf2,qrain,qgraup,                                              &
    n_drop_out,                                                                &
    aerosol,                                                                   &
    q,                                                                         &
    vfall, vfall2, vfall_rain, vfall_graup,                                    &
    frac_ice_above, cttemp, rainfrac, precfrac_k, precfrac_fall,               &
    l_cosp_lsp, wtrac_mp,                                                      &
    !Compressed super arrays
    spr, wtrac_mp_cpr)

  !-----------------------------------------------------------------------
  !  Deallocate arrays
  !-----------------------------------------------------------------------

  deallocate( wtrac_mp_cpr )
  deallocate( spr )
  deallocate( spl )

end do ! Loop over points
!$OMP end do

!$OMP end PARALLEL


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_ppnc


!-----------------------------------------------------------------------
!  Gather variables, i.e initialise the data
!-----------------------------------------------------------------------

subroutine ls_ppnc_gather( level, ix, ip, i1,                                  &
  lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                                 &
  cf,cfl,cff,                                                                  &
  qcf,qcl,tnuc_new, t, qcf2,qrain,qgraup,                                      &
  n_drop_tpr, n_drop_out,                                                      &
  aerosol,                                                                     &
  hmteff, zb,                                                                  &
!--------------------------------------------------------
! Variables passed from layer above
!--------------------------------------------------------
  q, p_theta_levels, r_theta_levels,                                           &
  deltaz, rhodz_dry, rhodz_moist,                                              &
  rhc_row_length, rhc_rows, bland, rhcrit, f_arr, cos_theta_latitude,          &
  vfall, vfall2, vfall_rain, vfall_graup,                                      &
  frac_ice_above, cttemp, rainfrac, rainfrac_impr, precfrac_k, precfrac_fall,  &
  uk, vk, ukp1, vkp1, onelevel, wtrac_mp,                                      &
  !Compressed super arrays
  spr,spl,wtrac_mp_cpr)


implicit none

integer, intent(in) ::                                                         &
 level,                                                                        &
              ! level number
 ip,                                                                           &
              ! Number of points
 i1,                                                                           &
              ! Starting point
 rhc_row_length,                                                               &
 rhc_rows,                                                                     &
 ix (  tdims%i_len * tdims%j_len , 2 )
              ! Gather/scatter index

real(kind=real_umphys), intent(in) ::                                          &
  cf (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                   &
              ! Cloud fraction.
  cfl(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                   &
              ! Cloud liquid fraction.
  cff(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                   &
              ! Cloud ice fraction.
  p_theta_levels( tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end ),                                 &
              ! Pressure on theta levels (Pa)
  r_theta_levels( tdims_l%i_start:tdims_l%i_end,                               &
                  tdims_l%j_start:tdims_l%j_end,                               &
                  tdims_l%k_start:tdims_l%k_end ),                             &
  deltaz( tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end ),                                         &
              ! Thickness of layer (m)
  rhodz_dry(tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end ),                                       &
              ! Dry air density
              ! * layer thickness (kg m-2)
  rhodz_moist(tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end ),                                     &
              ! Moist air density
              ! * layer thickness (kg m-2)
  rhcrit(rhc_row_length,rhc_rows),                                             &
              ! Critical humidity for cloud formation.
  f_arr(3, tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),              &
  cos_theta_latitude(pdims_s%i_start:pdims_s%i_end,                            &
                     pdims_s%j_start:pdims_s%j_end)

logical, intent(in) ::                                                         &
  bland(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end )
              ! Land/sea mask

real(kind=real_umphys), intent(in) ::                                          &
  q(tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end),                                                &
              ! Specific humidity (kg water/kg air).
  qcf(tdims%i_start:tdims%i_end,                                               &
      tdims%j_start:tdims%j_end),                                              &
              ! Cloud ice (kg per kg air).
  qcl(tdims%i_start:tdims%i_end,                                               &
      tdims%j_start:tdims%j_end),                                              &
              ! Cloud liquid water (kg per kg air).
  tnuc_new(tdims%i_start:tdims%i_end,                                          &
      tdims%j_start:tdims%j_end),                                              &
              !tnuc_new as function of dust (deg cel)
  qcf2(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end),                                             &
              ! Cloud ice2 (kg per kg air).
  qrain(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end),                                            &
              ! Rain water (kg per kg air).
  qgraup(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Graupel water (kg per kg air).
  t(tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end),                                                &
              ! Temperature (K).
  aerosol(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end),                                          &
              ! Aerosol (K).
  lsrain(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Surface rainfall rate (kg m^-2 s^-1).
  lssnow(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Surface snowfall rate (kg m^-2 s^-1).
  lssnow2(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end),                                          &
              ! Layer snowfall rate (kg m^-2 s^-1).
  lsgraup(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end),                                          &
              ! Layer graupelfall rate (kg m^-2 s^-1)
  droplet_flux(tdims%i_start:tdims%i_end,                                      &
               tdims%j_start:tdims%j_end),                                     &
              ! Water droplet flux / kg m^-2 s^-1
  cttemp(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Ice cloud top temperature (K)
  rainfrac(tdims%i_start:tdims%i_end,                                          &
           tdims%j_start:tdims%j_end),                                         &
              ! Rain fraction.
  precfrac_k(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end),                                       &
              ! Precipitation fraction at theta-level k
  precfrac_fall(tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
              ! Falling precipitation fraction at rho-levels
  frac_ice_above(tdims%i_start:tdims%i_end,                                    &
                 tdims%j_start:tdims%j_end),                                   &
              ! Ice fraction from layer above
  vfall(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end),                                            &
              ! Fall velocity of ice (m per
  vfall2(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Fall vel. of rain (m/s)
  vfall_rain(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end),                                       &
              ! Fall vel. of rain (m/s)
  vfall_graup(tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end)
              ! Fall vel. of graupel (m/s)

real(kind=real_umphys), intent(in) ::                                          &
  rainfrac_impr(tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end)
              ! Improved rain fraction

real(kind=real_umphys),intent(in) ::                                           &
  ! U and V wind at levels k and k+1 on P grid, so defined using pdims.
  uk  (pdims%i_start : pdims%i_end,                                            &
       pdims%j_start : pdims%j_end),                                           &
  ukp1(pdims%i_start : pdims%i_end,                                            &
       pdims%j_start : pdims%j_end),                                           &
  vk  (pdims%i_start : pdims%i_end,                                            &
       pdims%j_start : pdims%j_end),                                           &
  vkp1(pdims%i_start : pdims%i_end,                                            &
       pdims%j_start : pdims%j_end)

! Level we are interested in for r_theta_level
integer, intent(in) :: onelevel

real(kind=real_umphys),intent(in) ::                                           &
  hmteff(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end)
real(kind=real_umphys),intent(in) ::                                           &
  zb(tdims%i_start:tdims%i_end,                                                &
     tdims%j_start:tdims%j_end)

real(kind=real_umphys),intent(in) ::                                           &
  n_drop_tpr(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end )
              ! Tapered droplet number
real(kind=real_umphys), intent(in) ::                                          &
  n_drop_out( tdims%i_start : tdims%i_end,                                     &
                  tdims%j_start : tdims%j_end )
              ! Droplet number from autoconversion

! Water tracer working arrays
type(mp_wtrac_type), intent(in) :: wtrac_mp(n_wtrac)

! Compressed super arrays ---------<------
real(kind=real_lsprec), intent(out) :: spr(ip, lsp_num_real)
logical, intent(out) :: spl(ip, lsp_num_log)

! Water tracer compressed structure
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

! Local variables

integer :: multrhc
              ! Zero if (rhc_row_length*rhc_rows) le 1, else 1

integer :: i,ii,ij,                                                            &
              ! Loop counters: I - horizontal field index;
           irhi,irhj
              !  IRHI,IRHJ-indices for RHcrit.

real(kind=jprb)             :: zhook_handle
character(len=*), parameter :: RoutineName='LS_PPNC_GATHER'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Multiple RHC rows?

if ( (rhc_row_length * rhc_rows)  >   1) then
  multrhc = 1
else
  multrhc = 0
end if

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do i=1, ip
  spr(i,tnuc_new_i) = zero
  spr(i,qcf2_i)   = zero
  spr(i,qrain_i)  = zero
  spr(i,qgraup_i) = zero
  spr(i,aero_i)   = zero

  ii   =  ix(i+i1-1,1)
  ij   =  ix(i+i1-1,2)
  irhi = (multrhc * (ii - 1)) + 1
  irhj = (multrhc * (ij - 1)) + 1

  spr(i,lsrain_i)          = real(lsrain(ii,ij), kind=real_lsprec)
  spr(i,lssnow_i)          = real(lssnow(ii,ij), kind=real_lsprec)
  spr(i,lssnow2_i)         = real(lssnow2(ii,ij), kind=real_lsprec)
  spr(i,lsgraup_i)         = real(lsgraup(ii,ij), kind=real_lsprec)
  spr(i,droplet_flux_i)    = real(droplet_flux(ii,ij), kind=real_lsprec)
  spr(i,cttemp_i)          = real(cttemp(ii,ij), kind=real_lsprec)
  spr(i,rainfrac_i)        = real(rainfrac(ii,ij), kind=real_lsprec)
  spr(i,rainfrac_impr_i)   = real(rainfrac_impr(ii,ij), kind=real_lsprec)
  if ( l_mcr_precfrac ) then
    spr(i,precfrac_k_i)      = real(precfrac_k(ii,ij), kind=real_lsprec)
    spr(i,precfrac_fall_i)   = real(precfrac_fall(ii,ij), kind=real_lsprec)
  else
    spr(i,precfrac_k_i)      = zero
    spr(i,precfrac_fall_i)   = zero
  end if
  spr(i,frac_ice_above_i)  = real(frac_ice_above(ii,ij), kind=real_lsprec)
  spr(i,p_i)               = real(p_theta_levels(ii,ij), kind=real_lsprec)
  spr(i,deltaz_i)          = real(deltaz(ii,ij), kind=real_lsprec)
  spr(i,uk_i)              = real(uk(ii,ij), kind=real_lsprec)
  spr(i,vk_i)              = real(vk(ii,ij), kind=real_lsprec)
  spr(i,ukp1_i)            = real(ukp1(ii,ij), kind=real_lsprec)
  spr(i,vkp1_i)            = real(vkp1(ii,ij), kind=real_lsprec)
  spr(i,r_theta_levels_i)  = real(r_theta_levels(ii,ij,onelevel),              &
                                  kind=real_lsprec)
  spr(i,r_theta_surf_i)    = real(r_theta_levels(ii,ij,0), kind=real_lsprec)
  spr(i,g_cos_theta_latitude_i)  = real(cos_theta_latitude(ii,ij),             &
                                        kind=real_lsprec)
  spr(i,rhodz_dry_i)   = real(rhodz_dry(ii,ij), kind=real_lsprec)
  spr(i,rhodz_moist_i) = real(rhodz_moist(ii,ij), kind=real_lsprec)
  spl(i,bland_i)       = bland(ii,ij) !This is a logical
  spr(i,cf_i)          = real(cf(ii,ij), kind=real_lsprec)
  spr(i,cfl_i)         = real(cfl(ii,ij), kind=real_lsprec)
  spr(i,cff_i)         = real(cff(ii,ij), kind=real_lsprec)
  spr(i,qcf_i)         = real(qcf(ii,ij), kind=real_lsprec)
  spr(i,qcl_i)         = real(qcl(ii,ij), kind=real_lsprec)
  if (l_progn_tnuc) then
    spr(i,tnuc_new_i)= real(tnuc_new(ii,ij), kind=real_lsprec)
  else
    spr(i,tnuc_new_i) = zero
  end if
  spr(i,q_i)           = real(q(ii,ij), kind=real_lsprec)
  spr(i,t_i)           = real(t(ii,ij), kind=real_lsprec)
  spr(i,n_drop_tpr_i)  = real(n_drop_tpr(ii, ij), kind=real_lsprec)

  if (l_mcr_qcf2) then
    spr(i,qcf2_i)   = real(qcf2(ii,ij), kind=real_lsprec)
  else
    spr(i,qcf2_i) = zero
  end if

  if (l_mcr_qrain) then
    spr(i,qrain_i)  = real(qrain(ii,ij), kind=real_lsprec)
  else
    spr(i,qrain_i) = zero
  end if


  if (l_mcr_qgraup) then
    spr(i,qgraup_i) = real(qgraup(ii,ij), kind=real_lsprec)
  else
    spr(i,qgraup_i) = zero
  end if

  if (l_murk_advect) then
    spr(i,aero_i)   = real(aerosol(ii,ij), kind=real_lsprec)
  else
    spr(i,aero_i) = zero
  end if

  !--------------------------------------------------------------------
  ! Pass through transfer diagnostics. Note that the 3D arrays are all
  ! initialized to zero in microphys_ctl, so the 1D versions below
  ! should just pick up zeros on the first microphysics iteration and
  ! if there are any subsequent iterations, then they should pick up
  ! the value of the diagnostic from previous iterations. Increments
  ! to each diagnostic are then added in the various transfer routines
  !--------------------------------------------------------------------

  if (l_psdep_diag) then
    spr(i,psdep_i) = real(psdep(ii,ij, level), kind=real_lsprec)
  else
    spr(i,psdep_i) = zero
  end if
  if (l_psaut_diag) then
    spr(i,psaut_i) = real(psaut(ii,ij, level), kind=real_lsprec)
  else
    spr(i,psaut_i) = zero
  end if
  if (l_psacw_diag) then
    spr(i,psacw_i) = real(psacw(ii,ij, level), kind=real_lsprec)
  else
    spr(i,psacw_i) = zero
  end if
  if (l_psacr_diag) then
    spr(i,psacr_i) = real(psacr(ii,ij, level), kind=real_lsprec)
  else
    spr(i,psacr_i) = zero
  end if
  if (l_psaci_diag) then
    spr(i,psaci_i) = real(psaci(ii,ij, level), kind=real_lsprec)
  else
    spr(i,psaci_i) = zero
  end if
  if (l_psmlt_diag) then
    spr(i,psmlt_i) = real(psmlt(ii,ij, level), kind=real_lsprec)
  else
    spr(i,psmlt_i) = zero
  end if

  if (l_psmltevp_diag) then
    spr(i,psmltevp_i) = real(psmltevp(ii,ij, level), kind=real_lsprec)
  else
    spr(i,psmltevp_i) = zero
  end if
  if (l_praut_diag) then
    spr(i,praut_i)    = real(praut(ii,ij, level), kind=real_lsprec)
  else
    spr(i,praut_i) = zero
  end if
  if (l_pracw_diag) then
    spr(i,pracw_i)    = real(pracw(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pracw_i) = zero
  end if
  if (l_prevp_diag) then
    spr(i,prevp_i)    = real(prevp(ii,ij, level), kind=real_lsprec)
  else
    spr(i,prevp_i) = zero
  end if
  if (l_pgaut_diag) then
    spr(i,pgaut_i)    = real(pgaut(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pgaut_i) = zero
  end if
  if (l_pgacw_diag) then
    spr(i,pgacw_i)    = real(pgacw(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pgacw_i) = zero
  end if
  if (l_pgacs_diag) then
    spr(i,pgacs_i)    = real(pgacs(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pgacs_i) = zero
  end if
  if (l_pgmlt_diag) then
    spr(i,pgmlt_i)    = real(pgmlt(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pgmlt_i) = zero
  end if
  if (l_pifrw_diag) then
    spr(i,pifrw_i)    = real(pifrw(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pifrw_i) = zero
  end if
  if (l_piprm_diag) then
    spr(i,piprm_i)    = real(piprm(ii,ij, level), kind=real_lsprec)
  else
    spr(i,piprm_i) = zero
  end if
  if (l_piprr_diag) then
    spr(i,piprr_i)    = real(piprr(ii,ij, level), kind=real_lsprec)
  else
    spr(i,piprr_i) = zero
  end if
  if (l_pidep_diag) then
    spr(i,pidep_i)    = real(pidep(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pidep_i) = zero
  end if
  if (l_piacw_diag) then
    spr(i,piacw_i)    = real(piacw(ii,ij, level), kind=real_lsprec)
  else
    spr(i,piacw_i) = zero
  end if
  if (l_piacr_diag) then
    spr(i,piacr_i)    = real(piacr(ii,ij, level), kind=real_lsprec)
  else
    spr(i,piacr_i)  = zero
  end if
  if (l_pimlt_diag) then
    spr(i,pimlt_i)    = real(pimlt(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pimlt_i)  = zero
  end if

  if (l_pimltevp_diag) then
    spr(i,pimltevp_i) = real(pimltevp(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pimltevp_i) = zero
  end if
  if (l_pifall_diag) then
    spr(i,pifall_i)   = real(pifall(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pifall_i)  = zero
  end if
  if (l_psfall_diag) then
    spr(i,psfall_i)   = real(psfall(ii,ij, level), kind=real_lsprec)
  else
    spr(i,psfall_i) = zero
  end if
  if (l_prfall_diag) then
    spr(i,prfall_i)   = real(prfall(ii,ij, level), kind=real_lsprec)
  else
    spr(i,prfall_i)  = zero
  end if
  if (l_pgfall_diag) then
    spr(i,pgfall_i)   = real(pgfall(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pgfall_i) = zero
  end if
  if (l_plset_diag) then
    spr(i,plset_i)    = real(plset(ii,ij, level), kind=real_lsprec)
  else
    spr(i,plset_i) = zero
  end if
  if (l_plevpset_diag) then
    spr(i,plevpset_i) = real(plevpset(ii,ij, level), kind=real_lsprec)
  else
    spr(i,plevpset_i) = zero
  end if
  if (l_pifrr_diag) then
    spr(i,pifrr_i)    = real(pifrr(ii,ij, level), kind=real_lsprec)
  else
    spr(i,pifrr_i) = zero
  end if

  if (l_sfwater_diag) then
    spr(i,sfwater_i) = real(sfwater(ii,ij, level), kind=real_lsprec)
  else
    spr(i,sfwater_i) = zero
  end if
  if (l_sfrain_diag) then
    spr(i,sfrain_i) = real(sfrain(ii,ij, level), kind=real_lsprec)
  else
    spr(i,sfrain_i) = zero
  end if
  if (l_sfsnow_diag) then
    spr(i,sfsnow_i) = real(sfsnow(ii,ij, level), kind=real_lsprec)
  else
    spr(i,sfsnow_i) = zero
  end if

  spr(i,n_drop_out_i)  = real(n_drop_out(ii,ij), kind=real_lsprec)
  spr(i,vfall_i)       = real(vfall(ii,ij), kind=real_lsprec)
  spr(i,vfall2_i)      = real(vfall2(ii,ij), kind=real_lsprec)
  spr(i,vfall_rain_i)  = real(vfall_rain(ii,ij), kind=real_lsprec)
  spr(i,vfall_graup_i) = real(vfall_graup(ii,ij), kind=real_lsprec)
  spr(i,rhc_i)         = real(rhcrit(irhi,irhj), kind=real_lsprec)

  spr(i,hmteff_i)     = real(hmteff(ii,ij), kind=real_lsprec)
  spr(i,zb_i)         = real(zb(ii,ij), kind=real_lsprec)

  if (l_fsd_generator) then
    spr(i,f_arr1_i) = real(f_arr(1,ii,ij), kind=real_lsprec)
    spr(i,f_arr2_i) = real(f_arr(2,ii,ij), kind=real_lsprec)
    spr(i,f_arr3_i) = real(f_arr(3,ii,ij), kind=real_lsprec)
  else
    spr(i,f_arr1_i) = zero
    spr(i,f_arr2_i) = zero
    spr(i,f_arr3_i) = zero
  end if

end do

! Set up water tracer super array equivalent
if (l_wtrac) then
  call ls_ppnc_gather_wtrac(ip, i1, ix, level, wtrac_mp, wtrac_mp_cpr)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_ppnc_gather


!-----------------------------------------------------------------------
!  Scatter variables, i.e post-process the data
!
! Default:    variables are changed from real_lsprec to real_umphys
! Exceptions: diagnostics in mphys_diags_mod which are always real_lsprec
!             (note USEed in at module level)
!-----------------------------------------------------------------------

subroutine ls_ppnc_scatter( level, ix, ip, i1,                                 &
  lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                                 &
  cf,cfl,cff,                                                                  &
  qcf,qcl,t, qcf2,qrain,qgraup,                                                &
  n_drop_out,                                                                  &
  aerosol,                                                                     &
!--------------------------------------------------------
! Variables passed from layer above
!--------------------------------------------------------
  q,                                                                           &
  vfall, vfall2, vfall_rain, vfall_graup,                                      &
  frac_ice_above, cttemp, rainfrac, precfrac_k, precfrac_fall,                 &
!--------------------------------------------------------
! COSP flag
!--------------------------------------------------------
  l_cosp_lsp, wtrac_mp,                                                        &
  !Compressed super arrays
  spr, wtrac_mp_cpr)

implicit none

integer, intent(in) ::                                                         &
 level,                                                                        &
              ! level number
 ip,                                                                           &
              ! Number of points
 i1,                                                                           &
              ! Starting point
 ix (  tdims%i_len * tdims%j_len , 2 )
              ! Gather/scatter index

real(kind=real_umphys), intent(in out) ::                                      &
  cf (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                   &
              ! Cloud fraction.
  cfl(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                   &
              ! Cloud liquid fraction.
  cff(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
              ! Cloud ice fraction.

logical, intent(in) ::  l_cosp_lsp
              ! Switch for COSP LS diagnostics

real(kind=real_umphys), intent(in out) ::                                      &
  q(tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end),                                                &
              ! Specific humidity (kg water/kg air).
  qcf(tdims%i_start:tdims%i_end,                                               &
      tdims%j_start:tdims%j_end),                                              &
              ! Cloud ice (kg per kg air).
  qcl(tdims%i_start:tdims%i_end,                                               &
      tdims%j_start:tdims%j_end),                                              &
              ! Cloud liquid water (kg per kg air).
  qcf2(tdims%i_start:tdims%i_end,                                              &
       tdims%j_start:tdims%j_end),                                             &
              ! Cloud ice2 (kg per kg air).
  qrain(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end),                                            &
              ! Rain water (kg per kg air).
  qgraup(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Graupel water (kg per kg air).
  t(tdims%i_start:tdims%i_end,                                                 &
    tdims%j_start:tdims%j_end),                                                &
              ! Temperature (K).
  aerosol(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end),                                          &
              ! Aerosol (K).
  lsrain(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Surface rainfall rate (kg m^-2 s^-1).
  lssnow(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Surface snowfall rate (kg m^-2 s^-1).
  lssnow2(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end),                                          &
              ! Layer snowfall rate (kg m^-2 s^-1).
  lsgraup(tdims%i_start:tdims%i_end,                                           &
          tdims%j_start:tdims%j_end),                                          &
              ! Layer graupelfall rate (kg m^-2 s^-1)
  droplet_flux(tdims%i_start:tdims%i_end,                                      &
               tdims%j_start:tdims%j_end),                                     &
              ! Water droplet flux / kg m^-2 s^-1
  cttemp(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Ice cloud top temperature (K)
  rainfrac(tdims%i_start:tdims%i_end,                                          &
           tdims%j_start:tdims%j_end),                                         &
              ! Rain fraction.
  precfrac_k(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end),                                       &
              ! Precipitation fraction at theta-level k
  precfrac_fall(tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
              ! Falling precipitation fraction at rho-levels
  frac_ice_above(tdims%i_start:tdims%i_end,                                    &
                 tdims%j_start:tdims%j_end),                                   &
              ! Ice fraction from layer above
  vfall(tdims%i_start:tdims%i_end,                                             &
        tdims%j_start:tdims%j_end),                                            &
              ! Fall velocity of ice (m per
  vfall2(tdims%i_start:tdims%i_end,                                            &
         tdims%j_start:tdims%j_end),                                           &
              ! Fall vel. of rain (m/s)
  vfall_rain(tdims%i_start:tdims%i_end,                                        &
             tdims%j_start:tdims%j_end),                                       &
              ! Fall vel. of rain (m/s)
  vfall_graup(tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end)
              ! Fall vel. of graupel (m/s)

real(kind=real_umphys), intent(in out) ::                                      &
  n_drop_out( tdims%i_start : tdims%i_end,                                     &
                  tdims%j_start : tdims%j_end )
              ! Droplet number from autoconversion

type(mp_wtrac_type), intent(in out) :: wtrac_mp(n_wtrac)

! Compressed super arrays ---------------
real(kind=real_lsprec), intent(in) :: spr(ip, lsp_num_real)

! Water tracer compressed structure
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

! Local variables

integer :: i,ii,ij
              ! Loop counters: I - horizontal field index;


real(kind=jprb)             :: zhook_handle
character(len=*), parameter :: RoutineName='LS_PPNC_SCATTER'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! assumption that all pairs ii,ij generated for each i are distinct
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do i=1, ip

  ii         = ix(i+i1-1,1)
  ij         = ix(i+i1-1,2)

  t(ii,ij)   = real(spr(i,t_i), kind=real_umphys)
  q(ii,ij)   = real(spr(i,q_i), kind=real_umphys)
  qcf(ii,ij) = real(spr(i,qcf_i), kind=real_umphys)
  qcl(ii,ij) = real(spr(i,qcl_i), kind=real_umphys)

  if (l_mcr_qcf2) then
    qcf2(ii,ij) = real(spr(i,qcf2_i), kind=real_umphys)

  end if

  if (l_mcr_qrain) then
    qrain(ii,ij) = real(spr(i,qrain_i), kind=real_umphys)
  end if

  if (l_mcr_qgraup) then
    qgraup(ii,ij) = real(spr(i,qgraup_i), kind=real_umphys)
  end if

  if (l_murk_advect) then
    aerosol(ii,ij) = real(spr(i,aero_i), kind=real_umphys)
  end if

  n_drop_out(ii,ij)     = real(spr(i,n_drop_out_i), kind=real_umphys)
  lsrain(ii,ij)         = real(spr(i,lsrain_i), kind=real_umphys)
  lssnow(ii,ij)         = real(spr(i,lssnow_i), kind=real_umphys)
  lssnow2(ii,ij)        = real(spr(i,lssnow2_i), kind=real_umphys)
  lsgraup(ii,ij)        = real(spr(i,lsgraup_i), kind=real_umphys)
  droplet_flux(ii,ij)   = real(spr(i,droplet_flux_i), kind=real_umphys)
  cttemp(ii,ij)         = real(spr(i,cttemp_i), kind=real_umphys)
  rainfrac(ii,ij)       = real(spr(i,rainfrac_i), kind=real_umphys)
  if ( l_mcr_precfrac ) then
    precfrac_k(ii,ij)     = real(spr(i,precfrac_k_i), kind=real_umphys)
    precfrac_fall(ii,ij)  = real(spr(i,precfrac_fall_i), kind=real_umphys)
  end if
  frac_ice_above(ii,ij) = real(spr(i,frac_ice_above_i), kind=real_umphys)

  if (i_cld_vn == i_cld_pc2) then
    cff(ii,ij)  = real(spr(i,cff_i), kind=real_umphys)
    cfl(ii,ij)  = real(spr(i,cfl_i), kind=real_umphys)
    cf(ii,ij)   = real(spr(i,cf_i), kind=real_umphys)
  end if  ! i_cld_pc2

  vfall(ii,ij)       = real(spr(i,vfall_i), kind=real_umphys)
  vfall2(ii,ij)      = real(spr(i,vfall2_i), kind=real_umphys)
  vfall_rain(ii,ij)  = real(spr(i,vfall_rain_i), kind=real_umphys)
  vfall_graup(ii,ij) = real(spr(i,vfall_graup_i), kind=real_umphys)

  ! Only store process rates in array for diagnostic
  ! if a particular diagnostic is requested,
  ! otherwise overwriting will occur
  ! (space for the 3D array in microphys_ctl is only allocated
  ! if the diagnostic is active, to save memory)
  if (l_aggfr_diag .or. l_cosp_lsp) then
    frac_agg(ii,ij, level) = real(spr(i,frac_agg_i), kind=real_lsprec)
  end if
  if (l_psdep_diag) then
    psdep(ii,ij, level) = real(spr(i,psdep_i), kind=real_lsprec)
  end if
  if (l_psaut_diag) then
    psaut(ii,ij, level) = real(spr(i,psaut_i), kind=real_lsprec)
  end if
  if (l_psacw_diag) then
    psacw(ii,ij, level) = real(spr(i,psacw_i), kind=real_lsprec)
  end if
  if (l_psacr_diag) then
    psacr(ii,ij, level) = real(spr(i,psacr_i), kind=real_lsprec)
  end if
  if (l_psaci_diag) then
    psaci(ii,ij, level) = real(spr(i,psaci_i), kind=real_lsprec)
  end if
  if (l_psmlt_diag) then
    psmlt(ii,ij, level) = real(spr(i,psmlt_i), kind=real_lsprec)
  end if
  if (l_psmltevp_diag) then
    psmltevp(ii,ij, level) = real(spr(i,psmltevp_i), kind=real_lsprec)
  end if
  if (l_praut_diag) then
    praut(ii,ij, level) = real(spr(i,praut_i), kind=real_lsprec)
  end if
  if (l_pracw_diag) then
    pracw(ii,ij, level) = real(spr(i,pracw_i), kind=real_lsprec)
  end if
  if (l_prevp_diag) then
    prevp(ii,ij, level) = real(spr(i,prevp_i), kind=real_lsprec)
  end if
  if (l_pgaut_diag) then
    pgaut(ii,ij, level) = real(spr(i,pgaut_i), kind=real_lsprec)
  end if
  if (l_pgacw_diag) then
    pgacw(ii,ij, level) = real(spr(i,pgacw_i), kind=real_lsprec)
  end if
  if (l_pgacs_diag) then
    pgacs(ii,ij, level) = real(spr(i,pgacs_i), kind=real_lsprec)
  end if
  if (l_pgmlt_diag) then
    pgmlt(ii,ij, level) = real(spr(i,pgmlt_i), kind=real_lsprec)
  end if
  if (l_pifrw_diag) then
    pifrw(ii,ij, level) = real(spr(i,pifrw_i), kind=real_lsprec)
  end if
  if (l_piprm_diag) then
    piprm(ii,ij, level) = real(spr(i,piprm_i), kind=real_lsprec)
  end if
  if (l_piprr_diag) then
    piprr(ii,ij, level) = real(spr(i,piprr_i), kind=real_lsprec)
  end if
  if (l_pidep_diag) then
    pidep(ii,ij, level) = real(spr(i,pidep_i), kind=real_lsprec)
  end if
  if (l_piacw_diag) then
    piacw(ii,ij, level) = real(spr(i,piacw_i), kind=real_lsprec)
  end if
  if (l_piacr_diag) then
    piacr(ii,ij, level) = real(spr(i,piacr_i), kind=real_lsprec)
  end if
  if (l_pimlt_diag) then
    pimlt(ii,ij, level) = real(spr(i,pimlt_i), kind=real_lsprec)
  end if
  if (l_pimltevp_diag) then
    pimltevp(ii,ij, level) = real(spr(i,pimltevp_i), kind=real_lsprec)
  end if
  if (l_pifall_diag) then
    pifall(ii,ij, level) = real(spr(i,pifall_i), kind=real_lsprec)
  end if
  if (l_psfall_diag) then
    psfall(ii,ij, level) = real(spr(i,psfall_i), kind=real_lsprec)
  end if
  if (l_prfall_diag) then
    prfall(ii,ij, level) = real(spr(i,prfall_i), kind=real_lsprec)
  end if
  if (l_pgfall_diag) then
    pgfall(ii,ij, level) = real(spr(i,pgfall_i), kind=real_lsprec)
  end if
  if (l_plset_diag) then
    plset(ii,ij, level) = real(spr(i,plset_i), kind=real_lsprec)
  end if
  if (l_plevpset_diag) then
    plevpset(ii,ij, level) = real(spr(i,plevpset_i), kind=real_lsprec)
  end if
  if (l_pifrr_diag) then
    pifrr(ii, ij, level) = real(spr(i,pifrr_i), kind=real_lsprec)
  end if
  if (l_vtbranch_diag) then
    vtbranch_flag(ii, ij, level) = real(spr(i,vtbranch_flag_i),                &
                                        kind=real_lsprec)
  end if
  if (l_vm_cry_diag) then
    vm_cry(ii, ij, level) = real(spr(i,vm_cry_i), kind=real_lsprec)
  end if
  if (l_vm_agg_diag) then
    vm_agg(ii, ij, level) = real(spr(i,vm_agg_i), kind=real_lsprec)
  end if
  if (l_vm_used_diag) then
    vm_used(ii, ij, level) = real(spr(i,vm_used_i), kind=real_lsprec)
  end if

  ! Seeder feeder diagnostics
  if (l_sfwater_diag) then
    sfwater(ii,ij, level) = real(spr(i,sfwater_i), kind=real_lsprec)
  end if
  if (l_sfrain_diag) then
    sfrain(ii,ij, level) = real(spr(i,sfrain_i), kind=real_lsprec)
  end if
  if (l_sfsnow_diag) then
    sfsnow(ii,ij, level) = real(spr(i,sfsnow_i), kind=real_lsprec)
  end if

  ! Radar reflectivity diagnostics
  if (l_ref_diag) then
    dbz_tot(ii, ij, level) = real(spr(i,dbz_tot_i), kind=real_lsprec)
    dbz_g(ii, ij, level)   = real(spr(i,dbz_g_i), kind=real_lsprec)
    dbz_i(ii, ij, level)   = real(spr(i,dbz_i_i), kind=real_lsprec)
    dbz_i2(ii, ij, level)  = real(spr(i,dbz_i2_i), kind=real_lsprec)
    dbz_r(ii, ij, level)   = real(spr(i,dbz_r_i), kind=real_lsprec)
    dbz_l(ii, ij, level)   = real(spr(i,dbz_l_i), kind=real_lsprec)
  end if

end do

! Scatter water tracer fields and deallocate compressed arrays
if (l_wtrac) then
  call ls_ppnc_scatter_wtrac(ip, i1, ix, level, wtrac_mp_cpr, wtrac_mp)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ls_ppnc_scatter


end module ls_ppnc_mod
