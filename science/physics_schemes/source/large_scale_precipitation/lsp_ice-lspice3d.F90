! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Precipitation microphysics calculations.
! Subroutine Interface:
module lsp_ice_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_ICE_MOD'

contains

subroutine lsp_ice(                                                            &
  p,deltaz,rhodz_dry,rhodz_moist,                                              &
  points, rhcpt, hmteff,zb,                                                    &
  qcf,qcl,tnuc_new,q,qcf2,qrain,qgraup, n_drop_tpr, n_drop_out,                &
  rainrate, vf_rain, snow_agg, vf_agg,                                         &
  snow_cry, vf_cry, grauprate, vf_graup, droplet_flux,                         &
  frac_ice_above, frac_agg, cttemp,                                            &
  rainfrac, rainfrac_impr, precfrac_k, precfrac_fall,                          &
  t,cfkeep,cfliqkeep,cficekeep,bland,                                          &
  psdep,psaut,psacw,psacr,psaci,psmlt,psmltevp,                                &
  praut,pracw,prevp,                                                           &
  pgaut,pgacw,pgacs,pgmlt,                                                     &
  pifrw,pifrr,piprm,piprr,pidep,piacw,piacr,pimlt,pimltevp,                    &
  pifall,psfall,prfall,pgfall,plset,plevpset,                                  &
  dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_l, dbz_r,                                 &
  sfwater, sfrain, sfsnow,                                                     &
  niters_mp, uk, vk, ukp1, vkp1,                                               &
  r_theta_levels_c, fv_cos_theta_latitude_c,                                   &
  r_theta_surf_c,                                                              &
  f_arr1, f_arr2, f_arr3,                                                      &
  vm_cry, vm_agg, vtbranch_flag, vm_used, wtrac_mp_cpr                         &
      )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod,           only: lcrcp, lfrcp, timestep, m0, t_scaling,         &
                                qcf0, timestep_mp, zero, one, small_number

use um_types,             only: real_lsprec

  ! General atmosphere modules- logicals and integers
use gen_phys_inputs_mod,  only: l_mr_physics

! Temporary bug-fixes
use science_fixes_mod,    only: l_fix_mcr_frac_ice

! Microphysics modules
use mphys_inputs_mod,     only: l_psd, l_mcr_qcf2, l_mcr_qrain,                &
                                l_mcr_qgraup, l_shape_rime,                    &
                                l_orograin,l_orogrime,                         &
                                not_generic_size_dist, graupel_option,         &
                                gr_field_psd, gr_srcols,                       &
                                sediment_loc,                                  &
                                all_sed_start, fall_end, all_sed_end,          &
                                rain_sed_end, warm_sed_end, i_mcr_iter,        &
                                i_mcr_iter_none,                               &
                                l_mcr_precfrac, l_subgrid_graupel_frac,        &
                                l_proc_fluxes,                                 &
                                i_update_precfrac, i_homog_areas, i_sg_correl

use mphys_bypass_mod,     only: l_crystals, l_ref_diag

! Dr Hook Modules
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

use lsp_accretion_mod,    only: lsp_accretion
use lsp_autoc_mod,        only: lsp_autoc
use lsp_capture_mod,      only: lsp_capture
use lsp_collection_mod,   only: lsp_collection
use lsp_deposition_mod,   only: lsp_deposition
use lsp_evap_mod,         only: lsp_evap
use lsp_evap_snow_mod,    only: lsp_evap_snow
use lsp_fall_mod,         only: lsp_fall, lsp_fall_graupel, lsp_fall_ice,      &
                                lsp_fall_rain
use lsp_fall_precfrac_mod,only: lsp_fall_precfrac
use lsp_combine_precfrac_mod,   only: lsp_combine_precfrac
use lsp_update_precfrac_mod,    only: lsp_update_precfrac
use lsp_graup_autoc_mod,  only: lsp_graup_autoc
use lsp_init_mod,         only: lsp_init
use lsp_melting_mod,      only: lsp_melting
use lsp_nucleation_mod,   only: lsp_nucleation
use lsp_riming_mod,       only: lsp_riming, lsp_riming_sphere,                 &
                                lsp_riming_graupel
use lsp_settle_mod,       only: lsp_settle
use lsp_snow_autoc_mod,   only: lsp_snow_autoc
use lsp_subgrid_mod,      only: lsp_subgrid
use lsp_tidy_mod,         only: lsp_tidy
use mphys_reflec_mod,     only: mphys_reflec
use lsp_orogwater_mod,    only: lsp_orogwater

! Water tracers
use free_tracers_inputs_mod,    only: l_wtrac, n_wtrac
use wtrac_mphys_mod,            only: mp_cpr_wtrac_type, mp_cpr_old_wtrac_type
use lsp_init_wtrac_mod,         only: lsp_init_wtrac
use lsp_settle_wtrac_mod,       only: lsp_settle_wtrac
use lsp_sedim_eulexp_wtrac_mod, only: lsp_sedim_eulexp_wtrac
use lsp_nucleation_wtrac_mod,   only: lsp_nucleation_wtrac
use lsp_deposition_wtrac_mod,   only: lsp_deposition_wtrac
use lsp_gen_wtrac_mod,          only: lsp_gen_wtrac

implicit none

! Description:
!   Updates ice, liquid and vapour contents, temperature and
!   cloud fractions as a result of microphysical processes.

! Method:
!   Calculates transfers of water between vapour, ice/snow,
!   cloud liquid water, rain and graupel.

!   Processes included are:
!   - Fall of hydrometeor into and out of the layer (sedimentation)
!   - Homogeneous and heterogeneous nucleation of ice;
!   - Deposition and sublimation of ice/snow;
!   - Autoconversion of ice->snow, snow->graupel, liquid->rain
!   - Collection processes
!   - Melting and evaporation

!   This is described in Unified Model Documentation Paper 26.

!   There are a number of different options for prognostic
!   hydrometeor variables. Each is independent of the others.
!   - Second prognostic cloud ice variables
!      Active if l_mcr_qcf2=.True. (in CNTLATM namelist)
!      The code supports the use of a second cloud ice prognostic
!      variable so that both cloud ice aggregates (QCF/QCF_AGG)
!      and cloud ice pristine crystals (QCF2/QCF_CRY) can be
!      represented and advected separately.
!      If False, then there is a diagnostic split within each level
!      at each timestep into ice crystals and snow aggregates.
!   - Prognostic rain
!      Active if l_mcr_qrain=.True. (in CNTLATM namelist)
!      If False, then rain is treated diagnostically.
!   - Prognostic graupel
!      Active if l_mcr_qgraup=.True. (in CNTLATM namelist)
!      If False, then graupel is not represented.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

! Declarations:

! Subroutine arguments

integer, intent(in) ::                                                         &
  points,                                                                      &
!         Number of points to be processed.
  niters_mp
!         Number of iterations of microphysics


logical, intent(in) ::                                                         &
  bland(points)
!         Land/sea mask

real (kind=real_lsprec), intent(in) ::                                         &
  p(points),                                                                   &
!         Air pressure at this level (Pa).
    deltaz(points),                                                            &
                           ! Thickness of layer (m)
    rhodz_dry(points),                                                         &
                           ! Dry air density * layer thickness (kg m-2)
    rhodz_moist(points),                                                       &
                           ! Moist air density * layer thick.  (kg m-2)
    hmteff(points),                                                            &
!         Effective mountain height forcing ascent for seeder feeder
    zb(points),                                                                &
!         Blocked layer depth for seeder feeder
    rhcpt(points),                                                             &
!         Critical relative humidity of all points for cloud formation.
    n_drop_tpr(points),                                                        &
!         droplet number (/ m 3) for autoconversion and settling
    f_arr1(points),                                                            &
    f_arr2(points),                                                            &
    f_arr3(points)

real (kind=real_lsprec), intent(in out) ::                                     &
  cfliqkeep(points),                                                           &
!         Liquid cloud fraction in this layer (no units).
    cficekeep(points),                                                         &
!         Frozen cloud fraction in this layer (no units).
    cfkeep(points),                                                            &
!         Total cloud fraction in this layer (no units).
    q(points),                                                                 &
!         Specific humidity at this level (kg water per kg air).
    qcf(points),                                                               &
!         Cloud ice (kg water per kg air).
    qcl(points),                                                               &
!         Cloud liquid water (kg water per kg air).
    qcf2(points),                                                              &
!         Second cloud ice (kg water per kg air).
    qrain(points),                                                             &
!         Rain (kg water per kg air).
    qgraup(points),                                                            &
!         Graupel (kg water per kg air).
    t(points)
!         Temperature at this level (K).

!     Hydrometeor flux between layers (kg m-2 s-1)
!     On input: Hydrometeor flux entering this layer from above.
!     On output: Hydrometeor flux leaving this layer.
!     Note: If only one ice prognostic is active (l_mcr_qcf2=.F.)
!     then SNOW_AGG contains all the ice/snow and SNOW_CRY is zero.
real (kind=real_lsprec), intent(in) ::                                         &
    tnuc_new(points)
!tnuc as function of dust (deg cel)
real (kind=real_lsprec), intent(in out) ::                                     &
  snow_agg(points),                                                            &
                            ! snow aggregates
  snow_cry(points),                                                            &
                            ! ice crystals
  rainrate(points),                                                            &
                            ! rain
  grauprate(points),                                                           &
                            ! graupel
  droplet_flux(points)
                            ! droplets

!     Hydrometeor mass-weighted fall velocities (m s-1)
!     On input: Fall velocity of hydrometeor entering layer.
!     On Output: Fall velocity of hydrometeor leaving layer.
!     Note: If only one ice prognostic is active, then only
!     VF_AGG is used.

real (kind=real_lsprec), intent(in out) ::                                     &
  vf_agg(points),                                                              &
                         ! snow aggregates
  vf_cry(points),                                                              &
                         ! ice crystals
  vf_rain(points),                                                             &
                         ! rain
  vf_graup(points)   ! graupel

!     Cloud/precipitation sub-grid fraction variables

real (kind=real_lsprec), intent(in out) ::                                     &
  cttemp(points),                                                              &
!         Ice cloud top temperature (K)
    rainfrac(points),                                                          &
!         Rain fraction (no units)
    precfrac_k(points),                                                        &
!         Prognostic precipitation fraction at theta-level k
!         (fraction associated with the prognostic rain and graupel fields).
    precfrac_fall(points),                                                     &
!         Precipitation fraction associated with the fall-flux of precip
!         and graupel;
!         in:  fraction falling in from the level above.
!         out: fraction falling out to the level below.
    frac_ice_above(points),                                                    &
!         in:  Fraction of ice falling from the layer above (no units)
!         out: Fraction of ice falling down to the next level
    frac_agg(points)
!         Aggregate fraction

real (kind=real_lsprec), intent(in) :: rainfrac_impr(points)

! Microphysical process rate diagnostics
real (kind=real_lsprec), intent(in out) ::                                     &
  psdep(points),                                                               &
                     ! Deposition of vapour to snow aggregates
  psaut(points),                                                               &
                     ! Autoconversion of aggregates from crystals
  psacw(points),                                                               &
                     ! Accretion of liq. water by snow aggregates
  psacr(points),                                                               &
                     ! Collection of rain by snow aggregates
  psaci(points),                                                               &
                     ! Collection of ice crystals by aggregates
  psmlt(points),                                                               &
                     ! Melting of snow aggregates
  psmltevp(points)  ! Evaporation of melting aggregates
real (kind=real_lsprec), intent(in out) ::                                     &
  praut(points),                                                               &
                     ! Autoconversion of cloud drops to rain
  pracw(points),                                                               &
                     ! Accretion of liq. water by rain
  prevp(points)  ! Evaporation of rain
real (kind=real_lsprec), intent(in out) ::                                     &
  pgaut(points),                                                               &
                     ! Autoconversion of graupel from aggregates
  pgacw(points),                                                               &
                     ! Accretion of liq. water by graupel
  pgacs(points),                                                               &
                     ! Collection of snow aggregates by graupel
  pgmlt(points)  ! Melting of graupel
real (kind=real_lsprec), intent(in out) ::                                     &
  pifrw(points),                                                               &
                     ! Homogeneous freezing nucleation
  pifrr(points),                                                               &
                     ! Homogeneous freezing of rain
  piprm(points),                                                               &
                     ! Heterogeneous (primary) nucleation
  piprr(points),                                                               &
                     ! Heterogeneous nucleation of rain
  pidep(points),                                                               &
                     ! Deposition of vapour to ice crystals
  piacw(points),                                                               &
                     ! Accretion of liq. water by ice crystals
  piacr(points),                                                               &
                     ! Collection of rain by ice crystals
  pimlt(points),                                                               &
                     ! Melting of ice crystals
  pimltevp(points)  ! Evaporation of melting ice crystals
real (kind=real_lsprec), intent(in out) ::                                     &
 pifall(points),                                                               &
                     ! Sedimentation of ice crystals
 psfall(points),                                                               &
                     ! Sedimentation of aggregates
 prfall(points),                                                               &
                     ! Sedimentation of rain
 pgfall(points) ! Sedimentation of graupel
real (kind=real_lsprec), intent(in out) ::                                     &
 plset(points),                                                                &
                     ! Droplet settling of liquid water
 plevpset(points) ! Evaporated settled droplets

! Seeder feeder water and precip rate diagnostics
real(kind=real_lsprec), intent(in out) ::                                      &
 sfwater(points),                                                              &
    ! Subgrid orographic water mixing ratio
 sfrain(points),                                                               &
    ! Subgrid orographic rain creation rate
 sfsnow(points)
    ! Subgrid orographic snow creation rate

real (kind=real_lsprec), intent(in) ::                                         &
  uk  (points),                                                                &
                     ! U wind at level k
  vk  (points),                                                                &
                     ! V wind at level k
  ukp1(points),                                                                &
                     ! U wind at level k+1
  vkp1(points),                                                                &
                     ! V wind at level k+1
  r_theta_levels_c(points),                                                    &
                     ! Distance from centre of the Earth
  r_theta_surf_c(points),                                                      &
                     ! ...and near surface (k=0)
  fv_cos_theta_latitude_c(points)
                     ! Finite volume cosine of latitude.

real (kind=real_lsprec), intent(in out) :: n_drop_out(points)
                    ! output droplet number from autoconversion

! Variables for defining different fallspeed relations for crystals
! and aggregates when the generic psd is used and l_diff_icevt = .true.
real (kind=real_lsprec), intent(in out) ::                                     &
 vm_agg(points),                                                               &
    ! mass-weighted mean fallspeed using aggregate fallspeed relation
 vm_cry(points),                                                               &
    ! mass-weighted mean fallspeed using crystal fallspeed relation
 vtbranch_flag(points),                                                        &
    ! 0=crystal vt-D relation; 1=aggregate vt-D relation
 vm_used(points)
    ! mass-weigthed fallspeed used

real (kind=real_lsprec), intent(out) :: dbz_tot(points)
    ! Total reflectivity (dBZ)
real (kind=real_lsprec), intent(out) :: dbz_g(points)
    ! Graupel reflectivity (dBZ)
real (kind=real_lsprec), intent(out) :: dbz_i(points)
    ! Ice Agg. reflectivity (dBZ)
real (kind=real_lsprec), intent(out) :: dbz_i2(points)
    ! Ice Cry. reflectivity (dBZ)
real (kind=real_lsprec), intent(out) :: dbz_l(points)
    ! Cloud liquid reflectivity (dBZ)
real (kind=real_lsprec), intent(out) :: dbz_r(points)
    ! Rain reflectivity (dBZ)

! Water tracer fields
! (Note, water tracer code only written for
!   l_mcr_qcf2 = l_mcr_graup = l_crystals = F)
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

!  Local scalars and dynamic arrays
integer ::  i     !  Loop counter (horizontal field index).
integer ::  i_wt  !  Loop counter (water tracer index)

real (kind=real_lsprec) :: lsrcp

real (kind=real_lsprec) ::                                                     &
  qs(points),                                                                  &
!         Saturated sp humidity for (T,p) in layer (kg kg-1)
    qsl(points)
!         Saturated sp humidity for (T,p) in layer
!         wrt water at all temps (kg kg-1)

!     Cumulative fall out of hydrometeor within iterations (kg m-2 s-1)
real (kind=real_lsprec) ::                                                     &
  snowt_agg(points),                                                           &
  snowt_cry(points),                                                           &
  rainratet(points),                                                           &
  graupratet(points)

! Subgrid orographic water MR
real (kind=real_lsprec) ::                                                     &
        ql_orog(points),                                                       &
        cf_orog(points),                                                       &
        ql_orog1(points),                                                      &
        t1a(points),                                                           &
        t1b(points),                                                           &
        qrain1a(points),                                                       &
        qrain1b(points),                                                       &
        qsnow1a(points),                                                       &
        qsnow1b(points),                                                       &
        pracw1a(points),                                                       &
        pracw1b(points),                                                       &
        psacw1a(points),                                                       &
        psacw1b(points),                                                       &
        dqrain(points),                                                        &
        dqsnow(points),                                                        &
        rain_liq_orog(points),                                                 &
        rain_mix_orog(points),                                                 &
        area_mix_orog(points),                                                 &
        area_liq_orog(points),                                                 &
        altitude(points)

logical ::  l_enh_rime(points) !whether riming is enhanced


real (kind=real_lsprec) ::                                                     &
  rho(points),                                                                 &
!         Density of air in the layer (kg m-3).
    rhor(points),                                                              &
!         1.0/RHO to speed up calculations (kg-1 m3).
    esi(points),                                                               &
!         saturation vapour pressure (wrt ice below zero Celsius)(Pa)
    esw(points)
!         saturation vapour pressure (wrt water at all temperatures)(Pa)

real (kind=real_lsprec) ::                                                     &
  cfice(points),                                                               &
!         fraction of ice inferred for the microphysics (no units).
    cficei(points),                                                            &
!         inverse of CFICE (no units)
    rainfraci(points),                                                         &
!         inverse of rain fraction
    cfliq(points),                                                             &
!         liquid cloud fraction for the microphysics
    cf(points),                                                                &
!         total cloud fraction for the microphysics (no units)
    dhi(points),                                                               &
!         CFL limit (s m-1)
    dhir(points)
!         1.0/DHI (m s-1)

real (kind=real_lsprec) ::                                                     &
    corr(points),                                                              &
!         density correction for fall speed (no units)
    rocor(points)
!         density correction for fall speed (no units)

real (kind=real_lsprec) ::                                                     &
    corr2(points),                                                             &
!         Temperature correction of viscosity etc. (no units)
    tcg(points),                                                               &
!         Temperature Factor for aggregate size distribution (no units)
    tcgi(points),                                                              &
!         Inverse of TCG (no units)
    tcgc(points),                                                              &
!         Temperature Factor for crystal size distribution (no units)
    tcgci(points),                                                             &
!         Inverse of TCGC (no units)
    tcgg(points),                                                              &
!         Temperature Factor for graupel size distribution (no units)
    tcggi(points)
!         Inverse of TCGC (no units)

real (kind=real_lsprec) ::                                                     &
  area_liq(points),                                                            &
!         Liquid only area of gridbox (no units)
    area_mix(points),                                                          &
!         Mixed phase area of gridbox (no units)
    area_ice(points),                                                          &
!         Ice only area of gridbox (no units)
    area_clear(points),                                                        &
!         Cloud free area of gridbox (no units)
    areamix_over_cfliq(points),                                                &
!         area_mix / cfliq (no units) - for perturbation sensitivity.
    rain_liq(points),                                                          &
!         Overlap fraction of gridbox between rain/graupel and liquid cloud
    rain_mix(points),                                                          &
!         Overlap fraction of gridbox between rain/graupel and mixed phase cloud
    rain_ice(points),                                                          &
!         Overlap fraction of gridbox between rain/graupel and ice cloud
    rain_clear(points),                                                        &
!         Overlap fraction of gridbox between rain/graupel and no cloud
    q_ice(points),                                                             &
!         Vapour content in the ice only part of the grid box (kg kg-1)
    q_clear(points),                                                           &
!         Vapour content in the cloud free part of the grid box(kg kg-1)
    qcf_agg(points),                                                           &
!         QCF in the form of aggregates (kg kg-1)
    qcf_cry(points)
!         QCF in the form of crystals (kg kg-1)

! Total amount of ice (crystals+aggregates)
real (kind=real_lsprec) :: qcf_tot(points)

! Fraction of ice falling down to the next level
real (kind=real_lsprec) :: frac_ice_fall(points)

real (kind=real_lsprec) ::                                                     &
  lheat_correc_liq(points),                                                    &
!         Reduction factor in evaporation limits because of latent heat
    lheat_correc_ice(points),                                                  &
!         Reduction factor in evaporation limits because of latent heat
    q_ice_1(points), q_ice_2(points),                                          &
    area_ice_1(points), area_ice_2(points),                                    &
!         Subgrid splitting for deposition term
    qcft(points)   ! Holds total ice content qcf_cry + qcf_agg

! Cloud and rain fraction transfer rate diagnostics
real (kind=real_lsprec) ::                                                     &
 cf_transfer_diag(points),                                                     &
                                ! Dummy to receive cf increment
 cfl_transfer_diag(points),                                                    &
                                ! Dummy to receive cfl increments
 cff_transfer_diag(points),                                                    &
                                ! Dummy to receive cff increments
 rf_transfer_diag(points)   ! Dummy to receive rainfrac increment

real (kind=real_lsprec) :: one_over_tsi        ! 1.0/(timestep_mp*iterations)

! Graupel autoconversion process rates
real (kind=real_lsprec) :: graut_psdep(points),                                &
                            ! Snow deposition rate for graupel autoconversion
                           graut_psacw(points),                                &
                            ! Aggregate-Snow riming rate for graupel
                            ! autoconversion
                           psacw_prev(points)
                            ! Value of aggregate-snow riming rate from a
                            ! previous calculation (used to differentiate
                            ! between psacw values on this substep from those
                            ! on an earlier substep of the microphysics.

! Dummy variable for lsp_deposition and mphys_reflec argument lists
real (kind=real_lsprec) :: dummy(points)

! Precipitation mass increment within each sub-region of the precipitation
! fraction (bits that overlap liquid-cloud, ice-cloud, clear-sky).
! Used to update the prognostic precipitation fraction.
real (kind=real_lsprec) :: dqprec_liq(points)
real (kind=real_lsprec) :: dqprec_mix(points)
real (kind=real_lsprec) :: dqprec_ice(points)
real (kind=real_lsprec) :: dqprec_clear(points)

! Area and precip increment within region where rain or graupel is generated
! in air that didn't yet contain any rain or graupel (so outside the above
! 4 fractions)
real (kind=real_lsprec) :: rain_new(points)
real (kind=real_lsprec) :: dqprec_new(points)

! Work variables for rain-mass increment calculations from
! orographic seeder-feeder enhancement
real (kind=real_lsprec) :: qprec(points)
real (kind=real_lsprec) :: frac
real (kind=real_lsprec) :: tmp
logical :: l_update_dqprec

! Fraction of ice, rain and graupel passing through the current model-level
! that are not falling out during the current sub-step
real (kind=real_lsprec) :: cry_nofall(points)
real (kind=real_lsprec) :: agg_nofall(points)
real (kind=real_lsprec) :: rain_nofall(points)
real (kind=real_lsprec) :: graup_nofall(points)

! Variables for calls to lsp_collection
logical ::                                                                     &
 l_use_area,                                                                   &
!       Use the ice partition to calculate transfer rates rather than
!       assuming a uniform distribution through the gridbox
   l_no_t_check
!       Do not check that the temperature is below zero degrees Celsius

logical ::                                                                     &
 l_use_agg_vt(points)              ! For each point defines which
! branch of the fallspeed relation is used. .true.=aggregate params.
! Set in lsp_init

integer ::                                                                     &
 ice_type1, ice_type2
!       Category of ice involved in collision process:
!       0 - crystals; 1 - aggregates; 3 - graupel.

! Structure to store old normal water fields - used for water tracer
! calculations
type(mp_cpr_old_wtrac_type) :: wtrac_mp_cpr_old

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_ICE'

!- End of header

! ======================================================================
!       Initialize variables and calculate microphysical quantities
! ======================================================================

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call lsp_init(points, timestep,                                                &
              timestep_mp,                                                     &
              t, p, cttemp,                                                    &
              deltaz, rhodz_dry, rhodz_moist, rho, rhor,                       &
              dhi, dhir,                                                       &
              qcf, qcf2, qrain, qgraup, rainrate,                              &
              qcf_agg, qcf_cry,                                                &
              qcf_tot, frac_agg,                                               &
              qs, qsl, esi, esw,                                               &
              cf_transfer_diag, cfl_transfer_diag,                             &
              cff_transfer_diag, rf_transfer_diag,                             &
              snow_cry, snow_agg, snowt_cry, snowt_agg,                        &
              rainratet, graupratet,                                           &
              lheat_correc_liq, lheat_correc_ice,                              &
              corr, corr2, rocor,                                              &
              agg_nofall, cry_nofall, rain_nofall, graup_nofall,               &
              tcg, tcgi, tcgc, tcgci, tcgg, tcggi,                             &
              cficekeep,                                                       &
              vm_cry, vm_agg, l_use_agg_vt, vm_used,                           &
              graut_psdep, graut_psacw, psacw_prev,                            &
              rainfrac, precfrac_k )

if (l_wtrac) then
  call lsp_init_wtrac(points, timestep, rhodz_dry, rhodz_moist,                &
                      q, qcl, qcf, qrain, t, droplet_flux, wtrac_mp_cpr_old,   &
                      wtrac_mp_cpr)
end if

! Calculate one_over_tsi once here rather than multiple times later
one_over_tsi = one/(timestep_mp*niters_mp)

! Latent heat of sublimation / cp (K).
lsrcp=lcrcp+lfrcp

! ======================================================================
!       Start iterations
! ======================================================================

do i = 1, points
       !------------------------------------------------------------
       ! Set fallspeed branch flag for diagnostics
       !------------------------------------------------------------
  if ( l_use_agg_vt(i) ) then
    vtbranch_flag(i) = one   ! Fallspeed setting for aggregates
  else
    vtbranch_flag(i) = zero   ! Fallspeed setting for crystals
  end if

  ! -------------------------------------------------------------
  ! Check that ice cloud fraction is sensible.
  ! -------------------------------------------------------------

  cfice(i)  = max( cficekeep(i), 0.001_real_lsprec )
  cficei(i) = one / cfice(i)
  cf(i)     = cfkeep(i)
  cfliq(i)  = cfliqkeep(i)
  cf(i)     = min( max(cf(i),cfice(i)) ,(cfice(i)+cfliq(i)) )

       ! -------------------------------------------------------------
       ! Calculate overlaps of liquid, ice and rain fractions
       ! -------------------------------------------------------------

  area_ice(i)   = max(cf(i)-cfliq(i),zero)
  area_clear(i) = max(one-cf(i),zero)

       ! -------------------------------------------------------------
       ! Update ice cloud top temperature if no ice falling in
       ! -------------------------------------------------------------

  if (snow_cry(i)+snow_agg(i)  <=  zero) then
    cttemp(i) = t(i)
  end if

end do

if ( l_fix_mcr_frac_ice ) then
  ! The copy of ice-cloud fraction limited to be above 0.001 above
  ! is used to calculate the increments to cff in various parts
  ! of the code. But this causes problems because if the actual
  ! ice cloud fraction is smaller than 0.001, the increments can
  ! remove nearly all the ice cloud fraction whilst leaving a
  ! significant ice mass. The implied in-cloud ice qcf/cff then
  ! becomes stupidly large.
  ! Under this bug-fix switch, avoid this by resetting the
  ! prognostic ice fraction consistent with the copies
  ! used to compute the increments
  do i = 1, points
    if ( qcf_agg(i) > zero  .or. qcf_cry(i) > zero .or.                        &
         snow_agg(i) > zero .or. snow_cry(i) > zero ) then
      ! Only reset where either there is already ice mass present
      ! or ice is falling in from above
      cficekeep(i) = cfice(i)
      cfkeep(i)    = cf(i)
    end if
  end do
end if

! ======================================================================
!        Droplet settling (if at start of microphysics)
! ======================================================================

if ( sediment_loc == all_sed_start .or. sediment_loc == fall_end               &
     .or. sediment_loc == rain_sed_end ) then

  call lsp_settle(points, one_over_tsi,                                        &
                  q, qcl, t, droplet_flux, bland,                              &
                  cfliq, rho, rhor, corr2, lcrcp,                              &
                  dhi, dhir,                                                   &
                  n_drop_tpr,                                                  &
                  plset, plevpset                                              &
                 )

  ! Water tracers droplet settling
  if (l_wtrac) then
    call lsp_settle_wtrac(points, cfliq, rhor, dhi, q, qcl, t,                 &
                          droplet_flux, wtrac_mp_cpr_old, wtrac_mp_cpr)
  end if

end if ! sediment_loc

! ======================================================================
!  Sedimentation of ice, rain and graupel (if at start of microphysics)
! ======================================================================

if ( sediment_loc == all_sed_start ) then

  ! This is the historic sedimentation location used by the UM microphysics

  ! Fall of prognostic precipitation fraction
  if ( l_mcr_precfrac ) then
    ! Need to call lsp_fall_precfrac just before the rain fall is done.
    ! Note: if using sub-grid fraction for graupel as well, then
    ! the graupel and rain fall need to be done together, so the
    ! options to do them in different places can't be used.
    call lsp_fall_precfrac( points, qrain, qgraup, rainrate, grauprate,        &
                            dhi, rhor, precfrac_k, precfrac_fall )

  end if

  call lsp_fall(points,                                                        &
                qcf_cry, qcf_agg, frac_agg, qrain, qgraup, t,                  &
                snow_agg, snow_cry, rainrate, grauprate,                       &
                snowt_agg, snowt_cry, rainratet, graupratet,                   &
                vf_agg, vf_cry, vf_rain, vf_graup,                             &
                area_clear, area_ice, cfice, cficei,                           &
                frac_ice_above, frac_ice_fall, cfkeep, cfliqkeep, cficekeep,   &
                rho, rhor, tcgi, tcgci,                                        &
                corr, dhi, dhir, rainfrac,                                     &
                pifall, psfall, prfall, pgfall,                                &
                one_over_tsi,                                                  &
                cf_transfer_diag, cff_transfer_diag,                           &
                uk, vk, ukp1, vkp1,                                            &
                r_theta_levels_c, fv_cos_theta_latitude_c,                     &
                l_use_agg_vt,                                                  &
                vm_used                                                        &
               )

  if (l_wtrac) then
    call lsp_sedim_eulexp_wtrac(points, rho, rhor, dhi, dhir, t,               &
                       vf_agg, qcf_agg, 'ice', wtrac_mp_cpr_old, wtrac_mp_cpr)

    if (l_mcr_qrain) then
      call lsp_sedim_eulexp_wtrac(points, rho, rhor, dhi, dhir, t,             &
                       vf_rain, qrain, 'rain', wtrac_mp_cpr_old, wtrac_mp_cpr)

    end if   ! l_mcr_qrain

  end if  ! l_wtrac

else if (sediment_loc == rain_sed_end .or. sediment_loc == warm_sed_end ) then

   ! Keep ice species sedimentation here at the start of the microphysics,
   ! but move the sedimentation of the rain part of lsp_fall to later on.

   ! Keep terms in the same order that they were in in lsp_fall, so ice
   ! first.

  call lsp_fall_ice(points, rho, rhor, t, qcf_agg, qcf_cry, cfice,             &
                    cficei, cfliqkeep, vf_agg, vf_cry, corr, tcgi, tcgci,      &
                    snow_agg, snow_cry, vm_used, dhi, dhir, frac_agg,          &
                    snowt_agg, snowt_cry, psfall, pifall,                      &
                    frac_ice_above, frac_ice_fall,                             &
                    cf_transfer_diag, cff_transfer_diag,                       &
                    cfkeep, cficekeep, area_clear, area_ice,                   &
                    ukp1, uk, vkp1, vk, one_over_tsi, r_theta_levels_c,        &
                    fv_cos_theta_latitude_c, l_use_agg_vt)

  ! Water tracers ice fall
  if (l_wtrac) then
    call lsp_sedim_eulexp_wtrac(points, rho, rhor, dhi, dhir, t,               &
                       vf_agg, qcf_agg, 'ice', wtrac_mp_cpr_old, wtrac_mp_cpr)
  end if

  ! Next, perform the graupel sedimentation, if prognostic graupel
  ! is in use.

  if ( l_mcr_qgraup ) then
    call lsp_fall_graupel(points, qgraup, rho, rhor, corr, dhi, dhir,          &
                          grauprate, graupratet, vf_graup, one_over_tsi,       &
                          rainfrac, pgfall )
  end if ! l_mcr_qgraup

end if ! sediment_loc

if ( l_proc_fluxes ) then
  ! Add falling out precip fluxes back onto the prognostics,
  ! but save the fraction that is not falling out
  ! (note: these fractions are initialised to 1.0 in lsp_init, so
  ! crucially they get left as 1.0 if l_proc_fluxes==.false., or if there
  ! is no fall-out flux).

  ! Ice crystals
  do i = 1, points
    tmp = snowt_cry(i) * dhi(i) * rhor(i)
    if ( tmp > zero ) then
      cry_nofall(i) = qcf_cry(i) / ( qcf_cry(i) + tmp )
      qcf_cry(i) = qcf_cry(i) + tmp
    end if
  end do
  ! Ice aggregates
  do i = 1, points
    tmp = snowt_agg(i) * dhi(i) * rhor(i)
    if ( tmp > zero ) then
      agg_nofall(i) = qcf_agg(i) / ( qcf_agg(i) + tmp )
      qcf_agg(i) = qcf_agg(i) + tmp
    end if
  end do
  ! Rain
  if ( l_mcr_qrain ) then
    do i = 1, points
      tmp = rainratet(i) * dhi(i) * rhor(i)
      if ( tmp > zero ) then
        rain_nofall(i) = qrain(i) / ( qrain(i) + tmp )
        qrain(i) = qrain(i) + tmp
      end if
    end do
  end if
  ! Graupel
  if ( l_mcr_qgraup ) then
    do i = 1, points
      tmp = graupratet(i) * dhi(i) * rhor(i)
      if ( tmp > zero ) then
        graup_nofall(i) = qgraup(i) / ( qgraup(i) + tmp )
        qgraup(i) = qgraup(i) + tmp
      end if
    end do
  end if

end if  ! ( l_proc_fluxes )

! ======================================================================
!        Subgrid-scale set-up calculations and tidy-ups
! ======================================================================
call lsp_subgrid(points,                                                       &
                q, qcf_cry, qcf_agg, qcf_tot, t,                               &
                qsl, qs, snow_cry, snow_agg, cry_nofall, agg_nofall,           &
                q_ice, q_clear, q_ice_1, q_ice_2,                              &
                area_liq,area_mix,area_ice,area_clear,                         &
                area_ice_1, area_ice_2,                                        &
                areamix_over_cfliq,                                            &
                rain_liq,rain_mix,rain_ice,rain_clear,                         &
                cf, cfliq, cfice, cficei,                                      &
                cfkeep, cficekeep, rainfrac, rainfraci, precfrac_k, rain_new,  &
                dqprec_liq, dqprec_mix, dqprec_ice, dqprec_clear, dqprec_new,  &
                lsrcp, rhcpt, wtrac_mp_cpr, wtrac_mp_cpr_old                   &
               )

! ======================================================================
!        HOMOGENEOUS (PIFRW) and HETEROGENEOUS NUCLEATION (PIPRM)
! ======================================================================
if (l_crystals) then
       ! Call nucleation with the crystals ice category (qcf_cry)
  call lsp_nucleation(points, timestep,                                        &
                q, qcl, tnuc_new, qrain, qcf_cry, qgraup, t,                   &
                qs, qsl,                                                       &
                cfliq,                                                         &
                area_liq, area_mix,                                            &
                cfkeep, cfliqkeep, cficekeep, rainfrac,                        &
                rain_liq, rain_mix, rain_ice, rain_clear, rainfraci,           &
                rho, rhor, lheat_correc_ice,                                   &
                corr, dhir, rain_nofall,                                       &
                lfrcp, lsrcp,                                                  &
                piprm, piprr, pifrw, pifrr, one_over_tsi,                      &
                cf_transfer_diag, cfl_transfer_diag,                           &
                cff_transfer_diag, rf_transfer_diag,                           &
                dqprec_liq, dqprec_mix, dqprec_ice, dqprec_clear,              &
                wtrac_mp_cpr_old                                               &
               )
else
       ! Call nucleation with the only ice category (qcf_agg)
  call lsp_nucleation(points, timestep,                                        &
                q, qcl, tnuc_new, qrain, qcf_agg, qgraup, t,                   &
                qs, qsl,                                                       &
                cfliq,                                                         &
                area_liq, area_mix,                                            &
                cfkeep, cfliqkeep, cficekeep, rainfrac,                        &
                rain_liq, rain_mix, rain_ice, rain_clear, rainfraci,           &
                rho, rhor, lheat_correc_ice,                                   &
                corr, dhir, rain_nofall,                                       &
                lfrcp, lsrcp,                                                  &
                piprm, piprr, pifrw, pifrr, one_over_tsi,                      &
                cf_transfer_diag, cfl_transfer_diag,                           &
                cff_transfer_diag, rf_transfer_diag,                           &
                dqprec_liq, dqprec_mix, dqprec_ice, dqprec_clear,              &
                wtrac_mp_cpr_old                                               &
               )

  if (l_wtrac) then
    call lsp_nucleation_wtrac(points, q, qcl, qcf_agg, qrain, t,               &
                              wtrac_mp_cpr_old, wtrac_mp_cpr)
  end if

end if  ! l_crystals

! ======================================================================
!             DEPOSITION/SUBLIMATION OF ICE CRYSTALS (PIDEP)
! ======================================================================
if (l_crystals) then

  do i = 1, points
       ! Calculate total ice content
    qcft(i) = qcf_cry(i) + qcf_agg(i)
  end do

  call lsp_deposition(points, timestep_mp,                                     &
                  q, qcl, qcf_cry, qcft, t, p,                                 &
                  q_ice_1, q_ice_2,                                            &
                  area_ice_1, area_ice_2,                                      &
                  esi, qs, qsl,                                                &
                  area_mix, cfliq, cfice, cficei,                              &
                  areamix_over_cfliq,                                          &
                  cfkeep, cfliqkeep, cficekeep,                                &
                  rho, tcgc, tcgci,                                            &
                  corr2, rocor, lheat_correc_ice, cry_nofall,                  &
                  lfrcp, lsrcp, 0,                                             &
                  not_generic_size_dist,                                       &
                  pidep, dummy, one_over_tsi,                                  &
                  cf_transfer_diag, cfl_transfer_diag,                         &
                  cff_transfer_diag,                                           &
                  l_use_agg_vt, wtrac_mp_cpr_old                               &
                 )

end if  ! l_crystals
! ======================================================================
!           DEPOSITION/SUBLIMATION OF SNOW AGGREGATES (PSDEP)
! ======================================================================
do i = 1, points
       ! Calculate total ice content
  qcft(i) = qcf_cry(i) + qcf_agg(i)
end do

call lsp_deposition(points, timestep_mp,                                       &
                q, qcl, qcf_agg, qcft, t, p,                                   &
                q_ice_1, q_ice_2,                                              &
                area_ice_1, area_ice_2,                                        &
                esi, qs, qsl,                                                  &
                area_mix, cfliq, cfice, cficei,                                &
                areamix_over_cfliq,                                            &
                cfkeep, cfliqkeep, cficekeep,                                  &
                rho, tcg, tcgi,                                                &
                corr2, rocor, lheat_correc_ice, agg_nofall,                    &
                lfrcp, lsrcp, 1,                                               &
                l_psd,                                                         &
                psdep, graut_psdep, one_over_tsi,                              &
                cf_transfer_diag, cfl_transfer_diag,                           &
                cff_transfer_diag,                                             &
                l_use_agg_vt, wtrac_mp_cpr_old                                 &
               )

if (l_wtrac) then
  call lsp_deposition_wtrac(points, q, qcl, qcf_agg, t, wtrac_mp_cpr_old,      &
                            wtrac_mp_cpr)
end if

if (l_mcr_qcf2) then  ! Note that l_mcr_qcf2==true implies
                           ! that l_crystals==true. It is *intended*
                           ! in parametrization development that
                           ! l_mcr_qcf2==true implies l_psd==false.
  ! ======================================================================
  !          AUTOCONVERSION OF ICE CRYSTALS to AGGREGATES (PSAUT)
  ! ======================================================================
  call lsp_snow_autoc(points, timestep_mp,                                     &
                qcf_cry, qcf_agg, t, cttemp,                                   &
                m0, t_scaling, cry_nofall, agg_nofall,                         &
                qcf0,                                                          &
                psaut, one_over_tsi                                            &
               )

  ! ======================================================================
  !           COLLECTION OF CRYSTALS BY AGGREGATES (PSACI)
  ! ======================================================================
  l_use_area=.false.
  l_no_t_check=.true.
  ice_type1=1  ! aggregates
  ice_type2=0  ! crystals
  call lsp_collection(points, timestep_mp,                                     &
              qcf_agg, qcf_cry, t,                                             &
              area_mix, area_ice, cficei, cficei,                              &
                    rho, rhor, m0, tcg, tcgi, tcgc, tcgci, corr,               &
                    agg_nofall, cry_nofall, ice_type1,ice_type2,               &
                    not_generic_size_dist,                                     &
                    l_use_area, l_no_t_check,                                  &
                    psaci, one_over_tsi,                                       &
                    l_use_agg_vt                                               &
                   )

end if  ! l_mcr_qcf2

! ======================================================================
!              RIMING OF ICE CRYSTALS BY CLOUD WATER (PIACW)
! ======================================================================
if (l_crystals) then
  if (l_shape_rime) then
    call lsp_riming(points, timestep_mp,                                       &
                    qcl, qcf_cry, t,                                           &
                    area_liq, area_mix, cfliq, cficei,                         &
                    rho, m0, tcgc, tcgci, corr, cry_nofall, lfrcp, 0,          &
                    not_generic_size_dist,                                     &
                    piacw, one_over_tsi,                                       &
                    l_use_agg_vt,                                              &
                    wtrac_mp_cpr_old                                           &
                         )
  else
    call lsp_riming_sphere(points, timestep_mp,                                &
                           qcl, qcf_cry, t,                                    &
                           area_liq, area_mix, cfliq, cficei,                  &
                           rho, m0, tcgc, tcgci, corr, cry_nofall, lfrcp, 0,   &
                           not_generic_size_dist,                              &
                           piacw, one_over_tsi,                                &
                           l_use_agg_vt                                        &
                          )
  end if  ! l_shape_rime
end if  ! l_crystals

! ======================================================================
!             RIMING OF SNOW AGGREGATES BY CLOUD WATER (PSACW)
! ======================================================================

if ( i_mcr_iter /= i_mcr_iter_none ) then
  ! There is more than one iteration of the large-scale precipitation scheme
  ! and therefore the value of pscaw may include contributions from a
  ! previous substep. Therefore, we need to work out the value of psacw
  ! and this point and save it, so we can determine the value on this
  ! substep of the microphysics.
  do i = 1, points
    psacw_prev(i) = psacw(i)
  end do
end if ! i_mcr_iter

! Enhance riming using orographic water if required
if (l_orograin .and. l_orogrime) then

  !   Initialise all arrays used for riming enhancement CHANGE
  do i = 1, points
    ql_orog(i) = zero
    cf_orog(i) = zero
    qsnow1a(i) = zero
    qsnow1b(i) = zero
    psacw1a(i) = zero
    psacw1b(i) = zero
    t1a(i) = zero
    t1b(i) = zero
    area_mix_orog(i) = zero
    area_liq_orog(i) = zero
    dqsnow(i) = zero
    altitude(i) = zero
    l_enh_rime(i) = .false.
  end do

  !   Calculate orographic water MR and cloud fraction
  call lsp_orogwater( points, hmteff, zb,                                      &
    r_theta_levels_c,r_theta_surf_c, fv_cos_theta_latitude_c,                  &
    p, t, q, qsl, esw, qcl,                                                    &
    ql_orog, cf_orog )

  !   Loop through points to initialize temp variables for riming calls
  do i = 1, points

    !     Save qrain and accr rate into temp variables for orog call
    qsnow1a(i) = qcf_agg(i)
    qsnow1b(i) = qcf_agg(i)
    t1a(i) = t(i)
    t1b(i) = t(i)
    psacw1a(i) = psacw(i)
    psacw1b(i) = psacw(i)
    !     Save ql_orog for accretion enhancement test
    ql_orog1(i) = ql_orog(i)
    !     Save into diagnostic for subgrid orographic water
    sfwater(i) = sfwater(i) + (ql_orog(i) / niters_mp)

    !     MINIMUM overlap of orographic cloud with ICE (as for resolved)
    area_mix_orog(i) = max(zero, (cf_orog(i) + cfice(i) - one) )
    area_liq_orog(i) = max(zero, (cf_orog(i) - area_mix_orog(i)) )

  end do ! Loop through all cloudy points


  if (l_shape_rime) then

    call lsp_riming(points, timestep_mp,                                       &
                    qcl, qcf_agg, t,                                           &
                    area_liq, area_mix, cfliq, cficei,                         &
                    rho, m0, tcg, tcgi, corr, agg_nofall, lfrcp, 1,            &
                    l_psd,                                                     &
                    psacw, one_over_tsi,                                       &
                    l_use_agg_vt,                                              &
                    wtrac_mp_cpr_old                                           &
                       )

    if (l_wtrac) then
      call lsp_gen_wtrac(points, qcl, qcf_agg, t, 'Rim',                       &
                         wtrac_mp_cpr_old%qcl, wtrac_mp_cpr_old%qcf,           &
                         wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,         &
                         wtrac_mp_cpr)
    end if

    call lsp_riming(points, timestep_mp,                                       &
                    ql_orog, qsnow1b, t1b,                                     &
                    area_liq_orog, area_mix_orog, cf_orog, cficei,             &
                    rho, m0, tcg, tcgi, corr, agg_nofall, lfrcp, 1,            &
                    l_psd,                                                     &
                    psacw1b, one_over_tsi,                                     &
                    l_use_agg_vt,                                              &
                    wtrac_mp_cpr_old                                           &
                       )

  else

    call lsp_riming_sphere(points, timestep_mp,                                &
                     qcl, qcf_agg, t,                                          &
                     area_liq, area_mix, cfliq, cficei,                        &
                     rho, m0, tcg, tcgi, corr, agg_nofall, lfrcp, 1,           &
                     l_psd,                                                    &
                     psacw, one_over_tsi,                                      &
                     l_use_agg_vt                                              &

                       )


    call lsp_riming_sphere(points, timestep_mp,                                &
                     ql_orog, qsnow1b, t1b,                                    &
                     area_liq_orog, area_mix_orog, cf_orog, cficei,            &
                     rho, m0, tcg, tcgi, corr, agg_nofall, lfrcp, 1,           &
                     l_psd,                                                    &
                     psacw1b, one_over_tsi,                                    &
                     l_use_agg_vt                                              &
                       )

  end if

  !   Apply riming enhancement to actual model variables
  do i = 1, points

    !     Amount of orog water accreted (ensure positive)
    dqsnow(i) = max(zero, qsnow1b(i) - qsnow1a(i) )

    if (dqsnow(i)>zero .and. ql_orog1(i)>zero) then

      l_enh_rime(i) = .true.
      qcf_agg(i) = qcf_agg(i) + dqsnow(i)
      q(i) = q(i) - dqsnow(i)
      if (l_wtrac) wtrac_mp_cpr_old%qchange(i) = dqsnow(i)

      !       Add LH for cond+freezing of rimed orog water
      t(i) = t(i) + (dqsnow(i) * lsrcp)

      !       Add mass transfer
      psacw(i) = psacw(i) + (psacw1b(i) - psacw1a(i))

      ! Seeder feeder snow production rate / kg kg-1 s-1
      sfsnow(i) = sfsnow(i) + dqsnow(i) * one_over_tsi

    end if

  end do ! Loop through all cloudy points

  if (l_wtrac) then
    call lsp_gen_wtrac(points, q, qcf_agg, t, 'Rim orog',                      &
                       wtrac_mp_cpr_old%q, wtrac_mp_cpr_old%qcf,               &
                       wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,           &
                       wtrac_mp_cpr)
  end if
else  ! seeder feeder off

  if (l_shape_rime) then
    call lsp_riming(points, timestep_mp,                                       &
                 qcl, qcf_agg, t,                                              &
                 area_liq, area_mix, cfliq, cficei,                            &
                 rho, m0, tcg, tcgi, corr, agg_nofall, lfrcp, 1,               &
                 l_psd,                                                        &
                 psacw, one_over_tsi,                                          &
                 l_use_agg_vt,                                                 &
                 wtrac_mp_cpr_old                                              &
                   )

    if (l_wtrac) then
      call lsp_gen_wtrac(points, qcl, qcf_agg, t, 'Rim',                       &
                         wtrac_mp_cpr_old%qcl, wtrac_mp_cpr_old%qcf,           &
                         wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,         &
                         wtrac_mp_cpr)
    end if

  else
    call lsp_riming_sphere(points, timestep_mp,                                &
                        qcl, qcf_agg, t,                                       &
                        area_liq, area_mix, cfliq, cficei,                     &
                        rho, m0, tcg, tcgi, corr, agg_nofall, lfrcp, 1,        &
                        l_psd,                                                 &
                        psacw, one_over_tsi,                                   &
                        l_use_agg_vt                                           &
                       )

  end if

end if ! l_orograin / l_orogrime

if ( i_mcr_iter == i_mcr_iter_none ) then
  ! No substepping/iterations of the microphysics scheme.
  ! Set graut_psacw to be the same as psacw
  do i = 1, points
    graut_psacw(i) = psacw(i)
  end do

else ! i_mcr_iter
  ! Substepping of the microphysics scheme is taking place.
  ! Need to subtract the value of psacw from the previous
  ! iteration, ensuring it is kept positive.
  do i = 1, points
    graut_psacw(i) = max(zero, (psacw(i) - psacw_prev(i)))
  end do

end if ! i_mcr_iter

  ! ======================================================================
  !             end OF RIMING OF SNOW AGGREGATES BY CLOUD WATER (PSACW)
  ! ======================================================================


if ( l_mcr_qgraup ) then
  ! ======================================================================
  !             AUTOCONVERSION OF SNOW to GRAUPEL  (PGAUT)
  ! ======================================================================

  call lsp_graup_autoc(points, timestep_mp,                                    &
                qrain, qcf_agg, qgraup, t, rho,                                &
                graut_psacw, graut_psdep, pgaut, agg_nofall, one_over_tsi,     &
                area_mix, area_mix_orog, l_enh_rime,                           &
                dqprec_mix, dqprec_ice, dqprec_new,                            &
                rain_mix, rain_ice, rain_new, precfrac_k )

  ! ======================================================================
  !              RIMING OF GRAUPEL BY CLOUD WATER (PGACW)
  ! ======================================================================
             ! Graupel is not included in ice cloud fraction
  if ( l_subgrid_graupel_frac ) then
    ! Riming accounting for sub-grid fraction of precip, including graupel
    ! Call wrapper routine which sets up the right area fraction variables
    ! to account for the graupel fraction (lsp_riming_sphere is called inside)
    call lsp_riming_graupel( points, timestep_mp, one_over_tsi,                &
                             qrain, qcl, qgraup, t, cfliq,                     &
                             rho, tcgg, tcggi, corr, pgacw,                    &
                             rainfraci, rain_liq, rain_mix, graup_nofall,      &
                             dqprec_liq, dqprec_mix, precfrac_k )
  else
    ! Riming with no sub-grid fraction of graupel
    ! (passing in whole mixed-phase cloud fraction; I think this
    !  is actually wrong, because if graupel has no sub-grid fraction
    !  it must overlap fully with the liquid cloud, so should pass in
    !  cfliq instead of area_mix.
    !  Also should pass in 1.0 instead of cficei)
    ! i.e. this is assuming that graupel exists only in the ice cloud fraction.
    ! Also this call wrongly passes in l_use_agg_vt for the ice-cloud,
    ! which has nothing to do with the graupel properties!
    call lsp_riming_sphere(points, timestep_mp,                                &
                           qcl, qgraup, t,                                     &
                           area_liq, area_mix, cfliq, cficei,                  &
                           rho, m0, tcgg, tcggi, corr, graup_nofall, lfrcp, 3, &
                           not_generic_size_dist,                              &
                           pgacw, one_over_tsi,                                &
                           l_use_agg_vt                                        &
                          )
  end if

  ! ======================================================================
  !              COLLECTION OF SNOW BY GRAUPEL (PGACS)
  ! ======================================================================
          ! For Graupel collecting snow the collision should result in a
          ! shattering of the ice and not a collection. For the moment,
          ! we shall not include this term when there is extra graupel
          ! production from snow-rain collisions. This should maintain
          ! bit-reproducibility with the present graupel scheme.
          !
          ! The break-up of ice aggregates to crystals will be dealt with
          ! at a later release.

  if (  graupel_option == gr_field_psd .or. graupel_option == gr_srcols ) then

    ! Set diagnostic to zero to show there is no transfer due to this process

    pgacs(:) = zero

  else

    ! Perform collection process as the standard UM graupel scheme

    ! Set flags input to lsp_collection,
    ! for collection of aggregates by graupel:
    l_no_t_check=.false.
    ice_type1=3  ! graupel
    ice_type2=1  ! aggregates


    if ( l_subgrid_graupel_frac ) then
      ! Using sub-grid graupel fraction

      l_use_area=.true.

      call lsp_collection(points, timestep_mp,                                 &
                          qgraup, qcf_agg, t,                                  &
                          rain_mix, rain_ice, rainfraci, cficei,               &
                          rho, rhor, m0, tcgg, tcggi,                          &
                          tcg, tcgi, corr, graup_nofall, agg_nofall,           &
                          ice_type1, ice_type2,                                &
                          l_psd,                                               &
                          l_use_area, l_no_t_check,                            &
                          pgacs, one_over_tsi,                                 &
                          l_use_agg_vt,                                        &
                          dqprec_mix=dqprec_mix, dqprec_ice=dqprec_ice,        &
                          qrain=qrain, precfrac_k=precfrac_k                   &
                         )

    else
      ! No sub-grid graupel fraction;
      ! note: what's actually done here is not consistent with the
      ! graupel being homogeneous over the grid-box as you'd expect;
      ! qgraup still gets scaled by cficei inside lsp_collection,
      ! so its really assuming graupel only exists inside the ice
      ! cloud fraction.  Setting of l_use_area=.false. just means
      ! the increment doesn't get scaled down by cfice, so it ends
      ! up inconsistent.

      l_use_area=.false.

      call lsp_collection(points, timestep_mp,                                 &
                          qgraup, qcf_agg, t,                                  &
                          area_mix, area_ice, cficei, cficei,                  &
                          rho, rhor, m0, tcgg, tcggi,                          &
                          tcg, tcgi, corr, graup_nofall, agg_nofall,           &
                          ice_type1, ice_type2,                                &
                          l_psd,                                               &
                          l_use_area, l_no_t_check,                            &
                          pgacs, one_over_tsi,                                 &
                          l_use_agg_vt                                         &
                         )

    end if  ! ( l_subgrid_graupel_frac )

  end if ! graupel_option

end if  ! l_mcr_graup

! ======================================================================
!               COLLECTION OF RAIN BY ICE CRYSTALS (PIACR)
! ======================================================================
if (l_crystals) then
  call lsp_capture(points, timestep_mp,                                        &
                   qcf_cry, qrain, qgraup, t, cficei,                          &
                   rainfrac,rain_liq,rain_mix,rain_ice,rain_clear,             &
                   rho, rhor, m0, tcgc, tcgci,                                 &
                   corr, dhir, cry_nofall, rain_nofall,                        &
                   lfrcp, 0,                                                   &
                   not_generic_size_dist,                                      &
                   piacr, one_over_tsi,                                        &
                   l_use_agg_vt,                                               &
                   dqprec_mix, dqprec_ice, precfrac_k,                         &
                   wtrac_mp_cpr_old                                            &
                        )
end if

! ======================================================================
!               COLLECTION OF RAIN BY SNOW AGGREGATES (PSACR)
! ======================================================================
call lsp_capture(points, timestep_mp,                                          &
                 qcf_agg, qrain, qgraup, t, cficei,                            &
                 rainfrac,rain_liq,rain_mix,rain_ice,rain_clear,               &
                 rho, rhor, m0, tcg, tcgi,                                     &
                 corr, dhir, agg_nofall, rain_nofall,                          &
                 lfrcp, 1,                                                     &
                 l_psd,                                                        &
                 psacr, one_over_tsi,                                          &
                 l_use_agg_vt,                                                 &
                 dqprec_mix, dqprec_ice, precfrac_k,                           &
                 wtrac_mp_cpr_old                                              &
                    )

if (l_wtrac) then
  call lsp_gen_wtrac(points, qrain, qcf_agg, t, 'Capture',                     &
                     wtrac_mp_cpr_old%qrain, wtrac_mp_cpr_old%qcf,             &
                     wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,             &
                     wtrac_mp_cpr)
end if

! ======================================================================
!                 EVAPORATION OF MELTING ICE CRYSTALS (PIMLTEVP)
! ======================================================================
if (l_crystals) then

  do i = 1, points
       ! Calculate total ice content
    qcft(i) = qcf_cry(i) + qcf_agg(i)
  end do

  call lsp_evap_snow(points, timestep_mp,                                      &
                     q, q_ice, qcf_cry, qcft, t, p, esw, qsl,                  &
                     area_ice, cficei, cfkeep, cficekeep,                      &
                     rho, tcgc, tcgci,                                         &
                     corr2, rocor, cry_nofall, lheat_correc_liq,               &
                     lsrcp, 0,                                                 &
                     not_generic_size_dist, pimltevp, one_over_tsi,            &
                     cf_transfer_diag, cff_transfer_diag,                      &
                     l_use_agg_vt,                                             &
                     wtrac_mp_cpr_old                                          &
                    )
end if  ! l_crystals

! ======================================================================
!                 EVAPORATION OF MELTING SNOW AGGREGATES (PSMLTEVP)
! ======================================================================
do i = 1, points
       ! Calculate total ice content
  qcft(i) = qcf_cry(i) + qcf_agg(i)
end do

call lsp_evap_snow(points, timestep_mp,                                        &
                   q, q_ice, qcf_agg, qcft, t, p, esw, qsl,                    &
                   area_ice, cficei, cfkeep, cficekeep,                        &
                   rho, tcg, tcgi,                                             &
                   corr2, rocor, agg_nofall, lheat_correc_liq, lsrcp, 1,       &
                   l_psd, psmltevp, one_over_tsi,                              &
                   cf_transfer_diag,cff_transfer_diag,                         &
                   l_use_agg_vt,                                               &
                   wtrac_mp_cpr_old                                            &
                  )

if (l_wtrac) then
  call lsp_gen_wtrac(points, qcf_agg, q, t, 'Evapsnow',                        &
                     wtrac_mp_cpr_old%qcf, wtrac_mp_cpr_old%q,                 &
                     wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,             &
                     wtrac_mp_cpr)
end if

     ! Perform the melting steps here

! ======================================================================
!                    MELTING OF ICE CRYSTALS (PIMLT)
! ======================================================================
if (l_crystals) then

  do i = 1, points
     ! Calculate total ice content
    qcft(i) = qcf_cry(i) + qcf_agg(i)
  end do

  call lsp_melting(points, timestep_mp,                                        &
                   q, q_ice, qgraup, qcf_cry, qcft, qrain, qsl, t, p,          &
                   area_liq, area_mix, area_ice, area_clear,                   &
                   cficei, frac_ice_above,                                     &
                   rainfrac, rain_liq, rain_mix,                               &
                   rain_ice, rain_clear, rain_new, cfkeep, cficekeep,          &
                   rho, rhor, m0, tcg, tcgi, corr2, rocor, cry_nofall,         &
                   lfrcp, 0,                                                   &
                   not_generic_size_dist,                                      &
                   pimlt, one_over_tsi,                                        &
                   cf_transfer_diag, cff_transfer_diag,                        &
                   rf_transfer_diag,                                           &
                   l_use_agg_vt,                                               &
                   dqprec_ice, dqprec_mix, dqprec_new, precfrac_k,             &
                   wtrac_mp_cpr_old                                            &
                  )
end if  ! l_crystals

! ======================================================================
!                    MELTING OF SNOW AGGREGATES (PSMLT)
! ======================================================================
do i = 1, points
     ! Calculate total ice content
  qcft(i) = qcf_cry(i) + qcf_agg(i)
end do

call lsp_melting(points, timestep_mp,                                          &
                 q, q_ice, qgraup, qcf_agg, qcft, qrain, qsl, t, p,            &
                 area_liq, area_mix, area_ice, area_clear,                     &
                 cficei, frac_ice_above,                                       &
                 rainfrac, rain_liq, rain_mix,                                 &
                 rain_ice, rain_clear, rain_new, cfkeep, cficekeep,            &
                 rho, rhor, m0, tcgc, tcgci,                                   &
                 corr2, rocor, agg_nofall,                                     &
                 lfrcp, 1,                                                     &
                 l_psd,                                                        &
                 psmlt, one_over_tsi,                                          &
                 cf_transfer_diag, cff_transfer_diag,                          &
                 rf_transfer_diag,                                             &
                 l_use_agg_vt,                                                 &
                 dqprec_ice, dqprec_mix, dqprec_new, precfrac_k,               &
                 wtrac_mp_cpr_old                                              &
                )

if (l_wtrac) then
  call lsp_gen_wtrac(points, qcf_agg, qrain, t, 'Melt',                        &
                     wtrac_mp_cpr_old%qcf, wtrac_mp_cpr_old%qrain,             &
                     wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,             &
                     wtrac_mp_cpr)
end if

! ======================================================================
!                    MELTING OF GRAUPEL (PGMLT)
! ======================================================================
if (l_mcr_qgraup) then
   ! Graupel does not update cloud fractions so there is no need
   ! to update qcft (it is not used)
  if ( l_subgrid_graupel_frac ) then
    ! Using sub-grid graupel fraction; pass in rainfrac (which has been
    ! set equal to the prognostic precipitation fraction,
    ! which includes graupel) instead of cfice
    call lsp_melting(points, timestep_mp,                                      &
                     q, q_ice, qgraup, qgraup, qcft, qrain, qsl, t, p,         &
                     area_liq, area_mix, area_ice, area_clear,                 &
                     rainfraci, frac_ice_above,                                &
                     rainfrac, rain_liq, rain_mix,                             &
                     rain_ice, rain_clear, rain_new, cfkeep, cficekeep,        &
                     rho, rhor, m0, tcgg, tcggi, corr2, rocor, graup_nofall,   &
                     lfrcp, 3,                                                 &
                     not_generic_size_dist,                                    &
                     pgmlt, one_over_tsi,                                      &
                     cf_transfer_diag, cff_transfer_diag,                      &
                     rf_transfer_diag,                                         &
                     l_use_agg_vt,                                             &
                     dqprec_ice, dqprec_mix, dqprec_new, precfrac_k,           &
                     wtrac_mp_cpr_old                                          &
                    )
  else
    ! No sub-grid graupel fraction; passing in cfice for fraction but
    ! the checking on ice_type=3 inside should prevent it from being used.
    call lsp_melting(points, timestep_mp,                                      &
                     q, q_ice, qgraup, qgraup, qcft, qrain, qsl, t, p,         &
                     area_liq, area_mix, area_ice, area_clear,                 &
                     cficei, frac_ice_above,                                   &
                     rainfrac, rain_liq, rain_mix,                             &
                     rain_ice, rain_clear, rain_new, cfkeep, cficekeep,        &
                     rho, rhor, m0, tcgg, tcggi, corr2, rocor, graup_nofall,   &
                     lfrcp, 3,                                                 &
                     not_generic_size_dist,                                    &
                     pgmlt, one_over_tsi,                                      &
                     cf_transfer_diag, cff_transfer_diag,                      &
                     rf_transfer_diag,                                         &
                     l_use_agg_vt,                                             &
                     dqprec_ice, dqprec_mix, dqprec_new, precfrac_k,           &
                     wtrac_mp_cpr_old                                          &
                    )
  end if  ! ( l_subgrid_graupel_frac )
end if  ! l_mcr_qgraup

! ======================================================================
!                   EVAPORATION OF RAINDROPS (PREVP)
! ======================================================================
call lsp_evap(points, timestep_mp, p, q, qrain, t,                             &
              qgraup, q_ice, q_clear,                                          &
              rainfrac, rain_liq, rain_mix,                                    &
              rain_ice, rain_clear,                                            &
              rho, corr, corr2, rocor,                                         &
              dhir, lcrcp, lheat_correc_liq,                                   &
              qsl, esw, rain_nofall,                                           &
              prevp, rf_transfer_diag, one_over_tsi,                           &
              dqprec_ice, dqprec_clear, precfrac_k,                            &
              wtrac_mp_cpr_old                                                 &
             )

if (l_wtrac) then
  call lsp_gen_wtrac(points, qrain, q, t, 'Evaprain',                          &
                     wtrac_mp_cpr_old%qrain, wtrac_mp_cpr_old%q,               &
                     wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,             &
                     wtrac_mp_cpr)
end if

! ======================================================================
!             ACCRETION OF CLOUD DROPLETS ON RAINDROPS (PRACW)
! ======================================================================
!  Enhanced by subgrid orographic motions
! ================================================================

if (l_orograin) then

  !   Calculate orographic water MR and cloud fraction
  call lsp_orogwater( points, hmteff, zb,                                      &
    r_theta_levels_c,r_theta_surf_c, fv_cos_theta_latitude_c,                  &
    p, t, q, qsl, esw, qcl,                                                    &
    ql_orog, cf_orog )

  !   Loop through points to initialize temp variables for accr calls
  do i = 1, points

    !     Save qrain into temp variables for orog call
    qrain1a(i) = qrain(i)
    qrain1b(i) = qrain(i)
    pracw1a(i) = pracw(i)
    pracw1b(i) = pracw(i)
    !     Save ql_orog for outputting afterwards
    ql_orog1(i) = ql_orog(i)

    !     Overlap with RAIN - assume max value (only use sum)
    rain_liq_orog(i) = max(zero, min( cf_orog(i), rainfrac(i) ) )
    rain_mix_orog(i) =  zero

  end do

  l_update_dqprec = .true.
  call lsp_accretion(points, timestep_mp, qgraup, qcl, qrain,                  &
                    cfliq, rainfrac, rain_liq, rain_mix,                       &
                    rho, corr, dhir, rain_nofall,                              &
                    pracw, one_over_tsi,                                       &
                    r_theta_levels_c, fv_cos_theta_latitude_c,                 &
                    f_arr1, f_arr2, f_arr3,                                    &
                    l_update_dqprec, dqprec_liq, dqprec_mix, precfrac_k,       &
                    wtrac_mp_cpr_old                                           &
                        )

  if (l_wtrac) then
    call lsp_gen_wtrac(points, qcl, qrain, t, 'Accretion',                     &
                       wtrac_mp_cpr_old%qcl, wtrac_mp_cpr_old%qrain,           &
                       wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,           &
                       wtrac_mp_cpr)
  end if

  l_update_dqprec = .false.
  call lsp_accretion(points, timestep_mp, qgraup, ql_orog, qrain1b,            &
                   cf_orog, rainfrac, rain_liq_orog, rain_mix_orog,            &
                   rho, corr, dhir, rain_nofall,                               &
                   pracw1b, one_over_tsi,                                      &
                   r_theta_levels_c, fv_cos_theta_latitude_c,                  &
                   f_arr1, f_arr2, f_arr3,                                     &
                   l_update_dqprec, dqprec_liq, dqprec_mix, precfrac_k,        &
                   wtrac_mp_cpr_old                                            &
                        )


  !   Apply accretion enhancement to actual model variables
  !   ...(but only if riming wasnt enhanced at that point).
  do i = 1, points

    !     Subgrid orog water MR actually accreted by rain
    dqrain(i) = max(zero, qrain1b(i) - qrain1a(i) )

    if (dqrain(i)>zero .and. ql_orog1(i)>zero .and.                            &
          .not. l_enh_rime(i) ) then

      !       Add orographic accretion to that from resolved cloud
      qrain(i) = qrain(i) + dqrain(i)
      q(i) = q(i) - dqrain(i)
      t(i) = t(i) + (dqrain(i) * lcrcp)
      if (l_wtrac) wtrac_mp_cpr_old%qchange(i) = dqrain(i)

      !       Add mass transfer
      pracw(i) = pracw(i) +  (pracw1b(i) - pracw1a(i))

      ! Seeder feeder rain production rate / kg kg-1 s-1
      sfrain(i) = sfrain(i) + dqrain(i) * one_over_tsi

    else
      dqrain(i) = zero
    end if

  end do

  if (l_wtrac) then
    call lsp_gen_wtrac(points, q, qrain, t, 'Accr orog',                       &
                       wtrac_mp_cpr_old%q, wtrac_mp_cpr_old%qrain,             &
                       wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,           &
                       wtrac_mp_cpr)
  end if

  if ( l_mcr_precfrac ) then
    ! If using prognostic precip fraction, update increments within
    ! the orographic liquid cloud region...
    if ( i_update_precfrac == i_homog_areas ) then
      do i = 1, points
        tmp = one / max( rain_liq_orog(i), small_number )

        ! Assign the increment between the various sub-partitions of
        ! the rain-fraction, assuming maximal overlap with first
        ! liquid-only cloud, then mixed-phase cloud, then ice-cloud...
        frac = rain_liq_orog(i)
        if ( frac > zero ) then
          dqprec_liq(i) = dqprec_liq(i) + dqrain(i)*min(rain_liq(i),frac) * tmp
          frac = frac - rain_liq(i)
        end if
        if ( frac > zero ) then
          dqprec_mix(i) = dqprec_mix(i) + dqrain(i)*min(rain_mix(i),frac) * tmp
          frac = frac - rain_mix(i)
        end if
        if ( frac > zero ) then
          dqprec_ice(i) = dqprec_ice(i) + dqrain(i)*min(rain_ice(i),frac) * tmp
          frac = frac - rain_ice(i)
        end if
        if ( frac > zero ) then
          dqprec_clear(i)=dqprec_clear(i)+dqrain(i)*min(rain_clear(i),frac)*tmp
          frac = frac - rain_clear(i)
        end if

        ! Note: we must have rain_liq_orog <= rainfrac, so can't end up with
        ! any of the increment dqrain leftover unasigned to the dqprec stores.

      end do
    else if ( i_update_precfrac == i_sg_correl ) then

      ! Set precip mass before the increment was applied
      ! (ignoring negative values)
      do i = 1, points
        qprec(i) = max( qrain(i) - dqrain(i), zero )
        dqrain(i) = max( qrain(i), zero ) - qprec(i)
      end do
      if ( l_subgrid_graupel_frac ) then
        do i = 1, points
          qprec(i) = qprec(i) + max( qgraup(i), zero )
        end do
      end if
      ! Calculate combined fraction of the existing precip and increment
      ! (precfrac_k is updated)
      call lsp_combine_precfrac( points,                                       &
                                 qprec, dqrain, precfrac_k, rain_liq_orog )

    end if  ! ( i_update_precfrac )
  end if  ! ( l_mcr_precfrac )

else  ! seeder feeder off

  l_update_dqprec = .true.
  call lsp_accretion(points, timestep_mp, qgraup, qcl, qrain,                  &
                     cfliq, rainfrac, rain_liq, rain_mix,                      &
                     rho, corr, dhir, rain_nofall,                             &
                     pracw, one_over_tsi,                                      &
                    r_theta_levels_c, fv_cos_theta_latitude_c,                 &
                    f_arr1, f_arr2, f_arr3,                                    &
                    l_update_dqprec, dqprec_liq, dqprec_mix, precfrac_k,       &
                    wtrac_mp_cpr_old                                           &
                        )

  if (l_wtrac) then
    call lsp_gen_wtrac(points, qcl, qrain, t, 'Accretion',                     &
                       wtrac_mp_cpr_old%qcl, wtrac_mp_cpr_old%qrain,           &
                       wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,           &
                       wtrac_mp_cpr)
  end if

end if


! ======================================================================
!              AUTOCONVERSION OF CLOUD LIQUID to RAIN (PRAUT)
! ======================================================================
call lsp_autoc(points, timestep_mp, qgraup, qcl, qrain, t, p,                  &
               cfliq, rhcpt,                                                   &
               area_liq, area_mix, area_ice, rainfrac,                         &
               rain_liq, rain_mix, rain_ice, rain_clear, rain_new,             &
               rho, rhor, corr2, lcrcp,                                        &
               one_over_tsi, praut, rf_transfer_diag,                          &
               n_drop_tpr, n_drop_out,                                         &
               r_theta_levels_c, fv_cos_theta_latitude_c,                      &
               f_arr1, f_arr2, f_arr3,                                         &
               dqprec_liq, dqprec_mix, dqprec_new, precfrac_k,                 &
               wtrac_mp_cpr_old                                                &
                 )

if (l_wtrac) then
  call lsp_gen_wtrac(points, qcl, qrain, t, 'Autoc',                           &
                     wtrac_mp_cpr_old%qcl, wtrac_mp_cpr_old%qrain,             &
                     wtrac_mp_cpr_old%t, wtrac_mp_cpr_old%qchange,             &
                     wtrac_mp_cpr)
end if

! ======================================================================
!             Update prognostic precip fraction
! ======================================================================
if ( l_mcr_precfrac ) then
  if ( i_update_precfrac == i_homog_areas) then
    ! We have now stored the increments to qrain (and optionally including
    ! qgraup as well) attributable to each of the sub-grid regions,
    ! in the 5 dqprec stores.
    ! We now calculate the new in-region precip mixing-ratio
    ! in each region.  Based on these, we compute a new effective precip
    ! fraction such that the in-region mean value of qp is the
    ! precip-mass-weighted mean over the regions, <qp*qp> / <qp>
    call lsp_update_precfrac( points,                                          &
                              rain_liq, rain_mix, rain_ice,                    &
                              rain_clear, rain_new,                            &
                              dqprec_liq, dqprec_mix, dqprec_ice,              &
                              dqprec_clear, dqprec_new,                        &
                              qrain, qgraup, precfrac_k )
  end if
  ! Reset rainfrac consistent with the updated precfrac
  do i = 1, points
    rainfrac(i) = precfrac_k(i)
  end do
end if

! ======================================================================
!             Droplet settling (if at end of microphysics)
! ======================================================================

if ( sediment_loc == all_sed_end .or. sediment_loc == warm_sed_end ) then

  ! In this case, droplet settling needs to occur at the end of the
  ! microphysics, but before lsp_fall to match the same order as they are
  ! performed above.

  call lsp_settle(points, one_over_tsi,                                        &
                  q, qcl, t, droplet_flux, bland,                              &
                  cfliq, rho, rhor, corr2, lcrcp,                              &
                  dhi, dhir,                                                   &
                  n_drop_tpr,                                                  &
                  plset, plevpset                                              &
                 )

  ! Water tracers droplet settling
  if (l_wtrac) then
    call lsp_settle_wtrac(points, cfliq, rhor, dhi, q, qcl, t,                 &
                          droplet_flux, wtrac_mp_cpr_old, wtrac_mp_cpr)
  end if

end if ! sediment_loc

! ======================================================================
!  Sedimentation of ice, rain and graupel (if at end of microphysics)
! ======================================================================

! Fall of prognostic precipitation fraction
if ( l_mcr_precfrac .and. ( .not. sediment_loc == all_sed_start ) ) then
  ! All options for sediment_loc other than all_sed_start have the rain
  ! fall here at the end.  Need to call lsp_fall_precfrac just before
  ! the rain fall is done.
  call lsp_fall_precfrac( points, qrain, qgraup, rainrate, grauprate,          &
                          dhi, rhor, precfrac_k, precfrac_fall )
end if

if ( sediment_loc == fall_end .or. sediment_loc == all_sed_end ) then

  ! Option to move fall of all main hydrometeor species to the end of the
  ! microphysics.

  call lsp_fall(points,                                                        &
                qcf_cry, qcf_agg, frac_agg, qrain, qgraup, t,                  &
                snow_agg, snow_cry, rainrate, grauprate,                       &
                snowt_agg, snowt_cry, rainratet, graupratet,                   &
                vf_agg, vf_cry, vf_rain, vf_graup,                             &
                area_clear, area_ice, cfice, cficei,                           &
                frac_ice_above, frac_ice_fall, cfkeep, cfliqkeep, cficekeep,   &
                rho, rhor, tcgi, tcgci,                                        &
                corr, dhi, dhir, rainfrac,                                     &
                pifall, psfall, prfall, pgfall,                                &
                one_over_tsi,                                                  &
                cf_transfer_diag, cff_transfer_diag,                           &
                uk, vk, ukp1, vkp1,                                            &
                r_theta_levels_c, fv_cos_theta_latitude_c,                     &
                l_use_agg_vt,                                                  &
                vm_used                                                        &
               )

  ! Water tracers ice and rain
  if (l_wtrac) then
    call lsp_sedim_eulexp_wtrac(points, rho, rhor, dhi, dhir, t,               &
                      vf_agg, qcf_agg, 'ice', wtrac_mp_cpr_old, wtrac_mp_cpr)

    if (l_mcr_qrain) then
      call lsp_sedim_eulexp_wtrac(points, rho, rhor, dhi, dhir, t,             &
                      vf_rain, qrain, 'rain', wtrac_mp_cpr_old, wtrac_mp_cpr)
    end if
  end if

else if ( sediment_loc == rain_sed_end   .or.                                  &
          sediment_loc == warm_sed_end ) then

  ! In these cases, the sedimentation of rain alone (or of warm species)
  ! is performed at the end of the microphysics.
  ! If sediment_loc == warm_sed_end, then lsp_settle is called immediately
  ! above. But in both cases, we need to do the fall of rain here, right
  ! at the end of the microphysics.

  call lsp_fall_rain(points, qrain, rho, rhor, corr, dhi, dhir,                &
                     rainrate, rainratet, vf_rain, one_over_tsi,               &
                     rainfrac, prfall )

  ! Water tracers rain
  if (l_wtrac) then
    if (l_mcr_qrain) then
      call lsp_sedim_eulexp_wtrac(points, rho, rhor, dhi, dhir, t,             &
                     vf_rain, qrain, 'rain', wtrac_mp_cpr_old, wtrac_mp_cpr)
    end if
  end if

end if ! sediment_loc

if ( l_proc_fluxes ) then
  ! If we have applied process rates to fluxes, then the fall-out fluxes
  ! have been included in the updated qcf_agg, qcf_cry, qrain and qgraup...

  ! Repartition ice, rain and graupel masses into an amount left on the current
  ! model-level and a flux falling out, assuming they keep the same ratio
  ! they had just after the sedimentation calculation.

  ! Ice crystals
  do i = 1, points
    snowt_cry(i) = (one - cry_nofall(i)) * qcf_cry(i)                          &
                                         / ( dhi(i) * rhor(i) )
    qcf_cry(i) = cry_nofall(i) * qcf_cry(i)
  end do
  ! Ice aggregates
  do i = 1, points
    snowt_agg(i) = (one - agg_nofall(i)) * qcf_agg(i)                          &
                                         / ( dhi(i) * rhor(i) )
    qcf_agg(i) = agg_nofall(i) * qcf_agg(i)
  end do
  ! Rain
  if ( l_mcr_qrain ) then
    do i = 1, points
      rainratet(i) = (one - rain_nofall(i)) * qrain(i)                         &
                                            / ( dhi(i) * rhor(i) )
      qrain(i) = rain_nofall(i) * qrain(i)
    end do
  end if
  ! Graupel
  if ( l_mcr_qgraup ) then
    do i = 1, points
      graupratet(i) = (one - graup_nofall(i)) * qgraup(i)                      &
                                              / ( dhi(i) * rhor(i) )
      qgraup(i) = graup_nofall(i) * qgraup(i)
    end do
  end if

  ! Set the falling-out ice-cloud and precipitation fractions equal to their
  ! updated values at level k, since the fall-out flux has been handled
  ! as part of the mass on level k until now...
  do i = 1, points
    frac_ice_fall(i) = cficekeep(i)
  end do
  if ( l_mcr_precfrac ) then
    do i = 1, points
      precfrac_fall(i) = precfrac_k(i)
    end do
  end if

end if  ! ( l_proc_fluxes )

! ======================================================================


!                  COPY ICE/SNOW VARIABLES and FLUXES
!                        to OUTPUT VARIABLES


! ======================================================================
! Copy contents of ice/snow variables and fluxes to output variables
! to fall into next layer down
! ----------------------------------------------------------------------

do i = 1, points

  if (l_mcr_qcf2) then ! two ice prognostics
    qcf(i)  = qcf_agg(i)
    qcf2(i) = qcf_cry(i)
    snow_cry(i) = snowt_cry(i)
    snow_agg(i) = snowt_agg(i)
  else ! only one ice prognostic, put all snow in to snow_agg
    qcf(i) = qcf_cry(i) + qcf_agg(i)
    snow_cry(i) = zero  ! Redundant variable
    snow_agg(i) = snowt_cry(i) + snowt_agg(i)
  end if ! on l_mcr_qcf2

  if (l_mcr_qrain)  rainrate(i) = rainratet(i)
  if (l_mcr_qgraup) grauprate(i) = graupratet(i)

end do ! Points

! Update water tracer fluxes
if (l_wtrac) then
  do i_wt = 1, n_wtrac
    do i = 1, points
      wtrac_mp_cpr(i_wt)%lssnow(i) = wtrac_mp_cpr(i_wt)%snowratet(i)
    end do
    if (l_mcr_qrain) then
      do i = 1, points
        wtrac_mp_cpr(i_wt)%lsrain(i) =  wtrac_mp_cpr(i_wt)%rainratet(i)
      end do
    end if
  end do
end if

! ======================================================================
!              NUMERICAL TIDYING UP OF SMALL VALUES
! ======================================================================
call lsp_tidy(points, one_over_tsi,                                            &
             q, qcl, qcf, qcf2, qrain, qgraup, t,                              &
             area_liq, area_mix, area_ice,                                     &
             cfice, cficei, cfkeep, cfliqkeep, cficekeep,                      &
             rainfrac, rain_liq, rain_mix, rain_ice, rain_clear,               &
             q_ice, qs, qsl, snow_agg, snow_cry, rainrate, grauprate,          &
             rho, rhor, p,                                                     &
             cttemp,dhi,dhir,frac_ice_fall,                                    &
             lcrcp, lfrcp, lsrcp,                                              &
             psdep, pidep, psmlt, pimlt, prevp,                                &
             cf_transfer_diag,cfl_transfer_diag,                               &
             cff_transfer_diag, rf_transfer_diag,                              &
             precfrac_k, precfrac_fall,                                        &
             wtrac_mp_cpr, wtrac_mp_cpr_old                                    &
             )

      !------------------------------------------------
      ! Now update fraction of ice in layer above
      ! for next layer down
      !------------------------------------------------
if ( l_fix_mcr_frac_ice ) then
  ! Bug fix for the fraction of ice falling down to the next level;
  ! needs to be set using values intercepted in lsp_fall_ice,
  ! where it is consistent with the fall-out flux of ice that is
  ! transfered down to the next level.
  do i = 1, points
    frac_ice_above(i) = frac_ice_fall(i)
  end do
else
  ! Without this fix, ice fraction falling down to next level gets set
  ! with values of cff near the end of this routine, after other
  ! processes have modified the ice fraction at level k so that it
  ! is no longer consistent with the fall-flux.
  do i = 1, points
    frac_ice_above(i)=cficekeep(i)
  end do ! Points
end if



!========================================================================
!                   RADAR REFLECTIVITY DIAGNOSTICS
!========================================================================

if (l_ref_diag) then
  if ( l_mcr_precfrac ) then
    ! If using prognostic precip fraction, pass that in for the rain and
    ! graupel fractions
    call mphys_reflec(points, rho,t, qgraup, qcf, qcf_cry, qrain, qcl,         &
                      n_drop_tpr, cficekeep, cfliqkeep, precfrac_k, precfrac_k,&
                      tcg, tcgc, dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_l, dbz_r )
  else
    ! Not using prognostic fraction for rain/graupel; pass in rainfrac_impr
    ! (the fraction for graupel will not be used if l_mcr_precfrac is off)
    call mphys_reflec(points, rho,t, qgraup, qcf, qcf_cry, qrain, qcl,         &
                      n_drop_tpr, cficekeep, cfliqkeep, rainfrac_impr, dummy,  &
                      tcg, tcgc, dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_l, dbz_r )
  end if
end if

! ======================================================================
!       if DIAGNOSTIC RAIN, convert MASS (kg/kg) to FLUX (kg/m2/s)
! ======================================================================
if (.not. l_mcr_qrain) then

  if (l_mr_physics) then

    do i = 1, points

      ! Use mixing ratio formulation
      rainrate(i) = qrain(i) * rhodz_dry(i) / timestep

    end do ! points

    if (l_wtrac) then
      do i_wt  = 1, n_wtrac
        do i = 1, points
          wtrac_mp_cpr(i_wt)%lsrain(i) = qrain(i) * rhodz_dry(i) / timestep
        end do
      end do
    end if

  else ! l_mr_physics

    do i = 1, points

      ! Use specific humidity formulation
      rainrate(i) = qrain(i) * rhodz_moist(i) / timestep

    end do ! points

    if (l_wtrac) then
      do i_wt  = 1, n_wtrac
        do i = 1, points
          wtrac_mp_cpr(i_wt)%lsrain(i) = qrain(i) * rhodz_moist(i) / timestep
        end do
      end do
    end if

  end if  ! l_mr_physics

end if ! on prognostic rain mixing ratio

! ----------------------------------------------------------------------
!   End of the LSP_ICE subroutine
! ----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_ice
end module lsp_ice_mod
