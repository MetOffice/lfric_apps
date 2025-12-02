! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Deep convection scheme

module deep_conv_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Deep convection scheme
!   works on points diagnosed as deep in subroutine CONV_DIAG.
!   Called by GLUE_CONV.
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


character(len=*), parameter, private :: ModuleName='DEEP_CONV_6A_MOD'

contains

subroutine deep_conv_6a(nbl,nlev,ntra,n_wtrac,n_cca_lev,n_dp,trlev,            &
                       bland, delthvu,                                         &
                       exner_rho,                                              &
                       exner_layer_centres,                                    &
                       exner_layer_boundaries,                                 &
                       l_q_interact,                                           &
                       l_tracer,ntml,ntpar,                                    &
                       pstar,p_layer_centres,                                  &
                       p_layer_boundaries,                                     &
                       z_theta, z_rho,                                         &
                       r_theta, r_rho,                                         &
                       rho_theta, rho,                                         &
                       rho_dry_theta, rho_dry,                                 &
                       r2rho_th, r2rho,                                        &
                       dr_across_th, dr_across_rh,                             &
                       conv_prog_precip,                                       &
                       q,th,                                                   &
                       timestep,                                               &
                       u,v,w,uw0,vw0,w_max,wstar,qsat_lcl,                     &
                       entrain_coef,zlcl_uv,tnuc_nlcl,                         &
                       freeze_lev,recip_pstar,qse,                             &
                       l_scm_convss_dg,                                        &
                       g_ccp, h_ccp, ccp_strength,                             &

                       ! InOut
                       cf_frozen,cf_liquid,qcf,                                &
                       qcl,tracer,wtrac_e,w2p,conv_prog_flx,                   &
                       scm_convss_dg,                                          &

                       ! Out
                       cape_out,cclwp,ccw,cca,                                 &
                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,             &
                       dqclbydt,dthbydt,                                       &
                       dubydt,dvbydt,dtrabydt,                                 &
                       detrain_up,detrain_dwn,                                 &
                       entrain_up,entrain_dwn,                                 &
                       iccb,icct,lcca,                                         &
                       lcbase,lctop,rain,snow,                                 &
                       rain_3d, snow_3d, up_flux, up_flux_half,                &
                       dwn_flux,uw_deep,vw_deep,kterm,tcw,cca_2d,              &
                       ind_cape_reduced,cape_ts_used,cfl_limited,              &
                       ind_deep,dt_dd,dq_dd, area_ud,                          &
                       error_point)

use planet_constants_mod, only:                                                &
    r, kappa, pref, repsilon, c_virtual, g

use cv_run_mod, only:                                                          &
    l_mom, l_eman_dd, cldbase_opt_dp, cape_min,                                &
    w_cape_limit, cape_timescale, deep_cmt_opt,                                &
    cca2d_dp_opt,                                                              &
    l_anvil, l_cv_conserve_check,                                              &
    l_cmt_heating, l_snow_rain, l_fcape,                                       &
    l_conv_prog_flx, tau_conv_prog_flx, cape_ts_min, cape_ts_max,              &
    cnv_cold_pools, l_ccp_parcel_dp, ccp_buoyancy

use cv_param_mod, only:                                                        &
    total_condensed_water, srf_precip, avg_mass_flux,                          &
    a_land, a_sea, b_land, b_sea, w_cca,                                       &
    thpixs_deep, qpixs_deep, c_mass, wcape_fac,                                &
    max_dp_thpert, min_dp_thpert, max_dp_qpert_fac,                            &
    a_cape, b_cape, a_cb, b_cb, c_cb,                                          &
    term_undil, ccp_off, ccp_h_coef

use cv_dependent_switch_mod, only:                                             &
    dp_on, dp_new_termc

use cv_stash_flg_mod, only:                                                    &
    flg_up_flx,  flg_up_flx_half, flg_dwn_flx,                                 &
    flg_entr_up, flg_detr_up, flg_detr_up, flg_detr_dwn,                       &
    flg_entr_dwn, flg_uw_dp, flg_vw_dp, flg_mf_deep,                           &
    flg_area_ud

use scm_convss_dg_mod, only: scm_convss_dg_type

use water_constants_mod, only: lc, lf, tm

use conversions_mod, only: rsec_per_day

use yomhook,     only: lhook, dr_hook
use parkind1,    only: jprb, jpim
use ereport_mod, only: ereport
use umPrintMgr,  only: umPrint, umMessage, newline,                            &
                       printstatus, prstatus_normal
use lift_par_6a_mod,       only: lift_par_6a
use lift_undil_par_6a_mod, only: lift_undil_par_6a
use convec2_6a_mod,        only: convec2_6a
use water_loading_mod,     only: water_loading
use cor_engy_6a_mod,       only: cor_engy_6a
use mix_ipert_6a_mod,      only: mix_ipert_6a
use cmt_heating_mod,       only: cmt_heating
use layer_cn_6a_mod,       only: layer_cn_6a, deep
use eman_dd_rev_mod,       only: eman_dd_rev

use calc_3d_cca_mod,          only: calc_3d_cca
use cmt_mass_mod,             only: cmt_mass
use dd_all_call_6a_mod,       only: dd_all_call_6a
use deep_cmt_incr_mod,        only: deep_cmt_incr
use deep_grad_stress_mod,     only: deep_grad_stress
use deep_ngrad_stress_mod,    only: deep_ngrad_stress
use deep_turb_cmt_mod,        only: deep_turb_cmt
use eman_dd_mod,              only: eman_dd
use evap_bcb_nodd_all_6a_mod, only: evap_bcb_nodd_all_6a
use flag_wet_mod,             only: flag_wet
use wtrac_conv_mod,           only: l_wtrac_conv,                              &
                                    conv_e_wtrac_type, conv_p_wtrac_type,      &
                                    wtrac_alloc_conv_p, wtrac_dealloc_conv_p
use water_tracers_mod,        only: wtrac_info
use wtrac_calc_ratio_mod,     only: wtrac_calc_ratio_fn
use correct_small_q_conv_mod, only: correct_small_q_conv

use qsat_mod, only: qsat

implicit none

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent in:


integer, intent(in) :: nbl      ! No. of boundary layer levels

integer, intent(in) :: nlev     ! No. of model layers

integer, intent(in) :: ntra     ! No. of tracer fields

integer, intent(in) :: n_wtrac  ! No. of water tracers

integer, intent(in) :: n_cca_lev! No. of convective cloud
                                ! amount levels (1 for 2D,
                                               ! nlevs for 3D)

integer, intent(in) :: n_dp     ! No. of deep convection points

integer, intent(in) :: trlev    ! No. of model levels on which
                                ! tracers are included

logical, intent(in) :: bland(n_dp) ! Land/sea mask

real(kind=real_umphys), intent(in)    :: delthvu(n_dp)
                                     ! a measure of CAPE used to cal wcld

real(kind=real_umphys), intent(in)    :: exner_rho(n_dp,nlev)
                                             ! Exner on rho levels

real(kind=real_umphys), intent(in)    :: exner_layer_centres(n_dp,0:nlev) !Exner

real(kind=real_umphys), intent(in)    :: exner_layer_boundaries(n_dp,0:nlev)
                                ! Exner at half level above
                                ! exner_layer_centres

logical, intent(in) :: l_q_interact ! Switch allows overwriting
                                    ! parcel variables when
                                    ! calculating condensate incr.

logical, intent(in) :: l_tracer ! Switch for inclusion of tracers

integer, intent(in) :: ntml(n_dp) ! Top level of surface mixed
                                  ! layer defined relative to
                                  ! theta,q grid

integer, intent(in) :: ntpar(n_dp) ! Top level of initial parcel
                                   ! ascent in BL scheme defined
                                   ! relative to theta,q grid

real(kind=real_umphys), intent(in)    :: pstar(n_dp) ! Surface pressure (Pa)

real(kind=real_umphys), intent(in)    :: p_layer_centres(n_dp,0:nlev)
                                                    ! Pressure (Pa)

real(kind=real_umphys), intent(in)    :: p_layer_boundaries(n_dp,0:nlev)
                                                       ! Pressure
                                                       ! at half level above
                                                       ! p_layer_centres (Pa)

real(kind=real_umphys), intent(in) :: z_theta(n_dp,nlev)
                                            ! height of theta levels (m)
real(kind=real_umphys), intent(in) :: z_rho(n_dp,nlev)
                                            ! height of rho levels (m)
real(kind=real_umphys), intent(in) :: r_theta(n_dp,0:nlev)
                                            ! radius of theta levels (m)
real(kind=real_umphys), intent(in) :: r_rho(n_dp,nlev)
                                            ! radius of rho levels (m)
real(kind=real_umphys), intent(in) :: rho_theta(n_dp,nlev)
                                            ! density for theta lev (kg/m3)
real(kind=real_umphys), intent(in) :: rho(n_dp,nlev)
                                            ! density for rho lev (kg/m3)
real(kind=real_umphys), intent(in) ::                                          &
  rho_dry_theta(n_dp,nlev)  & ! dry density for theta lev (kg/m3)
 ,rho_dry(n_dp,nlev)          ! dry density for rho lev (kg/m3)

real(kind=real_umphys), intent(in) :: r2rho_th(n_dp,nlev)
                                            ! radius**2 density for
                                            ! theta lev (kg/m)
real(kind=real_umphys), intent(in) :: r2rho(n_dp,nlev)
                                            ! radius**2 density for
                                            ! rho lev (kg/m)
real(kind=real_umphys), intent(in) :: dr_across_th(n_dp,nlev)
                                            ! thickness of theta levels (m)
real(kind=real_umphys), intent(in) :: dr_across_rh(n_dp,nlev)
                                            ! thickness of rho levels (m)
real(kind=real_umphys), intent(in) :: conv_prog_precip(n_dp,nlev)
                                                ! Surface precipitation based
                                                ! 3d convective prognostic in
                                                ! kg/m2/s
real(kind=real_umphys), intent(in)    :: q(n_dp,nlev)
                                    ! Model mixing ratio (kg/kg)

real(kind=real_umphys), intent(in)    :: th(n_dp,nlev)
                                     !Model potential temperature (K)

real(kind=real_umphys), intent(in)    :: timestep ! Model timestep (s)

real(kind=real_umphys), intent(in)    ::                                       &
  u(n_dp,nlev)         & ! Model u field (m/s)
 ,v(n_dp,nlev)         & ! Model v field (m/s)
 ,w(n_dp,nlev)           ! Model w field (m/s)

real(kind=real_umphys), intent(in)    :: uw0(n_dp)
                                 ! U-comp of surface stress (N/m2)

real(kind=real_umphys), intent(in)    :: vw0(n_dp)
                                 ! V-comp of surface stress (N/m2)

real(kind=real_umphys), intent(in)    :: wstar(n_dp)
                                   ! Convective velocity scale (m/s)

real(kind=real_umphys), intent(in)    :: w_max(n_dp) ! max w in column
                                   !for use in scale dependent cape timescale

real(kind=real_umphys), intent(in)    :: entrain_coef(n_dp)
                                          ! entrainment coefficient

real(kind=real_umphys), intent(in)    :: qsat_lcl(n_dp)
                                      ! qsat at cloud base (kg/kg)

real(kind=real_umphys), intent(in)    :: zlcl_uv(n_dp)
                                     !Lifting condensation level
                                     ! defined for the uv grid (m)
real(kind=real_umphys), intent(in)    :: tnuc_nlcl(n_dp)
                                     !nucleation temperature as function of
                                     !dust indexed using nlcl(deg cel)
integer, intent(in) :: freeze_lev(n_dp) ! Level index for freezing level

real(kind=real_umphys), intent(in) :: recip_pstar(n_dp)
                                      ! Reciprocal of pstar array

real(kind=real_umphys), intent(in) :: qse(n_dp,nlev)
                                   ! Saturation mixing ratio of
                                   ! cloud environment (kg/kg)

real(kind=real_umphys), intent(in) ::   g_ccp(n_dp)                            &
                                          ! Cold_pool reduced gravity
                    , h_ccp(n_dp) &       ! Cold-pool depth
                    , ccp_strength(n_dp)  ! Cold-pool strength

! Arguments with intent INOUT:


real(kind=real_umphys), intent(in out) :: cf_frozen(n_dp,nlev)
                                             ! Frozen water cloud volume ( )

real(kind=real_umphys), intent(in out) :: cf_liquid(n_dp,nlev)
                                             ! Liq water cloud volume ( )

real(kind=real_umphys), intent(in out) :: qcf(n_dp,nlev)
                                       ! Ice condensate mix ratio (kg/kg)

real(kind=real_umphys), intent(in out) :: qcl(n_dp,nlev)
                                       ! Liq condensate mix ratio (kg/kg)

real(kind=real_umphys), intent(in out) :: tracer(n_dp,trlev,ntra)
                                                !Model tracer fields (kg/kg)

real(kind=real_umphys), intent(in out) :: w2p(n_dp,nlev)
                                       ! (Parcel vertical velocity)^2 [(m/s)^2]

real(kind=real_umphys), intent(in out) :: conv_prog_flx(n_dp,nlev)
                                                 ! Mass flux convective
                                                 ! prognostic in Pa/s

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                               ! Structure containing water
                                               ! tracer fields

! Structure containing SCM convection sub-step diagnostics
! (needs intent inout as contains allocatable arrays that need to
! retain their allocated status on input as well as output)
type(scm_convss_dg_type), intent(in out) :: scm_convss_dg( n_dp )
! Flag for SCM convection sub-step diagnostics
logical, intent(in) :: l_scm_convss_dg


! Arguments with intent out:

real(kind=real_umphys), intent(out) :: cape_out(n_dp)
                                    ! Saved convective available
                                    ! potential energy for diagnostic
                                    ! output (J/kg)

real(kind=real_umphys), intent(out) :: cclwp(n_dp)
                                    ! Condensed water path (kg/m^2)

real(kind=real_umphys), intent(out) :: ccw(n_dp,nlev)
                                ! Convective cloud liquid water
                                ! on model levels (kg/kg)

real(kind=real_umphys), intent(out) :: cca(n_dp,n_cca_lev)
                                ! Convective cloud amount
                                ! on model levels (0-1)

real(kind=real_umphys), intent(out) ::                                         &
  dbcfbydt(n_dp,nlev)  & ! Increments to total cld volume due to convection(/s)
 ,dcffbydt(n_dp,nlev)  & ! Increments to ice cloud volume due to convection (/s)
 ,dcflbydt(n_dp,nlev)  & ! Increments to liq cloud volume due to convection (/s)
 ,dqbydt(n_dp,nlev)    & ! Increments to q due to convection (kg/kg/s)
 ,dqcfbydt(n_dp,nlev)  & ! Increments to ice
                         ! condensate due to convection(kg/kg/s)
 ,dqclbydt(n_dp,nlev)  & ! Increments to liq condensate due to convection
                         ! (kg/kg/s)
 ,dthbydt(n_dp,nlev)   & ! Increments to potential temp. due to convection (K/s)
 ,dubydt(n_dp,nlev)    & ! Increments to U due to CMT (m/s2)
 ,dvbydt(n_dp,nlev)      ! Increments to V due to CMT (m/s2)

real(kind=real_umphys), intent(out) ::                                         &
  dtrabydt(n_dp,nlev,ntra)   !Increment to tracer due to convection (kg/kg/s)

real(kind=real_umphys), intent(out) ::                                         &
  detrain_up(n_dp,nlev)  & ! Fractional detrainment rate into updraughts
                           ! (Pa/s)
 ,detrain_dwn(n_dp,nlev) & ! Fractional detrainment rate into downdraughts
                           ! (Pa/s)
 ,entrain_up(n_dp,nlev)  & ! Fractional entrainment rate into updraughts
                           ! (Pa/s)
 ,entrain_dwn(n_dp,nlev)   ! Fractional entrainment rate into downdraughts
                           ! (Pa/s)

integer, intent(out) :: iccb(n_dp) ! Convective cloud base level

integer, intent(out) :: icct(n_dp) ! Convective cloud top level

real(kind=real_umphys), intent(out) :: lcca(n_dp) ! Lowest conv. cloud amt. (%)

integer, intent(out) :: lcbase(n_dp) ! Lowest conv. cloud base level

integer, intent(out) :: lctop(n_dp) ! Lowest conv. cloud top level

real(kind=real_umphys), intent(out) :: rain(n_dp)
                                ! Surface convective rainfall (kg/m2/s)

real(kind=real_umphys), intent(out) :: snow(n_dp)
                                ! Surface convective snowfall (kg/m2/s)

real(kind=real_umphys), intent(out) :: rain_3d(n_dp,nlev)
                                        ! Convective rainfall flux (kg/m2/s)

real(kind=real_umphys), intent(out) :: snow_3d(n_dp,nlev)
                                        ! Convective snowfall flux (kg/m2/s)

real(kind=real_umphys), intent(out) :: up_flux(n_dp,nlev)
                                        ! Updraught mass flux (Pa/s)

real(kind=real_umphys), intent(out) :: up_flux_half(n_dp,nlev)
                                             ! Updraught mass flux on
                                             ! half levels (Pa/s)
real(kind=real_umphys), intent(out) :: dwn_flux(n_dp,nlev)
                                         ! Downdraught mass flux (Pa/s)

real(kind=real_umphys), intent(out) ::                                         &
  uw_deep(n_dp,nlev)  & ! X-comp. of stress from deep convection (kg/m/s2)
 ,vw_deep(n_dp,nlev)    ! Y-comp. of stress from deep convection (kg/m/s2)

integer, intent(out) :: kterm(n_dp) ! Level at which deep
                                    ! convection terminates,
                                    ! required by mid level scheme
real(kind=real_umphys), intent(out) :: tcw(n_dp)
                                 ! Total condensed water(kg/m2/s)
                                 ! required by mid-level CCA cal.

real(kind=real_umphys), intent(out) :: cca_2d(n_dp)
                                  ! 2D convective cloud amount (%)

real(kind=real_umphys), intent(out) ::                                         &
  ind_cape_reduced(n_dp)  & ! 1.0 - if CAPE reduced applies to several
                            ! CAPE options
 ,cape_ts_used(n_dp)      & ! cape timescale used for deep convection (s)
 ,cfl_limited(n_dp)       & ! Indicator of CFL limited convection
 ,ind_deep(n_dp)          & ! 1.0 if real deep convection else 0.0
 ,dt_dd(n_dp,nlev)        & ! dT/dt from DD and evap below cloud base (K/s)
 ,dq_dd(n_dp,nlev)        & ! dq/dt from DD and evap below cloud base (kg/kg/s)
 ,area_ud(n_dp,nlev)        ! fractional area of updraughts

integer, intent(out) :: error_point     ! 0 no error
                                        ! location of problem deep point

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

integer :: index1(n_dp),index2(n_dp)

integer :: ncposs               ! No. of points which may convect

integer :: nconv                ! No. of convecting points

real(kind=real_umphys) :: amdetk(n_dp)
                                ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

real(kind=real_umphys) :: b_calc                  ! Coefficient in thpert calc.

real(kind=real_umphys) :: c_calc                  ! Coefficient in thpert calc.

real(kind=real_umphys) :: cape(n_dp)
                                ! Convective available potential
                                ! energy (J/kg)
real(kind=real_umphys) :: fcape(n_dp)
                                ! Convective available potential
                                ! energy weighted by f.det profile (J/kg)

real(kind=real_umphys) :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

real(kind=real_umphys) :: dcpbydt(n_dp)
                                ! Rate of change of cape (J/kg/s)

real(kind=real_umphys) :: depth(n_dp)
                                ! Depth of convective cloud (m)

real(kind=real_umphys) :: ekp14(n_dp)             ! Entrainment coefficients at
                                ! level k+1/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

real(kind=real_umphys) :: ekp34(n_dp)             ! Entrainment coefficients at
                                ! level k+3/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

real(kind=real_umphys) :: exk(n_dp)               ! Exner ratio at layer k

real(kind=real_umphys) :: exkp1(n_dp)             ! Exner ratio at layer k+1

real(kind=real_umphys) :: flxmax(n_dp)            ! Maximum initial convective
                                ! mass flux (Pa/s)

real(kind=real_umphys) :: flx_init(n_dp)
                                ! Initial mass flux at cloud base
                                ! (Pa/s)

real(kind=real_umphys) :: flx_init_new(n_dp)
                                ! flx_init scaled to destroy cape
                                ! over timescale cape_timescale (Pa/s)

real(kind=real_umphys) :: flxmax_init(n_dp)
                                ! Maximum possible initial mass
                                ! flux (limited to the mass in
                                ! the initial convecting layer
                                ! in Pa/s)

real(kind=real_umphys) :: flxbyrho_int(n_dp)
                                ! mass weighted vertical integral of
                                ! mass flux/rho used in cca_2d calculation
                                ! (m/s)

real(kind=real_umphys) :: max_cfl(n_dp)
                                ! Max cfl ratio over a convecting
                                ! layer

real(kind=real_umphys) :: decay_amount
                                ! decay fraction for time-smoothed mass flux

real(kind=real_umphys) :: tmp_conv_prog_flx(n_dp,nlev)
                                     ! Copy of Mass flux convective
                                ! prognostic in Pa/s

real(kind=real_umphys) :: p_lcl(n_dp)             ! Pressure at LCL (Pa)

real(kind=real_umphys) :: precip(n_dp,nlev)
                                ! Amount of precip from each layer
                                ! from each layer (kg/m/s)

real(kind=real_umphys) :: pk(n_dp)
                                ! Pressure at midpoint of layer
                                ! k (Pa)

real(kind=real_umphys) :: pkp1(n_dp)
                                ! Pressure at midpoint of layer
                                ! k+1 (Pa)

real(kind=real_umphys) :: delpk(n_dp)
                                ! Pressure difference over layer
                                ! k (Pa)

real(kind=real_umphys) :: delpkp1(n_dp)
                                ! Pressure difference over layer
                                ! k+1 (Pa)

real(kind=real_umphys) :: delpkp12(n_dp)          ! Pressure difference between
                                ! layers k and k+1 (Pa)

real(kind=real_umphys) :: delp_uv_k(n_dp)
                                ! Pressure difference across uv
                                ! layer k (Pa)

real(kind=real_umphys) :: delp_uv_kp1(n_dp)
                                ! Pressure difference across uv
                                ! layer k+1 (Pa)

real(kind=real_umphys) :: q_lcl(n_dp)             ! Mixing ratio at LCL (kg/kg)

real(kind=real_umphys) :: qse_lcl(n_dp)           ! Saturated q at LCL (kg/kg)

real(kind=real_umphys) :: rhum(n_dp)              ! Dummy relative humidity
                                ! (only used on shallow points)

real(kind=real_umphys) :: t_lcl(n_dp)             ! Temperature at LCL (K)

real(kind=real_umphys) :: th_lcl(n_dp)            ! Theta at LCL (K)

real(kind=real_umphys) :: thv_pert(n_dp)
                                ! Theta_v parcel perturbation (K)

real(kind=real_umphys) :: thpert(n_dp)
                                ! Theta parcel perturbation (K)

real(kind=real_umphys) :: qpert(n_dp)
                                ! q parcel perturbation (kg/kg)

real(kind=real_umphys) :: wsc_o_mb(n_dp)          ! Dummy argument for layer_cn
                                ! Convective velocity scale divided
                                ! by cloud base mass flux mb

integer :: start_lev(n_dp)      ! Convection initiation level

logical :: bgmk(n_dp)           ! Mask for points where parcel in
                                ! layer k is saturated
logical :: bgmk_term(n_dp)      ! Mask for points where parcel in
                                ! layer k is saturated at termination

logical :: blatent(n_dp)        ! Mask for points where latent heat has
                                ! been released

logical :: bwater(n_dp,2:nlev)  ! Mask for points at which
                                ! condensate is liquid

logical :: blowst(n_dp)         ! Dummy variable indicating low
                                ! enough stability for convection
                                ! to occur

logical :: bterm(n_dp)          ! Mask for points which have
                                ! stopped convecting

logical :: bconv(n_dp)          ! Mask for points at which
                                ! convection is occurring

logical :: bcposs(n_dp)         ! Mask for points passing
                                ! initial stability test


! Parcel variables


real(kind=real_umphys) :: qpi   ! Initial parcel mixing ratio
                                !(kg/kg)

real(kind=real_umphys) :: qp(n_dp,nlev)           ! Parcel mixing ratio (kg/kg)

real(kind=real_umphys) :: thpi  ! Initial parcel potential temp. (K)

real(kind=real_umphys) :: thp(n_dp,nlev)          ! Parcel potential temp (K)

real(kind=real_umphys) :: trap(n_dp,nlev,ntra)
                                ! Tracer content of parcel (kg/kg)

real(kind=real_umphys) :: flx(n_dp,nlev)          ! Parcel massflux (Pa/s)

real(kind=real_umphys) :: xsbmin_v(n_dp,nlev)
                                ! Minmum parcel buoyancy excess

real(kind=real_umphys) :: thpixs_v(n_dp,nlev)     ! Theta parcel excess (K)

real(kind=real_umphys) :: qpixs_v(n_dp,nlev)      ! Q parcel excess(kg/kg)

! Water tracers
type(conv_p_wtrac_type) :: wtrac_p(n_wtrac)
                                ! Structure containing water
                                ! tracer fields relating to the parcel
real(kind=real_umphys), allocatable :: qpixs_v_wtrac(:,:,:)
                                ! Water tracer q parcel excess (kg/kg)
real(kind=real_umphys), allocatable :: ratio_wtrac(:,:,:)
                                ! Ratio of water tracer to normal water

! Undilute parcel variables
real(kind=real_umphys) :: qu(n_dp,nlev)
                                ! Undilute Parcel mixing ratio (kg/kg)
real(kind=real_umphys) :: thu(n_dp,nlev)
                                ! Undilute Parcel potential temp (K)


! PC2

real(kind=real_umphys) :: qclp(n_dp,nlev)
                                ! Parcel liquid condensated mixing
                                ! ratio in layer k (kg/kg)

real(kind=real_umphys) :: qcfp(n_dp,nlev)
                                ! Parcel frozen condensated mixing
                                ! ratio in layer k (kg/kg)


! Parameters
real(kind=real_umphys), parameter :: cfl_limit = 1.0 ! Max CFL ratio allowed
real(kind=real_umphys), parameter :: minflx = tiny(flx_init_new)
                                                ! minimum allowable
                                                ! initial mass flux

! CMT variables  - those used depend on scheme

! Required by Gregory-Kershaw scheme operating in plume calculation

integer ::                                                                     &
 nstart(n_dp)        ! Level for start of plume

real(kind=real_umphys) ::                                                      &
 eflux_u_ud(n_dp)  & ! Vertical eddy flux of momentum due to UD at
                     !  top of layer (Pa m/s)
,eflux_v_ud(n_dp)  & ! Vertical eddy flux of momentum due to UD at
                     ! bottom of layer (Pa m/s)
,up(n_dp,nlev)     & ! Parcel U (m/s)
,vp(n_dp,nlev)     & ! Parcel V (m/s)
,zsurf(n_dp)         ! Height of start of plume = 0.1*zlcl

! Required by Turbulence base scheme called after plume calculation

integer ::                                                                     &
 nlcl_uv(n_dp)         & ! Level index for LCL
,ntop_uv(n_dp)         & ! Level index for top of layer
,n_0degc(n_dp)         & ! Level index for zero degrees
,cu_term(n_dp)         & ! Indicies for CMT subs
,cu_tend(n_dp)           ! Indicies for CMT subs

real(kind=real_umphys) ::                                                      &
 mass_dwn(nlev,n_dp)   & ! Downdraught mass flux (Pa/s)
,p_uv(nlev,n_dp)       & ! Pressure of model level (Pa)
,phalf_uv(nlev,n_dp)   & ! Pressure of half level (Pa)
,plcl_uv(n_dp)         & ! Pressure at LCL (Pa)
,ptop_uv(n_dp)         & ! Pressure at top of cloud layer (Pa)
,p_0degc_uv(n_dp)      & ! Pressure of zero degree level (Pa)
,rho_uv(nlev,n_dp)     & ! Density on uv level (kg/m3)
,visc(nlev,n_dp)       & ! CMT eddy viscosity (m2/s)
,uw(nlev,n_dp)         & ! U- comp stress profile (N/m2)
                         ! (units vary through calls)
,vw(nlev,n_dp)         & ! V-comp stress profile (N/m2)
,uw_base(nlev,n_dp)    & ! Cloud base U stress (N/m2)
,vw_base(nlev,n_dp)    & ! Cloud base V stress (N/m2)
,ue_p(nlev,n_dp)       & ! Environment U profile (m/s)
,ve_p(nlev,n_dp)         ! Environment V profile (m/s)

real(kind=real_umphys) :: exk_temp                ! Temporary exner

! Required by all version of CMT

real(kind=real_umphys) :: flxkp12(nlev,n_dp)
                                ! Mass flux on half level (Pa/s)

real(kind=real_umphys) :: mb(n_dp)
                                ! Cloud base mass flux (m/s or Pa/s depending
                                ! on location in code)

logical :: l_mom_gk             ! true if Gregory-Kershaw CMT required
logical :: l_mom_gk_stable      ! true if stabilized Gregory-Kershaw CMT
                                ! (different from the original)

! Cape scaling/closure variables

integer :: det_lev(n_dp)        ! Level at which split final
                                ! detrainment last occurred

integer :: nterm                ! No. of points where conv.
                                ! has terminated

real(kind=real_umphys) :: tempnum
                                ! Temporary variable for storage


! Arrays used by various CAPE closures
real(kind=real_umphys) ::                                                      &
  scale_f(n_dp)         & ! scale factor
 ,cape_ts_new(n_dp)     & ! Used as variable in RH-based closure
 ,relh(n_dp)            & ! RH integral (average when convection terminates)
 ,rh_mean(n_dp)         & ! mean RH over cloud on termination
 ,wls_mean(n_dp)        & ! wls - sum over convecting levels
 ,mass_mean(n_dp)       & ! mass - sum over convecting levels
 ,dptot(n_dp)             ! Delta P integral

real(kind=real_umphys) :: deltaktot(n_dp)
                                ! Integrated forced detrainment

! Original downdraught scheme variables

integer :: npossdd              ! Max. no. of downdraughts
                                ! possible

integer :: nnodd                ! No. of downdraughts not possible

integer :: index_possdd(n_dp)   ! Index of downdraughts possible

integer :: index_nodd(n_dp)     ! Index of downdraughts not
                                ! possible
integer :: kmax_term            ! maximum termination level + 1

real(kind=real_umphys) :: deltap_cld
                                ! Pressure thickness of convective
                                ! cloud (Pa)

logical ::                                                                     &
  bgmkp1(n_dp)

real(kind=real_umphys) ::                                                      &
   qpkp1(n_dp)                                                                 &
  ,qclpkp1(n_dp)                                                               &
  ,qcfpkp1(n_dp)                                                               &
  ,thpkp1(n_dp)

real(kind=real_umphys) :: Qlkp1(n_dp)
real(kind=real_umphys) :: Qfkp1(n_dp)
real(kind=real_umphys) :: Frezkp1(n_dp)

real(kind=real_umphys) :: watldek(n_dp)
real(kind=real_umphys) :: watldpk(n_dp)
real(kind=real_umphys) :: watldekp1(n_dp)
real(kind=real_umphys) :: watldpkp1(n_dp)


! Local compressed arrays


real(kind=real_umphys) :: thrk(n_dp)
                             ! potential temperature of forced detrained air

real(kind=real_umphys) :: qrk(n_dp)
                             ! specific humidity of forced detrained air

real(kind=real_umphys) :: deltak(n_dp)         ! Parcel forced detrainment rate
                             ! in layer k multiplied by
                             ! appropriate layer thickness

real(kind=real_umphys) :: flxkp12_t(n_dp)      ! Half level mass flux (Pa/s)

real(kind=real_umphys) :: qpk(n_dp)
real(kind=real_umphys) :: thpk(n_dp)

real(kind=real_umphys) :: rbuoyk(n_dp)       ! Par. buoyancy at k (K)
real(kind=real_umphys) :: rbuoykp1(n_dp)    ! Par. buoyancy at k+1 (K)
real(kind=real_umphys) :: rbuoyukp1(n_dp)   ! undilute Par. buoy at k+1 (K)

logical :: b_nodd(n_dp)   ! points with no downdraught
logical :: b_dd(n_dp)     ! points with downdraught on termination

real(kind=real_umphys) ::                                                      &
  rh_test    &   ! critical RH value for convective closure option
                 ! cldbase_opt_dp == 6
                 ! (RH-based CAPE timescale, timestep limited, reduced by w)
 ,rh_fac         ! factor for calculation in the above closure

! required by CMT

real(kind=real_umphys) ::                                                      &
  wcld(n_dp)    &  ! Convective veloicty scale
 ,zlcl(n_dp)       ! lifting condensation level

real(kind=real_umphys), parameter :: qmin = 1.0e-8 ! Global minimum allowed Q

! local versions of CAPE timescale max/min
real(kind=real_umphys) :: cape_ts_min_use
real(kind=real_umphys) :: cape_ts_max_use

integer :: warning   !  local integer (-ve) to pass to ereport

! Loop counters

integer :: i,j,k,ktra,kt,kk,i_wt

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEEP_CONV_6A'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise error_point


error_point=0

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

! Initialise logicals

do i = 1,n_dp
  blowst(i)    = .true.
  bconv(i)     = .false.
  bcposs(i)    = .false.
  b_nodd(i)    = .false.
  b_dd(i)      = .false.
  blatent(i)   = .false.
  bgmk_term(i) = .false.
end do

do i = 1,n_dp
  ind_cape_reduced(i) = 0.0
  cape_ts_used(i)     = 0.0
  cfl_limited(i)      = 0.0
  kterm(i)            = 0
  ind_deep(i)         = 0.0
  scale_f(i)          = 0.0
  start_lev(i)        = ntml(i)
end do

!-----------------------------------------------------------------------
! 2.1  Initialise parcel properties and increment arrays
!-----------------------------------------------------------------------

!initialise parcel values over all levels
do k = 1, nlev
  do i = 1, n_dp
    qp(i,k)     = 0.0
    thp(i,k)    = 0.0
    qclp(i,k)   = 0.0
    qcfp(i,k)   = 0.0
    flx(i,k)    = 0.0
    precip(i,k) = 0.0
    ccw(i,k)    = 0.0
    area_ud(i,k) = 0.0
  end do
end do

if (dp_new_termc == term_undil) then
  do k = 1, nlev
    do i = 1, n_dp
      qu(i,k)     = 0.0
      thu(i,k)    = 0.0
    end do
  end do
end if

if (l_tracer) then
  do ktra = 1,ntra
    do k=1,nlev
      do i = 1,n_dp
        trap(i,k,ktra) = 0.0
      end do
    end do
  end do
end if

! Allocate water tracer arrays and initialise

call wtrac_alloc_conv_p(n_dp, nlev, n_wtrac, wtrac_p)

if (l_wtrac_conv) then

  do i_wt = 1, n_wtrac
    do k = 1, nlev
      do i = 1, n_dp
        wtrac_p(i_wt)%q(i,k)      = 0.0
        wtrac_p(i_wt)%qcl(i,k)    = 0.0
        wtrac_p(i_wt)%qcf(i,k)    = 0.0
        wtrac_p(i_wt)%precip(i,k) = 0.0
      end do
    end do
  end do
end if  !l_wtrac_conv

do k = 1,nlev
  do i = 1,n_dp
    dthbydt(i,k)  = 0.0
    dqbydt(i,k)   = 0.0
    dqclbydt(i,k) = 0.0
    dqcfbydt(i,k) = 0.0
    dbcfbydt(i,k) = 0.0
    dcflbydt(i,k) = 0.0
    dcffbydt(i,k) = 0.0
    dt_dd(i,k)    = 0.0
    dq_dd(i,k)    = 0.0
  end do
end do

if (l_mom) then
  do k = 1,nlev
    do i = 1,n_dp
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    end do
  end do
end if  ! L_mom

if (l_tracer) then
  do ktra = 1,ntra
    do k = 1,nlev
      do i = 1,n_dp
        dtrabydt(i,k,ktra) = 0.0
      end do
    end do
  end do
end if  ! L_tracer

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do k = 1, nlev
      do i = 1, n_dp
        wtrac_e(i_wt)%dqbydt(i,k)   = 0.0
        wtrac_e(i_wt)%dqclbydt(i,k) = 0.0
        wtrac_e(i_wt)%dqcfbydt(i,k) = 0.0
      end do
    end do
  end do
end if  ! l_wtrac_conv

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------
if (flg_up_flx .or. flg_mf_deep) then
  do k = 1,nlev
    do i = 1,n_dp
      up_flux(i,k)      = 0.0
    end do
  end do
end if
if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,n_dp
      up_flux_half(i,k) = 0.0
    end do
  end do
end if
if (flg_dwn_flx) then
  do k = 1,nlev
    do i = 1,n_dp
      dwn_flux(i,k)     = 0.0
    end do
  end do
end if
if (flg_entr_up) then
  do k = 1,nlev
    do i = 1,n_dp
      entrain_up(i,k)   = 0.0
    end do
  end do
end if
if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,n_dp
      detrain_up(i,k)   = 0.0
    end do
  end do
end if
if (flg_entr_dwn) then
  do k = 1,nlev
    do i = 1,n_dp
      entrain_dwn(i,k)  = 0.0
    end do
  end do
end if
if (flg_detr_dwn) then
  do k = 1,nlev
    do i = 1,n_dp
      detrain_dwn(i,k)  = 0.0
    end do
  end do
end if

if (l_mom) then
  if (flg_uw_dp) then
    do k = 1,nlev
      do i = 1,n_dp
        uw_deep(i,k)    = 0.0
      end do
    end do
  end if
  if (flg_vw_dp) then
    do k = 1,nlev
      do i = 1,n_dp
        vw_deep(i,k)    = 0.0
      end do
    end do
  end if
end if  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
do i = 1,n_dp
  cca_2d(i) = 0.0
  iccb(i)   = 0
  icct(i)   = 0
  tcw(i)    = 0.0
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
end do

do k= 1,n_cca_lev
  do i= 1,n_dp
    cca(i,k) = 0.0
  end do
end do


do i = 1,n_dp
  !-----------------------------------------------------------------------
  ! 2.5  Initialise diagnostics for scaling and closure calculations
  !-----------------------------------------------------------------------
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)         = 0.0
  fcape(i)        = 0.0
  cape_out(i)     = 0.0
  dcpbydt(i)      = 0.0
  max_cfl(i)      = 0.0
  det_lev(i)      = 0
  relh(i)         = 0.0
  dptot(i)        = 0.0
  deltaktot(i)    = 0.0
  wls_mean(i)     = 0.0
  mass_mean(i)    = 0.0

  !-----------------------------------------------------------------------
  ! 2.6  Initialise eddy flux arrays for updraught
  !-----------------------------------------------------------------------
  eflux_u_ud(i)   = 0.0
  eflux_v_ud(i)   = 0.0

  !-----------------------------------------------------------------------
  ! 2.7  Initialise surface precipitation arrays
  !-----------------------------------------------------------------------
  rain(i)         = 0.0
  snow(i)         = 0.0

  !-----------------------------------------------------------------------
  ! Initialise dummy variables that are not used by deep
  !-----------------------------------------------------------------------
  wsc_o_mb(i)     = 0.0
end do

! Initialise surface precip arrays for water tracers
if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i = 1,n_dp
      wtrac_e(i_wt)%rain(i) = 0.0
      wtrac_e(i_wt)%snow(i) = 0.0
    end do
  end do
end if



!Initialise adaptive entrainment variables
!initialise to level 2 'cos that's where parcel lift starts from
do i = 1,n_dp
  thpk(i) = th(i,2)
  qpk(i)  = q(i,2)
end do

! Initialise the cape timescale, limited by deep convection timestep.
if (cape_ts_min < timestep) then
  cape_ts_min_use = timestep
else
  cape_ts_min_use = cape_ts_min
end if
if (cape_ts_max < timestep) then
  cape_ts_max_use = timestep
else
  cape_ts_max_use = cape_ts_max
end if

!-----------------------------------------------------------------------
! Calculate parcel perturbations
!-----------------------------------------------------------------------

! Calculate XSBMIN and THPIXS constants based on layer thickness (Pa)
do k = 1,nlev-1
  do i = 1,n_dp
    xsbmin_v(i,k) = min( ((p_layer_centres(i,k) -                              &
              p_layer_centres(i,k+1))/5000.0),1.0) *0.2

    thpixs_v(i,k) = min( ((p_layer_centres(i,k) -                              &
              p_layer_centres(i,k+1))/5000.0),1.0) * thpixs_deep

    qpixs_v(i,k)  = qpixs_deep
  end do
end do  ! nlev

if (l_wtrac_conv) then

  allocate(qpixs_v_wtrac(n_dp,nlev,n_wtrac))
  allocate(ratio_wtrac(n_dp,nlev,n_wtrac))

  do i_wt = 1, n_wtrac
    do k = 1,nlev-1
      do i = 1,n_dp
        ! Calculate water tracer content of parcel q perturbation
        ratio_wtrac(i,k,i_wt) = wtrac_calc_ratio_fn(i_wt,                      &
                                wtrac_e(i_wt)%q(i,k), q(i,k))
        qpixs_v_wtrac(i,k,i_wt) = ratio_wtrac(i,k,i_wt) * qpixs_v(i,k)
      end do
    end do  ! nlev
  end do    ! n_wtrac
end if ! l_wtrac_conv


! Calculate cloud base mass flux
if ( cldbase_opt_dp == 8 .or. cldbase_opt_dp == 9 ) then
  ! Closure options that use the surface-flux based mass-flux at cloud-base

  ! Use value from CRM fits to deep simulations
  ! gives a value different to c_mass
  ! mb in units of m/s at this point.
  do i = 1,n_dp
    mb(i) = b_cb*a_cb * wstar(i)
  end do

else ! original value
  do i = 1,n_dp
    mb(i) = c_mass * wstar(i)
  end do
end if

! Define the LCL at the half level above ntml. Find environmental
! T at p_lcl by approximating theta there with
! th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is tunable.
! Similarly for q.

do i = 1,n_dp
  k=ntml(i)
  p_lcl(i)  = p_layer_boundaries(i,k)
  th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
  t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
  q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
end do

! Calculate saturation mixing ratio at LCL
call qsat(qse_lcl,t_lcl,p_lcl,n_dp)

!-----------------------------------------------------------------------
! Initialize arrays required for Convective Momentum Transport(CMT)
!-----------------------------------------------------------------------
if (l_mom) then
  select case (deep_cmt_opt)

  case (2,6)         ! Gregory-Kershaw CMT

    ! need level near surface for initial parcel U & V values
    ! zsurf = 0.1*z_lcl

    do i = 1,n_dp
      zsurf(i)  = 0.1*z_rho(i,ntml(i))
    end do
    do k=nlev-1,1,-1
      do i = 1,n_dp
        if (zsurf(i) <= z_theta(i,k)) then
          nstart(i) = k
        end if
      end do
    end do
    l_mom_gk = .true.

    ! Initialise winds for Gregory-Kershaw parcel calculation
    do k=1,nlev
      do i = 1,n_dp
        up(i,k) = 0.0
        vp(i,k) = 0.0
      end do
    end do

    if (deep_cmt_opt==6) then  ! Stabilized version
      l_mom_gk_stable = .true.
    else                      ! Original version
      l_mom_gk_stable = .false.
    end if

  case DEFAULT    ! (0/1/5) Alan Grant's eddy viscosity based CMT

    ! Note: In terms of array indices p and phalf follow the convention
    !       used in the boundary layer scheme. phalf(k,*) refers to the
    !       lower boundary of uv layer k. This follows the convention for
    !       um UM4.5 and before

    !       Also note that p_layer_boundaries(0) and p_layer_centres(0)
    !       = pstar, so p_uv(k,1) and phalf_uv(k,1) will be equal.

    !       Because of the definition of nlcl, the pressure of the top of
    !       the mixed layer is phalf_uv(nlcl,*)


    k=1
    do i = 1,n_dp
      p_uv(k,i)     = p_layer_boundaries(i,k-1)
      phalf_uv(k,i) = p_layer_centres(i,k-1)
      ue_p(k,i)     = u(i,k)
      ve_p(k,i)     = v(i,k)
      flxkp12(k,i)  = 0.0
      nlcl_uv(i)    = ntml(i) + 1
      n_0degc(i)    = freeze_lev(i)
    end do

    do i = 1,n_dp
      do k = 2,nlev
        p_uv(k,i)     = p_layer_boundaries(i,k-1)
        phalf_uv(k,i) = p_layer_centres(i,k-1)
        ue_p(k,i)     = u(i,k)
        ve_p(k,i)     = v(i,k)
        flxkp12(k,i)  = 0.0
        exk_temp      = (p_uv(k,i)/pref)**kappa
        rho_uv(k,i)   = 2.0 * p_uv(k,i) / (r * exk_temp *                      &
                        (th(i,k-1) + th(i,k)))
      end do
      plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
      p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
      rho_uv(1,i)     = rho_uv(2,i)
    end do

    l_mom_gk = .false.
    l_mom_gk_stable = .false.

  end select      ! test on deep_cmt_opt

else

  ! Initialise variable
  l_mom_gk = .false.
  l_mom_gk_stable = .false.

end if     !L_mom

! Calculate theta and q perturbation (perturbation is based on
! environment buoyancy gradient)
! First calculate thv_pert
do i = 1,n_dp

  thv_pert(i) = -0.5 * (th(i,ntml(i)+1)                                        &
              * (1.0+c_virtual * q(i,ntml(i)+1))                               &
              - th(i,ntml(i)) * (1.0 + c_virtual                               &
              * q(i,ntml(i)))) + 0.5

end do ! n_dp
! Test for inclusion of additional cold-pool perturbation
if ((cnv_cold_pools > ccp_off) .and. l_ccp_parcel_dp) then
  do i = 1,n_dp

    thv_pert(i) = thv_pert(i) +                                                &
            ccp_buoyancy * g_ccp(i) *                                          &
            exp(-ccp_h_coef*z_theta(i,ntml(i))/(1.0+h_ccp(i)))

  end do ! n_dp
end if
! Reset thpixs and qpixs at ntml
do i = 1,n_dp

  if (t_lcl(i) >  tm) then
    dq_sat_env  = repsilon * lc * qse_lcl(i)                                   &
                / (r * t_lcl(i) * t_lcl(i))
  else
    dq_sat_env  = repsilon * (lc+lf) * qse_lcl(i)                              &
                / (r * t_lcl(i) * t_lcl(i))
  end if

  b_calc      = t_lcl(i) * c_virtual * dq_sat_env + 1.0                        &
              + c_virtual * qse_lcl(i)

  c_calc    = th_lcl(i) * c_virtual * (qse_lcl(i)                              &
            - q_lcl(i)) - thv_pert(i)

  thpert(i) = max(min(-c_calc / b_calc, max_dp_thpert),                        &
              min_dp_thpert) ! ignore term in thpert**2

  thpixs_v(i,ntml(i)) = thpert(i)

  qpert(i)  = max(min(qse_lcl(i) + ((p_lcl(i) / pref)                          &
            **kappa) * thpert(i) * dq_sat_env                                  &
            - q_lcl(i),                                                        &
              max_dp_qpert_fac * qse_lcl(i)),0.0)

  qpixs_v(i,ntml(i))  = qpert(i)

end do ! n_dp

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i = 1,n_dp
      wtrac_p(i_wt)%qpert(i)        = ratio_wtrac(i,ntml(i),i_wt) * qpert(i)
      qpixs_v_wtrac(i,ntml(i),i_wt) = wtrac_p(i_wt)%qpert(i)
    end do
  end do
end if ! l_wtrac_conv

! Set bwater=.true. on points where water will condense rather than
! ice.
call flag_wet(n_dp,n_dp,nlev,th,exner_layer_centres,bwater)

!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------

do k = 2,nlev-1

  do i = 1,n_dp
    bterm(i)  = .false.
  end do

  !-----------------------------------------------------------------------
  ! Initialise environment variables
  ! NB These variable are only used by layer_cn.
  !-----------------------------------------------------------------------
  do i = 1,n_dp
    !Note that unlike p_layer_boundaries, where k indexing is offset
    !by one compared to the dynamics numbering, z retains the numbering
    !convention for dynamics variables i.e. for theta levels, k->k
    !and for rho levels k+1/2 -> k+1
    rhum(i)   = q(i,k) / qse(i,k)
  end do

  !-----------------------------------------------------------------------
  ! Initialise parcel properties (theta,q,tracer,momentum) if convection
  ! is not occurring at level k and has not convected in column before
  !-----------------------------------------------------------------------
  do i = 1,n_dp
    if ( .not. bconv(i) .and. det_lev(i) == 0) then
      bgmk(i)     = .false.
      depth(i)    = 0.0
      thpi        = th(i,k) + thpixs_v(i,k)
      thp(i,k)    = thpi
      thpk(i)     = thp(i,k)
      qpi         = q(i,k)  + qpixs_v(i,k)
      qp(i,k)     = qpi
      qpk(i)      = qp(i,k)
      if (dp_new_termc == term_undil) then
        thu(i,k)  = thpi
        qu(i,k)   = qpi
      end if
      if (l_q_interact) then
        qclp(i,k) = qcl(i,k)
        qcfp(i,k) = qcf(i,k)
      end if
      if (l_mom_gk) then  ! Gregory Kershaw CMT
        ! Set initial parcel values at cloud
        ! base to values of near surface winds
        up(i,k)   = u(i,nstart(i))
        vp(i,k)   = v(i,nstart(i))
      end if
    end if
  end do  ! n_dp
  if (l_tracer) then
    do ktra=1,ntra
      do i = 1,n_dp
        if ( .not. bconv(i)) then
          trap(i,k,ktra)  = tracer(i,k,ktra)
        end if  !not bconv
      end do
    end do
  end if

  ! Water tracers
  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i = 1,n_dp
        if ( .not. bconv(i) .and. det_lev(i) == 0) then
          wtrac_p(i_wt)%q(i,k) = wtrac_e(i_wt)%q(i,k) + qpixs_v_wtrac(i,k,i_wt)

          if (l_q_interact) then
            wtrac_p(i_wt)%qcl(i,k) = wtrac_e(i_wt)%qcl(i,k)
            wtrac_p(i_wt)%qcf(i,k) = wtrac_e(i_wt)%qcf(i,k)
          else
            wtrac_p(i_wt)%qcl(i,k) = 0.0
            wtrac_p(i_wt)%qcf(i,k) = 0.0
          end if
        end if
      end do
    end do
  end if   ! l_wtrac_conv

  !-----------------------------------------------------------------------
  ! 3.1  Calculate layer dependent constants (pressure,
  !      layer thickness, entrainment coefficients, detrainment
  !      coefficients)
  !-----------------------------------------------------------------------

  call layer_cn_6a(k, n_dp, nlev,                                              &
                   ntml, ntpar, start_lev,                                     &
                   exner_layer_centres,                                        &
                   p_layer_boundaries, p_layer_centres,                        &
                   z_rho,                                                      &
                   conv_prog_precip,                                           &
                   recip_pstar, entrain_coef, rhum,                            &
                   ccp_strength,                                               &
                   z_theta(:,k), z_rho(:,k+1), z_theta(:,k+1),                 &
                   wsc_o_mb, qsat_lcl, w_max,                                  &
                   deep,                                                       &
                   bconv,                                                      &
                   ! Out
                   pk, pkp1, exk, exkp1,                                       &
                   delpk, delpkp12, delpkp1,                                   &
                   delp_uv_k, delp_uv_kp1,                                     &
                   ekp14, ekp34, amdetk                                        &
                   )

  ! Maximum initial convective mass flux
  do i = 1,n_dp
    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)
  end do  ! n_dp

  !-----------------------------------------------------------------------
  ! Initial test to check if convection is possible in layer k
  !-----------------------------------------------------------------------
  ! Convection is possible if
  ! - the point was convecting (bconv = .T.) and did not terminate
  !   in the previous layer
  ! - or if at the top level of the surface mixed layer (k = ntml)
  do i = 1,n_dp
    bcposs(i) = bconv(i) .or. k  ==  ntml(i)
  end do  ! n_dp

  ! Calculate number of points which may convect (ncposs) and
  ! set compression indices (index1)
  ncposs = 0
  do i = 1,n_dp
    if (bcposs(i)) then
      ncposs          = ncposs + 1
      index1(ncposs)  = i
    end if
  end do

  !-----------------------------------------------------------------------
  ! 3.2  Lift parcel from layer k to layer k+1
  !-----------------------------------------------------------------------
  if (ncposs > 0) then

    call lift_par_6a(k, n_dp, n_wtrac, th(:,k), th(:,k+1),                     &
                q(:,k), q(:,k+1), qcl(:,k), qcf(:,k),                          &
                qcl(:,k+1), qcf(:,k+1),                                        &
                pkp1(:), exkp1(:),                                             &
                thp(:,k), qp(:,k), qclp(:,k),qcfp(:,k),                        &
                ekp14(:), ekp34(:),                                            &
                l_q_interact, bwater(:,k+1), wtrac_e,                          &
                ! Inout
                wtrac_p,                                                       &
                ! Out
                bgmkp1, thpkp1, qpkp1,                                         &
                qclpkp1, qcfpkp1,                                              &
                Qlkp1, Qfkp1, Frezkp1, index1, ncposs)

    !-----------------------------------------------------------------------
    ! Calculate the water loading for level k and k+1
    !-----------------------------------------------------------------------
    call water_loading(n_dp, qcl(:,k), qcf(:,k), qclp(:,k), qcfp(:,k),         &
                     watldek, watldpk, index1, ncposs)

    call water_loading(n_dp, qcl(:,k+1), qcf(:,k), qclpkp1, qcfpkp1,           &
                     watldekp1, watldpkp1, index1, ncposs)


    if (dp_new_termc == term_undil) then
      !---------------------------------------------------------------------
      ! Calculate undilute parcel properties
      !---------------------------------------------------------------------
      call lift_undil_par_6a(n_dp,                                             &
                  pkp1(:), exkp1(:),                                           &
                  thu(:,k), qu(:,k),                                           &
                  Frezkp1(:), bwater(:,k+1),                                   &
                  !Out
                  thu(:,k+1), qu(:,k+1), index1, ncposs)
    end if

    !-----------------------------------------------------------------------
    ! Test if convection is starting from layer k
    !-----------------------------------------------------------------------

!DIR$ IVDEP
    do i = 1,ncposs ! Loop over points which may convect

      ! Calculate buoyancy (virt. potential temp.) of parcel in layer k and k+1
      rbuoyk(index1(i))   = thp(index1(i),k)                                   &
                    * (1.0 + c_virtual *qp(index1(i),k)                        &
                    - watldpk(index1(i)))                                      &
                    - th(index1(i),k) * (1.0 + c_virtual *q(index1(i),k)       &
                    - watldek(index1(i)))

      rbuoykp1(index1(i)) = thpkp1(index1(i))                                  &
                    * (1.0 + c_virtual *qpkp1(index1(i))                       &
                    - watldpkp1(index1(i)))                                    &
                    - th(index1(i),k+1) * (1.0 + c_virtual * q(index1(i),k+1)  &
                    - watldekp1(index1(i)))

      if (dp_new_termc == term_undil) then
        !undilute parcel buoyancy
        rbuoyukp1(index1(i))= thu(index1(i),k+1)                               &
                      * (1.0 + c_virtual *qu(index1(i), k+1)                   &
                      - watldpkp1(index1(i)))                                  &
                      - th(index1(i),k+1) * (1.0 + c_virtual *q(index1(i),k+1) &
                      - watldekp1(index1(i)))
      end if

      ! Allow parcel to convect from ntml.
      if (k  ==  ntml(index1(i))) then
        bconv(index1(i))  = .true.  ! convection active
        blowst(index1(i)) = .true.  ! convection initialised in layer

        ! Set parcel mass flux
        flx(index1(i),k)  = mb(index1(i)) * g                                  &
                          * p_layer_centres(index1(i),k)                       &
                          / ( r * thp(index1(i),k)                             &
                          * (p_layer_centres(index1(i),k)/ pref)**kappa )

        ! If mass flux out of the initial layer is greater than the mass flux
        ! of the layer over the timestep then limit mass flux to mass of layer.
        if (flx(index1(i),k) > flxmax(index1(i))) then
          flx(index1(i),k) = flxmax(index1(i))
        end if

        ! Store diagnostics linked to initial convective mass flux for
        ! calculation of final closure.
        flx_init(index1(i))    = flx(index1(i),k)
        flxmax_init(index1(i)) = flxmax(index1(i))

      else     ! k=ntml test
        blowst(index1(i)) = .false. ! convection not initialise in layer
      end if   ! k=ntml test


      ! Reset threshold for forced detrainment to the initial
      ! (positive or negative) buoyancy (limit positive buoy.
      ! threshold to XSBMIN fn(delta P)), only for first 5 levels of lift

      if (k  >=  ntml(index1(i)) .and.                                         &
          k  <=  ntml(index1(i)) + 4) then

        xsbmin_v(index1(i),k) = min ( xsbmin_v(index1(i),k),                   &
                                      thv_pert(index1(i)))
      end if

    end do  !ncposs
  end if    !ncposs>0

  ! Calculate number of points which are convecting  (nconv)
  ! set compression indices (index2).
  nconv = 0
  do i = 1,ncposs
    if (bconv(index1(i))) then
      nconv         = nconv + 1
      index2(nconv) = index1(i)
    end if
  end do

  !-----------------------------------------------------------------------
  ! Second compression to form arrays of length nconv to be passed
  ! into CONVEC2.
  ! NB This process is used to update some single level variables that are
  ! defined on level k and kp1 by using the full field variables.
  ! NB The order in which the variables are compressed
  ! is the same as the argument list for CONVEC2.
  !-----------------------------------------------------------------------
  if (nconv  >   0) then

    !-----------------------------------------------------------------------
    ! 3.3  Calculate the rest of the parcel ascent  and the effect of
    !      convection on the large-scale atmosphere.

    !      Subroutine CONVEC2

    !      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
    !-----------------------------------------------------------------------

    ! (Note, if using water tracers, the water tracer arrays passed into
    !  convec2 must have size 'npnts', i.e. the 2nd argument in the call)

    call convec2_6a  (k, n_dp, n_dp, nlev, ntra, n_wtrac, trlev,               &
                      dp_on, dp_new_termc, start_lev, timestep,                &
                      pk, pkp1, delpk,                                         &
                      delpkp1, delp_uv_k, delp_uv_kp1,                         &
                      exk, exkp1,                                              &
                      th(:,k), th(:,k+1), q(:,k), q(:,k+1),                    &
                      qcl(:,k), qcl(:,k+1), qcf(:,k), qcf(:,k+1),              &
                      qse(:,k), qse(:,k+1),                                    &
                      cf_liquid(:,k), cf_liquid(:,k+1),                        &
                      cf_frozen(:,k),  cf_frozen(:,k+1),                       &
                      thp(:,k), qp(:,k), qclp(:,k), qcfp(:,k),                 &
                      rbuoyk, rbuoykp1, rbuoyukp1,                             &
                      watldek, watldekp1, watldpk, watldpkp1,                  &
                      Qlkp1, Qfkp1, Frezkp1,                                   &
                      ekp14, ekp34, amdetk, flx(:,k), flx_init,                &
                      u(:,k), u(:,k+1), v(:,k), v(:,k+1),                      &
                      up(:,k), vp(:,k),                                        &
                      tracer,                                                  &
                      z_theta(:,k), z_rho(:,k+1),                              &
                      z_theta(:,k+1),                                          &
                      l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,       &
                      bgmk, bgmkp1, bwater(:,k),                               &
                      bwater(:,k+1), blowst, bland,                            &

                      ! In/out
                      lcbase, lctop,                                           &
                      thpkp1, qpkp1, qclpkp1, qcfpkp1,                         &
                      dthbydt(:,k), dqbydt(:,k), dqclbydt(:,k), dqcfbydt(:,k), &
                      tcw, depth, cclwp, lcca,                                 &
                      cape, fcape, dcpbydt,                                    &
                      relh, dptot, deltaktot, max_cfl,                         &
                      eflux_u_ud, eflux_v_ud,                                  &
                      dubydt(:,k), dvbydt(:,k),                                &
                      dtrabydt, trap, wtrac_e, wtrac_p,                        &
                      w2p(:,k), bterm, blatent, xsbmin_v(:,k),                 &

                      ! Out
                      iccb, icct,                                              &
                      dcflbydt(:,k), dcffbydt(:,k), dbcfbydt(:,k),             &
                      dthbydt(:,k+1), dqbydt(:,k+1),                           &
                      dqclbydt(:,k+1), dqcfbydt(:,k+1),                        &
                      dcflbydt(:,k+1), dcffbydt(:,k+1), dbcfbydt(:,k+1),       &
                      precip(:,k+1), thrk, qrk, deltak,                        &
                      flxkp12_t, flx(:,k+1),                                   &
                      cca_2d, ccw(:,k+1),                                      &
                      up(:,k+1), vp(:,k+1),                                    &
                      dubydt(:,k+1), dvbydt(:,k+1),                            &
                      w2p(:,k+1),tnuc_nlcl,                                    &

                      ! Indirect indexing
                      index2,nconv)

  end if ! nconv > 0

  !-----------------------------------------------------------------------
  ! Decompression of compressed variables coming out of
  ! of CONVEC2.
  ! NB The order in which the variables are decompressed
  ! is the same as the the argument list for CONVEC2.
  !-----------------------------------------------------------------------
  do i = 1,n_dp
    depth(i)      = 0.0
    bgmk(i)       = .false.
  end do

  if (nconv  >   0) then
    !Decompression for intent(in)
!DIR$ IVDEP
    do i = 1,nconv
      j = index2(i)
      bgmk(j)         = bgmkp1(j)
    end do
    !Decompression for intent(INOUT)
!DIR$ IVDEP
    do i = 1,nconv
      j = index2(i)
      thp(j,k+1)      = thpkp1(j)
      qp(j,k+1)       = qpkp1(j)
      qclp(j,k+1)     = qclpkp1(j)
      qcfp(j,k+1)     = qcfpkp1(j)
    end do
    if (l_mom) then    ! needed for all versions
      do i = 1,nconv
        j = index2(i)
        flxkp12(k,j)  = flxkp12_t(j)
      end do
    end if  ! L_mom
  end if      ! nconv > 0

  !-----------------------------------------------------------------------
  ! 3.4  Cape and CFL scaling - adjust initial mass flux so that cape is
  !      removed by convection over timescale cape_timescale.
  !-----------------------------------------------------------------------

  ! Set up integer nterm which is the total number of points where
  ! convection has terminated.


  nterm = 0
  do i = 1,n_dp
    if (bterm(i)) then
      nterm       = nterm + 1
      bgmk_term(i)= bgmk(i)
      rh_mean(i)  = relh(i)/dptot(i)
      kterm(i)    = k
      ! If convection has terminated write cape to diagnostic output
      ! variable (cape_out).
      if (l_fcape) then
        cape_out(i) = fcape(i)/deltaktot(i)
      else
        cape_out(i) = cape(i)
      end if
      cape(i)     = 0.0
      fcape(i)    = 0.0
      bconv(i)    = .false.
      det_lev(i)  = k+1
      do kk=ntml(i),k
        wls_mean(i)  = wls_mean(i)  + w(i,kk)*r2rho_th(i,kk)*dr_across_th(i,kk)
        mass_mean(i) = mass_mean(i) + r2rho_th(i,kk)*dr_across_th(i,kk)
      end do
      wls_mean(i) =  wls_mean(i)/mass_mean(i)
      ! Need to initialise mass-flux at half-levels `flxkp12(det_lev(i),i)`
      ! as we may apply closure scalings at this model level later.
      flxkp12(k+1,i) = 0.0
    end if
  end do

  !-----------------------------------------------------------------------
  ! Write out entrainment, detrainment and half-level mass flux diagnostics.
  ! They will be scaled by the full level mass flux outside
  ! of the level loop
  !-----------------------------------------------------------------------
  ! Calculate fractional entrainment rate for level k.
  if (flg_entr_up) then
!DIR$ IVDEP
    do i = 1,nconv
      j = index2(i)
      entrain_up(j,k) = (1.0 - deltak(j))                                      &
               * (1.0 - amdetk(j)) * (ekp14(j) + ekp34(j)                      &
               * (1.0 + ekp14(j)))
    end do
  end if

  ! Calculate fractional detrainment rate for level k
  if (flg_detr_up) then
!DIR$ IVDEP
    do i = 1,nconv
      j = index2(i)
      detrain_up(j,k) = -(amdetk(j)                                            &
                      + deltak(j) * (1.0 - amdetk(j)))
    end do
  end if

  ! Calculate the half level mass flux for level k
  ! Only the scaling factor between full level and half levels is calculated
  ! here. This is scaled by the full level mass flux outside the level loop
  if (flg_up_flx_half) then
!DIR$ IVDEP
    do i =1,nconv
      j = index2(i)
      up_flux_half(j,k) = (1.0 - deltak(j))                                    &
                * (1.0 - amdetk(j)) * (1.0 + ekp14(j))
    end do
  end if

  ! Copy SCM profile diagnostics for which profiles are not stored elsewhere
  if (l_scm_convss_dg) then
!DIR$ IVDEP
    do i = 1, nconv
      j = index2(i)

      ! Save buoyancy profiles at level k
      scm_convss_dg(j) % par_thetav_excess(k) = rbuoyk(j)
      scm_convss_dg(j) % par_thetav(k)                                         &
           = thp(j,k) * ( 1.0 + c_virtual*qp(j,k) - watldpk(j) )
      scm_convss_dg(j) % env_thetav(k)                                         &
           = th(j,k) * ( 1.0 + c_virtual*q(j,k) - watldek(j) )

      ! If this is the termination level, also save these at k+1
      if ( bterm(j) ) then
        scm_convss_dg(j) % par_thetav_excess(k+1) = rbuoykp1(j)
        scm_convss_dg(j) % par_thetav(k+1)                                     &
           = thpkp1(j) * ( 1.0 + c_virtual*qpkp1(j) - watldpkp1(j) )
        scm_convss_dg(j) % env_thetav(k+1)                                     &
           = th(j,k+1) * ( 1.0 + c_virtual*q(j,k+1) - watldekp1(j) )
      end if

      ! Entrainment and mixing detrainment profiles
      scm_convss_dg(j) % ekp14(k)        = ekp14(j)
      scm_convss_dg(j) % ekp34(k)        = ekp34(j)
      scm_convss_dg(j) % amdetk(k)       = amdetk(j)

      ! Save adaptive detrainment diagnostics at level k / k+1
      scm_convss_dg(j) % deltak(k)       = deltak(j)
      scm_convss_dg(j) % rbuoy_star(k+1) = rbuoykp1(j)
      scm_convss_dg(j) % xsbmin(k+1)     = xsbmin_v(j,k)
      if (k==ntml(j)) then
        scm_convss_dg(j) % rbuoy_star(k) = rbuoyk(j)
        scm_convss_dg(j) % xsbmin(k)     = rbuoyk(j)
      end if
      scm_convss_dg(j) % thrk(k)         = thrk(j)
      scm_convss_dg(j) % qrk(k)          = qrk(j)
      ! buoyancy excess of detrained air, assuming same waterloading as parcel
      if ( thrk(j) > 0.0 ) then
        scm_convss_dg(j) % thvrk_excess(k)                                     &
            = thrk(j) * ( 1.0 + c_virtual*qrk(j) - watldpk(j) )                &
            - th(j,k) * ( 1.0 + c_virtual*q(j,k) - watldek(j) )
      end if

      ! Save the updraft mass-flux before closure scaling
      scm_convss_dg(j) % up_flx_guess(k) = flx(j,k)

    end do
  end if

  !-----------------------------------------------------------------------
  ! 3.6  End of main loop over levels
  !-----------------------------------------------------------------------
end do

if (l_wtrac_conv) then
  deallocate(ratio_wtrac)
  deallocate(qpixs_v_wtrac)
end if

! If used, copy convection profile diagnostics
if (l_scm_convss_dg) then

  ! Copy profile diagnostics
  do k = 1, nlev
    do i = 1, n_dp
      if ( k>=ntml(i) .and. k<=det_lev(i) ) then
        scm_convss_dg(i) % par_theta(k) = thp(i,k)
        scm_convss_dg(i) % par_q(k)     = qp(i,k)
        scm_convss_dg(i) % par_qcl(k)   = qclp(i,k)
        scm_convss_dg(i) % par_qcf(k)   = qcfp(i,k)
      end if
      scm_convss_dg(i) % env_theta(k) = th(i,k)
      scm_convss_dg(i) % env_q(k)     = q(i,k)
      scm_convss_dg(i) % env_qcl(k)   = qcl(i,k)
      scm_convss_dg(i) % env_qcf(k)   = qcf(i,k)
    end do
  end do

end if


!-----------------------------------------------------------------------
! 4.0 Choice of cloud base closure option
!-----------------------------------------------------------------------
select case ( cldbase_opt_dp )

  ! Default 4a convection scheme - RH-based CAPE closure
case ( 0 )

  do i = 1,n_dp
    if (dcpbydt(i)  >   0.0) then
      cape_ts_new(i) =                                                         &
                min(max(900.0*(1.0 - rh_mean(i))/0.1,60.0),cape_timescale)
      if (cape_ts_new(i) < cape_timescale) then
        ind_cape_reduced(i) = 1.0
      end if
    end if
  end do

  ! Modified 4a convection scheme - RH-based CAPE closure
  ! timescale limited to timestep
case ( 1 )

  do i = 1,n_dp
    if (dcpbydt(i)  >   0.0) then
      cape_ts_new(i) =                                                         &
                  min(max(cape_timescale*(1.0-rh_mean(i))/0.4,timestep)        &
                ,cape_timescale)
      if (cape_ts_new(i) < cape_timescale) then
        ind_cape_reduced(i) = 1.0
      end if
    end if  ! dcpbydt > 0
  end do

  ! Fixed cape timescale
case ( 2 )

  do i = 1,n_dp
    if (dcpbydt(i)  >   0.0) then
      cape_ts_new(i) =  cape_timescale
    end if
  end do

  ! w based cape closure; if w_cape_limit > 1000. reverts to cape_timescale
case ( 3 )
  ! Initialise array to cape timescale and then alter as required.
  cape_ts_new(:) =  cape_timescale

  if ( w_cape_limit < 1000.0 ) then
    !  This section includes test on w_max
    do i = 1,n_dp
      if ( dcpbydt(i) > 0.0 ) then
        ! new denominator introduced at vn6.6
        if ( w_max(i) > w_cape_limit ) then
          cape_ts_new(i) =   cape_timescale * w_cape_limit/                    &
                    (w_cape_limit+ (w_max(i)-w_cape_limit)*wcape_fac)
          ! set indicator that CAPE reduced
          ind_cape_reduced(i) = 1.0
        end if !  w_max(i) > w_cape_limit
      end if  ! dcpbydt > 0
    end do
  end if  ! w_cape_limit

  ! w based cape closure with grid-box area scaling
case ( 4 )

  if ( w_cape_limit < 1000.0 ) then
    !  This section includes test on w_max
    do i = 1,n_dp
      if ( dcpbydt(i) > 0.0 ) then
        if ( w_max(i) > w_cape_limit ) then
          cape_ts_new(i) = cape_timescale * w_cape_limit/w_max(i)
        else
          cape_ts_new(i) = cape_timescale * cape_out(i) / cape_min +           &
                             cape_timescale * exp(-cape_out(i) / cape_min)
        end if !  w_max(i) > w_cape_limit
      end if  ! dcpbydt > 0
    end do
  else
    do i = 1,n_dp
      if ( dcpbydt(i) > 0.0 ) then
        cape_ts_new(i) = cape_timescale * cape_out(i) / cape_min +             &
                         cape_timescale * exp( - cape_out(i) /cape_min)
      end if  ! dcpbydt > 0
    end do
  end if  ! w_cape_limit

  ! w based cape closure (experimental option; not used);
case ( 5 )

  if ( w_cape_limit < 1000.0 ) then
    !  This section includes test on w_max
    do i = 1,n_dp
      if ( dcpbydt(i) > 0.0 ) then
        if ( w_max(i) > w_cape_limit ) then
          cape_ts_new(i) =   cape_timescale * w_cape_limit/w_max(i)
        else
          if ( rh_mean(i) >= 0.75 ) then
            cape_ts_new(i) = cape_timescale *( 0.2373 / (rh_mean(i))**5)
            ind_cape_reduced(i) = 1.0
          else
            cape_ts_new(i) = cape_timescale
          end if ! rh_mean(i) >= 0.75
        end if   !  w_max(i) > w_cape_limit
      end if  ! dcpbydt > 0
    end do
  else
    do i = 1,n_dp
      if ( dcpbydt(i) > 0.0 ) then
        if ( rh_mean(i) >= 0.75 ) then
          cape_ts_new(i) = cape_timescale *( 0.2373 / (rh_mean(i))**5)
        else
          cape_ts_new(i) = cape_timescale
        end if ! rh_mean(i) >= 0.75
      end if   ! dcpbydt > 0
    end do
  end if  ! w_cape_limit

  ! RH and w based CAPE option
  ! Expects a sensible w_cape_limit or will do nothing
case ( 6 )

  rh_test = 0.60         ! critical RH value
  rh_fac  = 1.0/ (1.0 - rh_test)

  if ( w_cape_limit < 1000.0 ) then
    !  This section includes test on w_max
    do i = 1,n_dp
      if ( dcpbydt(i) > 0.0 ) then
        ! work out any reduction in cape_timescale due to RH
        ! linearly falls to 1/2 given cape_timescale for RH above rh_test
        if ( rh_mean(i) >= rh_test ) then
          cape_ts_new(i) = cape_timescale*0.5*                                 &
                             (1.0+(1.0-rh_mean(i))*rh_fac)
          ind_cape_reduced(i) =1.0
        else
          cape_ts_new(i) = cape_timescale
        end if
        ! Further reduction if w_max above critical value
        if ( w_max(i) > w_cape_limit ) then
          cape_ts_new(i) = cape_ts_new(i) * w_cape_limit/                      &
                          (w_cape_limit + (w_max(i)-w_cape_limit)*wcape_fac)
          ind_cape_reduced(i) =1.0
        end if
        ! Limit CAPE timescale to convective timestep
        cape_ts_new(i)  = max(cape_ts_new(i), timestep)
      end if  ! dcpbydt > 0
    end do
  end if  ! w_cape_limit

  ! CAPE timescale dependent on large-scale w,
  ! or
  ! CAPE timescale dependent on large-scale w but with a lower limit of a
  !   cloud base mass flux of a_cb*b_cb*wstar (applied later)
case ( 7, 8 )

  do i = 1,n_dp
    if ( dcpbydt(i) > 0.0 ) then
      if (wls_mean(i) > 0.0) then
        cape_ts_new(i) = a_cape *(wls_mean(i)**b_cape)
      else    ! set to maximum value
        cape_ts_new(i) = cape_ts_max_use
      end if

      ! Limit CAPE timescale within the min/max values:
      cape_ts_new(i)  = min(cape_ts_new(i), cape_ts_max_use)
      cape_ts_new(i)  = max(cape_ts_new(i), cape_ts_min_use)
    end if  ! dcpbydt > 0
  end do

  ! Boundary layer and large-scale vertical velocity closure
  ! NO CAPE element to this closure
  ! Works out required new cloud base mass flux directly
case ( 9 )

  do i = 1,n_dp
    if ( dcpbydt(i) > 0.0 ) then
      if (wls_mean(i) > 0.0) then
        ! Closure based on large-scale vertical velocity and wstar
        ! Currently flx_init(i) = a_cb*b_cb*wstar*rho_cb
        ! Instead of a fixed sigma_cb = b_cb we want to use
        ! sigma_cb = b_cb+c_cb*wls_mean
        flx_init_new(i) = flx_init(i)*(1.0+c_cb*wls_mean(i)/b_cb)
      else    ! Assume a fixed fractional area
        flx_init_new(i) = flx_init(i)
      end if
    end if  ! dcpbydt > 0
  end do

end select        ! cldbase_opt_dp


!---------------------------------------------------------------------------
! Calculate closure scaling to apply to convective tendencies etc
!---------------------------------------------------------------------------
! If using some sort of CAPE-closure (all closure options except the bl / w one
if (cldbase_opt_dp /= 9) then
  do i = 1,n_dp
    if (dcpbydt(i) > 0.0) then
      ! Calculate new mass-flux at cloud-base, so-as to remove
      ! CAPE over time cape_ts_new(i)
      flx_init_new(i)   = flx_init(i) * cape_out(i) /                          &
                          (cape_ts_new(i)*dcpbydt(i))
      ! Apply limit for the mass-flux out of the originating layer
      if (flx_init_new(i) > flxmax_init(i)) then
        flx_init_new(i) = flxmax_init(i)
      end if
      ! Scale max_cfl with cape scale
      max_cfl(i)        = max_cfl(i) * flx_init_new(i) / flx_init(i)
    else
      ! If convective tendencies don't erode CAPE, set mass-flux to zero
      flx_init_new(i) = 0.0
    end if  ! dcpbydt > 0
  end do
end if


! If using surface-limited w-based CAPE timescale:
if (cldbase_opt_dp == 8) then
  do i = 1,n_dp
    if ( flx_init_new(i) > 0.0 .and. flx_init(i) > flx_init_new(i) ) then

      ! Use original closure i.e. surface based, if the surface-based
      ! closure implies a larger mass-flux than the CAPE-based closure

      ! First Scale max_cfl with cape scale
      max_cfl(i)        = max_cfl(i) * flx_init(i) / flx_init_new(i)

      ! And then apply the surface based closure
      flx_init_new(i)   = flx_init(i)

      ! Apply limit for the mass-flux out of the originating layer
      if (flx_init_new(i) > flxmax_init(i)) then
        flx_init_new(i) = flxmax_init(i)
      end if

    end if
  end do
end if

! Time smoothed mass flux for closure
if (l_conv_prog_flx) then

  decay_amount = timestep/tau_conv_prog_flx

  do i = 1,n_dp
    ! Calculate scaling.
    scale_f(i)      = flx_init_new(i) / flx_init(i)
  end do

  ! Take a copy of the smoothed mass flux in case of failed convection
  ! And update the smoothed mass flux
  do kt = 1, nlev
    do i = 1,n_dp
      tmp_conv_prog_flx(i,kt) = conv_prog_flx(i,kt)
      conv_prog_flx(i,kt)     = conv_prog_flx(i,kt) +                          &
                                decay_amount * scale_f(i) * flx(i,kt)
    end do
  end do

  ! And update flx_init_new to its time-smoothed value and update
  ! max_cfl appropriately
  do i = 1,n_dp
    max_cfl(i)      = max_cfl(i) * conv_prog_flx(i,start_lev(i)) /             &
                      flx_init_new(i)
    flx_init_new(i) = conv_prog_flx(i,start_lev(i))

    if (flx_init_new(i) > flxmax_init(i)) then
      max_cfl(i)      = max_cfl(i) * flxmax_init(i)/flx_init_new(i)
      flx_init_new(i) = flxmax_init(i)
    end if
  end do

end if

! Work out scaled mass flux needed to keep cfl ratio below limit.
! This applies whether CAPE or surface based closure.
do i = 1,n_dp
  ! Scale max_cfl by the timestep
  max_cfl(i)        = max_cfl(i) * timestep

  ! If the max CFL in the profile exceeds the limit, reduce the mass-flux
  ! accordingly
  if (max_cfl(i) > cfl_limit) then
    flx_init_new(i) = flx_init_new(i) * cfl_limit / max_cfl(i)
    cfl_limited(i)  = 1.0       ! flag that deep convection is CFL limited
  end if

  ! Also re-apply separately calculated CFL-limit for the initiating layer.
  if (flx_init_new(i)  >   flxmax_init(i)) then
    flx_init_new(i) = flxmax_init(i)
  end if

  max_cfl(i) = 0.0
end do

! Compute the final closure scaling and the effective CAPE timescale diagnostic
do i = 1,n_dp
  ! Calculate the final scaling.
  scale_f(i)      = flx_init_new(i) / flx_init(i)

  ! set flx_init to the new value to provide the real initial mass
  ! flux in all conditions
  flx_init(i)     = flx_init_new(i)

  if (flx_init_new(i) > 0.0) then
    ! Work out effective CAPE timescale being used taking into account
    ! all restrictions being applied to cloud base mass flux (diagnostic).
    cape_ts_used(i) = cape_out(i)/(scale_f(i) * dcpbydt(i))
  end if
end do


!---------------------------------------------------------------------------
! Test to catch "failed convection" events
!---------------------------------------------------------------------------

! First, save SCM diagnostic of whether the convection failed, and if so, why:
if ( l_scm_convss_dg ) then
  do i = 1,n_dp

    ! Convection failed due to the ascent being too shallow or cloud-free
    if ( icct(i)-iccb(i) <= 3 .or. ( .not. blatent(i) ) ) then
      scm_convss_dg(i) % status_deep = 1

      ! Convection failed because the closure set the mass-flux to zero
    else if ( flx_init_new(i)  <=  minflx ) then
      scm_convss_dg(i) % status_deep = 2

      ! Real deep convection occurred!
    else
      scm_convss_dg(i) % status_deep = 3

    end if

  end do
end if

! Do the actual test for failed convection, and reset outputs to zero
! where convection has failed
do i = 1,n_dp
  if ( (flx_init_new(i)  <=  minflx)                                           &
                  .or. ( (icct(i)-iccb(i)) <= 3 )                              &
                  .or. ( .not. blatent(i) ) ) then
    ! False convection diagnosed if:
    ! - the new initial mass flux is less than zero
    ! - or the convecting layer it too thin
    ! - or the parcel has never released latent heat
    ! 3d variables are reset below by setting scale_f to zero.
    flx_init(i)     = 0.0
    flx_init_new(i) = 0.0
    mb(i)           = 0.0
    scale_f(i)      = 0.0
    cca_2d(i)       = 0.0
    iccb(i)         = 0
    icct(i)         = 0
    tcw(i)          = 0.0
    cclwp(i)        = 0.0
    lcca(i)         = 0.0
    lctop(i)        = 0
    lcbase(i)       = 0
    ! kterm should be reset to zero if convection fails but because it is used
    ! for reasons of efficiency in 3d variable scaling it is reset later
    ! instead here
    !   kterm(i)        = 0
    cape_ts_used(i) = 0.0     ! Not real deep so set to zero
    ind_cape_reduced(i) = 0   ! Ensure not set to 1 as not real deep
  else
    ! True convection
    ind_deep(i)     = 1.0  ! real deep event indicator
  end if
end do


! Maximum termination level
kmax_term = 2
do i = 1,n_dp
  if (kterm(i)+1 >  kmax_term) then
    kmax_term = kterm(i)+1
  end if
end do

if (kmax_term > nlev) then
  kmax_term = nlev
end if
!-----------------------------------------------------------------------
! Apply cape and cfl scaling
!-----------------------------------------------------------------------
do kt = 2, kmax_term
  do i = 1,n_dp
    if (kt  >=  ntml(i) .and. kt <= kterm(i)+1 ) then

      dthbydt(i,kt)   = dthbydt(i,kt)  * scale_f(i)
      dqbydt(i,kt)    = dqbydt(i,kt)   * scale_f(i)
      dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
      dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
      dcflbydt(i,kt)  = dcflbydt(i,kt) * scale_f(i)
      dcffbydt(i,kt)  = dcffbydt(i,kt) * scale_f(i)
      dbcfbydt(i,kt)  = dbcfbydt(i,kt) * scale_f(i)

      if (l_mom_gk) then
        dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
        dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
      end if
      if (l_mom) then     ! required for all versions
        flxkp12(kt,i) = flxkp12(kt,i) * scale_f(i)
      end if
      if (l_tracer) then
        do ktra = 1,ntra
          dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
        end do
      end if
      if (l_wtrac_conv) then
        do i_wt = 1, n_wtrac
          wtrac_e(i_wt)%dqbydt(i,kt)   =                                       &
                wtrac_e(i_wt)%dqbydt(i,kt)   * scale_f(i)
          wtrac_e(i_wt)%dqclbydt(i,kt) =                                       &
                wtrac_e(i_wt)%dqclbydt(i,kt) * scale_f(i)
          wtrac_e(i_wt)%dqcfbydt(i,kt) =                                       &
                wtrac_e(i_wt)%dqcfbydt(i,kt) * scale_f(i)
          wtrac_p(i_wt)%precip(i,kt)   =                                       &
                wtrac_p(i_wt)%precip(i,kt)   * scale_f(i)
        end do
      end if
      flx(i,kt)    = flx(i,kt) * scale_f(i)
      precip(i,kt) = precip(i,kt) * scale_f(i)

    end if !kt>ntml and flx_init_new>0
  end do  ! i loop
end do  ! kt loop

!-----------------------------------------------------------------------
! Scale cloud fraction
! Additional check on scale_f needed because of logarithmic dependence
!-----------------------------------------------------------------------
do i = 1,n_dp
  if (scale_f(i) > 0.0) then
    cca_2d(i)  = cca_2d(i) + 0.06 * log(scale_f(i))

    ! Check scaled cloud fraction not smaller than minimum value
    ! (2.0E-5) or greater than unity.
    cca_2d(i) = max(2.0e-5, cca_2d(i))
    cca_2d(i) = min(1.0e+0, cca_2d(i))
  end if
end do      ! i

!-----------------------------------------------------------------------
! Reset kterm value now that increments etc for all levels have been set
! to zero in the cases of failed convection
!-----------------------------------------------------------------------
do i = 1,n_dp
  if (scale_f(i) == 0.0) then
    kterm(i) = 0
  end if
end do      ! i

!-----------------------------------------------------------------------
! Correct time smoothed mass flux in case of failed convection
!-----------------------------------------------------------------------
if (l_conv_prog_flx) then
  ! Update the smoothed mass flux
  do kt = 1, nlev
    do i = 1,n_dp
      ! Test for failed convection
      if (kterm(i) == 0) then
        conv_prog_flx(i,kt) = tmp_conv_prog_flx(i,kt)
      end if
    end do
  end do
end if


!-----------------------------------------------------------------------
! Write out updraught massflux diagnostics and scale the
! entrainment and detrainment diagnostics by the mass flux.
!-----------------------------------------------------------------------
if (flg_up_flx .or. flg_mf_deep) then
  do k = 1,nlev
    do i = 1,n_dp
      up_flux(i,k) = flx(i,k)
    end do
  end do
end if

if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,n_dp
      up_flux_half(i,k) = up_flux_half(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_entr_up) then
  do k = 1,nlev
    do i = 1,n_dp
      entrain_up(i,k) = entrain_up(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,n_dp
      detrain_up(i,k) = detrain_up(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_area_ud) then  ! updraught area
  do k = 1,nlev
    do i = 1,n_dp
      if (flx(i,k) > 0.0) then
        ! expect to get a valid wup
        ! updraught core area from wup and flx

        if (w2p(i,k) > 0.0) then
          area_ud(i,k) =  flx(i,k)/(g*sqrt(w2p(i,k))* rho_theta(i,k))
        else
          ! The calculation of w2p is not producing a sensible value
          ! Assume wcb =1m/s and calculate a value for cloud base as value
          ! equired for new DD scheme.
          if (k == iccb(i)+1) then
            area_ud(i,k) =  flx(i,k)/(g* rho_theta(i,k))
          end if
        end if
      end if
    end do
  end do
end if

!-----------------------------------------------------------------------
! 5.0  Mixes the increments from the initial parcel perturbation throughout
!      the subcloud layer
!-----------------------------------------------------------------------

call mix_ipert_6a(n_dp, nlev, nbl, n_wtrac, ntml,                              &
             p_layer_boundaries, exner_layer_centres,                          &
             dthbydt, dqbydt, wtrac_e, flx_init,                               &
             thpert, qpert, wtrac_p)

!-----------------------------------------------------------------------
! 6.0 Down draughts  - now 2 options
!                     (a) Emanuel downdraught scheme
!                     (b) Original mass flux code
!                     (c) New mass flux scheme
! Note the level at which deep convection terminates has been stored
! in the above updraught loop as Kterm.
!-----------------------------------------------------------------------

if (l_eman_dd) then

  !-----------------------------------------------------------------------
  ! (a) Emanuel downdraught scheme
  !-----------------------------------------------------------------------

  ! Work out maximum termination level
  kmax_term = 2
  do i = 1,n_dp
    if (kterm(i) >  kmax_term) then
      kmax_term = kterm(i)
    end if
  end do

  if (l_snow_rain) then   ! revised Emanuel DD

    call eman_dd_rev (n_dp,kmax_term,nlev,trlev,ntra,                          &
                  kterm,l_tracer,                                              &
                  exner_layer_centres,                                         &
                  p_layer_centres, p_layer_boundaries,                         &
                  scale_f, th, q, qse, tracer, precip,                         &
                  dthbydt, dqbydt, dtrabydt,                                   &
                  rain, snow ,dwn_flux, dt_dd, dq_dd)

  else

    call eman_dd (n_dp,kmax_term,nlev,trlev,ntra                               &
,                      kterm,l_tracer                                          &
,                      exner_layer_centres                                     &
,                      p_layer_centres, p_layer_boundaries                     &
,                      timestep, th, q, qse, tracer, precip                    &
,                      dthbydt, dqbydt, dtrabydt                               &
,                      rain, snow ,dwn_flux, dt_dd, dq_dd                      &
                     )
  end if

else

  !-----------------------------------------------------------------------
  ! (b) Original mass flux downdraughts & evaporation scheme
  !-----------------------------------------------------------------------
  ! 6.1  Downdraught calculation - on all points where convection is
  !      terminating.
  !      UM Documentation Paper 27
  !-----------------------------------------------------------------------

  do i = 1,n_dp
    if (kterm(i) >= ntml(i) .and. flx_init_new(i) >0.0) then
      k = kterm(i)
      tempnum = 0.0
      if (iccb(i)  >   0) then
        deltap_cld = p_layer_centres(i,iccb(i))                                &
                         - p_layer_centres(i,k)
        do kt = iccb(i), k+1
          tempnum = tempnum + precip(i,kt)
        end do
      else
        deltap_cld = 0.0
      end if


      ! Downdraughts possible if pressure thickness of convective
      ! cloud (deltap_cld) is greater than 15000Pa, the point is saturated
      ! and the precip. in the layer is greater than a threshold
      ! value (1E-12).
      ! Set logical for use later
      if (deltap_cld  >   15000.0 .and. bgmk_term(i) .and.                     &
                             tempnum  >   1e-12) then
        b_dd(i) = .true.
      else
        b_nodd(i) = .true.
      end if

    end if   ! test on whether deep really happened
  end do
  !-----------------------------------------------------------------------
  ! 6.2  Downdraught calculation - on all points where convection is
  !      terminating.
  !      UM Documentation Paper 27
  !-----------------------------------------------------------------------

  npossdd = 0
  do i = 1,n_dp
    if (b_dd(i)) then
      npossdd = npossdd +1
      index_possdd(npossdd) = i
    end if
  end do

  if (npossdd  >   0) then

    ! Work out maximum termination level
    kmax_term = 2
    do i = 1,npossdd
      if (kterm(index_possdd(i)) >  kmax_term) then
        kmax_term = kterm(index_possdd(i))
      end if
    end do

    call dd_all_call_6a (n_dp,npossdd,kmax_term,nlev,trlev,ntra, n_wtrac       &
,                      kterm, iccb, icct, index_possdd, l_tracer               &
,                      bwater(1,2)                                             &
,                      exner_layer_centres,exner_layer_boundaries              &
,                      p_layer_centres, p_layer_boundaries,pstar               &
,                      recip_pstar,timestep , cca_2d                           &
,                      thp, qp, th, q, qse, trap,tracer, flx,precip            &
,                      dthbydt, dqbydt, dtrabydt                               &
,                      rain, snow , rain_3d, snow_3d                           &
,                      wtrac_p, wtrac_e                                        &
,                      dwn_flux, entrain_dwn, detrain_dwn, dt_dd, dq_dd)


  end if

  !-----------------------------------------------------------------------
  ! 6.3 Surface precipitation calculation for terminating points with
  !     no downdraught (moved outside level loop) ie do this calculation
  !     on all points at the end.
  !-----------------------------------------------------------------------
  ! Points where no downdraught possible
  nnodd = 0
  do i = 1,n_dp

    if (b_nodd(i)) then
      nnodd = nnodd +1
      index_nodd(nnodd) = i
    end if
  end do

  if (nnodd  >   0) then

    ! Work out maximum termination level
    kmax_term = 2
    do i = 1,nnodd
      if (kterm(index_nodd(i)) >  kmax_term) then
        kmax_term = kterm(index_nodd(i))
      end if
    end do
    ! Only add 1 if kmax_term is less than model levels (which should be
    ! true).
    if (kmax_term  <  nlev ) then
      kmax_term = kmax_term + 1
    end if

    ! Surface precipitation calculation


    call evap_bcb_nodd_all_6a(n_dp,nnodd,n_wtrac,kmax_term,kterm               &
,                      iccb, index_nodd, bwater(1,2)                           &
,                      exner_layer_centres                                     &
,                      p_layer_centres, p_layer_boundaries                     &
,                      timestep , cca_2d, th, q, qse, precip                   &
,                      dthbydt, dqbydt                                         &
,                      rain, snow, rain_3d, snow_3d                            &
,                      dt_dd, dq_dd, wtrac_p, wtrac_e)

  end if

end if        ! test on down draught type

! Deallocate water tracer parcel fields
call wtrac_dealloc_conv_p(n_wtrac, wtrac_p)

!-----------------------------------------------------------------------
! 7.0  Convective Momentum Transport (if L_mom = .true.)
!-----------------------------------------------------------------------

if (l_mom) then

  select case (deep_cmt_opt)
  case (2,6)     ! Gregory-Kershaw deep CMT

    ! Do nothing here as calculated in parcel ascent earlier

  case DEFAULT   ! (0/1/5) Alan Grant's Eddy viscosity CMT

    ! altered to use kterm instead of ntpar
    do i = 1,n_dp
      if (kterm(i) >=  nlcl_uv(i)) then
        ntop_uv(i)    = kterm(i) + 1
      else     ! case where deep convection fails
        ! I think in this case the cloud base mass flux will be zero so
        ! there will be no CMT. (The value will not matter)
        ntop_uv(i)    = ntpar(i) + 1
      end if

      ptop_uv(i)    = phalf_uv(ntop_uv(i),i)
    end do

    nterm = 0

    ! Set cloud base mass flux equal to mass flux at half level below
    ! the LCL. mb now completely reset and units changed to Pa/s

    do i = 1,n_dp
      if (mb(i) > 0.0 .and. kterm(i) < nlev -1) then
        nterm = nterm + 1
        cu_term(nterm) = i
        cu_tend(nterm) = i
        mb(i) = flxkp12(nlcl_uv(i),i)
        do j = 1,nlev
          flxkp12(j,i) = 0.0
        end do
      else if (kterm(i) == nlev - 1) then
        ! Problem deep convection has gone to the top of the model
        ! Return location of deep problem point plus profiles
        if (printstatus >= prstatus_normal) then
          warning = -1
          write(umMessage,'(A)')                                               &
            "  PROBLEM: deep convection has gone to the top of the model."
          call umPrint(umMessage,src='deep_conv_6a')
          write(umMessage,'(a23,i6,a7,i6,a7,i6)')                              &
              ' Deep convection point ',i,                                     &
              ' kterm ',kterm(i),' nlcl  ',nlcl_uv(i)
          call umPrint(umMessage,src='deep_conv_6a')
          write(umMessage,'(A6,4A26)') 'k','Theta','q','qcl','qcf'
          call umPrint(umMessage,src='deep_conv_6a')
          do k=1,nlev
            write(umMessage,'(I6,4G26.18)') k,th(i,k),q(i,k),qcl(i,k),qcf(i,k)
            call umPrint(umMessage,src='deep_conv_6a')
          end do
          call ereport( 'deep_conv_6a', warning,                               &
            "Convection has gone to the top of the model."    //newline//      &
            "This usually means there is junk in the inputs"  //newline//      &
            "to the convection scheme.  Check the profiles at"//newline//      &
            "the grid-point where this has occured (these are"//newline//      &
            "printed in the run output file)." )
        end if
        error_point = i
        if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,        &
                                zhook_handle)
        return

      end if

      ! initialise output arrays as lower level subroutines don't set all
      ! values

      do j = 1,nlev
        uw(j,i)       =0.0
        vw(j,i)       =0.0
        uw_base(j,i)  =0.0
        vw_base(j,i)  =0.0
        visc(j,i)     =0.0
        mass_dwn(j,i) =0.0
      end do

    end do

    if (nterm  >   0) then
      call cmt_mass(n_dp, n_dp, nlev, nterm, cu_term,                          &
              kterm, cu_tend, nlcl_uv, ntop_uv,                                &
              mb, p_0degc_uv, plcl_uv, ptop_uv, phalf_uv,                      &
            ! Output arguments
              flxkp12 ,mass_dwn, visc)

      call deep_grad_stress(n_dp,nlev,nlcl_uv,ntop_uv,                         &
                            nterm,cu_term,                                     &
                            ue_p,ve_p,visc,phalf_uv,p_uv,timestep,             &
                            ! Output
                            uw,vw)


      call deep_ngrad_stress(n_dp,n_dp,n_dp,nterm,nlev,                        &
                             nlcl_uv,ntop_uv,cu_term,cu_tend,cu_tend,          &
                             pstar,uw0,vw0,zlcl_uv,ue_p,ve_p,                  &
                             flxkp12,p_uv,phalf_uv,rho_uv,timestep,            &
                             ! Input/output
                             uw,vw,                                            &
                             ! Output
                             uw_base,vw_base,uw_deep,vw_deep)


      call deep_cmt_incr(n_dp,n_dp,n_dp,nlev,nterm,                            &
                         nlcl_uv,ntop_uv,cu_term,cu_tend,                      &
                         zlcl_uv,phalf_uv,                                     &
                         uw,vw,                                                &
                         ! Output
                         dubydt,dvbydt)

    end if  ! nterm > 0

  case (3,4)        ! New Turbulence scheme using heights

    nterm = 0   ! count of number of deep points which actually convected
    do i = 1,n_dp
      if (kterm(i)  >=  nlcl_uv(i)) then
        nterm = nterm + 1
        cu_term(nterm) = i

        ! Use CAPE scaled mass flux as initial mass flux rather than CRM derived
        ! value to be consistent with thermodynamic part of convection.
        ! mb now completely reset, units now Pa/s
        mb(i) = flxkp12(nlcl_uv(i),i)
        zlcl(i) = z_rho(i,ntml(i))

        ! Cloud velocity scale - derived from CRM simulations
        !             wcld = (C_mass*wstar*CAPE)**(1/3)

        wcld(i) = (delthvu(i) * c_mass * wstar(i) * g / (th(i,ntml(i))         &
              * (1.0 + c_virtual * q(i,ntml(i)))))**0.3333
      end if
    end do

    if (nterm > 0) then
      call deep_turb_cmt (n_dp, nterm, nlev, deep_cmt_opt,                     &
                    ntml, kterm,cu_term,freeze_lev,                            &
                    timestep,                                                  &
                    uw0, vw0, mb, wcld, wstar ,zlcl,                           &
                    flx,                                                       &
                    z_rho, z_theta,rho,rho_theta,                              &
                    r2rho, r2rho_th, dr_across_th, dr_across_rh,               &
                    u, v,                                                      &
                    dubydt, dvbydt, uw_deep, vw_deep)

    end if           ! nterm > 0

  end select           ! deep_cmt_opt

end if  ! L_mom

!-----------------------------------------------------------------------
! 7.1  Add the dissipative heating from the CMT to the theta increment
!-----------------------------------------------------------------------
if (l_mom .and. l_cmt_heating) then
  call cmt_heating(n_dp, nlev,                                                 &
                   z_theta, z_rho, exner_layer_centres,                        &
                   u, v, dubydt, dvbydt,                                       &
                   ! Out
                   dthbydt)
end if


!-----------------------------------------------------------------------
! 8.0  Energy (and optionally water) correction calculation
!-----------------------------------------------------------------------
do i = 1,n_dp
  index1(i) = i
end do

if (l_cv_conserve_check) then
  call cor_engy_6a(n_dp, n_dp, nlev, n_wtrac, index1, timestep,                &
                   p_layer_boundaries, exner_layer_centres,                    &
                   exner_rho, r_theta, r_rho, rho_dry,                         &
                   r2rho, rho_theta, rho_dry_theta,                            &
                   dr_across_rh,                                               &
                   dubydt, dvbydt, dqclbydt, dqcfbydt,                         &
                   rain, snow, q, qcl, qcf, u, v,                              &
                   !In/Out
                   dthbydt, dqbydt, wtrac_e)
end if

!-----------------------------------------------------------------------
! 9.0  Correct negative/very small humidities
!-----------------------------------------------------------------------

call correct_small_q_conv(n_dp, n_dp, nlev, n_wtrac, index1, timestep, qmin,   &
                          r2rho_th, dr_across_th,  q, 'water', 'Deep',         &
                          dqbydt, wtrac_e=wtrac_e)

! Then use same method to correct any remaining negative/v small values of
! water tracer vapour

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    call correct_small_q_conv(n_dp, n_dp, nlev, n_wtrac, index1, timestep,     &
                              wtrac_info(i_wt)%qlimit, r2rho_th, dr_across_th, &
                              wtrac_e(i_wt)%q, 'wtrac', 'Deep',                &
                              wtrac_e(i_wt)%dqbydt)
  end do
end if

!-----------------------------------------------------------------------
! 10.0  3D - Convective cloud amount assumed 3d required ie L_3d_cca
!      is true in old code
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------------
!      CCRad - Calculate CCA for Deep levels only
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------
! 10.1 Calculate CCA_2D of Deep Cloud
!-----------------------------------------------------------------
select case (cca2d_dp_opt)
case (srf_precip)
  do i=1, n_dp
    if (iccb(i) > 0) then ! Deep convection was successful

      ! Use determination of CCA_2D based on surface precip. rate
      ! from deep cloud.

      ! NOTE: at present the a_ and b_ parameters for LAND and SEA
      !       are equal, so there will be no difference between
      !       land and sea points.

      if ((rain(i) + snow(i)) > 0.0) then
        if (bland(i)) then
          ! Land point
          tempnum = a_land + b_land                                            &
                  * log(rsec_per_day * (rain(i)+snow(i)))
        else
          ! Sea point

          tempnum = a_sea + b_sea                                              &
                  * log(rsec_per_day * (rain(i)+snow(i)))
        end if

        cca_2d(i) = max(2.0e-5, tempnum)

        ! Grab lowest cca value before any tuning occurs
        ! This will overwrite lcca in ni_conv_ctl only if neither
        ! shallow or deep have occured.
        lcca(i) = cca_2d(i)

      end if
    end if    ! iccb
  end do      ! i (n_dp)

case (avg_mass_flux)
  !Initialise the integrated flx/rho
  do i=1, n_dp
    flxbyrho_int(i) = 0.0
    mass_mean(i)    = 0.0
  end do

  !calculate the integrated flx/rho and recalc convecting layer mass
  do k = 1,nlev
    do i= 1, n_dp
      if (flx(i,k) > 0.0) then
        flxbyrho_int(i) = flxbyrho_int(i) + flx(i,k)/g  * dr_across_th(i,k)
        mass_mean(i)    = mass_mean(i) + rho_theta(i,k) * dr_across_th(i,k)
      end if
    end do
  end do

  do i=1, n_dp
    if (iccb(i) > 0) then ! Deep convection was successful
      !-----------------------------------------------------------
      ! CCA_2D is derived from the mass averaged mass flux
      !-----------------------------------------------------------

      cca_2d(i) = flxbyrho_int(i) / (w_cca * mass_mean(i))

      !Limit to avoid unphysical values
      cca_2d(i) = max(min(cca_2d(i), 1.0),0.0)

      ! Grab lowest cca value before any tuning occurs
      ! This will overwrite lcca in ni_conv_ctl only if neither
      ! shallow or deep have occured.
      lcca(i) = cca_2d(i)

    end if    ! iccb
  end do      ! i (n_dp)

case (total_condensed_water)
  ! cca_2d_dp left unchanged from code, which is based on
  ! TCW (Total Condensed Water) (This is a rate)

end select

! l_dcpl_cld4pc2 is set to true for 5a scheme, so
! CCRad tuning knobs applied in glue_conv

!---------------------------------------------------------------------
! 10.2 Apply CCA_2D to 3d cloud profile
!---------------------------------------------------------------------

if (l_anvil) then

  ! Apply anvil scheme to deep cloud
  call calc_3d_cca                                                             &
    ( n_dp, n_dp, nlev, n_cca_lev, nbl, iccb, icct                             &
    , p_layer_boundaries, freeze_lev, cca_2d, cca, z_theta, z_rho )

  ! NOTE: iccb, icct are layer centres (theta levels) at this
  !        point.

else

  ! Apply cca_2d to all levels from deep base to deep top
  do i=1, n_dp
    ! Only copy across to cca if both iccb(i) and icct(i) have
    ! been set and are within cca array bounds
    if ((iccb(i) > 0) .and. (icct(i) > 0) .and.                                &
        (iccb(i) <= n_cca_lev)) then
      do k=iccb(i), min(icct(i), n_cca_lev)
        cca(i,k) = cca_2d(i)
      end do
    end if
  end do
end if      ! l_anvil

! l_dcpl_cld4pc2 is set to true for 5a scheme, so
! CCRad tuning knobs applied in glue_conv


!-----------------------------------------------------------------------
! Final SCM convection sub-step diagnostics
!-----------------------------------------------------------------------
if ( l_scm_convss_dg ) then

  ! 3-D diagnostics
  do k = 1, nlev
    do i = 1, n_dp
      if ( k>=ntml(i) .and. k<=det_lev(i) ) then
        ! Save the final updraft mass-flux
        scm_convss_dg(i) % up_flx(k) = flx(i,k)
      end if
    end do
  end do

  ! 2-D diagnostics
  do i = 1, n_dp
    scm_convss_dg(i) % precip_deep = rain(i) + snow(i)
  end do

end if


!-----------------------------------------------------------------------
! 11.0  End Subroutine
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine deep_conv_6a
end module deep_conv_6a_mod
