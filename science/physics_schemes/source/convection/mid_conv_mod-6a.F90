! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Mid-level convection scheme

module mid_conv_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Mid level convection scheme
!
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


character(len=*), parameter, private :: ModuleName='MID_CONV_6A_MOD'

contains

subroutine  mid_conv_6a(nbl,nlev,ntra,n_wtrac,n_cca_lev,npnts,trlev,           &
                       bland, w_max,                                           &
                       exner_rho,                                              &
                       exner_layer_centres,                                    &
                       exner_layer_boundaries,                                 &
                       l_q_interact,                                           &
                       l_tracer,midtrig,ntml,ntpar,freeze_lev,                 &
                       pstar,p_layer_centres,p_layer_boundaries,               &
                       z_theta, z_rho,                                         &
                       r_theta, r_rho,                                         &
                       rho_theta, rho,                                         &
                       rho_dry_theta, rho_dry,                                 &
                       r2rho_th, r2rho,                                        &
                       dr_across_th, dr_across_rh,                             &
                       conv_prog_precip,                                       &
                       q,th,                                                   &
                       timestep,ftl,fqt,u,v,w,recip_pstar,qse,                 &
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
                       lcbase,lctop,rain,snow,rain_3d,snow_3d,                 &
                       up_flux,up_flux_half,                                   &
                       dwn_flux,l_mid_all,cca_2d,                              &
                       uw_mid,vw_mid, cfl_limited, dt_dd, dq_dd,               &
                       area_ud,                                                &
                       error_point,tnuc_new                                    &
                       )


use water_constants_mod, only: lc

use cv_derived_constants_mod, only: ls

use cv_run_mod, only:                                                          &
    l_mom, l_eman_dd, cldbase_opt_md, cape_min,                                &
    w_cape_limit, cape_timescale, mid_cmt_opt,                                 &
    mid_cnv_pmin, cca2d_md_opt,                                                &
    l_anvil, l_cv_conserve_check,                                              &
    l_cmt_heating, tice, l_fcape, l_prog_pert,                                 &
    l_conv_prog_flx, tau_conv_prog_flx, cape_ts_min, cape_ts_max,              &
    cnv_cold_pools, l_ccp_parcel_md, ccp_buoyancy,                             &
    md_pert_opt, efrac, max_pert_scale, thpixs_mid

use cv_param_mod, only:                                                        &
    total_condensed_water, srf_precip, avg_mass_flux,                          &
    a_land, a_sea, b_land, b_sea, w_cca, xsbmin,                               &
    qpixs_mid,                                                                 &
    mparb, c_mid, d_mid, wcape_fac, delthst,                                   &
    a_cape, b_cape, term_undil, ccp_off, ccp_h_coef,                           &
    md_pert_orig, md_pert_efrac, md_pert_bowen, refqsat,                       &
    min_pert_scale, logp_min, logp_max

use cv_dependent_switch_mod, only:                                             &
    md_on, md_new_termc

use mphys_inputs_mod, only: l_progn_tnuc

use cv_stash_flg_mod, only:                                                    &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                                  &
    flg_entr_up, flg_detr_up, flg_entr_dwn, flg_detr_dwn,                      &
    flg_mf_midlev, flg_area_ud

use scm_convss_dg_mod, only: scm_convss_dg_type

use conversions_mod, only: rsec_per_day

use planet_constants_mod, only: cp, c_virtual, g

use yomhook,    only: lhook, dr_hook
use parkind1,   only: jprb, jpim
use umPrintMgr, only: umPrint, ummessage, printstatus, prstatus_normal

! subroutines
use lift_par_6a_mod,       only: lift_par_6a
use lift_undil_par_6a_mod, only: lift_undil_par_6a
use convec2_6a_mod,        only: convec2_6a
use water_loading_mod,     only: water_loading
use cor_engy_6a_mod,       only: cor_engy_6a
use cmt_heating_mod,       only: cmt_heating
use layer_cn_6a_mod,       only: layer_cn_6a, midlevel

use calc_3d_cca_mod,      only: calc_3d_cca
use dd_call_6a_mod,       only: dd_call_6a
use eman_cex_mod,         only: eman_cex
use evap_bcb_nodd_6a_mod, only: evap_bcb_nodd_6a
use flag_wet_mod,         only: flag_wet
use mid_conv_dif_cmt_mod, only: mid_conv_dif_cmt

use wtrac_conv_mod,           only: l_wtrac_conv,                              &
                                    conv_e_wtrac_type, conv_p_wtrac_type,      &
                                    wtrac_alloc_conv_p, wtrac_dealloc_conv_p
use correct_small_q_conv_mod, only: correct_small_q_conv
use water_tracers_mod,        only: wtrac_info
use wtrac_calc_ratio_mod,     only: wtrac_calc_ratio_fn

implicit none

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent in:

integer, intent(in) ::                                                         &
  nbl                  & ! No. of boundary layer levels
 ,nlev                 & ! No. of model layers
 ,ntra                 & ! No. of tracer fields
 ,n_wtrac              & ! No. of water tracer fields
 ,n_cca_lev            & ! No. of convective cloud amount levels (1 for 2D,
                         ! nlevs for 3D)
 ,npnts                & ! No. of deep convection points
 ,trlev                  ! No. of model levels on which tracers are included

logical, intent(in) :: bland(npnts) ! Land/sea mask

real(kind=real_umphys), intent(in)    :: exner_rho(npnts,nlev)
                                              ! Exner on rho levels

real(kind=real_umphys), intent(in)    ::                                       &
  exner_layer_centres(npnts,0:nlev)    & ! Exner
 ,exner_layer_boundaries(npnts,0:nlev)   ! Exner at half level above
                                         ! exner_layer_centres

logical, intent(in) ::                                                         &
  l_q_interact          & ! Switch allows overwriting parcel variables when
                          ! calculating condensate incr.
 ,l_tracer                ! Switch for inclusion of tracers

integer, intent(in) ::                                                         &
  midtrig(npnts)       & ! Lowest trigger levelfor convection
 ,ntml(npnts)          & ! Top level of surface mixed layer defined relative to
                         ! theta,q grid
 ,ntpar(npnts)         & ! Top level of initial parcel ascent in BL scheme
                         ! defined relative to theta,q grid
 ,freeze_lev(npnts)      ! freezing level

real(kind=real_umphys), intent(in)    ::                                       &
  pstar(npnts)                     & ! Surface pressure (Pa)
 ,p_layer_centres(npnts,0:nlev)    & ! Pressure (Pa)
 ,p_layer_boundaries(npnts,0:nlev)   ! Pressure at half level above
                                     ! p_layer_centres (Pa)

real(kind=real_umphys), intent(in) :: z_theta(npnts,nlev)
                                             ! height of theta levels (m)
real(kind=real_umphys), intent(in) :: z_rho(npnts,nlev)
                                             ! height of rho levels (m)
real(kind=real_umphys), intent(in) :: r_theta(npnts,0:nlev)
                                             ! radius of theta levels (m)
real(kind=real_umphys), intent(in) :: r_rho(npnts,nlev)
                                             ! radius of rho levels (m)
real(kind=real_umphys), intent(in) :: rho_theta(npnts,nlev)
                                             ! density for theta lev (kg/m3)
real(kind=real_umphys), intent(in) :: rho(npnts,nlev)
                                             ! density for rho lev (kg/m3)
real(kind=real_umphys), intent(in) ::                                          &
  rho_dry_theta(npnts,nlev) & ! dry density on theta levels (kg/m3)
 ,rho_dry(npnts,nlev)         ! dry density on rho levels (kg/m3)

real(kind=real_umphys), intent(in) :: r2rho_th(npnts,nlev)
                                             ! radius**2 density for
                                             ! theta lev (kg/m)
real(kind=real_umphys), intent(in) :: r2rho(npnts,nlev)
                                             ! radius**2 density for
                                             ! rho lev (kg/m)
real(kind=real_umphys), intent(in) :: dr_across_th(npnts,nlev)
                                             ! thickness of theta levels (m)
real(kind=real_umphys), intent(in) :: dr_across_rh(npnts,nlev)
                                             ! thickness of rho levels (m)
real(kind=real_umphys), intent(in) :: conv_prog_precip(npnts,nlev)
                                                ! Surface precipitation based
                                                ! 3d convective prognostic in
                                                ! kg/m2/s
real(kind=real_umphys), intent(in) ::  ftl(npnts)
                                        ! Surface sensible heat flux divided
                                        ! by cp (K kg/m2/s)
real(kind=real_umphys), intent(in) ::  fqt(npnts)
                                        ! Surface total water flux (kg/m2/s)

real(kind=real_umphys), intent(in)    ::                                       &
  q(npnts,nlev)         & ! Model mixing ratio (kg/kg)
 ,th(npnts,nlev)        & ! Model potential temperature (K)
 ,timestep              & ! Model timestep (s)
 ,u(npnts,nlev)         & ! Model u field (m/s)
 ,v(npnts,nlev)         & ! Model v field (m/s)
 ,w(npnts,nlev)         & ! Model w field (m/s)
 ,w_max(npnts)          & ! max w in column
                          ! for use in scale dependent cape timescale
 ,recip_pstar(npnts)    & ! Reciprocal of pstar array
 ,qse(npnts,nlev)         ! Saturation mixing ratio of
                          ! cloud environment (kg/kg)

real(kind=real_umphys), intent(in)    ::   g_ccp(npnts)                        &
                                              ! Cold_pool reduced gravity
                       , h_ccp(npnts) &       ! Cold-pool depth
                       , ccp_strength(npnts)  ! Cold-pool strength

! Arguments with intent INOUT:

real(kind=real_umphys), intent(in out) ::                                      &
  cf_frozen(npnts,nlev)  & ! Frozen water cloud volume ( )
 ,cf_liquid(npnts,nlev)  & ! Liq water cloud volume ( )
 ,qcf(npnts,nlev)        & ! Ice condensate mix ratio (kg/kg)
 ,qcl(npnts,nlev)          ! Liq condensate mix ratio (kg/kg)

real(kind=real_umphys), intent(in out) :: tracer(npnts,trlev,ntra)
                                                 !Model tracer fields (kg/kg)

real(kind=real_umphys), intent(in out) :: w2p(npnts,nlev) ! (Parcel vertical velocity)^2 [(m/s)^2]

real(kind=real_umphys), intent(in out) :: conv_prog_flx(npnts,nlev)
                                                 ! Mass flux convective
                                                ! prognostic in Pa/s

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                               ! Structure containing water
                                               ! tracer fields

! Structure containing SCM convection sub-step diagnotics
! (needs intent inout as contains allocatable arrays that need to
! retain their allocated status on input as well as output)
type(scm_convss_dg_type), intent(in out) :: scm_convss_dg( npnts )
! Flag for SCM convection sub-step diagnostics
logical, intent(in) :: l_scm_convss_dg


! Arguments with intent out:

real(kind=real_umphys), intent(out) ::                                         &
  cape_out(npnts)      & ! Saved convective available
                         ! potential energy for diagnostic output (J/kg)
 ,cclwp(npnts)         & ! Condensed water path (kg/m^2)
 ,ccw(npnts,nlev)      & ! Convective cloud liquid water
                         ! on model levels (kg/kg)
 ,cca(npnts,n_cca_lev) & ! Convective cloud amount on model levels (0-1)
 ,dbcfbydt(npnts,nlev) & ! Increments to total cld volume due to convection(/s)
 ,dcffbydt(npnts,nlev) & ! Increments to ice cloud volume due to convection(/s)
 ,dcflbydt(npnts,nlev) & ! Increments to liq cloud volume due to convection (/s)
 ,dqbydt(npnts,nlev)   & ! Increments to q due to convection (kg/kg/s)
 ,dqcfbydt(npnts,nlev) & ! Increments to ice condensate due to convection
                         ! (kg/kg/s)
 ,dqclbydt(npnts,nlev) & ! Increments to liq condensate due to convection
                         ! (kg/kg/s)
 ,dthbydt(npnts,nlev)  & ! Increments to potential temp. due to convection (K/s)
 ,dubydt(npnts,nlev) & ! Increments to U due to CMT (m/s2)
 ,dvbydt(npnts,nlev)   ! Increments to V due to CMT (m/s2)

real(kind=real_umphys), intent(out) ::                                         &
  dtrabydt(npnts,nlev,ntra)   !Increment to tracer convection (kg/kg/s)

real(kind=real_umphys), intent(out) ::                                         &
  detrain_up(npnts,nlev)   & ! Fractional detrainment rate into updraughts
                             ! (Pa/s)
 ,detrain_dwn(npnts,nlev)  & ! Fractional detrainment rate into downdraughts
                             ! (Pa/s)
 ,entrain_up(npnts,nlev)   & ! Fractional entrainment rate into updraughts
                             ! (Pa/s)
 ,entrain_dwn(npnts,nlev)    ! Fractional entrainment rate into downdraughts
                             ! (Pa/s)

integer, intent(out) ::                                                        &
  iccb(npnts)           & ! Convective cloud base level
 ,icct(npnts)             ! Convective cloud top level

real(kind=real_umphys), intent(out) :: lcca(npnts) ! Lowest conv. cloud amt. (%)

integer, intent(out) ::                                                        &
  lcbase(npnts)         & ! Lowest conv. cloud base level
 ,lctop(npnts)            ! Lowest conv. cloud top level

real(kind=real_umphys), intent(out) ::                                         &
  rain(npnts)         & ! Surface convective rainfall (kg/m2/s)
 ,snow(npnts)         & ! Surface convective snowfall (kg/m2/s)
 ,rain_3d(npnts,nlev) & ! Convective rainfall flux (kg/m2/s)
 ,snow_3d(npnts,nlev) & ! Convective snowfall flux (kg/m2/s)
 ,up_flux(npnts,nlev) & ! Updraught mass flux (Pa/s)
 ,up_flux_half(npnts,nlev) & ! Updraught mass flux on half levels(Pa/s)
 ,dwn_flux(npnts,nlev)  ! Downdraught mass flux (Pa/s)

logical, intent(out) ::                                                        &
  l_mid_all(npnts)        ! Points where mid level convection
                          ! occurs at some level
real(kind=real_umphys),intent(in out) ::                                       &
  cca_2d(npnts)         !2D convective cloud amount(%)

! CMT diagnostics
real(kind=real_umphys), intent(out) ::                                         &
  uw_mid(npnts,nlev)      & ! U component of stress from mid-level convection
                            ! (kg/m/s2)
 ,vw_mid(npnts,nlev)        ! V component of stress from mid-level convection
                            ! (kg/m/s2)
real(kind=real_umphys), intent(out) ::                                         &
  cfl_limited(npnts)        ! Inidicator of CFL limited mid-level convection

real(kind=real_umphys), intent(out) ::                                         &
  dt_dd(npnts,nlev)       & ! dT/dt from DD & evap below cloud base (K/s)
 ,dq_dd(npnts,nlev)       & ! dq/dt from DD & evap below cloud base (kg/kg/s)
 ,area_ud(npnts,nlev)       ! fractional area of updraughts

integer, intent(out) ::                                                        &
  error_point               ! 0 no problem, > 0 problem point

real(kind=real_umphys), intent(in) :: tnuc_new(npnts,nlev)
                                         ! 3D ice nucleation temperature as
                                         ! function of dust

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

! Height above surface of model levels ...
real(kind=real_umphys) :: zkm1(npnts)             ! ... k-1   [m]

real(kind=real_umphys) :: entrain_coef(npnts)
                                ! entrainment coefficients unset

integer :: index1(npnts),index2(npnts)

integer :: ncposs               ! No. of points which may convect

integer :: nconv                ! No. of convecting points

real(kind=real_umphys) :: amdetk(npnts)
                                ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

real(kind=real_umphys) :: cape(npnts)
                                ! Convective available potential
                                ! energy (J/kg)
real(kind=real_umphys) :: fcape(npnts)
                                ! Convective available potential
                                ! energy weighted by f.det profile (J/kg)

real(kind=real_umphys) :: dcpbydt(npnts)
                                ! Rate of change of cape (J/kg/s)

real(kind=real_umphys) :: depth(npnts)
                                ! Depth of convective cloud (m)

real(kind=real_umphys) :: eminds(npnts)
                                ! Minimum buoyancy for convection
                                ! to initiate from level k
                                ! (Kelvin)

real(kind=real_umphys) :: ekp14(npnts)            ! Entrainment coefficients at
                                ! level k+1/4 multiplied by
                                ! appropriate layer thickness
                                !(Dimensionless)

real(kind=real_umphys) :: ekp34(npnts)            ! Entrainment coefficients at
                                ! level k+3/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

real(kind=real_umphys) :: exk(npnts)              ! Exner ratio at layer k

real(kind=real_umphys) :: exkp1(npnts)            ! Exner ratio at layer k+1

real(kind=real_umphys) :: flxmax(npnts)           ! Maximum initial convective
                                ! mass flux (Pa/s)

real(kind=real_umphys) :: flx_init(npnts)
                                ! Initial mass flux at cloud base (Pa/s)

real(kind=real_umphys) :: flx_init_new(npnts)
                                ! flx_init scaled to destroy cape
                                ! over timescale cape_timescale (Pa/s)

real(kind=real_umphys) :: flxmax_init(npnts)
                                ! Maximum possible initial mass
                                ! flux (limited to the mass in
                                ! the initial convecting layer in Pa/s)

real(kind=real_umphys) :: flxbyrho_int(npnts)
                                ! mass weighted vertical integral of
                                ! mass flux/rho used in cca_2d calculation
                                ! (m/s)

real(kind=real_umphys) :: tcw(npnts)
                                ! Total condensed water(kg/m2/s)

real(kind=real_umphys) :: tot_mass_mean(npnts)
                                ! total mass of all convecting layers (kg/m2)

real(kind=real_umphys) :: max_cfl(npnts)
                                ! Max cfl ratio over a convecting
                                ! layer

real(kind=real_umphys) :: decay_amount
                                ! decay fraction for time-smoothed mass flux

real(kind=real_umphys) :: tmp_conv_prog_flx(npnts,nlev)
                                      ! Copy of Mass flux convective
                                ! prognostic in Pa/s

real(kind=real_umphys) :: precip(npnts,nlev)
                                ! Amount of precip from each layer
                                ! from each layer (kg/m/s)

real(kind=real_umphys) :: pk(npnts)
                                ! Pressure at midpoint of layer k (Pa)

real(kind=real_umphys) :: pkp1(npnts)
                                ! Pressure at midpoint of layer k+1 (Pa)

real(kind=real_umphys) :: delpk(npnts)
                                ! Pressure difference over layer k (Pa)

real(kind=real_umphys) :: delpkp1(npnts)
                                ! Pressure difference over layer k+1 (Pa)

real(kind=real_umphys) :: delpkp12(npnts)         ! Pressure difference between
                                ! layers k and k+1 (Pa)

real(kind=real_umphys) :: delp_uv_k(npnts)
                                ! Pressure difference across uv layer k (Pa)

real(kind=real_umphys) :: delp_uv_kp1(npnts)
                                ! Pressure difference across uv layer k+1 (Pa)

real(kind=real_umphys) :: rhum(npnts)             ! Dummy relative humidity
                                ! (only used on shallow points)

real(kind=real_umphys) :: wsc_o_mb(npnts)         ! Dummy argument for layer_cn
                                ! Convective velocity scale divided
                                ! by cloud base mass flux mb

logical :: bgmk(npnts)          ! Mask for points where parcel in
                                ! layer k is saturated
logical :: blatent(npnts)       ! Mask for points where latent heat has
                                ! been released

logical :: bwater(npnts,nlev)   ! Mask for points at which
                                ! condensate is liquid

logical :: bterm(npnts)         ! Mask for points which have
                                ! stopped convecting

logical :: bconv(npnts)         ! Mask for points at which
                                ! convection is occurring

logical :: bcposs(npnts)        ! Mask for points passing
                                ! initial stability test

logical :: binit(npnts)         ! Mask for points initiating

! Parcel variables

real(kind=real_umphys) :: qpi   ! Initial parcel mixing ratio (kg/kg)

real(kind=real_umphys) :: qp(npnts,nlev)          ! Parcel mixing ratio (kg/kg)

real(kind=real_umphys) :: thpi  ! Initial parcel potential temp.(K)

real(kind=real_umphys) :: thp(npnts,nlev)         ! Parcel potential temp (K)

real(kind=real_umphys) :: up(npnts,nlev)          ! Parcel U (m/s)

real(kind=real_umphys) :: vp(npnts,nlev)          ! Parcel V (m/s)

real(kind=real_umphys) :: trap(npnts,nlev,ntra)
                                ! Tracer content of parcel (kg/kg)

real(kind=real_umphys) :: flx(npnts,nlev)         ! Parcel massflux (Pa/s)

real(kind=real_umphys) :: xsbmin_v(npnts,nlev)
                                ! Minmum parcel buoyancy excess

real(kind=real_umphys) :: thpixs_v(npnts,nlev)    ! Theta parcel excess (K)

real(kind=real_umphys) :: qpixs_v(npnts,nlev)     ! Q parcel excess(kg/kg)

! Water tracers
type(conv_p_wtrac_type) :: wtrac_p(n_wtrac)
                                ! Structure containing water
                                ! tracer fields relating to the parcel
real(kind=real_umphys), allocatable :: qpixs_v_wtrac(:,:,:)
                                ! Water tracer q parcel excess (kg/kg)
real(kind=real_umphys) :: ratio_wt
                                !  Ratio of water tracer to normal water

real(kind=real_umphys) :: dpmin(npnts)
                                ! work array for parcel excess cal

real(kind=real_umphys) :: pert_scale(npnts,nlev)
                                ! Scaling of initial perturbation due to
                                ! precipitation based prognostic

real(kind=real_umphys) :: qclp(npnts,nlev)
                                 ! Parcel liquid condensate mixing
                                 ! ratio in layer k (kg/kg)

real(kind=real_umphys) :: qcfp(npnts,nlev)
                                 ! Parcel frozen condensate mixing
                                 ! ratio in layer k (kg/kg)

real(kind=real_umphys) :: el     ! Latent heat of vapour to liquid or ice
                                 ! (J/kg)
real(kind=real_umphys) :: denom  ! Denominator used in parcel perturbation
real(kind=real_umphys) :: tmp_efrac
                                 ! Temporary evaporative fraction used in
                                 ! parcel perturbation calculation
real(kind=real_umphys) :: wght_efrac
                                 ! Weighting between namelist efrac and
                                 ! surface flux based efrac


! Undilute parcel variables
real(kind=real_umphys) :: qu(npnts,nlev)
                                ! Undilute Parcel mixing ratio (kg/kg)
real(kind=real_umphys) :: thu(npnts,nlev)
                                ! Undilute Parcel potential temp (K)

! Parameters

real(kind=real_umphys), parameter :: cfl_limit = 1.0 ! Max CFL ratio allowed
real(kind=real_umphys), parameter :: minflx = tiny(flx_init_new)
                                                ! minimum allowable
                                                ! initial mass flux

! CMT variables

integer :: kterm(npnts)         ! Level index for termination of

real(kind=real_umphys) :: eflux_u_ud(npnts)
                                ! Vertical eddy flux of momentum
                                ! due to UD at top of layer (Pa m/s2)

real(kind=real_umphys) :: eflux_v_ud(npnts)
                                ! Vertical eddy flux of momentum
                                ! due to UD at bottom of layer (Pa m/s2)

logical :: l_mom_gk             ! true if Gregory-Kershaw CMT
logical :: l_mom_gk_stable      ! true if stabilized Gregory-Kershaw CMT
                                ! (different from the original)

! Cape scaling/closure variables

integer :: start_lev(npnts)     ! Convection initiation level

integer :: det_lev(npnts)       ! Level at which split final
                                ! detrainment last occurred

integer :: nterm                ! No. of points where conv.
                                ! has terminated

integer :: index_nterm(npnts)   ! Index for points where conv.
                                ! has terminated

real(kind=real_umphys) :: tempnum
                                ! Temporary variable for storage

real(kind=real_umphys) :: scale_f(npnts)          ! scaling factor
real(kind=real_umphys) :: scale_low(npnts)
                                ! scaling factor lowest mid-level

real(kind=real_umphys) :: cape_ts_new(npnts)      ! Used as variable in RH-based
                                ! closure

real(kind=real_umphys) :: relh(npnts)             ! RH integral (average when
                                ! convection terminates)

real(kind=real_umphys) ::                                                      &
  wls_mean(npnts)             & ! large-scale w over convecting column
 ,mass_mean(npnts)              ! mass over convecting column

real(kind=real_umphys) :: dptot(npnts)            ! Delta P integral
real(kind=real_umphys) :: deltaktot(npnts)
                                ! Integrated forced detrainment
! Original downdraught scheme variables

integer :: npossdd              ! Max. no. of downdraughts
                                ! possible

integer :: nnodd                ! No. of downdraughts not possible

integer :: index_possdd(npnts)  ! Index of downdraughts possible

integer :: index_nodd(npnts)    ! Index of downdraughts not
                                ! possible

real(kind=real_umphys) :: deltap_cld
                                ! Pressure thickness of convective
                                ! cloud (Pa)

! Arrays required by Emanuel downdraughts

integer ::                                                                     &
  kterm_mid(npnts)  & ! termination level of highest mid level
, kterm_max

logical ::                                                                     &
  bgmkp1(npnts)

real(kind=real_umphys) ::                                                      &
   qpkp1(npnts)                                                                &
  ,qclpkp1(npnts)                                                              &
  ,qcfpkp1(npnts)                                                              &
  ,thpkp1(npnts)

real(kind=real_umphys) :: Qlkp1(npnts)
real(kind=real_umphys) :: Qfkp1(npnts)
real(kind=real_umphys) :: Frezkp1(npnts)

real(kind=real_umphys) :: watldek(npnts)
real(kind=real_umphys) :: watldpk(npnts)
real(kind=real_umphys) :: watldekp1(npnts)
real(kind=real_umphys) :: watldpkp1(npnts)


! Local arrays

real(kind=real_umphys) :: thrk(npnts)
                             ! potential temperature of forced detrained air

real(kind=real_umphys) :: qrk(npnts)
                             ! specific humidity of forced detrained air

real(kind=real_umphys) :: deltak(npnts)        ! Parcel forced detrainment rate
                             ! in layer k multiplied by
                             ! appropriate layer thickness

real(kind=real_umphys) :: flxkp12_t(npnts)     ! Half level mass flux (Pa/s)

real(kind=real_umphys) :: rbuoyk(npnts)       ! Par. buoyancy at k (K)
real(kind=real_umphys) :: rbuoykp1(npnts)    ! Par. buoyancy at k+1 (K)
real(kind=real_umphys) :: rbuoyukp1(npnts)   ! undilute Par. buoy at k+1 (K)


real(kind=real_umphys) :: qsat_lcl(npnts)     !qsat at the LCL

real(kind=real_umphys) :: tnuc_nlcl(npnts)
                          ! Ice nucleation temperature at parcel starting level


!===============================================================
! CCRad Variables local variables
!===============================================================
integer          :: n_mdcld     ! Number of gridpoints with
                                ! single mid-level cloud banks
                                !
integer          :: n_mdcld_mult   ! Number of gridpoints with
                                   ! multiple mid-level cloud
                                   ! banks
                                   !

!===============================================================
! Allocatable arrays, because we do not know how many gridpoints
! will contain mid-level until we test for them.
!===============================================================

integer                            :: dum1(npnts)
integer                            :: dum2(npnts)

integer, allocatable :: mdcldi(:)
                                ! INDICES in full array of
                                ! gridpoints single mid-level
                                ! cloud banks
integer, allocatable :: mdcldi_mult(:)
                                ! index in full array of
                                ! gridpoints multiple mid-level
                                ! cloud banks

!===============================================================
! Compressed arrays for gridpoints which have one bank of
! mid-level cloud and multiple banks of mid-level cloud.
! Requires the use of compressed allocatable arrays because the
! location and number of gridpoints with mid-level has not been
! diagnosed yet.
!===============================================================

!===============================================================
! For multiple mid-level cloud
!===============================================================
integer, allocatable ::  iccb_md_c(:)
integer, allocatable ::  icct_md_c(:)
integer, allocatable ::  freeze_lev_md_c(:)
real(kind=real_umphys), allocatable ::  cca_2d_md_c(:)
real(kind=real_umphys), allocatable ::  cca_md_c(:,:)
real(kind=real_umphys), allocatable ::  ccw_md_c(:,:)
real(kind=real_umphys), allocatable ::  z_theta_md_c(:,:)
real(kind=real_umphys), allocatable ::  z_rho_md_c(:,:)
real(kind=real_umphys), allocatable ::  p_lyr_bnds_md_c(:,:)

!===============================================================
! For multiple mid-level cloud
!===============================================================
integer, allocatable ::  iccb_md_mult_c(:)
integer, allocatable ::  icct_md_mult_c(:)
integer, allocatable ::  freeze_lev_md_mult_c(:)
real(kind=real_umphys), allocatable ::  cca_2d_md_mult_c(:)
real(kind=real_umphys), allocatable ::  cca_md_mult_c(:,:)
real(kind=real_umphys), allocatable ::  ccw_md_mult_c(:,:)
real(kind=real_umphys), allocatable ::  z_theta_md_mult_c(:,:)
real(kind=real_umphys), allocatable ::  z_rho_md_mult_c(:,:)
real(kind=real_umphys), allocatable ::  p_lyr_bnds_md_mult_c(:,:)

!===============================================================
! End CCRad Variables local variables
!===============================================================

!   required by check on -ve q
real(kind=real_umphys), parameter :: qmin = 1.0e-8 ! Global minimum allowed Q

integer ::                                                                     &
 nmid             & ! total number of points where mid-level
                    ! convection occurs
,nmax_layers        ! Maximum number of allow mid-level layers

real(kind=real_umphys) ::                                                      &
  rh_test    &   ! critical RH value for cldbase_opt_md=6
 ,rh_fac         ! factor for cldbase_opt_md=6 calculation

! local versions of CAPE timescale max/min
real(kind=real_umphys) :: cape_ts_min_use
real(kind=real_umphys) :: cape_ts_max_use

! Loop counters
integer :: i,j,k,ktra,kt,kk,i_wt

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MID_CONV_6A'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
error_point = 0

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

! Initialise logicals

do i = 1,npnts
  bconv(i)     = .false.
  l_mid_all(i) = .false.
  blatent(i)   = .false.
end do

!-----------------------------------------------------------------------
! Decide whether Gregory-Kershaw CMT scheme is to be used if l_mom is true.
! The Gregory-Kershaw scheme requires calculations in the main plume ascent
! loop whereas the alternative diffusive scheme is called after the plume
! calculation.

if (mid_cmt_opt == 1 ) then
  l_mom_gk = .false.       ! Use diffusive scheme
else
  l_mom_gk = l_mom         ! Use Gregory-Kershaw scheme
end if

if (l_mom .and. mid_cmt_opt == 2 ) then  ! Stabilized Gregory-Kershaw scheme
  l_mom_gk_stable = .true.
else                                     ! Original Gregory-Kershaw scheme
  l_mom_gk_stable = .false.
end if

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!-----------------------------------------------------------------------
do i = 1,npnts
  kterm(i)       = 0
  kterm_mid(i)   = 0
  scale_low(i)   = 1.0   ! scaling factor for DD
end do


!initialise parcel values over all levels
do k = 1, nlev
  do i = 1, npnts
    qp(i,k)     = 0.0
    thp(i,k)    = 0.0
    qclp(i,k)   = 0.0
    qcfp(i,k)   = 0.0
    flx(i,k)    = 0.0
    precip(i,k) = 0.0
    ccw(i,k)    = 0.0
  end do
end do

if (md_new_termc == term_undil) then
  do k = 1, nlev
    do i = 1, npnts
      qu(i,k)   = 0.0
      thu(i,k)  = 0.0
    end do
  end do
end if

if (l_mom_gk) then
  do k=1,nlev
    do i = 1,npnts
      up(i,k) = 0.0
      vp(i,k) = 0.0
    end do
  end do
end if

if (l_tracer) then
  do ktra = 1,ntra
    do k=1,nlev
      do i = 1,npnts
        trap(i,k,ktra) = 0.0
      end do
    end do
  end do
end if

! Allocate water tracer arrays and initialise

call wtrac_alloc_conv_p(npnts, nlev, n_wtrac, wtrac_p)

if (l_wtrac_conv) then

  do i_wt = 1, n_wtrac
    do k = 1, nlev
      do i = 1, npnts
        wtrac_p(i_wt)%q(i,k)      = 0.0
        wtrac_p(i_wt)%qcl(i,k)    = 0.0
        wtrac_p(i_wt)%qcf(i,k)    = 0.0
        wtrac_p(i_wt)%precip(i,k) = 0.0
      end do
    end do
  end do
end if  !l_wtrac_conv

do k=1,nlev
  do i=1, npnts
    dthbydt(i,k)  = 0.0
    dqbydt(i,k)   = 0.0
    dqclbydt(i,k) = 0.0
    dqcfbydt(i,k) = 0.0
    dbcfbydt(i,k) = 0.0
    dcflbydt(i,k) = 0.0
    dcffbydt(i,k) = 0.0
    dt_dd(i,k)=0.0
    dq_dd(i,k)=0.0
    area_ud(i,k) = 0.0
  end do
end do

if (l_mom) then
  do k = 1,nlev
    do i = 1,npnts
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    end do
  end do
  do k = 1,nlev
    do i = 1,npnts
      uw_mid(i,k) = 0.0
      vw_mid(i,k) = 0.0
    end do
  end do
end if  ! L_mom

if (l_tracer) then
  do ktra = 1,ntra
    do k = 1,nlev
      do i = 1,npnts
        dtrabydt(i,k,ktra) = 0.0
      end do
    end do
  end do
end if  ! L_tracer

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do k = 1, nlev
      do i = 1, npnts
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
if (flg_up_flx .or. flg_mf_midlev) then
  do k = 1,nlev
    do i = 1,npnts
      up_flux(i,k)      = 0.0
    end do
  end do
end if
if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,npnts
      up_flux_half(i,k) = 0.0
    end do
  end do
end if
if (flg_dwn_flx) then
  do k = 1,nlev
    do i = 1,npnts
      dwn_flux(i,k)     = 0.0
    end do
  end do
end if
! Now required in all cases
do k = 1,nlev
  do i = 1,npnts
    entrain_up(i,k)   = 0.0
  end do
end do
if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,npnts
      detrain_up(i,k)   = 0.0
    end do
  end do
end if
if (flg_entr_dwn) then
  do k = 1,nlev
    do i = 1,npnts
      entrain_dwn(i,k)  = 0.0
    end do
  end do
end if
if (flg_detr_dwn) then
  do k = 1,nlev
    do i = 1,npnts
      detrain_dwn(i,k)  = 0.0
    end do
  end do
end if

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
! Zeroes values to remove mid-level dependence on shallow and deep
! cloud
do i = 1,npnts
  iccb(i)   = 0
  icct(i)   = 0
  tcw(i)    = 0.0
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
end do

do k= 1,n_cca_lev
  do i= 1,npnts
    cca(i,k) = 0.0
  end do
end do


do i = 1,npnts
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
  cfl_limited(i)  = 0.0

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
  !      Initialise qsat at lcl for a parcel lifted from level 1.
  !      This will be equal to level 1 humidity.
  !-----------------------------------------------------------------------
  qsat_lcl(i)     = q(i,1)

  !-----------------------------------------------------------------------
  ! Initialise dummy variables that are not used by mid-level
  !-----------------------------------------------------------------------
  entrain_coef(i) = -99.0
  wsc_o_mb(i)     = 0.0
end do

! Initialise surface precip arrays for water tracers
if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i = 1,npnts
      wtrac_e(i_wt)%rain(i)   = 0.0
      wtrac_e(i_wt)%snow(i)   = 0.0
    end do
  end do
end if

! Initialise the cape timescale, limited by mid convection timestep.
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
! Set bwater=.true. on points where water will condense rather than
! ice.
!-----------------------------------------------------------------------

call flag_wet(npnts,npnts,nlev,th,exner_layer_centres,bwater(1,2))

! Set bwater at level 1 which is a special case not handled by flag_wet
do i = 1, npnts
  bwater(i,1) = th(i,1)*exner_layer_centres(i,1) > tice
end do

!-----------------------------------------------------------------------
! Enhance initial perturbations using precip based convective prognostic
!-----------------------------------------------------------------------
if (l_prog_pert) then
  denom = 1.0/(logp_max - logp_min)
  do k = 1,nlev-1
    do i = 1,npnts
      pert_scale(i,k) = denom * (max_pert_scale-min_pert_scale) *              &
                        (log10(conv_prog_precip(i,k) * refqsat/qse(i,k)) -     &
                        logp_min) + min_pert_scale
      pert_scale(i,k) = min(max(pert_scale(i,k), min_pert_scale),              &
                        max_pert_scale)
    end do
  end do  ! nlev
else
  do k = 1,nlev-1
    do i = 1,npnts
      pert_scale(i,k) = 1.0
    end do
  end do  ! nlev
end if


!-----------------------------------------------------------------------
! Calculate parcel perturbations
!-----------------------------------------------------------------------
select case (md_pert_opt)
case (md_pert_orig)
  ! Calculate xsbmin_v, thpixs_v and qpixs_v constants based on layer
  ! thickness (Pa).
  do k = 1,nlev-1
    do i = 1,npnts
      dpmin(i) = min( ((p_layer_centres(i,k) -                                 &
                         p_layer_centres(i,k+1))/5000.0),1.0)

      xsbmin_v(i,k) = dpmin(i) * thpixs_mid
      thpixs_v(i,k) = dpmin(i) * thpixs_mid
      qpixs_v(i,k)  = qpixs_mid
    end do
  end do  ! nlev
case (md_pert_efrac)
  ! Calculate xsbmin_v, thpixs_v and qpixs_v constants based on a proscribed
  ! ratio between latent heat flux and total heat flux at cloud base
  do k = 1,nlev-1
    do i = 1,npnts
      !Pressure thickness scaling
      dpmin(i) = min( ((p_layer_centres(i,k) -                                 &
                         p_layer_centres(i,k+1))/5000.0),1.0)

      !Calculate the buoyancy (thetav) perturbation
      xsbmin_v(i,k) = pert_scale(i,k) * dpmin(i) * thpixs_mid

      !Set latent heat dependent on phase
      if (bwater(i,k)) then
        el = lc
      else
        el = ls
      end if

      !Calculate the humidity perturbation
      denom = 1.0/(el*(1.0-efrac) +                                            &
              efrac*c_virtual*cp*th(i,k)*exner_layer_centres(i,k))
      qpixs_v(i,k)  = efrac*cp*exner_layer_centres(i,k)*xsbmin_v(i,k)*denom

      !Limit to humidity perturbation to avoid
      !(1) supersaturation of parcel (unless environment is supersaturated),
      !(2) perturbation being greater than 20% RH,
      !(3) being larger than 80% of environment humidity and
      !(4) perturbation being negative.
      qpixs_v(i,k)  = max(min(qpixs_v(i,k), qse(i,k)-q(i,k), 0.2*qse(i,k),     &
                          0.8*q(i,k)),0.0)

      !And update the theta perturbation such that the total buoyancy
      !perturbation is equal xsbmin_v
      thpixs_v(i,k) = xsbmin_v(i,k) - c_virtual*th(i,k)*qpixs_v(i,k)
    end do
  end do  ! nlev
case (md_pert_bowen)
  ! Calculate xsbmin_v, thpixs_v and qpixs_v constants based on the surface
  ! evaporative fraction if near the surface tapering to efrac away from surface
  do k = 1,nlev-1
    do i = 1,npnts
      !Pressure thickness scaling
      dpmin(i) = min( ((p_layer_centres(i,k) -                                 &
                         p_layer_centres(i,k+1))/5000.0),1.0)

      !Calculate the buoyancy (thetav) perturbation
      xsbmin_v(i,k) = pert_scale(i,k) *  dpmin(i) * thpixs_mid

      !Set latent heat dependent on phase
      if (bwater(i,k)) then
        el = lc
      else
        el = ls
      end if

      !Calculate the surface evaporative fraction
      if ((cp*ftl(i) + el*fqt(i)) >= 0.0) then
        !positive surface flux
        tmp_efrac = el*fqt(i)/(cp*ftl(i) + el*fqt(i))
      else
        tmp_efrac = efrac
      end if

      !Calculate z_ntml dependent weighting
      wght_efrac = min(1.0, exp(1.0-z_theta(i,k)/z_theta(i,ntml(i))))

      !Update the evaporative fraction to be applied at cumulus cloud base
      tmp_efrac = wght_efrac*tmp_efrac + (1.0-wght_efrac)*efrac

      !Calculate the humidity perturbation
      denom = 1.0/(el*(1.0-tmp_efrac) +                                        &
              tmp_efrac*c_virtual*cp*th(i,k)*exner_layer_centres(i,k))
      qpixs_v(i,k)  = efrac*cp*exner_layer_centres(i,k)*xsbmin_v(i,k)*denom

      !Limit to humidity perturbation to avoid
      !(1) supersaturation of parcel (unless environment is supersaturated),
      !(2) perturbation being greater than 20% RH,
      !(3) being larger than 80% of environment humidity and
      !(4) perturbation being negative.
      qpixs_v(i,k)  = max(min(qpixs_v(i,k), qse(i,k)-q(i,k), 0.2*qse(i,k),     &
                          0.8*q(i,k)),0.0)

      !And update the theta perturbation such that the total buoyancy
      !perturbation is equal xsbmin_v
      thpixs_v(i,k) = xsbmin_v(i,k) - c_virtual*th(i,k)*qpixs_v(i,k)
    end do
  end do  ! nlev
end select

if (l_wtrac_conv) then

  allocate(qpixs_v_wtrac(npnts,nlev,n_wtrac))

  do i_wt = 1, n_wtrac
    do k = 1,nlev-1
      do i = 1,npnts
        ! Calculate water tracer content of parcel q perturbation
        ratio_wt = wtrac_calc_ratio_fn(i_wt, wtrac_e(i_wt)%q(i,k), q(i,k))
        qpixs_v_wtrac(i,k,i_wt) = ratio_wt * qpixs_v(i,k)
      end do
    end do  ! nlev
  end do    ! n_wtrac
end if ! l_wtrac_conv

!-----------------------------------------------------------------------
! Additional cold-pool perturbation to theta.
!-----------------------------------------------------------------------
if ((cnv_cold_pools > ccp_off) .and. l_ccp_parcel_md) then
  do k = 1,nlev-1
    do i = 1,npnts
      thpixs_v(i,k) = thpixs_v(i,k) + dpmin(i) * ccp_buoyancy*g_ccp(i) *       &
                      exp(-ccp_h_coef*z_theta(i,k)/(1.0+h_ccp(i)))
    end do
  end do  ! nlev
end if


!-----------------------------------------------------------------------
! Main loop over levels
!-----------------------------------------------------------------------
do k = 1,nlev-1

  do i = 1,npnts
    bterm(i)   = .false.
  end do

  !-----------------------------------------------------------------------
  ! Initialise environment variables
  ! NB These variable are only used by layer_cn.
  !-----------------------------------------------------------------------
  do i = 1,npnts
    !Note that unlike p_layer_boundaries, where k indexing is offset
    !by one compared to the dynamics numbering, z retains the numbering
    !convention for dynamics variables i.e. for theta levels, k->k
    !and for rho levels k+1/2 -> k+1
    rhum(i)   = q(i,k) / qse(i,k)
  end do

  if (k == 1) then
    do i = 1,npnts
      zkm1(i) = 0.0
    end do
  else
    do i = 1,npnts
      zkm1(i) = z_theta (i,k-1)
    end do
  end if

  !-----------------------------------------------------------------------
  ! Set initial parcel properties (theta,q,tracer,momentum) if convection
  ! is not occurring at level k
  !-----------------------------------------------------------------------
  do i = 1,npnts
    if ( .not. bconv(i)) then
      bgmk(i)     = .false.
      depth(i)    = 0.0
      thpi        = th(i,k) + thpixs_v(i,k)
      thp(i,k)    = thpi
      qpi         = q(i,k)  + qpixs_v(i,k)
      qp(i,k)     = qpi
      start_lev(i)= k
      if (md_new_termc == term_undil) then
        thu(i,k)  = thpi
        qu(i,k)   = qpi
      end if
      if (l_q_interact) then
        qclp(i,k) = qcl(i,k)
        qcfp(i,k) = qcf(i,k)
      end if
      if (l_mom_gk) then
        up(i,k)   = u(i,k)
        vp(i,k)   = v(i,k)
      end if
    end if  !not bconv
  end do
  if (l_tracer) then
    do ktra=1,ntra
      do i = 1,npnts
        if ( .not. bconv(i)) then
          trap(i,k,ktra)  = tracer(i,k,ktra)
        end if  !not bconv
      end do    ! npnts
    end do
  end if

  ! Water tracers
  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      do i = 1,npnts
        if ( .not. bconv(i) ) then

          wtrac_p(i_wt)%q(i,k) = wtrac_e(i_wt)%q(i,k) + qpixs_v_wtrac(i,k,i_wt)

          if (l_q_interact) then
            wtrac_p(i_wt)%qcl(i,k) = wtrac_e(i_wt)%qcl(i,k)
            wtrac_p(i_wt)%qcf(i,k) = wtrac_e(i_wt)%qcf(i,k)
          end if
        end if
      end do
    end do
  end if   ! l_wtrac_conv

  ! Initialise binit i.e. set as no convection initialised in layer
  ! at start of this levels calculations
  do i = 1,npnts
    binit(i) = .false.
  end do

  !-----------------------------------------------------------------------
  ! 3.1  Calculate layer dependent constants (pressure,
  !      layer thickness, entrainment coefficients, detrainment
  !      coefficients)
  !-----------------------------------------------------------------------

  call layer_cn_6a(k, npnts, nlev,                                             &
                   ntml, ntpar, start_lev,                                     &
                   exner_layer_centres,                                        &
                   p_layer_boundaries, p_layer_centres,                        &
                   z_rho,                                                      &
                   conv_prog_precip,                                           &
                   recip_pstar, entrain_coef, rhum,                            &
                   ccp_strength,                                               &
                   z_theta(:,k), z_rho(:,k+1), z_theta(:,k+1),                 &
                   wsc_o_mb, qsat_lcl, w_max,                                  &
                   midlevel,                                                   &
                   bconv,                                                      &
                   ! Out
                   pk, pkp1, exk, exkp1,                                       &
                   delpk, delpkp12, delpkp1,                                   &
                   delp_uv_k, delp_uv_kp1,                                     &
                   ekp14, ekp34, amdetk                                        &
                   )

  do i = 1,npnts
    ! Maximum initial convective mass flux
    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)

    ! Minimum buoyancy for convection to start from level k
    eminds(i) = mparb * delpkp12(i) * recip_pstar(i)
  end do

  !-----------------------------------------------------------------------
  ! Initial test to check if convection is possible in layer k
  !-----------------------------------------------------------------------
  ! Convection is possible from layer k to k+1 if it is already convecting
  ! (bconv) or the all the following conditions are met:
  ! 1. it is at the triggering level (midtrig) or above. midtrig is set in
  !    glue_conv
  ! 2. it is not too high up (p>mid_cnv_pmin)
  ! 3. level k of the environment is unstable or only just stable
  !    (delthst) with respect to level k+1
  !-----------------------------------------------------------------------
  do i = 1,npnts
    bcposs(i) =   bconv(i)                                  .or.               &
                ( k >= midtrig(i)                           .and.              &
                  p_layer_centres(i,k) > mid_cnv_pmin       .and.              &
                ( th(i,k)   * (1.0 + c_virtual * q(i,k)) + delthst  +          &
                  max(0.0,(q(i,k)-qse(i,k+1))) * lc/(cp * exkp1(i)) >          &
                  th(i,k+1) * (1.0 + c_virtual * q(i,k+1)) ) )
  end do  ! npnts

  ! Calculate number of points which may convect (ncposs) and
  ! set compression indices (index1)
  ncposs = 0
  do i = 1,npnts
    if (bcposs(i)) then
      ncposs          = ncposs + 1
      index1(ncposs)  = i
    end if
  end do

  !-----------------------------------------------------------------------
  ! 3.2  Lift parcel from layer k to layer k+1
  !-----------------------------------------------------------------------
  if (ncposs > 0) then

    call lift_par_6a (k, npnts, n_wtrac, th(:,k), th(:,k+1),                   &
                q(:,k), q(:,k+1), qcl(:,k), qcf(:,k),                          &
                qcl(:,k+1), qcf(:,k+1),                                        &
                pkp1(:), exkp1(:),                                             &
                thp(:,k), qp(:,k), qclp(:,k),qcfp(:,k),                        &
                ekp14(:), ekp34(:),                                            &
                l_q_interact, bwater(:,k+1), wtrac_e,                          &
                !InOut
                wtrac_p,                                                       &
                !Out
                bgmkp1, thpkp1, qpkp1,                                         &
                qclpkp1, qcfpkp1,                                              &
                Qlkp1, Qfkp1, Frezkp1, index1, ncposs)


    !-----------------------------------------------------------------------
    ! Calculate the water loading for level k and k+1
    !-----------------------------------------------------------------------
    call water_loading(npnts, qcl(:,k), qcf(:,k), qclp(:,k), qcfp(:,k),        &
                     watldek, watldpk, index1, ncposs)

    call water_loading(npnts, qcl(:,k+1), qcf(:,k), qclpkp1, qcfpkp1,          &
                     watldekp1, watldpkp1, index1, ncposs)


    if (md_new_termc == term_undil) then
      !---------------------------------------------------------------------
      ! Calculate undilute parcel properties
      !---------------------------------------------------------------------
      call lift_undil_par_6a(npnts,                                            &
                  pkp1(:), exkp1(:),                                           &
                  thu(:,k), qu(:,k),                                           &
                  Frezkp1(:), bwater(:,k+1),                                   &
                  !Out
                  thu(:,k+1), qu(:,k+1), index1, ncposs)
    end if

    !-----------------------------------------------------------------------
    ! Test if convection is starting from layer k
    ! And set initial mass flux
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

      if (md_new_termc == term_undil) then
        !undilute parcel buoyancy
        rbuoyukp1(index1(i))= thu(index1(i),k+1)                               &
                      * (1.0 + c_virtual *qu(index1(i), k+1)                   &
                      - watldpkp1(index1(i)))                                  &
                      - th(index1(i),k+1) * (1.0 + c_virtual *q(index1(i),k+1) &
                      - watldekp1(index1(i)))
      end if

    end do

!DIR$ IVDEP
    do i = 1,ncposs ! Loop over points which may convect
      !-----------------------------------------------------------------------
      ! Allow parcel to start convection at midtrig or above if
      ! 1. it does not overlap with a lower event
      ! 2. it is not already convecting
      ! 3. it is buoyant at the next level.
      !-----------------------------------------------------------------------
      if (k >= midtrig(index1(i))   .and.                                      &
          k > det_lev(index1(i))    .and.                                      &
          .not. bconv(index1(i))    .and.                                      &
          rbuoykp1(index1(i))  > (eminds(index1(i)) + xsbmin)) then
        start_lev(index1(i))  =  k
        bconv(index1(i))      = .true.
        binit(index1(i))      = .true.

        ! Set parcel mass flux (UM Documentation paper 27, section 1.5)
        flx(index1(i),k)   = 1.0e-3*pstar(index1(i)) * (d_mid + c_mid *        &
                      pstar(index1(i)) * ((rbuoykp1(index1(i)) - xsbmin)       &
                      / delpkp12(index1(i))))

        ! If mass flux out of the initial layer is greater than the mass flux
        ! of the layer over the timestep then limit mass flux to mass of layer.
        if (flx(index1(i),k) > flxmax(index1(i))) then
          flx(index1(i),k) = flxmax(index1(i))
        end if

        ! Apply the parcel perturbation increments
        dthbydt(index1(i),k) = dthbydt(index1(i),k)                            &
                               - flx(index1(i),k)/delpk(index1(i))             &
                               *(thp(index1(i),k)-th(index1(i),k))
        dqbydt(index1(i),k)  = dqbydt(index1(i),k)                             &
                               - flx(index1(i),k)/delpk(index1(i))             &
                               *(qp(index1(i),k)-q(index1(i),k))

        if (l_wtrac_conv) then
          do i_wt = 1, n_wtrac
            wtrac_e(i_wt)%dqbydt(index1(i),k) =                                &
                  wtrac_e(i_wt)%dqbydt(index1(i),k)                            &
                  - flx(index1(i),k)/delpk(index1(i))                          &
                  *(wtrac_p(i_wt)%q(index1(i),k)-wtrac_e(i_wt)%q(index1(i),k))
          end do
        end if

        ! Store diagnostics linked to initial convective mass flux for
        ! calculation of final closure.
        flx_init(index1(i))   = flx(index1(i),k)
        flxmax_init(index1(i))= flxmax(index1(i))

        if (l_progn_tnuc) then
          ! Set ice nucleation temperature at starting level
          tnuc_nlcl(index1(i)) = tnuc_new(index1(i),k)
        end if

      end if

    end do  !ncposs loop

  end if  !ncposs>0

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
    !
    !      Subroutine CONVEC2
    !
    !      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
    !-----------------------------------------------------------------------

    ! (Note, if using water tracers, the water tracer arrays passed into
    !  convec2 must have size 'npnts', i.e. the 2nd argument in the call)

    call convec2_6a  (k, npnts, npnts, nlev, ntra, n_wtrac, trlev,             &
                      md_on, md_new_termc, start_lev, timestep,                &
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
                      z_theta(:,k), z_rho(:,k+1), z_theta(:,k+1),              &
                      l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,       &
                      bgmk, bgmkp1, bwater(:,k),                               &
                      bwater(:,k+1), binit, bland,                             &

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
  do i = 1,npnts
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
  end if     ! nconv >0


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
      if (k==start_lev(j)) then
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
  ! 3.4  Cape and CFL scaling - adjust initial mass flux so that cape is
  !      removed by convection over timescale cape_timescale.
  !-----------------------------------------------------------------------

  ! Set up integer nterm which is the total number of points where
  ! convection has terminated.
  ! Index to full array (npnts) with index_nterm

  nterm = 0
  do i = 1,npnts
    if (bterm(i)) then
      nterm = nterm + 1
      index_nterm(nterm) = i
    end if
  end do

  if (l_fcape .and. nterm > 0) then
    !If l_fcape then replace cape with f.det profile weighted cape
    do j = 1,nterm
      i = index_nterm(j)
      cape(i) = fcape(i)/deltaktot(i)
    end do
  end if

  if (nterm >  0) then

    select case (cldbase_opt_md)

      ! Default 4a convection scheme - RH-based CAPE closure
    case (0)

      ! NEC compiler directive
      !CDIR NODEP
      do j = 1,nterm
        i = index_nterm(j)
        if (dcpbydt(i)  >   0.0) then
          cape_ts_new(i) =                                                     &
                  min(max(900.0*(1.0 - relh(i)/dptot(i))/0.1,60.0)             &
                                                      ,cape_timescale)
        end if
      end do

      ! Modified 4a convection scheme - RH-based CAPE closure
      ! timescale limited to timestep
    case (1)

      ! NEC compiler directive
      !CDIR NODEP
      do j = 1,nterm
        i = index_nterm(j)
        if (dcpbydt(i)  >   0.0) then
          ! Only use RH cape if thickness of convective layer is greater
          ! than 150hPa
          if (p_layer_centres(i,start_lev(i))-                                 &
               p_layer_centres(i,k)  >   15000.0) then

            cape_ts_new(i) =                                                   &
                    min(max(cape_timescale*(1.0-relh(i)/dptot(i))/0.4          &
                  ,timestep)                                                   &
                  ,cape_timescale)

          else
            cape_ts_new(i) = cape_timescale
          end if
        end if  ! dcpbydt > 0
      end do

      ! Fixed cape timescale
    case (2)

      ! NEC compiler directive
      !CDIR NODEP
      do j = 1,nterm
        i = index_nterm(j)
        if (dcpbydt(i)  >   0.0) then
          cape_ts_new(i) =  cape_timescale
        end if
      end do

      ! w based cape closure; if w_cape_limit > 1000. reverts to cape_timescale
    case (3)

      ! Initialise array to cape timescale and then alter as required.
      cape_ts_new(:) =  cape_timescale

      if ( w_cape_limit < 1000.0 ) then
        !  This section includes test on w_max

        !CDIR NODEP
        do j = 1, nterm
          i = index_nterm(j)
          if ( dcpbydt(i) > 0.0 ) then
            ! new denominator introduced at vn6.6
            if ( w_max(i) > w_cape_limit ) then
              cape_ts_new(i) =   cape_timescale * w_cape_limit/                &
                      (w_cape_limit+ (w_max(i)-w_cape_limit)*wcape_fac)
            else
              cape_ts_new(i) = cape_timescale
            end if !  w_max(i) > w_cape_limit
          end if  ! dcpbydt > 0
        end do
      end if  ! w_cape_limit

      ! Option 4 - a w based cape closure; Grid-box area scaled CAPE closure
    case (4)

      if ( w_cape_limit < 1000.0 ) then
        !  This section includes test on w_max

        !CDIR NODEP
        do j = 1, nterm
          i = index_nterm(j)
          if ( dcpbydt(i) > 0.0 ) then
            ! new denominator introduced at vn6.6
            if ( w_max(i) > w_cape_limit ) then
              cape_ts_new(i) =   cape_timescale * w_cape_limit/w_max(i)
            else
              cape_ts_new(i) = cape_timescale * cape(i) / cape_min +           &
                               cape_timescale * exp(-cape(i) / cape_min)
            end if !  w_max(i) > w_cape_limit
          end if  ! dcpbydt > 0
        end do
      else
        do j = 1, nterm
          i = index_nterm(j)
          if ( dcpbydt(i) > 0.0 ) then

            cape_ts_new(i) = cape_timescale * cape(i) / cape_min +             &
                             cape_timescale * exp( - cape(i) /cape_min)

          end if  ! dcpbydt > 0
        end do
      end if  ! w_cape_limit

      ! Option 5 - a w based cape closure;
    case (5)

      if ( w_cape_limit < 1000.0 ) then
        !  This section includes test on w_max

        !CDIR NODEP
        do j = 1, nterm
          i = index_nterm(j)
          if ( dcpbydt(i) > 0.0 ) then
            if ( w_max(i) > w_cape_limit ) then
              cape_ts_new(i) =   cape_timescale * w_cape_limit/w_max(i)
            else
              if ( relh(i) / dptot(i) >= 0.75 ) then
                cape_ts_new(i) = cape_timescale *                              &
                                    ( 0.2373 / (relh(i) / dptot(i))**5)
              else
                cape_ts_new(i) = cape_timescale
              end if ! relh(i) / dptot(i) >= 0.75
            end if !  w_max(i) > w_cape_limit
          end if  ! dcpbydt > 0
        end do
      else
        do j = 1, nterm
          i = index_nterm(j)
          if ( dcpbydt(i) > 0.0 ) then
            if ( relh(i) / dptot(i) >= 0.75 ) then
              cape_ts_new(i) = cape_timescale *                                &
                                  ( 0.2373 / (relh(i) / dptot(i))**5)
            else
              cape_ts_new(i) = cape_timescale
            end if ! relh(i) / dptot(i) >= 0.75
          end if  ! dcpbydt > 0
        end do
      end if  ! w_cape_limit

      ! RH and w based CAPE option
      ! Expects a sensible w_cape_limit or will do nothing
    case (6)

      rh_test = 0.60         ! critical RH value
      rh_fac  = 1.0/ (1.0 - rh_test)

      if ( w_cape_limit < 1000.0 ) then
        !  This section includes test on w_max
        !CDIR NODEP
        do j = 1, nterm
          i = index_nterm(j)
          if ( dcpbydt(i) > 0.0 ) then
            ! work out any reduction in cape_timescale due to RH
            ! linearly falls to 1/2 given cape_timescale for RH above rh_test
            if ( relh(i) / dptot(i) >= rh_test ) then
              cape_ts_new(i) = cape_timescale*0.5*                             &
                               (1.0+(1.0-(relh(i) / dptot(i)))*rh_fac)
            else
              cape_ts_new(i) = cape_timescale
            end if
            ! Further reduction if w_max above critical value
            if ( w_max(i) > w_cape_limit ) then
              cape_ts_new(i) = cape_ts_new(i) * w_cape_limit/                  &
                       (w_cape_limit + (w_max(i)-w_cape_limit)*wcape_fac)

            end if
            ! Limit CAPE timescale to convective timestep
            cape_ts_new(i) = max(cape_ts_new(i), timestep)

          end if  ! dcpbydt > 0
        end do
      end if  ! w_cape_limit

      ! CAPE based on large-scale w over convecting column

    case (7)

      do j = 1, nterm
        i = index_nterm(j)
        if ( dcpbydt(i) > 0.0 ) then
          mass_mean(i) = 0.0
          wls_mean(i)  = 0.0
          do kk=start_lev(i),k
            wls_mean(i) = wls_mean(i) +                                        &
                              w(i,kk)*r2rho_th(i,kk)*dr_across_th(i,kk)
            mass_mean(i) = mass_mean(i) + r2rho_th(i,kk)*dr_across_th(i,kk)
          end do
          wls_mean(i) =  wls_mean(i)/mass_mean(i)
          if (wls_mean(i) > 0.0) then
            cape_ts_new(i) = a_cape *(wls_mean(i)**b_cape)
          else    ! set to maximum value
            cape_ts_new(i) = cape_ts_max_use
          end if
          ! Limit CAPE timescale within min/max values:
          cape_ts_new(i) = min(cape_ts_new(i), cape_ts_max_use)
          cape_ts_new(i) = max(cape_ts_new(i), cape_ts_min_use)
        end if  ! dcpbydt > 0
      end do

    end select        ! cldbase_opt_md

    ! Use new cape timescale calculated above

    ! NEC compiler directive
    !CDIR NODEP
    do j = 1,nterm
      i = index_nterm(j)
      if (dcpbydt(i)  >   0.0) then
        flx_init_new(i) = flx_init(i)*cape(i)/(cape_ts_new(i)*dcpbydt(i))

        if (flx_init_new(i)  >   flxmax_init(i)) then
          flx_init_new(i) = flxmax_init(i)
        end if

        ! Scale max_cfl with cape scale

        max_cfl(i) = max_cfl(i) * flx_init_new(i) / flx_init(i)
      else
        flx_init_new(i) = 0.0
      end if  ! dcpbydt > 0

    end do

    ! Time smoothed mass flux for closure
    if (l_conv_prog_flx) then

      decay_amount = timestep/tau_conv_prog_flx

      do j = 1,nterm
        i = index_nterm(j)
        ! Calculate scaling.
        scale_f(i)      = flx_init_new(i) / flx_init(i)
      end do

      ! Take a copy of the smoothed mass flux in case of failed convection
      ! And update the smoothed mass flux
      do kt = 1, nlev
        do j = 1,nterm
          i = index_nterm(j)
          if (kt  >=  start_lev(i)) then
            tmp_conv_prog_flx(i,kt) = conv_prog_flx(i,kt)
            conv_prog_flx(i,kt)     = conv_prog_flx(i,kt) +                    &
                                  decay_amount * scale_f(i) * flx(i,kt)
          end if
        end do
      end do

      ! And update flx_init_new to its time-smoothed value
      do j = 1,nterm
        i = index_nterm(j)
        max_cfl(i)      = max_cfl(i) * conv_prog_flx(i,start_lev(i)) /         &
                          flx_init_new(i)
        flx_init_new(i) = conv_prog_flx(i,start_lev(i))

        if (flx_init_new(i) > flxmax_init(i)) then
          max_cfl(i)      = max_cfl(i) * flxmax_init(i)/flx_init_new(i)
          flx_init_new(i) = flxmax_init(i)
        end if
      end do

    end if


    ! Work out scaled mass flux needed to keep cfl ratio below limit.
    ! L_CAPE assumed to be true.

    do j = 1,nterm
      i = index_nterm(j)
      max_cfl(i) = max_cfl(i) * timestep

      if (max_cfl(i)  >   cfl_limit) then
        flx_init_new(i) = flx_init_new(i) * cfl_limit / max_cfl(i)
        cfl_limited(i) = 1.0
      else
        flx_init_new(i) = flx_init_new(i)
      end if

      if (flx_init_new(i)  >   flxmax_init(i)) then
        flx_init_new(i) = flxmax_init(i)
      end if
      max_cfl(i) = 0.0

      ! Scale cloud fraction

      if (flx_init_new(i)  >   0.0) then
        scale_f(i) = flx_init_new(i) / flx_init(i)
        cca_2d(i) = cca_2d(i) + 0.06 * log(scale_f(i))

        ! Check scaled cloud fraction not smaller than minimum value
        ! (2.0E-5) or greater than unity.

        cca_2d(i) = max(2.0e-5,cca_2d(i))
        if (cca_2d(i)  >   1.0) then
          cca_2d(i) = 1.0
        end if
      end if

    end do  ! nterm

    !-----------------------------------------------------------------------
    ! Check for false/true convection and reset variables appropriately
    !-----------------------------------------------------------------------

    ! First, save SCM diagnostic of whether convection failed, and if so, why:
    if ( l_scm_convss_dg ) then
      do j = 1,nterm
        i = index_nterm(j)

        ! Set mid-convection status; note that we don't want to set the
        ! status to failed if mid-level convection actually succeeded
        ! in a lower layer.  So test pre-existing status before resetting.

        ! Convection failed due to the ascent being too shallow or cloud-free
        if ( icct(i)-iccb(i) <= 3 .or. ( .not. blatent(i) ) ) then
          if ( scm_convss_dg(i) % status_mid < 1 )                             &
               scm_convss_dg(i) % status_mid = 1

          ! Convection failed because the closure set the mass-flux to zero
        else if ( flx_init_new(i)  <=  minflx ) then
          if ( scm_convss_dg(i) % status_mid < 2 )                             &
               scm_convss_dg(i) % status_mid = 2

          ! Real mid-level convection occurred!
        else
          if ( scm_convss_dg(i) % status_mid < 3 )                             &
               scm_convss_dg(i) % status_mid = 3

        end if

      end do
    end if

    ! Do the actual test for failed convection...
    do j = 1,nterm
      i = index_nterm(j)
      if ( (flx_init_new(i)  <=  minflx)                                       &
                  .or. ( (icct(i)-iccb(i)) <= 3 )                              &
                  .or. ( .not. blatent(i) ) ) then
        ! False convection diagnosed if:
        ! - the new initial mass flux is less than zero
        ! - or the convecting layer it too thin
        ! - or the parcel has never released latent heat
        ! 3d variables are reset below by setting scale_f to zero.
        flx_init(i)     = 0.0
        flx_init_new(i) = 0.0
        scale_f(i)      = 0.0
        cca_2d(i)       = 0.0
        iccb(i)         = 0
        icct(i)         = 0
        tcw(i)          = 0.0
        cclwp(i)        = 0.0
        lcca(i)         = 0.0
        lctop(i)        = 0
        lcbase(i)       = 0
        kterm_mid(i)    = 0
        kterm(i)        = 0
        det_lev(i)      = 0
      else
        ! True convection
        l_mid_all(i)    = .true.
        kterm_mid(i)    = k
        kterm(i)        = k
        det_lev(i)      = k+1
        ! Set lowest mid-level scaling factor if not set
        if (scale_low(i) == 1.0) then
          scale_low(i) = scale_f(i)
        end if
      end if
    end do  ! nterm

    !-----------------------------------------------------------------------
    ! Apply closure and cfl scaling
    !-----------------------------------------------------------------------
    do kt = 1, k+1
      do j = 1,nterm
        i = index_nterm(j)
        if (kt  >=  start_lev(i)) then
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

          flx(i,kt)    = flx(i,kt)    * scale_f(i)
          precip(i,kt) = precip(i,kt) * scale_f(i)

        end if ! kt>start_lev and flx_init_new >0
      end do  ! j loop
    end do  ! kt loop

    if (l_tracer) then
      do ktra = 1,ntra
        do kt = 2, k+1
          do j = 1,nterm
            i = index_nterm(j)
            if (kt  >=  start_lev(i)) then
              dtrabydt(i,kt,ktra) =dtrabydt(i,kt,ktra)*scale_f(i)
            end if ! kt>start_lev and flx_init_new >0
          end do  ! j loop
        end do  ! kt loop
      end do
    end if

    if (l_wtrac_conv) then
      do i_wt = 1,n_wtrac
        do kt = 1, k+1
          do j = 1,nterm
            i = index_nterm(j)
            if (kt  >=  start_lev(i)) then
              wtrac_e(i_wt)%dqbydt(i,kt)   =                                   &
                      wtrac_e(i_wt)%dqbydt(i,kt)   * scale_f(i)
              wtrac_e(i_wt)%dqclbydt(i,kt) =                                   &
                      wtrac_e(i_wt)%dqclbydt(i,kt) * scale_f(i)
              wtrac_e(i_wt)%dqcfbydt(i,kt) =                                   &
                      wtrac_e(i_wt)%dqcfbydt(i,kt) * scale_f(i)
              wtrac_p(i_wt)%precip(i,kt)    =                                  &
                      wtrac_p(i_wt)%precip(i,kt)   * scale_f(i)

            end if ! kt>start_lev
          end do  ! j loop
        end do  ! kt loop
      end do
    end if

    ! Correct time smoothed mass flux in case of failed convection
    if (l_conv_prog_flx) then
      ! Update the smoothed mass flux
      do kt = 1, nlev
        do j = 1,nterm
          i = index_nterm(j)
          ! Test for failed convection
          if (kterm(i) == 0) then
            if (kt  >=  start_lev(i)) then
              conv_prog_flx(i,kt) = tmp_conv_prog_flx(i,kt)
            end if
          end if
        end do
      end do
    end if

    !-----------------------------------------------------------------------
    ! 3.5  Original downdraught calculation - on all points where convection
    !      is terminating.
    !      UM Documentation Paper 27
    !-----------------------------------------------------------------------

    if (l_eman_dd) then

      ! save termination level

      do j=1,nterm
        i = index_nterm(j)
        kterm_mid(i) = kterm(i)
      end do

    else        ! original down draughts

      npossdd = 0
      nnodd = 0
      do j = 1,nterm
        i = index_nterm(j)
        tempnum = 0.0
        if (iccb(i)  >   0) then
          deltap_cld = p_layer_centres(i,iccb(i))                              &
                             - p_layer_centres(i,k)
          do kt = iccb(i), k+1
            tempnum = tempnum + precip(i,kt)
          end do
        else
          deltap_cld = 0.0
        end if

        ! Downdraughts possible if pressure thickness of convective
        ! cloud (deltap_cld) is greater than 15000Pa, the point is saturated
        ! and the precip. in the layer is greater than a threashold
        ! value (1E-12).

        if (deltap_cld  >   15000.0 .and. bgmk(i) .and.                        &
                                   tempnum  >   1e-12) then
          npossdd = npossdd + 1
          index_possdd(npossdd) = i
        else
          nnodd = nnodd + 1
          index_nodd(nnodd) = i
        end if
      end do  ! nterm loop

      ! If some downdraughts are possible (npossdd > 0) then call
      ! downdraught code

      if (npossdd  >   0) then

        call dd_call_6a(npnts, npossdd, k, nlev, trlev, ntra, n_wtrac          &
          ,iccb, icct, index_possdd                                            &
          ,l_tracer                                                            &
          ,bwater(1,2)                                                         &
          ,exner_layer_centres, exner_layer_boundaries                         &
          ,p_layer_centres, p_layer_boundaries, pstar                          &
          ,recip_pstar, timestep , cca_2d                                      &
          ,thp(1,1), qp(1,1), th(1,1), q(1,1),qse                              &
          ,trap, tracer                                                        &
          ,flx(1,1), precip(1,1)                                               &
          ,dthbydt(1,1), dqbydt(1,1), dtrabydt                                 &
          ,rain, snow ,rain_3d, snow_3d                                        &
          ,wtrac_p, wtrac_e                                                    &
          ,dwn_flux, entrain_dwn, detrain_dwn, dt_dd, dq_dd)

      end if

      ! Surface precipitation calculation for points where downdraught not
      ! possible.

      if (nnodd  >   0) then
        call evap_bcb_nodd_6a (npnts,nnodd,n_wtrac,k,iccb,index_nodd,          &
                     bwater(1,2),timestep,                                     &
                     p_layer_centres,p_layer_boundaries,                       &
                     exner_layer_centres,                                      &
                     th, q, qse, cca_2d,                                       &
                     precip, dthbydt, dqbydt,                                  &
                     rain, snow, rain_3d, snow_3d,                             &
                     dt_dd, dq_dd, wtrac_p, wtrac_e)

      end if

    end if  ! Emanuel/new DD scheme test


    ! If convection has terminated write cape to diagnostic output
    ! variable (cape_out). Note if more than one lot of mid-level
    ! convection or if mid above deep or shallow then stores highest
    ! convection results only.
    ! Zero integrals as convection terminated.
    ! NB convection may start again at the
    ! same locations at higher levels in the atmosphere.
    ! NEC compiler directive
    !CDIR NODEP
    do j = 1,nterm
      i=index_nterm(j)
      cape_out(i)   = cape(i)
      dcpbydt(i)    = 0.0
      cape(i)       = 0.0
      bconv(i)      = .false.
      bterm(i)      = .false.     ! reset bterm to false
      start_lev(i)  = 0
      relh(i)       = 0.0      ! needed for RH CAPE
      dptot(i)      = 0.0      !
      bgmk(i)       = .false.
      blatent(i)    = .false.
    end do

  end if        ! nterm >0

  !-----------------------------------------------------------------------
  ! Write out entrainment, detrainment and half-level mass flux diagnostics.
  ! They will be scaled by the full level mass flux outside
  ! of the level loop
  !-----------------------------------------------------------------------
  ! Calculate fractional entrainment rate for level k.
  if (flg_entr_up .or. mid_cmt_opt == 1) then
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

  !-----------------------------------------------------------------------
  ! 3.6  End of main loop over levels
  !-----------------------------------------------------------------------
end do

! Deallocate water tracer fields
if (l_wtrac_conv) then
  deallocate(qpixs_v_wtrac)
end if
call wtrac_dealloc_conv_p(n_wtrac, wtrac_p)

! If used, copy convection profile diagnostics
if (l_scm_convss_dg) then

  ! Copy profile diagnostics
  do k = 1, nlev
    do i = 1, npnts
      ! Take care not to overwrite parcel properties of any shallow
      ! or deep convection beneath the mid-level convection.
      if ( k >= midtrig(i) .and. scm_convss_dg(i) % par_thetav(k)>0.0 ) then
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
! Write out updraught massflux diagnostics and scale the
! entrainment and detrainment diagnostics by the mass flux.
!-----------------------------------------------------------------------
if (flg_up_flx .or. flg_mf_midlev) then
  do k = 1,nlev
    do i = 1,npnts
      up_flux(i,k) = flx(i,k)
    end do
  end do
end if

if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,npnts
      up_flux_half(i,k) = up_flux_half(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_entr_up) then
  do k = 1,nlev
    do i = 1,npnts
      entrain_up(i,k) = entrain_up(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,npnts
      detrain_up(i,k) = detrain_up(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_area_ud) then
  do k = 1,nlev
    do i = 1,npnts
      if (flx(i,k) > 0.0) then
        ! expect to get a valid wup
        ! updraught core area from wup and flx

        if (w2p(i,k) > 0.0) then
          area_ud(i,k) =  flx(i,k)/(g*sqrt(w2p(i,k))* rho_theta(i,k))
        else
          ! The calculation of w2p is not producing a sensible value
          ! Assume wcb =1m/s and calculate a value for cloud base
          if (k == iccb(i)+1) then
            area_ud(i,k) =  flx(i,k)/(g* rho_theta(i,k))
          end if
        end if
      end if
    end do
  end do
end if

! Has mid-level convection gone to the top of the model ?

do i=1,npnts
  if (kterm_mid(i) == nlev-1 .and. l_mid_all(i) ) then
    if (printstatus >= prstatus_normal) then
      write(umMessage,'(a41,i8)') 'WARNING: mid_conv has reached nlev-1 at ',i
      call umPrint(umMessage,src='mid_conv_mod-6a')
      ! Extra info which may prove useful to debugging problem
      ! The move to umPrint causes problem for these write statements as
      ! nlev varies and no format is appropriate to print all of each
      ! profile. Redoing to loop over level printing several fields at once
      write(umMessage,'(a7,g16.8)')    ' rain  ',rain(i)
      call umPrint(umMessage,src='mid_conv_mod-6a')
      write(umMessage,'(a7,g16.8)')    ' snow  ',snow(i)
      call umPrint(umMessage,src='mid_conv_mod-6a')
      write(umMessage,'(A6,11A26)') 'k','Theta','q','qcl','qcf',               &
         'up mass flux','dq','dqcl','dqcf','dth','du','dv'
      call umPrint(umMessage,src='mid_conv_mod-6a')
      do k=1,nlev
        write(umMessage,'(I6,11G26.18)') k,th(i,k),q(i,k),qcl(i,k),qcf(i,k),   &
            flx(i,k),dqbydt(i,k),dqclbydt(i,k),dqcfbydt(i,k),dthbydt(i,k),     &
            dubydt(i,k),dvbydt(i,k)
        call umPrint(umMessage,src='mid_conv_mod-6a')
      end do
    end if

    error_point = i
    if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    return
  end if
end do


! (Was copied out of scaling if test to ensure these limits at all
! times, not just when cca_2d is scaled)

! Check scaled cloud fraction not smaller than minimum value
! (2.0E-5) or greater than unity.

do i=1, npnts
  if (cca_2d(i) > 0.0) then
    cca_2d(i) = max(2.0e-5, cca_2d(i))
    if (cca_2d(i) > 1.0) then
      cca_2d(i) = 1.0
    end if
  end if
end do

!-----------------------------------------------------------------------
! How many columns have mid-level convection
!-----------------------------------------------------------------------

nmid = 0
do i=1,npnts
  if (l_mid_all(i)) then
    nmid = nmid +1
    index1(nmid) = i
  end if
end do



!-----------------------------------------------------------------------
! Emanuel down draughts - one call for all points with mid-level conv
!-----------------------------------------------------------------------

if (l_eman_dd) then

  ! how many points have mid level convection?
  ! What is the maximum termination level ?

  kterm_max = 2

  do i=1,npnts
    if (l_mid_all(i)) then

      if (kterm_mid(i) >  kterm_max) then
        kterm_max = kterm_mid(i)
      end if

    end if
  end do

  ! Call routine to compress to just those points with convection
  ! and call Emanuel downdraughts then expand back to full grid.

  call eman_cex(npnts, nmid, kterm_max, nlev, trlev ,ntra                      &
,                  kterm_mid, l_mid_all, l_tracer                              &
,                  exner_layer_centres,exner_layer_boundaries                  &
,                  p_layer_centres, p_layer_boundaries                         &
,                  timestep, scale_low, th, q, qse, tracer                     &
,                  precip, dthbydt, dqbydt, dtrabydt                           &
,                  rain, snow, dwn_flux, dt_dd, dq_dd)

end if

!-----------------------------------------------------------------------
! 4.0  Diffusive Convective momentum tranport option
!-----------------------------------------------------------------------

if (l_mom .and. mid_cmt_opt == 1) then

  ! Maximum number of mid-level layers, required for holding cloud bases
  ! and top of each layer. Note each layer must be at least 2 levels and have at
  ! least one level between. The bottom model level is not used and don't expect
  ! convection in any stratospheric levels so trying estimating by dividing
  ! by 6. This should give a number bigger than required but not excessive.

  nmax_layers = nlev/6

  ! Only call if there are points with mid_level convection.
  if (nmid > 0) then

    call mid_conv_dif_cmt(npnts, nmid, nlev, nmax_layers,                      &
                       index1,                                                 &
                       timestep,                                               &
                       u, v, r_theta, r_rho,                                   &
                       z_theta, z_rho, rho,                                    &
                       flx, entrain_up,                                        &
                       dubydt, dvbydt, uw_mid, vw_mid )

  end if  ! test on mid points
end if  ! test on mid_cmt_opt 1


!-----------------------------------------------------------------------
! 4.1  Add the dissipative heating from the CMT to the theta increment
!-----------------------------------------------------------------------
if (l_mom .and. l_cmt_heating) then
  call cmt_heating(npnts, nlev,                                                &
                   z_theta, z_rho, exner_layer_centres,                        &
                   u, v, dubydt, dvbydt,                                       &
                   ! Out
                   dthbydt)
end if


! ---------------------------------------------------------------------
! 5.0  Energy (and optionally water) correction calculation
! ---------------------------------------------------------------------



if (nmid >  0) then

  if (l_cv_conserve_check) then
    call cor_engy_6a(npnts, nmid, nlev, n_wtrac, index1, timestep,             &
                     p_layer_boundaries, exner_layer_centres,                  &
                     exner_rho, r_theta, r_rho, rho_dry,                       &
                     r2rho, rho_theta, rho_dry_theta,                          &
                     dr_across_rh,                                             &
                     dubydt, dvbydt, dqclbydt, dqcfbydt,                       &
                     rain, snow, q, qcl, qcf, u, v,                            &
                     !In/Out
                     dthbydt, dqbydt, wtrac_e)
  end if

  !-----------------------------------------------------------------------
  ! 6.0  Correct negative/very small humidities
  !-----------------------------------------------------------------------

  call correct_small_q_conv(nmid, npnts, nlev, n_wtrac, index1, timestep,      &
                            qmin, r2rho_th, dr_across_th,  q, 'water',         &
                            'Mid-level', dqbydt, wtrac_e=wtrac_e)

  ! Then use same method to correct any remaining negative/v small values
  ! of water tracer vapour

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
      call correct_small_q_conv(nmid, npnts, nlev, n_wtrac, index1, timestep,  &
                                wtrac_info(i_wt)%qlimit, r2rho_th,             &
                                dr_across_th,                                  &
                                wtrac_e(i_wt)%q, 'wtrac', 'Mid-level',         &
                                wtrac_e(i_wt)%dqbydt)
    end do
  end if

end if    ! test on nmid

!-----------------------------------------------------------------------------
! 7.1 Calculate CCA fow MID-CONVECTION levels only
!-----------------------------------------------------------------------------

n_mdcld      = 0    ! Initialise no of points with mid-cld
n_mdcld_mult = 0    ! Initialise no of points with multiple mid-cld

do i=1, npnts
  if ((l_mid_all(i)) .and. (lcbase(i) /= 0)) then

    ! Mid-level cloud present
    ! Create index of gridpoints with single/multiple mid-level
    ! cloud

    if (lcbase(i) /= iccb(i)) then

      n_mdcld_mult       = n_mdcld_mult + 1
      dum2(n_mdcld_mult) = i

    else if (lcbase(i) == iccb(i)) then

      n_mdcld       = n_mdcld + 1
      dum1(n_mdcld) = i
    end if
  end if
end do

!---------------------------------------------------------------
! Calculate CCA_2D of Mid Cloud
!---------------------------------------------------------------
select case (cca2d_md_opt)
case (srf_precip)
  do i=1, npnts
    if ((l_mid_all(i)) .and. (lcbase(i) /= 0)) then
      !-----------------------------------------------------------
      ! CCA_2D based on surface precipitation rate from mid
      ! cloud
      !-----------------------------------------------------------

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
        ! shallow or deep have occured.  This is under a switch in
        ! the 4a scheme.
        !
        ! NOTE: Downdraughts & Evaporation still being fed cca_2d
        !       derived from TCW, This issue may require further
        !       investigation.
        lcca(i) = cca_2d(i)

      end if
    end if    ! mid-level cloud present
  end do      ! npnts

case (avg_mass_flux)
  !Initialise the average flx/rho
  do i=1, npnts
    flxbyrho_int(i)   = 0.0
    tot_mass_mean(i)  = 0.0
  end do

  !calculate the average flx/rho and convecting layer mass
  do k = 1,nlev
    do i= 1, npnts
      if (flx(i,k) > 0.0) then
        flxbyrho_int(i)  = flxbyrho_int(i)  + flx(i,k)/g     * dr_across_th(i,k)
        tot_mass_mean(i) = tot_mass_mean(i) + rho_theta(i,k) * dr_across_th(i,k)
      end if
    end do
  end do

  do i=1, npnts
    if (iccb(i) > 0) then ! mid convection was successful
      !-----------------------------------------------------------
      ! CCA_2D is derived from the mass averaged mass flux
      !-----------------------------------------------------------

      cca_2d(i) = flxbyrho_int(i) / (w_cca*tot_mass_mean(i))

      !Limit to avoid unphysical values
      cca_2d(i) = max(min(cca_2d(i), 1.0),0.0)

      ! Grab lowest cca value before any tuning occurs
      ! This will overwrite lcca in ni_conv_ctl only if neither
      ! shallow or deep have occured.
      lcca(i) = cca_2d(i)

    end if      ! iccb
  end do      ! i (npnts

case (total_condensed_water)
    ! cca_2d left unchanged from code, which is based on
    ! TCW (Total Condensed Water) (This is a rate)

end select      ! cca2d_md_opt

! l_dcpl_cld4pc2 is set to true for 5a scheme, so
! CCRad tuning knobs applied in glue_conv

!-------------------------------------------------------------------
! 7.2 Apply CCA_3D
!-------------------------------------------------------------------
! A value of cca_2d has been calculated, though we need to know
! which levels to apply it to. So we identify which layers require
! CCA based on layers (above shallow/deep, if they have occurred)
! with non-zero ccw.
!-------------------------------------------------------------------

!===================================================================
! Single mid-level cloud bank gridpoints
!===================================================================

if (n_mdcld > 0) then

  ! Resize compressed arrays for single Mid-level events

  allocate (mdcldi(n_mdcld))

  ! Now we know the number and indexes of gridpoints with
  ! single/multiple mid-level cloud, we can now assign the
  ! compressed mid-level arrays to pass to the anvil scheme.

  allocate (iccb_md_c       (n_mdcld       ))
  allocate (icct_md_c       (n_mdcld       ))
  allocate (freeze_lev_md_c (n_mdcld       ))
  allocate (cca_2d_md_c     (n_mdcld       ))
  allocate (cca_md_c        (n_mdcld,  n_cca_lev))
  allocate (ccw_md_c        (n_mdcld,  nlev))
  allocate (z_theta_md_c    (n_mdcld,  nlev))
  allocate (z_rho_md_c      (n_mdcld,  nlev))
  allocate (p_lyr_bnds_md_c (n_mdcld,0:nlev))

  do i=1, n_mdcld
    mdcldi          (i) = dum1(i)
    iccb_md_c       (i) = iccb       (mdcldi(i))
    icct_md_c       (i) = icct       (mdcldi(i))
    cca_2d_md_c     (i) = cca_2d     (mdcldi(i))
    freeze_lev_md_c (i) = freeze_lev (mdcldi(i))
  end do

  do k=1, nlev
    do i=1, n_mdcld
      ccw_md_c        (i,k) = ccw     (mdcldi(i),k)
      z_theta_md_c    (i,k) = z_theta (mdcldi(i),k)
      z_rho_md_c      (i,k) = z_rho   (mdcldi(i),k)
      p_lyr_bnds_md_c (i,k) =                                                  &
                    p_layer_boundaries(mdcldi(i),k)
    end do
  end do
  do i=1, n_mdcld
    p_lyr_bnds_md_c(i,0) = p_layer_boundaries(mdcldi(i),0)
  end do

  do k=1, n_cca_lev
    do i=1, n_mdcld
      cca_md_c(i,k) = 0.0
    end do
  end do

  !-----------------------------------------------------------------
  ! End Resizing compressed arrays for single Mid-level events.
  ! There is only one mid-level cloud bank if a mid-event has
  ! occurred, so does not require checking to see if there are
  ! multiple mid-level cloud banks.
  !-----------------------------------------------------------------

  if (l_anvil) then
    call calc_3d_cca                                                           &
      ( n_mdcld, n_mdcld, nlev, n_cca_lev, nbl, iccb_md_c, icct_md_c           &
      , p_lyr_bnds_md_c, freeze_lev_md_c, cca_2d_md_c, cca_md_c                &
      , z_theta_md_c, z_rho_md_c )
  else

    !---------------------------------------------------------------
    ! Do not use anvil scheme and apply cca_2d_md_c to all cloud
    ! levels between mid-level base and top
    !---------------------------------------------------------------
    do k=1, n_cca_lev
      do i=1, n_mdcld
        if ( k >= iccb_md_c(i) .and.                                           &
             k <= icct_md_c(i) ) then
          cca_md_c(i,k) = cca_2d_md_c(i)
        end if
      end do      ! i (n_mdcld)
    end do      ! k (n_cca_lev)
  end if     ! l_anvil

  !-----------------------------------------------------------------
  ! Merge cca_md/ccw_md to full cca array and scale ccw_md_c by
  ! ccw_md_knob
  !-----------------------------------------------------------------

  do k=1, n_cca_lev
    do i=1, n_mdcld
      cca(mdcldi(i),k) = cca(mdcldi(i),k) + cca_md_c(i,k)
    end do      ! i (n_mdcld)
  end do      ! k (n_cca_lev)

  ! l_dcpl_cld4pc2 is set to true for 5a scheme, so
  ! CCRad tuning knobs applied in glue_conv

  ! Deallocate compressed single mid-level cloud arrays

  deallocate (mdcldi)
  deallocate (iccb_md_c)
  deallocate (icct_md_c)
  deallocate (freeze_lev_md_c)
  deallocate (cca_2d_md_c)
  deallocate (cca_md_c)
  deallocate (ccw_md_c)
  deallocate (z_theta_md_c)
  deallocate (z_rho_md_c)
  deallocate (p_lyr_bnds_md_c)

end if ! n_mdcld > 0

!===================================================================
! Multiple mid-level cloud bank gridpoints
!===================================================================

if (n_mdcld_mult > 0) then
  ! Resize index array for gridpoints with multiple Mid-level cloud
  ! banks and apply indices

  allocate(mdcldi_mult(n_mdcld_mult))

  !-----------------------------------------------------------------
  ! Allocate arrays with multiple-mid level cloud gridpoints
  !-----------------------------------------------------------------
  allocate (iccb_md_mult_c       (n_mdcld_mult       ))
  allocate (icct_md_mult_c       (n_mdcld_mult       ))
  allocate (freeze_lev_md_mult_c (n_mdcld_mult       ))
  allocate (cca_2d_md_mult_c     (n_mdcld_mult       ))

  ! Lets allocate them with levels first since it will then be continuous in
  ! memory later on when we call calc_3d_cca.
  allocate (cca_md_mult_c        ( n_cca_lev, n_mdcld_mult))
  allocate (ccw_md_mult_c        (  nlev, n_mdcld_mult))
  allocate (z_theta_md_mult_c    (  nlev, n_mdcld_mult))
  allocate (z_rho_md_mult_c      (  nlev, n_mdcld_mult))
  allocate (p_lyr_bnds_md_mult_c (0:nlev, n_mdcld_mult))

  do i=1, n_mdcld_mult
    mdcldi_mult          (i) = dum2(i)
    iccb_md_mult_c       (i) = 0
    icct_md_mult_c       (i) = 0
    freeze_lev_md_mult_c (i) = freeze_lev (mdcldi_mult(i))
    cca_2d_md_mult_c     (i) = cca_2d     (mdcldi_mult(i))

    do k=1, nlev
      ccw_md_mult_c       (k,i) = ccw     (mdcldi_mult(i),k)
      z_theta_md_mult_c   (k,i) = z_theta (mdcldi_mult(i),k)
      z_rho_md_mult_c     (k,i) = z_rho   (mdcldi_mult(i),k)
      p_lyr_bnds_md_mult_c(k,i) =                                              &
                        p_layer_boundaries(mdcldi_mult(i),k)
    end do
    p_lyr_bnds_md_mult_c(0,i) = p_layer_boundaries(mdcldi_mult(i),0)

    do k=1, n_cca_lev
      cca_md_mult_c(k,i) = 0.0
    end do
  end do     ! i (n_mdcld_mult)
  ! Breaking loop to help with optimisation
  do i=1, n_mdcld_mult
    ! For multiple mid-level cloud banks, increment up model levels
    ! through compressed ccw_md_c array. Enter CAL3DCCA on locating
    ! cloud/bases of multiple clouds.

    do k=2, nlev

      !-------------------------------------------------------------
      ! Check for cloud base
      !-------------------------------------------------------------
      if ((ccw_md_mult_c(k,i) >  0.0) .and.                                    &
          (iccb_md_mult_c(i)  == 0)   .and.                                    &
          (icct_md_mult_c(i)  == 0)) then
        iccb_md_mult_c(i) = k
      end if

      !-------------------------------------------------------------
      ! Check for cloud top
      !-------------------------------------------------------------
      if ((ccw_md_mult_c(k,i)   <= 0.0) .and.                                  &
          (ccw_md_mult_c(k-1,i) >  0.0) .and.                                  &
          (iccb_md_mult_c(i)    /= 0)   .and.                                  &
          (icct_md_mult_c(i)    == 0)) then
        icct_md_mult_c(i) = k-1
      end if

      !-------------------------------------------------------------
      ! Check for for anvil if both a cloud base and top found
      !-------------------------------------------------------------
      if (iccb_md_mult_c(i) /= 0 .and.                                         &
          icct_md_mult_c(i) /= 0) then

        ! Apply CCA to vertically continuous cloud
        if (l_anvil) then
          ! Since we are performing this on each individual point we can
          ! pass down the array with level information in the first rank
          ! instead of the second rank.  This reduces the amount of temporary
          ! data which is required otherwise it has to step through each
          ! array by number of points to get the all level information for it.
          call calc_3d_cca                                                     &
            ( 1, 1, nlev, n_cca_lev, nbl, iccb_md_mult_c(i)                    &
            , icct_md_mult_c(i), p_lyr_bnds_md_mult_c(:,i)                     &
            , freeze_lev_md_mult_c(i), cca_2d_md_mult_c(i)                     &
            , cca_md_mult_c(:,i), z_theta_md_mult_c(:,i)                       &
            , z_rho_md_mult_c(:,i) )
        else

          ! Copy cca_2d_md_c to all levels from cloud base to cloud
          ! top of located cloud

          do j=iccb_md_mult_c(i), icct_md_mult_c(i)
            cca_md_mult_c(j,i) = cca_2d_md_mult_c(i)
          end do      ! j
        end if      ! l_anvil

        iccb_md_mult_c(i) = 0
        icct_md_mult_c(i) = 0

      end if      ! Test on iccb_md_c(i) and icct_md_c(i)
    end do      ! k (nlev)
  end do     ! i (n_mdcld_mult)
  ! Breaking loop to help with optimisation
  do i=1, n_mdcld_mult

    !-----------------------------------------------------------------
    ! Merge cca_md/ccw_md to cca full array and scale ccw_md_c by
    ! ccw_md_knob
    !-----------------------------------------------------------------

    do k=1, n_cca_lev
      cca(mdcldi_mult(i),k) = cca(mdcldi_mult(i),k)                            &
                            + cca_md_mult_c(k,i)
    end do      ! k (n_cca_lev)
  end do      ! i (n_mdcld_mult)

  ! l_dcpl_cld4pc2 is set to true for 5a scheme, so
  ! CCRad tuning knobs applied in glue_conv

  !-----------------------------------------------------------------
  ! Deallocate arrays
  !-----------------------------------------------------------------
  deallocate (mdcldi_mult)
  deallocate (iccb_md_mult_c)
  deallocate (icct_md_mult_c)
  deallocate (freeze_lev_md_mult_c)
  deallocate (cca_2d_md_mult_c)
  deallocate (cca_md_mult_c)
  deallocate (ccw_md_mult_c)
  deallocate (z_theta_md_mult_c)
  deallocate (z_rho_md_mult_c)
  deallocate (p_lyr_bnds_md_mult_c)

end if      ! n_mdcld_mult


!-----------------------------------------------------------------------
! Final SCM convection sub-step diagnostics
!-----------------------------------------------------------------------
if ( l_scm_convss_dg ) then

  ! 3-D diagnostics
  do k = 1, nlev
    do i = 1, npnts
      if ( k>=midtrig(i) ) then
        ! Save the final updraft mass-flux
        scm_convss_dg(i) % up_flx(k) = flx(i,k)
      end if
    end do
  end do

  ! 2-D diagnostics
  do i = 1, npnts
    scm_convss_dg(i) % precip_mid = rain(i) + snow(i)
  end do

end if


!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mid_conv_6a
end module mid_conv_6a_mod
