! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Shallow convection scheme

module shallow_conv_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Shallow convection scheme
!   Works only on points diagnosed as shallow in subroutine CONV_DIAG
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


character(len=*), parameter, private :: ModuleName='SHALLOW_CONV_6A_MOD'

contains

subroutine shallow_conv_6a(nbl,nlev,ntra,n_wtrac,n_cca_lev,n_sh,trlev,         &
                       bland,delthvu,                                          &
                       exner_rho,                                              &
                       exner_layer_centres,                                    &
                       exner_layer_boundaries,                                 &
                       l_q_interact,                                           &
                       l_tracer, ntml, ntpar,                                  &
                       pstar,p_layer_centres,                                  &
                       p_layer_boundaries,                                     &
                       z_theta, z_rho,                                         &
                       r_theta, r_rho,                                         &
                       rho_theta,                                              &
                       rho_dry_theta, rho_dry,                                 &
                       r2rho_th, r2rho,                                        &
                       dr_across_th, dr_across_rh,                             &
                       conv_prog_precip,                                       &
                       q,th,timestep,u,v,uw0,vw0,                              &
                       wstar,wthvs,entrain_coef,delta_smag,                    &
                       zlcl_uv,tnuc_nlcl,ztop_uv,freeze_lev,recip_pstar,qse,   &
                       l_scm_convss_dg,                                        &
                       ccp_strength,                                           &

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
                       rain_3d, snow_3d,                                       &
                       up_flux, up_flux_half,                                  &
                       dwn_flux,uw_shall,vw_shall,tcw,cca_2d,                  &
                       kterm,                                                  &
                       ind_shallow, dt_dd, dq_dd, area_ud )

use planet_constants_mod, only:                                                &
    cp, r, kappa, pref, repsilon, c_virtual, g

use water_constants_mod, only: lc, lf, tm

use cv_run_mod, only:                                                          &
    l_mom, sh_pert_opt,                                                        &
    cca2d_sh_opt, c_mass_sh,                                                   &
    l_cv_conserve_check, l_eman_dd,                                            &
    l_snow_rain, l_cmt_heating, cldbase_opt_sh, cape_timescale,                &
    l_conv_prog_flx, tau_conv_prog_flx

use cv_param_mod, only:                                                        &
    total_condensed_water, grant_lock, grant_lock_over,                        &
    thpixs_shallow, qpixs_shallow, c_mass,                                     &
    max_sh_thpert, min_sh_thpert, max_sh_qpert_fac, beta_cu,                   &
    sh_wstar_closure, sh_grey_closure,                                         &
    term_undil

use cv_dependent_switch_mod, only:                                             &
    sh_on, sh_new_termc

use cv_stash_flg_mod, only:                                                    &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                                  &
    flg_entr_up, flg_detr_up, flg_entr_dwn, flg_detr_dwn,                      &
    flg_uw_shall, flg_vw_shall, flg_mf_shall,                                  &
    flg_area_ud

use scm_convss_dg_mod, only: scm_convss_dg_type

use bl_option_mod, only:                                                       &
    kprof_cu, off, max_cu_depth, klcl_entr, rlinfac, linear0

use yomhook,    only: lhook, dr_hook
use parkind1,   only: jprb, jpim

! Subroutines
use lift_par_6a_mod,       only: lift_par_6a
use lift_undil_par_6a_mod, only: lift_undil_par_6a
use convec2_6a_mod,        only: convec2_6a
use water_loading_mod,     only: water_loading
use cor_engy_6a_mod,       only: cor_engy_6a
use mix_ipert_6a_mod,      only: mix_ipert_6a
use cmt_heating_mod,       only: cmt_heating
use layer_cn_6a_mod,       only: layer_cn_6a, shallow
use eman_dd_rev_mod,       only: eman_dd_rev

use dd_all_call_6a_mod,       only: dd_all_call_6a
use evap_bcb_nodd_all_6a_mod, only: evap_bcb_nodd_all_6a
use flag_wet_mod,             only: flag_wet
use shallow_base_stress_mod,  only: shallow_base_stress
use shallow_cmt_incr_mod,     only: shallow_cmt_incr
use shallow_grad_stress_mod,  only: shallow_grad_stress

use qsat_mod, only: qsat

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
 ,n_sh                 & ! No. of shallow convection points
 ,trlev                  ! No. of model levels on which tracers are included

logical, intent(in) :: bland(n_sh) ! Land/sea mask

real(kind=real_umphys), intent(in)    :: delthvu(n_sh)
                                     !Integral of undilute parcel
                                     ! buoyancy over convective cloud
                                     ! layer (Kelvin m)

real(kind=real_umphys), intent(in)    :: exner_rho(n_sh,nlev)
                                             ! Exner on rho levels

real(kind=real_umphys), intent(in)    ::                                       &
  exner_layer_centres(n_sh,0:nlev)   & ! Exner
 ,exner_layer_boundaries(n_sh,0:nlev)  ! Exner at half level above
                                       ! exner_layer_centres

logical, intent(in) ::                                                         &
  l_q_interact         & ! Switch allows overwriting parcel variables when
                         ! calculating condensate incr.
 ,l_tracer               ! Switch for inclusion of tracers

integer, intent(in) ::                                                         &
  ntml(n_sh)           & ! Top level of surface mixed layer defined relative to
                         ! theta,q grid
 ,ntpar(n_sh)            ! Top level of initial parcel ascent in BL scheme
                         ! defined relative to theta,q grid

real(kind=real_umphys), intent(in)    ::                                       &
  pstar(n_sh)                   & ! Surface pressure (Pa)
 ,p_layer_centres(n_sh,0:nlev)  & ! Pressure (Pa)
 ,p_layer_boundaries(n_sh,0:nlev) ! Pressure at half level above
                                  ! p_layer_centres (Pa)

! Note heights passed in but not currently used - will be
! required by new turbulence based scheme therefore been added to
! arguement list ready for future developments.

real(kind=real_umphys), intent(in) :: z_theta(n_sh,nlev)
                                            ! height of theta levels (m)
real(kind=real_umphys), intent(in) :: z_rho(n_sh,nlev)
                                            ! height of rho levels (m)
real(kind=real_umphys), intent(in) :: r_theta(n_sh,0:nlev)
                                            ! radius of theta levels (m)
real(kind=real_umphys), intent(in) :: r_rho(n_sh,nlev)
                                            ! radius of rho levels (m)
real(kind=real_umphys), intent(in)    ::                                       &
  rho_theta(n_sh,nlev)      & ! wet density for theta lev (kg/m3)
 ,rho_dry_theta(n_sh,nlev)  & ! dry density on theta levels (kg/m3)
 ,rho_dry(n_sh,nlev)          ! dry density on rho levels (kg/m3)

real(kind=real_umphys), intent(in) :: r2rho_th(n_sh,nlev)
                                            ! radius**2 density for
                                            ! theta lev (kg/m)
real(kind=real_umphys), intent(in) :: r2rho(n_sh,nlev)
                                            ! radius**2 density for
                                            ! rho lev (kg/m)
real(kind=real_umphys), intent(in) :: dr_across_th(n_sh,nlev)
                                            ! thickness of theta levels (m)
real(kind=real_umphys), intent(in) :: dr_across_rh(n_sh,nlev)
                                            ! thickness of rho levels (m)
real(kind=real_umphys), intent(in) :: conv_prog_precip(n_sh,nlev)
                                                ! Surface precipitation based
                                                ! 3d convective prognostic in
                                                ! kg/m2/s
real(kind=real_umphys), intent(in)    :: q(n_sh,nlev)
                                    ! Model mixing ratio (kg/kg)

real(kind=real_umphys), intent(in)    :: th(n_sh,nlev)
                                     ! Model potential temperature (K)

real(kind=real_umphys), intent(in)    :: timestep    ! Model timestep (s)

real(kind=real_umphys), intent(in)    ::                                       &
  u(n_sh,nlev)         & ! Model U field (m/s)
 ,v(n_sh,nlev)           ! Model V field (m/s)

real(kind=real_umphys), intent(in)    ::                                       &
  uw0(n_sh)            & ! U-comp of surface stress (N/m2)
 ,vw0(n_sh)            & ! V-comp of surface stress (N/m2)
 ,wstar(n_sh)          & ! Convective velocity scale (m/s)
 ,wthvs(n_sh)          & ! Surface flux of THV (Pa m/s2)
 ,entrain_coef(n_sh)   & ! Entrainment coefficients
 ,delta_smag(n_sh)     & ! grid size (m)
 ,zlcl_uv(n_sh)        & ! Lifting condensation level defined for the uv
                         ! grid (m)
 ,tnuc_nlcl(n_sh)      & ! nucleation temperature as function of dust
                         ! indexed using nlcl (deg cel)
 ,ztop_uv(n_sh)          ! Top of cloud layer defined for the uv grid (m)

integer, intent(in) :: freeze_lev(n_sh) ! Level index for freezing level

real(kind=real_umphys), intent(in) ::                                          &
  recip_pstar(n_sh)    & ! Reciprocal of pstar array
 ,qse(n_sh,nlev)         ! Saturation mixing ratio of cloud environment (kg/kg)

real(kind=real_umphys), intent(in) ::                                          &
  ccp_strength(n_sh)     ! Cold-pool strength

! Arguments with intent INOUT:

real(kind=real_umphys), intent(in out) ::                                      &
  cf_frozen(n_sh,nlev)   & ! Frozen water cloud volume ( )
 ,cf_liquid(n_sh,nlev)   & ! Liq water cloud volume ( )
 ,qcf(n_sh,nlev)         & ! Ice condensate mix ratio (kg/kg)
 ,qcl(n_sh,nlev)         & ! Liq condensate mix ratio (kg/kg)
 ,tracer(n_sh,trlev,ntra)  ! Model tracer fields (kg/kg)

real(kind=real_umphys), intent(in out) ::                                      &
  w2p(n_sh,nlev)           ! (Parcel vertical velocity)^2, [(m/s)^2]

real(kind=real_umphys), intent(in out) :: conv_prog_flx(n_sh,nlev)
                                                 ! Mass flux convective
                                                ! prognostic in Pa/s

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                               ! Structure containing water
                                               ! tracer fields

! Structure containing SCM convection sub-step diagnostics
! (needs intent inout as contains allocatable arrays that need to
! retain their allocated status on input as well as output)
type(scm_convss_dg_type), intent(in out) :: scm_convss_dg( n_sh )
! Flag for SCM convection sub-step diagnostics
logical, intent(in) :: l_scm_convss_dg


! Arguments with intent out:

real(kind=real_umphys), intent(out) ::                                         &
  cape_out(n_sh)     & ! Saved convective available potential energy for
                       ! diagnostic output (J/kg)
 ,cclwp(n_sh)        & ! Condensed water path (kg/m2)
 ,ccw(n_sh,nlev)     & ! Convective cloud liquid water on model levels (kg/kg)
 ,cca(n_sh,n_cca_lev)  ! Convective cloud amount on model levels (fraction)

real(kind=real_umphys), intent(out) ::                                         &
  dbcfbydt(n_sh,nlev)  & ! Increments to total cld volume due to convection(/s)
 ,dcffbydt(n_sh,nlev)  & ! Increments to ice cloud volume due to convection (/s)
 ,dcflbydt(n_sh,nlev)  & ! Increments to liq cloud volume due to convection (/s)
 ,dqbydt(n_sh,nlev)    & ! Increments to q due to convection (kg/kg/s)
 ,dqcfbydt(n_sh,nlev)  & ! Increments to ice condensate due to convection
                         ! (kg/kg/s)
 ,dqclbydt(n_sh,nlev)  & ! Increments to liq condensate due to convection
                         ! (kg/kg/s)
 ,dthbydt(n_sh,nlev)   & ! Increments to potential temp. due to convection (K/s)
 ,dubydt(n_sh,nlev)    & ! Increments to U due to CMT (m/s2)
 ,dvbydt(n_sh,nlev)      ! Increments to V due to CMT (m/s2)

real(kind=real_umphys), intent(out) ::                                         &
  dtrabydt(n_sh,nlev,ntra)  ! Increment to tracer due to convection (kg/kg/s)


real(kind=real_umphys), intent(out) ::                                         &
  detrain_up(n_sh,nlev)  & ! Fractional detrainment rate into updraughts (Pa/s)
 ,detrain_dwn(n_sh,nlev) & ! Fractional detrainment rate into downdraughts
                           ! (Pa/s)
 ,entrain_up(n_sh,nlev)  & ! Fractional entrainment rate into updraughts (Pa/s)
 ,entrain_dwn(n_sh,nlev)   ! Fractional entrainment rate into downdraughts
                           ! (Pa/s)

integer, intent(out) ::                                                        &
  iccb(n_sh)            & ! Convective cloud base level
 ,icct(n_sh)              ! Convective cloud top level

real(kind=real_umphys), intent(out) :: lcca(n_sh) ! Lowest conv. cloud amt. (%)

integer, intent(out) ::                                                        &
  lcbase(n_sh)          & ! Lowest conv. cloud base level
 ,lctop(n_sh)             ! Lowest conv. cloud top level

real(kind=real_umphys), intent(out) :: rain(n_sh)
                                ! Surface convective rainfall (kg/m2/s)

real(kind=real_umphys), intent(out) :: snow(n_sh)
                                ! Surface convective snowfall (kg/m2/s)

real(kind=real_umphys), intent(out) :: rain_3d(n_sh,nlev)
                                        ! Convective rainfall flux (kg/m2/s)

real(kind=real_umphys), intent(out) :: snow_3d(n_sh,nlev)
                                        ! Convective snowfall flux (kg/m2/s)

real(kind=real_umphys), intent(out) :: up_flux(n_sh,nlev)
                                        ! Updraught mass flux (Pa/s)

real(kind=real_umphys), intent(out) :: up_flux_half(n_sh,nlev)
                                             ! Updraught mass flux on
                                             !half levels(Pa/s)

real(kind=real_umphys), intent(out) :: dwn_flux(n_sh,nlev)
                                         ! Downdraught mass flux (Pa/s)

real(kind=real_umphys), intent(out) :: uw_shall(n_sh,nlev) ! X-comp. of stress
                                         ! from shallow convection (kg/m/s2)

real(kind=real_umphys), intent(out) :: vw_shall(n_sh,nlev) ! Y-comp. of stress
                                         ! from shallow convection (kg/m/s2)

real(kind=real_umphys), intent(out) :: tcw(n_sh)
                                ! Total condensed water(kg/m2/s)

real(kind=real_umphys), intent(out) :: cca_2d(n_sh)
                                  ! 2D convective cloud amount (%)

real(kind=real_umphys), intent(out) :: ind_shallow(n_sh)
                                        ! 1.0 if real shallow convection
                                        ! else 0.0

integer, intent(out) :: kterm(n_sh) ! termination level for shallow
                                    ! convection
! Downdraught and evap below cloud base
real(kind=real_umphys), intent(out) ::                                         &
  dt_dd(n_sh,nlev)        & ! dT/dt from DD and evap below cloud base (K/s)
 ,dq_dd(n_sh,nlev)        & ! dq/dt from DD and evap below cloud base (kg/kg/s)
 ,area_ud(n_sh,nlev)        ! fractional updraught area

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

integer :: index1(n_sh),index2(n_sh)

integer :: ncposs               ! No. of points which may convect

integer :: nconv                ! No. of convecting points

real(kind=real_umphys) :: amdetk(n_sh)
                                ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

real(kind=real_umphys) :: b_calc                  ! Coefficient in thpert calc.

real(kind=real_umphys) :: c_calc                  ! Coefficient in thpert calc.

real(kind=real_umphys) :: cape(n_sh)
                                ! Convective available potential
                                ! energy (J/kg)
real(kind=real_umphys) :: fcape(n_sh)
                                ! Convective available potential
                                ! energy weighted by f.det profile (J/kg)

real(kind=real_umphys) :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

real(kind=real_umphys) :: dcpbydt(n_sh)
                                ! Rate of change of cape (J/kg/s)

real(kind=real_umphys) :: depth(n_sh)
                                ! Depth of convective cloud (m)

real(kind=real_umphys) :: ekp14(n_sh)             ! Entrainment coefficients at
                                ! level k+1/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

real(kind=real_umphys) :: ekp34(n_sh)             ! Entrainment coefficients at
                                ! level k+3/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

real(kind=real_umphys) :: exk(n_sh)               ! Exner ratio at layer k

real(kind=real_umphys) :: exkp1(n_sh)             ! Exner ratio at layer k+1

real(kind=real_umphys) :: flxmax(n_sh)            ! Maximum initial convective
                                ! mass flux (Pa/s)

real(kind=real_umphys) :: flx_init(n_sh)
                                ! Initial mass flux at cloud base
                                ! (Pa/s)

real(kind=real_umphys) :: flx_init_new(n_sh)      ! flx_init scaled (Pa/s)

real(kind=real_umphys) :: flxmax_init(n_sh)
                                ! Maximum possible initial mass
                                ! flux (limited to the mass in
                                ! the initial convecting layer
                                ! in Pa/s)

real(kind=real_umphys) :: max_cfl(n_sh)
                                ! Max cfl ratio over a convecting
                                ! layer

real(kind=real_umphys) :: decay_amount
                                ! decay fraction for time-smoothed mass flux

real(kind=real_umphys) :: p_lcl(n_sh)             ! Pressure at LCL (Pa)

real(kind=real_umphys) :: precip(n_sh,nlev)
                                ! Amount of precip from each layer
                                ! from each layer (kg/m/s)

real(kind=real_umphys) :: pk(n_sh)
                                ! Pressure at midpoint of layer
                                ! k (Pa)

real(kind=real_umphys) :: pkp1(n_sh)
                                ! Pressure at midpoint of layer
                                ! k+1 (Pa)

real(kind=real_umphys) :: delpk(n_sh)
                                ! Pressure difference over layer
                                ! k (Pa)

real(kind=real_umphys) :: delpkp1(n_sh)
                                ! Pressure difference over layer
                                ! k+1 (Pa)

real(kind=real_umphys) :: delpkp12(n_sh)          ! Pressure difference between
                                ! layers k and k+1 (Pa)

real(kind=real_umphys) :: delp_uv_k(n_sh)
                                ! Pressure difference across uv
                                ! layer k (Pa)

real(kind=real_umphys) :: delp_uv_kp1(n_sh)
                                ! Pressure difference across uv
                                ! layer k+1 (Pa)

real(kind=real_umphys) :: q_lcl(n_sh)             ! Mixing ratio at LCL (kg/kg)

real(kind=real_umphys) :: qse_lcl(n_sh)           ! Saturated q at LCL (kg/kg)

real(kind=real_umphys) :: rhum(n_sh)              ! Dummy relative humidity
                                ! (only used on shallow points)

real(kind=real_umphys) :: t_lcl(n_sh)             ! Temperature at LCL (K)

real(kind=real_umphys) :: th_lcl(n_sh)            ! Theta at LCL (K)

real(kind=real_umphys) :: dthv_ma                 ! Moist adiabtic change in thv
                                ! from ntml to ntml+1 (K)

real(kind=real_umphys) :: thv_pert(n_sh)
                                ! Theta_v parcel pertubation (K)

real(kind=real_umphys) :: thpert(n_sh)            ! Theta parcel pertubation (K)

real(kind=real_umphys) :: qpert(n_sh)             ! q parcel pertubation (kg/kg)

real(kind=real_umphys) :: rho_k ! density on level k

integer :: start_lev(n_sh)      ! Convection initiation level

real(kind=real_umphys) :: wsc(n_sh)
                                ! Convective velocity scale (m/s)

real(kind=real_umphys) :: wsc_o_mb(n_sh)
                                ! Convective velocity scale divided
                                ! by cloud base mass flux mb
real(kind=real_umphys) :: w_max(n_sh)
                                ! dummy variable for maximum w in column (m/s)

logical :: bgmk(n_sh)           ! Mask for points where parcel in
                                ! layer k is saturated

logical :: blatent(n_sh)        ! Mask for points where latent heat has
                                ! been released

logical :: bwater(n_sh,2:nlev)  ! Mask for points at which
                                ! condensate is liquid

logical :: blowst(n_sh)         ! Dummy variable indicating low
                                ! enough stability for convection
                                ! to occur

logical :: bterm(n_sh)          ! Mask for points which have
                                ! stopped convecting

logical :: bconv(n_sh)          ! Mask for points at which
                                ! convection is occurring

logical :: bcposs(n_sh)         ! Mask for points passing
                                ! initial stability test


! Parcel variables


real(kind=real_umphys) :: qpi   ! Initial parcel mixing ratio
                                !(kg/kg)

real(kind=real_umphys) :: qp(n_sh,nlev)           ! Parcel mixing ratio (kg/kg)

real(kind=real_umphys) :: thpi  ! Initial parcel potential temp.
                                !(K)

real(kind=real_umphys) :: thp(n_sh,nlev)          ! Parcel potential temp (K)

real(kind=real_umphys) :: up(n_sh,nlev)           ! Parcel U (m/s)

real(kind=real_umphys) :: vp(n_sh,nlev)           ! Parcel V  (m/s)

real(kind=real_umphys) :: trap(n_sh,nlev,ntra)    ! Tracer content of parcel
                                ! (kg/kg)

real(kind=real_umphys) :: flx(n_sh,nlev)          ! Parcel massflux (Pa/s)

real(kind=real_umphys) :: xsbmin_v(n_sh,nlev)
                                ! Minmum parcel buoyancy excess

real(kind=real_umphys) :: thpixs_v(n_sh,nlev)     ! Theta parcel excess (K)

real(kind=real_umphys) :: qpixs_v(n_sh,nlev)      ! Q parcel excess(kg/kg)

type(conv_p_wtrac_type) :: wtrac_p(n_wtrac)
                                ! Structure containing water
                                ! tracer fields relating to the parcel
real(kind=real_umphys), allocatable :: qpixs_v_wtrac(:,:,:)
                                ! Water tracer q parcel excess (kg/kg)
real(kind=real_umphys), allocatable :: ratio_wtrac(:,:,:)
                                ! Ratio of water tracer to normal water

real(kind=real_umphys) :: qclp(n_sh,nlev)
                                ! Parcel liquid condensated mixing
                                ! ratio in layer k (kg/kg)

real(kind=real_umphys) :: qcfp(n_sh,nlev)
                                ! Parcel frozen condensated mixing
                                ! ratio in layer k (kg/kg)

! Undilute parcel variables
real(kind=real_umphys) :: qu(n_sh,nlev)
                                ! Undilute Parcel mixing ratio (kg/kg)
real(kind=real_umphys) :: thu(n_sh,nlev)
                                ! Undilute Parcel potential temp (K)

! Parameters

real(kind=real_umphys), parameter :: cfl_limit = 1.0 ! Max CFL ratio allowed
real(kind=real_umphys), parameter :: minflx = tiny(flx_init_new)
                                                ! minimum allowable
                                                ! initial mass flux

! CMT variables

integer :: nlcl_uv(n_sh)        ! Level index for LCL

integer :: ntop_uv(n_sh)        ! Level index for top of layer

integer :: n_0degc(n_sh)        ! Level index for zero degrees

integer :: cu_term(n_sh),cu_tend(n_sh) !Indicies for CMT subs

real(kind=real_umphys) :: exk_temp                ! Temporary exner

real(kind=real_umphys) :: eflux_u_ud(n_sh)
                                ! Vertical eddy flux of momentum
                                ! due to UD at top of layer
                                ! (Pa m/s2)

real(kind=real_umphys) :: eflux_v_ud(n_sh)
                                ! Vertical eddy flux of momentum
                                ! due to UD at bottom of layer
                                ! (Pa m/s2)

real(kind=real_umphys) :: mb(n_sh)                ! Cloud base mass flux (Pa/s)

real(kind=real_umphys) :: p_uv(nlev,n_sh)         ! Pressure of model level (Pa)

real(kind=real_umphys) :: phalf_uv(nlev,n_sh)     ! Pressure of half level (Pa)

real(kind=real_umphys) :: plcl_uv(n_sh)           ! Pressure at LCL (Pa)

real(kind=real_umphys) :: ptop_uv(n_sh)
                                ! Pressure at top of cloud layer
                                ! (Pa)

real(kind=real_umphys) :: rho_uv(nlev,n_sh)       ! Density on uv level (kg/m3)

real(kind=real_umphys) :: uw(nlev,n_sh)
                                ! U- comp stress profile (N/m2)
                                ! (units change through calls)

real(kind=real_umphys) :: ue_p(nlev,n_sh)         ! Environment U profile (m/s)

real(kind=real_umphys) :: vw(nlev,n_sh)           ! V-comp stress profile (N/m2)

real(kind=real_umphys) :: ve_p(nlev,n_sh)         ! Environment V profile (m/s)

real(kind=real_umphys) :: zcld(n_sh)              ! Depth of cloud layer (m)

logical :: l_mom_gk             ! true if Gregory-Kershaw CMT
logical :: l_mom_gk_stable      ! true if stabilized Gregory-Kershaw CMT
                                ! (different from the original)

! CFL scaling variables


integer :: det_lev(n_sh)        ! Level at which split final
                                ! detrainment last occurred

integer :: nterm                ! No. of points where conv.
                                ! has terminated

integer :: index_nterm(n_sh)    ! Index for points where conv.
                                ! has terminated

real(kind=real_umphys) :: tempnum
                                ! Temporary variable for storage

real(kind=real_umphys) :: scale_f(n_sh)           ! store scaling factor

real(kind=real_umphys) :: weight_param(n_sh)      ! Weighting factor

real(kind=real_umphys) :: deltaktot(n_sh)
                                ! Integrated forced detrainment

! Original downdraught scheme variables
integer :: nnodd                ! No. of downdraughts not possible

integer :: index_nodd(n_sh)     ! Index of downdraughts not
                                ! possible
integer :: npossdd              ! No. downdraughts possible

integer :: index_possdd(n_sh)   ! Index of downdraughts

integer :: kmax_term            ! maximum termination level + 1

real(kind=real_umphys) :: deltap_cld
                                ! pressure thickness of convective
                                ! cloud (Pa)

! Limit nlev loop to those levels actually required using ntpar
! diagnosed in conv_diag

integer :: ntpar_max            ! max ntpar value

! parameter for qmin checks
real(kind=real_umphys), parameter :: qmin = 1.0e-8 ! Global minimum allowed Q

logical ::                                                                     &
  bgmkp1(n_sh)

real(kind=real_umphys) ::                                                      &
   qpkp1(n_sh)                                                                 &
  ,qclpkp1(n_sh)                                                               &
  ,qcfpkp1(n_sh)                                                               &
  ,thpkp1(n_sh)

real(kind=real_umphys) :: Qlkp1(n_sh)
real(kind=real_umphys) :: Qfkp1(n_sh)
real(kind=real_umphys) :: Frezkp1(n_sh)

real(kind=real_umphys) :: watldek(n_sh)
real(kind=real_umphys) :: watldpk(n_sh)
real(kind=real_umphys) :: watldekp1(n_sh)
real(kind=real_umphys) :: watldpkp1(n_sh)

real(kind=real_umphys) :: dummy1(n_sh)
real(kind=real_umphys) :: dummy2(n_sh)

! Local compressed arrays

real(kind=real_umphys) :: thrk(n_sh)
                              ! potential temperature of forced detrained air

real(kind=real_umphys) :: qrk(n_sh)
                              ! specific humidity of forced detrained air

real(kind=real_umphys) :: deltak(n_sh)          ! Parcel forced detrainment rate
                              ! in layer k multiplied by
                              ! appropriate layer thickness

real(kind=real_umphys) :: flxkp12_t(n_sh)       ! Half level mass flux (Pa/s)

real(kind=real_umphys) :: rbuoyk(n_sh)      ! Par. buoyancy at k (K)
real(kind=real_umphys) :: rbuoykp1(n_sh)    ! Par. buoyancy at k+1 (K)
real(kind=real_umphys) :: rbuoyukp1(n_sh)   ! undilute Par. buoy at k+1 (K)

real(kind=real_umphys) :: qsat_lcl(n_sh)         ! not used
real(kind=real_umphys) :: z_scale                ! cloud depth scale

logical :: b_nodd(n_sh)   ! points with no downdraught
logical :: b_dd(n_sh)     ! points with downdraught on termination

!===============================================================
! CCRad Variables local variables
!===============================================================

real(kind=real_umphys)   :: overlap_fac(n_sh)  ! Factor designed to improve
                             ! shallow Cu cover by allowing
                             ! for non-vertical clouds.

real(kind=real_umphys)   :: zpr         ! method (BL Fluxes) only if
                      !   l_ccrad = T .and. cca2d_sh_opt  = 1

!===============================================================
! End CCRad Variables local variables
!===============================================================

! Loop counters

integer :: i,i2,j,k,ktra,kt,i_wt

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SHALLOW_CONV_6A'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

! Initialise logicals

do i = 1,n_sh
  blowst(i)    = .true.
  bconv(i)     = .false.
  bcposs(i)    = .false.
  b_nodd(i)    = .false.
  b_dd(i)      = .false.
  blatent(i)   = .false.
end do

do i = 1, n_sh
  kterm(i)        = 0
  ind_shallow(i)  = 0.0
  start_lev(i)    = ntml(i)
end do

l_mom_gk = .false.       ! not Gregory-Kershaw CMT
l_mom_gk_stable = .false. !
                         ! Shallow code uses turbulence based CMT

!-----------------------------------------------------------------------
! 2.1  Initialise parcel properties and increment arrays
!-----------------------------------------------------------------------

!initialise parcel values over all levels
do k = 1, nlev
  do i = 1, n_sh
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

if (sh_new_termc == term_undil) then
  do k = 1, nlev
    do i = 1, n_sh
      qu(i,k)     = 0.0
      thu(i,k)    = 0.0
    end do
  end do
end if

if (l_mom_gk) then
  do k=1,nlev
    do i = 1,n_sh
      up(i,k) = 0.0
      vp(i,k) = 0.0
    end do
  end do
end if

if (l_tracer) then
  do ktra = 1,ntra
    do k=1,nlev
      do i = 1,n_sh
        trap(i,k,ktra) = 0.0
      end do
    end do
  end do
end if

! Allocate water tracer arrays and initialise

call wtrac_alloc_conv_p(n_sh, nlev, n_wtrac, wtrac_p)

if (l_wtrac_conv) then

  do i_wt = 1, n_wtrac
    do k = 1, nlev
      do i = 1, n_sh
        wtrac_p(i_wt)%q(i,k)      = 0.0
        wtrac_p(i_wt)%qcl(i,k)    = 0.0
        wtrac_p(i_wt)%qcf(i,k)    = 0.0
        wtrac_p(i_wt)%precip(i,k) = 0.0
      end do
    end do
  end do
end if  !l_wtrac_conv

do k = 1,nlev
  do i = 1,n_sh
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
    do i = 1,n_sh
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    end do
  end do
end if  ! L_mom

if (l_tracer) then
  do ktra = 1,ntra
    do k = 1,nlev
      do i = 1,n_sh
        dtrabydt(i,k,ktra) = 0.0
      end do
    end do
  end do
end if  ! L_tracer

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do k = 1, nlev
      do i = 1, n_sh
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
if (flg_up_flx .or. flg_mf_shall) then
  do k = 1,nlev
    do i = 1,n_sh
      up_flux(i,k)      = 0.0
    end do
  end do
end if
if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,n_sh
      up_flux_half(i,k) = 0.0
    end do
  end do
end if
if (flg_dwn_flx) then
  do k = 1,nlev
    do i = 1,n_sh
      dwn_flux(i,k)     = 0.0
    end do
  end do
end if
if (flg_entr_up) then
  do k = 1,nlev
    do i = 1,n_sh
      entrain_up(i,k)   = 0.0
    end do
  end do
end if
if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,n_sh
      detrain_up(i,k)   = 0.0
    end do
  end do
end if
if (flg_entr_dwn) then
  do k = 1,nlev
    do i = 1,n_sh
      entrain_dwn(i,k)   = 0.0
    end do
  end do
end if
if (flg_detr_dwn) then
  do k = 1,nlev
    do i = 1,n_sh
      detrain_dwn(i,k)  = 0.0
    end do
  end do
end if
if (l_mom) then
  if (flg_uw_shall) then
    do k = 1,nlev
      do i = 1,n_sh
        uw_shall(i,k)   = 0.0
      end do
    end do
  end if
  if (flg_vw_shall) then
    do k = 1,nlev
      do i = 1,n_sh
        vw_shall(i,k)   = 0.0
      end do
    end do
  end if
end if  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
do i = 1,n_sh
  cca_2d(i) = 0.0
  iccb(i)   = 0
  icct(i)   = 0
  tcw(i)    = 0.0
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
end do

do k = 1,n_cca_lev
  do i = 1,n_sh
    cca(i,k) = 0.0
  end do
end do

do i = 1,n_sh
  !-----------------------------------------------------------------------
  ! 2.5  Initialise diagnostics for scaling calculations
  !-----------------------------------------------------------------------
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)         = 0.0
  fcape(i)        = 0.0
  cape_out(i)     = 0.0
  dcpbydt(i)      = 0.0
  max_cfl(i)      = 0.0
  det_lev(i)      = 0
  deltaktot(i)    = 0.0

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
end do

! Initialise surface precip arrays for water tracers
if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i = 1,n_sh
      wtrac_e(i_wt)%rain(i) = 0.0
      wtrac_e(i_wt)%snow(i) = 0.0
    end do
  end do
end if


! Note: In terms of array indices p and phalf follow the convention
!       used in the boundary layer scheme. phalf(k,*) refers to the
!       lower boundary of uv layer k. This follows the convention for
!       um UM4.5 and before
!
!       Also note that p_layer_boundaries(0) and p_layer_centres(0)
!       = pstar, so p_uv(k,1) and phalf_uv(k,1) will be equal.
!
!       Because of the definition of nlcl, the pressure of the top of
!       the mixed layer is phalf_uv(nlcl,*)

if (l_mom) then

  ! Initialize arrays required for Convective Momentum Transport(CMT)

  k=1
  do i = 1,n_sh
    p_uv(k,i)     = p_layer_boundaries(i,k-1)
    phalf_uv(k,i) = p_layer_centres(i,k-1)
    ue_p(k,i)     = u(i,k)
    ve_p(k,i)     = v(i,k)
  end do

  do i = 1,n_sh
    nlcl_uv(i)    = ntml(i) + 1
    ntop_uv(i)    = ntpar(i) + 1
    n_0degc(i)    = freeze_lev(i)
  end do

  do i = 1,n_sh
    do k = 2,nlev
      p_uv(k,i)     = p_layer_boundaries(i,k-1)
      phalf_uv(k,i) = p_layer_centres(i,k-1)
      ue_p(k,i)     = u(i,k)
      ve_p(k,i)     = v(i,k)
      exk_temp      = (p_uv(k,i)/pref)**kappa
      rho_uv(k,i)   = 2.0 * p_uv(k,i) / (r * exk_temp * (th(i,k-1) + th(i,k)))
    end do
    plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
    ptop_uv(i)      = phalf_uv(ntop_uv(i),i)
    rho_uv(1,i)     = rho_uv(2,i)
  end do
end if     !L_mom

!-----------------------------------------------------------------------
! Calculate parcel perturbations
!-----------------------------------------------------------------------

! Calculate XSBMIN and THPIXS constants based on layer thickness (Pa)
do k = 1,nlev-1
  do i = 1,n_sh
    xsbmin_v(i,k) = min( ((p_layer_centres(i,k) -                              &
              p_layer_centres(i,k+1))/5000.0),1.0) *0.2

    thpixs_v(i,k) = min( ((p_layer_centres(i,k) -                              &
              p_layer_centres(i,k+1))/5000.0),1.0) * thpixs_shallow

    qpixs_v(i,k)  = qpixs_shallow
  end do
end do  ! nlev

if (l_wtrac_conv) then

  allocate(qpixs_v_wtrac(n_sh,nlev,n_wtrac))
  allocate(ratio_wtrac(n_sh,nlev,n_wtrac))

  do i_wt = 1, n_wtrac
    do k = 1,nlev-1
      do i = 1,n_sh
        ! Calculate water tracer content of parcel q perturbation
        ratio_wtrac(i,k,i_wt) = wtrac_calc_ratio_fn(i_wt,                      &
                                wtrac_e(i_wt)%q(i,k), q(i,k))
        qpixs_v_wtrac(i,k,i_wt) = ratio_wtrac(i,k,i_wt) * qpixs_v(i,k)
      end do
    end do  ! nlev
  end do    ! n_wtrac
end if ! l_wtrac_conv


! Calculate convective velocity scale and cloud base mass flux
do i = 1,n_sh
  wsc(i) = (delthvu(i) * c_mass * wstar(i) * g / (th(i,ntml(i))                &
               * (1.0 + c_virtual * q(i,ntml(i)))))**0.3333
  mb(i)  = c_mass_sh * wstar(i)
  zcld(i) = ztop_uv(i) - zlcl_uv(i)
  wsc_o_mb(i) = wsc(i)/mb(i)
  weight_param(i) = 1.0
end do

! Define the LCL

if ( sh_pert_opt == 0) then

  !         Define the LCL at the half level above ntml. Find
  !         environmental T at p_lcl by approximating theta there with
  !         th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is
  !         tunable.  Similarly for q.

  do i = 1,n_sh
    k =ntml(i)
    p_lcl(i)  = p_layer_boundaries(i,k)
    th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
    t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
    q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
  end do

else  ! Sh_pert_Opt = 1

  !         Define the LCL at the half level above ntml. Find
  !         environmental T at p_lcl by approximating theta there with
  !         th(i,k) Similarly for q.

  do i = 1,n_sh
    k =ntml(i)
    p_lcl(i)  = p_layer_boundaries(i,k)
    th_lcl(i) = th(i,k)
    t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
    q_lcl(i)  = q(i,k)
  end do

end if

! Calculate saturation mixing ratio at LCL
call qsat(qse_lcl,t_lcl,p_lcl,n_sh)

! Calculate theta and q perturbation (perturbation is based on
! environment buoyancy gradient)
! Reset thpixs and qpixs at ntml
if ( sh_pert_opt == 0) then
  do i = 1,n_sh

    k = ntml(i)

    if (t_lcl(i) >  tm) then
      dq_sat_env = repsilon * lc * qse_lcl(i)                                  &
                      / (r * t_lcl(i) * t_lcl(i))

      ! Estimate of moist adiabatic lapse rate

      dthv_ma    = ( (lc/cp) - (1.0+c_virtual)*th(i,k) )*                      &
                   dq_sat_env*(g/cp)/(1.0+(lc/cp)*dq_sat_env)
    else
      dq_sat_env = repsilon * (lc+lf) * qse_lcl(i)                             &
                      / (r * t_lcl(i) * t_lcl(i))

      ! Estimate of moist adiabatic lapse rate (in K/m)

      dthv_ma    = ( ((lc+lf)/cp) - (1.0+c_virtual)*th(i,k) )*                 &
                   dq_sat_env*(g/cp)/(1.0+((lc+lf)/cp)*dq_sat_env)
    end if

    b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0                         &
                        + c_virtual * qse_lcl(i)

    ! Calculate theta_v perturbation:

    if (kprof_cu == off) then
      thv_pert(i) = -0.17 * wthvs(i) / mb(i)                                   &
                    + (th(i,k+1) * (1.0 + c_virtual                            &
                    * q(i,k+1)) - th(i,k)                                      &
                    * (1.0 + c_virtual * q(i,k)))
    else
      ! "entrainment flux" at LCL given by BL scheme
      thv_pert(i) =  (th(i,k+1) * (1.0 + c_virtual                             &
                    * q(i,k+1)) - th(i,k)                                      &
                    * (1.0 + c_virtual * q(i,k)))
    end if

    c_calc   = th_lcl(i) * c_virtual * (qse_lcl(i) - q_lcl(i))                 &
                         - thv_pert(i)

    thpert(i) = -c_calc / b_calc  ! ignore term in thpert**2

    thpixs_v(i,k) = thpert(i)

    qpert(i)  = qse_lcl(i) + ((p_lcl(i) / pref)                                &
                           **kappa) * thpert(i) * dq_sat_env                   &
                           - q_lcl(i)

    qpixs_v(i,ntml(i))  = qpert(i)

  end do !n_sh

else ! Sh_pert_opt = 1

  do i = 1,n_sh

    k = ntml(i)

    if (t_lcl(i) >  tm) then
      dq_sat_env = repsilon * lc * qse_lcl(i)                                  &
                      / (r * t_lcl(i) * t_lcl(i))

      !         Estimate of moist adiabatic lapse rate

      dthv_ma    = ( (lc/cp) - (1.0+c_virtual)*th(i,k) )*                      &
                   dq_sat_env*(g/cp)/(1.0+(lc/cp)*dq_sat_env)
    else
      dq_sat_env = repsilon * (lc+lf) * qse_lcl(i)                             &
                      / (r * t_lcl(i) * t_lcl(i))

      !         Estimate of moist adiabatic lapse rate (in K/m)

      dthv_ma    = ( ((lc+lf)/cp) - (1.0+c_virtual)*th(i,k) )*                 &
                   dq_sat_env*(g/cp)/(1.0+((lc+lf)/cp)*dq_sat_env)
    end if

    b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0                         &
                        + c_virtual * qse_lcl(i)


    !       Calculate theta_v perturbation:

    !       First convert moist adiabatic lapse rate to thv difference
    !       between levels k and k+1

    rho_k   = p_layer_centres(i,k) /                                           &
            (r * th(i,k) * (p_layer_centres(i,k)/ pref)**kappa)
    dthv_ma = -dthv_ma*                                                        &
            (p_layer_centres(i,k+1)-p_layer_centres(i,k)) /                    &
            (rho_k*g)

    !       Make perturbation relative to a target lapse rate (namely
    !       0.6*dthv_ma, which is approximately what is seen in LES)

    if (kprof_cu == off) then
      thv_pert(i) = -0.17 * wthvs(i) / mb(i) +  0.6*dthv_ma                    &
                - (  th(i,k+1)*(1.0 + c_virtual*q(i,k+1))                      &
                - th(i,k)  *(1.0 + c_virtual*q(i,k))  )
    else
      ! "entrainment flux" at LCL done by BL scheme
      thv_pert(i) = 0.6*dthv_ma                                                &
                    - (  th(i,k+1)*(1.0 + c_virtual*q(i,k+1))                  &
                    - th(i,k)  *(1.0 + c_virtual*q(i,k))  )

    end if

    !       limit thv_pert to physically sensible values

    thv_pert(i) = max(min(thv_pert(i), max_sh_thpert),                         &
                  min_sh_thpert)

    c_calc      = th_lcl(i) * c_virtual * (qse_lcl(i) - q_lcl(i))              &
                - thv_pert(i)

    thpert(i)   = max(min(-c_calc / b_calc, max_sh_thpert),                    &
                  min_sh_thpert)  ! ignore term in thpert**2

    thpixs_v(i,k) = thpert(i)

    qpert(i)    = max(min(qse_lcl(i) + ((p_lcl(i) / pref)                      &
                **kappa) * thpert(i) * dq_sat_env                              &
                - q_lcl(i),                                                    &
                  max_sh_qpert_fac * qse_lcl(i)),0.0)

    qpixs_v(i,ntml(i))  = qpert(i)

  end do !n_sh

end if

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    do i = 1,n_sh
      wtrac_p(i_wt)%qpert(i)        = ratio_wtrac(i,ntml(i),i_wt) * qpert(i)
      qpixs_v_wtrac(i,ntml(i),i_wt) = wtrac_p(i_wt)%qpert(i)
    end do
  end do
end if ! l_wtrac_conv


! Set bwater=.true. on points where water will condense rather than
! ice.
call flag_wet(n_sh,n_sh,nlev,th,exner_layer_centres,bwater)


!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------
! To reduce cost limit level loop to NTPAR_MAX maximum NTPAR value.
! NTPAR is the top of the parcel ascent for shallow convection.

if (sh_on == 1 .or. cldbase_opt_sh == sh_grey_closure) then
                     ! Adaptive forced detrainment or grey shallow param
                     ! No limit on convection top
  ntpar_max = nlev-3 ! What is a sensible value to have here?

else                 ! Top limited

  ntpar_max=0
  do i = 1,n_sh
    if (ntpar(i)+1 >  ntpar_max) then
      ntpar_max=ntpar(i)+1
    end if
  end do
  ! Ensure that ntpar_max does not exceed nlev-1
  ntpar_max = min(ntpar_max, nlev-1)
end if

do k = 2,ntpar_max  !loop over model levels

  do i = 1,n_sh
    bterm(i)  = .false.
  end do

  !-----------------------------------------------------------------------
  ! Initialise environment variables
  ! NB These variable are only used by layer_cn.
  !-----------------------------------------------------------------------
  do i = 1,n_sh
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
  do i = 1,n_sh
    if ( .not. bconv(i) .and. det_lev(i) == 0) then
      bgmk(i)     = .false.
      depth(i)    = 0.0
      thpi        = th(i,k) + thpixs_v(i,k)
      thp(i,k)    = thpi
      qpi         = q(i,k)  + qpixs_v(i,k)
      qp(i,k)     = qpi
      if (sh_new_termc == term_undil) then
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
    end if
  end do  ! n_sh
  if (l_tracer) then
    do ktra=1,ntra
      do i = 1,n_sh
        if ( .not. bconv(i)) then
          trap(i,k,ktra)  = tracer(i,k,ktra)
        end if  !not bconv
      end do
    end do
  end if

  ! Water tracers
  if (l_wtrac_conv) then

    do i_wt = 1, n_wtrac
      do i = 1,n_sh
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

  call layer_cn_6a(k, n_sh, nlev,                                              &
                   ntml, ntpar, start_lev,                                     &
                   exner_layer_centres,                                        &
                   p_layer_boundaries, p_layer_centres,                        &
                   z_rho,                                                      &
                   conv_prog_precip,                                           &
                   recip_pstar, entrain_coef, rhum,                            &
                   ccp_strength,                                               &
                   z_theta(:,k), z_rho(:,k+1), z_theta(:,k+1),                 &
                   wsc_o_mb, qsat_lcl, w_max,                                  &
                   shallow,                                                    &
                   bconv,                                                      &
                   ! Out
                   pk, pkp1, exk, exkp1,                                       &
                   delpk, delpkp12, delpkp1,                                   &
                   delp_uv_k, delp_uv_kp1,                                     &
                   ekp14, ekp34, amdetk                                        &
                   )

  ! Maximum initial convective mass flux
  do i = 1,n_sh
    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)
  end do

  !-----------------------------------------------------------------------
  ! Initial test to check if convection is possible in layer k
  !-----------------------------------------------------------------------
  ! Convection is possible if
  ! - the point was convecting (bconv = .T.) and did not terminate
  !   in the previous layer
  ! - or if at the top level of the surface mixed layer (k = ntml)
  do i = 1,n_sh
    bcposs(i) = bconv(i) .or. k  ==  ntml(i)
  end do  ! n_sh

  ! Calculate number of points which may convect (ncposs) and
  ! set compression indices (index1)
  ncposs = 0
  do i = 1,n_sh
    if (bcposs(i)) then
      ncposs          = ncposs + 1
      index1(ncposs)  = i
    end if
  end do

  !-----------------------------------------------------------------------
  ! 3.2  Lift parcel from layer k to layer k+1
  !-----------------------------------------------------------------------
  if (ncposs > 0) then

    call lift_par_6a (k, n_sh, n_wtrac, th(:,k), th(:,k+1),                    &
                q(:,k), q(:,k+1), qcl(:,k), qcf(:,k),                          &
                qcl(:,k+1), qcf(:,k+1),                                        &
                 pkp1(:), exkp1(:),                                            &
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
    call water_loading(n_sh, qcl(:,k), qcf(:,k), qclp(:,k), qcfp(:,k),         &
                     watldek, watldpk, index1, ncposs)

    call water_loading(n_sh, qcl(:,k+1), qcf(:,k), qclpkp1, qcfpkp1,           &
                     watldekp1, watldpkp1, index1, ncposs)


    if (sh_new_termc == term_undil) then
      !---------------------------------------------------------------------
      ! Calculate undilute parcel properties
      !---------------------------------------------------------------------
      call lift_undil_par_6a(n_sh,                                             &
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

      if (sh_new_termc == term_undil) then
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

        ! Store diagnostics linked to initial convective mass flux for
        ! calculation of final closure.
        flx_init(index1(i))    = flx(index1(i),k)
        flxmax_init(index1(i)) = flxmax(index1(i))

      else
        blowst(index1(i)) = .false. ! convection not initialise in layer
      end if

      ! Reset threshold for forced detrainment to the initial
      ! (potentially negative) buoyancy
      xsbmin_v(index1(i),k) = thv_pert(index1(i))

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
    !Compression for intent(in)
!DIR$ IVDEP
    do i = 1,nconv
      dummy1(index2(i)) = 0.0
      dummy2(index2(i)) = 0.0
    end do

    ! Force shallow convection to stop at the parcel top from conv_diag (ntpar)
    ! unless using adaptive forced detrainment or grez-zone shallow.
    if (sh_on == 0 .and. cldbase_opt_sh == sh_wstar_closure) then
!DIR$ IVDEP
      do i = 1,nconv
        if (k  ==  ntpar(index2(i))) then
          bterm(index2(i)) = .true.
        end if
      end do
    end if

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

    call convec2_6a  (k, n_sh, n_sh, nlev, ntra, n_wtrac, trlev,               &
                      sh_on, sh_new_termc, start_lev, timestep,                &
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
                      dummy1, dummy2, deltaktot, max_cfl,                      &
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
  do i = 1,n_sh
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
  end if      ! nconv > 0

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
  ! 3.4  CFL scaling and grey zone closures
  !-----------------------------------------------------------------------
  ! Set up integer nterm which is the total number of points where
  ! convection has terminated.
  ! Index to full array (n_sh) with index_nterm

  nterm = 0
  do i = 1,n_sh
    if (bterm(i)) then
      nterm = nterm + 1
      index_nterm(nterm) = i
    end if
  end do

  if (nterm >  0) then

    ! Work out scaled mass flux needed to keep cfl ratio below limit.
    ! First scale by grey zone parametrizations, if requested
    ! Then to keep cfl ratio below limit.
    ! Note L_CAPE not applied to shallow convection

    do j = 1,nterm
      i = index_nterm(j)

      flx_init_new(i) = flx_init(i)
      if (cldbase_opt_sh == sh_grey_closure) then
        !-----------------------------------------------------------------------
        ! So far the "shallow" closure has been used.
        ! If cloud depth is deep (>4km) then rescale with CAPE.
        ! Smoothly match between shallow and deep closure as cloud layer deepens
        ! from 1.5 to 4km.
        !-----------------------------------------------------------------------
        ! Calculate cloud layer depth relative to typical expected shallow cloud
        ! depth of 1.5km, so transition is from zscale between 0 amd 1
        !-----------------------------------------------------------------------
        z_scale = z_theta(i,k)-z_theta(i,ntml(i))
        z_scale = ( z_scale - 1500.0 )/(4000.0 - 1500.0)

        if (dcpbydt(i) > 0.0 .and. z_scale > 0.0 ) then

          !         ! Reset flx_init_new using deep closure but ensuring
          !         ! greater than the original shallow closure
          flx_init_new(i) = flx_init_new(i)*max( 1.0,                          &
                            cape(i)/(cape_timescale*dcpbydt(i)) )

          if (z_scale < 1.0) then
            !           ! Take linear sum of shallow (flux_init) and
            !           ! deep (flux_init_new) closures
                        ! [zscale varies between 0 and 1]
            flx_init_new(i)= z_scale * flx_init_new(i) +                       &
                       (1.0-z_scale) * flx_init(i)
          end if

        end if  ! dcpbydt > 0
        !-----------------------------------------------------------
        ! Calculate Honnert et al style subgrid weighting as a
        ! function of the cloud-top height to grid size ratio
        !-----------------------------------------------------------
        weight_param(i) = 1.0 -                                                &
              tanh( beta_cu*z_theta(i,k)/delta_smag(i)) * max( 0.0,            &
           min( 1.0, (linear0-delta_smag(i)/z_theta(i,k))*rlinfac ) )
        ! Scale flx_init using parametrization weighting
        flx_init_new(i) = weight_param(i)*flx_init_new(i)
      end if

      max_cfl(i) = max_cfl(i) * timestep
      if (max_cfl(i)  >   cfl_limit) then
        flx_init_new(i) = flx_init_new(i) * cfl_limit / max_cfl(i)
      end if

      if (flx_init_new(i)  >   flxmax_init(i)) then
        flx_init_new(i) = flxmax_init(i)
      end if
      max_cfl(i) = 0.0
    end do      ! j (nterm)

    !
    ! Scale cloud fraction
    !
    do j = 1,nterm
      i = index_nterm(j)

      if (flx_init_new(i) > 0.0) then
        scale_f(i) = flx_init_new(i) / flx_init(i)
        cca_2d(i)  = cca_2d(i) + 0.06 * log(scale_f(i))

        ! set the flx_init to the new value to provide the real initial mass
        ! flux under all conditions
        flx_init(i) = flx_init_new(i)

      end if

      ! Check scaled cloud fraction not smaller than minimum value
      ! (2.0E-5) or greater than unity.
      !
      ! (Was moved out of scaling if test to ensure these limits
      ! at all times, not just when cca_2d is scaled)
      !
      cca_2d(i) = max(2.0e-5, cca_2d(i))
      cca_2d(i) = min(1.0e+0, cca_2d(i))

    end do      ! j (nterm)

    !-----------------------------------------------------------------------
    ! Check for false/true convection and reset variables appropriately
    !-----------------------------------------------------------------------

    ! First, save SCM diagnostic of whether convection failed, and if so, why:
    if ( l_scm_convss_dg ) then
      do j = 1,nterm
        i = index_nterm(j)

        ! Convection failed due to the ascent being too shallow or cloud-free
        if ( icct(i)-iccb(i) <= 3 .or. ( .not. blatent(i) ) ) then
          scm_convss_dg(i) % status_shallow = 1

          ! Convection failed because the closure set the mass-flux to zero
        else if ( flx_init_new(i)  <=  minflx ) then
          scm_convss_dg(i) % status_shallow = 2

          ! Real shallow convection occurred!
        else
          scm_convss_dg(i) % status_shallow = 3

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
        kterm(i)        = 0
      else
        ! True convection
        kterm(i)        = k
        ind_shallow(i)  = 1.0
      end if
    end do  ! nterm

    !-----------------------------------------------------------------------
    ! Apply cfl scaling
    !-----------------------------------------------------------------------
    do kt = 2, k+1
      do j = 1,nterm
        i = index_nterm(j)
        if (kt  >=  ntml(i)) then
          dthbydt(i,kt)   = dthbydt(i,kt)  * scale_f(i)
          dqbydt(i,kt)    = dqbydt(i,kt)   * scale_f(i)
          dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
          dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
          dcflbydt(i,kt)  = dcflbydt(i,kt) * scale_f(i)
          dcffbydt(i,kt)  = dcffbydt(i,kt) * scale_f(i)
          dbcfbydt(i,kt)  = dbcfbydt(i,kt) * scale_f(i)
          if (l_mom) then
            dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
            dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
          end if
          if (l_tracer) then
            do ktra = 1,ntra
              dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
            end do
          end if
          if (l_wtrac_conv) then
            do i_wt = 1, n_wtrac
              wtrac_e(i_wt)%dqbydt(i,kt)   =                                   &
                         wtrac_e(i_wt)%dqbydt(i,kt)   * scale_f(i)
              wtrac_e(i_wt)%dqclbydt(i,kt) =                                   &
                         wtrac_e(i_wt)%dqclbydt(i,kt) * scale_f(i)
              wtrac_e(i_wt)%dqcfbydt(i,kt) =                                   &
                         wtrac_e(i_wt)%dqcfbydt(i,kt) * scale_f(i)
              wtrac_p(i_wt)%precip(i,kt) =                                     &
                         wtrac_p(i_wt)%precip(i,kt)   * scale_f(i)
            end do
          end if

          flx(i,kt)    = flx(i,kt)    *  scale_f(i)
          precip(i,kt) = precip(i,kt) *  scale_f(i)

        end if !kt >ntml and flx_init_new >0
      end do  ! j loop
    end do  ! kt loop

    do j = 1,nterm
      i = index_nterm(j)
    end do  ! nterm loop

    !-----------------------------------------------------------------------
    ! 3.5  Downdraught calculation - on all points where convection is
    !      terminating. Downdraughts are possible for some deeper shallow
    !      convection.
    !
    !      Subroutine DD_CALL
    !
    !      UM Documentation Paper 27, part 2
    !
    !-----------------------------------------------------------------------

    npossdd=0
    nnodd = 0

    do i = 1,nterm
      i2=index_nterm(i)
      tempnum=0.0
      if (iccb(i2) >  0) then
        deltap_cld=p_layer_centres(i2,iccb(i2)) -p_layer_centres(i2,k)
        do kt=iccb(i2),k+1
          tempnum=tempnum+precip(i2,kt)
        end do
      else
        deltap_cld = 0.0
      end if

      ! Set logical to determine if downdraughts are allowed or not.
      if (deltap_cld >  15000.0 .and. bgmk(i2) .and. tempnum >  1e-12) then
        b_dd(i2) = .true.
      else
        b_nodd(i2) = .true.
      end if
    end do  ! nterm loop

    ! If convection has terminated write cape to diagnostic output
    ! variable (cape_out).

    do j = 1,nterm
      i=index_nterm(j)
      cape_out(i) = cape(i)
      dcpbydt(i) = 0.0
      cape(i) = 0.0
      bconv(i) = .false.
      det_lev(i)= k+1 ! Set final detrainment level (but not used).
    end do

  end if  ! nterm > 0

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
                * (1.0 - amdetk(j)) * (ekp14(j) + ekp34(j)                     &
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

if (l_wtrac_conv) then
  deallocate(ratio_wtrac)
  deallocate(qpixs_v_wtrac)
end if

!-----------------------------------------------------------------------
! Time smoothed mass flux
!-----------------------------------------------------------------------
if (l_conv_prog_flx) then
  decay_amount = timestep/tau_conv_prog_flx
  ! Update the smoothed mass flux
  do k = 1, nlev
    do i = 1,n_sh
      conv_prog_flx(i,k) = conv_prog_flx(i,k) + decay_amount * flx(i,k)
    end do
  end do
end if


! If used, copy convection profile diagnostics
if (l_scm_convss_dg) then

  ! Copy profile diagnostics
  do k = 1, nlev
    do i = 1, n_sh
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
! Write out updraught massflux diagnostics and scale the
! entrainment and detrainment diagnostics by the mass flux.
!-----------------------------------------------------------------------
if (flg_up_flx .or. flg_mf_shall) then
  do k = 1,nlev
    do i = 1,n_sh
      up_flux(i,k) = flx(i,k)
    end do
  end do
end if

if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,n_sh
      up_flux_half(i,k) = up_flux_half(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_entr_up) then
  do k = 1,nlev
    do i = 1,n_sh
      entrain_up(i,k) = entrain_up(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,n_sh
      detrain_up(i,k) = detrain_up(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_area_ud) then
  do k = 1,nlev
    do i = 1,n_sh
      if (flx(i,k) > 0.0) then
        ! expect to get a valid wup
        ! updraught core area from wup and flx
        if (w2p(i,k) > 0.0) then
          ! This seems to give cloud areas where the area is smaller higher
          ! in the convective core
          area_ud(i,k) =  flx(i,k)/(g*sqrt(w2p(i,k))* rho_theta(i,k))
          if (area_ud(i,k) >= 1.0) then
            ! Something very wrong try
            area_ud(i,k) = flx(i,k)/(g* rho_theta(i,k))
          end if
        else
          ! The calculation of w2p is not producing a sensible value
          ! for all shallow ascents.
          ! What do I want to do here? Assumption that at cloud base
          ! wup ~ 1.0m/s use this as may give something useable?
          ! Taking this approach gives cloud top areas greater than
          ! cloud base areas which causes problems with my DD area
          ! parametrization. So may be better to only do this at cloud base.
          if (k == iccb(i)+1) then
            area_ud(i,k) =  flx(i,k)/(g* rho_theta(i,k))
          end if
        end if
      else
        area_ud(i,k)    = 0.0
      end if
    end do
  end do
end if
!-----------------------------------------------------------------------
! 4.0  Mixes the increments from the initial parcel perturbation throughout
!      the subcloud layer
!-----------------------------------------------------------------------

call mix_ipert_6a(n_sh, nlev, nbl, n_wtrac, ntml,                              &
                  p_layer_boundaries, exner_layer_centres,                     &
                  dthbydt, dqbydt, wtrac_e, flx_init,                          &
                  thpert, qpert, wtrac_p)

!-----------------------------------------------------------------------
! 5.0 All shallow convection will terminate at some level. This level
!     has been stored in the main level loop.
!     The convection will either have a down draught or none will be
!     possible.
!     Downdraughts & evaporation - now 3 options
!                     (a) Emanuel downdraught scheme
!                     (b) Original mass flux code
!                     (c) New mass flux scheme
!-----------------------------------------------------------------------

if (l_eman_dd .and. l_snow_rain) then
  ! (a) Use revised Emanuel downdraughts for DD and evaporation below cloud
  ! base

  ! Work out maximum termination level
  kmax_term = 2
  do i = 1,n_sh
    if (kterm(i) >  kmax_term) then
      kmax_term = kterm(i)
    end if
  end do

  call eman_dd_rev (n_sh,kmax_term,nlev,trlev,ntra,                            &
                  kterm,l_tracer,                                              &
                  exner_layer_centres,                                         &
                  p_layer_centres, p_layer_boundaries,                         &
                  scale_f, th, q, qse, tracer, precip,                         &
                  dthbydt, dqbydt, dtrabydt,                                   &
                  rain, snow ,dwn_flux, dt_dd, dq_dd)

else ! Standard downdraughts and evaporation below cloud base

  !-----------------------------------------------------------------------
  ! (b) Original downdraught calculation - on all points where convection is
  !      terminating.
  !
  !      Subroutine DD_ALL_CALL
  !
  !      UM Documentation Paper 27
  !
  !-----------------------------------------------------------------------

  npossdd = 0
  do i = 1,n_sh
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

    call dd_all_call_6a(n_sh,npossdd,kmax_term,nlev,trlev,ntra,n_wtrac         &
                  , kterm,iccb,icct,index_possdd,l_tracer                      &
                  , bwater(1,2),exner_layer_centres                            &
                  , exner_layer_boundaries,p_layer_centres                     &
                  , p_layer_boundaries,pstar,recip_pstar,timestep              &
                  , cca_2d,thp,qp,th,q,qse,trap,tracer,flx,precip              &
                  , dthbydt,dqbydt,dtrabydt,rain,snow,rain_3d,snow_3d          &
                  , wtrac_p,wtrac_e                                            &
                  , dwn_flux,entrain_dwn,detrain_dwn,dt_dd,dq_dd)
  end if

  !-----------------------------------------------------------------------
  ! 5.2 Surface precipitation calculation for terminating points with
  !     no downdraught (moved outside level loop) ie do this calculation
  !     on all points at the end.
  !-----------------------------------------------------------------------
  ! Points where no downdraught possible
  nnodd = 0
  do i = 1,n_sh
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

    call evap_bcb_nodd_all_6a(n_sh,nnodd,n_wtrac,kmax_term,kterm               &
,                        iccb, index_nodd, bwater(1,2)                         &
,                        exner_layer_centres                                   &
,                        p_layer_centres, p_layer_boundaries                   &
,                        timestep , cca_2d, th, q, qse, precip                 &
,                        dthbydt, dqbydt                                       &
,                        rain, snow, rain_3d, snow_3d                          &
,                        dt_dd, dq_dd, wtrac_p, wtrac_e )

  end if

end if ! not Emanuel DD

! Deallocate water tracer parcel fields
call wtrac_dealloc_conv_p(n_wtrac, wtrac_p)

!-----------------------------------------------------------------------
! 6.0  Convective Momentum Transport (if L_mom = .T.)
!-----------------------------------------------------------------------

if (l_mom) then

  if (cldbase_opt_sh == sh_grey_closure) then
    !   ! use convection parcel top, rather than adiabatic
    do i = 1,n_sh
      ntop_uv(i) = kterm(i) + 1
      ptop_uv(i) = phalf_uv(ntop_uv(i),i)
      zcld(i)    = z_rho(i,ntop_uv(i)) - zlcl_uv(i)
      if ( (kprof_cu == klcl_entr .and. zcld(i) < max_cu_depth) .or.           &
           zcld(i) < 0.1* zlcl_uv(i) ) then
        ! Require non-zero cloud depth so test on top of BL mixing
        ! (max_cu_depth, or 1.1*zlcl as a sensible default)
        ind_shallow(i) = 0.0
      end if
    end do
  end if

  nterm = 0
  do i = 1, n_sh
    ! Check that cloud base mass flux is non-zero.
    ! If convection failed then mb=0
    ! If cldbase_opt_sh=sh_grey_closure then cloud depth could have been
    ! zero so test on ind_shallow too (which is a real equal to 0 or 1
    ! so test on 0.5)
    if ( mb(i) > 0.0 .and. (ind_shallow(i) > 0.5 .or.                          &
                    cldbase_opt_sh == sh_wstar_closure) ) then
      nterm = nterm + 1
      cu_term(nterm) = i
      cu_tend(nterm) = i
    end if
  end do

  if (nterm  >   0) then

    call shallow_grad_stress(n_sh,n_sh,nterm,nlev,cu_term,                     &
                             nlcl_uv,ntop_uv,mb,wsc,wstar,zcld,                &
                             plcl_uv,ptop_uv,p_uv,phalf_uv,                    &
                             rho_uv,ue_p,ve_p,timestep,                        &
                             weight_param,                                     &
                             ! in
                             uw,vw)

    call shallow_base_stress(n_sh,n_sh,n_sh,nlev,nterm,cu_term,                &
                             cu_tend,nlcl_uv,ntop_uv,mb,wsc,                   &
                             zlcl_uv,uw0,vw0,plcl_uv,                          &
                             ptop_uv,ue_p,ve_p,phalf_uv,p_uv,                  &
                             rho_uv,timestep,weight_param,                     &
                             ! INOUT
                             uw,vw,                                            &
                             ! out
                             uw_shall,vw_shall)

    call shallow_cmt_incr(n_sh,n_sh,n_sh,nlev,nterm,cu_term,                   &
                          cu_tend,nlcl_uv,ntop_uv,uw,vw,phalf_uv,              &
                          rho_uv,zlcl_uv,                                      &
                          !out
                          dubydt,dvbydt)

  end if  ! nterm>0
end if ! L_mom

!-----------------------------------------------------------------------
! 6.1  Add the dissipative heating from the CMT to the theta increment
!-----------------------------------------------------------------------
if (l_mom .and. l_cmt_heating) then
  call cmt_heating(n_sh, nlev,                                                 &
                   z_theta, z_rho, exner_layer_centres,                        &
                   u, v, dubydt, dvbydt,                                       &
                   ! Out
                   dthbydt)
end if

!-----------------------------------------------------------------------
! 7.0  Energy (and optionally water) correction calculation
!-----------------------------------------------------------------------
do i = 1,n_sh
  index1(i) = i
end do

if (l_cv_conserve_check) then
  call cor_engy_6a(n_sh, n_sh, nlev, n_wtrac, index1, timestep,                &
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
! 8.0  Correct negative/very small humidities
!-----------------------------------------------------------------------

call correct_small_q_conv(n_sh, n_sh, nlev, n_wtrac, index1, timestep, qmin,   &
                          r2rho_th, dr_across_th,  q, 'water', 'Shallow',      &
                          dqbydt, wtrac_e=wtrac_e)

! Then use same method to correct any remaining negative/v small values of
! water tracer vapour

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    call correct_small_q_conv(n_sh, n_sh, nlev, n_wtrac, index1, timestep,     &
                              wtrac_info(i_wt)%qlimit, r2rho_th, dr_across_th, &
                              wtrac_e(i_wt)%q, 'wtrac', 'Shallow',             &
                              wtrac_e(i_wt)%dqbydt)
  end do
end if

!-------------------------------------------------------------------------
! 9.1 CCRad - Calculate CCA fow shallow levels only
!-------------------------------------------------------------------------

do i=1, n_sh

  overlap_fac(i) = 0.0

  if (iccb(i) /= 0) then ! Shallow convection occured

    overlap_fac(i) = 1.0

    if (cca2d_sh_opt == grant_lock_over) then
      !---------------------------------------------------------------
      ! Grant and Lock (2004) LES show mb/wsc nicely scales the cloud
      ! fraction profiles but not the TCA.  Also the UM overlap
      ! assumption in radiation is maximal.  This implies significant
      ! underestimate of TCA. So include a further parametrization of
      ! Cu "overlap", based again on the LES of Grant and Lock (2004).
      ! This increases cca_2d proportional to the ratio of the cloud
      ! to sub-cloud layer depths.  In order to preserve the grid-box
      ! cloud water, ccw will be divided by the same factor.
      !---------------------------------------------------------------

      overlap_fac(i) = 2.0 * (  z_rho(i,ntpar(i)+1)                            &
                              - z_rho(i, ntml(i)+1) )                          &
                     / z_rho(i,ntml(i)+1)
    end if     ! cca2d_sh_opt


  end if     ! iccb
end do     ! n_sh

!---------------------------------------------------------------
! 9.11 Calculate CCA
!---------------------------------------------------------------
select case (cca2d_sh_opt)
case (grant_lock, grant_lock_over)

  do i=1, n_sh
    if (iccb(i) /= 0) then ! Shallow convection occured
      tempnum = 2.0*mb(i)/wsc(i)

      cca_2d(i) = max(2.0e-5, tempnum)
      !------------------------------------------------------
      ! Include grey-zone weighting
      !------------------------------------------------------
      if (cldbase_opt_sh == sh_grey_closure)                                   &
            cca_2d(i) = cca_2d(i)*weight_param(i)

      cca_2d(i) = min( 1.0, cca_2d(i))

      ! Will be used by name, grab lowest cca_2d before any
      ! Tuning knobs applied
      lcca(i) = cca_2d(i)
    end if     ! iccb
  end do     ! n_sh


case (total_condensed_water)
    ! cca_2d is left unchanged from that calculated in the
    ! code, which is based on TCW (Total Condensed Water)
    ! (TCW is a rate)

end select


if (cca2d_sh_opt == grant_lock_over) then

  do i=1, n_sh

    overlap_fac(i) = max( 0.5, overlap_fac(i) )
    overlap_fac(i) = min( 5.0, overlap_fac(i) )

    if (overlap_fac(i)*cca_2d(i) > 0.99) then
      overlap_fac(i) = 0.99/cca_2d(i)
    end if
  end do      ! i (n_sh)

end if ! cca2d_sh_opt


!-------------------------------------------------------------------
! 9.12 Fill cca with cca_2d where non-zero ccw
!-------------------------------------------------------------------
do k=1, nlev
  do i=1, n_sh
    if (iccb(i) /= 0) then ! Shallow convection occured

      if (ccw(i,k) > 0.0) then

        zpr = (z_rho(i,k)          - z_rho(i,ntml(i)+1))                       &
            / (z_rho(i,ntpar(i)+1) - z_rho(i,ntml(i)+1))

        ! Apply Shape-function
        !
        ! Apply overlap_fac to cca, also preserving grid-box water
        ! by dividing ccw by overlap_fac, at least at cloud-base

        ccw(i,k)  = ccw(i,k)/overlap_fac(i)
        zpr       = min(1.0,zpr)
        cca(i,k)  = overlap_fac(i)*cca_2d(i)                                   &
                       * 0.25*( 1.0 + 3.0*exp(-5.0*zpr) )

      end if       ! ccw
    end if       ! iccb
  end do       ! i (n_sh)
end do       ! k (nlev)


!-----------------------------------------------------------------------
! Final SCM convection sub-step diagnostics
!-----------------------------------------------------------------------
if ( l_scm_convss_dg ) then

  ! 3-D diagnostics
  do k = 1, nlev
    do i = 1, n_sh
      if ( k>=ntml(i) .and. k<=det_lev(i) ) then
        ! Save the final updraft mass-flux
        scm_convss_dg(i) % up_flx(k) = flx(i,k)
      end if
    end do
  end do

  ! 2-D diagnostics
  do i = 1, n_sh
    scm_convss_dg(i) % precip_shallow = rain(i) + snow(i)
  end do

end if


!-----------------------------------------------------------------------
! 10.0  End Subroutine
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine shallow_conv_6a
end module shallow_conv_6a_mod
