! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Congestus convection scheme

module congest_conv_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Congetus convection scheme
!   works on points diagnosed as congestus in subroutine CONV_DIAG.
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


character(len=*), parameter, private :: ModuleName='CONGEST_CONV_6A_MOD'

contains

subroutine congest_conv_6a(nbl,nlev,ntra,n_cca_lev,n_cg,trlev,                 &
                       bland,delthvu,                                          &
                       exner_rho,                                              &
                       exner_layer_centres,                                    &
                       exner_layer_boundaries,                                 &
                       l_q_interact,                                           &
                       l_tracer,ntml,ntpar,                                    &
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
                       wstar,wthvs,entrain_coef,                               &
                       zlcl_uv,tnuc_nlcl,ztop_uv,freeze_lev,                   &
                       recip_pstar,qse,                                        &
                       ccp_strength,                                           &

                       ! InOut
                       bulk_cf,cf_frozen,cf_liquid,qcf,                        &
                       qcl,tracer,w2p,conv_prog_flx,                           &

                       ! Out
                       cape_out,cclwp,ccw,cca,                                 &
                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,             &
                       dqclbydt,dthbydt,                                       &
                       dubydt,dvbydt,dtrabydt,                                 &
                       detrain_up,detrain_dwn,                                 &
                       entrain_up,entrain_dwn,                                 &
                       iccb,icct,lcca,                                         &
                       lcbase,lctop,rain,snow,                                 &
                       rain_3d,snow_3d, up_flux, up_flux_half,                 &
                       dwn_flux,uw_shall,vw_shall,kterm,                       &
                       tcw,cca_2d, dt_dd, dq_dd, area_ud )

use planet_constants_mod, only:                                                &
    r, pref, kappa, repsilon, c_virtual, g

use cv_run_mod, only:                                                          &
    l_mom,                                                                     &
    cca2d_sh_opt,                                                              &
    l_cv_conserve_check,                                                       &
    l_cmt_heating, l_conv_prog_flx, tau_conv_prog_flx

use cv_param_mod, only:                                                        &
    total_condensed_water, grant_lock,                                         &
    thpixs_shallow, qpixs_shallow, c_mass,                                     &
    term_undil

use cv_dependent_switch_mod, only:                                             &
    cg_on, cg_new_termc

use cv_stash_flg_mod, only:                                                    &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                                  &
    flg_entr_up, flg_entr_dwn,                                                 &
    flg_detr_up, flg_detr_dwn, flg_uw_shall, flg_vw_shall,                     &
    flg_w_eqn, flg_area_ud

use water_constants_mod, only: lc, lf, tm

use yomhook,    only: lhook, dr_hook
use parkind1,   only: jprb, jpim

! subroutines
use lift_par_6a_mod,       only: lift_par_6a
use lift_undil_par_6a_mod, only: lift_undil_par_6a
use convec2_6a_mod,        only: convec2_6a
use water_loading_mod,     only: water_loading
use cor_engy_6a_mod,       only: cor_engy_6a
use mix_ipert_6a_mod,      only: mix_ipert_6a
use cmt_heating_mod,       only: cmt_heating
use layer_cn_6a_mod,       only: layer_cn_6a, congestus

use dd_all_call_6a_mod,       only: dd_all_call_6a
use evap_bcb_nodd_all_6a_mod, only: evap_bcb_nodd_all_6a
use flag_wet_mod,             only: flag_wet
use shallow_base_stress_mod,  only: shallow_base_stress
use shallow_cmt_incr_mod,     only: shallow_cmt_incr
use shallow_grad_stress_mod,  only: shallow_grad_stress

use qsat_mod, only: qsat
use correct_small_q_conv_mod, only: correct_small_q_conv

use wtrac_conv_mod, only: conv_e_wtrac_type, conv_p_wtrac_type,                &
                          wtrac_alloc_conv_e, wtrac_dealloc_conv_e,            &
                          wtrac_alloc_conv_p, wtrac_dealloc_conv_p


implicit none

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent in:

integer, intent(in) ::                                                         &
  nbl                  & ! No. of boundary layer levels
 ,nlev                 & ! No. of model layers
 ,ntra                 & ! No. of tracer fields
 ,n_cca_lev            & ! No. of convective cloud amount levels (1 for 2D,
                         ! nlevs for 3D)
 ,n_cg                 & ! No. of congestus convection points
 ,trlev                  ! No. of model levels on which tracers are included

logical, intent(in) :: bland(n_cg) ! Land/sea mask

real(kind=real_umphys), intent(in)    :: delthvu(n_cg)
                                     !Integral of undilute parcel
                                     ! buoyancy over convective cloud
                                     ! layer (Kelvin m)

real(kind=real_umphys), intent(in)    :: exner_rho(n_cg,nlev)
                                             ! Exner on rho levels

real(kind=real_umphys), intent(in)    ::                                       &
  exner_layer_centres(n_cg,0:nlev)    & ! Exner
 ,exner_layer_boundaries(n_cg,0:nlev)   ! Exner at half level above
                                        ! exner_layer_centres

logical, intent(in) ::                                                         &
  l_q_interact         & ! Switch allows overwriting parcel variables
                         ! when calculating condensate incr.
 ,l_tracer               ! Switch for inclusion of tracers

integer, intent(in) ::                                                         &
  ntml(n_cg)           & ! Top level of surface mixed layer defined relative to
                         ! theta,q grid
 ,ntpar(n_cg)            ! Top level of initial parcel ascent in BL scheme
                         ! defined relative to theta, q grid

real(kind=real_umphys), intent(in)    ::                                       &
  pstar(n_cg)                     & ! Surface pressure (Pa)
 ,p_layer_centres(n_cg,0:nlev)    & ! Pressure (Pa)
 ,p_layer_boundaries(n_cg,0:nlev)   ! Pressure at half level above
                                    ! p_layer_centres (Pa)

! Note heights passed in but not currently used - will be
! required by new turbulence based scheme therefore been added to
! arguement list ready for future developments.

real(kind=real_umphys), intent(in) :: z_theta(n_cg,nlev)
                                            ! height of theta levels (m)
real(kind=real_umphys), intent(in) :: z_rho(n_cg,nlev)
                                            ! height of rho levels (m)
real(kind=real_umphys), intent(in) :: r_theta(n_cg,0:nlev)
                                            ! radius of theta levels (m)
real(kind=real_umphys), intent(in) :: r_rho(n_cg,nlev)
                                            ! radius of rho levels (m)
real(kind=real_umphys), intent(in) :: rho_theta(n_cg,nlev)
                                            ! density for theta lev (kg/m3)
real(kind=real_umphys), intent(in) ::                                          &
  rho_dry_theta(n_cg,nlev) & ! dry density on theta levels (kg/m3)
 ,rho_dry(n_cg,nlev)         ! dry density on rho levels (kg/m3)

real(kind=real_umphys), intent(in) :: r2rho_th(n_cg,nlev)
                                            ! radius**2 density for
                                            ! theta lev (kg/m)
real(kind=real_umphys), intent(in) :: r2rho(n_cg,nlev)
                                            ! radius**2 density for
                                            ! rho lev (kg/m)
real(kind=real_umphys), intent(in) :: dr_across_th(n_cg,nlev)
                                            ! thickness of theta levels (m)
real(kind=real_umphys), intent(in) :: dr_across_rh(n_cg,nlev)
                                            ! thickness of rho levels (m)
real(kind=real_umphys), intent(in) :: conv_prog_precip(n_cg,nlev)
                                                ! Surface precipitation based
                                                ! 3d convective prognostic in
                                                ! kg/m2/s
real(kind=real_umphys), intent(in)  ::                                         &
  q(n_cg,nlev)        & ! Model mixing ratio (kg/kg)
 ,th(n_cg,nlev)       & ! Model potential temperature (K)
 ,timestep            & ! Model timestep (s)
 ,u(n_cg,nlev)        & ! Model U field (m/s)
 ,v(n_cg,nlev)        & ! Model V field (m/s)
 ,uw0(n_cg)           & ! U-comp of surface stress (N/m2)
 ,vw0(n_cg)           & ! V-comp of surface stress (N/m2)
 ,wstar(n_cg)         & ! Convective velocity scale (m/s)
 ,wthvs(n_cg)         & ! Surface flux of THV (Pa m/s2)
 ,entrain_coef(n_cg)  & ! entrianment coefficients
 ,zlcl_uv(n_cg)       & ! Lifting condensation level defined for the uv grid (m)
 ,tnuc_nlcl(n_cg)     & ! nucleation temperature as function of dust
                        ! indexed using nlcl(deg cel)
 ,ztop_uv(n_cg)         ! Top of cloud layer defined for the uv grid (m)

integer, intent(in) :: freeze_lev(n_cg) ! Level index for freezing level

real(kind=real_umphys), intent(in) :: recip_pstar(n_cg)
                                       ! Reciprocal of pstar array

real(kind=real_umphys), intent(in) :: qse(n_cg,nlev)
                                   ! Saturation mixing ratio of
                                   ! cloud environment (kg/kg)

real(kind=real_umphys), intent(in)    ::                                       &
  ccp_strength(n_cg)  ! Cold-pool strength

! Arguments with intent INOUT:

real(kind=real_umphys), intent(in out) ::                                      &
  bulk_cf(n_cg,nlev)     & ! Bulk total cloud volume ( )
 ,cf_frozen(n_cg,nlev)   & ! Frozen water cloud volume ( )
 ,cf_liquid(n_cg,nlev)   & ! Liq water cloud volume ( )
 ,qcf(n_cg,nlev)         & ! Ice condensate mix ratio (kg/kg)
 ,qcl(n_cg,nlev)         & ! Liq condensate mix ratio (kg/kg)
 ,tracer(n_cg,trlev,ntra)  ! Model tracer fields(kg/kg)

real(kind=real_umphys), intent(in out) ::                                      &
  w2p(n_cg,nlev)           ! (Parcel vertical velocity)^2 [(m/s)^2]

real(kind=real_umphys), intent(in out) :: conv_prog_flx(n_cg,nlev)
                                                 ! Mass flux convective
                                                ! prognostic in Pa/s

! Arguments with intent out:

real(kind=real_umphys), intent(out) ::                                         &
  cape_out(n_cg)     & ! Saved convective available potential energy for
                       ! diagnostic output (J/kg)
 ,cclwp(n_cg)        & ! Condensed water path (k/m2)
 ,ccw(n_cg,nlev)     & ! Convective cloud liquid water on model levels (g/kg)
 ,cca(n_cg,n_cca_lev)  ! Convective cloud amount on model levels (fraction)

real(kind=real_umphys), intent(out) ::                                         &
  dbcfbydt(n_cg,nlev) & ! Increments to total cld volume due to convection(/s)
 ,dcffbydt(n_cg,nlev) & ! Increments to ice cloud volume due to convection (/s)
 ,dcflbydt(n_cg,nlev) & ! Increments to liq cloud volume due to convection (/s)
 ,dqbydt(n_cg,nlev)   & ! Increments to q due to convection (kg/kg/s)
 ,dqcfbydt(n_cg,nlev) & ! Increments to ice condensate due to conv (kg/kg/s)
 ,dqclbydt(n_cg,nlev) & ! Increments to liq condensate due to conv (kg/kg/s)
 ,dthbydt(n_cg,nlev)  & ! Increments to potential temp. due to convection (K/s)
 ,dubydt(n_cg,nlev)   & ! Increments to U due to CMT (m/s2)
 ,dvbydt(n_cg,nlev)     ! Increments to V due to CMT (m/s2)

real(kind=real_umphys), intent(out) ::                                         &
  dtrabydt(n_cg,nlev,ntra)   !Increment to tracer due to convection (kg/kg/s)

real(kind=real_umphys), intent(out) ::                                         &
  detrain_up(n_cg,nlev)  & ! Fractional detrainment rate into updraughts (Pa/s)
 ,detrain_dwn(n_cg,nlev) & ! Fractional detrainment rate into
                           ! downdraughts (Pa/s)
 ,entrain_up(n_cg,nlev)  & ! Fractional entrainment rate into updraughts (Pa/s)
 ,entrain_dwn(n_cg,nlev)   ! Fractional entrainment rate into
                           ! downdraughts (Pa/s)

integer, intent(out) :: iccb(n_cg) ! Convective cloud base level

integer, intent(out) :: icct(n_cg) ! Convective cloud top level

real(kind=real_umphys), intent(out) :: lcca(n_cg) ! Lowest conv. cloud amt. (%)


integer, intent(out) :: lcbase(n_cg) ! Lowest conv. cloud base level

integer, intent(out) :: lctop(n_cg) ! Lowest conv. cloud top level

real(kind=real_umphys), intent(out) ::                                         &
  rain(n_cg)           & ! Surface convective rainfall (kg/m2/s)
 ,snow(n_cg)           & ! Surface convective snowfall (kg/m2/s)
 ,rain_3d(n_cg,nlev)   & ! rainfall flux (kg/m2/s)
 ,snow_3d(n_cg,nlev)     ! snowfall flux (kg/m2/s)

real(kind=real_umphys), intent(out) ::                                         &
  up_flux(n_cg,nlev)   & ! Updraught mass flux (Pa/s)
 ,up_flux_half(n_cg,nlev)& ! Updraught mass flux on half levels (Pa/s)
 ,dwn_flux(n_cg,nlev)  & ! Downdraught mass flux (Pa/s)
 ,uw_shall(n_cg,nlev)  & ! X-comp. of stress from shallow convection (kg/m/s2)
 ,vw_shall(n_cg,nlev)    ! Y-comp. of stress from shallow convection (kg/m/s2)

integer, intent(out) :: kterm(n_cg) ! termination level

real(kind=real_umphys), intent(out) ::                                         &
  tcw(n_cg)            & ! Total condensed water(kg/m2/s)
 ,cca_2d(n_cg)           ! 2D convective cloud amount (%)

! Downdraught and evap below cloud base
real(kind=real_umphys), intent(out) ::                                         &
  dt_dd(n_cg,nlev)        & ! dT/dt from DD and evap below cloud base (K/s)
 ,dq_dd(n_cg,nlev)        & ! dq/dt from DD and evap below cloud base (kg/kg/s)
 ,area_ud(n_cg,nlev)        ! fractional area of updraughts
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

! Indices of work points.
! This routine will use a special case of indirect indexing where all
! the point are contiguous.
integer :: idx_contig(n_cg)

! Height above surface of model levels ...
real(kind=real_umphys) :: zkm1(n_cg)              ! ... k-1   [m]
real(kind=real_umphys) :: zk(n_cg)                ! ... k     [m]
real(kind=real_umphys) :: zkp12(n_cg)             ! ... k+1/2 [m]
real(kind=real_umphys) :: zkp1(n_cg)              ! ... k+1   [m]

integer :: index1(n_cg),index2(n_cg)

integer :: ncposs               ! No. of points which may convect

integer :: nconv                ! No. of convecting points

real(kind=real_umphys) :: amdetk(n_cg)
                                ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

real(kind=real_umphys) :: b_calc                  ! Coefficient in thpert calc.

real(kind=real_umphys) :: c_calc                  ! Coefficient in thpert calc.

real(kind=real_umphys) :: cape(n_cg)
                                ! Convective available potential
                                ! energy (J/kg)
real(kind=real_umphys) :: fcape(n_cg)
                                ! Convective available potential
                                ! energy weighted by f.det profile (J/kg)

real(kind=real_umphys) :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

real(kind=real_umphys) :: dcpbydt(n_cg)
                                ! Rate of change of cape (J/kg/s)

real(kind=real_umphys) :: depth(n_cg)
                                ! Depth of convective cloud (m)

real(kind=real_umphys) ::                                                      &
  ekp14(n_cg)  & ! Entrainment coefficients at level k+1/4 multiplied by
                 ! appropriate layer thickness (dimensionless)
 ,ekp34(n_cg)    ! Entrainment coefficients at level k+3/4 multiplied by
                 ! appropriate layer thickness (dimensionless)

real(kind=real_umphys) ::                                                      &
  exk(n_cg)          & ! Exner ratio at layer k
 ,exkp1(n_cg)        & ! Exner ratio at layer k+1
 ,flxmax(n_cg)       & ! Maximum initial convective mass flux (Pa/s)
 ,flx_init(n_cg)     & ! Initial mass flux at cloud base (Pa/s)
 ,flx_init_new(n_cg) & ! flx_init scaled (Pa/s)
 ,flxmax_init(n_cg)  & ! Maximum possible initial mass flux (limited to the
                       ! mass in the initial convecting layer in Pa/s)
 ,max_cfl(n_cg)      & ! Max cfl ratio over a convecting layer
 ,p_lcl(n_cg)        & ! Pressure at LCL (Pa)
 ,precip(n_cg,nlev)    ! Amount of precip from each layer (kg/m/s)

real(kind=real_umphys) :: decay_amount
                       ! decay fraction for time-smoothed mass flux

real(kind=real_umphys) ::                                                      &
  pk(n_cg)           & ! Pressure at midpoint of layer k (Pa)
 ,pkp1(n_cg)         & ! Pressure at midpoint of layer k+1 (Pa)
 ,delpk(n_cg)        & ! Pressure difference over layer k (Pa)
 ,delpkp1(n_cg)      & ! Pressure difference over layer k+1 (Pa)
 ,delpkp12(n_cg)     & ! Pressure difference between layers k and k+1 (Pa)
 ,delp_uv_k(n_cg)    & ! Pressure difference across uv layer k (Pa)
 ,delp_uv_kp1(n_cg)    ! Pressure difference across uv layer k+1 (Pa)

real(kind=real_umphys) ::                                                      &
  q_lcl(n_cg)        & ! Mixing ratio at LCL (kg/kg)
 ,qse_lcl(n_cg)      & ! Saturated q at LCL (kg/kg)
 ,rhum(n_cg)         & ! Dummy relative humidity (only used on shallow points)
 ,t_lcl(n_cg)        & ! Temperature at LCL (K)
 ,th_lcl(n_cg)       & ! Theta at LCL (K)
 ,thv_pert(n_cg)     & ! Theta_v parcel pertubation (K)
 ,thpert(n_cg)       & ! Theta parcel pertubation (K)
 ,qpert(n_cg)          ! q parcel pertubation (kg/kg)

integer :: start_lev(n_cg)      ! Convection initiation level
integer :: start_lev_c2(n_cg)   ! Compressed convection initiation level

real(kind=real_umphys) :: wsc(n_cg)
                                ! Convective velocity scale (m/s)

real(kind=real_umphys) :: wsc_o_mb(n_cg)
                                ! Convective velocity scale /mb
real(kind=real_umphys) :: w_max(n_cg)
                                ! dummy variable for maximum w in column (m/s)


logical ::                                                                     &
  bgmk(n_cg)         & ! Mask for points where parcel in layer k is saturated
 ,bwater(n_cg,2:nlev)& ! Mask for points at which condensate is liquid
 ,blowst(n_cg)       & ! Dummy variable indicating low enough stability for
                       ! convection to occur
 ,bterm(n_cg)        & ! Mask for points which have stopped convecting
 ,bconv(n_cg)        & ! Mask for points at which convection is occurring
 ,bcposs(n_cg)         ! Mask for points passing initial stability test

logical :: blatent(n_cg)! Mask for points where latent heat has
                        ! been released

! Parcel variables

real(kind=real_umphys) :: qpi(n_cg)
                                ! Initial parcel mixing ratio (kg/kg)

real(kind=real_umphys) :: qp(n_cg,nlev)           ! Parcel mixing ratio (kg/kg)

real(kind=real_umphys) :: thpi(n_cg)
                                ! Initial parcel potential temp.(K)

real(kind=real_umphys) :: thp(n_cg,nlev)          ! Parcel potential temp (K)

real(kind=real_umphys) :: up(n_cg,nlev)           ! Parcel U (m/s)

real(kind=real_umphys) :: vp(n_cg,nlev)           ! Parcel V  (m/s)

real(kind=real_umphys) :: trap(n_cg,nlev,ntra)
                                ! Tracer content of parcel (kg/kg)

real(kind=real_umphys) :: expi(n_cg)
                                ! Initial parcel exner pressure

real(kind=real_umphys) :: flx(n_cg,nlev)          ! Parcel massflux (Pa/s)

real(kind=real_umphys) :: xsbmin_v(n_cg,nlev)
                                ! Minmum parcel buoyancy excess

real(kind=real_umphys) :: thpixs_v(n_cg,nlev)     ! Theta parcel excess (K)

real(kind=real_umphys) :: qpixs_v(n_cg,nlev)      ! Q parcel excess(kg/kg)

real(kind=real_umphys) :: qclp(n_cg,nlev)
                                ! Parcel liquid condensated mixing
                                ! ratio in layer k (kg/kg)

real(kind=real_umphys) :: qcfp(n_cg,nlev)
                                ! Parcel frozen condensated mixing
                                ! ratio in layer k (kg/kg)

! Undilute parcel variables
real(kind=real_umphys) :: qu(n_cg,nlev)
                                ! Undilute Parcel mixing ratio (kg/kg)
real(kind=real_umphys) :: thu(n_cg,nlev)
                                ! Undilute Parcel potential temp (K)

! Parameters

real(kind=real_umphys), parameter :: cfl_limit = 1.0 ! Max CFL ratio allowed
real(kind=real_umphys), parameter :: minflx = tiny(flx_init_new)
                                                ! minimum allowable
                                                ! initial mass flux

! CMT variables

integer ::                                                                     &
  nlcl_uv(n_cg)      & ! Level index for LCL
 ,ntop_uv(n_cg)      & ! Level index for top of layer
 ,n_0degc(n_cg)      & ! Level index for zero degrees
 ,cu_term(n_cg)      & ! Indicies for CMT subs
 ,cu_tend(n_cg)        !

real(kind=real_umphys) ::                                                      &
  exk_temp             & ! Temporary exner
 ,eflux_u_ud(n_cg)     & ! Vertical eddy flux of momentum due to UD at
                         ! top of layer (Pa m/s2)
 ,eflux_v_ud(n_cg)     & ! Vertical eddy flux of momentum due to UD at
                         ! bottom of layer (Pa m/s2)
 ,mb(n_cg)               ! Cloud base mass flux (Pa/s)

real(kind=real_umphys) :: p_uv(nlev,n_cg)         ! Pressure of model level (Pa)

real(kind=real_umphys) :: phalf_uv(nlev,n_cg)     ! Pressure of half level (Pa)

real(kind=real_umphys) :: plcl_uv(n_cg)           ! Pressure at LCL (Pa)

real(kind=real_umphys) :: ptop_uv(n_cg)
                                ! Pressure at top of cloud layer (Pa)

real(kind=real_umphys) :: p_0degc_uv(n_cg)
                                ! Pressure of zero degree level (Pa)

real(kind=real_umphys) :: rho_uv(nlev,n_cg)       ! Density on uv level (kg/m3)

real(kind=real_umphys) :: uw(nlev,n_cg)
                                ! U- comp stress profile (N/m2)
                                ! (units change through calls)

real(kind=real_umphys) :: ue_p(nlev,n_cg)         ! Environment U profile (m/s)

real(kind=real_umphys) :: vw(nlev,n_cg)           ! V-comp stress profile (N/m2)

real(kind=real_umphys) :: ve_p(nlev,n_cg)         ! Environment V profile (m/s)

real(kind=real_umphys) :: zcld(n_cg)              ! Depth of cloud layer (m)

logical :: l_mom_gk     ! Set to true if using Gregory-Kershaw CMT scheme
logical :: l_mom_gk_stable  ! Set to true if using stabilized
                            ! Gregory-Kershaw CMT scheme
                            ! (different from the original)

! CFL scaling variables

integer ::                                                                     &
  det_lev(n_cg)      & ! Level at which split final detrainment last occurred
 ,nterm              & ! No. of points where conv. has terminated
 ,index_nterm(n_cg)    ! Index for points where conv. has terminated

real(kind=real_umphys) ::                                                      &
  tempnum             & ! Temporary variable for storage
 ,scale_f(n_cg)       & ! store scaling factor
 ,weight_param(n_cg)    ! Weighting factor

real(kind=real_umphys) :: deltaktot(n_cg) ! Integrated forced detrainment

integer ::                                                                     &
  nnodd               & ! No. of downdraughts not possible
 ,index_nodd(n_cg)    & ! Index of downdraughts not possible
 ,npossdd             & ! No. downdraughts possible
 ,index_possdd(n_cg)  & ! Index of downdraughts
 ,kmax_term             ! maximum termination level + 1

real(kind=real_umphys) :: deltap_cld
                        ! pressure thickness of convective cloud (Pa)

! Limit nlev loop to those levels actually required using ntpar
! diagnosed in conv_diag

integer :: ntpar_max           ! max ntpar value

! parameters etc for qmin checks

real(kind=real_umphys), parameter :: qmin = 1.0e-8 ! Global minimum allowed Q

! Local compressed arrays

logical ::                                                                     &
  bgmkp1_c(n_cg)  & ! Mask for points where parcel in layer k+1 is saturated
 ,bgmkp1_c2(n_cg) & ! Mask for points where parcel in layer k+1 is saturated
 ,bwk_c(n_cg)     & ! bwater mask in layer k
 ,bwk_c2(n_cg)    & ! bwater mask in layer k
 ,bwkp1_c(n_cg)   & ! bwater mask in layer k+1
 ,bwkp1_c2(n_cg)    ! bwater mask in layer k+1

logical :: blatent_c2(n_cg)     ! Mask for points where latent heat has
                                ! been released

real(kind=real_umphys) :: thrk_c2(n_cg)
                                ! potential temperature of forced detrained air

real(kind=real_umphys) :: qrk_c2(n_cg)
                                ! specific humidity of forced detrained air

real(kind=real_umphys) :: deltak_c2(n_cg)
                                ! Parcel forced detrainment rate
                                ! in layer k multiplied by
                                ! appropriate layer thickness

real(kind=real_umphys) ::                                                      &
  dqek_c2(n_cg)    & ! Increment to q due to convection in layer k (kg/kg)
 ,dqekp1_c2(n_cg)  & ! Increment to q due to convection in layer k+1 (kg/kg)
 ,dthek_c2(n_cg)   & ! Increment to potential temp. due to conv in layer k
 ,dthekp1_c2(n_cg)   ! Increment to potential temp. due to conv in layer k+1

real(kind=real_umphys) ::                                                      &
  dtrae_c2(n_cg,nlev,ntra) ! Incr to model tracer due to conv. (kg/kg/s)

real(kind=real_umphys) ::                                                      &
  duek_c2(n_cg)     & ! Increment to model U in layer k due to CMT (m/s2)
 ,duekp1_c2(n_cg)   & ! Increment to model U in layer k+1 due to CMT (m/s2)
 ,dvek_c2(n_cg)     & ! Increment to model V in layer k due to CMT (m/s2)
 ,dvekp1_c2(n_cg)     ! Increment to model V in layer k+1 due to CMT (m/s2)

real(kind=real_umphys) ::                                                      &
  flxk_c(n_cg)    & ! Parcel mass flux in layer k (Pa/s)
 ,flxk_c2(n_cg)   & ! Parcel mass flux in layer k (Pa/s)
 ,flxkp12_c2(n_cg)& ! Half level mass flux (Pa/s)
 ,flxkp1_c2(n_cg)   ! Parcel mass flux in layer k+1 (Pa/s)

real(kind=real_umphys) :: flx_init_c2(n_cg)
                          ! Initial parcal mass flux at cloud base (Pa/s)

real(kind=real_umphys) ::                                                      &
  prekp1_c2(n_cg) ! Precip from parcel as it rises from layer k to k+1 (kg/m2/s)

real(kind=real_umphys) ::                                                      &
  qpk_c(n_cg)     & ! Parcel mixing ratio in layer k (kg/kg)
 ,qpk_c2(n_cg)    & ! Parcel mixing ratio in layer k (kg/kg)
 ,qpkp1_c(n_cg)   & ! Parcel mixing ratio in layer k+1 (kg/kg)
 ,qpkp1_c2(n_cg)  & ! Parcel mixing ratio in layer k+1 (kg/kg)
 ,qek_c(n_cg)     & ! Env. mixing ratio in layer k (kg/kg)
 ,qek_c2(n_cg)    & ! Env. mixing ratio in layer k (kg/kg)
 ,qek(n_cg)       & ! Env. mixing ratio in layer k (kg/kg)
 ,qekp1_c(n_cg)   & ! Env. mixing ratio in layer k+1 (kgkg-1)
 ,qekp1_c2(n_cg)  & ! Env. mixing ratio in layer k+1 (kgkg-1)
 ,qekp1(n_cg)     & ! Env. mixing ratio in layer k+1 (kgkg-1)
 ,qsek_c2(n_cg)   & ! Saturation mixing ratio of cld. env. in layer k (kg/kg)
 ,qsekp1_c2(n_cg)   ! Saturation mixing ratio of cld. env. in layer k+1 (kg/kg)

real(kind=real_umphys) :: quk_c(n_cg)
                            ! Undilute parcel mixing ratio in layer k (kg/kg)
real(kind=real_umphys) :: qukp1_c(n_cg)
                            ! Undilute parcel mixing ratio in layer k+1 (kg/kg)

real(kind=real_umphys) ::                                                      &
  thek_c(n_cg)    & ! Env. potential temp in layer k (K)
 ,thek_c2(n_cg)   & ! Env. potential temp in layer k (K)
 ,thek(n_cg)      & ! Env. potential temp in layer k (K)
 ,thekp1_c(n_cg)  & ! Env. potential temp i in layer k+1 (K)
 ,thekp1_c2(n_cg) & ! Env. potential temp i in layer k+1 (K)
 ,thekp1(n_cg)    & ! Env. potential temp i in layer k+1 (K)
 ,thpk_c(n_cg)    & ! Parcel potential temp in layer k (K)
 ,thpk_c2(n_cg)   & ! Parcel potential temp in layer k (K)
 ,thpkp1_c(n_cg)  & ! Parcel potential temp in layer k+1 (K)
 ,thpkp1_c2(n_cg)   ! Parcel potential temp in layer k+1 (K)

real(kind=real_umphys) :: thuk_c(n_cg)
                             ! Undilute parcel pot temp in layer k (K)
real(kind=real_umphys) :: thukp1_c(n_cg)
                             ! Undilute parcel pot temp in layer k+1 (K)

real(kind=real_umphys) :: trae_c2(n_cg,nlev,ntra)  ! Env. Tracer content (kg/kg)
real(kind=real_umphys) :: trap_c2(n_cg,nlev,ntra)
                                ! Parcel Tracer content (kg/kg)

real(kind=real_umphys) :: rbuoyk_c(n_cg), rbuoyk_c2(n_cg)
                                              ! Par. buoyancy at k (K)
real(kind=real_umphys) :: rbuoykp1_c(n_cg),rbuoykp1_c2(n_cg)
                                              ! Par. buoyancy at k+1 (K)
real(kind=real_umphys) :: rbuoyukp1_c(n_cg),rbuoyukp1_c2(n_cg)
                                              ! undilute Par. buoy at k+1 (K)

real(kind=real_umphys) :: watldek_c(n_cg), watldek_c2(n_cg)
                                              ! Env. water loading
                                              ! in layer k (kg/kg)
real(kind=real_umphys) :: watldpk_c(n_cg), watldpk_c2(n_cg)
                                              ! Par. water loading
                                              ! in layer k (kg/kg)
real(kind=real_umphys) :: watldekp1_c(n_cg), watldekp1_c2(n_cg)
                                              ! Env. water loading
                                              ! in layer k+1 (kg/kg)
real(kind=real_umphys) :: watldpkp1_c(n_cg), watldpkp1_c2(n_cg)
                                              ! Par. water loading
                                              ! in layer k+1 (kg/kg)

real(kind=real_umphys) :: Qlkp1_c(n_cg),                                       &
                             ! Amount of condensation to liquid water
        Qlkp1_c2(n_cg)       ! in the parcel (kg/kg)
real(kind=real_umphys) :: Qfkp1_c(n_cg),                                       &
                             ! Amount of deposition to ice water
        Qfkp1_c2(n_cg)       ! in the parcel (kg/kg)
real(kind=real_umphys) :: Frezkp1_c(n_cg), &   ! Amount of freezing from liquid
        Frezkp1_c2(n_cg)     ! to frozen water in the parcel (kg/kg)

real(kind=real_umphys) :: uek_c2(n_cg)    ! Env. U in layer k (m/s)
real(kind=real_umphys) :: uekp1_c2(n_cg)  ! Env. U in layer k+1 (m/s)
real(kind=real_umphys) :: vek_c2(n_cg)    ! Env. V in layer k (m/s)
real(kind=real_umphys) :: vekp1_c2(n_cg)  ! Env. V in layer k+1 (m/s)
real(kind=real_umphys) :: upk_c2(n_cg)    ! Parcel U in layer k (m/s)
real(kind=real_umphys) :: upkp1_c2(n_cg)  ! Parcel U in layer k+1 (m/s)
real(kind=real_umphys) :: vpk_c2(n_cg)    ! Parcel V in layer k (m/s)
real(kind=real_umphys) :: vpkp1_c2(n_cg)  ! Parcel V in layer k+1 (m/s)


!====================================================
! Local compressed arrays for calculation of
! parcel vertical velocity

! Compressed arrays for (Parcel vertical velocity)^2 on ...
real(kind=real_umphys) :: w2p_km1_c2 (n_cg)   ! ...layer centre k-1 [(m/s)^2]
real(kind=real_umphys) :: w2p_k_c2   (n_cg)   ! ...layer centre k   [(m/s)^2]
real(kind=real_umphys) :: w2p_kp1_c2 (n_cg)   ! ...layer centre k+1 [(m/s)^2]

! Compressed arrays for height above surface of model level ...
real(kind=real_umphys) :: zkm1_c2    (n_cg)   ! ...k-1   [m]
real(kind=real_umphys) :: zk_c2      (n_cg)   ! ...k     [m]
real(kind=real_umphys) :: zkp12_c2   (n_cg)   ! ...k+1/2 [m]
real(kind=real_umphys) :: zkp1_c2    (n_cg)   ! ...k+1   [m]
!====================================================

! PC2 compression arrays

real(kind=real_umphys) ::                                                      &
  qclek_c(n_cg)    & ! Env liquid condensate mixing ratio in layer k (kg/kg)
 ,qclek_c2(n_cg)   & ! Env liquid condensate mixing ratio in layer k (kg/kg)
 ,qclekp1_c(n_cg)  & ! Env liquid condensate mixing ratio in layer k+1 (kg/kg)
 ,qclekp1_c2(n_cg) & ! Env liquid condensate mixing ratio in layer k+1 (kg/kg)
 ,qcfek_c(n_cg)    & ! Env frozen condensate mixing ratio in layer k (kg/kg)
 ,qcfek_c2(n_cg)   & ! Env frozen condensate mixing ratio in layer k (kg/kg)
 ,qcfekp1_c(n_cg)  & ! Env frozen condensate mixing ratio in layer k+1 (kg/kg)
 ,qcfekp1_c2(n_cg) & ! Env frozen condensate mixing ratio in layer k+1 (kg/kg)
 ,qclpk_c(n_cg)    & ! Parcel liquid condensate mixing ratio in layer k (kg/kg)
 ,qclpk_c2(n_cg)   & ! Parcel liquid condensate mixing ratio in layer k (kg/kg)
 ,qclpkp1_c(n_cg)  & ! Parcel liquid condensate mixing ratio in layer k+1(kg/kg)
 ,qclpkp1_c2(n_cg) & ! Parcel liquid condensate mixing ratio in layer k+1(kg/kg)
 ,qcfpk_c(n_cg)    & ! Parcel frozen condensate mixing ratio in layer k (kg/kg)
 ,qcfpk_c2(n_cg)   & ! Parcel frozen condensate mixing ratio in layer k (kg/kg)
 ,qcfpkp1_c(n_cg)  & ! Parcel frozen condensate mixing ratio in layer k+1(kg/kg)
 ,qcfpkp1_c2(n_cg)   ! Parcel frozen condensate mixing ratio in layer k+1(kg/kg)

real(kind=real_umphys) ::                                                      &
  cflek_c2(n_cg)   & ! Environment liquid water cloud volume ( )
 ,cflekp1_c2(n_cg) & ! Environment liquid water cloud volume ( )
 ,cffek_c2(n_cg)   & ! Environment frozen water cloud volume ( )
 ,cffekp1_c2(n_cg) & ! Environment frozen water cloud volume ( )
 ,bcfek_c2(n_cg)   & ! Environment bulk total cloud volume ( )
 ,bcfekp1_c2(n_cg)   ! Environment bulk total cloud volume ( )

real(kind=real_umphys) ::                                                      &
  dqclek_c2(n_cg)  & ! Environment increments to liquid condensate mixing
 ,dqclekp1_c2(n_cg)& ! ratio to convection (kg/kg/s)
 ,dqcfek_c2(n_cg)  & ! Environment increments to frozen condensate mixing
 ,dqcfekp1_c2(n_cg)  !  ratio to convection (kg/kg/s)

real(kind=real_umphys) ::                                                      &
  dcflek_c2(n_cg)   & ! Environment increments to liquid water cloud
 ,dcflekp1_c2(n_cg) & ! volume due to convection (/s)
 ,dcffek_c2(n_cg)   & ! Environment increments to frozen water cloud
 ,dcffekp1_c2(n_cg) & ! volume due to convection (/s)
 ,dbcfek_c2(n_cg)   & ! Environment increments to bulk total cloud
 ,dbcfekp1_c2(n_cg)   ! volume due to convection (/s)

real(kind=real_umphys) :: amdetk_c2(n_cg)
logical :: bgmk_c2(n_cg)
logical :: bland_c2(n_cg)
logical :: blowst_c2(n_cg)
logical :: bterm_c2(n_cg)
real(kind=real_umphys) :: cape_c2(n_cg)
real(kind=real_umphys) :: fcape_c2(n_cg)
real(kind=real_umphys) :: cca_2d_c2(n_cg)
real(kind=real_umphys) :: cclwp_c2(n_cg)
real(kind=real_umphys) :: ccwkp1_c2(n_cg)
real(kind=real_umphys) :: dcpbydt_c2(n_cg)
real(kind=real_umphys) :: delpk_c2(n_cg)
real(kind=real_umphys) :: delpkp12_c2(n_cg)
real(kind=real_umphys) :: delpkp1_c2(n_cg)
real(kind=real_umphys) :: delp_uv_k_c2(n_cg)
real(kind=real_umphys) :: delp_uv_kp1_c2(n_cg)
real(kind=real_umphys) :: depth_c2(n_cg)
real(kind=real_umphys) :: dptot_c2(n_cg)
real(kind=real_umphys) :: deltaktot_c2(n_cg)
real(kind=real_umphys) :: eflux_u_ud_c2(n_cg)
real(kind=real_umphys) :: eflux_v_ud_c2(n_cg)
real(kind=real_umphys) :: ekp14_c(n_cg),ekp14_c2(n_cg)
real(kind=real_umphys) :: ekp34_c(n_cg),ekp34_c2(n_cg)
real(kind=real_umphys) :: exk_c(n_cg), exk_c2(n_cg)
real(kind=real_umphys) :: exkp1_c(n_cg),exkp1_c2(n_cg)
real(kind=real_umphys) :: expi_c2(n_cg)
integer :: icct_c2(n_cg)
integer :: iccb_c2(n_cg)
integer :: lctop_c2(n_cg)
integer :: lcbase_c2(n_cg)
real(kind=real_umphys) :: lcca_c2(n_cg)
real(kind=real_umphys) :: max_cfl_c2(n_cg)
real(kind=real_umphys) :: pk_c(n_cg),pk_c2(n_cg)
real(kind=real_umphys) :: pkp1_c(n_cg),pkp1_c2(n_cg)
real(kind=real_umphys) :: pstar_c2(n_cg)
real(kind=real_umphys) :: qpi_c2(n_cg)
real(kind=real_umphys) :: relh_c2(n_cg)
real(kind=real_umphys) :: tcw_c2(n_cg)
real(kind=real_umphys) :: thpi_c2(n_cg)
real(kind=real_umphys) :: xsbmin_v_c2(n_cg)
real(kind=real_umphys) :: qsat_lcl(n_cg)    ! not used
logical :: b_nodd(n_cg)   ! points with no downdraught
logical :: b_dd(n_cg)     ! points with downdraught on termination

!===============================================================
! CCRad Variables local variables As SHALLOW
!===============================================================

real(kind=real_umphys)   :: overlap_fac(n_cg)  ! Factor designed to improve
                             ! shallow Cu cover by allowing
                             ! for non-vertical clouds.

real(kind=real_umphys)   :: zpr         ! method (BL Fluxes) only if
                      !   l_ccrad = T .and. cca2d_sh_opt  = 1

!===============================================================
! End CCRad Variables local variables
!===============================================================

! Dummy variables for water tracer fields as water tracers have not been
! coded in congestus convection routine

integer, parameter :: n_wtrac_d = 1

type(conv_e_wtrac_type) :: wtrac_e(n_wtrac_d)
type(conv_p_wtrac_type) :: wtrac_p(n_wtrac_d)

! Loop counters

integer :: i,i2,j,k,ktra,kt

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CONGEST_CONV_6A'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

! Initialise logicals

do i = 1,n_cg
  blowst(i)    = .true.
  bterm(i)     = .false.
  bconv(i)     = .false.
  bcposs(i)    = .false.
  b_nodd(i)    = .false.
  b_dd(i)      = .false.
  blatent(i)   = .false.
end do

l_mom_gk = .false.    ! not using Gregory-Kershaw scheme
l_mom_gk_stable = .false.

!-----------------------------------------------------------------------
! 2.1  Initialise parcel properties and increment arrays
!-----------------------------------------------------------------------

!initialise parcel values over all levels
do k = 1, nlev
  do i = 1, n_cg
    qp(i,k)     = 0.0
    thp(i,k)    = 0.0
    qclp(i,k)   = 0.0
    qcfp(i,k)   = 0.0
    flx(i,k)    = 0.0
    precip(i,k) = 0.0
    ccw(i,k)    = 0.0
  end do
end do

if (cg_new_termc == term_undil) then
  do k = 1, nlev
    do i = 1, n_cg
      qu(i,k)     = 0.0
      thu(i,k)    = 0.0
    end do
  end do
end if

if (l_mom_gk) then
  do k=1,nlev
    do i = 1,n_cg
      up(i,k) = 0.0
      vp(i,k) = 0.0
    end do
  end do
end if

if (l_tracer) then
  do ktra = 1,ntra
    do k=1,nlev
      do i = 1,n_cg
        trap(i,k,ktra) = 0.0
      end do
    end do
  end do
end if

do k = 1,nlev
  do i = 1,n_cg
    dthbydt(i,k)  = 0.0
    dqbydt(i,k)   = 0.0
    dqclbydt(i,k) = 0.0
    dqcfbydt(i,k) = 0.0
    dbcfbydt(i,k) = 0.0
    dcflbydt(i,k) = 0.0
    dcffbydt(i,k) = 0.0
  end do
end do

if (l_mom) then
  do k = 1,nlev
    do i = 1,n_cg
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    end do
  end do
end if  ! L_mom

if (l_tracer) then
  do ktra = 1,ntra
    do k = 1,nlev
      do i = 1,n_cg
        dtrabydt(i,k,ktra) = 0.0
      end do
    end do
  end do
end if  ! L_tracer

!Initialise the termination level
do i =1, n_cg
  kterm(i) = 0
end do

! Allocate water tracer arrays - these are just dummy fields
call wtrac_alloc_conv_e(1,1,1,wtrac_e)
call wtrac_alloc_conv_p(1,1,1,wtrac_p)

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------
if (flg_up_flx) then
  do k = 1,nlev
    do i = 1,n_cg
      up_flux(i,k)      = 0.0
    end do
  end do
end if
if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,n_cg
      up_flux_half(i,k) = 0.0
    end do
  end do
end if
if (flg_dwn_flx) then
  do k = 1,nlev
    do i = 1,n_cg
      dwn_flux(i,k)     = 0.0
    end do
  end do
end if
if (flg_entr_up) then
  do k = 1,nlev
    do i = 1,n_cg
      entrain_up(i,k)   = 0.0
    end do
  end do
end if
if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,n_cg
      detrain_up(i,k)   = 0.0
    end do
  end do
end if
if (flg_entr_dwn) then
  do k = 1,nlev
    do i = 1,n_cg
      entrain_dwn(i,k)  = 0.0
    end do
  end do
end if
if (flg_detr_dwn) then
  do k = 1,nlev
    do i = 1,n_cg
      detrain_dwn(i,k)  = 0.0
    end do
  end do
end if
if (l_mom) then
  if (flg_uw_shall) then
    do k = 1,nlev
      do i = 1,n_cg
        uw_shall(i,k)   = 0.0
      end do
    end do
  end if
  if (flg_vw_shall) then
    do k = 1,nlev
      do i = 1,n_cg
        vw_shall(i,k)   = 0.0
      end do
    end do
  end if
end if  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
do i = 1,n_cg
  cca_2d(i) = 0.0
  iccb(i)   = 0
  icct(i)   = 0
  tcw(i)    = 0.0
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
end do

do k= 1, n_cca_lev
  do i= 1, n_cg
    cca(i,k) = 0.0
  end do
end do


do i = 1,n_cg
  start_lev(i)        = ntml(i)

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
  eflux_u_ud(i) = 0.0
  eflux_v_ud(i) = 0.0

  !-----------------------------------------------------------------------
  ! 2.7  Initialise surface precipitation arrays
  !-----------------------------------------------------------------------
  rain(i) = 0.0
  snow(i) = 0.0
end do

!-----------------------------------------------------------------------
! Calculate parcel perturbations
!-----------------------------------------------------------------------

! Calculate XSBMIN and THPIXS constants based on layer thickness (Pa)
do k = 1,nlev-1
  do i = 1,n_cg
    xsbmin_v(i,k) = min( ((p_layer_centres(i,k) -                              &
              p_layer_centres(i,k+1))/5000.0),1.0) *0.2

    thpixs_v(i,k) = min( ((p_layer_centres(i,k) -                              &
              p_layer_centres(i,k+1))/5000.0),1.0) * thpixs_shallow

    qpixs_v(i,k)  = qpixs_shallow
  end do
end do  ! nlev

! Calculate convective velocity scale and cloud base mass flux
do i = 1,n_cg
  wsc(i) = (delthvu(i) * c_mass * wstar(i) * g / (th(i,ntml(i))                &
             * (1.0 + c_virtual * q(i,ntml(i)))))**0.3333
  mb(i)  = c_mass * wstar(i)
  zcld(i) = ztop_uv(i) - zlcl_uv(i)
  wsc_o_mb(i) = wsc(i)/mb(i)
  weight_param(i) = 1.0              ! no grey scale changes
end do

! Define the LCL at the half level above ntml. Find environmental
! T at p_lcl by approximating theta there with
! th(i,k) + constant*(th(i,k+1)-th(i,k)) where constant is tunable.
! Similarly for q.
do i = 1,n_cg
  k =ntml(i)
  p_lcl(i)  = p_layer_boundaries(i,k)
  th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
  t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
  q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
end do

! Calculate saturation mixing ratio at LCL

call qsat(qse_lcl,t_lcl,p_lcl,n_cg)

! Calculate theta and q perturbation (perturbation is based on
! environment buoyancy gradient)
! Re-set thpixs and qpixs at ntml
do i = 1,n_cg
  if (t_lcl(i) >  tm) then
    dq_sat_env = repsilon * lc * qse_lcl(i) / (r * t_lcl(i) * t_lcl(i))
  else
    dq_sat_env = repsilon * (lc+lf) * qse_lcl(i) / (r * t_lcl(i) * t_lcl(i))
  end if

  b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0                           &
               + c_virtual * qse_lcl(i)

  thv_pert(i) = -0.17 * wthvs(i) / mb(i) +                                     &
              ( th(i,ntml(i)+1) *(1.0 + c_virtual*q(i,ntml(i)+1))              &
              - th(i,ntml(i))  * (1.0 + c_virtual*q(i,ntml(i))) )

  c_calc   = th_lcl(i) * c_virtual * (qse_lcl(i) - q_lcl(i)) - thv_pert(i)

  thpert(i)   = -c_calc / b_calc  !ignore term in THPERT**2

  thpixs_v(i,ntml(i)) = thpert(i)

  qpert(i)  = qse_lcl(i) + ((p_lcl(i) / pref)**kappa)                          &
                             * thpert(i) * dq_sat_env - q_lcl(i)

  qpixs_v(i,ntml(i))  = qpert(i)

end do !n_cg


! Set bwater=.true. on points where water will condense rather than
! ice.
call flag_wet(n_cg,n_cg,nlev,th,exner_layer_centres,bwater)


!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------
! To reduce cost limit level loop to NTPAR_MAX maximum NTPAR value.
! NTPAR is the top of the parcel ascent for shallow convection.

if (cg_on == 1) then   ! Adaptive forced detrainment on
                       ! No limit on convection top
  ntpar_max = nlev-3   ! What is a sensible value to have here?

else                   ! Top limited
  ntpar_max=0
  do i = 1,n_cg
    if (ntpar(i)+1 > ntpar_max) then
      ntpar_max=ntpar(i)+1
    end if
  end do
  ! Ensure that ntpar_max does not exceed nlev-1
  ntpar_max = min(ntpar_max, nlev-1)
end if

do k = 2,ntpar_max  !loop over model levels
  !-----------------------------------------------------------------------
  ! Initialise environment variables
  ! NB These variable are only used by layer_cn.
  !-----------------------------------------------------------------------
  do i = 1,n_cg
    thek(i)   = th(i,k)
    qek(i)    = q(i,k)
    thekp1(i) = th(i,k+1)
    qekp1(i)  = q(i,k+1)
    !Note that unlike p_layer_boundaries, where k indexing is offset
    !by one compared to the dynamics numbering, z retains the numbering
    !convention for dynamics variables i.e. for theta levels, k->k
    !and for rho levels k+1/2 -> k+1
    zkm1(i)   = z_theta (i,k-1)
    zk(i)     = z_theta(i,k)
    zkp12(i)  = z_rho(i, k+1)
    zkp1(i)   = z_theta(i, k+1)
    rhum(i)   = q(i,k) / qse(i,k)
  end do

  !-----------------------------------------------------------------------
  ! Initialise parcel properties (theta,q,tracer,momentum) if convection
  ! is not occurring at level k and has not convected in column before
  !-----------------------------------------------------------------------
  do i = 1,n_cg
    if ( .not. bconv(i) .and. det_lev(i) == 0) then
      expi(i)     = exner_layer_centres(i,k)
      bgmk(i)     = .false.
      depth(i)    = 0.0
      thpi(i)     = th(i,k) + thpixs_v(i,k)
      thp(i,k)    = thpi(i)
      qpi(i)      = q(i,k)  + qpixs_v(i,k)
      qp(i,k)     = qpi(i)
      if (cg_new_termc == term_undil) then
        thu(i,k)  = thpi(i)
        qu(i,k)   = qpi(i)
      end if
      if (l_q_interact) then
        qclp(i,k) = qcl(i,k)
        qcfp(i,k) = qcf(i,k)
      else
        qclp(i,k) = 0.0
        qcfp(i,k) = 0.0
      end if
      if (l_mom_gk) then
        up(i,k)   = u(i,k)
        vp(i,k)   = v(i,k)
      end if
    end if
  end do  ! n_cg
  if (l_tracer) then
    do ktra=1,ntra
      do i = 1,n_cg
        if ( .not. bconv(i)) then
          trap(i,k,ktra)  = tracer(i,k,ktra)
        end if  !not bconv
      end do
    end do
  end if

  !-----------------------------------------------------------------------
  ! 3.1  Calculate layer dependent constants (pressure,
  !      layer thickness, entrainment coefficients, detrainment
  !      coefficients)
  !-----------------------------------------------------------------------

  call layer_cn_6a(k, n_cg, nlev,                                              &
                   ntml, ntpar, start_lev,                                     &
                   exner_layer_centres,                                        &
                   p_layer_boundaries, p_layer_centres,                        &
                   z_rho,                                                      &
                   conv_prog_precip,                                           &
                   recip_pstar, entrain_coef, rhum,                            &
                   ccp_strength,                                               &
                   zk, zkp12, zkp1,                                            &
                   wsc_o_mb, qsat_lcl, w_max,                                  &
                   congestus,                                                  &
                   bconv,                                                      &
                   ! Out
                   pk, pkp1, exk, exkp1,                                       &
                   delpk, delpkp12, delpkp1,                                   &
                   delp_uv_k, delp_uv_kp1,                                     &
                   ekp14, ekp34, amdetk                                        &
                   )

  ! Maximum initial convective mass flux
  do i = 1,n_cg
    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)
  end do

  !-----------------------------------------------------------------------
  ! Initial test to check if convection is possible in layer k
  !-----------------------------------------------------------------------
  ! Convection is possible if
  ! - the point was convecting (bconv = .T.) and did not terminate
  !   in the previous layer
  ! - or if at the top level of the surface mixed layer (k = ntml)
  do i = 1,n_cg
    bcposs(i) = bconv(i) .or. k  ==  ntml(i)
  end do  ! n_cg

  ! Calculate number of points which may convect (ncposs) and
  ! set compression indices (index1)
  ncposs = 0
  do i = 1,n_cg
    idx_contig(i) = i
    if (bcposs(i)) then
      ncposs          = ncposs + 1
      index1(ncposs)  = i
    end if
  end do

  !-----------------------------------------------------------------------
  ! Compress points where convection may occur
  ! NB This process is used to update some single level variables that are
  ! defined on level k and kp1 by using the full field variables.
  ! NB The order in which the variables are compressed
  ! is the same as the argument list for LIFT_PAR
  !-----------------------------------------------------------------------
  if (ncposs  >   0) then
    do i = 1,ncposs
      !intent(in) for lift_par/water_loading
      thek_c(i)     = th(index1(i),k)
      thekp1_c(i)   = th(index1(i),k+1)
      qek_c(i)      = q(index1(i),k)
      qekp1_c(i)    = q(index1(i),k+1)
      qclek_c(i)    = qcl(index1(i),k)
      qcfek_c(i)    = qcf(index1(i),k)
      qclekp1_c(i)  = qcl(index1(i),k+1)
      qcfekp1_c(i)  = qcf(index1(i),k+1)
      pk_c(i)       = pk(index1(i))
      pkp1_c(i)     = pkp1(index1(i))
      exkp1_c(i)    = exkp1(index1(i))
      thpk_c(i)     = thp(index1(i),k)
      qpk_c(i)      = qp(index1(i),k)
      qclpk_c(i)    = qclp(index1(i),k)
      qcfpk_c(i)    = qcfp(index1(i),k)
      ekp14_c(i)    = ekp14(index1(i))
      ekp34_c(i)    = ekp34(index1(i))
      bwk_c(i)      = bwater(index1(i),k)
      bwkp1_c(i)    = bwater(index1(i),k+1)
      !intent(in) for water_loading only
      exk_c(i)      = exk(index1(i))
    end do
  end if  ! ncposs>0

  if (cg_new_termc == term_undil .and. ncposs > 0) then
    do i = 1,ncposs
      !intent(in) for undilute lift_par
      thuk_c(i)     = thu(index1(i),k)
      quk_c(i)      = qu(index1(i),k)
    end do
  end if

  !-----------------------------------------------------------------------
  ! 3.2  Lift parcel from layer k to layer k+1
  !-----------------------------------------------------------------------
  if (ncposs > 0) then

    call lift_par_6a (k, ncposs, n_wtrac_d, thek_c, thekp1_c,                  &
                qek_c, qekp1_c, qclek_c, qcfek_c,                              &
                qclekp1_c, qcfekp1_c,                                          &
                pkp1_c, exkp1_c,                                               &
                thpk_c, qpk_c, qclpk_c, qcfpk_c,                               &
                ekp14_c, ekp34_c,                                              &
                l_q_interact, bwkp1_c, wtrac_e,                                &
                !InOut
                wtrac_p,                                                       &
                !Out
                bgmkp1_c, thpkp1_c, qpkp1_c,                                   &
                qclpkp1_c, qcfpkp1_c,                                          &
                Qlkp1_c, Qfkp1_c, Frezkp1_c, idx_contig, ncposs)

    !-----------------------------------------------------------------------
    ! Calculate the water loading for level k and k+1
    !-----------------------------------------------------------------------
    call water_loading(ncposs, qclek_c, qcfek_c, qclpk_c, qcfpk_c,             &
                     watldek_c, watldpk_c, idx_contig, ncposs)

    call water_loading(ncposs, qclekp1_c, qcfekp1_c, qclpkp1_c, qcfpkp1_c,     &
                     watldekp1_c, watldpkp1_c, idx_contig, ncposs)

    if (cg_new_termc == term_undil) then
      !---------------------------------------------------------------------
      ! Calculate undilute parcel properties
      !---------------------------------------------------------------------
      call lift_undil_par_6a(ncposs,                                           &
                  pkp1_c, exkp1_c,                                             &
                  thuk_c, quk_c,                                               &
                  Frezkp1_c, bwkp1_c,                                          &
                  !Out
                  thukp1_c, qukp1_c, idx_contig, ncposs)

      ! Update undilute profile
      do i = 1,ncposs ! Loop over points which may convect
        thu(index1(i),k+1) = thukp1_c(i)
        qu(index1(i),k+1)  = qukp1_c(i)
      end do
    end if

    ! NEC compiler directive
    !CDIR NODEP

    !-----------------------------------------------------------------------
    ! Test if convection is starting from layer k
    !-----------------------------------------------------------------------
    do i = 1,ncposs ! Loop over points which may convect

      ! Calculate buoyancy (virt. potential temp.) of parcel in layer k and k+1
      rbuoyk_c(i)   = thpk_c(i) * (1.0 + c_virtual *qpk_c(i)                   &
                    - watldpk_c(i))                                            &
                    - thek_c(i) * (1.0 + c_virtual *qek_c(i)                   &
                    - watldek_c(i))

      rbuoykp1_c(i) = thpkp1_c(i) * (1.0 + c_virtual *qpkp1_c(i)               &
                    - watldpkp1_c(i))                                          &
                    - thekp1_c(i) * (1.0 + c_virtual *qekp1_c(i)               &
                    - watldekp1_c(i))

      if (cg_new_termc == term_undil) then
        !undilute parcel buoyancy
        rbuoyukp1_c(i)= thukp1_c(i) * (1.0 + c_virtual *qukp1_c(i)             &
                      - watldpkp1_c(i))                                        &
                      - thekp1_c(i) * (1.0 + c_virtual *qekp1_c(i)             &
                      - watldekp1_c(i))
      end if

      ! Allow parcel to convect from ntml.
      if (k  ==  ntml(index1(i))) then
        bconv(index1(i))  = .true.  ! convection active
        blowst(index1(i)) = .true.  ! convection initialised in layer

        ! Set parcel mass flux
        flxk_c(i)         = mb(index1(i)) * g                                  &
                          * p_layer_centres(index1(i),k)                       &
                          / ( r * thpk_c(i)                                    &
                          * (p_layer_centres(index1(i),k)/ pref)**kappa )

        ! Write compressed mass flux back to full array
        flx(index1(i),k)  = flxk_c(i)

        ! Store diagnostics linked to initial convective mass flux for
        ! calculation of final closure.
        flx_init(index1(i))    = flxk_c(i)
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
      index2(nconv) = i
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
    do i = 1,nconv
      start_lev_c2(i)   = ntml(index1(index2(i)))
      pstar_c2(i)       = pstar(index1(index2(i)))
      pk_c2(i)          = pk_c(index2(i))
      pkp1_c2(i)        = pkp1_c(index2(i))
      delpk_c2(i)       = delpk(index1(index2(i)))
      delpkp1_c2(i)     = delpkp1(index1(index2(i)))
      delpkp12_c2(i)    = delpkp12(index1(index2(i)))
      delp_uv_k_c2(i)   = delp_uv_k(index1(index2(i)))
      delp_uv_kp1_c2(i) = delp_uv_kp1(index1(index2(i)))
      exk_c2(i)         = exk_c(index2(i))
      exkp1_c2(i)       = exkp1_c(index2(i))
      thek_c2(i)        = thek_c(index2(i))
      thekp1_c2(i)      = thekp1_c(index2(i))
      qek_c2(i)         = qek_c(index2(i))
      qekp1_c2(i)       = qekp1_c(index2(i))
      qclek_c2(i)       = qclek_c(index2(i))
      qclekp1_c2(i)     = qclekp1_c(index2(i))
      qcfek_c2(i)       = qcfek_c(index2(i))
      qcfekp1_c2(i)     = qcfekp1_c(index2(i))
      qsek_c2(i)        = qse(index1(index2(i)),k)
      qsekp1_c2(i)      = qse(index1(index2(i)),k+1)
      cflek_c2(i)       = cf_liquid(index1(index2(i)),k)
      cflekp1_c2(i)     = cf_liquid(index1(index2(i)),k+1)
      cffek_c2(i)       = cf_frozen(index1(index2(i)),k)
      cffekp1_c2(i)     = cf_frozen(index1(index2(i)),k+1)
      bcfek_c2(i)       = bulk_cf(index1(index2(i)),k)
      bcfekp1_c2(i)     = bulk_cf(index1(index2(i)),k+1)
      thpk_c2(i)        = thpk_c(index2(i))
      qpk_c2(i)         = qpk_c(index2(i))
      qclpk_c2(i)       = qclpk_c(index2(i))
      qcfpk_c2(i)       = qcfpk_c(index2(i))
      thpi_c2(i)        = thpi(index1(index2(i)))
      qpi_c2(i)         = qpi(index1(index2(i)))
      expi_c2(i)        = expi(index1(index2(i)))
      rbuoyk_c2(i)      = rbuoyk_c(index2(i))
      rbuoykp1_c2(i)    = rbuoykp1_c(index2(i))
      rbuoyukp1_c2(i)   = rbuoyukp1_c(index2(i))
      xsbmin_v_c2(i)    = xsbmin_v(index1(index2(i)),k)
      watldek_c2(i)     = watldekp1_c(index2(i))
      watldekp1_c2(i)   = watldekp1_c(index2(i))
      watldpk_c2(i)     = watldpkp1_c(index2(i))
      watldpkp1_c2(i)   = watldpkp1_c(index2(i))
      Qlkp1_c2(i)       = Qlkp1_c(index2(i))
      Qfkp1_c2(i)       = Qfkp1_c(index2(i))
      Frezkp1_c2(i)     = Frezkp1_c(index2(i))
      ekp14_c2(i)       = ekp14_c(index2(i))
      ekp34_c2(i)       = ekp34_c(index2(i))
      amdetk_c2(i)      = amdetk(index1(index2(i)))
      flxk_c2(i)        = flx(index1(index2(i)),k)
      flx_init_c2(i)    = flx_init(index1(index2(i)))
    end do
    if (l_mom_gk) then
      do i = 1,nconv
        uek_c2(i)       = u(index1(index2(i)),k)
        uekp1_c2(i)     = u(index1(index2(i)),k+1)
        vek_c2(i)       = v(index1(index2(i)),k)
        vekp1_c2(i)     = v(index1(index2(i)),k+1)
        upk_c2(i)       = up(index1(index2(i)),k)
        vpk_c2(i)       = vp(index1(index2(i)),k)
      end do
    end if
    if (l_tracer) then
      do ktra = 1,ntra
        do i = 1,nconv
          trae_c2(i,k,ktra)   = tracer(index1(index2(i)),k,ktra)
          trae_c2(i,k+1,ktra) = tracer(index1(index2(i)),k+1,ktra)
          trap_c2(i,k,ktra)  = trap(index1(index2(i)),k,ktra)
        end do
      end do
    end if
    if (flg_w_eqn) then
      do i=1, nconv
        zkm1_c2(i)      = zkm1  (index1(index2(i)))
        zk_c2(i)        = zk    (index1(index2(i)))
        zkp12_c2(i)     = zkp12 (index1(index2(i)))
        zkp1_c2(i)      = zkp1  (index1(index2(i)))
        w2p_km1_c2(i)   = w2p(index1(index2(i)),k-1)
      end do
    end if
    do i = 1,nconv
      bgmk_c2(i)        = bgmk(index1(index2(i)))
      bgmkp1_c2(i)      = bgmkp1_c(index2(i))
      bwk_c2(i)         = bwk_c(index2(i))
      bwkp1_c2(i)       = bwkp1_c(index2(i))
      blowst_c2(i)      = blowst(index1(index2(i)))
      bland_c2(i)       = bland(index1(index2(i)))
    end do
    !Compression for intent(INOUT)
    do i = 1,nconv
      lcbase_c2(i)      = lcbase(index1(index2(i)))
      lctop_c2(i)       = lctop(index1(index2(i)))
      thpkp1_c2(i)      = thpkp1_c(index2(i))
      qpkp1_c2(i)       = qpkp1_c(index2(i))
      qclpkp1_c2(i)     = qclpkp1_c(index2(i))
      qcfpkp1_c2(i)     = qcfpkp1_c(index2(i))
      dthek_c2(i)       = dthbydt(index1(index2(i)),k)
      dqek_c2(i)        = dqbydt(index1(index2(i)),k)
      dqclek_c2(i)      = dqclbydt(index1(index2(i)),k)
      dqcfek_c2(i)      = dqcfbydt(index1(index2(i)),k)
      tcw_c2(i)         = tcw(index1(index2(i)))
      depth_c2(i)       = depth(index1(index2(i)))
      cclwp_c2(i)       = cclwp(index1(index2(i)))
      lcca_c2(i)        = lcca(index1(index2(i)))
      cape_c2(i)        = cape(index1(index2(i)))
      fcape_c2(i)       = fcape(index1(index2(i)))
      dcpbydt_c2(i)     = dcpbydt(index1(index2(i)))
      relh_c2(i)        = 0.0 ! dummy variable
      dptot_c2(i)       = 0.0 ! dummy variable
      deltaktot_c2(i)   = deltaktot(index1(index2(i)))
      max_cfl_c2(i)     = max_cfl(index1(index2(i)))
    end do
    if (l_mom_gk) then
      do i = 1,nconv
        eflux_u_ud_c2(i)= eflux_u_ud(index1(index2(i)))
        eflux_v_ud_c2(i)= eflux_v_ud(index1(index2(i)))
        duek_c2(i)      = dubydt(index1(index2(i)),k)
        dvek_c2(i)      = dvbydt(index1(index2(i)),k)
      end do
    end if
    if (l_tracer) then
      do ktra = 1,ntra
        do i = 1,nconv
          dtrae_c2(i,k,ktra) = dtrabydt(index1(index2(i)),k,ktra)
        end do
      end do
    end if
    if (flg_w_eqn) then
      do i=1, nconv
        w2p_k_c2(i)     = w2p(index1(index2(i)),k)
      end do
    end if
    do i = 1,nconv
      bterm_c2(i)       = .false.
      blatent_c2(i)     = blatent(index1(index2(i)))
    end do
    !Compression for intent(out)
    ! Most intent(out) variables do not need to be initialised because they
    ! are always set in CONVEC2 or if they are not set in CONVEC2 then they
    ! are not used. However, several of the cloud variables do need to be
    ! initialised.
    do i = 1,nconv
      iccb_c2(i)        = iccb(index1(index2(i)))
      icct_c2(i)        = icct(index1(index2(i)))
      cca_2d_c2(i)      = cca_2d(index1(index2(i)))
    end do

    ! Force congestus convection to stop at the parcel top from
    ! conv_diag (ntpar) unless using adaptive forced detrainment.
    if (cg_on == 0) then
      do i = 1,nconv
        if (k  ==  ntpar(index1(index2(i)))) then
          bterm_c2(i) = .true.
        end if
      end do

    end if
    !-----------------------------------------------------------------------
    ! 3.3  Calculate the rest of the parcel ascent  and the effect of
    !      convection on the large-scale atmosphere.

    !      Subroutine CONVEC2

    !      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
    !-----------------------------------------------------------------------

    call convec2_6a  (k, nconv, n_cg, nlev, ntra, n_wtrac_d, nlev,             &
                      cg_on, cg_new_termc, start_lev_c2, timestep,             &
                      pk_c2, pkp1_c2, delpk_c2,                                &
                      delpkp1_c2, delp_uv_k_c2, delp_uv_kp1_c2,                &
                      exk_c2, exkp1_c2,                                        &
                      thek_c2, thekp1_c2, qek_c2, qekp1_c2,                    &
                      qclek_c2, qclekp1_c2, qcfek_c2, qcfekp1_c2,              &
                      qsek_c2, qsekp1_c2,                                      &
                      cflek_c2, cflekp1_c2,  cffek_c2,  cffekp1_c2,            &
                      thpk_c2, qpk_c2, qclpk_c2, qcfpk_c2,                     &
                      rbuoyk_c2, rbuoykp1_c2, rbuoyukp1_c2,                    &
                      watldek_c2, watldekp1_c2, watldpk_c2, watldpkp1_c2,      &
                      Qlkp1_c2, Qfkp1_c2, Frezkp1_c2,                          &
                      ekp14_c2, ekp34_c2, amdetk_c2, flxk_c2, flx_init_c2,     &
                      uek_c2, uekp1_c2, vek_c2, vekp1_c2,                      &
                      upk_c2, vpk_c2,                                          &
                      trae_c2,                                                 &
                      zk_c2, zkp12_c2, zkp1_c2,                                &
                      l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,       &
                      bgmk_c2, bgmkp1_c2, bwk_c2,                              &
                      bwkp1_c2, blowst_c2, bland_c2,                           &

                      ! In/out
                      lcbase_c2, lctop_c2,                                     &
                      thpkp1_c2, qpkp1_c2, qclpkp1_c2, qcfpkp1_c2,             &
                      dthek_c2, dqek_c2, dqclek_c2, dqcfek_c2,                 &
                      tcw_c2, depth_c2, cclwp_c2, lcca_c2,                     &
                      cape_c2, fcape_c2, dcpbydt_c2,                           &
                      relh_c2, dptot_c2, deltaktot_c2, max_cfl_c2,             &
                      eflux_u_ud_c2, eflux_v_ud_c2,                            &
                      duek_c2, dvek_c2,                                        &
                      dtrae_c2, trap_c2,                                       &
                      wtrac_e, wtrac_p,                                        &
                      w2p_k_c2, bterm_c2, blatent_c2, xsbmin_v_c2,             &

                      ! Out
                      iccb_c2, icct_c2,                                        &
                      dcflek_c2, dcffek_c2, dbcfek_c2,                         &
                      dthekp1_c2, dqekp1_c2, dqclekp1_c2, dqcfekp1_c2,         &
                      dcflekp1_c2, dcffekp1_c2, dbcfekp1_c2,                   &
                      prekp1_c2, thrk_c2, qrk_c2, deltak_c2,                   &
                      flxkp12_c2, flxkp1_c2,                                   &
                      cca_2d_c2, ccwkp1_c2,                                    &
                      upkp1_c2, vpkp1_c2,                                      &
                      duekp1_c2, dvekp1_c2,                                    &
                      w2p_kp1_c2,tnuc_nlcl,                                    &

                      !Indirect indexing
                      idx_contig,nconv )

  end if ! nconv > 0

  !-----------------------------------------------------------------------
  ! Decompression of compressed variables coming out of
  ! of CONVEC2.
  ! NB The order in which the variables are decompressed
  ! is the same as the the argument list for CONVEC2.
  !-----------------------------------------------------------------------
  do i = 1,n_cg
    depth(i)      = 0.0
    bgmk(i)       = .false.
    bterm(i)      = .false.
  end do

  if (nconv  >   0) then
    !Decompression for intent(in)
    do i = 1,nconv
      bgmk(index1(index2(i)))         = bgmkp1_c2(i)
    end do
    !Decompression for intent(INOUT)
    do i = 1,nconv
      lcbase(index1(index2(i)))       = lcbase_c2(i)
      lctop(index1(index2(i)))        = lctop_c2(i)
      thp(index1(index2(i)),k+1)      = thpkp1_c2(i)
      qp(index1(index2(i)),k+1)       = qpkp1_c2(i)
      qclp(index1(index2(i)),k+1)     = qclpkp1_c2(i)
      qcfp(index1(index2(i)),k+1)     = qcfpkp1_c2(i)
      dthbydt(index1(index2(i)),k)    = dthek_c2(i)
      dqbydt(index1(index2(i)),k)     = dqek_c2(i)
      dqclbydt(index1(index2(i)),k)   = dqclek_c2(i)
      dqcfbydt(index1(index2(i)),k)   = dqcfek_c2(i)
      tcw(index1(index2(i)))          = tcw_c2(i)
      depth(index1(index2(i)))        = depth_c2(i)
      cclwp(index1(index2(i)))        = cclwp_c2(i)
      lcca(index1(index2(i)))         = lcca_c2(i)
      cape(index1(index2(i)))         = cape_c2(i)
      fcape(index1(index2(i)))        = fcape_c2(i)
      dcpbydt(index1(index2(i)))      = dcpbydt_c2(i)
      !dummy                          = relh_c2(i)
      !dummy                          = dptot_c2(i)
      deltaktot(index1(index2(i)))    = deltaktot_c2(i)
      max_cfl(index1(index2(i)))      = max_cfl_c2(i)
    end do
    if (l_mom_gk) then
      do i = 1,nconv
        eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
        eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
        dubydt(index1(index2(i)),k)   = duek_c2(i)
        dvbydt(index1(index2(i)),k)   = dvek_c2(i)
      end do
    end if
    if (l_tracer) then
      do i = 1,nconv
        do ktra = 1,ntra
          dtrabydt(index1(index2(i)),k,ktra)    = dtrae_c2(i,k,ktra)
        end do
      end do
    end if
    if (flg_w_eqn) then
      do i=1, nconv
        w2p(index1(index2(i)),k)      = w2p_k_c2(i)
      end do
    end if
    do i = 1,nconv
      bterm(index1(index2(i)))        = bterm_c2(i)
      blatent(index1(index2(i)))      = blatent_c2(i)
    end do
    !Decompression for intent(out)
    do i = 1,nconv
      iccb(index1(index2(i)))         = iccb_c2(i)
      icct(index1(index2(i)))         = icct_c2(i)
      dcflbydt(index1(index2(i)),k)   = dcflek_c2(i)
      dcffbydt(index1(index2(i)),k)   = dcffek_c2(i)
      dbcfbydt(index1(index2(i)),k)   = dbcfek_c2(i)
      dthbydt(index1(index2(i)),k+1)  = dthekp1_c2(i)
      dqbydt(index1(index2(i)),k+1)   = dqekp1_c2(i)
      dqclbydt(index1(index2(i)),k+1) = dqclekp1_c2(i)
      dqcfbydt(index1(index2(i)),k+1) = dqcfekp1_c2(i)
      dcflbydt(index1(index2(i)),k+1) = dcflekp1_c2(i)
      dcffbydt(index1(index2(i)),k+1) = dcffekp1_c2(i)
      dbcfbydt(index1(index2(i)),k+1) = dbcfekp1_c2(i)
      precip(index1(index2(i)),k+1)   = prekp1_c2(i)
      !dummy                          = deltak_c2(i)
      !dummy                          = flxkp12_c2(i)
      flx(index1(index2(i)),k+1)      = flxkp1_c2(i)
      cca_2d(index1(index2(i)))       = cca_2d_c2(i)
      ccw(index1(index2(i)),k+1)      = ccwkp1_c2(i)
    end do
    if (l_mom_gk) then
      do i = 1,nconv
        up(index1(index2(i)),k+1)     = upkp1_c2(i)
        vp(index1(index2(i)),k+1)     = vpkp1_c2(i)
        dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
        dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
      end do
    end if
    if (l_tracer) then
      do i = 1,nconv
        do ktra = 1,ntra
          trap(index1(index2(i)),k+1,ktra)      = trap_c2(i,k+1,ktra)
          dtrabydt(index1(index2(i)),k+1,ktra)  = dtrae_c2(i,k+1,ktra)
        end do
      end do
    end if
    if (flg_w_eqn) then
      do i=1, nconv
        w2p(index1(index2(i)),k+1)    = w2p_kp1_c2(i)
      end do
    end if

  end if      ! nconv > 0

  !-----------------------------------------------------------------------
  ! 3.4  CFL scaling
  !-----------------------------------------------------------------------

  ! Set up integer nterm which is the total number of points where
  ! convection has terminated.
  ! Index to full array (n_cg) with index_nterm


  nterm = 0
  do i = 1,n_cg
    if (bterm(i)) then
      nterm = nterm + 1
      index_nterm(nterm) = i
    end if
  end do

  if (nterm >  0) then


    ! Work out scaled mass flux needed to keep cfl ratio below limit.
    ! Note L_CAPE not applied to shallow convection

    do j = 1,nterm
      i = index_nterm(j)

      max_cfl(i) = max_cfl(i) * timestep
      if (max_cfl(i)  >   cfl_limit) then
        flx_init_new(i) = flx_init(i) * cfl_limit / max_cfl(i)
      else
        flx_init_new(i) = flx_init(i)
      end if

      if (flx_init_new(i)  >   flxmax_init(i)) then
        flx_init_new(i) = flxmax_init(i)
      end if
      max_cfl(i) = 0.0
    end do      ! j (nterm)

    ! Scale cloud fraction

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
          if (l_mom_gk) then
            dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
            dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
          end if
          if (l_tracer) then
            do ktra = 1,ntra
              dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
            end do
          end if

          flx(i,kt)    = flx(i,kt)    *  scale_f(i)
          precip(i,kt) = precip(i,kt) *  scale_f(i)

        end if !kt >ntml and flx_init_new >0
      end do  ! j loop
    end do  ! kt loop

    !-----------------------------------------------------------------------
    ! 3.5  Setup dowdnruaght arrays used by original scheme b_dd & b_nodd
    !-----------------------------------------------------------------------
    do i = 1,nterm
      i2=index_nterm(i)
      tempnum=0.0
      if (iccb(i2) >  0) then
        deltap_cld=p_layer_centres(i2,iccb(i2))-p_layer_centres(i2,k)

        do kt=iccb(i2),k+1
          tempnum=tempnum+precip(i2,kt)
        end do
      else
        deltap_cld = 0.0
      end if

      ! Set logical to determine if downdraughts are allowed or not.
      if (deltap_cld >  15000.0 .and. bgmk(i2)                                 &
                              .and. tempnum >  1e-12) then
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
      dcpbydt(i)  = 0.0
      cape(i)     = 0.0
      bconv(i)    = .false.
      det_lev(i)  = k+1 ! Set final detrainment level (but not used).
    end do

  end if  ! nterm > 0

  !-----------------------------------------------------------------------
  ! Write out entrainment, detrainment and half-level mass flux diagnostics.
  ! They will be scaled by the full level mass flux outside
  ! of the level loop
  !-----------------------------------------------------------------------
  ! Calculate fractional entrainment rate for level k.
  if (flg_entr_up) then
    do i = 1,nconv
      entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i))                   &
               * (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i)             &
               * (1.0 + ekp14_c2(i)))
    end do
  end if

  ! Calculate fractional detrainment rate for level k
  if (flg_detr_up) then
    do i = 1,nconv
      detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)                         &
                      + deltak_c2(i) * (1.0 - amdetk_c2(i)))
    end do
  end if

  ! Calculate the half level mass flux for level k
  ! Only the scaling factor between full level and half levels is calculated
  ! here. This is scaled by the full level mass flux outside the level loop
  if (flg_up_flx_half) then
    do i = 1,nconv
      up_flux_half(index1(index2(i)),k) = (1.0 - deltak_c2(i))                 &
                * (1.0 - amdetk_c2(i)) * (1.0 + ekp14_c2(i))
    end do
  end if

  !-----------------------------------------------------------------------
  ! 3.6  End of main loop over levels
  !-----------------------------------------------------------------------
end do

!-----------------------------------------------------------------------
! Time smoothed mass flux
!-----------------------------------------------------------------------
if (l_conv_prog_flx) then
  decay_amount = timestep/tau_conv_prog_flx
  ! Update the smoothed mass flux
  do k = 1, nlev
    do i = 1, n_cg
      conv_prog_flx(i,k) = conv_prog_flx(i,k) + decay_amount * flx(i,k)
    end do
  end do
end if


!-----------------------------------------------------------------------
! Write out updraught massflux diagnostics and scale the
! entrainment and detrainment diagnostics by the mass flux.
!-----------------------------------------------------------------------
if (flg_up_flx) then
  do k = 1,nlev
    do i = 1,n_cg
      up_flux(i,k) = flx(i,k)
    end do
  end do
end if

if (flg_up_flx_half) then
  do k = 1,nlev
    do i = 1,n_cg
      up_flux_half(i,k) = up_flux_half(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_entr_up) then
  do k = 1,nlev
    do i = 1,n_cg
      entrain_up(i,k) = entrain_up(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_detr_up) then
  do k = 1,nlev
    do i = 1,n_cg
      detrain_up(i,k) = detrain_up(i,k) * flx(i,k)
    end do
  end do
end if

if (flg_area_ud) then
  do k = 1,nlev
    do i = 1,n_cg
      ! updraught core area from wup and flx
      if (w2p(i,k) > 0.0) then
        area_ud(i,k) =  flx(i,k)/(g*sqrt(w2p(i,k))* rho_theta(i,k))
      else
        area_ud(i,k)    = 0.0
      end if
    end do
  end do
end if
!-----------------------------------------------------------------------
! 4.0  Mixing of the convective increments in the boundary
!      layer.
!-----------------------------------------------------------------------

!      Mixes the increments from the initial parcel perturbation throughout
!      the subcloud layer if this option is selected.

call mix_ipert_6a (n_cg, nlev, nbl, n_wtrac_d, ntml,                           &
                p_layer_boundaries, exner_layer_centres,                       &
                dthbydt, dqbydt, wtrac_e, flx_init,                            &
                thpert, qpert, wtrac_p)

!-----------------------------------------------------------------------
! 5.0 All congestus convection will terminate at some level. This level
!     has been stored in the main level loop.
!
! Orignal downdraught calculation - on all points where convection is
! terminating.
!
!      Subroutine DD_ALL_CALL
!
!      UM Documentation Paper 27
!-----------------------------------------------------------------------

npossdd = 0
do i = 1,n_cg
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

  call dd_all_call_6a (n_cg,npossdd,kmax_term,nlev,trlev,ntra,n_wtrac_d        &
,                      kterm, iccb, icct, index_possdd, l_tracer               &
,                      bwater(1,2)                                             &
,                      exner_layer_centres,exner_layer_boundaries              &
,                      p_layer_centres, p_layer_boundaries,pstar               &
,                      recip_pstar,timestep , cca_2d                           &
,                      thp, qp, th, q, qse, trap,tracer, flx,precip            &
,                      dthbydt, dqbydt, dtrabydt                               &
,                      rain, snow, rain_3d, snow_3d                            &
,                      wtrac_p, wtrac_e                                        &
,                      dwn_flux, entrain_dwn, detrain_dwn                      &
,                      dt_dd, dq_dd)

end if

!-----------------------------------------------------------------------
! 5.2 Surface precipitation calculation for terminating points with
!     no downdraught (moved outside level loop) ie do this calculation
!     on all points at the end.
!-----------------------------------------------------------------------
! Points where no downdraught possible
nnodd = 0
do i = 1,n_cg

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


  call evap_bcb_nodd_all_6a (n_cg,nnodd,n_wtrac_d,kmax_term,kterm              &
                     , iccb,index_nodd,bwater(1,2)                             &
                     , exner_layer_centres                                     &
                     , p_layer_centres,p_layer_boundaries                      &
                     , timestep, cca_2d, th, q, qse, precip                    &
                     , dthbydt, dqbydt                                         &
                     , rain, snow, rain_3d, snow_3d                            &
                     , dt_dd, dq_dd, wtrac_p, wtrac_e)

end if

! Deallocate water tracer arrays - these are just dummy fields
call wtrac_dealloc_conv_p(1,wtrac_p)
call wtrac_dealloc_conv_e(1,wtrac_e)


!-----------------------------------------------------------------------
! 6.0  Convective Momentum Transport (if L_mom = .T.)
!-----------------------------------------------------------------------

if (l_mom) then


  ! Initialize arrays required for Convective Momentum Transport(CMT)

  k=1
  do i = 1,n_cg
    p_uv(k,i)     = p_layer_boundaries(i,k-1)
    phalf_uv(k,i) = p_layer_centres(i,k-1)
    ue_p(k,i)     = u(i,k)
    ve_p(k,i)     = v(i,k)
  end do

  do i = 1,n_cg
    nlcl_uv(i)    = ntml(i) + 1
    n_0degc(i)    = freeze_lev(i)
  end do

  do i = 1,n_cg
    do k = 2,nlev
      p_uv(k,i)     = p_layer_boundaries(i,k-1)
      phalf_uv(k,i) = p_layer_centres(i,k-1)
      ue_p(k,i)     = u(i,k)
      ve_p(k,i)     = v(i,k)
      exk_temp      = (p_uv(k,i)/pref)**kappa
      rho_uv(k,i)   = 2.0 * p_uv(k,i) / (r * exk_temp *                        &
                      (th(i,k-1) + th(i,k)))
    end do
    plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
    p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
    rho_uv(1,i)     = rho_uv(2,i)
  end do

  ! altered to use kterm instead of ntpar
  do i=1,n_cg
    if (kterm(i) >= nlcl_uv(i)) then
      ntop_uv(i) = kterm(i) +1
    else     ! case where congestus convection fails
      ! I think in this case the mass flux will be zero so no CMT
      ntop_uv(i) = ntpar(i) + 1
    end if
    ptop_uv(i) = phalf_uv(ntop_uv(i),i)
  end do

  ! Calculate CMT for required points

  nterm = 0

  do i = 1, n_cg
    ! Check that cloud base mass flux is non-zero.
    ! If convection failed then mb=0
    if (mb(i) > 0.0) then
      nterm = nterm + 1
      cu_term(nterm) = i
      cu_tend(nterm) = i
    end if
  end do

  ! Note using shallow CMT assumptions but may be using top from kterm
  ! May not be a sensible choice.

  if (nterm  >   0) then

    call shallow_grad_stress(n_cg,n_cg,nterm,nlev,cu_term,                     &
                             nlcl_uv,ntop_uv,mb,wsc,wstar,zcld,                &
                             plcl_uv,ptop_uv,p_uv,phalf_uv,                    &
                             rho_uv,ue_p,ve_p,timestep,                        &
                             weight_param,                                     &
                             ! in
                             uw,vw)

    call shallow_base_stress(n_cg,n_cg,n_cg,nlev,nterm,cu_term,                &
                             cu_tend,nlcl_uv,ntop_uv,mb,wsc,                   &
                             zlcl_uv,uw0,vw0,plcl_uv,                          &
                             ptop_uv,ue_p,ve_p,phalf_uv,p_uv,                  &
                             rho_uv,timestep,weight_param,                     &
                             ! INOUT
                             uw,vw,                                            &
                             ! out
                             uw_shall,vw_shall)


    call shallow_cmt_incr(n_cg,n_cg,n_cg,nlev,nterm,cu_term,                   &
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
  call cmt_heating(n_cg, nlev,                                                 &
                   z_theta, z_rho, exner_layer_centres,                        &
                   u, v, dubydt, dvbydt,                                       &
                   ! Out
                   dthbydt)
end if

!-----------------------------------------------------------------------
! 7.0  Energy (and optionally water) correction calculation
!-----------------------------------------------------------------------
do i = 1,n_cg
  index1(i) = i
end do

if (l_cv_conserve_check) then
  call cor_engy_6a(n_cg, n_cg, nlev, n_wtrac_d, index1, timestep,              &
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

call correct_small_q_conv(n_cg, n_cg, nlev, n_wtrac_d, index1, timestep, qmin, &
                          r2rho_th, dr_across_th,  q, 'water', 'Congest',      &
                          dqbydt, wtrac_e=wtrac_e)

!-----------------------------------------------------------------------
! 9.0  Calculate convective cloud amount on model levels - no anvils
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------------
! 9.1 CCRad - Calculate CCA fow shallow levels only
!-----------------------------------------------------------------------------

do i=1, n_cg

  overlap_fac(i) = 0.0

  if (iccb(i) /= 0) then ! Shallow convection occured

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

    overlap_fac(i) = 2.0*zcld(i) / z_rho(i,ntml(i)+1)

  end if     ! iccb
end do     ! n_cg

!---------------------------------------------------------------
! 9.11 Calculate CCA   - use shallow values
!---------------------------------------------------------------
select case (cca2d_sh_opt)
case (grant_lock)

  do i=1, n_cg
    if (iccb(i) /= 0) then ! Shallow convection occured
      cca_2d(i) = 2.0*mb(i)/wsc(i)

      ! Grab lowest cca value before any tuning occurs
      ! This will overwrite lcca in ni_conv_ctl only if neither
      ! shallow or deep have occured.  This is under a switch in
      ! the 4a scheme.
      !
      ! NOTE: Downdraughts & Evaporation still being fed cca_2d
      !       derived from TCW, This issue may require further
      !       investigation.
      lcca(i) = cca_2d(i)

    end if     ! iccb
  end do     ! n_cg

case (total_condensed_water)
    ! cca_2d is left unchanged from that calculated in the
    ! code, which is based on TCW (Total Condensed Water)
    ! (TCW is a rate)

end select


do i=1, n_cg

  overlap_fac(i) = max( 0.5, overlap_fac(i) )
  overlap_fac(i) = min( 5.0, overlap_fac(i) )

  if (overlap_fac(i)*cca_2d(i) > 0.99) then
    overlap_fac(i) = 0.99/cca_2d(i)
  end if
end do      ! i (n_cg)


!-------------------------------------------------------------------
! 9.12 Fill cca with cca_2d where non-zero ccw
!-------------------------------------------------------------------
do k=1, nlev
  do i=1, n_cg
    if (iccb(i) /= 0) then ! Shallow convection occured

      if (ccw(i,k) > 0.0) then
        zpr = (z_rho(i,k) - z_rho(i,ntml(i)+1)) / zcld(i)

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
  end do       ! i (n_cg)
end do       ! k (nlev)

!-----------------------------------------------------------------------
! 10.0  End Subroutine
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine congest_conv_6a
end module congest_conv_6a_mod
