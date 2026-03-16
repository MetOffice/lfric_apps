! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module glue_conv_6a_mod

use ereport_mod, only: ereport
use umPrintMgr,  only: umPrint, ummessage, printstatus, prstatus_normal
use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='GLUE_CONV_6A_MOD'

contains
!
subroutine glue_conv_6a(ncells,seg_size,nlev,n_wtrac,nbl,                      &
                     call_number,seg_num                                       &
   ,                 th,q,qcl,qcf                                              &
   ,                 q_wtrac, qcl_wtrac, qcf_wtrac                             &
   ,                 cf_liquid,cf_frozen,bulk_cf,pstar                         &
   ,                 bland,u,v,w                                               &
   ,                 tracer, dthbydt, dqbydt, dqclbydt, dqcfbydt               &
   ,                 dcflbydt, dcffbydt, dbcfbydt, dubydt, dvbydt              &
   ,                 dqbydt_wtrac, dqclbydt_wtrac, dqcfbydt_wtrac              &
   ,                 rain, snow, rain_3d, snow_3d, rain_wtrac, snow_wtrac      &
   ,                 cca0_dp, cca0_md, cca0_sh                                 &
   , cca0, iccb0, icct0, cclwp0, ccw0, lcbase0, lctop, lcca                    &
   , cca,  iccb,  icct,  cclwp,  ccw,  lcbase,  cca_2d, freeze_lev             &
   ,                 deep_cfl_limited,mid_cfl_limited                          &
   ,                 l_mid_all,kterm_deep,kterm_shall                          &
   ,                 precip_deep,precip_shall,precip_mid,precip_cong           &
   ,                 wstar_dn,wstar_up,mb1,mb2,kterm_congest                   &
   ,                 uw0,vw0,w_max,zlcl,zlcl_uv,tnuc_new,tnuc_nlcl             &
   ,                 ztop_uv                                                   &
   ,                 entrain_coef,conv_prog_precip,conv_prog_flx               &
   ,                 deep_flag,past_conv_ht,cape_out                           &
   ,                 n_dp,n_cg, n_sh, n_md                                     &
   ,                 r_rho,r_theta,rho, rho_theta, rho_dry                     &
   ,                 rho_dry_theta, delta_smag                                 &
   ,                 exner_rho_levels                                          &
   ,                 exner_layer_boundaries                                    &
   ,                 exner_layer_centres                                       &
   ,                 p_layer_boundaries                                        &
   ,                 p_layer_centres                                           &
   ,                 z_theta, z_rho                                            &
   ,                 timestep,t1_sd,q1_sd                                      &
   ,                 ntml,ntpar,conv_type,l_shallow_bl                         &
   ,                 l_congestus, l_mid                                        &
   ,                 cumulus_bl,wstar,wthvs,delthvu_bl,ql_ad                   &
   ,                 qsat_lcl, ftl, fqt                                        &
   ,                 l_tracer, ntra, trlev, n_cca_lev                          &
   ,                 l_calc_dxek, l_q_interact                                 &
   ,                 up_flux_half, up_flux, dwn_flux                           &
   ,                 entrain_up,detrain_up,entrain_dwn,detrain_dwn             &
   ,                 uw_deep,vw_deep,uw_shall,vw_shall,uw_mid,vw_mid           &
   ,                 wqt_flux_sh,wthetal_flux_sh                               &
   ,                 wthetav_flux_sh,wql_flux_sh                               &
   ,                 mf_deep,mf_congest,mf_shall,mf_midlev                     &
   ,                 dt_deep,dt_congest,dt_shall,dt_midlev                     &
   ,                 dq_deep,dq_congest,dq_shall,dq_midlev                     &
   ,                 du_deep,du_congest,du_shall,du_midlev                     &
   ,                 dv_deep,dv_congest,dv_shall,dv_midlev                     &
   ,                 ind_cape_reduced, cape_ts_used, ind_deep                  &
   ,                 ind_shall,w2p, dt_dd, dq_dd, du_dd, dv_dd                 &
   ,                 area_ud, area_dd                                          &
   ,                 scm_convss_dg, l_scm_convss_dg                            &
   ,                 g_ccp, h_ccp, ccp_strength                                &
   )
  ! Purpose:
  !  Gather-scatter routine for deep and shallow convection points.
  !  Interface to deep, shallow and mid level convection
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !  Language: F90
  !  This code is written to UMDP3 v6 programming standards
  !
use scm_convss_dg_mod,     only: scm_convss_dg_type

use cv_run_mod,  only:                                                         &
   l_mom, iconv_shallow, iconv_congestus, iconv_mid,                           &
   iconv_deep,                                                                 &
   cca_sh_knob,  cca_md_knob,   cca_dp_knob,                                   &
   ccw_sh_knob,  ccw_md_knob,   ccw_dp_knob,                                   &
   l_conv_hist,                                                                &
   l_conv_prog_flx, l_conv_prog_precip,                                        &
   midtrig_opt

use cv_stash_flg_mod, only:                                                    &
   flg_up_flx, flg_up_flx_half, flg_dwn_flx,                                   &
   flg_entr_up, flg_entr_dwn, flg_detr_up, flg_detr_dwn,                       &
   flg_uw_dp, flg_vw_dp, flg_uw_shall, flg_vw_shall,                           &
   flg_uw_mid, flg_vw_mid, flg_w_eqn,                                          &
   flg_wqt_flux, flg_wthetal_flux,                                             &
   flg_wthetav_flux, flg_wql_flux,                                             &
   flg_mf_deep, flg_mf_congest, flg_mf_shall, flg_mf_midlev,                   &
   flg_dt_deep, flg_dt_congest, flg_dt_shall, flg_dt_midlev,                   &
   flg_dq_deep, flg_dq_congest, flg_dq_shall, flg_dq_midlev,                   &
   flg_du_deep, flg_du_congest, flg_du_shall, flg_du_midlev,                   &
   flg_dv_deep, flg_dv_congest, flg_dv_shall, flg_dv_midlev,                   &
   flg_area_ud

use cv_param_mod, only:                                                        &
   mtrig_ntmlplus2, mtrig_ntmlplus1, mtrig_ntml, mtrig_surface


use cv_hist_constants_mod, only:                                               &
   decay_period

use tcs_warm_mod, only:                                                        &
   tcs_warm


use water_constants_mod, only: tm
use shallow_conv_6a_mod, only: shallow_conv_6a
use congest_conv_6a_mod, only: congest_conv_6a
use deep_conv_6a_mod,    only: deep_conv_6a
use mid_conv_6a_mod,     only: mid_conv_6a
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim
use rad_input_mod,       only: l_cca_dp_prog, l_cca_md_prog,                   &
                               l_cca_sh_prog

use model_domain_mod,      only: model_type, mt_single_column, mt_lfric

use wtrac_conv_mod,         only: l_wtrac_conv, conv_e_wtrac_type,             &
                                   wtrac_alloc_conv_e, wtrac_dealloc_conv_e
use wtrac_gather_conv_mod,  only: wtrac_gather_conv
use wtrac_scatter_conv_mod, only: wtrac_scatter_conv, wtrac_scatter_conv_mid

! Initialisation routines
use init_conv6a_mod,       only: init_6a_rain, init_6a_w2p, init_6a_npnts,     &
                                 init_6a_moist

use errormessagelength_mod, only: errormessagelength


use qsat_mod, only: qsat, qsat_mix

use gen_phys_inputs_mod, only: l_mr_physics
use mphys_inputs_mod,    only: l_progn_tnuc
implicit none
! Model Constants required


!------------------------------------------------------------------
! Subroutine Arguments
!------------------------------------------------------------------
!
! Arguments with intent in:
!
integer, intent(in) :: ncells ! No of points in a full field
!(note all multi-dimensional fields
!  passed to routine MUST be
!   dimensioned with ncells)

integer, intent(in) :: seg_size    ! No. of points in segment

! NOTE - loops over points should be over seg_size or less NEVER over ncells.

integer, intent(in) :: nlev        ! No. of model layers used in convection
integer, intent(in) :: n_wtrac     ! No. of water tracers,
                                   !   set to 1 if l_wtrac_conv=F
integer, intent(in) :: ntml(seg_size) ! Top level of surface mixed
! layer defined relative to
! theta,q grid

integer, intent(in) :: ntpar(seg_size) ! Top level of initial parcel
! ascent in BL scheme defined
! relative to theta,q grid

integer, intent(in) :: conv_type(seg_size)
! Integer index describing convective type:
!    0=no convection
!    1=non-precipitating shallow
!    2=drizzling shallow
!    3=warm congestus
!    ...
!    8=deep convection

integer, intent(in) :: n_cca_lev! No. of convective cloud
! amount levels (1 for 2D,
! nlevs for 3D)

integer, intent(in) :: nbl      ! No. of boundary layer levels

integer, intent(in) :: call_number ! Current sweep of convection

integer, intent(in) :: seg_num  ! Segment number

integer, intent(in) :: ntra     ! No. of tracer fields

integer, intent(in) :: trlev    ! No. of model levels on which
! tracers are included

integer, intent(in) :: n_dp     ! Number of deep points

integer, intent(in) :: n_cg     ! Number of congestus points

integer, intent(in) :: n_sh     ! Number of shallow points

integer, intent(in) :: n_md     ! Number of mid points

real(kind=real_umphys), intent(in) :: delthvu_bl(seg_size)
                                      !Integral of undilute parcel
! buoyancy over convective cloud
! layer (Kelvin m)

real(kind=real_umphys), intent(in) ::                                          &
    ql_ad(seg_size)      & ! adiabatic liquid water content at inversion (kg/kg)
   , qsat_lcl(seg_size)  & ! qsat at cloud base (kg/kg)
   , ftl(seg_size)       & ! Surface sensible heat flux divided by cp from BL
                        ! (K kg/m2/s) i.e. rho*w'tl'
   , fqt(seg_size)       & ! Total water flux from surface (kg/m2/s)
                        ! i.e. rho*w'qT'
   , delta_smag(seg_size)  ! grid size used in smagorinsky length scale (m)


real(kind=real_umphys), intent(in) :: r_rho(ncells,nlev)
                                            ! radius rho levels(m)
real(kind=real_umphys), intent(in) :: r_theta(ncells,0:nlev)
                                              ! theta levels (m)

real(kind=real_umphys), intent(in) ::                                          &
  rho(ncells,nlev)            & ! Wet density on rho levels (kg/m3)
 ,rho_theta(ncells,nlev)      & ! Wet denisty on theta levels (kg/m3)
 ,rho_dry(ncells,nlev)        & ! dry density on rho levels (kg/m3)
 ,rho_dry_theta(ncells,nlev)    ! dry density on theta levels (kg/m3)

real(kind=real_umphys), intent(in) :: exner_rho_levels(ncells,nlev)
                                                      !Exner on rho levels
real(kind=real_umphys), intent(in) :: exner_layer_centres(ncells,0:nlev)
                                                         !Exner

real(kind=real_umphys), intent(in) :: exner_layer_boundaries(ncells,0:nlev)
                                                            !Exner
! at half level above
! exner_layer_centres

real(kind=real_umphys), intent(in) :: pstar(seg_size) ! Surface pressure (Pa)

real(kind=real_umphys), intent(in) :: p_layer_centres(ncells,0:nlev)
                                                     !Pressure(Pa)


real(kind=real_umphys), intent(in) :: p_layer_boundaries(ncells,0:nlev)
                                                        ! Pressure
! at half level above
! p_layer_centres (Pa)

real(kind=real_umphys), intent(in) :: z_theta(ncells,nlev) ! height of theta
! levels above surface (m)

real(kind=real_umphys), intent(in) :: z_rho(ncells,nlev)
                                            ! height of rho levels
! above surface (m)


real(kind=real_umphys), intent(in) :: t1_sd(seg_size) ! Standard deviation of
! turbulent fluctuations of
! layer 1 temp. (K)

real(kind=real_umphys), intent(in) :: q1_sd(seg_size) ! Standard deviation of
! turbulent fluctuations of
! layer 1 q (kg/kg)

real(kind=real_umphys), intent(in) :: th(ncells,nlev) !Model potential
! temperature (K)

real(kind=real_umphys), intent(in) :: q(ncells,nlev)
                                     ! Model water vapour (kg/kg)

real(kind=real_umphys), intent(in) :: timestep    ! Model timestep (s)

real(kind=real_umphys), intent(in) :: q_wtrac(ncells,nlev,n_wtrac)
                                                  ! Water tracer vapour (kg/kg)
real(kind=real_umphys), intent(in) :: qcl_wtrac(ncells,nlev,n_wtrac)
                                                  ! Water tracer liquid
                                                  ! condensate (kg/kg)
real(kind=real_umphys), intent(in) :: qcf_wtrac(ncells,nlev,n_wtrac)
                                                  ! Water tracer ice
                                                  ! condensate (kg/kg)

real(kind=real_umphys), intent(in) :: uw0(seg_size)
                               ! U-comp of surface stress(N/m2)

real(kind=real_umphys), intent(in) :: vw0(seg_size)
                               ! V-comp of surface stress(N/m2)

real(kind=real_umphys), intent(in) :: w_max(seg_size) ! max w in column

real(kind=real_umphys), intent(in) :: u(ncells,nlev) !Model U field (m/s)

real(kind=real_umphys), intent(in) :: v(ncells,nlev) !Model V field (m/s)

real(kind=real_umphys), intent(in) :: w(ncells,nlev) !Model W field (m/s)

real(kind=real_umphys), intent(in) :: wstar(seg_size) ! Convective velocity scale
! (m/s)

real(kind=real_umphys), intent(in) :: wthvs(seg_size)
                                 ! Surface flux of THV  (Pa m/s2)

real(kind=real_umphys), intent(in) :: zlcl(seg_size)
                                ! Lifting condensation level accurate
! height (m) not a model level.

real(kind=real_umphys), intent(in) :: zlcl_uv(seg_size)
                                   ! Lifting condensation level
! defined for the uv grid (m)
real(kind=real_umphys), intent(in) :: tnuc_new(ncells,nlev)
                                            ! 3D ice nucleation temperature
                                            ! as function of dust
real(kind=real_umphys), intent(in) :: tnuc_nlcl(seg_size)
                                       !nucleation temperature as function
                                       !of dust indexed using nlcl
real(kind=real_umphys), intent(in) :: ztop_uv(seg_size) ! Top of cloud layer
! defined for the uv
! grid (m)

real(kind=real_umphys), intent(in) ::                                          &
   entrain_coef(seg_size)     ! entrainment coefficients

real(kind=real_umphys), intent(in) ::                                          &
  g_ccp(seg_size)         & ! cold-pool reduced gravity
 ,h_ccp(seg_size)         & ! cold-pool depth
 ,ccp_strength(seg_size)    ! cold-pool strength

real(kind=real_umphys), intent(in) :: conv_prog_precip(ncells,nlev)
                                                    ! Surface precipitation
                                                    ! based 3d convective
                                                    ! prognostic in kg/m2/s
real(kind=real_umphys), intent(in out) :: conv_prog_flx(ncells,nlev)
                                                     ! Mass flux convective
                                                    ! prognostic in Pa/s
! History prognostics only in use if l_conv_hist = .true.

real(kind=real_umphys), intent(in out) ::                                      &
   deep_flag(seg_size)     & ! 0.0-1.0, 1. indicates deep last time step
   ,past_conv_ht(seg_size)   ! convective height (m)

logical, intent(in) :: l_tracer ! Switch for inclusion of tracers

logical, intent(in) :: l_shallow_bl(seg_size) ! Shallow cumulus
! indicator

logical, intent(in) :: l_congestus(seg_size) ! congestus cumulus

logical, intent(in) :: l_mid(seg_size)       ! possible mid-level cnv

logical, intent(in) :: cumulus_bl(seg_size) ! Cumulus indicator

logical, intent(in) :: l_calc_dxek ! Switch for calculation of
! condensate increment

logical, intent(in) :: l_q_interact ! Switch allows overwriting
! parcel variables when
! calculating condensate incr.

logical, intent(in) :: bland(seg_size) ! Land/sea mask

!
! Arguments with intent INOUT:
!

! NOTE - All moist variables passed down to this routine are :
! specific humidities if l_mr_physics = .false.
! mixing ratios       if l_mr_physics = .true.

real(kind=real_umphys), intent(in out) :: qcl(ncells,nlev)
                                           ! Liq condensate (kg/kg)

real(kind=real_umphys), intent(in out) :: qcf(ncells,nlev)
                                           ! Ice condensate (kg/kg)

real(kind=real_umphys), intent(in out) :: cf_liquid(ncells,nlev)
! Liq water cloud volume (fraction)

real(kind=real_umphys), intent(in out) :: cf_frozen(ncells,nlev)
! Frozen water cloud volume (fraction?)

real(kind=real_umphys), intent(in out) :: bulk_cf(ncells,nlev)
                                               ! Bulk total cloud
                                              ! volume ( )

real(kind=real_umphys), intent(in out) :: tracer(ncells,trlev,ntra)
                                                    !Model tracer
                                              ! fields (kg/kg)

real(kind=real_umphys), intent(in out) ::                                      &
  w2p(ncells,nlev)        ! (Parcel vertical velocity)^2 [(m/s)^2]

!
! Arguments with intent out:
!

real(kind=real_umphys), intent(out) :: dqclbydt(ncells,nlev)
                                             ! Increments to liq
! condensate due to convection
! (kg/kg/s)

real(kind=real_umphys), intent(out) :: dqcfbydt(ncells,nlev)
                                             ! Increments to ice
! condensate due to convection
! (kg/kg/s)

real(kind=real_umphys), intent(out) :: dcflbydt(ncells,nlev)
                                             ! Increments to liq
! cloud volume due to convection
! (/s)

real(kind=real_umphys), intent(out) :: dcffbydt(ncells,nlev)
                                             ! Increments to ice
! cloud volume due to convection
! (/s)

real(kind=real_umphys), intent(out) :: dbcfbydt(ncells,nlev) ! Increments to
! total cld volume due to
! convection(/s)

real(kind=real_umphys), intent(out) :: dthbydt(ncells,nlev) ! Increments to
! potential temp. due to convection (K/s)

real(kind=real_umphys), intent(out) :: dqbydt(ncells,nlev)
                                           ! Increments to q due
! to convection (kg/kg/s)

real(kind=real_umphys), intent(out) :: dubydt(ncells,nlev+1)
                                             ! Increments to U due
! to CMT (m/s2)

real(kind=real_umphys), intent(out) :: dvbydt(ncells,nlev+1)
                                             ! Increments to V due
! to CMT (m/s2)

real(kind=real_umphys), intent(out) :: dqbydt_wtrac(ncells,nlev,n_wtrac)
                                            ! Increments to q_wtrac due
                                            ! to convection (kg/kg/s)
real(kind=real_umphys), intent(out) :: dqclbydt_wtrac(ncells,nlev,n_wtrac)
                                            ! Increments to qcl_wtrac
                                            ! due to convection (kg/kg/s)
real(kind=real_umphys), intent(out) :: dqcfbydt_wtrac(ncells,nlev,n_wtrac)
                                            ! Increments to qcf_wtrac
                                            ! due to convection (kg/kg/s)

real(kind=real_umphys), intent(out) :: rain(seg_size) ! Surface convective rainfall
! (kg/m2/s)

real(kind=real_umphys), intent(out) :: snow(seg_size) ! Surface convective snowfall
! (kg/m2/s)

real(kind=real_umphys), intent(in out) :: rain_3d(ncells,nlev)
                                            ! convective rainfall flux
! (kg/m2/s)

real(kind=real_umphys), intent(in out) :: snow_3d(ncells,nlev)
                                            ! convective snowfall flux
! (kg/m2/s)

real(kind=real_umphys), intent(out) :: rain_wtrac(ncells,n_wtrac)
                                            ! Surface convective water
                                            ! tracer rainfall (kg/m2/s)

real(kind=real_umphys), intent(out) :: snow_wtrac(ncells,n_wtrac)
                                            ! Surface convective water
                                            ! tracer snowfall (kg/m2/s)


! Section 5 Convective Cloud properties
real(kind=real_umphys), intent(out) ::                                         &
  cca(ncells,n_cca_lev)  &! Cnv. cld amount (0-1)
, ccw(ncells,nlev)       &! Cnv. in-cld liquid water (kg/kg)
, cclwp(seg_size)             &! Condensed water path (kg/m^2)
, lcca(seg_size)               ! Lowest cnv. cld amt. (0-1)

integer, intent(out) ::                                                        &
  iccb(seg_size)              &! Cnv. cld base level
, icct(seg_size)              &! Cnv. cld top level
, lcbase(seg_size)            &! Lowest cnv. cld base level
, lctop(seg_size)              ! Lowest cnv. cld top level


! Section 0 Convective Cloud properties for Radiative impacts
real(kind=real_umphys) ::                                                      &
  cca0(ncells,n_cca_lev)    &! Cnv. cld amount (0-1)
, cca0_sh(ncells,n_cca_lev) &! Shallow cnv. cld amount (0-1)
, cca0_md(ncells,n_cca_lev) &! Mid-level cnv. cld amount (0-1)
, cca0_dp(ncells,n_cca_lev) &! Deep cnv. cld amount (0-1)
, ccw0(ncells,nlev)         &! Cnv. in-cld liquid water (kg/kg)
, cclwp0(seg_size)                ! Cond. water path (kg/m^2)

integer, intent(out) ::                                                        &
  iccb0(seg_size)             &! Cnv. cld base level
, icct0(seg_size)             &! Cnv. cld top level
, lcbase0(seg_size)            ! Lowest cnv. cld base level


integer, intent(out) :: freeze_lev(seg_size) !index for freezing lev

integer, intent(out) :: kterm_deep(seg_size) ! index deep conv
integer, intent(out) :: kterm_shall(seg_size) ! level for shallow termination

logical, intent(out) :: l_mid_all(seg_size)  ! on exit true if mid level
                                          ! convection has triggered

real(kind=real_umphys), intent(out) ::                                         &
  deep_cfl_limited(seg_size) & !  indicator for cfl limited deep conv
 ,mid_cfl_limited(seg_size)    !  indicator for cfl limited mid conv

real(kind=real_umphys), intent(out) :: precip_deep(seg_size)
                                        ! deep precip (kg/m2/s)

real(kind=real_umphys), intent(out) :: precip_shall(seg_size)
                                         ! shallow precip(kg/m2/s)

real(kind=real_umphys), intent(out) :: precip_mid(seg_size) ! mid precip (kg/m2/s)
real(kind=real_umphys), intent(out) :: precip_cong(seg_size)
                                        ! congest precip (kg/m2/s)

real(kind=real_umphys), intent(out) ::                                         &
   wstar_dn(seg_size)        & ! subcloud layer convective velocity scale(m/s)
 , wstar_up(seg_size)        & ! cumulus layer convective velocity scale (m/s)
 , mb1(seg_size)             & ! cloud base mass flux from wstar_dn (m/s)
 , mb2(seg_size)               ! cloud base mass flux for cloud layer (m/s)

integer, intent(out) :: kterm_congest(seg_size) ! termination level
!                                                      for congestus


!
! Meaning of this diagnostic depends on scheme
!  Plume model - updraught mass flux  (Pa/s)
!  turbulence model - mass flux (not exactly updraught) (m/s)

real(kind=real_umphys), intent(out) :: up_flux(ncells,nlev) ! mass flux

real(kind=real_umphys), intent(out) :: up_flux_half(ncells,nlev)
                                                 !mass flux on rho
!dummy variable not used in turbulence

!
! Diagnostics with no meaning for turbulence based schemes
!
real(kind=real_umphys), intent(out) ::                                         &
  dwn_flux(ncells,nlev)    & ! Downdraught mass flux (Pa/s)

 ,entrain_up(ncells,nlev)  & ! Fractional entrainment rate into
                               ! updraughts (Pa/s)
 ,detrain_up(ncells,nlev)  & ! Fractional detrainment rate into
                               ! updraughts (Pa/s)
 ,entrain_dwn(ncells,nlev) & ! Fractional entrainment rate into
                               ! downdraughts (Pa/s)
 ,detrain_dwn(ncells,nlev)   ! Fractional detrainment rate into
                               ! downdraughts (Pa/s)
!
! Diagnostics relating to momentum fluxes
!
real(kind=real_umphys), intent(out) ::                                         &
   uw_deep(ncells,nlev)  & ! X-comp. of stress from deep convection
                             !(kg/m/s2)
 , vw_deep(ncells,nlev)  & ! Y-comp. of stress from deep convection
                             !(kg/m/s2)
 , uw_shall(ncells,nlev) & ! X-comp. of stress from shallow
                             ! convection (kg/m/s2)
 , vw_shall(ncells,nlev) & ! Y-comp. of stress from shallow
                             ! convection (kg/m/s2)
 , uw_mid(ncells,nlev)   & ! U comp of stress from mid convection (kg/m/s2)
 , vw_mid(ncells,nlev)     ! V comp of stress from mid convection (kg/m/s2)

real(kind=real_umphys), intent(out) :: cape_out(seg_size)
                                     ! Saved convective available
! potential energy for diagnostic
! output (Jkg-1)

! Fluxes from turbulence based convection schemes

real(kind=real_umphys), intent(out) ::                                         &
   wqt_flux_sh(ncells,nlev)      & ! w'qt' flux (m/s kg/kg)
 , wthetal_flux_sh(ncells,nlev)  & ! w'thetal' flux  (m/s K)
 , wthetav_flux_sh(ncells,nlev)  & ! w'thetav' flux  (m/s K)
 , wql_flux_sh(ncells,nlev)      & ! w'ql' flux  (m/s kg/kg)

 , mf_deep(ncells,nlev)          & ! mass flux deep
 , mf_congest(ncells,nlev)       & ! mass flux congestus
 , mf_shall(ncells,nlev)         & ! mass flux shallow
 , mf_midlev(ncells,nlev)        & ! mass flux mid-lev

 , dt_deep(ncells,nlev)          & ! dt increment deep   (K/s)
 , dt_congest(ncells,nlev)       & ! dt increment congestus (K/s)
 , dt_shall(ncells,nlev)         & ! dt increment shallow (K/s)
 , dt_midlev(ncells,nlev)        & ! dt increment mid-level (K/s)

 , dq_deep(ncells,nlev)          & ! dq increment deep (kg/kg/s)
 , dq_congest(ncells,nlev)       & ! dq increment congestus (kg/kg/s)
 , dq_shall(ncells,nlev)         & ! dq increment shallow (kg/kg/s)
 , dq_midlev(ncells,nlev)        & ! dq increment mid-level (kg/kg/s)

 , du_deep(ncells,nlev+1)        & ! du increment deep (m/s)
 , du_congest(ncells,nlev+1)     & ! du increment congestus (m/s)
 , du_shall(ncells,nlev+1)       & ! du increment shallow (m/s)
 , du_midlev(ncells,nlev+1)      & ! du increment mid-level (m/s)

 , dv_deep(ncells,nlev+1)        & ! dv increment deep (m/s)
 , dv_congest(ncells,nlev+1)     & ! dv increment congestus (m/s)
 , dv_shall(ncells,nlev+1)       & ! dv increment shallow (m/s)
 , dv_midlev(ncells,nlev+1)        ! dv increment mid-level (m/s)


real(kind=real_umphys), intent(out) ::                                         &
   ind_cape_reduced(seg_size)   & ! indicates cape timescale reduced
   ,cape_ts_used(seg_size)      & ! cape timescale for deep (s)
   ,ind_deep(seg_size)          & ! indicator of real deep
   ,ind_shall(seg_size)           ! indicator of real shallow

real(kind=real_umphys), intent(out) ::                                         &
   dt_dd(ncells,nlev)      & ! dT/dt from DD & evap below cloud base (K/s)
  ,dq_dd(ncells,nlev)      & ! dq/dt from DD & evap below cloud base (kg/kg/s)
  ,du_dd(ncells,nlev)      & ! du/dt from DD (m/s2)
  ,dv_dd(ncells,nlev)      & ! dv/dt from DD (m/s2)
  ,area_ud(ncells,nlev)    & ! updraught fractional area
  ,area_dd(ncells,nlev)      ! downdraught fractional area

! Structure containing SCM convection sub-step diagnostics
! (needs intent inout as contains allocatable arrays that need to
! retain their allocated status on input as well as output)
type(scm_convss_dg_type), intent(in out) :: scm_convss_dg( seg_size )
! Flag for SCM convection sub-step diagnostics
logical, intent(in) :: l_scm_convss_dg


! local required in move mixing to glue

real(kind=real_umphys) :: dtrabydt(seg_size,nlev,ntra) ! Increment to tracer due to
! convection (kg/kg/s)

! Structure containing compressed water fields (used for deep, shallow and
! mid-level convection)
type(conv_e_wtrac_type) :: wtrac_e(n_wtrac)

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to DEEP convection scheme.
! Arrays are identified by underscore DP (_dp) and are of length
! n_dp where n_dp is the number of points diagnosed as deep in the
! boundary layer diagnosis routine. For full desriptions of variables
! see above.
!-----------------------------------------------------------------------
integer :: error_point      ! location of problem deep point

integer :: dpi(n_dp)        ! index for deep points in full grid

integer ::                                                                     &
   ntml_dp(n_dp)                                                               &
   , ntpar_dp(n_dp)

real(kind=real_umphys) ::                                                      &
   pstar_dp(n_dp)                                                              &
   , recip_pstar_dp(n_dp)                                                      &
   , t1_sd_dp(n_dp)                                                            &
   , q1_sd_dp(n_dp)                                                            &
   , uw0_dp(n_dp)                                                              &
   , vw0_dp(n_dp)                                                              &
   , zlcl_uv_dp(n_dp)                                                          &
   , tnuc_nlcl_dp(n_dp)                                                        &
   !similar to tnuc_nlcl. To be passed to deep convection
   , wstar_dp(n_dp)                                                            &
   , delthvu_dp(n_dp)                                                          &
   , entrain_coef_dp(n_dp)                                                     &
   , qsat_lcl_dp(n_dp)

! cold-pool variables
real(kind=real_umphys) ::                                                      &
   ccp_strength_dp(n_dp)                                                       &
   , g_ccp_dp(n_dp)                                                            &
   , h_ccp_dp(n_dp)

! Added for deep turbulence scheme

real(kind=real_umphys) ::                                                      &
   p_layer_centres_dp(n_dp,0:nlev)                                             &
   , p_layer_boundaries_dp(n_dp,0:nlev)                                        &
   , exner_rho_levels_dp(n_dp,nlev)                                            &
   , exner_layer_centres_dp(n_dp,0:nlev)                                       &
   , exner_layer_boundaries_dp(n_dp,0:nlev)                                    &
   , r2rho_dp(n_dp,nlev)                                                       &
   , r2rho_th_dp(n_dp,nlev)                                                    &
   , rho_dp(n_dp,nlev)                                                         &
   , rho_theta_dp(n_dp,nlev)                                                   &
   , dr_across_th_dp(n_dp,nlev)                                                &
   , dr_across_rh_dp(n_dp,nlev)                                                &
   , z_theta_dp(n_dp,nlev)                                                     &
   , z_rho_dp(n_dp,nlev)                                                       &
   , r_rho_dp(n_dp,nlev)                                                       &
   , r_theta_dp(n_dp,0:nlev)

! Surface precipitation based 3d convective prognostic in kg/m2/s
real(kind=real_umphys) :: conv_prog_precip_dp(n_dp,nlev)
! Mass flux convective prognostic in Pa/s
real(kind=real_umphys) :: conv_prog_flx_dp(n_dp,nlev)


logical :: bland_dp(n_dp)

real(kind=real_umphys) ::                                                      &
   u_dp(n_dp,nlev)                                                             &
   , v_dp(n_dp,nlev)                                                           &
   , w_dp(n_dp,nlev)                                                           &
   , th_dp(n_dp,nlev)                                                          &
   , q_dp(n_dp,nlev)                                                           &
   , qse_dp(n_dp,nlev)

real(kind=real_umphys) :: tracer_dp(n_dp,trlev,ntra)
! increments
real(kind=real_umphys) ::                                                      &
   dthbydt_dp(n_dp,nlev)                                                       &
   , dqbydt_dp(n_dp,nlev)                                                      &
   , dubydt_dp(n_dp,nlev)                                                      &
   , dvbydt_dp(n_dp,nlev)

! output variables
real(kind=real_umphys) ::                                                      &
   rain_dp(n_dp)                                                               &
   , snow_dp(n_dp)                                                             &
   , rain_3d_dp(n_dp,nlev)                                                     &
   , snow_3d_dp(n_dp,nlev)                                                     &
   , tcw_dp(n_dp)                                                              &
   , cclwp_dp(n_dp)                                                            &
   , lcca_dp(n_dp)                                                             &
   , cape_out_dp(n_dp)                                                         &
   , mb_dp(n_dp)

integer ::                                                                     &
   iccb_dp(n_dp)                                                               &
   , icct_dp(n_dp)                                                             &
   , lcbase_dp(n_dp)                                                           &
   , lctop_dp(n_dp)                                                            &
   , freeze_lev_dp(n_dp)

real(kind=real_umphys) ::                                                      &
   ccw_dp(n_dp,nlev)                                                           &
   , up_flux_dp(n_dp,nlev)                                                     &
   , up_flux_half_dp(n_dp,nlev)                                                &
   , dwn_flux_dp(n_dp,nlev)                                                    &
   , entrain_up_dp(n_dp,nlev)                                                  &
   , detrain_up_dp(n_dp,nlev)                                                  &
   , entrain_dwn_dp(n_dp,nlev)                                                 &
   , detrain_dwn_dp(n_dp,nlev)                                                 &
   , uw_deep_dp(n_dp,nlev)                                                     &
   , vw_deep_dp(n_dp,nlev)                                                     &
   , w_max_dp(n_dp)

!PC2
real(kind=real_umphys) ::                                                      &
   qcl_dp(n_dp,nlev)                                                           &
   , qcf_dp(n_dp,nlev)                                                         &
   , cf_liquid_dp(n_dp,nlev)                                                   &
   , cf_frozen_dp(n_dp,nlev)                                                   &
   , dqclbydt_dp(n_dp,nlev)                                                    &
   , dqcfbydt_dp(n_dp,nlev)                                                    &
   , dcflbydt_dp(n_dp,nlev)                                                    &
   , dcffbydt_dp(n_dp,nlev)                                                    &
   , dbcfbydt_dp(n_dp,nlev)

real(kind=real_umphys) :: dtrabydt_dp(n_dp,nlev,ntra)
                                    ! Increment to tracer due to
                                    ! convection (kg/kg/s)
integer :: kterm_dp(n_dp)    ! required by mid level scheme

real(kind=real_umphys) ::                                                      &
   cca_2d_dp(n_dp)             & ! required by mid level scheme
   , cca_dp(n_dp,n_cca_lev)    & ! 3d CCA
   , ind_cape_reduced_dp(n_dp) & ! indicates reduced cape timescale
   , cape_ts_used_dp(n_dp)     & ! cape timescale deep
   , cfl_limited_dp(n_dp)      & ! CFL limited conv
   , ind_deep_dp(n_dp)           ! indicator of real deep

! Temporary arrays contianing copies of fields for deep, shallow, congestus
! or mid-level convection

real(kind=real_umphys), allocatable ::                                         &
  dt_dd_temp(:,:)           & ! dt from DD & evap below cloud from scheme
 ,dq_dd_temp(:,:)           & ! dq from DD & evap below cloud from scheme
 ,area_ud_temp(:,:)           ! updraught fractional area

real(kind=real_umphys), allocatable ::                                         &
  rho_dry_temp(:,:)         & ! dry density on rho levels (kg/m3)
 ,rho_dry_theta_temp(:,:)     ! dry density on theta levels (kg/m3)

! Version of the structure containing SCM convection sub-step diagnostics,
! to be compressed for passing into deep, mid or shallow
type(scm_convss_dg_type), allocatable :: scm_convss_dg_c(:)


!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to congestus convection scheme.
! Arrays are identified by underscore cg (_cg) and are of length
! n_cg where n_cg is the number of points diagnosed as congestus in the
! boundary layer diagnosis routine.For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

integer :: cgi(n_cg)     ! index for congestus points in full grid

integer ::                                                                     &
   ntml_cg(n_cg)                                                               &
   , ntpar_cg(n_cg)

real(kind=real_umphys) ::                                                      &
   pstar_cg(n_cg)                                                              &
   , recip_pstar_cg(n_cg)                                                      &
   , delthvu_cg(n_cg)                                                          &
   , uw0_cg(n_cg)                                                              &
   , vw0_cg(n_cg)                                                              &
   , wstar_cg(n_cg)                                                            &
   , wthvs_cg(n_cg)                                                            &
   , zlcl_uv_cg(n_cg)                                                          &
   , tnuc_nlcl_cg(n_cg)                                                        &
   !similar to tnuc_nlcl. To be passed to congestus convection
   , ztop_uv_cg(n_cg)                                                          &
   , entrain_coef_cg(n_cg)

real(kind=real_umphys) :: ccp_strength_cg(n_cg) ! cold-pool strength

real(kind=real_umphys) ::                                                      &
   p_layer_centres_cg(n_cg,0:nlev)                                             &
   , p_layer_boundaries_cg(n_cg,0:nlev)                                        &
   , exner_rho_levels_cg(n_cg,nlev)                                            &
   , exner_layer_centres_cg(n_cg,0:nlev)                                       &
   , exner_layer_boundaries_cg(n_cg,0:nlev)                                    &
   , z_theta_cg(n_cg,nlev)                                                     &
   , z_rho_cg(n_cg,nlev)                                                       &
   , u_cg(n_cg,nlev)                                                           &
   , v_cg(n_cg,nlev)                                                           &
   , th_cg(n_cg,nlev)                                                          &
   , q_cg(n_cg,nlev)                                                           &
   , ccw_cg(n_cg,nlev)                                                         &
   , qse_cg(n_cg,nlev)                                                         &
   , r_rho_cg(n_cg,nlev),r_theta_cg(n_cg,0:nlev)                               &
   , r2rho_th_cg(n_cg,nlev)                                                    &
                            ! radius**2 density theta lev (kg/m)
   , r2rho_cg(n_cg,nlev)                                                       &
                            ! radius**2 density rho lev (kg/m)
   , dr_across_th_cg(n_cg,nlev)                                                &
                            ! thickness of theta levels (m)
   , dr_across_rh_cg(n_cg,nlev)                                                &
                            ! thickness of rho levels (m)
   , rho_theta_cg(n_cg,nlev)        ! rho on theta levels (kg/m3)

! Surface precipitation based 3d convective prognostic in kg/m2/s
real(kind=real_umphys) :: conv_prog_precip_cg(n_cg,nlev)
! Mass flux convective prognostic in Pa/s
real(kind=real_umphys) :: conv_prog_flx_cg(n_cg,nlev)


logical :: bland_cg(n_cg)

! tracers
real(kind=real_umphys) :: tracer_cg(n_cg,trlev,ntra)                           &
   , dtrabydt_cg(n_cg,nlev,ntra)

real(kind=real_umphys) ::                                                      &
   dthbydt_cg(n_cg,nlev)                                                       &
   , dqbydt_cg(n_cg,nlev)                                                      &
   , dubydt_cg(n_cg,nlev)                                                      &
   , dvbydt_cg(n_cg,nlev)

real(kind=real_umphys) ::                                                      &
   rain_cg(n_cg)                                                               &
   , snow_cg(n_cg)                                                             &
   , rain_3d_cg(n_cg,nlev)                                                     &
   , snow_3d_cg(n_cg,nlev)                                                     &
   , tcw_cg(n_cg)                                                              &
   , cclwp_cg(n_cg)                                                            &
   , lcca_cg(n_cg)                                                             &
   , cape_out_cg(n_cg)

real(kind=real_umphys) ::                                                      &
   cca_2d_cg(n_cg)         & ! required by mid level scheme
   , cca_cg(n_cg,n_cca_lev)  ! 3d CCA

integer ::                                                                     &
   iccb_cg(n_cg)                                                               &
   , icct_cg(n_cg)                                                             &
   , lcbase_cg(n_cg)                                                           &
   , lctop_cg(n_cg)                                                            &
   , freeze_lev_cg(n_cg)                                                       &
   , kterm_cg(n_cg)    ! required by mid level scheme

! diagnostics
real(kind=real_umphys) ::                                                      &
   up_flux_cg(n_cg,nlev)                                                       &
   , up_flux_half_cg(n_cg,nlev)                                                &
   , dwn_flux_cg(n_cg,nlev)                                                    &
   , entrain_up_cg(n_cg,nlev)                                                  &
   , detrain_up_cg(n_cg,nlev)                                                  &
   , entrain_dwn_cg(n_cg,nlev)                                                 &
   , detrain_dwn_cg(n_cg,nlev)                                                 &
   , uw_shall_cg(n_cg,nlev)                                                    &
   , vw_shall_cg(n_cg,nlev)

!PC2
real(kind=real_umphys) ::                                                      &
   qcl_cg(n_cg,nlev)                                                           &
   , qcf_cg(n_cg,nlev)                                                         &
   , cf_liquid_cg(n_cg,nlev)                                                   &
   , cf_frozen_cg(n_cg,nlev)                                                   &
   , bulk_cf_cg(n_cg,nlev)                                                     &
   , dqclbydt_cg(n_cg,nlev)                                                    &
   , dqcfbydt_cg(n_cg,nlev)                                                    &
   , dcflbydt_cg(n_cg,nlev)                                                    &
   , dcffbydt_cg(n_cg,nlev)                                                    &
   , dbcfbydt_cg(n_cg,nlev)


!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to SHALLOW convection scheme.
! Arrays are identified by underscore SH (_sh) and are of length
! n_sh where n_sh is the number of points diagnosed as shallow in the
! boundary layer diagnosis routine.For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

integer :: shi(n_sh)      ! index for shallow points in full grid

! compressed inputs

integer :: ntml_sh(n_sh)                                                       &
   , ntpar_sh(n_sh)                                                            &
   , kterm_sh(n_sh)  ! May be needed in adaptive option 4 or 6

integer ::                                                                     &
   conv_type_sh(n_sh)
! Integer index describing convective type:
!    1=non-precipitating shallow
!    2=drizzling shallow
!    3=warm congestus
! This array is defined on points diagnosed as
! shallow.


real(kind=real_umphys) ::                                                      &
   pstar_sh(n_sh)                                                              &
   , recip_pstar_sh(n_sh)                                                      &
   , delthvu_sh(n_sh)                                                          &
   , ql_ad_sh(n_sh)                                                            &
   , uw0_sh(n_sh)                                                              &
   , vw0_sh(n_sh)                                                              &
   , wstar_sh(n_sh)                                                            &
   , wthvs_sh(n_sh)                                                            &
   , zlcl_sh(n_sh)                                                             &
   , zlcl_uv_sh(n_sh)                                                          &
   , tnuc_nlcl_sh(n_sh)                                                        &
   !similar to tnuc_nlcl. To be passed to shallow convection
   , ztop_uv_sh(n_sh)                                                          &
   , entrain_coef_sh(n_sh)                                                     &
   , delta_smag_sh(n_sh)                                                       &
   , ind_shall_sh(n_sh)          ! indicator of real shallow

real(kind=real_umphys) :: ccp_strength_sh(n_sh) ! cold-pool strength

real(kind=real_umphys) ::                                                      &
   p_layer_centres_sh(n_sh,0:nlev)                                             &
   , p_layer_boundaries_sh(n_sh,0:nlev)                                        &
   , exner_rho_levels_sh(n_sh,nlev)                                            &
   , exner_layer_centres_sh(n_sh,0:nlev)                                       &
   , exner_layer_boundaries_sh(n_sh,0:nlev)                                    &
   , z_theta_sh(n_sh,nlev), z_rho_sh(n_sh,nlev)                                &
   , rho_sh(n_sh,nlev)                                                         &
   , r_rho_sh(n_sh,nlev),r_theta_sh(n_sh,0:nlev)                               &
   , r2rho_th_sh(n_sh,nlev)                                                    &
                            ! radius**2 density theta lev (kg/m)
   , r2rho_sh(n_sh,nlev)                                                       &
                            ! radius**2 density rho lev (kg/m)
   , dr_across_th_sh(n_sh,nlev)                                                &
                            ! thickness of theta levels (m)
   , dr_across_rh_sh(n_sh,nlev)                                                &
                            ! thickness of rho levels (m)
   , rho_theta_sh(n_sh,nlev)        ! rho on theta levels (kg/m3)

! Surface precipitation based 3d convective prognostic in kg/m2/s
real(kind=real_umphys) :: conv_prog_precip_sh(n_sh,nlev)
! Mass flux convective prognostic in Pa/s
real(kind=real_umphys) :: conv_prog_flx_sh(n_sh,nlev)


logical :: bland_sh(n_sh)

real(kind=real_umphys) ::                                                      &
   u_sh(n_sh,nlev)                                                             &
   , v_sh(n_sh,nlev)                                                           &
   , w_sh(n_sh,nlev)                                                           &
   , th_sh(n_sh,nlev)                                                          &
   , q_sh(n_sh,nlev)                                                           &
   , qse_sh(n_sh,nlev)

! output
real(kind=real_umphys) ::                                                      &
   dthbydt_sh(n_sh,nlev)                                                       &
   , dqbydt_sh(n_sh,nlev)                                                      &
   , dubydt_sh(n_sh,nlev)                                                      &
   , dvbydt_sh(n_sh,nlev)

real(kind=real_umphys) ::                                                      &
   rain_sh(n_sh)                                                               &
   , snow_sh(n_sh)                                                             &
   , rain_3d_sh(n_sh,nlev)                                                     &
   , snow_3d_sh(n_sh,nlev)                                                     &
   , tcw_sh(n_sh)                                                              &
   , cclwp_sh(n_sh)                                                            &
   , lcca_sh(n_sh)                                                             &
   , cape_out_sh(n_sh)                                                         &
   , cca_2d_sh(n_sh)                                                           &
   , cca_sh(n_sh,n_cca_lev) ! required by mid level scheme
                            ! 3d CCA

!split these out for lfric ifdef
real(kind=real_umphys) ::                                                      &
   wstar_up_sh(n_sh)                                                           &
   , mb1_sh(n_sh)                                                              &
   , mb2_sh(n_sh)


integer :: iccb_sh(n_sh)                                                       &
   , icct_sh(n_sh)                                                             &
   , lcbase_sh(n_sh)                                                           &
   , lctop_sh(n_sh)                                                            &
   , freeze_lev_sh(n_sh)

real(kind=real_umphys) ::                                                      &
   ccw_sh(n_sh,nlev)                                                           &
   , up_flux_sh(n_sh,nlev)                                                     &
   , up_flux_half_sh(n_sh,nlev)                                                &
   , dwn_flux_sh(n_sh,nlev)                                                    &
   , entrain_up_sh(n_sh,nlev)                                                  &
   , detrain_up_sh(n_sh,nlev)                                                  &
   , entrain_dwn_sh(n_sh,nlev)                                                 &
   , detrain_dwn_sh(n_sh,nlev)                                                 &
   , uw_shall_sh(n_sh,nlev)                                                    &
   , vw_shall_sh(n_sh,nlev)

! PC2  input
real(kind=real_umphys) ::                                                      &
   qcl_sh(n_sh,nlev)                                                           &
   , qcf_sh(n_sh,nlev)                                                         &
   , cf_liquid_sh(n_sh,nlev)                                                   &
   , cf_frozen_sh(n_sh,nlev)
! PC2 output
real(kind=real_umphys) ::                                                      &
   dqclbydt_sh(n_sh,nlev)                                                      &
   , dqcfbydt_sh(n_sh,nlev)                                                    &
   , dcflbydt_sh(n_sh,nlev)                                                    &
   , dcffbydt_sh(n_sh,nlev)                                                    &
   , dbcfbydt_sh(n_sh,nlev)

! Tracers in and out
real(kind=real_umphys) :: tracer_sh(n_sh,trlev,ntra)
real(kind=real_umphys) :: dtrabydt_sh(n_sh,nlev,ntra)
                                    ! Increment to tracer due to
                                    ! convection (kg/kg/s)

! Local arrays for (Parcel vertical velocity)^2 [(m/s)^2] on ...
real(kind=real_umphys) :: w2p_sh(n_sh,nlev) ! ...shallow   convection points
real(kind=real_umphys) :: w2p_cg(n_cg,nlev) ! ...congestus convection points
real(kind=real_umphys) :: w2p_md(n_md,nlev) ! ...mid-level convection points
real(kind=real_umphys) :: w2p_dp(n_dp,nlev) ! ...deep      convection points

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to MID-LEVEL convection scheme.
! Arrays are identified by underscore MD (_md) and are of length
! seg_size where seg_size is the total number of points in the grid (since
! mid-level convection may occur on any point previously diagnosed as
! shallow or deep).
!-----------------------------------------------------------------------


integer :: midtrig(n_md)   ! Level at which mid level convection
                           ! may start
integer :: mdi(n_md)       ! index for mid points in full grid

integer ::                                                                     &
   ntml_md(n_md)                                                               &
   , ntpar_md(n_md)

real(kind=real_umphys) ::                                                      &
   pstar_md(n_md)                                                              &
   , recip_pstar_md(n_md)

! cold-pool variables
real(kind=real_umphys) ::                                                      &
   ccp_strength_md(n_md)                                                       &
   , g_ccp_md(n_md)                                                            &
   , h_ccp_md(n_md)


logical :: bland_md(n_md)
logical :: l_mid_md(n_md)     ! true if mid level convection
                              ! compressed version of l_mid_all
real(kind=real_umphys) ::  ftl_md(n_md)
                              ! Compressed surface sensible heat flux/cp
real(kind=real_umphys) ::  fqt_md(n_md)
                              ! Compressed surface total water flux

real(kind=real_umphys) ::                                                      &
   cca_2d_md(n_md)         &! required by mid level scheme
   ,cca_md(n_md,n_cca_lev)  ! 3d CCA

real(kind=real_umphys) ::                                                      &
   p_layer_centres_md(n_md,0:nlev)                                             &
   , p_layer_boundaries_md(n_md,0:nlev)                                        &
   , exner_rho_levels_md(n_md,nlev)                                            &
   , exner_layer_centres_md(n_md,0:nlev)                                       &
   , exner_layer_boundaries_md(n_md,0:nlev)                                    &
   , z_theta_md(n_md,nlev), z_rho_md(n_md,nlev)                                &
   , r_theta_md(n_md,0:nlev), r_rho_md(n_md,nlev)                              &
   , r2rho_th_md(n_md,nlev)                                                    &
   , r2rho_md(n_md,nlev)                                                       &
   , rho_md(n_md,nlev)                                                         &
   , rho_theta_md(n_md,nlev)                                                   &
   , dr_across_th_md(n_md,nlev)                                                &
   , dr_across_rh_md(n_md,nlev)                                                &
   , u_md(n_md,nlev)                                                           &
   , v_md(n_md,nlev)                                                           &
   , w_md(n_md,nlev)                                                           &
   , th_md(n_md,nlev)                                                          &
   , q_md(n_md,nlev)

! Surface precipitation based 3d convective prognostic in kg/m2/s
real(kind=real_umphys) :: conv_prog_precip_md(n_md,nlev)
! Mass flux convective prognostic in Pa/s
real(kind=real_umphys) :: conv_prog_flx_md(n_md,nlev)

! output variables
real(kind=real_umphys) ::                                                      &
   dthbydt_md(n_md,nlev)                                                       &
   , dqbydt_md(n_md,nlev)                                                      &
   , dubydt_md(n_md,nlev)                                                      &
   , dvbydt_md(n_md,nlev)

real(kind=real_umphys) :: rain_md(n_md)                                        &
   , snow_md(n_md)                                                             &
   , rain_3d_md(n_md,nlev)                                                     &
   , snow_3d_md(n_md,nlev)                                                     &
   , uw_mid_md(n_md,nlev)                                                      &
   , vw_mid_md(n_md,nlev)                                                      &
   , cclwp_md(n_md)                                                            &
   , ccw_md(n_md,nlev)                                                         &
   , lcca_md(n_md)                                                             &
   , cape_out_md(n_md)                                                         &
   , cfl_limited_md(n_md)

integer ::                                                                     &
   iccb_md(n_md)                                                               &
   , icct_md(n_md)                                                             &
   , lcbase_md(n_md)                                                           &
   , lctop_md(n_md)                                                            &
   , freeze_lev_md(n_md)

real(kind=real_umphys) ::                                                      &
   up_flux_md(n_md,nlev)                                                       &
   , up_flux_half_md(n_md,nlev)                                                &
   , dwn_flux_md(n_md,nlev)                                                    &
   , entrain_up_md(n_md,nlev)                                                  &
   , detrain_up_md(n_md,nlev)                                                  &
   , entrain_dwn_md(n_md,nlev)                                                 &
   , detrain_dwn_md(n_md,nlev)                                                 &
   , w_max_md(n_md)

! PC2
real(kind=real_umphys) ::                                                      &
   qcl_md(n_md,nlev)                                                           &
   , qcf_md(n_md,nlev)                                                         &
   , cf_liquid_md(n_md,nlev)                                                   &
   , cf_frozen_md(n_md,nlev)                                                   &
   , qse_md(n_md,nlev)

real(kind=real_umphys) ::                                                      &
   dqclbydt_md(n_md,nlev)                                                      &
   , dqcfbydt_md(n_md,nlev)                                                    &
   , dcflbydt_md(n_md,nlev)                                                    &
   , dcffbydt_md(n_md,nlev)                                                    &
   , dbcfbydt_md(n_md,nlev)

! tracers
real(kind=real_umphys) :: tracer_md(n_md,trlev,ntra)
real(kind=real_umphys) :: dtrabydt_md(n_md,nlev,ntra)
                                    ! Increment to tracer due to
                                    ! convection (kg/kg/s)

real(kind=real_umphys) :: tnuc_new_md(n_md,nlev)
! ice nucleation temperature to be passed to mid convection
integer :: errorstatus       ! error status

! added for deep flag
logical ::                                                                     &
    l_deep           ! indicates deep convection

real(kind=real_umphys) ::                                                      &
    decay_amount  &  ! fraction to reduce deep flag by
   ,cld_depth        ! cloud depth

! should be in a include file
real(kind=real_umphys), parameter :: deep_depth=2500.0
                                         ! Minimum depth to be classed
                                          ! as a deep cloud

!-----------------------------------------------------------------------
! Mixing ratio / specific humidity conversions
!-----------------------------------------------------------------------
!  The full model holds q, qcl and qcf (high resolution versions can also
!  hold other moist variables.) The current mass flux scheme without
!  PC2 ignores the presence of qcl and qcf.
!  This assumption is also true of the new shallow turbulence based scheme.
!
!  Let m indicate mixing ratio and q specific then
!
!  Now taking account of all moist variables in the dynamics, we use start
!  of timestep fields (at time t) and calculate mt(t) from:
!      1 + mt(t) = rho_wet/rho_dry
!
!  Not PC2
!      dqcl = 0.0,      dqcf = 0.0
!      dmcl = 0.0,      dmcf = 0.0
!
! In the following, the left column is for scheme's working in m but
! inputs are q (l_mr_physics=F); the right is for scheme's working in q
! but inputs are m (l_mr_physics=T)
!
! Before scheme
!  mv(t) = qv(t)*(1+mt(t))          or qv(t) = mv(t)/(1+mt(t))
!
! How to convert increments
!  dqv = dmv/(1+mt(t))              or dmv = dqv*(1+mt(t))
!
! where
!  dmv= mv(t+dt) -  mv(t)           or dqv = qv(t+dt)- qv(t)
!
!
! Conversion of qsat
! ------------------
! rsat - mixing ratio saturation value
!
! qsat=rsat/(1+rsat)      and rsat=qsat(1-qsat)
!
!-----------------------------------------------------------------------
!  Arrays required for mixing ratio to specific conversions

real(kind=real_umphys) :: mt(seg_size,nlev)    ! total water mixing ratio (kg/kg)
                                            ! rhowet/rhodry = 1 + mt
real(kind=real_umphys) :: denom             ! 1/denominator

!-----------------------------------------------------------------------
!  Arrays required by 3d CCA calculation

real(kind=real_umphys), intent(in out) :: cca_2d(seg_size)
                                             ! 2d convective cloud (Section 5)

! Saturation mixing ratio calulation and 1/pstar

real(kind=real_umphys) ::                                                      &
  qse_mix(seg_size,nlev)     ! Saturation mixing ratio of cloud
                          ! cloud environment (kg/kg) only used if
                          ! turbulence schemes in use
real(kind=real_umphys) ::                                                      &
   recip_pstar(seg_size)                                                          &
                            ! Reciprocal of pstar array
   , qse(seg_size,nlev)                                                           &
                            ! Saturation specific humidity of cloud
                            ! cloud environment (kg/kg)
   , pt(seg_size)                                                                 &
                            ! Temporary store for P in calc. of sat.
                            ! value. (Pa)
   , tt(seg_size)                                                                 &
                            ! Temporary store for T in calc.
                            ! of saturation value. (K)
   , ttkm1(seg_size)           ! Temporary store for T in layer
                            ! k-1 for use in freezing level
                            ! calc. for anvil. (K)

!  arrays required by grid calculations

real(kind=real_umphys) ::                                                      &
   r2rho_th(seg_size,nlev)      & ! radius**2 density theta lev (kg/m)
   ,r2rho(seg_size,nlev)         & ! radius**2 density rho lev (kg/m)
   ,dr_across_th(seg_size,nlev)  & ! thickness of theta levels (m)
   ,dr_across_rh(seg_size,nlev)    ! thickness of rho levels (m)

!  arrays required by for identifying convective points

integer :: index1(seg_size)
integer :: nconv_all         ! Total number of points convecting

! array required by tracers at end

real(kind=real_umphys) :: limited_step(seg_size)     ! Reduced step size for tracer
! mixing

real(kind=real_umphys) :: step_test(seg_size)
                                ! Work array used in reducing step

!
! Parameters
!

real(kind=real_umphys), parameter :: safety_margin = 1.0e-100
                                            ! Small number used in
! tracer step reduction

!---------------------------------------------------------------------
!
! Loop counters
!

integer :: i,j,k,ktra,i_wt

integer :: ierror                 ! error status

character(len=errormessagelength) :: cmessage    ! error message
character (len=*), parameter ::  RoutineName = 'GLUE_CONV_6A'


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialise output variables not initialised already
!-----------------------------------------------------------------------
!
! Initialise 3d rain variables

call init_6a_rain(seg_size, nlev, n_sh, n_dp, n_cg, n_md, rain_3d_sh, snow_3d_sh, &
                  rain_3d_dp, snow_3d_dp, rain_3d_cg, snow_3d_cg, rain_3d_md,  &
                  snow_3d_md, ind_cape_reduced, cape_ts_used, ind_deep,        &
                  ind_shall, kterm_deep, kterm_shall, kterm_congest )

!==========================================
! Initialise local arrays for w-eqn
!
if (flg_w_eqn) then

  call init_6a_w2p(nlev, n_sh, n_dp, n_cg, n_md, w2p_sh, w2p_cg, w2p_md,       &
                   w2p_dp )

end if
!
!==========================================

 ! Initialise various convection variables, outputs and diagnostics to zero
 ! at all points. Points not in the compression lists for the respective calls
 ! to DPCONV, SHCONV, MDCONV or CGCONV will remain zero.

call init_6a_npnts(seg_size, nlev, n_cca_lev, ncells, n_wtrac,                  &
                   dqbydt, dthbydt, dqclbydt, dqcfbydt, dcflbydt, dcffbydt,    &
                   dbcfbydt, dqbydt_wtrac, dqclbydt_wtrac,                     &
                   dqcfbydt_wtrac, ccw, cca,                                   &
                   dt_dd, dq_dd, du_dd, dv_dd, area_ud, area_dd, up_flux,      &
                   up_flux_half, dwn_flux,                                     &
                   entrain_up, detrain_up, entrain_dwn, detrain_dwn,           &
                   mf_deep, mf_congest,                                        &
                   mf_shall, mf_midlev, dt_deep, dt_congest, dt_shall,         &
                   dt_midlev, dq_deep, dq_congest, dq_shall, dq_midlev,        &
                   du_deep, du_congest, du_shall, du_midlev, dv_deep,          &
                   dv_congest, dv_shall, dv_midlev, wqt_flux_sh,               &
                   wql_flux_sh, wthetal_flux_sh, wthetav_flux_sh,              &
                   uw_deep, vw_deep, uw_shall, vw_shall, uw_mid, vw_mid )


! Initialisation of mt used for mixing/specific conversions using wet/dry rho.
! Also zero dubydt, dvbydt and dtrabydt arrays.

call init_6a_moist(seg_size, nlev, ncells, ntra, rho_theta, rho_dry_theta,      &
                   mt, dubydt, dvbydt, dtrabydt, l_tracer)
!
! Required to get same convective cloud as old scheme
! Need to pass values from deep and shallow to mid level.
! Currently conversion of 2d to 3d makes use of total 2d from
! shallow/deep and mid in a column. If in future the conversion
! does not work on a column total this could be removed.
!

!-----------------------------------------------------------------------
!  grid information
!-----------------------------------------------------------------------
! Calculate quantities involving density etc for future use

do k=1,nlev
  do i=1,seg_size
    dr_across_rh(i,k) = r_theta(i,k) - r_theta(i,k-1)
    r2rho(i,k)        = r_rho(i,k)*r_rho(i,k)*rho(i,k)
  end do
end do
!
! rho_theta only set for nlev-1
!
k=1     ! bottom theta level thicker
do i=1,seg_size
  dr_across_th(i,k) = r_rho(i,k+1) - r_theta(i,0)
  r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k)
end do

do k=2,nlev-1
  do i=1,seg_size
    dr_across_th(i,k) = r_rho(i,k+1) - r_rho(i,k)
    r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k)
  end do
end do

k=nlev     ! top layer  (hope not used ?
!             assume density as layer below)
do i=1,seg_size
  dr_across_th(i,k) = r_theta(i,nlev) - r_rho(i,k)
  r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k-1)
end do

!-----------------------------------------------------------------------
! 1.0 Section to calculate fields used by all convection types.
! Note cheaper to do here than in individual routines.
! Create saturation mixing ratio arrays  & calculate freeze level
!-----------------------------------------------------------------------

!
! Calculate 1/pstar and initialize freeze_lev array.
!

do i = 1,seg_size
  recip_pstar(i)=1.0 / pstar(i)
  freeze_lev(i) = 1
end do

!
! Loop over levels
!

do k = 1,nlev

  !
  ! Find freezing level
  !

  if (k  ==  1) then
    do i = 1,seg_size
      tt(i) = th(i,k) * exner_layer_centres(i,k)
      pt(i) = p_layer_centres(i,k)
      !
      ! Commented out as initialisation sets freeze_lev to 1.
      ! Code left incase altered in future.
      !
      !            If (tt(i)  <   TM) then
      !              freeze_lev(i) = k
      !            End If

    end do
  else
    do i = 1,seg_size
      ttkm1(i) = tt(i)
      tt(i) = th(i,k) * exner_layer_centres(i,k)
      pt(i) = p_layer_centres(i,k)
      if (tt(i)  <   tm .and. ttkm1(i)  >=  tm) then
        if (freeze_lev(i) == 1) then
          freeze_lev(i) = k
        end if
      end if
    end do
  end if

  !
  ! Calculate saturation specific humidity/mixing ratio  lq_mix=.false.
  !

  call qsat(qse(:,k),tt,pt,seg_size)

  if (iconv_deep == 2 .or. (iconv_shallow >= 2) ) then
    ! Using a turbulence scheme also require a mixing ratio version

    call qsat_mix(qse_mix(:,k),tt,pt,seg_size)

  end if
end do  ! nlev


! Setup compressed array structure for SCM convection sub-step
! diagnostics, if used

! Allocate the array of diagnostic structures.
! We will use the same array for deep, shallow and mid-level (just
! rezero the fields between each call).  Size of array needs to be
! max( n_dp, n_sh, n_mid ); then we only use elements 1:n_dp for
! deep etc.  This saves cost of repeatedly allocating / deallocating
! with different sizes for deep, shallow and mid.
! Note: this is allocated in all runs since the calls to convection
! routines below require the outer object to be indexed

if (max(n_dp, n_sh, n_md) > 0) then
  allocate( scm_convss_dg_c( max(n_dp, n_sh, n_md) ) )
end if


!-----------------------------------------------------------------------
! 1.0 DEEP Convection
! 1.1 Compress input variable arrays for deep convection scheme to
!     length n_dp (deep points only)
!-----------------------------------------------------------------------

if (n_dp  >   0 .and. iconv_deep >  0 ) then

  allocate(rho_dry_temp(n_dp,nlev))
  allocate(rho_dry_theta_temp(n_dp,nlev))

  j = 0
  do i = 1,seg_size
    if (cumulus_bl(i) .and. .not. l_shallow_bl(i) .and.                        &
       .not. l_congestus(i)) then
      j                        = j+1
      dpi(j)                   = i
    end if
  end do
  !
  ! In only variables
  !
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do j=1,n_dp
    bland_dp(j)           = bland(dpi(j))
    ntml_dp(j)            = ntml(dpi(j))
    ntpar_dp(j)           = ntpar(dpi(j))
    pstar_dp(j)           = pstar(dpi(j))
    recip_pstar_dp(j)     = recip_pstar(dpi(j))
    q1_sd_dp(j)           = q1_sd(dpi(j))
    t1_sd_dp(j)           = t1_sd(dpi(j))
    uw0_dp(j)             = uw0(dpi(j))
    vw0_dp(j)             = vw0(dpi(j))
    w_max_dp(j)           = w_max(dpi(j))
    wstar_dp(j)           = wstar(dpi(j))
    qsat_lcl_dp(j)        = qsat_lcl(dpi(j))
    zlcl_uv_dp(j)         = zlcl_uv(dpi(j))
    freeze_lev_dp(j)      = freeze_lev(dpi(j))
    delthvu_dp(j)         = delthvu_bl(dpi(j))
    entrain_coef_dp(j)    = entrain_coef(dpi(j))
    g_ccp_dp(j)           = g_ccp(dpi(j))
    h_ccp_dp(j)           = h_ccp(dpi(j))
    ccp_strength_dp(j)    = ccp_strength(dpi(j))
  end do

  do k = 0,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j=1,n_dp
      p_layer_centres_dp(j,k)    = p_layer_centres(dpi(j),k)
      p_layer_boundaries_dp(j,k) = p_layer_boundaries(dpi(j),k)

      exner_layer_centres_dp(j,k)    = exner_layer_centres(dpi(j),k)
      exner_layer_boundaries_dp(j,k) = exner_layer_boundaries(dpi(j),k)

      r_theta_dp(j,k)  = r_theta(dpi(j),k)
    end do
  end do

  do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j=1,n_dp
      u_dp(j,k)           = u(dpi(j),k)
      v_dp(j,k)           = v(dpi(j),k)
      w_dp(j,k)           = w(dpi(j),k)
      th_dp(j,k)          = th(dpi(j),k)
      exner_rho_levels_dp(j,k)  = exner_rho_levels(dpi(j),k)
      z_theta_dp(j,k)     = z_theta(dpi(j),k)
      z_rho_dp(j,k)       = z_rho(dpi(j),k)
      r_rho_dp(j,k)       = r_rho(dpi(j),k)
      r2rho_th_dp(j,k)     = r2rho_th(dpi(j),k)
      r2rho_dp(j,k)        = r2rho(dpi(j),k)
      rho_theta_dp(j,k)    = rho_theta(dpi(j),k)
      rho_dp(j,k)          = rho(dpi(j),k)
      rho_dry_theta_temp(j,k) = rho_dry_theta(dpi(j),k)
      rho_dry_temp(j,k)       = rho_dry(dpi(j),k)
      dr_across_th_dp(j,k) = dr_across_th(dpi(j),k)
      dr_across_rh_dp(j,k) = dr_across_rh(dpi(j),k)
    end do
  end do
  if (l_progn_tnuc) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j = 1, n_dp
      tnuc_nlcl_dp(j) = tnuc_nlcl(dpi(j))
    end do
  end if

  if (l_conv_prog_precip) then  ! Initialise precip based prognostic
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1, n_dp
        conv_prog_precip_dp(j,k) = conv_prog_precip(dpi(j),k)
      end do
    end do
  end if

  if (l_conv_prog_flx) then  ! Initialise mass flux based prognostic
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1, n_dp
        conv_prog_flx_dp(j,k) = conv_prog_flx(dpi(j),k)
      end do
    end do
  end if


  if (iconv_deep == 1) then  ! G-R mass flux scheme

    if (l_mr_physics) then  !  Conversion from m to q required

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_dp
          denom = 1.0 / ( 1.0 + mt(dpi(j),k) )
          q_dp(j,k)   = q(dpi(j),k)   * denom
          qse_dp(j,k) = qse(dpi(j),k) ! copy specific version
          qcl_dp(j,k) = qcl(dpi(j),k) * denom
          qcf_dp(j,k) = qcf(dpi(j),k) * denom
        end do
      end do

    else         ! input is specific humidity therefore no problems

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_dp
          q_dp(j,k)   = q(dpi(j),k)
          qse_dp(j,k) = qse(dpi(j),k) ! copy specific version
          qcl_dp(j,k) = qcl(dpi(j),k)
          qcf_dp(j,k) = qcf(dpi(j),k)
        end do
      end do

    end if      ! Test on l_mr_physics

  else        ! Future schemes to be code in mixing ratio

    if (l_mr_physics) then   ! Input as required

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_dp
          q_dp(j,k)   = q(dpi(j),k)
          qse_dp(j,k) = qse_mix(dpi(j),k) ! copy mixing ratio version
          qcl_dp(j,k) = qcl(dpi(j),k)
          qcf_dp(j,k) = qcf(dpi(j),k)
        end do
      end do

    else      !  Conversion from q to m required

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_dp
          q_dp(j,k)   = q(dpi(j),k)  *(1.0+mt(dpi(j),k))
          qcl_dp(j,k) = qcl(dpi(j),k)*(1.0+mt(dpi(j),k))
          qcf_dp(j,k) = qcf(dpi(j),k)*(1.0+mt(dpi(j),k))
          qse_dp(j,k) = qse_mix(dpi(j),k) ! copy mixing ratio version
        end do
      end do

    end if    ! Test on l_mr_physics

  end if   ! test on deep scheme

  ! Water tracers - allocate and set compressed dp arrays and convert between
  ! q and mixing ratio if required

  call wtrac_alloc_conv_e(n_dp, nlev, n_wtrac, wtrac_e)

  if (l_wtrac_conv) then

    call wtrac_gather_conv(ncells,seg_size,n_dp,nlev,n_wtrac, iconv_deep,       &
                           dpi, mt, q_wtrac, qcl_wtrac, qcf_wtrac, wtrac_e)

  end if

  !
  ! In/out variables
  !

  do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j=1,n_dp
      cf_liquid_dp(j,k)   = cf_liquid(dpi(j),k)
      cf_frozen_dp(j,k)   = cf_frozen(dpi(j),k)
    end do
  end do

  if (l_tracer) then
    do ktra = 1,ntra
      do k = 1,trlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j=1,n_dp
          tracer_dp(j,k,ktra)  = tracer(dpi(j),k,ktra)
        end do
      end do
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j=1,n_dp
          dtrabydt_dp(j,k,ktra)  = 0.0
        end do
      end do
    end do
  end if

  allocate(dt_dd_temp(n_dp,nlev))
  allocate(dq_dd_temp(n_dp,nlev))
  allocate(area_ud_temp(n_dp,nlev))


  !-----------------------------------------------------------------------
  ! 1.2 Call deep convection code
  !-----------------------------------------------------------------------

  if (iconv_deep == 1) then  ! 4A like Gregory Rowntree deep conv

    call deep_conv_6a(                                                         &
                            !in
       nbl,nlev,ntra,n_wtrac,n_cca_lev,n_dp,trlev,                             &
       bland_dp, delthvu_dp,                                                   &
       exner_rho_levels_dp,                                                    &
       exner_layer_centres_dp,                                                 &
       exner_layer_boundaries_dp,                                              &
       l_q_interact,                                                           &
       l_tracer, ntml_dp, ntpar_dp,                                            &
       pstar_dp,p_layer_centres_dp,                                            &
       p_layer_boundaries_dp,                                                  &
       z_theta_dp, z_rho_dp,                                                   &
       r_theta_dp, r_rho_dp,                                                   &
       rho_theta_dp, rho_dp,                                                   &
       rho_dry_theta_temp, rho_dry_temp,                                       &
       r2rho_th_dp, r2rho_dp,                                                  &
       dr_across_th_dp, dr_across_rh_dp,                                       &
       conv_prog_precip_dp,                                                    &
       q_dp,th_dp,timestep,                                                    &
       u_dp,v_dp,w_dp,                                                         &
       uw0_dp,vw0_dp,w_max_dp,wstar_dp,qsat_lcl_dp,                            &
       entrain_coef_dp,                                                        &
       zlcl_uv_dp,tnuc_nlcl_dp,freeze_lev_dp,                                  &
       recip_pstar_dp,qse_dp,                                                  &
       l_scm_convss_dg,                                                        &
       g_ccp_dp, h_ccp_dp, ccp_strength_dp,                                    &
                            !INOUT
       cf_frozen_dp,cf_liquid_dp,                                              &
       qcf_dp,qcl_dp,tracer_dp,wtrac_e,                                        &
       w2p_dp,conv_prog_flx_dp,scm_convss_dg_c(1:n_dp),                        &
                            !out
       cape_out_dp,cclwp_dp,ccw_dp,cca_dp,                                     &
       dbcfbydt_dp,dcffbydt_dp,dcflbydt_dp,                                    &
       dqbydt_dp,dqcfbydt_dp,dqclbydt_dp,                                      &
       dthbydt_dp,dubydt_dp,dvbydt_dp,dtrabydt_dp,                             &
       detrain_up_dp,detrain_dwn_dp,entrain_up_dp,                             &
       entrain_dwn_dp,                                                         &
       iccb_dp,icct_dp,                                                        &
       lcca_dp,lcbase_dp,lctop_dp,                                             &
       rain_dp,snow_dp, rain_3d_dp, snow_3d_dp,                                &
       up_flux_dp,up_flux_half_dp,                                             &
       dwn_flux_dp,uw_deep_dp,vw_deep_dp,kterm_dp,                             &
       tcw_dp,cca_2d_dp,                                                       &
       ind_cape_reduced_dp,                                                    &
       cape_ts_used_dp,cfl_limited_dp,ind_deep_dp,                             &
       dt_dd_temp, dq_dd_temp, area_ud_temp,                                   &
       error_point                                                             &
       )

    if (error_point /= 0) then
      errorstatus = 2   ! will cause model to fail
      write(cmessage,'(a37,i12,a8,i3,a9,i2)')                                  &
      'Deep conv went to model top at point ',                                 &
       dpi(error_point),' in seg ',seg_num,' on call ',call_number

      call ereport(routinename, errorstatus, cmessage )

    end if
  end if

  !-----------------------------------------------------------------------
  ! 1.3 Write data from deep convection points to full arrays
  !-----------------------------------------------------------------------
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i = 1,n_dp
    cape_out(dpi(i))           = cape_out_dp(i)
    cclwp(dpi(i))              = cclwp_dp(i)
    iccb(dpi(i))               = iccb_dp(i)
    icct(dpi(i))               = icct_dp(i)
    lcca(dpi(i))               = lcca_dp(i)
    lcbase(dpi(i))             = lcbase_dp(i)
    lctop(dpi(i))              = lctop_dp(i)

    rain(dpi(i))               = rain_dp(i)
    snow(dpi(i))               = snow_dp(i)
    precip_deep(dpi(i))        = rain_dp(i) + snow_dp(i)
    kterm_deep(dpi(i))         = kterm_dp(i)
    ind_cape_reduced(dpi(i))   = ind_cape_reduced_dp(i)
    cape_ts_used(dpi(i))       = cape_ts_used_dp(i)
    deep_cfl_limited(dpi(i))   = cfl_limited_dp(i)
    ind_deep(dpi(i))           = ind_deep_dp(i)
    cca_2d(dpi(i)) = cca_2d_dp(i)
  end do

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        rain_wtrac(dpi(i),i_wt) = wtrac_e(i_wt)%rain(i)
        snow_wtrac(dpi(i),i_wt) = wtrac_e(i_wt)%snow(i)
      end do
    end do
  end if

  if (l_mom) then
    do k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        dubydt(dpi(i),k)         = dubydt_dp(i,k)
        dvbydt(dpi(i),k)         = dvbydt_dp(i,k)
      end do
    end do
  end if

  if (flg_w_eqn) then
    do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_dp
        w2p(dpi(i),k) = w2p_dp(i,k)
      end do
    end do
  end if

  if (l_conv_prog_flx) then
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1, n_dp
        conv_prog_flx(dpi(i),k) = conv_prog_flx_dp(i,k)
      end do
    end do
  end if

  if (flg_area_ud) then
    do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_dp
        area_ud(dpi(i),k)        = area_ud_temp(i,k)
      end do
    end do
  end if

  do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_dp
      dthbydt(dpi(i),k)        = dthbydt_dp(i,k)
      dcflbydt(dpi(i),k)       = dcflbydt_dp(i,k)
      dcffbydt(dpi(i),k)       = dcffbydt_dp(i,k)
      dbcfbydt(dpi(i),k)       = dbcfbydt_dp(i,k)
      ccw(dpi(i),k)            = ccw_dp(i,k)
      dt_dd(dpi(i),k)          = dt_dd_temp(i,k)
    end do
  end do

  ! Water tracer (dqbydt_wtrac etc.) arrays are not decompressed here as this
  ! is done below in combination with any necessary q/mixing ratio conversions

  do k = 1,n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_dp
      cca(dpi(i),k) = cca_dp(i,k)
    end do
  end do

  ! Then scaling is done by CCRad, even if PC2 = .true.
  ! To zero CCA (section 0) for PC2 use the CCRAD knobs
  ! cca_knobs should be set to 0.0 in gui/namelist for PC2

  ! Use of Anvils will affect both section 0/5 diags if ccrad knobs
  ! are not set to 0.0, this may required further attention in future
  do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1, n_dp
      ccw0(dpi(i),k) = ccw_dp(i,k) * ccw_dp_knob
    end do
  end do

  do k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1, n_dp
      cca0(dpi(i),k) = cca_dp(i,k) * cca_dp_knob
    end do
  end do


!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_dp
    cclwp0 (dpi(i)) = cclwp  (dpi(i))
    iccb0  (dpi(i)) = iccb   (dpi(i))
    icct0  (dpi(i)) = icct   (dpi(i))
    lcbase0(dpi(i)) = lcbase (dpi(i))
  end do

  if (l_cca_dp_prog) then
    do k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_dp
        cca0_dp(dpi(i),k) = cca_dp(i,k)
      end do
    end do
  end if


  if (iconv_deep == 1) then  !G-R mass flux scheme
    if (l_mr_physics) then ! requires conversion from q incs to m incs

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_dp
          denom = 1.0 + mt(dpi(i),k)
          dqbydt(dpi(i),k)   = denom * dqbydt_dp(i,k)
          dqclbydt(dpi(i),k) = denom * dqclbydt_dp(i,k)
          dqcfbydt(dpi(i),k) = denom * dqcfbydt_dp(i,k)
          dq_dd(dpi(i),k)    = denom * dq_dd_temp(i,k)
        end do
      end do

    else   ! output is specific humidity therefore no problems

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_dp
          dqbydt(dpi(i),k)   = dqbydt_dp(i,k)
          dqclbydt(dpi(i),k) = dqclbydt_dp(i,k)
          dqcfbydt(dpi(i),k) = dqcfbydt_dp(i,k)
          dq_dd(dpi(i),k)    = dq_dd_temp(i,k)
        end do
      end do

    end if      ! Test on l_mr_physics

  else         ! Future schemes output mixing ratio incs

    if (l_mr_physics) then ! ok
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_dp
          dqbydt(dpi(i),k)   = dqbydt_dp(i,k)
          dqclbydt(dpi(i),k) = dqclbydt_dp(i,k)
          dqcfbydt(dpi(i),k) = dqcfbydt_dp(i,k)
          dq_dd(dpi(i),k)    = dq_dd_temp(i,k)
        end do
      end do

    else   ! output needs to be specific humidity

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_dp
          denom = 1.0 / ( 1.0 + mt(dpi(i),k) )
          dqbydt(dpi(i),k)   = denom * dqbydt_dp(i,k)
          dqclbydt(dpi(i),k) = denom * dqclbydt_dp(i,k)
          dqcfbydt(dpi(i),k) = denom * dqcfbydt_dp(i,k)
          dq_dd(dpi(i),k)    = denom * dq_dd_temp(i,k)
        end do
      end do

    end if      ! Test on l_mr_physics
  end if        ! test on deep scheme

  ! Convert water tracers between specific humidity and mixing ratios if
  ! required and decompress arrays
  if (l_wtrac_conv) then

    call wtrac_scatter_conv(ncells, seg_size, n_dp, nlev, n_wtrac,              &
                            iconv_deep, dpi, mt, wtrac_e,                      &
                            dqbydt_wtrac, dqclbydt_wtrac, dqcfbydt_wtrac)
  end if


  deallocate(area_ud_temp)
  deallocate(dq_dd_temp)
  deallocate(dt_dd_temp)
  deallocate(rho_dry_theta_temp)
  deallocate(rho_dry_temp)

  if (flg_up_flx) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        up_flux(dpi(i),k)        = up_flux_dp(i,k)
      end do
    end do
  end if
  if (flg_up_flx_half .and. iconv_deep == 1) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        up_flux_half(dpi(i),k)        = up_flux_half_dp(i,k)
      end do
    end do
  end if
  if (flg_dwn_flx) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        dwn_flux(dpi(i),k)       = dwn_flux_dp(i,k)
      end do
    end do
  end if
  if (flg_entr_up) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        entrain_up(dpi(i),k)     = entrain_up_dp(i,k)
      end do
    end do
  end if
  if (flg_detr_up) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        detrain_up(dpi(i),k)     = detrain_up_dp(i,k)
      end do
    end do
  end if
  if (flg_entr_dwn) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        entrain_dwn(dpi(i),k)    = entrain_dwn_dp(i,k)
      end do
    end do
  end if
  if (flg_detr_dwn) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        detrain_dwn(dpi(i),k)    = detrain_dwn_dp(i,k)
      end do
    end do
  end if
  if (flg_uw_dp) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        uw_deep(dpi(i),k)        = uw_deep_dp(i,k)
      end do
    end do
  end if
  if (flg_vw_dp) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        vw_deep(dpi(i),k)        = vw_deep_dp(i,k)
      end do
    end do
  end if
  if (flg_mf_deep) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        mf_deep(dpi(i),k)        = up_flux_dp(i,k)
      end do
    end do
  end if
  if (flg_dt_deep) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        dt_deep(dpi(i),k)        = dthbydt_dp(i,k)
      end do
    end do
  end if
  if (flg_dq_deep) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_dp
        dq_deep(dpi(i),k)        = dqbydt_dp(i,k)
      end do
    end do
  end if
  if (l_mom) then
    if (flg_du_deep) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_dp
          du_deep(dpi(i),k)        = dubydt_dp(i,k)
        end do
      end do
    end if
    if (flg_dv_deep) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_dp
          dv_deep(dpi(i),k)        = dvbydt_dp(i,k)
        end do
      end do
    end if
  end if

  if (l_tracer) then
    do ktra = 1,ntra
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1, n_dp
          dtrabydt(dpi(i),k,ktra)  = dtrabydt_dp(i,k,ktra)
        end do
      end do
    end do
  end if

  ! Merge dp 3d rain & snow profiles

  do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1, n_dp
      rain_3d(dpi(i),k) = rain_3d_dp(i,k)
      snow_3d(dpi(i),k) = snow_3d_dp(i,k)
    end do
  end do

  !----------------------------------------------------------------------
  ! Setting of deep convective history flag
  ! Only class as deep convection if depth of convection greater than
  ! a set value.
  !----------------------------------------------------------------------

  if (l_conv_hist) then

    decay_amount= timestep/decay_period

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_dp
      cld_depth = z_rho_dp(i,kterm_dp(i)+1) - z_rho_dp(i,ntml_dp(i)+1)

      ! Store depth of convection
      past_conv_ht(dpi(i)) = cld_depth

      if (cld_depth > deep_depth) then
        deep_flag(dpi(i)) = 1.0      ! deep convection occurred
      else   ! failed deep convection
        deep_flag(dpi(i)) = deep_flag(dpi(i)) - decay_amount
        deep_flag(dpi(i)) = max(deep_flag(dpi(i)),0.0)
      end if
    end do       ! loop over deep poins

  end if

  ! Deallocate water tracer deep convection arrays
  call wtrac_dealloc_conv_e(n_wtrac, wtrac_e)

end if ! test on n_dp >0


!-----------------------------------------------------------------------
! If not deep decay deep flag
!-----------------------------------------------------------------------

if (l_conv_hist) then

  decay_amount= timestep/decay_period

  do i = 1,seg_size
    l_deep=cumulus_bl(i) .and. .not. l_shallow_bl(i) .and.                     &
                                      .not. l_congestus(i)
    ! Decay flag ensuring 0.0 is the minimum allowed value.
    if (.not. l_deep) then
      deep_flag(i) = deep_flag(i) - decay_amount
      deep_flag(i) = max(deep_flag(i),0.0)

      ! Do I want to decay old convective depth here? At present NO
      !            past_conv_ht(i) = (1.0-decay_amount) * past_conv_ht(i)

    end if
  end do

end if    ! test on l_conv_hist

!-----------------------------------------------------------------------
! 2.0 SHALLOW convection
! 2.1 Compress input variable arrays for shallow convection scheme to
!     length n_sh (shallow points only)
!-----------------------------------------------------------------------

if (n_sh  >   0 .and. iconv_shallow >  0) then

  allocate(rho_dry_temp(n_sh,nlev))
  allocate(rho_dry_theta_temp(n_sh,nlev))

  j = 0
  do i = 1,seg_size
    if (cumulus_bl(i) .and. l_shallow_bl(i)) then
      j                        = j+1
      shi(j)                   = i
    end if
  end do
  !
  ! In only variables
  !
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do j = 1,n_sh
    bland_sh(j)           = bland(shi(j))
    delthvu_sh(j)         = delthvu_bl(shi(j))
    ql_ad_sh(j)           = ql_ad(shi(j))
    ntml_sh(j)            = ntml(shi(j))
    ntpar_sh(j)           = ntpar(shi(j))
    conv_type_sh(j)       = conv_type(shi(j))
    pstar_sh(j)           = pstar(shi(j))
    recip_pstar_sh(j)     = recip_pstar(shi(j))
    uw0_sh(j)             = uw0(shi(j))
    vw0_sh(j)             = vw0(shi(j))
    wstar_sh(j)           = wstar(shi(j))
    wthvs_sh(j)           = wthvs(shi(j))
    zlcl_uv_sh(j)         = zlcl_uv(shi(j))
    ztop_uv_sh(j)         = ztop_uv(shi(j))
    freeze_lev_sh(j)      = freeze_lev(shi(j))
    entrain_coef_sh(j)    = entrain_coef(shi(j))
    delta_smag_sh(j)      = delta_smag(shi(j))
    ccp_strength_sh(j)    = ccp_strength(shi(j))

    ! initialise to zero as not used in every option

    wstar_up_sh(j) = 0.0
    mb1_sh(j) = 0.0
    mb2_sh(j) = 0.0

  end do

  do k = 0,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j = 1,n_sh
      p_layer_centres_sh(j,k)     = p_layer_centres(shi(j),k)
      p_layer_boundaries_sh(j,k)  = p_layer_boundaries(shi(j),k)
      exner_layer_centres_sh(j,k) = exner_layer_centres(shi(j),k)
      exner_layer_boundaries_sh(j,k)                                           &
         = exner_layer_boundaries(shi(j),k)
      r_theta_sh(j,k)             = r_theta(shi(j),k)
    end do
  end do

  do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j = 1,n_sh
      u_sh(j,k)           = u(shi(j),k)
      v_sh(j,k)           = v(shi(j),k)
      w_sh(j,k)           = w(shi(j),k)
      th_sh(j,k)          = th(shi(j),k)
      exner_rho_levels_sh(j,k)  = exner_rho_levels(shi(j),k)
      z_theta_sh(j,k)     = z_theta(shi(j),k)
      z_rho_sh(j,k)       = z_rho(shi(j),k)
      rho_sh(j,k)         = rho(shi(j),k)
      rho_theta_sh(j,k)   = rho_theta(shi(j),k)
      r2rho_sh(j,k)       = r2rho(shi(j),k)
      r2rho_th_sh(j,k)    = r2rho_th(shi(j),k)
      rho_dry_theta_temp(j,k) = rho_dry_theta(shi(j),k)
      rho_dry_temp(j,k)       = rho_dry(shi(j),k)
      r_rho_sh(j,k)       = r_rho(shi(j),k)
      dr_across_rh_sh(j,k)       = dr_across_rh(shi(j),k)
      dr_across_th_sh(j,k)       = dr_across_th(shi(j),k)
    end do
  end do
  if (l_progn_tnuc) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j = 1, n_sh
      tnuc_nlcl_sh(j) = tnuc_nlcl(shi(j))
    end do
  end if
  if (l_conv_prog_precip) then  ! Initialise precip based prognostic
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1, n_sh
        conv_prog_precip_sh(j,k) = conv_prog_precip(shi(j),k)
      end do
    end do
  end if

  if (l_conv_prog_flx) then  ! Initialise mass flux based prognostic
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1, n_sh
        conv_prog_flx_sh(j,k) = conv_prog_flx(shi(j),k)
      end do
    end do
  end if

  !
  ! moisture  - input depends on scheme
  !
  if (iconv_shallow == 1) then     ! G-R scheme
    !
    ! G-R requires input of specific humidity
    !
    if (l_mr_physics) then

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_sh
          denom = 1.0 / ( 1.0 + mt(shi(j),k) )
          q_sh(j,k)   = q(shi(j),k)   * denom
          qse_sh(j,k) = qse(shi(j),k) ! copy specific value
          qcl_sh(j,k) = qcl(shi(j),k) * denom
          qcf_sh(j,k) = qcf(shi(j),k) * denom
        end do
      end do

    else       ! Input is specific humidity therefore no problems

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_sh
          q_sh(j,k)   = q(shi(j),k)
          qse_sh(j,k) = qse(shi(j),k) ! copy specific value
          qcl_sh(j,k) = qcl(shi(j),k)
          qcf_sh(j,k) = qcf(shi(j),k)
        end do
      end do

    end if     ! Test on l_mr_physics

  end if        ! test on scheme

  ! Water tracers - allocate compressed sh arrays and convert between
  ! q and mixing ratio if required

  call wtrac_alloc_conv_e(n_sh, nlev, n_wtrac, wtrac_e)

  if (l_wtrac_conv) then
    call wtrac_gather_conv(ncells,seg_size,n_sh,nlev,n_wtrac, iconv_shallow,    &
                           shi, mt, q_wtrac, qcl_wtrac, qcf_wtrac, wtrac_e)
  end if

  !
  ! In/out variables
  !

  do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j = 1,n_sh
      cf_liquid_sh(j,k)   = cf_liquid(shi(j),k)
      cf_frozen_sh(j,k)   = cf_frozen(shi(j),k)
    end do
  end do

  if (l_tracer) then
    do ktra = 1,ntra
      do k = 1,trlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_sh
          tracer_sh(j,k,ktra)  = tracer(shi(j),k,ktra)
        end do
      end do
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_sh
          dtrabydt_sh(j,k,ktra)  = 0.0
        end do
      end do
    end do
  end if

  allocate(dt_dd_temp(n_sh,nlev))
  allocate(dq_dd_temp(n_sh,nlev))
  allocate(area_ud_temp(n_sh,nlev))


  !-----------------------------------------------------------------------
  ! 2.2 Call shallow convection code
  !-----------------------------------------------------------------------

  if (iconv_shallow == 1) then !  G-R mass flux type scheme

    call shallow_conv_6a(                                                      &
                            !in
       nbl,nlev,ntra,n_wtrac,n_cca_lev,n_sh,trlev,                             &
       bland_sh,                                                               &
       delthvu_sh,                                                             &
       exner_rho_levels_sh,                                                    &
       exner_layer_centres_sh,                                                 &
       exner_layer_boundaries_sh,                                              &
       l_q_interact,l_tracer,ntml_sh,ntpar_sh,                                 &
       pstar_sh,p_layer_centres_sh,                                            &
       p_layer_boundaries_sh,                                                  &
       z_theta_sh, z_rho_sh,                                                   &
       r_theta_sh, r_rho_sh,                                                   &
       rho_theta_sh,                                                           &
       rho_dry_theta_temp, rho_dry_temp,                                       &
       r2rho_th_sh, r2rho_sh,                                                  &
       dr_across_th_sh, dr_across_rh_sh,                                       &
       conv_prog_precip_sh,                                                    &
       q_sh,th_sh,timestep,                                                    &
       u_sh,v_sh,uw0_sh,vw0_sh,wstar_sh,wthvs_sh,                              &
       entrain_coef_sh,delta_smag_sh,                                          &
       zlcl_uv_sh,tnuc_nlcl_sh,ztop_uv_sh,freeze_lev_sh,                       &
       recip_pstar_sh,qse_sh,                                                  &
       l_scm_convss_dg,                                                        &
       ccp_strength_sh,                                                        &
                            !INOUT
       cf_frozen_sh,cf_liquid_sh,                                              &
       qcf_sh,qcl_sh,tracer_sh,wtrac_e,                                        &
       w2p_sh,conv_prog_flx_sh,scm_convss_dg_c(1:n_sh),                        &
                            !out
       cape_out_sh,cclwp_sh,ccw_sh,cca_sh,                                     &
       dbcfbydt_sh,dcffbydt_sh,dcflbydt_sh,                                    &
       dqbydt_sh,dqcfbydt_sh,dqclbydt_sh,                                      &
       dthbydt_sh,dubydt_sh,dvbydt_sh,dtrabydt_sh,                             &
       detrain_up_sh,detrain_dwn_sh,entrain_up_sh,                             &
       entrain_dwn_sh,                                                         &
       iccb_sh,icct_sh,                                                        &
       lcca_sh,lcbase_sh,lctop_sh,                                             &
       rain_sh,snow_sh, rain_3d_sh, snow_3d_sh,                                &
       up_flux_sh,up_flux_half_sh,                                             &
       dwn_flux_sh,uw_shall_sh,vw_shall_sh,                                    &
       tcw_sh,cca_2d_sh,kterm_sh,ind_shall_sh,                                 &
       dt_dd_temp,dq_dd_temp,area_ud_temp )

  end if        ! test on shallow convection scheme

  !-----------------------------------------------------------------------
  ! 2.3 Write data from shallow convection points to full arrays
  !-----------------------------------------------------------------------

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i = 1,n_sh
    cape_out(shi(i))           = cape_out_sh(i)
    cclwp(shi(i))              = cclwp_sh(i)
    iccb(shi(i))               = iccb_sh(i)
    icct(shi(i))               = icct_sh(i)
    lcca(shi(i))               = lcca_sh(i)
    lcbase(shi(i))             = lcbase_sh(i)
    lctop(shi(i))              = lctop_sh(i)

    rain(shi(i))               = rain_sh(i)
    snow(shi(i))               = snow_sh(i)
    precip_shall(shi(i))       = rain_sh(i) + snow_sh(i)
    wstar_dn(shi(i))           = wstar_sh(i)  ! wstar_dn
    wstar_up(shi(i))           = wstar_up_sh(i)
    mb1(shi(i))                = mb1_sh(i)
    mb2(shi(i))                = mb2_sh(i)
    ind_shall(shi(i))          = ind_shall_sh(i)
    kterm_shall(shi(i))        = kterm_sh(i)
    cca_2d(shi(i)) = cca_2d_sh (i)
  end do

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        rain_wtrac(shi(i),i_wt) = wtrac_e(i_wt)%rain(i)
        snow_wtrac(shi(i),i_wt) = wtrac_e(i_wt)%snow(i)
      end do
    end do
  end if

  if (l_mom) then
    do k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        dubydt(shi(i),k)         = dubydt_sh(i,k)
        dvbydt(shi(i),k)         = dvbydt_sh(i,k)
      end do
    end do
  end if

  if (flg_w_eqn) then
    do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_sh
        w2p(shi(i),k) = w2p_sh(i,k)
      end do
    end do
  end if

  if (l_conv_prog_flx) then
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1, n_sh
        conv_prog_flx(shi(i),k) = conv_prog_flx_sh(i,k)
      end do
    end do
  end if

  if (flg_area_ud) then
    do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_sh
        area_ud(shi(i),k) = area_ud_temp(i,k)
      end do
    end do
  end if

  do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_sh
      dthbydt(shi(i),k)        = dthbydt_sh(i,k)
      dcflbydt(shi(i),k)       = dcflbydt_sh(i,k)
      dcffbydt(shi(i),k)       = dcffbydt_sh(i,k)
      dbcfbydt(shi(i),k)       = dbcfbydt_sh(i,k)
      ccw(shi(i),k)            = ccw_sh(i,k)
      dt_dd(shi(i),k)          = dt_dd_temp(i,k)
    end do
  end do

  do k = 1,n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_sh
      cca(shi(i),k) = cca_sh(i,k)
    end do
  end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_sh
    iccb0  (shi(i)) = iccb  (shi(i))
    icct0  (shi(i)) = icct  (shi(i))
    cclwp0 (shi(i)) = cclwp (shi(i))
    lcbase0(shi(i)) = lcbase(shi(i))
  end do

  do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1, n_sh
      ccw0(shi(i),k) = ccw_sh(i,k) * ccw_sh_knob
    end do
  end do

  do k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1, n_sh
      cca0(shi(i),k) = cca_sh(i,k) * cca_sh_knob
    end do
  end do

  if (l_cca_sh_prog) then
    do k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_sh
        cca0_sh(shi(i),k) = cca_sh(i,k)
      end do
    end do
  end if


  ! Set past convective depth BUT only if dump holds space for history
  ! prognostics

  if (l_conv_hist) then
    if (iconv_shallow == 1) then

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        cld_depth = z_rho_sh(i,kterm_sh(i)+1) - z_rho_sh(i,ntml_sh(i)+1)

        ! In some cases shallow fails and kterm remains set to zero so cld_depth
        ! will be negative.
        if ( cld_depth < 0.0 ) then
          cld_depth = 0.0
        end if
        past_conv_ht(shi(i)) = cld_depth

      end do
    else      ! shallow turbulence scheme so cloud ntpar
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        cld_depth = z_rho_sh(i,ntpar_sh(i)+1) - z_rho_sh(i,ntml_sh(i)+1)
        past_conv_ht(shi(i)) = cld_depth
      end do
    end if
  end if

    !
    ! moisture  - output depends on scheme
    !
  if (iconv_shallow == 1) then     ! G-R scheme
    !
    ! G-R requires input of specific humidity
    !
    if (l_mr_physics) then ! requires conversion

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          denom = 1.0 + mt(shi(i),k)
          dqbydt(shi(i),k)    = denom * dqbydt_sh(i,k)
          dqclbydt(shi(i),k)  = denom * dqclbydt_sh(i,k)
          dqcfbydt(shi(i),k)  = denom * dqcfbydt_sh(i,k)
          dq_dd(shi(i),k)     = denom * dq_dd_temp(i,k)
        end do
      end do

    else   ! output is specific humidity therefore no problems

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          dqbydt(shi(i),k)   = dqbydt_sh(i,k)
          dqclbydt(shi(i),k) = dqclbydt_sh(i,k)
          dqcfbydt(shi(i),k) = dqcfbydt_sh(i,k)
          dq_dd(shi(i),k)    = dq_dd_temp(i,k)
        end do
      end do

    end if     ! Test on l_mr_physics

  else
    !
    ! turbulence scheme outputs mixing ratios
    !
    if (l_mr_physics) then

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          dqbydt(shi(i),k)   = dqbydt_sh(i,k)
          dqclbydt(shi(i),k) = dqclbydt_sh(i,k)
          dqcfbydt(shi(i),k) = dqcfbydt_sh(i,k)
          dq_dd(shi(i),k)    = dq_dd_temp(i,k)
        end do
      end do

    else    ! Requires conversion of m increment to q

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          denom = 1.0 / ( 1.0 + mt(shi(i),k) )
          dqbydt(shi(i),k)   = denom * dqbydt_sh(i,k)
          dqclbydt(shi(i),k) = denom * dqclbydt_sh(i,k)
          dqcfbydt(shi(i),k) = denom * dqcfbydt_sh(i,k)
          dq_dd(shi(i),k)    = denom * dq_dd_temp(i,k)
        end do
      end do

    end if       ! Test on l_mr_physics

  end if        ! test on scheme

  ! Convert water tracers between specific humidity and mixing ratios if
  ! required and decompress arrays
  if (l_wtrac_conv) then

    call wtrac_scatter_conv(ncells,seg_size,n_sh,nlev,n_wtrac,                  &
                            iconv_shallow, shi, mt, wtrac_e,                   &
                            dqbydt_wtrac, dqclbydt_wtrac, dqcfbydt_wtrac)
  end if



  deallocate(area_ud_temp)
  deallocate(dq_dd_temp)
  deallocate(dt_dd_temp)
  deallocate(rho_dry_theta_temp)
  deallocate(rho_dry_temp)
  !
  ! diagnostics output
  !
  if (flg_up_flx) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        up_flux(shi(i),k)        = up_flux_sh(i,k)
      end do
    end do
  end if
  if (flg_up_flx_half .and. iconv_shallow == 1) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        up_flux_half(shi(i),k)        = up_flux_half_sh(i,k)
      end do
    end do
  end if
  ! Downdraught mass flux not available from shallow convection
  ! (it is calculated; why it isn't scattered back here is unknown!)
  !  if (flg_dwn_flx) then
  !    do k = 1,nlev
  !!DIR$ IVDEP
  !!DIR$ VECTOR ALWAYS
  !      do i = 1,n_sh
  !        dwn_flux(shi(i),k)       = dwn_flux_sh(i,k)
  !      end do
  !    end do
  !  end if
  if (flg_uw_shall) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        uw_shall(shi(i),k)       = uw_shall_sh(i,k)
      end do
    end do
  end if
  if (flg_vw_shall) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        vw_shall(shi(i),k)       = vw_shall_sh(i,k)
      end do
    end do
  end if


  if (flg_mf_shall) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        mf_shall(shi(i),k)        = up_flux_sh(i,k)
      end do
    end do
  end if
  if (iconv_shallow == 1) then
    if (flg_entr_up) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          entrain_up(shi(i),k)     = entrain_up_sh(i,k)
        end do
      end do
    end if
    if (flg_detr_up) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          detrain_up(shi(i),k)     = detrain_up_sh(i,k)
        end do
      end do
    end if
    if (flg_entr_dwn) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          entrain_dwn(shi(i),k)    = entrain_dwn_sh(i,k)
        end do
      end do
    end if
    if (flg_detr_dwn) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          detrain_dwn(shi(i),k)    = detrain_dwn_sh(i,k)
        end do
      end do
    end if
  end if  ! iconv_shallow equal to 1
  if (flg_dt_shall) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        dt_shall(shi(i),k)        = dthbydt_sh(i,k)
      end do
    end do
  end if
  if (flg_dq_shall) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_sh
        dq_shall(shi(i),k)        = dqbydt_sh(i,k)
      end do
    end do
  end if
  if (l_mom) then
    if (flg_du_shall) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          du_shall(shi(i),k)        = dubydt_sh(i,k)
        end do
      end do
    end if
    if (flg_dv_shall) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          dv_shall(shi(i),k)        = dvbydt_sh(i,k)
        end do
      end do
    end if
  end if

  if (l_tracer) then
    do ktra = 1,ntra
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_sh
          dtrabydt(shi(i),k,ktra)  = dtrabydt_sh(i,k,ktra)
        end do
      end do
    end do
  end if

  ! Merge sh 3d rain & snow profiles

  do k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1,n_sh
      rain_3d(shi(i),k) = rain_3d_sh(i,k)
      snow_3d(shi(i),k) = snow_3d_sh(i,k)
    end do
  end do

  ! Deallocate water tracer shallow convection arrays
  call wtrac_dealloc_conv_e(n_wtrac, wtrac_e)

end if ! n_sh > 0


!-----------------------------------------------------------------------
! 3.0 CONGESTUS Convection   - call depending on switches
!-----------------------------------------------------------------------

if (iconv_congestus ==  1) then

  if (n_cg  >   0) then   ! there must be congestus points

    allocate(rho_dry_temp(n_cg,nlev))
    allocate(rho_dry_theta_temp(n_cg,nlev))

    !-----------------------------------------------------------------------
    ! 3.1 Compress input variable arrays for congestus convection scheme to
    !     length n_cg (congestus points only)
    !-----------------------------------------------------------------------

    j = 0
    do i = 1,seg_size
      if (cumulus_bl(i) .and. l_congestus(i)) then
        j                        = j+1
        cgi(j)                   = i
      end if
    end do
    !
    ! In only variables
    !
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do j = 1,n_cg
      bland_cg(j)           = bland(cgi(j))
      delthvu_cg(j)         = delthvu_bl(cgi(j))
      ntml_cg(j)            = ntml(cgi(j))
      ntpar_cg(j)           = ntpar(cgi(j))
      pstar_cg(j)           = pstar(cgi(j))
      recip_pstar_cg(j)     = recip_pstar(cgi(j))
      uw0_cg(j)             = uw0(cgi(j))
      vw0_cg(j)             = vw0(cgi(j))
      wstar_cg(j)           = wstar(cgi(j))
      wthvs_cg(j)           = wthvs(cgi(j))
      zlcl_uv_cg(j)         = zlcl_uv(cgi(j))
      ztop_uv_cg(j)         = ztop_uv(cgi(j))
      freeze_lev_cg(j)      = freeze_lev(cgi(j))
      entrain_coef_cg(j)    = entrain_coef(cgi(j))
      kterm_cg(j) = 0     ! initialise array Will only get set
                          ! if call deep scheme
      ccp_strength_cg(j)    = ccp_strength(cgi(j))
    end do

    do k = 0,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1,n_cg
        p_layer_centres_cg(j,k)    = p_layer_centres(cgi(j),k)
        p_layer_boundaries_cg(j,k) = p_layer_boundaries(cgi(j),k)
        exner_layer_centres_cg(j,k)= exner_layer_centres(cgi(j),k)
        exner_layer_boundaries_cg(j,k)                                         &
           = exner_layer_boundaries(cgi(j),k)
        r_theta_cg(j,k)            = r_theta(cgi(j),k)
      end do
    end do

    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1,n_cg
        u_cg(j,k)           = u(cgi(j),k)
        v_cg(j,k)           = v(cgi(j),k)
        th_cg(j,k)          = th(cgi(j),k)
        exner_rho_levels_cg(j,k)  = exner_rho_levels(cgi(j),k)
        z_theta_cg(j,k)     = z_theta(cgi(j),k)
        z_rho_cg(j,k)       = z_rho(cgi(j),k)
        rho_theta_cg(j,k)   = rho_theta(cgi(j),k)
        rho_dry_theta_temp(j,k) = rho_dry_theta(cgi(j),k)
        rho_dry_temp(j,k)       = rho_dry(cgi(j),k)
        r2rho_cg(j,k)       = r2rho(cgi(j),k)
        r2rho_th_cg(j,k)    = r2rho_th(cgi(j),k)
        r_rho_cg(j,k)       = r_rho(cgi(j),k)
        dr_across_rh_cg(j,k)       = dr_across_rh(cgi(j),k)
        dr_across_th_cg(j,k)       = dr_across_th(cgi(j),k)
      end do
    end do

    if (l_progn_tnuc) then
      do j = 1, n_cg
        tnuc_nlcl_cg(j) = -10.0
      end do
    end if

    if (l_conv_prog_precip) then  ! Initialise precip based prognostic
      do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1, n_cg
          conv_prog_precip_cg(j,k) = conv_prog_precip(cgi(j),k)
        end do
      end do
    end if

    if (l_conv_prog_flx) then  ! Initialise mass flux based prognostic
      do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1, n_cg
          conv_prog_flx_cg(j,k) = conv_prog_flx(cgi(j),k)
        end do
      end do
    end if

    !
    ! moisture  - input depends on scheme
    !
    if (iconv_congestus == 1) then     ! G-R scheme
      !
      ! G-R requires input of specific humidity
      !
      if (l_mr_physics) then

        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do j = 1,n_cg
            denom = 1.0 / ( 1.0 + mt(cgi(j),k) )
            q_cg(j,k)   = q(cgi(j),k)   * denom
            qse_cg(j,k) = qse(cgi(j),k) ! copy specific value
            qcl_cg(j,k) = qcl(cgi(j),k) * denom
            qcf_cg(j,k) = qcf(cgi(j),k) * denom
          end do
        end do

      else   ! input is specific humidity therefore no problems

        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do j = 1,n_cg
            q_cg(j,k)   = q(cgi(j),k)
            qse_cg(j,k) = qse(cgi(j),k) ! copy specific value
            qcl_cg(j,k) = qcl(cgi(j),k)
            qcf_cg(j,k) = qcf(cgi(j),k)
          end do
        end do

      end if     ! l_mr_physics

    else         ! turbulence scheme requires mixing ratio
      !
      ! turbulence scheme requires mixing ratios as input
      !
      if (l_mr_physics) then

        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do j = 1,n_cg
            q_cg(j,k)   = q(cgi(j),k)
            qse_cg(j,k) = qse_mix(cgi(j),k) ! copy mixing ratio value
            qcl_cg(j,k) = qcl(cgi(j),k)
            qcf_cg(j,k) = qcf(cgi(j),k)
          end do
        end do

      else       ! Require conversion from q to m

        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do j = 1,n_cg

            q_cg(j,k) = q(cgi(j),k)*(1.0+mt(cgi(j),k))

            ! check for negative values
            if (q_cg(j,k) <  0.0) then
              if (printstatus >= prstatus_normal) then
                write(umMessage,'(a17,g26.18,a5,2i6)') 'problem q_mix -ve ',   &
                      q_cg(j,k),' j,k ',j,k
                call umPrint(umMessage,src='glue_conv-6a')
              end if
              q_cg(j,k) = 0.0
            end if
            qse_cg(j,k) = qse_mix(cgi(j),k) ! copy mixing ratio value
            ! No checks on qcl and qcf as not expected to be used or altered

            qcl_cg(j,k) = qcl(cgi(j),k)*(1.0+mt(cgi(j),k))
            qcf_cg(j,k) = qcf(cgi(j),k)*(1.0+mt(cgi(j),k))
          end do
        end do

      end if      ! Test on l_mr_physics

    end if        ! Test on scheme (iconv_congestus)

    !
    ! In/out variables
    !
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1,n_cg
        cf_liquid_cg(j,k)   = cf_liquid(cgi(j),k)
        cf_frozen_cg(j,k)   = cf_frozen(cgi(j),k)
        bulk_cf_cg(j,k)     = bulk_cf(cgi(j),k)
      end do
    end do

    if (l_tracer) then
      do ktra = 1,ntra
        do k = 1,trlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do j = 1,n_cg
            tracer_cg(j,k,ktra)  = tracer(cgi(j),k,ktra)
          end do
        end do
        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do j = 1,n_cg
            dtrabydt_cg(j,k,ktra)  = 0.0
          end do
        end do
      end do
    end if

    allocate(dt_dd_temp(n_cg,nlev))
    allocate(dq_dd_temp(n_cg,nlev))
    allocate(area_ud_temp(n_cg,nlev))

    !-----------------------------------------------------------------------
    ! 3.2 Call congestus convection code
    !-----------------------------------------------------------------------

    if (iconv_congestus == 1) then   ! G-R type scheme

      ! new mass flux routine (at present a copy of shallow scheme)

      call congest_conv_6a(nbl,nlev,ntra,n_cca_lev,n_cg,trlev                  &
         ,              bland_cg                                               &
         ,              delthvu_cg                                             &
         ,              exner_rho_levels_cg                                    &
         ,              exner_layer_centres_cg                                 &
         ,              exner_layer_boundaries_cg                              &
         ,              L_q_interact                                           &
         ,              L_tracer, ntml_cg, ntpar_cg                            &
         ,              pstar_cg,p_layer_centres_cg                            &
         ,              p_layer_boundaries_cg                                  &
         ,              z_theta_cg, z_rho_cg                                   &
         ,              r_theta_cg, r_rho_cg                                   &
         ,              rho_theta_cg                                           &
         ,              rho_dry_theta_temp, rho_dry_temp                       &
         ,              r2rho_th_cg, r2rho_cg                                  &
         ,              dr_across_th_cg, dr_across_rh_cg                       &
         ,              conv_prog_precip_cg                                    &
         ,              q_cg,th_cg,timestep                                    &
         ,              u_cg,v_cg,uw0_cg,vw0_cg,wstar_cg,wthvs_cg              &
         ,              entrain_coef_cg                                        &
         ,              zlcl_uv_cg,tnuc_nlcl_cg,ztop_uv_cg,freeze_lev_cg       &
         ,              recip_pstar_cg,qse_cg                                  &
         ,              ccp_strength_cg                                        &

                        ! InOut
         ,              bulk_cf_cg,cf_frozen_cg,cf_liquid_cg,qcf_cg            &
         ,              qcl_cg,tracer_cg,w2p_cg,conv_prog_flx_cg               &

                        ! Out
         ,              cape_out_cg,cclwp_cg,ccw_cg,cca_cg                     &
         ,              dbcfbydt_cg,dcffbydt_cg,dcflbydt_cg                    &
         ,              dqbydt_cg,dqcfbydt_cg,dqclbydt_cg,dthbydt_cg           &
         ,              dubydt_cg,dvbydt_cg,dtrabydt_cg                        &
         ,              detrain_up_cg,detrain_dwn_cg                           &
         ,              entrain_up_cg,entrain_dwn_cg                           &
         ,              iccb_cg,icct_cg,lcca_cg,lcbase_cg,lctop_cg             &
         ,              rain_cg,snow_cg,rain_3d_cg,snow_3d_cg                  &
         ,              up_flux_cg, up_flux_half_cg                            &
         ,              dwn_flux_cg,uw_shall_cg,vw_shall_cg,kterm_cg           &
         ,              tcw_cg,cca_2d_cg, dt_dd_temp,dq_dd_temp                &
         ,              area_ud_temp )


      !         Else If (iconv_congestus == 2) then  ! turbulence scheme
      !
      !            Call TURB_CONGESTUS_CONV(
      !
      !        Not available yet - still to be written
      !
    end if   ! test of type of congestus scheme

    !-----------------------------------------------------------------------
    ! 3.3 Write data from congestus convection points to full arrays
    !-----------------------------------------------------------------------
    !
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_cg
      cape_out(cgi(i))           = cape_out_cg(i)
      cclwp(cgi(i))              = cclwp_cg(i)
      iccb(cgi(i))               = iccb_cg(i)
      icct(cgi(i))               = icct_cg(i)
      lcca(cgi(i))               = lcca_cg(i)
      lcbase(cgi(i))             = lcbase_cg(i)
      lctop(cgi(i))              = lctop_cg(i)

      rain(cgi(i))               = rain_cg(i)
      snow(cgi(i))               = snow_cg(i)
      precip_cong(cgi(i))        = rain_cg(i) + snow_cg(i)
      wstar_dn(cgi(i))           = wstar_cg(i)  ! wstar_dn
      ! may be required if congestus not forced to stop at ntpar
      kterm_congest(cgi(i))      = kterm_cg(i)
      ! required at present to replicate old code
      cca_2d(cgi(i))             = cca_2d_cg(i)
    end do

    if (l_mom) then
      do k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_cg
          dubydt(cgi(i),k)         = dubydt_cg(i,k)
          dvbydt(cgi(i),k)         = dvbydt_cg(i,k)
        end do
      end do
    end if

    if (flg_w_eqn) then
      do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i=1, n_cg
          w2p(cgi(i),k) = w2p_cg(i,k)
        end do
      end do
    end if

    if (l_conv_prog_flx) then
      do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1, n_cg
          conv_prog_flx(cgi(i),k) = conv_prog_flx_cg(i,k)
        end do
      end do
    end if

    if (flg_area_ud) then
      do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i=1, n_cg
          area_ud(cgi(i),k) = area_ud_temp(i,k)
        end do
      end do
    end if

    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_cg
        dthbydt(cgi(i),k)        = dthbydt_cg(i,k)
        dcflbydt(cgi(i),k)       = dcflbydt_cg(i,k)
        dcffbydt(cgi(i),k)       = dcffbydt_cg(i,k)
        dbcfbydt(cgi(i),k)       = dbcfbydt_cg(i,k)
        ccw(cgi(i),k)            = ccw_cg(i,k)
        dt_dd(cgi(i),k)          = dt_dd_temp(i,k)
      end do
    end do

    do k = 1,n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_cg
        cca(cgi(i),k) = cca_cg(i,k)
      end do
    end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1, n_cg
      iccb0  (cgi(i)) = iccb  (cgi(i))
      icct0  (cgi(i)) = icct  (cgi(i))
      cclwp0 (cgi(i)) = cclwp (cgi(i))
      lcbase0(cgi(i)) = lcbase(cgi(i))
    end do

    do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_cg
        ccw0(cgi(i),k) = ccw_cg(i,k) * ccw_sh_knob
      end do
    end do

    do k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_cg
        cca0(cgi(i),k) = cca_cg(i,k) * cca_sh_knob
      end do
    end do

    !
    ! moisture  - output depends on scheme
    !
    if (iconv_congestus == 1) then     ! G-R scheme
      !
      ! Outputs specific humidity increments
      !
      if (l_mr_physics) then ! requires conversion

        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do i = 1,n_cg
            denom = 1.0 + mt(cgi(i),k)
            dqbydt(cgi(i),k)   = denom * dqbydt_cg(i,k)
            dqclbydt(cgi(i),k) = denom * dqclbydt_cg(i,k)
            dqcfbydt(cgi(i),k) = denom * dqcfbydt_cg(i,k)
            dq_dd(cgi(i),k)    = denom * dq_dd_temp(i,k)
          end do
        end do

      else   ! output is specific humidity therefore no problems

        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do i = 1,n_cg
            dqbydt(cgi(i),k)   = dqbydt_cg(i,k)
            dqclbydt(cgi(i),k) = dqclbydt_cg(i,k)
            dqcfbydt(cgi(i),k) = dqcfbydt_cg(i,k)
            dq_dd(cgi(i),k)    = dq_dd_temp(i,k)
          end do
        end do

      end if     ! l_mr_physics

    else         ! turbulence scheme
      !
      ! turbulence scheme outputs mixing ratios
      !
      if (l_mr_physics) then

        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do i = 1,n_cg
            dqbydt(cgi(i),k)   = dqbydt_cg(i,k)
            dqclbydt(cgi(i),k) = dqclbydt_cg(i,k)
            dqcfbydt(cgi(i),k) = dqcfbydt_cg(i,k)
            dq_dd(cgi(i),k)    = dq_dd_temp(i,k)
          end do
        end do

      else    ! Requires conversion of increment from m to q

        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do i = 1,n_cg
            denom = 1.0 / ( 1.0 + mt(cgi(i),k) )
            dqbydt(cgi(i),k)   = denom * dqbydt_cg(i,k)
            dqclbydt(cgi(i),k) = denom * dqclbydt_cg(i,k)
            dqcfbydt(cgi(i),k) = denom * dqcfbydt_cg(i,k)
            dq_dd(cgi(i),k)    = denom * dq_dd_temp(i,k)
          end do
        end do

      end if     ! Test on l_mr_physics

    end if        ! Test on scheme (iconv_congestus)

    deallocate(area_ud_temp)
    deallocate(dq_dd_temp)
    deallocate(dt_dd_temp)
    deallocate(rho_dry_theta_temp)
    deallocate(rho_dry_temp)

    ! Merge congestus 3d rain & snow profiles

    do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_cg
        rain_3d(cgi(i),k) = rain_3d_cg(i,k)
        snow_3d(cgi(i),k) = snow_3d_cg(i,k)
      end do
    end do

    if (flg_up_flx) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_cg
          up_flux(cgi(i),k)        = up_flux_cg(i,k)
        end do
      end do
    end if
    if (flg_up_flx_half) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_cg
          up_flux_half(cgi(i),k)        = up_flux_half_cg(i,k)
        end do
      end do
    end if

    if (flg_mf_congest) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_cg
          mf_congest(cgi(i),k)        = up_flux_cg(i,k)
        end do
      end do
    end if
    if (flg_dt_congest) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_cg
          dt_congest(cgi(i),k)        = dthbydt_cg(i,k)
        end do
      end do
    end if
    if (flg_dq_congest) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_cg
          dq_congest(cgi(i),k)        = dqbydt_cg(i,k)
        end do
      end do
    end if
    if (l_mom) then
      if (flg_du_congest) then
        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do i = 1,n_cg
            du_congest(cgi(i),k)        = dubydt_cg(i,k)
          end do
        end do
      end if
      if (flg_dv_congest) then
        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do i = 1,n_cg
            dv_congest(cgi(i),k)        = dvbydt_cg(i,k)
          end do
        end do
      end if
    end if
    if (l_tracer) then
      do ktra = 1,ntra
        do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do i = 1,n_cg
            dtrabydt(cgi(i),k,ktra)  = dtrabydt_cg(i,k,ktra)
          end do
        end do
      end do
    end if

  end if    ! number of congestus points > 0


  ! set past convective depth

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i = 1,n_cg
    cld_depth = z_rho_cg(i,kterm_cg(i)+1) - z_rho_cg(i,ntml_cg(i)+1)
    past_conv_ht(cgi(i)) = cld_depth
  end do

end if      ! iconv_congestus > 0

!-----------------------------------------------------------------------
! 4.0 MID-LEVEL Convection
! 4.1 Set lowest level that mid level convection can trigger
!-----------------------------------------------------------------------

if (iconv_mid >  0 .and. n_md > 0) then

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1,n_md
    cca_2d_md(i) = 0.0
  end do

  allocate(rho_dry_temp(n_md,nlev))
  allocate(rho_dry_theta_temp(n_md,nlev))
  !
  ! Create index of points where mid-level is possible
  !
  j = 0
  do i = 1,seg_size
    if (l_mid(i)) then
      j      = j+1
      mdi(j) = i
    end if
  end do

  select case (midtrig_opt)
  case (mtrig_ntmlplus2)
    do i = 1,n_md
      ! If there is no succcessful shallow, congestus or deep convection
      ! then start mid-level from 2 layers above BL top (ntml)
      ! If there are cumulus points then start mid-level from 2 layers
      ! above where shallow, congestus or deep have terminated.
      midtrig(i) = max( ntml(mdi(i))+2, kterm_shall(mdi(i))+2,                 &
                        kterm_congest(mdi(i))+2, kterm_deep(mdi(i))+2 )
    end do
  case (mtrig_ntmlplus1)
    do i = 1,n_md
      ! If there is no succcessful shallow, congestus or deep convection
      ! then start mid-level from layer above BL top (ntml)
      ! If there are cumulus points then start mid-level from layer
      ! above where shallow, congestus or deep has terminated.
      midtrig(i) = max( ntml(mdi(i))+1, kterm_shall(mdi(i))+1,                 &
                        kterm_congest(mdi(i))+1, kterm_deep(mdi(i))+1 )
    end do
  case (mtrig_ntml)
    do i = 1,n_md
      ! If there is no succcessful shallow, congestus or deep convection
      ! then start mid-level from layer BL top (ntml)
      ! If there are cumulus points then start mid-level from layer
      ! above where shallow, congestus or deep has terminated.
      midtrig(i) = max( ntml(mdi(i)), kterm_shall(mdi(i))+1,                   &
                        kterm_congest(mdi(i))+1, kterm_deep(mdi(i))+1 )
    end do
  case (mtrig_surface)
    do i = 1,n_md
      ! If there is no succcessful shallow, congestus or deep convection
      ! then start mid-level layer one.
      ! If there are cumulus points then start mid-level from layer
      ! above where shallow, congestus or deep has terminated.
      midtrig(i) = max( 1, kterm_shall(mdi(i))+1,                              &
                        kterm_congest(mdi(i))+1, kterm_deep(mdi(i))+1 )
    end do
  end select  !midtrig_opt

  !-----------------------------------------------------------------------
  ! 4.2 Copy all input arrays to arrays ending in _md for passing to
  !     mid-level scheme
  !-----------------------------------------------------------------------

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do j=1,n_md
    bland_md(j)           = bland(mdi(j))
    ntml_md(j)            = ntml(mdi(j))
    ntpar_md(j)           = ntpar(mdi(j))
    pstar_md(j)           = pstar(mdi(j))
    recip_pstar_md(j)     = recip_pstar(mdi(j))
    W_MAX_md(j)           = w_max(mdi(j))
    freeze_lev_md(j)      = freeze_lev(mdi(j))
    g_ccp_md(j)           = g_ccp(mdi(j))
    h_ccp_md(j)           = h_ccp(mdi(j))
    ccp_strength_md(j)    = ccp_strength(mdi(j))
    ftl_md(j)             = ftl(mdi(j))
    fqt_md(j)             = fqt(mdi(j))
  end do

  do k = 0,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_md
      p_layer_centres_md(i,k)        = p_layer_centres(mdi(i),k)
      p_layer_boundaries_md(i,k)     = p_layer_boundaries(mdi(i),k)
      exner_layer_centres_md(i,k)    = exner_layer_centres(mdi(i),k)
      exner_layer_boundaries_md(i,k) = exner_layer_boundaries(mdi(i),k)
      r_theta_md(i,k)                = r_theta(mdi(i),k)
    end do
  end do
  do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_md
      u_md(i,k)           = u(mdi(i),k)
      v_md(i,k)           = v(mdi(i),k)
      w_md(i,k)           = w(mdi(i),k)
      th_md(i,k)          = th(mdi(i),k)
      cf_liquid_md(i,k)   = cf_liquid(mdi(i),k)
      cf_frozen_md(i,k)   = cf_frozen(mdi(i),k)
      exner_rho_levels_md(i,k)  = exner_rho_levels(mdi(i),k)
      z_theta_md(i,k)     = z_theta(mdi(i),k)
      z_rho_md(i,k)       = z_rho(mdi(i),k)
      r_rho_md(i,k)       = r_rho(mdi(i),k)
      r2rho_th_md(i,k)     = r2rho_th(mdi(i),k)
      r2rho_md(i,k)        = r2rho(mdi(i),k)
      rho_theta_md(i,k)    = rho_theta(mdi(i),k)
      rho_dry_theta_temp(i,k) = rho_dry_theta(mdi(i),k)
      rho_dry_temp(i,k)       = rho_dry(mdi(i),k)
      rho_md(i,k)          = rho(mdi(i),k)
      dr_across_th_md(i,k) = dr_across_th(mdi(i),k)
      dr_across_rh_md(i,k) = dr_across_rh(mdi(i),k)
    end do
  end do

  if (l_progn_tnuc) then
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1, n_md
        tnuc_new_md(j,k) = tnuc_new(mdi(j),k)
      end do
    end do
  end if


  if (l_conv_prog_precip) then  ! Initialise precip based prognostic
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1, n_md
        conv_prog_precip_md(j,k) = conv_prog_precip(mdi(j),k)
      end do
    end do
  end if

  if (l_conv_prog_flx) then  ! Initialise mass flux based prognostic
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j = 1, n_md
        conv_prog_flx_md(j,k) = conv_prog_flx(mdi(j),k)
      end do
    end do
  end if


  do i = 1,n_md
    iccb_md(i)   = iccb(mdi(i))
    icct_md(i)   = icct(mdi(i))
  end do

  !
  ! moisture  - input depends on scheme
  !

  if (iconv_mid == 1) then     ! G-R scheme
    !
    ! G-R requires input of specific humidity
    !
    if (l_mr_physics) then

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_md
          denom = 1.0 / ( 1.0 + mt(mdi(j),k) )
          q_md(j,k)   = q(mdi(j),k)   * denom
          qse_md(j,k) = qse(mdi(j),k) ! copy specific value
          qcl_md(j,k) = qcl(mdi(j),k) * denom
          qcf_md(j,k) = qcf(mdi(j),k) * denom
        end do
      end do

    else   ! input is specific humidity therefore no problems

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_md
          q_md(j,k)   = q(mdi(j),k)
          qse_md(j,k) = qse(mdi(j),k)  ! copy specific value
          qcl_md(j,k) = qcl(mdi(j),k)
          qcf_md(j,k) = qcf(mdi(j),k)
        end do
      end do

    end if      ! Test on l_mr_physics

  else         ! future turbulence scheme requires mixing ratio
    !
    ! turbulence scheme requires mixing ratios as input
    !
    if (l_mr_physics) then

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_md
          q_md(j,k)   = q(mdi(j),k)
          qse_md(j,k) = qse_mix(mdi(j),k)  ! copy mixing ratio value
          qcl_md(j,k) = qcl(mdi(j),k)
          qcf_md(j,k) = qcf(mdi(j),k)
        end do
      end do

    else       ! Require conversion

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_md
          q_md(j,k)     = q(mdi(j),k)*(1.0+mt(mdi(j),k))
        end do
      end do

      if (printstatus >= prstatus_normal) then
        do k = 1,nlev
          do j = 1,n_md
            ! check for negative values
            if (q_md(j,k) <  0.0) then
              write(umMessage,'(a18,g26.18,a5,2i6)') 'problem q_mix -ve ',     &
                 q_md(j,k),' j,k ',j,k
              call umPrint(umMessage,src='glue_conv-6a')
            end if
          end do
        end do
      end if

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do j = 1,n_md
          if (q_md(j,k) <  0.0) then
            q_md(j,k) = 0.0
          end if
          qse_md(j,k) = qse_mix(mdi(j),k)  ! copy mixing ratio value

          ! Not sure whether qcl and qcf will be used or altered and don't
          ! know whether additional tests will be required if converting.

          qcl_md(j,k) = qcl(mdi(j),k)*(1.0+mt(mdi(j),k))
          qcf_md(j,k) = qcf(mdi(j),k)*(1.0+mt(mdi(j),k))
        end do
      end do

    end if     ! Test on l_mr_physics

  end if        ! test on scheme (iconv_mid)

  ! Water tracers - set up compressed md arrays and convert between
  ! q and mixing ratio if required

  call wtrac_alloc_conv_e(n_md, nlev, n_wtrac, wtrac_e)

  if (l_wtrac_conv) then
    call wtrac_gather_conv(ncells,seg_size,n_md,nlev,n_wtrac, iconv_mid,        &
                           mdi, mt, q_wtrac, qcl_wtrac, qcf_wtrac, wtrac_e)
  end if

  if (l_tracer) then
    do ktra = 1,ntra
      do k = 1,trlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          tracer_md(i,k,ktra)  = tracer(mdi(i),k,ktra)
        end do
      end do
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          dtrabydt_md(i,k,ktra)  = 0.0
        end do
      end do
    end do
  end if

  allocate(dt_dd_temp(n_md,nlev))
  allocate(dq_dd_temp(n_md,nlev))
  allocate(area_ud_temp(n_md,nlev))


  !-----------------------------------------------------------------------
  ! 4.3 Call mid-level convection code
  !-----------------------------------------------------------------------

  if (iconv_mid == 1) then ! Gregory-Rowntree mid level convection

    call mid_conv_6a(                                                          &
                  !in
                  nbl,nlev,ntra,n_wtrac,n_cca_lev,n_md,trlev,                  &
                  bland_md,w_max_md,                                           &
                  exner_rho_levels_md,                                         &
                  exner_layer_centres_md,                                      &
                  exner_layer_boundaries_md,                                   &
                  l_q_interact, l_tracer, midtrig, ntml_md,                    &
                  ntpar_md, freeze_lev_md,                                     &
                  pstar_md,p_layer_centres_md,                                 &
                  p_layer_boundaries_md,                                       &
                  z_theta_md, z_rho_md,                                        &
                  r_theta_md, r_rho_md,                                        &
                  rho_theta_md, rho_md,                                        &
                  rho_dry_theta_temp, rho_dry_temp,                            &
                  r2rho_th_md,  r2rho_md,                                      &
                  dr_across_th_md, dr_across_rh_md,                            &
                  conv_prog_precip_md,                                         &
                  q_md, th_md, timestep, ftl_md, fqt_md,                       &
                  u_md,v_md,w_md,recip_pstar_md,qse_md,                        &
                  l_scm_convss_dg,                                             &
                  g_ccp_md,h_ccp_md,ccp_strength_md,                           &
                  !INOUT
                  cf_frozen_md,cf_liquid_md,                                   &
                  qcf_md,qcl_md,tracer_md,wtrac_e,                             &
                  w2p_md,conv_prog_flx_md,scm_convss_dg_c(1:n_md),             &
                  !out
                  cape_out_md,cclwp_md,ccw_md,cca_md,                          &
                  dbcfbydt_md,dcffbydt_md,dcflbydt_md,                         &
                  dqbydt_md,dqcfbydt_md,dqclbydt_md,                           &
                  dthbydt_md,dubydt_md,dvbydt_md,dtrabydt_md,                  &
                  detrain_up_md,detrain_dwn_md,entrain_up_md,                  &
                  entrain_dwn_md,                                              &
                  iccb_md,icct_md,lcca_md,lcbase_md,lctop_md,                  &
                  rain_md,snow_md,rain_3d_md,snow_3d_md,                       &
                  up_flux_md,up_flux_half_md,                                  &
                  dwn_flux_md,l_mid_md,cca_2d_md,                              &
                  uw_mid_md,vw_mid_md, cfl_limited_md,                         &
                  dt_dd_temp, dq_dd_temp, area_ud_temp,                        &
                  error_point,tnuc_new_md )

    if (error_point /= 0) then
      errorstatus = 3   ! will cause model to fail
      write(cmessage,'(a47,i12,a16,i2)')                                       &
      'Mid conv went to the top of the model at point ',                       &
      mdi(error_point),' in seg on call ',call_number

      call ereport(routinename, errorstatus, cmessage )

    end if

  else        ! turbulence based alternative
    cmessage=' New mid-level scheme not available yet '
    ierror = 1      ! fatal error will cause model to stop
    call ereport( routinename, ierror, cmessage)
    !          call mid_turb_conv( ) ?
  end if

  !-----------------------------------------------------------------------
  ! 4.4 Write data from mid-level convection to full arrays
  !-----------------------------------------------------------------------

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i = 1,n_md
    !
    ! Cloud variables - only overwrite deep or shallow values if
    ! iccb_md  and icct_md > 0
    !
    if (iccb_md(i)  >   0 .and. icct_md(i)  >   0) then
      iccb(mdi(i))                    = iccb_md(i)
      icct(mdi(i))                    = icct_md(i)
    end if
    !
    ! Overwrite lowest cloud values only if cumulus = .F.
    !
    if (.not. cumulus_bl(mdi(i))) then
      lcca(mdi(i))                    = lcca_md(i)
      lcbase(mdi(i))                  = lcbase_md(i)
      lctop(mdi(i))                   = lctop_md(i)
    end if
  end do
  !
  ! Write remaining data to full arrays
  ! Split from loop with if tests to encourage optimisation
  !
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i = 1,n_md
    cape_out(mdi(i))                = cape_out(mdi(i)) + cape_out_md(i)
    cclwp(mdi(i))                   = cclwp(mdi(i))    + cclwp_md(i)
    rain(mdi(i))                    = rain(mdi(i))     + rain_md(i)
    snow(mdi(i))                    = snow(mdi(i))     + snow_md(i)
    precip_mid(mdi(i))              = rain_md(i)  + snow_md(i)
    mid_cfl_limited(mdi(i))         = cfl_limited_md(i)
    l_mid_all(mdi(i))       = l_mid_md(i)
    !===================================================================
    ! NOTE: At this point cca_2d_md is only equal
    !       to that from mid-level cloud on a given grid point.
    !===================================================================

  end do

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        rain_wtrac(mdi(i),i_wt) =                                              &
                     rain_wtrac(mdi(i),i_wt) + wtrac_e(i_wt)%rain(i)
        snow_wtrac(mdi(i),i_wt) =                                              &
                     snow_wtrac(mdi(i),i_wt) + wtrac_e(i_wt)%snow(i)
      end do
    end do
  end if

  ! Merge md 3d rain & snow profiles

  do k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1,n_md
      rain_3d(mdi(i),k) = rain_3d(mdi(i),k) + rain_3d_md(i,k)
      snow_3d(mdi(i),k) = snow_3d(mdi(i),k) + snow_3d_md(i,k)
    end do
  end do

  if (l_mom) then
    if (flg_uw_mid) then
      do k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i=1,n_md
          uw_mid(mdi(i),k)     = uw_mid(mdi(i),k) + uw_mid_md(i,k)
        end do
      end do
    end if
    if (flg_vw_mid) then
      do k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i=1,n_md
          vw_mid(mdi(i),k)     = vw_mid(mdi(i),k) + vw_mid_md(i,k)
        end do
      end do
    end if
    do k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        dubydt(mdi(i),k)            = dubydt(mdi(i),k) + dubydt_md(i,k)
        dvbydt(mdi(i),k)            = dvbydt(mdi(i),k) + dvbydt_md(i,k)
      end do
    end do
  end if

  do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_md
      dthbydt(mdi(i),k)     = dthbydt(mdi(i),k)  + dthbydt_md(i,k)
      dcflbydt(mdi(i),k)    = dcflbydt(mdi(i),k) + dcflbydt_md(i,k)
      dcffbydt(mdi(i),k)    = dcffbydt(mdi(i),k) + dcffbydt_md(i,k)
      dbcfbydt(mdi(i),k)    = dbcfbydt(mdi(i),k) + dbcfbydt_md(i,k)
      ccw(mdi(i),k)         = ccw(mdi(i),k)      + ccw_md(i,k)
      dt_dd(mdi(i),k)       = dt_dd(mdi(i),k)    + dt_dd_temp(i,k)
    end do
  end do

  do k = 1,n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i = 1,n_md
      cca(mdi(i),k) = cca(mdi(i),k) + cca_md(i,k)
    end do
  end do

  if (flg_w_eqn) then
    do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_md
        w2p(mdi(i),k) = w2p(mdi(i),k) + w2p_md(i,k)
      end do
    end do
  end if

  if (l_conv_prog_flx) then
    do k = 1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1, n_md
        conv_prog_flx(mdi(i),k) = conv_prog_flx_md(i,k)
      end do
    end do
  end if

  if (flg_area_ud) then
    do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_md
        area_ud(mdi(i),k) = area_ud(mdi(i),k) + area_ud_temp(i,k)
      end do
    end do
  end if

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_md
    cclwp0 (mdi(i)) = cclwp (mdi(i))
    iccb0  (mdi(i)) = iccb  (mdi(i))
    icct0  (mdi(i)) = icct  (mdi(i))
    lcbase0(mdi(i)) = lcbase(mdi(i))
  end do

  ! Then scaling is done by CCRad, even if PC2 = .true.
  ! To zero CCA (section 0) for PC2 use the CCRAD knobs
  ! cca_knobs should be set to 0.0 in gui/namelist for PC2

  ! Anvils if selected will affect both section 0/5 diags
  ! if ccrads are not set to 0.0, this may required further
  ! attention in future
  do k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1, n_md
      ccw0(mdi(i),k) = ccw0(mdi(i),k) + ccw_md(i,k) * ccw_md_knob
    end do
  end do

  do k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do i=1, n_md
      cca0(mdi(i),k) = cca0(mdi(i),k) + cca_md(i,k) * cca_md_knob
    end do
  end do

  ! Add on cca_2d_sh/cca_2d_cg/cca_2d_dp to cca_2d.
  ! Contributions fron cca_2d_sh and cca_2d_dp have already been
  ! applied to cca_sh/cca_cg/cca_dp in the _conv routines

  ! However, contributions from cca_2d_sh/cca_2d_cg/cca_2d_dp
  ! need to be included in the cca_2d diagnostic
  ! (which at this point only holds cca_2d_md)
  ! so as not to upset any Downstream products that use it

  ! cca_2d_md entered mid_conv as an empty array,
  ! i.e. no contribution from shallow/congestus/deep.
  ! Add any contribution from mid level to preserve
  ! CCA_2d diagnostic.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do i=1, n_md
    cca_2d(mdi(i)) = cca_2d(mdi(i)) + cca_2d_md(i)
  end do


  if (l_cca_md_prog) then
    do k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i=1, n_md
        cca0_md(mdi(i),k) = cca_md(i,k)
      end do
    end do
  end if

  !
  ! moisture  - output depends on scheme
  !
  if (iconv_mid == 1) then     ! G-R scheme
    !
    ! Output in specific humidity
    !
    if (l_mr_physics) then ! requires conversion

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          denom = 1.0 + mt(mdi(i),k)
          dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + denom * dqbydt_md(i,k)
          dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + denom * dqclbydt_md(i,k)
          dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + denom * dqcfbydt_md(i,k)
          dq_dd(mdi(i),k)    = dq_dd(mdi(i),k)    + denom * dq_dd_temp(i,k)
        end do
      end do

    else   ! output is specific humidity therefore no problems

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + dqbydt_md(i,k)
          dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + dqclbydt_md(i,k)
          dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + dqcfbydt_md(i,k)
          dq_dd(mdi(i),k)    = dq_dd(mdi(i),k)    + dq_dd_temp(i,k)
        end do
      end do

    end if     ! Test on scheme (iconv_mid)

  else         ! turbulence scheme requires mixing ratio
    !
    ! turbulence scheme outputs mixing ratios
    !
    if (l_mr_physics) then

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + dqbydt_md(i,k)
          dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + dqclbydt_md(i,k)
          dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + dqcfbydt_md(i,k)
          dq_dd(mdi(i),k)    = dq_dd(mdi(i),k)    + dq_dd_temp(i,k)
        end do
      end do

    else    ! Requires conversion of increments from mixing ratio to q

      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          denom = 1.0/(1.0+mt(mdi(i),k))
          dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + denom * dqbydt_md(i,k)
          dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + denom * dqclbydt_md(i,k)
          dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + denom * dqcfbydt_md(i,k)
          dq_dd(mdi(i),k)    = dq_dd(mdi(i),k)    + denom * dq_dd_temp(i,k)
        end do
      end do

    end if  ! Test on l_mr_physics

  end if        ! test on scheme (iconv_mid)

  ! Convert water tracers between specific humidity and mixing ratios if
  ! required and decompress arrays
  if (l_wtrac_conv) then

    call wtrac_scatter_conv_mid(ncells,seg_size,n_md,nlev,n_wtrac,              &
                                iconv_mid, mdi, mt, wtrac_e,                   &
                                dqbydt_wtrac, dqclbydt_wtrac, dqcfbydt_wtrac)
  end if


  deallocate(area_ud_temp)
  deallocate(dq_dd_temp)
  deallocate(dt_dd_temp)
  deallocate(rho_dry_theta_temp)
  deallocate(rho_dry_temp)

  if (flg_up_flx) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        up_flux(mdi(i),k)    = up_flux(mdi(i),k) + up_flux_md(i,k)
      end do
    end do
  end if
  if (flg_up_flx_half) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        up_flux_half(mdi(i),k)= up_flux_half(mdi(i),k)+up_flux_half_md(i,k)
      end do
    end do
  end if
  if (flg_dwn_flx) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        dwn_flux(mdi(i),k)    = dwn_flux(mdi(i),k) + dwn_flux_md(i,k)
      end do
    end do
  end if
  if (flg_entr_up) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        entrain_up(mdi(i),k)  = entrain_up(mdi(i),k)+ entrain_up_md(i,k)
      end do
    end do
  end if
  if (flg_detr_up) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        detrain_up(mdi(i),k)  = detrain_up(mdi(i),k) + detrain_up_md(i,k)
      end do
    end do
  end if
  if (flg_entr_dwn) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        entrain_dwn(mdi(i),k) = entrain_dwn(mdi(i),k) + entrain_dwn_md(i,k)
      end do
    end do
  end if
  if (flg_detr_dwn) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        detrain_dwn(mdi(i),k) = detrain_dwn(mdi(i),k) + detrain_dwn_md(i,k)
      end do
    end do
  end if
  if (flg_mf_midlev) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        mf_midlev(mdi(i),k)        = up_flux_md(i,k)
      end do
    end do
  end if
  if (flg_dt_midlev) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        dt_midlev(mdi(i),k)        = dthbydt_md(i,k)
      end do
    end do
  end if
  if (flg_dq_midlev) then
    do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do i = 1,n_md
        dq_midlev(mdi(i),k)        = dqbydt_md(i,k)
      end do
    end do
  end if
  if (l_mom) then
    if (flg_du_midlev) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          du_midlev(mdi(i),k)        = dubydt_md(i,k)
        end do
      end do
    end if
    if (flg_dv_midlev) then
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          dv_midlev(mdi(i),k)        = dvbydt_md(i,k)
        end do
      end do
    end if
  end if
  if (l_tracer) then
    do ktra = 1,ntra
      do k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        do i = 1,n_md
          dtrabydt(mdi(i),k,ktra) = dtrabydt(mdi(i),k,ktra)                    &
             + dtrabydt_md(i,k,ktra)
        end do
      end do
    end do
  end if

  ! Deallocate water tracer mid convection arrays
  call wtrac_dealloc_conv_e(n_wtrac, wtrac_e)

end if     ! test on iconv_mid

!-----------------------------------------------------------------------
! 5.0 Work out which points convection has occurred at.
!-----------------------------------------------------------------------

nconv_all=0
do i = 1,seg_size
  if (cumulus_bl(i) .or. l_mid_all(i)) then
    nconv_all = nconv_all + 1
    index1(nconv_all) = i
  end if
end do

!-----------------------------------------------------------------------
! 6.0  Update tracer field
!      More efficient to do here rather than in subroutines.
!      This seems to be an expensive bit of code on NEC.
!      Changed to operate on just convective points.
!-----------------------------------------------------------------------

if (l_tracer .and. (nconv_all >  0)) then

  !
  ! Adjust timestep to prevent any negative values invading the tracer
  ! fields (adjusted timestep is a function of geographical  location and
  ! tracer type.
  !
  do ktra=1,ntra
    ! initialise step to timestep
    do i = 1,nconv_all
      limited_step(i) =timestep
    end do

    do k=1,nlev
      do j = 1,nconv_all
        i=index1(j)
        ! negative increments  may have a problem
        if (dtrabydt(i,k,ktra)  <   0.0 ) then

          step_test(j) = (0.9999 * abs(tracer(i,k,ktra))) /                    &
             (abs(dtrabydt(i,k,ktra)) + safety_margin)

          if (step_test(j)   <   limited_step(j) ) then
            ! then increment is bigger than tracer and timestep
            ! needs to be reduced
            limited_step (j) = step_test(j)
          end if

        end if
      end do
    end do
    !
    ! Update tracer field using limited_step.
    !
    do k = 1,nlev
      ! NEC compiler directive
      !CDIR NODEP
      do j = 1,nconv_all
        i=index1(j)
        tracer(i,k,ktra) = tracer(i,k,ktra) +                                  &
           dtrabydt(i,k,ktra) * limited_step(j)
      end do
    end do

  end do  ! ktra loop
end if   ! L_tracer

! Added as something wrong with top level CCA
!  zeroing cca
do i=1,seg_size
  cca(i,n_cca_lev) = 0.0
  cca0(i,n_cca_lev) = 0.0
end do


! Deallocate the array of structures
if (max(n_dp, n_sh, n_md) > 0) then
  deallocate( scm_convss_dg_c )
end if


    !-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine glue_conv_6a

end module glue_conv_6a_mod
