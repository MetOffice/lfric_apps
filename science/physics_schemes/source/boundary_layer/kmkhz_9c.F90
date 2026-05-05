! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To calculate the non-local turbulent mixing
!           coefficients KM and KH

!  Programming standard: UMDP3

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module kmkhz_9c_mod

use um_types, only: rbl_eps, r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'KMKHZ_9C_MOD'
contains

subroutine kmkhz_9c (                                                          &
! in levels/switches
 bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                                       &
! in fields
 p,rho_wet_tq,rho_mix,rho_mix_tq,t,q,qcl,qcf,cf,qw,tl,dzl,rdz,z_tq,z_uv,       &
 rad_hr,micro_tends,                                                           &
 bt,bq,btm,bqm,dqsdt,btm_cld,bqm_cld,a_qs,a_qsm,a_dqsdtm,                      &
 v_s,fb_surf,rhostar_gb,ntpar,zh_prev,                                         &
 zhpar,z_lcl,zmaxb_for_dsc,zmaxt_for_dsc,l_shallow,wtrac_as,                   &
! INOUT fields
 ftl,fqw,zh,dzh,cumulus,ntml,w,etadot,t1_sd,q1_sd,wtrac_bl,                    &
! out fields
 rhokm,rhokh,rhokm_top,rhokh_top,zhsc,                                         &
 unstable,dsc,coupled,sml_disc_inv,dsc_disc_inv,                               &
 ntdsc,nbdsc,f_ngstress,tke_nl,grad_t_adj, grad_q_adj,                         &
 rhof2, rhofsc, ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb,                          &
 tothf_zh, tothf_zhsc, totqf_zh, totqf_zhsc,                                   &
 kent, we_lim, t_frac_tr, zrzi_tr,                                             &
 kent_dsc, we_lim_dsc, t_frac_dsc_tr, zrzi_dsc_tr,                             &
 k_plume                                                                       &
 )

use atm_fields_bounds_mod, only: pdims, tdims, tdims_l, pdims_s


use bl_diags_mod, only: strnewbldiag
use tuning_segments_mod, only:  bl_segment_size
use nlsizes_namelist_mod, only: model_levels
use bl_option_mod, only:                                                       &
    sc_cftol, entr_enhance_by_cu, Buoyrev_feedback, on, off, l_new_kcloudtop,  &
    sc_diag_opt, sc_diag_orig, sc_diag_cu_relax,                               &
    sc_diag_cu_rh_max, sc_diag_all_rh_max,                                     &
    a_grad_adj, max_t_grad, flux_grad, Locketal2000,                           &
    HoltBov1993, LockWhelan2006, entr_smooth_dec, entr_taper_zh,               &
    kprof_cu, l_use_sml_dsc_fixes, l_converge_ga,                              &
    bl_res_inv, cosine_inv_flux, target_inv_profile, pr_max, l_noice_in_turb,  &
    one_third, two_thirds, zero, one, one_half
use conversions_mod, only: pi => pi_bl
use cv_run_mod, only: l_param_conv
use gen_phys_inputs_mod, only: l_mr_physics
use level_heights_mod, only: eta_theta_levels
use model_domain_mod, only: model_type, mt_single_column
use missing_data_mod, only: rmdi
use planet_constants_mod, only: vkman => vkman_bl, r => rd_bl,                 &
     c_virtual => c_virtual_bl, g => g_bl, etar => etar_bl, grcp => grcp_bl,   &
     lcrcp => lcrcp_bl, lsrcp => lsrcp_bl
use s_scmop_mod,  only: default_streams,                                       &
     t_avg, d_bl, d_sl, scmdiag_bl
use timestep_mod, only: timestep
use water_constants_mod, only: tm => tm_bl

use qsat_mod, only: qsat, qsat_mix, qsat_wat, qsat_wat_mix

use excf_nl_9c_mod, only: excf_nl_9c

use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use wtrac_bl_mod,            only: bl_wtrac_type
use wtrac_atm_step_mod,      only: atm_step_wtrac_type
use kmkhz_9c_wtrac_mod,      only: calc_dqw_inv_wtrac, calc_fqw_inv_wtrac

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! in arguments
integer, intent(in) ::                                                         &
 bl_levels
                            ! in No. of atmospheric levels for
                            !    which boundary layer fluxes are
                            !    calculated.

!     Declaration of new BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

! Additional variables for SCM diagnostics which are dummy in full UM
integer, intent(in) ::                                                         &
  nSCMDpkgs             ! No of SCM diagnostics packages

logical, intent(in) ::                                                         &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

real(kind=r_bl), intent(in) ::                                                 &
 zmaxb_for_dsc,                                                                &
 zmaxt_for_dsc
                            ! in Max heights to look for DSC cloud
                            !    base and top

integer, intent(in) ::                                                         &
 ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! in Top level of parcel ascent.
                            !    Used in convection scheme.
                            !    NOTE: CAN BE > BL_LEVELS-1

logical, intent(in) ::                                                         &
 l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! in Flag to indicate shallow
                            !    convection

! Water tracer structure containing 'micro_tends' fields
type(atm_step_wtrac_type), intent(in) :: wtrac_as(n_wtrac)

real(kind=r_bl), intent(in) ::                                                 &
 bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                            ! in A buoyancy parameter for clear air
                            !    on p,T,q-levels (full levels).
 bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                            ! in A buoyancy parameter for clear air
                            !    on p,T,q-levels (full levels).
 bqm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                            ! in A buoyancy parameter for clear air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 btm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                            ! in A buoyancy parameter for clear air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 bqm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                            ! in A buoyancy parameter for cloudy air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 btm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                            ! in A buoyancy parameter for cloudy air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                            ! in Saturated lapse rate factor
                            !    on p,T,q-levels (full levels).
 a_qsm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                            ! in Saturated lapse rate factor
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 a_dqsdtm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          bl_levels),                                                          &
                            ! in Saturated lapse rate factor
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                            ! in Partial derivative of QSAT w.r.t.
                            !    temperature.
 p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,0:bl_levels),           &
                            ! in P(*,K) is pressure at full level k.
 qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                            ! in Total water content (kg per kg air).
 tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                            ! in Liquid/frozen water temperature (K).
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),             &
                            ! in Temperature (K).
 qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                            ! in Cloud ice (kg per kg air)
 qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                            ! in Cloud liquid water
 q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,                &
   tdims_l%k_start:bl_levels),                                                 &
                            ! in specific humidity
 cf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                            ! in Cloud fractions for boundary levs.
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                            ! in Z_tq(*,K) is the height of the
                            !    k-th full level above the surface.
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1),        &
                            ! in Z_uv(*,K) is the height of level
                            !       k-1/2 above the surface (m).
 dzl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                            ! in Layer depths (m).  DZL(,K) is the
                            !    distance from layer boundary K-1/2
                            !    to layer boundary K+1/2.  For K=1
                            !    the lower boundary is the surface.
 rdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                            ! in Reciprocal of distance between
                            !    full levels (m-1).  1/RDZ(,K) is
                            !    the vertical distance from level
                            !    K-1 to level K, except that for
                            !    K=1 it is the height of the
                            !    lowest atmospheric full level.
 rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels),                                                           &
                            ! in density on UV (ie. rho) levels,
                            !    used in RHOKH so dry density if
                            !    L_mr_physics is true
 rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                            ! in density on TQ (ie. theta) levels,
                            !    used in RHOKM so wet density
 rho_mix_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels)
                            ! in density on TQ (ie. theta) levels,
                            !    used in non-turb flux integration
                            !    so dry density if L_mr_physics is true

real(kind=r_bl), intent(in) ::                                                 &
 v_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                     &
                            ! in Surface friction velocity (m/s)
 fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                            ! in Surface buoyancy flux over density
                            !       (m^2/s^3).
 rhostar_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            ! in Surface air density in kg/m3
 z_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                            ! in Height of lifting condensation
                            !    level.
 zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                            ! in Height of top of NTPAR
                            !    NOTE: CAN BE ABOVE BL_LEVELS-1
 zh_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! in boundary layer height (m) from
                            !    previous timestep

real(kind=r_bl), intent(in) ::                                                 &
 rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2,bl_levels),                                                          &
                            ! in (LW,SW) radiative heating rates (K/s)
 micro_tends(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             2, bl_levels)
                            ! in Tendencies from microphysics
                            !    (TL, K/s; QW, kg/kg/s)

! INOUT arrays
integer, intent(in out) ::                                                     &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
    ! INOUT Number of model levels in the turbulently mixed layer.

logical, intent(in out) ::                                                     &
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! INOUT Flag for Cu in the bl

real(kind=r_bl), intent(in out) ::                                             &
 zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                      &
                            ! INOUT Boundary layer height (m).
 dzh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                     &
                            ! INOUT inversion thickness (m)
 ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
 fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                            ! INOUT "Explicit" fluxes of TL and QW
                            !       (rho*Km/s, rho*m/s)
                            !       in:  level 1 (surface flux)
                            !       out: entrainment-level flux
 t1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                            ! INOUT Standard Deviation of level 1
                            !    temperature (K).
 q1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                            ! INOUT Standard Deviation of level 1
                            !    specific humidity (kg/kg).
 w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,0:bl_levels),           &
                            ! INOUT Vertical velocity (m/s)
 etadot(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        0:bl_levels)
                            ! INOUT d(ETA)/dt

! Water tracer structure containing boundary layer fields
type(bl_wtrac_type), intent(in out) :: wtrac_bl(n_wtrac)

! out arrays
integer, intent(out) ::                                                        &
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                            ! out Top level for turb mixing in
                            !       cloud layer
 nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                            ! out Bottom level of any decoupled
                            !       turbulently mixed Sc layer
 sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                            ! out Flags for whether discontinuous
 dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                            ! out   inversions are diagnosed
 kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                            ! out Grid-levels of SML and DSC
 kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                            ! out   inversions (for tracer mixing)
 k_plume(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! out   Start grid-level for surface-driven plume

logical, intent(out) ::                                                        &
 unstable(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                            ! out Flag to indicate an unstable
                            !     surface layer.
 dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                     &
                            ! out Flag set if decoupled
                            !     stratocumulus layer found
 coupled(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! out Flag to indicate Sc layer weakly
                            !     decoupled (implies mixing at SML
                            !     top is through K profiles rather
                            !     than entrainment parametrization)

real(kind=r_bl), intent(out) ::                                                &
 rhokm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       2:bl_levels),                                                           &
                            ! out Non-local turbulent mixing
                            !     coefficient for momentum.
 rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
       2:bl_levels),                                                           &
                            ! out Non-local turbulent mixing
                            !     coefficient for scalars.
 rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
       2:bl_levels),                                                           &
                            ! out Top-down turbulent mixing
                            !     coefficient for momentum.
 rhokh_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
       2:bl_levels),                                                           &
                            ! out Top-down turbulent mixing
                            !     coefficient for scalars.
 tke_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
       2:bl_levels),                                                           &
                            ! out Non-local TKE diag (currently times rho)
 f_ngstress(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,       &
            2:bl_levels),                                                      &
                            ! out dimensionless function for
                            !     non-gradient stresses
 rhof2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
       2:bl_levels),                                                           &
 rhofsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
         2:bl_levels)
                            ! out f2 and fsc term shape profiles
                            !       multiplied by rho

real(kind=r_bl), intent(out) ::                                                &
  ft_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels+1),      &
                            ! out Non-turbulent heat and moisture
  fq_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels+1)
                            !       fluxes (rho*Km/s, rho*m/s)

real(kind=r_bl), intent(out) ::                                                &
 tothf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                            ! out Total heat fluxes at inversions
 tothf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            !     (rho*Km/s)
 totqf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                            ! out Total moisture fluxes at
 totqf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            !      inversions (rho*m/s)
 ft_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            ! out Non-turbulent heat and moisture
 fq_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            !      flux at the base of the DSC layer
 grad_t_adj(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            ! out Temperature gradient adjustment
                            !     for non-local mixing in unstable
                            !     turbulent boundary layer.
 grad_q_adj(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            ! out Humidity gradient adjustment
                            !     for non-local mixing in unstable
                            !     turbulent boundary layer.
 zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! out Cloud layer height (m).
    ! The following are used in tracer mixing.
    ! At 9B 3 elements were used - here just (i,j,2) is used.
real(kind=r_bl), intent(out) ::                                                &
 we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),                &
                            ! out rho*entrainment rate implied by
                            !     placing of subsidence
 zrzi_tr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),               &
                            ! out (z-z_base)/(z_i-z_base)
 t_frac_tr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),             &
                            ! out a fraction of the timestep
 we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),            &
                            ! out rho*entrainment rate implied by
                            !     placing of subsidence
 zrzi_dsc_tr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),           &
                            ! out (z-z_base)/(z_i-z_base)
 t_frac_dsc_tr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)
                            ! out a fraction of the timestep

!----------------------------------------------------------------------
!    Local and other symbolic constants :-

character(len=*), parameter ::  RoutineName = 'KMKHZ_9C'

real(kind=r_bl) ::                                                             &
        a_ga_hb93,a_ga_lw06,max_svl_grad,dfsw_frac,                            &
        ct_resid,svl_coup,svl_coup_max,fgf
parameter (                                                                    &
 a_ga_hb93=7.2_r_bl,                                                           &
 a_ga_lw06=10.0_r_bl,                                                          &
 max_svl_grad=1.0e-3_r_bl,                                                     &
                          ! maximum SVL gradient in a mixed layer
 ct_resid=200.0_r_bl,                                                          &
                          ! Parcel cloud-top residence time (in s)
 svl_coup_max=one,                                                             &
 svl_coup=one_half,                                                            &
                          ! Parameters controlling positioning of
                          ! surface-driven entrainment
 dfsw_frac = 0.35_r_bl,                                                        &
                          ! Fraction of SW flux difference to
                          ! contribute to net flux difference across
                          ! cloud top
 fgf=zero)                 ! Adiabatic gradient factor for ice

!  Define local storage.

!  (a) Workspace.

logical :: cloud_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                        ! Flag set when cloud base is reached.
           dsc_save(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                        ! Copy of DSC needed to indicate
                        !   decoupling diagnosed in EXCF_NL

real(kind=r_bl) ::                                                             &
        qs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),     &
                        ! Saturated sp humidity at pressure
                        !   and temperature of sucessive levels.
        cfl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                        ! Liquid cloud fraction.
        cff(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                        ! Frozen cloud fraction.
        dqcldz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels), &
                        ! Vertical gradient of in-cloud liquid cloud
                        !   water in a well-mixed layer.
        dqcfdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels), &
                        ! Vertical gradient of in-cloud frozen cloud
                        !   water in a well-mixed layer.
        sls_inc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                        ! SL and QW increments due to large-scale
        qls_inc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                        !    vertical advection (K s^-1, s^-1)
        df_over_cp(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end,bl_levels),                       &
                        ! Radiative flux change over layer / c_P
        dflw_over_cp(pdims%i_start:pdims%i_end,                                &
                     pdims%j_start:pdims%j_end,bl_levels),                     &
                        ! LW radiative flux change over layer / c_P
        dfsw_over_cp(pdims%i_start:pdims%i_end,                                &
                     pdims%j_start:pdims%j_end,bl_levels),                     &
                        ! SW radiative flux change over layer / c_P
        svl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                        ! Liquid/frozen water virtual temperature / CP
        sl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),     &
                        ! TL + G*Z/CP (K)
        z_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                        ! Z_TOP(*,K) is the height of
                        !   level k+1/2 above the surface.
        w_grad(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels)
                        ! Gradient of w

real(kind=r_bl) ::                                                             &
        db_ga_dry(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end,2:bl_levels),                      &
                        ! Out-of-cloud (DRY) and in-cloud buoyancy
        db_noga_dry(pdims%i_start:pdims%i_end,                                 &
                    pdims%j_start:pdims%j_end,2:bl_levels),                    &
                        !   jumps used in flux integral calculation.
        db_ga_cld(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end,2:bl_levels),                      &
                        !   GA terms include gradient adjustment
        db_noga_cld(pdims%i_start:pdims%i_end,                                 &
                    pdims%j_start:pdims%j_end,2:bl_levels)
                        !   arising from non-gradient fluxes. (m/s2)

!-----------------------------------------------------------------------
! The following fluxes, flux changes are in units of rho*Km/s for heat
! and rho*m/s for humidity
!-----------------------------------------------------------------------
real(kind=r_bl) ::                                                             &
 dfmic (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels,2),                                                          &
                                         ! Flux changes from microphys
 dfsubs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels,2),                                                          &
                                         !   and subsidence
 frad (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
        bl_levels+1),                                                          &
                                         ! Fluxes from net radiation,
 frad_lw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels+1),                                                         &
                                         !   LW,
 frad_sw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels+1),                                                         &
                                         !   SW,
 fmic (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
        bl_levels+1,2),                                                        &
                                         !   microphys and subsidence;
 fsubs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
        bl_levels+1,2),                                                        &
                                         !   for T, Q separately.
 svl_flux(bl_levels+1)                   ! svl flux across resolved inversion

real(kind=r_bl) ::                                                             &
        bflux_surf(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! Buoyancy flux at the surface.
        bflux_surf_sat(pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end),                             &
                        ! Saturated-air surface buoyancy flux
        db_top(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! Buoyancy jump at the top of the BL
        db_dsct(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
                        ! Buoyancy jump at the DSC layer top
        df_top_over_cp(pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end),                             &
                        ! Radiative flux change at cloud top / c_P
        df_dsct_over_cp(pdims%i_start:pdims%i_end,                             &
                        pdims%j_start:pdims%j_end),                            &
                        ! Radiative flux change at DSC top / CP
        env_svl_km1(pdims%i_start:pdims%i_end,                                 &
                    pdims%j_start:pdims%j_end),                                &
                        ! Density potential temperature layer K-1
        qcl_ic_top(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
        qcf_ic_top(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! In-cloud liquid and frozen water contents
                        !   at the top of the model layer
        bt_top(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! Buoyancy parameter at the top of the SML
        bt_dsct(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
                        !   and DSC
        btt_top(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
                        ! In-cloud buoyancy param at the top of the BL
        btt_dsct(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end),                                   &
                        !   and DSC
        btc_top(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
                        ! Cloud fraction weighted buoyancy parameter
        db_top_cld(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! In-cloud buoyancy jump at the top of the BL
        db_dsct_cld(pdims%i_start:pdims%i_end,                                 &
                    pdims%j_start:pdims%j_end),                                &
                        !   and DSC
        cld_factor(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! Fraction of grid box potentially giving
        cld_factor_dsc(pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end),                             &
                        !   evaporative entrainment, for SML and DSC
        chi_s_top(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end),                                  &
                        ! Mixing fraction of just saturated mixture
        chi_s_dsct(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        !   at top of the SML and DSC layer
        zeta_s(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! Non-cloudy fraction of mixing layer for
                        !   surface forced entrainment term.
        zeta_r(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! Non-cloudy fraction of mixing layer for
                        !   cloud-top radiative cooling entrainment
        zeta_r_dsc(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        !   term in SML and DSC layers
        zc(pdims%i_start:pdims%i_end,                                          &
           pdims%j_start:pdims%j_end),                                         &
                        ! Cloud depth (not cloud fraction weighted).
        zc_dsc(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        !   for SML and DSC layer (m)
        z_cld(pdims%i_start:pdims%i_end,                                       &
              pdims%j_start:pdims%j_end),                                      &
                        ! Cloud fraction weighted depth of cloud.
        z_cld_dsc(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end),                                  &
                        !   for SML and DSC layers (m)
        dscdepth(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end),                                   &
                        ! Depth of DSC layer (m)
        d_siems(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
        d_siems_dsc(pdims%i_start:pdims%i_end,                                 &
                    pdims%j_start:pdims%j_end),                                &
                        ! Siems (1990) et al. cloud-top entrainment
                        ! instability parm for SML and DSC inversions
        br_fback(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end),                                   &
        br_fback_dsc(pdims%i_start:pdims%i_end,                                &
                     pdims%j_start:pdims%j_end),                               &
                        ! Weight for degree of buoyancy reversal
                        ! feedback for SML and DSC inversions
        tv1_sd(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! Standard Deviation of level 1 Tv
        ft_nt_zh(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end),                                   &
                        ! FT_NT at ZH
        ft_nt_zhsc(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! FT_NT at ZHSC
        fq_nt_zh(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end),                                   &
                        ! FQ_NT at ZH
        fq_nt_zhsc(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! FQ_NT at ZHSC
        df_inv_sml(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! Radiative flux divergences
        df_inv_dsc(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        !   over inversion grid-level
        cf_sml(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! cloud fraction of SML
        cf_dsc(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! cloud fraction of DSC layer
        z_cf_base(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end),                                  &
                        ! cloud base height from cloud scheme
        z_ctop(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! cloud top height
        dqw_sml(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
                        ! QW and SL changes across SML disc inv
        dsl_sml(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &

        dqw_dsc(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
                        ! QW and SL changes across DSC disc inv
        dsl_dsc(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &

        rhokh_surf_ent(pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end),                             &
                        ! SML surf-driven entrainment KH
        rhokh_top_ent(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end),                              &
                        ! SML top-driven entrainment KH
        rhokh_dsct_ent(pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end),                             &
                        ! DSC top-driven entrainment KH
        zdsc_base(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end),                                  &
                        ! Height of base of K_top in DSC
        we_parm(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
                        ! Parametrised entrainment rates (m/s)
        we_dsc_parm(pdims%i_start:pdims%i_end,                                 &
                    pdims%j_start:pdims%j_end),                                &
                        !   for surf and DSC layers
        we_rho(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! rho*entrainment rate
        we_rho_dsc(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! rho*entrainment rate for DSC
        w_ls(pdims%i_start:pdims%i_end,                                        &
             pdims%j_start:pdims%j_end),                                       &
                        ! large-scale (subs) velocity
        w_ls_dsc(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end),                                   &
                        !   at subgrid inversion heights
        zh_np1(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
                        ! estimate of ZH at end of timestep
        zhsc_np1(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end),                                   &
                        ! estimate of ZHSC at end of timestep
        zh_frac(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end),                                    &
                        ! (ZH-ZHALF)/DZ
        zhsc_frac(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end),                                  &
                        ! (ZHSC-ZHALF)/DZ
        zrzi(pdims%i_start:pdims%i_end,                                        &
             pdims%j_start:pdims%j_end),                                       &
                        ! (z-z_base)/(z_i-z_base)
        zrzi_dsc(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end),                                   &
                        ! (z-z_base)/(z_i-z_base)
        t_frac(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end),                                     &
        t_frac_dsc(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! Fraction of timestep inversion is above
                        !   entr.t flux-level for SML and DSC layers
        svl_diff_frac(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end)
                        ! 1 - svl difference between ntdsc and ntml
                        ! divided by svl coupling threshold

real(kind=r_bl) ::                                                             &
 svl_plume,                                                                    &
                        ! SVL, SL and QW for a plume rising without
 sl_plume,                                                                     &
                        !   dilution from level 1.
 qw_plume

! Row arrays used to help vectorization of sub-sections 4. and 6.
real(kind=r_bl) ::                                                             &
dsldz(pdims%i_start:pdims%i_end),                                              &
                       ! Vertical gradient of SL in a well-mixed
                       ! layer.
cf_for_wb(pdims%i_start:pdims%i_end),                                          &
                       ! CF for use in wb calculation for decoupling
grad_t_adj_inv_rdz(pdims%i_start:pdims%i_end),                                 &
                       ! Temperature gradient adjustment times 1/RDZ
grad_q_adj_inv_rdz(pdims%i_start:pdims%i_end)
                       ! Humidity gradient adjustment times 1/RDZ

! For SCM diagnostics

integer ::  k_cloud_top(pdims%i_start:pdims%i_end,                             &
                        pdims%j_start:pdims%j_end),                            &
                        ! Level number of top of b.l. cloud.
            k_cloud_dsct(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end),                           &
                        ! Level number of top of dec. cloud.
            ntml_save(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end),                              &
                        ! Copy of NTML
            ntml_prev(pdims%i_start:pdims%i_end,                               &
                      pdims%j_start:pdims%j_end),                              &
                        ! NTML from previous timestep
            k_cbase(pdims%i_start:pdims%i_end,                                 &
                    pdims%j_start:pdims%j_end),                                &
                        ! grid-level above cloud-base
            dsc_removed(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                        ! Flag to indicate why EXCF_NL removed dsc layer

integer ::                                                                     &
 w_nonmono(pdims%i_start:pdims%i_end,                                          &
           pdims%j_start:pdims%j_end,bl_levels)
                        ! 0/1 flag for w being non-monotonic

! NEC vectorization
integer :: k_level(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end),                                 &
                        ! array to store level selection
           res_inv(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end)
                        ! Indicator of a resolved inversion

!  (b) Scalars.

integer ::                                                                     &
 k_cff
                        ! level counter for CFF
real(kind=r_bl) ::                                                             &
 virt_factor,                                                                  &
                       ! Temporary in calculation of buoyancy
                       ! parameters.
 dqw,                                                                          &
                       ! Total water content change across layer
                       ! interface.
 dsl,                                                                          &
                       ! Liquid/ice static energy change across
                       ! layer interface.
 dsl_ga,                                                                       &
                       ! As DSL but inc gradient adjustment
 dqw_ga,                                                                       &
                       ! As DQW but inc gradient adjustment
 dqcl,                                                                         &
                       ! Cloud liquid water change across layer
                       ! interface.
 dqcf,                                                                         &
                       ! Cloud frozen water change across layer
                       ! interface.
 q_vap_parc,                                                                   &
                       ! Vapour content of parcel
 q_liq_parc,                                                                   &
                       ! Condensed water content of parcel
 q_liq_env,                                                                    &
                       ! Condensed water content of environment
 t_parc,                                                                       &
                       ! Temperature of parcel
 t_dens_parc,                                                                  &
                       ! Density potential temperature of parcel
 t_dens_env,                                                                   &
                       ! Density potential temperature of
                       ! environment
 denv_bydz,                                                                    &
                       ! Gradient of density potential
                       ! temperature in environment
 dpar_bydz,                                                                    &
                       ! Gradient of density potential
                       ! temperature of parcel
 rho_dz,                                                                       &
                       ! rho*dz
 r_d_eta,                                                                      &
                       ! 1/(eta(k+1)-eta(k))
 svl_lapse,                                                                    &
                    ! Lapse rate of SVL above inversion (K/m)
 sl_lapse,                                                                     &
                    ! Lapse rate of SL above inversion (K/m)
 qw_lapse,                                                                     &
                    ! Lapse rate of QW above inversion (kg/kg/m)
 svl_lapse_base,                                                               &
                    ! Lapse rate of SVL above inversion (K/m)
 dsvl_top,                                                                     &
                    ! s_VL jump across inversion grid layer (K)
 tothf_efl,                                                                    &
                    ! total heat flux at entrainment flux grid-level
 totqf_efl,                                                                    &
                    ! Total QW flux at entrainment flux grid-level
 rhok_inv,                                                                     &
                    ! Estimate of rho*K within resolved inversion
 svl_lapse_rho,                                                                &
                    ! temporary in inversion flux calculation
 svl_target,                                                                   &
                    ! linear svl profile across resolved inversion
 recip_svl_lapse,                                                              &
                    ! one over the svl lapse rate
 Prandtl,                                                                      &
                    ! prandtl number
 ml_tend,                                                                      &
                    ! mixed layer tendency (d/dt)
 fa_tend,                                                                      &
                    ! free atmospheric tendency (d/dt)
 inv_tend,                                                                     &
                    ! limit on inversion grid-level tendency (d/dt)
 dflw_inv,                                                                     &
                    ! temporary in LW rad divergence calculation
 dfsw_inv,                                                                     &
                    ! temporary in SW rad divergence calculation
 z_rad_lim,                                                                    &
                    ! max height to search for peak LW cooling rate
 dz_disc_min,                                                                  &
                    ! smallest allowed DZ_DISC
 db_disc,                                                                      &
                    ! Temporary disc inversion buoyancy jump
 w_s_ent,                                                                      &
                    ! numerical (subsidence) entrainment rate
 dz_disc,                                                                      &
                    ! height of ZH below Z_uv(NTML+2)
 z_surf,                                                                       &
                    ! approx height of top of surface layer
 quad_a,                                                                       &
                    ! term `a' in quadratic solver for DZ_DISC
 quad_bm,                                                                      &
                    ! term `-b'in quadratic solver for DZ_DISC
 quad_c,                                                                       &
                    ! term `c' in quadratic solver for DZ_DISC
 w_m,                                                                          &
                    ! scaling velocity for momentum
 w_h,                                                                          &
                    ! scaling velocity for heat
 wstar3,                                                                       &
                    ! cube of convective velocity scale
 w_s_cubed,                                                                    &
                    ! convective velocity scale
 z_cbase,                                                                      &
                    ! cloud base height (m)
 zdsc_cbase,                                                                   &
                    ! DSC cloud base height (m)
 cfl_ml,cff_ml,                                                                &
                    ! liquid and frozen mixed layer fractions
 dfsw_top,                                                                     &
                    ! SW radiative flux change assoc with cloud-top
 c_ws,                                                                         &
                    ! Empirical constant multiplying Wstar
 c_tke,                                                                        &
                    ! Empirical constant in tke diagnostic
 pr_neut,                                                                      &
                    ! Neutral Prandtl number
 cu_depth_fac,                                                                 &
                    ! 0 to 1 factor related to cumulus depth
 w_curv,                                                                       &
                    ! curvature of w
 w_curv_nm,                                                                    &
 w_del_nm,                                                                     &
                    ! terms used to find where w is non-monotonic
 svl_diff           ! svl difference between ntdsc and ntml

real(kind=r_bl) :: rht_k, rht_kp1, rht_kp2
real(kind=r_bl) :: rht_max ( pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end )
real(kind=r_bl) :: interp, denom

integer ::                                                                     &
 i,                                                                            &
                 ! Loop counter (horizontal field index).
 j,                                                                            &
                 ! Offset counter in certain I loops.
 k,                                                                            &
                 ! Loop counter (vertical level index).
 kl,                                                                           &
                 ! K
 kp2,                                                                          &
                 ! K+2
 kp,km,                                                                        &
 kmax,                                                                         &
                 ! level of maximum
 k_rad_lim
                 ! limit on levels within which to search for
                 !   the max LW radiative cooling

integer :: i_wt   ! Water tracer counter

real(r_bl) :: w_var_inv ! vertical velocity variance at discontinuous inversions
real(r_bl) :: weight    ! Interpolation weight
real(r_bl) :: tke_nl_rh ! tke interpolated to rho-level
real(r_bl) :: delta_tke ! Increase in TKE at theta-levels due to inversion value


logical ::                                                                     &
 monotonic_inv,                                                                &
                 ! Flag that inversion grid-level properties are
                 ! monotonic, otherwise can't do subgrid
                 ! profile reconstruction (_DISC_INV flags)
 moisten(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! Indicator of whether inversion grid-level should
                 !   moisten this timestep (or dry)

! Flag for checking cloud-fraction at theta-level ntml+1
logical :: l_check_ntp1
integer :: p1

! Water tracer local arrays
integer, allocatable :: ntml_start(:,:)       ! Used to store ntml
integer, allocatable :: ntdsc_start(:,:)      ! Used to store ntdsc
logical, allocatable :: sml(:,:)              ! Dummy equivalent of 'dsc' field

! The following arrays are used to store details of the method used to solve
! for various quantities. These are used to update non-iso water tracers
! in a consistent way.
integer, allocatable :: dqw_sml_meth(:,:)
integer, allocatable :: dqw_dsc_meth(:,:)
integer, allocatable :: totqf_efl_meth1(:,:)
integer, allocatable :: totqf_efl_meth2(:,:)
logical, allocatable :: qw_lapse_zero_sml(:,:)
logical, allocatable :: qw_lapse_zero_dsc(:,:)

! Water tracer versions of water fields defined above
real(kind=r_bl), allocatable :: fsubs_wtrac(:,:,:,:)
real(kind=r_bl), allocatable :: fmic_wtrac(:,:,:,:)
real(kind=r_bl), allocatable :: dfmic_wtrac(:,:,:,:)
real(kind=r_bl), allocatable :: dfsubs_wtrac(:,:,:,:)
real(kind=r_bl), allocatable :: qls_inc_wtrac(:,:,:,:)
real(kind=r_bl), allocatable :: dqw_dsc_wtrac(:,:,:)
real(kind=r_bl), allocatable :: dqw_sml_wtrac(:,:,:)
real(kind=r_bl), allocatable :: fq_nt_zh_wtrac(:,:,:)
real(kind=r_bl), allocatable :: fq_nt_zhsc_wtrac(:,:,:)

real(kind=r_bl), allocatable :: z_uv_ntmlp1(:,:) ! Temporary copy of
                                                 ! z_uv at k=ntml+1

!Variables for cache-blocking
integer            :: jj         ! Block index

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate water tracer working arrays
if (l_wtrac) then
  allocate(ntml_start(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  allocate(ntdsc_start(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  allocate(fsubs_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                       bl_levels+1,n_wtrac))
  allocate(fmic_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
                      bl_levels+1,n_wtrac))
  allocate(dfmic_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                       bl_levels,n_wtrac))
  allocate(dfsubs_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
                        bl_levels,n_wtrac))
  allocate(qls_inc_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                         bl_levels,n_wtrac))
  allocate(dqw_dsc_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                         n_wtrac))
  allocate(dqw_sml_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                         n_wtrac))
  allocate(fq_nt_zh_wtrac(pdims%i_start:pdims%i_end,                           &
                          pdims%j_start:pdims%j_end, n_wtrac))
  allocate(fq_nt_zhsc_wtrac(pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end, n_wtrac))
  allocate(dqw_sml_meth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  allocate(dqw_dsc_meth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  allocate(totqf_efl_meth1                                                     &
                     (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  allocate(totqf_efl_meth2                                                     &
                     (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  allocate(qw_lapse_zero_sml                                                   &
                     (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  allocate(qw_lapse_zero_dsc                                                   &
                     (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
else
  ! These fields should always be allocated as they are passed to excf_nl_9c
  allocate(fq_nt_zh_wtrac(1,1,1))
  allocate(fq_nt_zhsc_wtrac(1,1,1))
end if

dz_disc_min = one_half * timestep * 1.0e-4_r_bl

!Start OpenMP parallel region

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP  private (i, j, jj, k, kp, kl, km, i_wt, w_curv_nm, w_del_nm, w_curv,    &
!$OMP  r_d_eta, rho_dz, z_surf, sl_plume, qw_plume,  q_liq_parc, q_liq_env,    &
!$OMP  t_parc, q_vap_parc, t_dens_parc, t_dens_env, dpar_bydz,  denv_bydz,     &
!$OMP  z_rad_lim,  k_rad_lim, dflw_inv, dfsw_inv, dfsw_top, svl_plume,         &
!$OMP  svl_diff, monotonic_inv, svl_lapse, svl_lapse_base,                     &
!$OMP  quad_a,  quad_bm, quad_c, dz_disc, dsvl_top, sl_lapse,  qw_lapse,       &
!$OMP  kp2, rht_k, rht_kp1, rht_kp2, interp, z_cbase )

!-----------------------------------------------------------------------
! Index to subroutine KMKHZ9C

! 1. Set up local variables, etc
! 2. Look for decoupled cloudy mixed-layer above SML top
! 3. Diagnose a discontinuous inversion structure.
! 4. Calculate the within-layer vertical gradients of cloud liquid
!      and frozen water
! 5. Calculate uniform mixed-layer cloud fractions and cloud depths
! 6. Calculate buoyancy flux factor used in the diagnosis of decoupling
! 7. Calculate inputs for the top of b.l. entrainment parametrization
! 8. Calculate the radiative flux change across cloud top
! 9. Calculate the non-turbulent fluxes at the layer boundaries.
! 10.Call subroutine EXCF_NL
!      - calculates parametrized entrainment rate, K profiles and
!        non-gradient flux/stress functions
! 11.Calculate "explicit" entrainment fluxes of SL and QW.

!-----------------------------------------------------------------------
! 1. Set up local variables, etc
!-----------------------------------------------------------------------
! 1.1 Calculate Z_TOP (top of levels) and NTML from previous timestep
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do j = jj, min((jj+bl_segment_size)-1,pdims%j_end)
    do i = pdims%i_start, pdims%i_end
      ntml_prev(i,j) = 1
    end do
  end do
  do k = 1, bl_levels-1
    do j = jj, min((jj+bl_segment_size)-1,pdims%j_end)
      do i = pdims%i_start, pdims%i_end
        z_top(i,j,k) = z_uv(i,j,k+1)
        !------------------------------------------------------------
        !find NTML from previous TS (for accurate gradient adjustment
        !of profiles - also note that NTML le BL_LEVELS-1)
        !------------------------------------------------------------
        if ( zh_prev(i,j)  >=  z_uv(i,j,k+1) ) ntml_prev(i,j)=k
      end do
    end do
  end do

  k = bl_levels
  do j = jj, min((jj+bl_segment_size)-1,pdims%j_end)
    do i = pdims%i_start, pdims%i_end
      z_top(i,j,k) = z_uv(i,j,k) + dzl(i,j,k)
    end do
  end do
end do ! jj
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 1.2 Calculate SVL: conserved buoyancy-like variable
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      sl(i,j,k)  = tl(i,j,k) + grcp * z_tq(i,j,k)
      svl(i,j,k) = sl(i,j,k) * ( one + c_virtual*qw(i,j,k) )
    end do
  end do
end do
!$OMP end do NOWAIT

!No halos
if (l_noice_in_turb) then
  ! use qsat_wat
  if ( l_mr_physics ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels
      call qsat_wat_mix(qs(:,:,k),t(:,:,k),p(:,:,k),pdims%i_end,pdims%j_end)
    end do
!$OMP end do
  else
!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels
      call qsat_wat(qs(:,:,k),t(:,:,k),p(:,:,k),pdims%i_end,pdims%j_end)
    end do
!$OMP end do
  end if ! l_mr_physics
else ! l_noice_in_turb
  if ( l_mr_physics ) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels
      call qsat_mix(qs(:,:,k),t(:,:,k),p(:,:,k),pdims%i_end,pdims%j_end)
    end do
!$OMP end do
  else
!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels
      call qsat(qs(:,:,k),t(:,:,k),p(:,:,k),pdims%i_end,pdims%j_end)
    end do
!$OMP end do
  end if ! l_mr_physics
end if ! test on l_noice_in_turb

!--------------------------------------------------------------------
! 1.3 Integrate non-turbulent increments to give flux profiles:
!     FT_NT, FQ_NT  are the flux profiles from non-turbulent processes
!                  (consisting of radiative FRAD, subsidence FSUBS and
!                   microphysical FMIC fluxes)
!--------------------------------------------------------------------
! For heat, units of rho * Km/s
! For humidity, units of rho * m/s
!----------------------------------
do k = 1, bl_levels
  km = max( 1, k-1 )
  kp = min( bl_levels, k+1 )
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end

    ! If the Intel Compiler vs 12 is used, a job fails with segmentation
    ! fault when the following loop is vectorised. Thus the following
    ! directive stop vectorization of this loop  when an Intel Compiler
    ! is used. Other compilers, for example Cray, should vectorise this
    ! loop automatically.
#if defined (IFORT_VERSION)
!DIR$ NOVECTOR
#endif

    do i = pdims%i_start, pdims%i_end
      w_grad(i,j,k) = (w(i,j,k)-w(i,j,km))*rdz(i,j,k)
      w_curv_nm = w(i,j,kp)-2.0_r_bl*w(i,j,k)+w(i,j,km)
      w_del_nm  = w(i,j,kp)-w(i,j,km)
      w_nonmono(i,j,k) = 0
      if ( abs(w_curv_nm) > abs(w_del_nm) ) w_nonmono(i,j,k) = 1
      sls_inc(i,j,k) = zero
      qls_inc(i,j,k) = zero
    end do ! i
  end do ! j
!$OMP end do NOWAIT
end do ! k

!$OMP BARRIER

!$OMP do SCHEDULE(DYNAMIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do k = 2, bl_levels-1
    do j = jj, min((jj+bl_segment_size)-1,pdims%j_end)
      do i = pdims%i_start, pdims%i_end

        if ( etadot(i,j,k)  < - tiny(one) .and.                                &
             etadot(i,j,k-1)< - tiny(one) ) then
            !-----------------------------------------------------------
            ! Only needed in subsidence regions
            ! Also don't attempt coupling with dynamics if w has
            ! significant vertical structure
            !-----------------------------------------------------------
          w_curv = (w_grad(i,j,k+1)-w_grad(i,j,k))/dzl(i,j,k)
          if (abs(w_curv) > 1.0e-6_r_bl .and. w_nonmono(i,j,k) == 1) then
              ! large curvature at a turning point
            sls_inc(i,j,k-1) = zero  ! need to make sure increments in
            qls_inc(i,j,k-1) = zero  ! level below are also set to zero
            etadot(i,j,k-1) = zero
            etadot(i,j,k)   = zero
            etadot(i,j,k+1) = zero
            w(i,j,k-1) = zero
            w(i,j,k)   = zero
            w(i,j,k+1) = zero
          else
            r_d_eta=1.0/(eta_theta_levels(k+1)-eta_theta_levels(k))
            sls_inc(i,j,k) = - etadot(i,j,k) * r_d_eta                         &
                                      * ( sl(i,j,k+1) - sl(i,j,k) )
            qls_inc(i,j,k) = - etadot(i,j,k) * r_d_eta                         &
                                      * ( qw(i,j,k+1) - qw(i,j,k) )
          end if  ! safe to calculate increments
        end if

      end do
    end do
  end do
end do !jj
!$OMP end do

! Repeat for necessary parts of last 2 loops for water tracers
if (l_wtrac) then
!$OMP do SCHEDULE(STATIC)
  do i_wt = 1, n_wtrac
    do k = 1, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          qls_inc_wtrac(i,j,k,i_wt) = zero
        end do ! i
      end do ! j
    end do ! k
  end do ! i_wt
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  do jj = pdims%j_start, pdims%j_end, bl_segment_size
    do k = 2, bl_levels-1
      do j = jj, min((jj+bl_segment_size)-1,pdims%j_end)
        do i = pdims%i_start, pdims%i_end

          if ( etadot(i,j,k)  < - tiny(one) .and.                              &
               etadot(i,j,k-1)< - tiny(one) ) then
            !-----------------------------------------------------------
            ! Only needed in subsidence regions
            ! Also don't attempt coupling with dynamics if w has
            ! significant vertical structure
            !-----------------------------------------------------------
            w_curv = (w_grad(i,j,k+1)-w_grad(i,j,k))/dzl(i,j,k)
            if (abs(w_curv) > 1.0e-6_r_bl .and. w_nonmono(i,j,k) == 1) then
              ! large curvature at a turning point
              do i_wt = 1, n_wtrac
                ! level below are also set to zero
                qls_inc_wtrac(i,j,k-1,i_wt) = zero
              end do
            else
              kp = k+1
              km = kp-1
              r_d_eta=one/(eta_theta_levels(kp)-eta_theta_levels(km))
              do i_wt = 1, n_wtrac
                qls_inc_wtrac(i,j,k,i_wt) = - etadot(i,j,k) * r_d_eta          &
                  * ( wtrac_bl(i_wt)%qw(i,j,kp) - wtrac_bl(i_wt)%qw(i,j,km) )
              end do
            end if  ! safe to calculate increments
          end if
        end do
      end do
    end do
  end do !jj
!$OMP end do
end if   ! l_wtrac

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
      ! Non-turbulent fluxes are defined relative to the surface
      ! so set them to zero at the surface
    frad(i,j,1)  = zero
    frad_lw(i,j,1) = zero
    frad_sw(i,j,1) = zero
    fsubs(i,j,1,1) = zero ! for heat
    fsubs(i,j,1,2) = zero ! for humidity
    fmic(i,j,1,1) = zero  ! for heat
    fmic(i,j,1,2) = zero  ! for humidity

    ! The two expressions are maintained as comment to highlight that
    ! ft_nt and fq_net need to be updated if there are any eventual
    ! changes in fmic/frad/fsubs at level 1.  A job fails with
    ! segmentation fault when the Intel compiler is used if the two
    ! arrays are not set to 0.0 directly at level 1.

    !ft_nt(i,j,1)  = frad(i,j,1) + fmic(i,j,1,1) + fsubs(i,j,1,1)
    !fq_nt(i,j,1)  =               fmic(i,j,1,2) + fsubs(i,j,1,2)

    ft_nt(i,j,1)  = zero
    fq_nt(i,j,1)  = zero

  end do
end do
!$OMP end do NOWAIT

! Repeat for necessary parts of last loop for water tracers
if (l_wtrac) then
!$OMP do SCHEDULE(STATIC)
  do i_wt = 1, n_wtrac
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        ! Non-turbulent fluxes are defined relative to the surface
        ! so set them to zero at the surface
        fsubs_wtrac(i,j,1,i_wt)     = zero
        fmic_wtrac(i,j,1,i_wt)      = zero
        wtrac_bl(i_wt)%fq_nt(i,j,1) = zero
      end do
    end do
  end do
!$OMP end do
end if

! This is the most computational expensive loop of the subroutine. The
! following parallelisation obtains lower times than the ones obtained
! by using the jj/bl_segment_size technique.  While it would be
! possible to remove the "df_over_cp", "dfmic", and "dfsubs" arrays,
! without them some jobs lose bit-comparability with the Intel
! compiler.

do k = 1, bl_levels
  kp = k+1
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      rho_dz = rho_mix_tq(i,j,k) * dzl(i,j,k)

      dflw_over_cp(i,j,k) = - rad_hr(i,j,1,k) * rho_dz
      dfsw_over_cp(i,j,k) = - rad_hr(i,j,2,k) * rho_dz
      df_over_cp(i,j,k)   = dflw_over_cp(i,j,k) + dfsw_over_cp(i,j,k)

      dfmic(i,j,k,1)  = - micro_tends(i,j,1,k) * rho_dz
      dfmic(i,j,k,2)  = - micro_tends(i,j,2,k) * rho_dz
      dfsubs(i,j,k,1) = - sls_inc(i,j,k)       * rho_dz
      dfsubs(i,j,k,2) = - qls_inc(i,j,k)       * rho_dz

      frad(i,j,kp)   = frad(i,j,k)    + df_over_cp(i,j,k)
      frad_lw(i,j,kp)= frad_lw(i,j,k) + dflw_over_cp(i,j,k)
      frad_sw(i,j,kp)= frad_sw(i,j,k) + dfsw_over_cp(i,j,k)
      fsubs(i,j,kp,1)= fsubs(i,j,k,1) + dfsubs(i,j,k,1)
      fsubs(i,j,kp,2)= fsubs(i,j,k,2) + dfsubs(i,j,k,2)
      fmic(i,j,kp,1) = fmic(i,j,k,1)  + dfmic(i,j,k,1)
      fmic(i,j,kp,2) = fmic(i,j,k,2)  + dfmic(i,j,k,2)

      ft_nt(i,j,kp)  = frad(i,j,kp) + fmic(i,j,kp,1) + fsubs(i,j,kp,1)
      fq_nt(i,j,kp)  =                fmic(i,j,kp,2) + fsubs(i,j,kp,2)
    end do ! i
  end do !j
!$OMP end do NOWAIT
end do ! k

! Repeat necessary parts of last loop for water tracer
if (l_wtrac) then
!$OMP do SCHEDULE(STATIC)
  do i_wt = 1, n_wtrac
    do k = 1, bl_levels
      do j = pdims%j_start, pdims%j_end
        kp = k+1
        do i = pdims%i_start, pdims%i_end
          rho_dz = rho_mix_tq(i,j,k) * dzl(i,j,k)
          dfmic_wtrac(i,j,k,i_wt)  =                                           &
                    - wtrac_as(i_wt)%micro_tends(i,j,k) * rho_dz
          dfsubs_wtrac(i,j,k,i_wt) = - qls_inc_wtrac(i,j,k,i_wt) * rho_dz

          fsubs_wtrac(i,j,kp,i_wt)=                                            &
                   fsubs_wtrac(i,j,k,i_wt) + dfsubs_wtrac(i,j,k,i_wt)
          fmic_wtrac(i,j,kp,i_wt) =                                            &
                   fmic_wtrac(i,j,k,i_wt)  + dfmic_wtrac(i,j,k,i_wt)
          wtrac_bl(i_wt)%fq_nt(i,j,kp)  =                                      &
                   fmic_wtrac(i,j,kp,i_wt) + fsubs_wtrac(i,j,kp,i_wt)
        end do ! i
      end do !j
    end do ! k
  end do   ! i_wt
!$OMP end do NOWAIT
end if   ! l_wtrac

!-----------------------------------------------------------------------
! 1.4 Set UNSTABLE flag and find first level above surface layer
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(DYNAMIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do j = jj, min(jj+bl_segment_size-1,pdims%j_end)
    do i = pdims%i_start, pdims%i_end
      unstable(i,j) = (fb_surf(i,j) >  zero)
      k_plume(i,j)  = -1
    end do
  end do

      !------------------------------------------------------------
      ! Find grid-level above top of surface layer, taken
      ! to be at a height, z_surf, given by:
      !       Z_SURF = 0.1*ZH_PREV
      ! Use ZH_prev since that will have determined the shape
      ! of the time-level n profiles.
      !------------------------------------------------------------

  do k = 1, bl_levels-1
    do j = jj, min(jj+bl_segment_size-1,pdims%j_end)
      do i = pdims%i_start, pdims%i_end

        if ( unstable(i,j) ) then

          z_surf = 0.1_r_bl * zh_prev(i,j)
          if ( z_tq(i,j,k) >= z_surf .and. k_plume(i,j) == -1 ) then
                 !reached z_surf
            k_plume(i,j)=k
          end if
          if ( svl(i,j,k+1) >= svl(i,j,k)                                      &
                  .and. k_plume(i,j) == -1 ) then
                 !reached inversion
            k_plume(i,j)=k
          end if

        end if
      end do
    end do
  end do

  do j = jj, min(jj+bl_segment_size-1,pdims%j_end)
    do i = pdims%i_start, pdims%i_end
      if (k_plume(i,j) == -1) k_plume(i,j)=1
    end do
  end do
end do ! jj
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 2.  Look for decoupled cloudy mixed-layer above SML top
!     Note that sc_diag_all_rh_max skips (2.1) and diagnoses a decoupled
!     layer in the same way for all regimes (2.2), irrespective of the
!     cumulus diagnosis.
!-----------------------------------------------------------------------
! 2.1  (if not CUMULUS: starting from level 3 and below 2.5km):
!      find cloud-base above SML inversion, ie. above NTML+1,
!      then cloud-top (ie. CF < SC_CFTOL)
!      and finally check that cloud is well-mixed.
!-----------------------------------------------------------------------
!      Initialise variables

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    cloud_base(i,j) = .false.
    dsc(i,j) = .false.
    coupled(i,j) = .false.
    zhsc(i,j)    = zero
    ntdsc(i,j)   = 0
  end do
end do
!$OMP end do

if ( .not. sc_diag_opt == sc_diag_all_rh_max ) then

  ! Use svl-gradient method to diagnose the Sc-top level at non-cumulus points

!$OMP do SCHEDULE(DYNAMIC)
  do jj = pdims%j_start, pdims%j_end, bl_segment_size
    do k = 3, bl_levels-1
      do j = jj, min(jj+bl_segment_size-1,pdims%j_end)
        do i = pdims%i_start, pdims%i_end

          !-----------------------------------------------------------------
          !..Find cloud-base (where cloud here means CF > SC_CFTOL)
          !-----------------------------------------------------------------

          if ( .not. cumulus(i,j) .and.                                        &
               z_tq(i,j,k) < zmaxb_for_dsc .and.                               &
               k  >   ntml(i,j)+1 .and. cf(i,j,k)  >   sc_cftol                &
                                .and. .not. cloud_base(i,j)                    &
            !                                  not yet found cloud-base
                                          .and. .not. dsc(i,j) ) then
            !                                  not yet found a Sc layer
            cloud_base(i,j) = .true.
          end if
          if ( cloud_base(i,j) .and. .not. dsc(i,j) .and.                      &
            !                 found cloud-base but not yet reached cloud-top
                         cf(i,j,k+1) < sc_cftol .and.                          &
                         z_tq(i,j,k) < zmaxt_for_dsc                           &
            !                 got to cloud-top below ZMAXT_FOR_DSC
                       ) then
            cloud_base(i,j) = .false.         ! reset CLOUD_BASE
            !-----------------------------------------------------------
            ! Look to see if at least top of cloud is well mixed:
            ! test SVL-gradient for top 2 pairs of levels, in case
            ! cloud top extends into the inversion.
            ! Parcel descent in Section 4.0 below will determine depth
            ! of mixed layer.
            !----------------------------------------------------------
            if ( (svl(i,j,k)-svl(i,j,k-1))                                     &
                         /(z_tq(i,j,k)-z_tq(i,j,k-1))                          &
                                                <   max_svl_grad ) then
              dsc(i,j) = .true.
              ntdsc(i,j) = k
              zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
            else if ( (svl(i,j,k-1)-svl(i,j,k-2))                              &
                          /(z_tq(i,j,k-1)-z_tq(i,j,k-2))                       &
                                              <   max_svl_grad ) then
              !---------------------------------------------------------
              ! Well-mixed layer with top at k-1 or k.  Check whether
              ! there is a buoyancy inversion between levels k-1 and k
              ! in a manner similar to the surface-driven plume: compare
              ! the buoyancy gradient between levels K-1 and K for an
              ! undiluted parcel and the environment
              !---------------------------------------------------------
              sl_plume = tl(i,j,k-1) + grcp * z_tq(i,j,k-1)
              qw_plume = qw(i,j,k-1)
              ! ------------------------------------------------------------
              ! calculate parcel water by linearising qsat about the
              ! environmental temperature.
              ! ------------------------------------------------------------
              if (t(i,j,k) >  tm) then
                q_liq_parc = max( zero, ( qw_plume - qs(i,j,k) -               &
                  dqsdt(i,j,k)*                                                &
                  ( sl_plume-grcp*z_tq(i,j,k)-t(i,j,k) )                       &
                                       ) *a_qs(i,j,k) )
                q_liq_env = max( zero, ( qw(i,j,k) - qs(i,j,k)                 &
                            -dqsdt(i,j,k)*( tl(i,j,k) - t(i,j,k) )             &
                                       ) *a_qs(i,j,k) )
                ! add on the difference in the environment's ql as
                ! calculated by the partial condensation scheme (using
                ! some RH_CRIT value) and what it would be if
                ! RH_CRIT=1. This then imitates partial condensation in
                ! the parcel.
                q_liq_parc = q_liq_parc + qcl(i,j,k) + qcf(i,j,k)              &
                               - q_liq_env
                t_parc = sl_plume - grcp * z_tq(i,j,k) +                       &
                                 lcrcp*q_liq_parc
              else
                q_liq_parc = max( zero, ( qw_plume - qs(i,j,k) -               &
                  dqsdt(i,j,k)*                                                &
                    ( sl_plume - grcp*z_tq(i,j,k)-t(i,j,k) )                   &
                                       ) *a_qs(i,j,k) )
                q_liq_env = max( zero, ( qw(i,j,k) - qs(i,j,k)                 &
                   -dqsdt(i,j,k)*( tl(i,j,k) - t(i,j,k) )                      &
                                       ) *a_qs(i,j,k) )
                ! add on difference in environment's ql between RH_CRIT and
                ! RH_CRIT=1
                q_liq_parc = q_liq_parc + qcl(i,j,k) + qcf(i,j,k)              &
                               - q_liq_env
                t_parc = sl_plume - grcp * z_tq(i,j,k) +                       &
                                 lsrcp*q_liq_parc
              end if
              q_vap_parc=qw_plume - q_liq_parc

              t_dens_parc=t_parc*(one+c_virtual*q_vap_parc-q_liq_parc)
              t_dens_env=t(i,j,k)*                                             &
                         (one+c_virtual*q(i,j,k)-qcl(i,j,k)-qcf(i,j,k))
              ! find vertical gradients in parcel and environment SVL
              ! (using values from level below (K-1))
              env_svl_km1(i,j) = t(i,j,k-1) * ( one+c_virtual*q(i,j,k-1)       &
                   -qcl(i,j,k-1)-qcf(i,j,k-1) ) + grcp*z_tq(i,j,k-1)
              dpar_bydz=(t_dens_parc+grcp*z_tq(i,j,k)-                         &
                          env_svl_km1(i,j)) /                                  &
                      (z_tq(i,j,k)-z_tq(i,j,k-1))
              denv_bydz=(t_dens_env+grcp*z_tq(i,j,k)-                          &
                          env_svl_km1(i,j))/                                   &
                      (z_tq(i,j,k)-z_tq(i,j,k-1))

              if ( denv_bydz >  1.25_r_bl*dpar_bydz ) then
                ! there is an inversion between levels K-1 and K
                if ( k  >=  ntml(i,j)+3 ) then
                  ! if NTDSC == NTML+1 then assume we're looking
                  ! at the same inversion and so don't set DSC
                  ntdsc(i,j) = k-1
                  zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
                  dsc(i,j) = .true.
                end if
              else
                ! no inversion between levels K-1 and K, assume there
                ! is an inversion between K and K+1 because of CF change
                ntdsc(i,j) = k
                zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
                dsc(i,j) = .true.
              end if
            end if
          end if
        end do ! i
      end do ! j
    end do ! k
  end do ! jj
!$OMP end do

end if  ! ( .not. sc_diag_opt == sc_diag_all_rh_max )

!-----------------------------------------------------------------------
! 2.2 If using sc_diag_all_rh_max, this is where a decoupled cloud layer
!       is diagnosed, irrespective of a cumulus diagnosis.
!     Otherwise, if the layer to ZHPAR is a cumulus layer capped by cloud and
!       an inversion, declare this layer a decoupled cloud layer and
!       set ZHSC and NTDSC accordingly.
!-----------------------------------------------------------------------
if ( sc_diag_opt == sc_diag_cu_rh_max .or.                                     &
     sc_diag_opt == sc_diag_all_rh_max ) then

  ! Options that diagnose the Sc-top using the max total-water RH method...

  ! j-loop outermost to allow parallelisation (k-loop is sequential)
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end

    do i = pdims%i_start, pdims%i_end
      ! Initialise max RH in the column to zero
      rht_max(i,j) = zero
    end do

    do k = 1, bl_levels-1
      do i = pdims%i_start, pdims%i_end
        if ( ( cumulus(i,j) .and. z_tq(i,j,k) > z_lcl(i,j) ) .or.              &
             ( sc_diag_opt == sc_diag_all_rh_max .and. (.not. cumulus(i,j))    &
               .and. k > ntml(i,j)+1 ) ) then
          ! If sc_diag_opt == sc_diag_cu_rh_max, only check cumulus points
          ! at heights above the LCL.
          ! If sc_diag_opt == sc_diag_all_rh_max, also check non-cumulus points
          ! at heights above ntml+1.

          ! Find max of RHt = qw/qsat(Tl) in the column
          ! (with extra check that value exceeds next level, to exclude
          !  points at bl_levels-1 where a larger value occurs at bl_levels)
          if ( cf(i,j,k) > sc_cftol .and. z_tq(i,j,k) < zmaxt_for_dsc ) then
            if ( l_mr_physics ) then
              call qsat_wat_mix( rht_k,   tl(i,j,k),   p(i,j,k) )
              call qsat_wat_mix( rht_kp1, tl(i,j,k+1), p(i,j,k+1) )
            else
              call qsat_wat( rht_k,   tl(i,j,k),   p(i,j,k) )
              call qsat_wat( rht_kp1, tl(i,j,k+1), p(i,j,k+1) )
            end if
            rht_k   = qw(i,j,k)   / rht_k
            rht_kp1 = qw(i,j,k+1) / rht_kp1
            if ( rht_k > rht_max(i,j) .and. rht_k > rht_kp1 ) then
              ntdsc(i,j) = k
              rht_max(i,j) = rht_k
            end if
          end if  ! ( cf(i,j,k) > sc_cftol etc )
        end if  ! ( cumulus(i,j) etc )
      end do  ! i = pdims%i_start, pdims%i_end
    end do  ! k = 1, bl_levels-1

    do i = pdims%i_start, pdims%i_end
      if ( ntdsc(i,j) > 0 .and. ( .not. dsc(i,j) ) ) then
        ! If we just found ntdsc above (exclude points where dsc already
        ! set true by finding ntdsc at non-cumulus points earlier under
        ! the option sc_diag_opt = sc_diag_cu_rh_max)...
        ! Set flag indicating we found a potential Sc layer top
        dsc(i,j) = .true.
        ! Now we assume that theta-level ntdsc is wholly within the Sc-layer,
        ! and that theta-level ntdsc+1 is composed of air with RH of
        ! level ntdsc up to height zhsc, and air with RH of level ntdsc+2
        ! above that height.  Interpolate to find the height of zhsc
        ! (fraction of theta-level ntdsc+1 that is below the Sc-top).
        k = ntdsc(i,j)
        kp2 = min(k+2, bl_levels)  ! Avoid out-of-bounds when k = bl_levels-1
        rht_k = rht_max(i,j)
        if ( l_mr_physics ) then
          call qsat_wat_mix( rht_kp1, tl(i,j,k+1), p(i,j,k+1) )
          call qsat_wat_mix( rht_kp2, tl(i,j,kp2), p(i,j,kp2) )
        else
          call qsat_wat( rht_kp1, tl(i,j,k+1), p(i,j,k+1) )
          call qsat_wat( rht_kp2, tl(i,j,kp2), p(i,j,kp2) )
        end if
        rht_kp1 = qw(i,j,k+1) / rht_kp1
        rht_kp2 = qw(i,j,kp2) / rht_kp2
        if ( rht_kp2 < rht_kp1 ) then
          ! RHt(k+1) lies between RHt(k) and RHt(k+2); compute fraction
          interp = (rht_kp1 - rht_kp2) / (rht_k - rht_kp2)
          zhsc(i,j) = (one-interp) * z_uv(i,j,k+1)                             &
                    +      interp  * z_uv(i,j,kp2)
        else
          ! Rht(k+1) is a local minimum; can't construct k+1 as a fraction
          ! of properties from k and k+2.  Set zhsc to what it would be in
          ! the limit RHt(k+2) = RHt(k+1)
          zhsc(i,j) = z_uv(i,j,k+1)
        end if
      end if ! ( ntdsc(i,j) > 0 )
    end do
  end do
!$OMP end do NOWAIT

else if ( sc_diag_opt == sc_diag_cu_relax ) then

  ! Diagnosed simply if significant cloud fraction at ZHPAR
  ! below the height threshold zmaxt_for_dsc

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      k = ntpar(i,j)
      if ( cumulus(i,j) .and. k < bl_levels  ) then
              ! cumulus layer within BL_LEVELS
        if ( z_tq(i,j,k) < zmaxt_for_dsc .and.                                 &
              ! cloud top below zmaxt_for_dsc
             ( max( cf(i,j,k-1),cf(i,j,k),cf(i,j,k+1) ) > sc_cftol )           &
              ! cloud-top sufficiently cloudy
            ) then
          dsc(i,j)  = .true.
          zhsc(i,j) = zhpar(i,j)
          ntdsc(i,j)= ntpar(i,j)
        end if
      end if
    end do
  end do
!$OMP end do NOWAIT

else if ( sc_diag_opt == sc_diag_orig ) then

  ! Original code, only diagnosed if shallow cu or not l_param_conv

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( (l_param_conv .and.                                                 &
              l_shallow(i,j) .and. ntpar(i,j)  <   bl_levels )                 &
              ! shallow cumulus layer within BL_LEVELS
           .or. (.not. l_param_conv .and.                                      &
                cumulus(i,j) .and. ntpar(i,j)  <   bl_levels ) ) then
              ! cumulus layer and inversion found
        if ( cf(i,j,ntpar(i,j))  >   sc_cftol  .or.                            &
             cf(i,j,ntpar(i,j)+1)  >   sc_cftol ) then
           ! cloudy
          dsc(i,j)  = .true.
          zhsc(i,j) = zhpar(i,j)
          ntdsc(i,j)= ntpar(i,j)
        end if
      end if
    end do
  end do
!$OMP end do NOWAIT

end if  ! test on sc_diag_opt

if ( l_use_sml_dsc_fixes ) then
  ! Need to override "NOWAIT" on the previous blocks if going in here:
!$OMP BARRIER
  ! If conv_diag has diagnosed a SML rising significantly above cloud-base,
  ! abort any DSC diagnosis higher-up.
  ! This is because at present, diagnosing an elevated DSC-layer prompts
  ! excf_nl_9c to test whether that DSC-layer should be recoupled with the SML,
  ! instead of testing whether the cloudy top of the SML should be decoupled.
  ! Skipping the decoupling test for the SML means we can force mixing up
  ! to cloud-top in cloudy boundary-layers which are really too stable to
  ! be well-mixed!
  ! The code currently only permits us to have one DSC-layer in the column,
  ! so if we have a cloudy boundary-layer, we need to reserve the one possible
  ! DSC for the case where the cloudy SML-top decouples to make a DSC.
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( unstable(i,j) .and. (.not. cumulus(i,j)) .and. dsc(i,j) ) then
        if ( zh(i,j) - z_lcl(i,j) >= 400.0  ) then
          ! Using cloud-layer depth >= 400m consistent with cumulus_test
          ! in conv_diag, as this indicates that conv_diag was prevented
          ! from diagnosing cumulus by its test on whether the qw profile
          ! looks well-mixed.  But it doesn't test on stability, so it can
          ! give a very misleading diagnosis of a well-mixed cloudy layer.
          ! In particular, it compares the cloud-layer qw gradient with the
          ! sub-cloud-layer qw gradient, giving a well-mixed diagnosis if
          ! they are similar even if both are large!
          dsc(i,j) = .false.
          ntdsc(i,j) = 0
        end if
      end if
    end do
  end do
!$OMP end do NOWAIT
end if  ! ( l_use_sml_dsc_fixes )

!-----------------------------------------------------------------------
! 2.3 Calculate the radiative flux changes across cloud top for the
!      stratocumulus layer and thence a first guess for the top-down
!      mixing depth of this layer, DSCDEPTH.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    k_cloud_dsct(i,j) = 0
    df_dsct_over_cp(i,j) = zero
    df_inv_dsc(i,j) = zero
  end do
end do
!$OMP end do

if (l_new_kcloudtop) then
  !---------------------------------------------------------------------
  ! improved method of finding the k_cloud_dsct, the top of the mixed
  ! layer as seen by radiation
  !---------------------------------------------------------------------
!$OMP do SCHEDULE(DYNAMIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
        !-------------------------------------------------------------
        ! Find k_cloud_dsct as the equivalent to ntdsc (the top level
        ! of the DSC) seen by radiation, which we take as the level two
        ! levels below the lowest level with free tropospheric radiative
        ! cooling.  This is done by finding the level with maximum LW
        ! cooling, below z_rad_lim and above the SML and 0.5*ZHSC
        ! (ie, restrict search to `close' to ZHSC)
        ! Necessary as radiation is not usually called every timestep.
        !-------------------------------------------------------------
      if ( dsc(i,j) .and. ntdsc(i,j)+2 <= bl_levels ) then
        z_rad_lim = max( z_tq(i,j,ntdsc(i,j)+2)+0.1_r_bl, 1.2_r_bl*zhsc(i,j) )

        k = ntml(i,j)+2
        do while (z_tq(i,j,k) < z_rad_lim .and. k < bl_levels)
          if ( z_tq(i,j,k) > one_half*zhsc(i,j)                                &
              .and. dflw_over_cp(i,j,k) > df_dsct_over_cp(i,j) ) then
            k_cloud_dsct(i,j) = k
            df_dsct_over_cp(i,j) = dflw_over_cp(i,j,k)
          end if
          k = k+1
        end do ! k
        ! Set K_CLOUD_DSCT to the level below if DF in the level
        ! above is less than 1.5 times the level above that
        ! (implies K_CLOUD_DSCT+1 is typical of free trop so
        !  K_CLOUD_DSCT must be inversion level, instead of ntdsc).
        ! DF in level K_CLOUD_DSCT+1 is then included as DF_INV_DSC
        ! (see below).
        k = k_cloud_dsct(i,j)
        if ( k >  1 .and. k < bl_levels -1 ) then
          if (dflw_over_cp(i,j,k+1) < 1.5_r_bl*dflw_over_cp(i,j,k+2))          &
              k_cloud_dsct(i,j) = k-1
        end if
      end if  ! DSC test separated out

    end do ! i
    !-----------------------------------------------------------------
    !  Find bottom grid-level (K_LEVEL) for cloud-top radiative flux
    !  divergence: higher of base of LW radiatively cooled layer,
    !  ZH and 0.5*ZHSC, since cooling must be in upper part of layer
    !  in order to generate turbulence.
    !-----------------------------------------------------------------
    do i = pdims%i_start, pdims%i_end
      k_level(i,j) = k_cloud_dsct(i,j)
      if ( k_cloud_dsct(i,j)  >   1 ) then
        k_rad_lim = ntml(i,j)+1
        k=k_cloud_dsct(i,j)-1
        kl=max(1,k)  ! only to avoid out-of-bounds compiler warning
        do while ( k  >   k_rad_lim                                            &
                .and. dflw_over_cp(i,j,kl)  >   zero                           &
                .and. z_tq(i,j,kl)  >   one_half*zhsc(i,j) )
          k_level(i,j) = k
          k = k-1
          kl=max(1,k)
        end do ! k
      end if
    end do ! i
  end do ! j
!$OMP end do

else

  ! original method of finding the k_cloud_dsct, the top of the mixed layer
  ! as seen by radiation, found to be resolution dependent

!$OMP do SCHEDULE(DYNAMIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      !-------------------------------------------------------------
      ! Find the layer with the greatest LW radiative flux jump in
      ! the upper half of the boundary layer and assume that this
      ! marks the top of the DSC layer.
      ! Necessary as radiation is not usually called every timestep.
      !-------------------------------------------------------------
      ! Limit the search to above the SML.
      k_rad_lim = ntml(i,j)+2

      do k = max(1,k_rad_lim), min(bl_levels,ntdsc(i,j)+2)

        if ( dsc(i,j) .and. z_tq(i,j,k) > one_half*zhsc(i,j)                   &
              .and. dflw_over_cp(i,j,k) > df_dsct_over_cp(i,j) ) then
          k_cloud_dsct(i,j) = k
            ! Set K_CLOUD_DSCT to the level below if its DF is greater
            ! than half the maximum.  DF in level K_CLOUD_DSCT+1 is then
            ! included as DF_INV_DSC below.
          if (dflw_over_cp(i,j,k-1)  >   one_half*dflw_over_cp(i,j,k))         &
             k_cloud_dsct(i,j) = k-1
          df_dsct_over_cp(i,j) = dflw_over_cp(i,j,k)
        end if

      end do ! k
    end do ! i

      !-----------------------------------------------------------------
      !  Find bottom grid-level (K_LEVEL) for cloud-top radiative flux
      !  divergence: higher of base of LW radiatively cooled layer,
      !  ZH and 0.5*ZHSC, since cooling must be in upper part of layer
      !  in order to generate turbulence.
      !-----------------------------------------------------------------

    do i = pdims%i_start, pdims%i_end
      k_level(i,j) = k_cloud_dsct(i,j)
      if ( k_cloud_dsct(i,j)  >   1 ) then
        k_rad_lim = ntml(i,j)+1
        k=k_cloud_dsct(i,j)-1
        kl=max(1,k)  ! only to avoid out-of-bounds compiler warning
        do while ( k  >   k_rad_lim                                            &
                  .and. dflw_over_cp(i,j,kl)  >   zero                         &
                  .and. z_tq(i,j,kl)  >   one_half*zhsc(i,j) )
          k_level(i,j) = k
          k = k-1
          kl=max(1,k)
        end do ! k
      end if
    end do ! i
  end do ! j
!$OMP end do

end if  ! test on l_new_kcloudtop

      !-----------------------------------------------------------------
      ! Calculate LW and SW flux divergences and combine into
      ! cloud-top turbulence forcing.
      ! Need to account for radiative divergence in cloud in inversion
      ! grid-level, DF_INV_DSC. Assume DF_OVER_CP(K_cloud_dsct+2) is
      ! representative of clear-air rad divergence and so subtract this
      ! `clear-air' part from the grid-level divergence.
      !-----------------------------------------------------------------

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    if ( k_cloud_dsct(i,j)  >   0 ) then
      dflw_inv = zero
      dfsw_inv = zero
      if ( k_cloud_dsct(i,j) <  bl_levels ) then
        k = k_cloud_dsct(i,j)+1
        if ( k  <   bl_levels ) then
          dflw_inv = dflw_over_cp(i,j,k)                                       &
                     - dflw_over_cp(i,j,k+1)                                   &
                            * dzl(i,j,k)/dzl(i,j,k+1)
          dfsw_inv = dfsw_over_cp(i,j,k)                                       &
                     - dfsw_over_cp(i,j,k+1)                                   &
                            * dzl(i,j,k)/dzl(i,j,k+1)
        else
          dflw_inv = dflw_over_cp(i,j,k)
          dfsw_inv = dfsw_over_cp(i,j,k)
        end if
        dflw_inv = max( dflw_inv, zero )
        dfsw_inv = min( dfsw_inv, zero )
      end if
      df_inv_dsc(i,j) = dflw_inv + dfsw_inv

      df_dsct_over_cp(i,j) = frad_lw(i,j,k_cloud_dsct(i,j)+1)                  &
                           - frad_lw(i,j,k_level(i,j))                         &
                           + dflw_inv

      dfsw_top = frad_sw(i,j,k_cloud_dsct(i,j)+1)                              &
               - frad_sw(i,j,k_level(i,j))                                     &
               + dfsw_inv

        !-----------------------------------------------------------
        ! Combine SW and LW cloud-top divergences into a net
        ! divergence by estimating SW flux divergence at a given
        ! LW divergence = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
        ! Empirically (from LEM data) a reasonable fit is found
        ! with A small and (1-exp{-A*kappa_sw/kappa_lw}) = dfsw_frac
        !-----------------------------------------------------------
      df_dsct_over_cp(i,j) = max( zero,                                        &
                    df_dsct_over_cp(i,j) + dfsw_frac * dfsw_top )
    end if
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
! 2.4 Set NBDSC, the bottom level of the DSC layer.
!     Note that this will only be used to give an estimate of the layer
!     depth, DSCDEPTH, used to calculate the entrainment
!     rate (the dependence is only weak), and that a more accurate
!     algorithm is subsequently used to determine the depth over which
!     the top-down mixing profiles will be applied.  If DSC is false,
!     DSCDEPTH = 0.  The plume descent here uses a radiative
!     perturbation to the cloud-layer SVL (use level NTDSC-1 in case
!     SVL is not yet well-mixed to NTDSC), based roughly
!     on a typical cloud-top residence time.  If the plume does not sink
!     and the cloud is decoupled from the surface (ie. above Alan
!     Grant's ZH), then it is assumed to be stable, ie. St rather than
!     Sc, and no mixing or entrainment is applied to it.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(DYNAMIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do j = jj, min(jj+bl_segment_size-1,pdims%j_end)
    do i = pdims%i_start, pdims%i_end
      nbdsc(i,j) = ntdsc(i,j)+1
      if (dsc(i,j)) then
        ! The depth of the radiatively-cooled layer tends to be less
        ! than O(50m) and so RAD_HR will be an underestimate of the
        ! cooling tendency there.  Compensate by multiplying by
        ! DZL/50. (~4) Recall that DF_OVER_CP(I,j,K) = RAD_HR *
        ! RHO_MIX_TQ * DZL Thus use cloud-top radiative forcing as
        ! follows:

        k = ntdsc(i,j)
        rho_dz = rho_mix_tq(i,j,k) * dzl(i,j,k)
        svl_plume=svl(i,j,k-1)                                                 &
           - ct_resid * dzl(i,j,k)*df_dsct_over_cp(i,j) / ( 50.0_r_bl*rho_dz )
      else
        svl_plume=zero
      end if
      do k = min(bl_levels-1, ntdsc(i,j)-1), 1, -1
        if (svl_plume  <   svl(i,j,k) ) then
          nbdsc(i,j) = k+1     ! marks lowest level within ML
        end if
      end do ! k
    end do ! i
  end do ! j
end do !jj
!$OMP end do
!----------------------------------------------------------------------
! 2.5 Tidy up variables associated with decoupled layer
!       NOTE that NTDSC ge 3 if non-zero
!----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ! Note that ZHSC-Z_UV(NTML+2) may = 0, so this test comes first!
    if (cumulus(i,j) .and. dsc(i,j))                                           &
                     nbdsc(i,j) = max( nbdsc(i,j), ntml(i,j)+2 )
    if ( ntdsc(i,j)  >=  1 ) then
      if ( nbdsc(i,j) <   ntdsc(i,j)+1 ) then
        if (sc_diag_opt==sc_diag_orig .or. sc_diag_opt==sc_diag_cu_relax) then
          ! Initial zhsc is normally on rho-level ntdsc+1
          dscdepth(i,j) = z_uv(i,j,ntdsc(i,j)+1) - z_uv(i,j,nbdsc(i,j))
        else
          ! Initial zhsc was interpolated between levels; use accurate zhsc
          dscdepth(i,j) = zhsc(i,j)              - z_uv(i,j,nbdsc(i,j))
        end if
      else
        !----------------------------------------------------------
        ! Indicates a layer of zero depth
        !----------------------------------------------------------
        if ( ( sc_diag_opt==sc_diag_orig .or. sc_diag_opt==sc_diag_cu_relax )  &
             .and. ntdsc(i,j) == ntpar(i,j) ) then
          !----------------------------------------------------------
          ! Indicates a Sc layer at the top of Cu: force mixing
          ! over single layer.
          !----------------------------------------------------------
          dscdepth(i,j) = dzl(i,j,ntdsc(i,j))
        else
          dsc(i,j)=.false.
          ntdsc(i,j)=0
          zhsc(i,j)=zero
          df_dsct_over_cp(i,j) = zero
          k_cloud_dsct(i,j) = 0
          df_inv_dsc(i,j)   = zero
          dscdepth(i,j) = zero
        end if
      end if
    else  ! ntdsc == 0, just to make sure
      dscdepth(i,j)=zero
      dsc(i,j)=.false.
      zhsc(i,j)=zero
      df_dsct_over_cp(i,j) = zero
      k_cloud_dsct(i,j) = 0
      df_inv_dsc(i,j)   = zero
    end if
  end do
end do
!$OMP end do NOWAIT

!----------------------------------------------------------------------
!2.6 If decoupled cloud-layer found test to see if it is, in fact,
!  only weakly decoupled from the surface mixed-layer:
!  if SVL difference between NTML and NTDSC is less than svl_coup (in K)
!  then assume there is still some coupling.  This will mean that
!  the surface-driven entrainment term will be applied at ZHSC, no
!  subgrid inversion or entrainment will be calculated for ZH and
!  ZHSC will be the length scale used in the entrainment inputs.
!  Note that for CUMULUS "surface-driven entrainment" will be done
!  by the convection scheme.
!----------------------------------------------------------------------

if ( entr_smooth_dec == on .or. entr_smooth_dec == entr_taper_zh ) then
  !-------------------------------------------------------------
  ! entr_smooth_dec:
  ! OFF - original method
  ! ON  - taper off surface terms to zero for svl_diff between
  !       svl_coup and svl_coup_max; also ignore cumulus diags
  !-------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      coupled(i,j)       = .false.
      svl_diff_frac(i,j) = zero   ! Fully coupled by default
      if ( dsc(i,j) ) then
        !-------------------------------------------------------------
        ! Calculate cloud to surface mixed layer SVL difference
        ! - avoid ntdsc as can be within base of inversion
        !-------------------------------------------------------------
        svl_diff = zero
        if ( ntdsc(i,j) >= 2 )                                                 &
                 svl_diff = svl(i,j,ntdsc(i,j)-1) - svl(i,j,ntml(i,j))
        if ( svl_diff  < svl_coup_max ) then
          coupled(i,j) = .true.
          svl_diff_frac(i,j) = one - max( zero,                                &
                               (svl_diff-svl_coup)/(svl_coup_max-svl_coup) )
                              ! to give 1 for svl_diff<svl_coup and
                              ! decrease linearly to 0 at svl_coup_max
          if ( entr_smooth_dec == entr_taper_zh ) then
            ! Adjust DSC depth smoothly from existing value to full height
            ! of the Sc-top as a function of coupling strength.
            dscdepth(i,j) =      svl_diff_frac(i,j)  * zhsc(i,j)               &
                          + (one-svl_diff_frac(i,j)) * dscdepth(i,j)
          else
            ! Not using smooth height option; jump fully to Sc-top if coupled
            dscdepth(i,j) = zhsc(i,j)
          end if
        end if
      end if  ! dsc test
    end do
  end do
!$OMP end do

else  ! entr_smooth_dec test

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      coupled(i,j)       = .false.
      svl_diff_frac(i,j) = zero   ! Fully coupled by default
      if ( dsc(i,j) .and. .not. cumulus(i,j) ) then
        !-----------------------------------------------------------
        ! Note this if test structure is required because if DSC is
        ! false then NTDSC = 0 and cannot be used to index SVL.
        !-----------------------------------------------------------
        if ( svl(i,j,ntdsc(i,j)) - svl(i,j,ntml(i,j)) < svl_coup )             &
               coupled(i,j) = .true.
      end if
    end do
  end do
!$OMP end do

end if  ! entr_smooth_dec test

! Store current values of ntml and ntdsc for water tracer use
if (l_wtrac) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ntml_start(i,j)  = ntml(i,j)
      ntdsc_start(i,j) = ntdsc(i,j)
    end do
  end do
!$OMP end do
end if

!-----------------------------------------------------------------------
! 3. Diagnose a discontinuous inversion structure:
!    - to this point in the code, ZH and ZHSC mark the half-level at
!      the base of the inversion
!    - now they will be interpolated into the level above assuming
!      SVL(NTML+1) is a volume average over a subgrid discontinuous
!      inversion structure
!    - discontinuous jumps of SL and QW (and thence buoyancy) can be
!      calculated and used to determine the entrainment rate
!    - parametrized grid-level fluxes at NTML,NTDSC can then be made
!      consistent with this assumed inversion structure
!-----------------------------------------------------------------------
       ! If any `problems' are encountered with this interpolation of ZH
       ! (such as ZH diagnosed at or below Z_UV(NTML+1)), then NTML
       ! is lowered a level and ZH is set fractionally below what has
       ! become Z_UV(NTML+2).  This distance is such that for a net
       ! dZH/dt of 1.E-4 m/s, ZH will be diagnosed as spending at least
       ! half the timestep in level NTML+2, leaving the growth only
       ! marginally affected.  Conversely, it allows a subsiding
       ! inversion to fall more readily.

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    if (l_wtrac) then
      dqw_sml_meth(i,j) = 0
      qw_lapse_zero_sml(i,j) = .false.
    end if

    sml_disc_inv(i,j) = 0  ! initialise flags to indicate whether a
    dsc_disc_inv(i,j) = 0  ! discontinuous inversion is diagnosed
    res_inv(i,j)      = 0  ! Flag for whether inversion is resolved
    if ( bl_res_inv /= off .and. .not. cumulus(i,j) .and.                      &
         dzh(i,j) > one .and. ntml(i,j)+2 <= bl_levels ) then
      if (zh(i,j)+dzh(i,j) > z_uv(i,j,ntml(i,j)+2) ) res_inv(i,j) = 1
    end if

    !..First interpolate to find ZH

    k = ntml(i,j)
    !..by default, keep ZH at the half-level where it was diagnosed
    !..initially and use grid-level jumps

    dsl_sml(i,j) = sl(i,j,k+1) - sl(i,j,k)
    dqw_sml(i,j) = qw(i,j,k+1) - qw(i,j,k)

    if ( .not. cumulus(i,j) .and. .not. coupled(i,j) .and.                     &
         res_inv(i,j) == 0 .and. k > 1 .and. k <= bl_levels-2 ) then

        !  Require SVL and SL to be monotonically increasing
        !  and QW to be simply monotonic
      monotonic_inv = ( svl(i,j,k+2) > svl(i,j,k+1) .and.                      &
                        svl(i,j,k+1) > svl(i,j,k) )                            &
                .and. ( sl(i,j,k+2) > sl(i,j,k+1) .and.                        &
                        sl(i,j,k+1) > sl(i,j,k) )                              &
                .and. ( ( qw(i,j,k+2) > qw(i,j,k+1) .and.                      &
                          qw(i,j,k+1) > qw(i,j,k) )                            &
                    .or. ( qw(i,j,k+2) < qw(i,j,k+1) .and.                     &
                          qw(i,j,k+1) < qw(i,j,k) ) )

      if ( monotonic_inv ) then

        if ( k  <=  bl_levels-3 ) then
          ! need to test for K+1 to K+2 gradient in case profile is
          ! concave (would mess up the inversion diagnosis so best
          ! just to ignore lapse)
          svl_lapse = max(zero,                                                &
                ( svl(i,j,k+3) - svl(i,j,k+2) ) * rdz(i,j,k+3)  )
          if ( svl_lapse  >                                                    &
                ( svl(i,j,k+2) - svl(i,j,k+1) ) * rdz(i,j,k+2) )               &
                svl_lapse = zero
        else
          svl_lapse = zero
        end if
        if ( k  >=  k_plume(i,j)+2 ) then
              ! Use mean mixed layer gradient (if resolved) to allow
              ! for stablisation by gradient-adjustment
              ! Ignore level K in case inversion is dropping
          svl_lapse_base = ( svl(i,j,k-1)-svl(i,j,k_plume(i,j)) )/             &
                        (z_tq(i,j,k-1)-z_tq(i,j,k_plume(i,j)))
          svl_lapse_base = max( zero, svl_lapse_base )
        else
          svl_lapse_base = zero
        end if

        quad_a  = one_half*( svl_lapse - svl_lapse_base )
        quad_bm = svl(i,j,k+2) - svl(i,j,k)                                    &
            - svl_lapse * ( z_tq(i,j,k+2)-z_uv(i,j,k+2) )                      &
            - svl_lapse_base * ( z_uv(i,j,k+1)-z_tq(i,j,k) +                   &
                                                    dzl(i,j,k+1) )
        quad_c  = dzl(i,j,k+1)*( svl(i,j,k+1) - svl(i,j,k) -                   &
            svl_lapse_base * (                                                 &
              z_uv(i,j,k+1)-z_tq(i,j,k) + one_half*dzl(i,j,k+1) ) )

        if ( quad_bm  >   zero ) then
          if ( quad_c  <=  zero) then
                ! SVL extrapolated from K to K+1 is greater than
                ! the level K+1 value - inversion needs to rise so
                ! place it as high as possible
            dz_disc = dz_disc_min
          else if ( quad_bm*quad_bm  >=  4.0_r_bl*quad_a*quad_c ) then
                ! solve equation for DZ_DISC...
            if ( abs(quad_a)  >=  rbl_eps ) then
                  !   ...quadratic if QUAD_A /= 0
              dz_disc = ( quad_bm - sqrt( quad_bm*quad_bm                      &
                                       - 4.0_r_bl*quad_a*quad_c )              &
                              ) / (2.0_r_bl*quad_a)
            else
                  !   ...linear if QUAD_A == 0
              dz_disc = quad_c / quad_bm
            end if
          else
            dz_disc = 99999.9_r_bl  ! large dummy value
          end if

          if ( dz_disc  >   0.9_r_bl * dzl(i,j,k+1) ) then
              ! ZH diagnosed very close to or below Z_UV(K+1):
            if ( svl(i,j,k)-svl(i,j,k-1)  >   zero) then
                  ! top of ML stably stratified so lower NTML but
                  ! set ZH only fractionally (DZ_DISC_MIN)
                  ! below the top of the inversion level.
              ntml(i,j) = ntml(i,j) - 1
              k=ntml(i,j)
              dz_disc = dz_disc_min
            else
                  ! top of ML well-mixed so don't lower the inversion
                  ! level but set ZH just (DZ_DISC_MIN) above the
                  ! half-level to allow the inversion to subside if
                  ! necessary.
              dz_disc = dzl(i,j,k+1) - dz_disc_min
            end if
          end if

        else
          !.. ignoring lapse rates
          dsvl_top = svl(i,j,k+2) - svl(i,j,k)
          dz_disc = dzl(i,j,k+1) *                                             &
                          (svl(i,j,k+1)-svl(i,j,k)) / dsvl_top
        end if

        zh(i,j) = z_uv(i,j,k+2) - dz_disc
        sml_disc_inv(i,j) = 1 ! set flag to indicate disc inv found

        !-----------------------------------------------------------
        !..Calculate SML inversion discontinuous jumps of SL and QW
        !-----------------------------------------------------------
                    ! Allow for lapse rate above inversion, if known
        dz_disc = z_tq(i,j,k+2) - zh(i,j)
        sl_lapse = zero
        qw_lapse = zero
        if ( k  <=  bl_levels-3 ) then
          sl_lapse = max( zero,                                                &
             ( sl(i,j,k+3) - sl(i,j,k+2) )*rdz(i,j,k+3) )
          qw_lapse = min( zero,                                                &
             ( qw(i,j,k+3) - qw(i,j,k+2) )*rdz(i,j,k+3) )

          if (l_wtrac) then      ! Store method
            if (qw_lapse >= 0.0) then
              qw_lapse_zero_sml(i,j) = .true.
            end if
          end if

        end if

          !-----------------
          ! First SL jump
          !-----------------
          ! Only reduce 2 level jump by at most half
        dsl_sml(i,j) = sl(i,j,k+2) - sl(i,j,k)
        dsl_sml(i,j) = dsl_sml(i,j) -                                          &
               min( one_half*dsl_sml(i,j), sl_lapse*dz_disc )
          !-----------------
          ! Next QW jump
          !-----------------
        if ( qw(i,j,k+2) < qw(i,j,k+1) .and.                                   &
             qw(i,j,k+1) < qw(i,j,k) ) then
            ! QW monotonically decreasing across inversion
            ! Only allow for QW lapse rate if both it and the
            ! 2 grid-level jump are negative (expected sign)
          dqw_sml(i,j) = qw(i,j,k+2) - qw(i,j,k)
          if (l_wtrac) dqw_sml_meth(i,j) = 1

          if ( dqw_sml(i,j) < zero ) then
            if (l_wtrac) then           ! Store method
              if (one_half*dqw_sml(i,j) > qw_lapse*dz_disc) then
                dqw_sml_meth(i,j) = 2
              else
                dqw_sml_meth(i,j) = 3
              end if
            end if
            dqw_sml(i,j) = dqw_sml(i,j) -                                      &
              max( one_half*dqw_sml(i,j), qw_lapse*dz_disc )
          end if
        else if ( qw(i,j,k+2) > qw(i,j,k+1) .and.                              &
                  qw(i,j,k+1) > qw(i,j,k) ) then
            ! QW monotonically increasing across inversion
            ! Suggests something unusual is going so not clear how
            ! to proceed, so currently leaving DQW as 2 level jump
          dqw_sml(i,j) = qw(i,j,k+2) - qw(i,j,k)
          if (l_wtrac) dqw_sml_meth(i,j) = 1
        end if

      end if  ! Monotonic inversion
    end if ! not cumulus and not at top of bl_levels

    !-------------------------------------------------------------------
    !..Second interpolate to find ZHSC
    !-------------------------------------------------------------------
    if (l_wtrac) then
      dqw_dsc_meth(i,j) = 0
      qw_lapse_zero_dsc(i,j) = .false.
    end if

    if ( dsc(i,j) ) then
      k = ntdsc(i,j)
      !..by default, keep ZHSC at the half-level where it was diagnosed
      !..initially and use grid-level jumps
      dsl_dsc(i,j) = sl(i,j,k+1) - sl(i,j,k)
      dqw_dsc(i,j) = qw(i,j,k+1) - qw(i,j,k)

      if ( k  <=  bl_levels-2 ) then

        !  Require SVL and SL to be monotonically increasing
        !  and QW to be simply monotonic
        monotonic_inv = ( svl(i,j,k+2) > svl(i,j,k+1) .and.                    &
                          svl(i,j,k+1) > svl(i,j,k) )                          &
                  .and. ( sl(i,j,k+2) > sl(i,j,k+1) .and.                      &
                          sl(i,j,k+1) > sl(i,j,k) )                            &
                  .and. ( ( qw(i,j,k+2) > qw(i,j,k+1) .and.                    &
                            qw(i,j,k+1) > qw(i,j,k) )                          &
                      .or. ( qw(i,j,k+2) < qw(i,j,k+1) .and.                   &
                            qw(i,j,k+1) < qw(i,j,k) ) )

        if ( monotonic_inv ) then

          if ( sc_diag_opt == sc_diag_cu_rh_max .or.                           &
               sc_diag_opt == sc_diag_all_rh_max ) then
            ! The initial zhsc can be between model-levels rather than exactly
            ! on a rho-level.
            ! Store height of base of DSC layer
            z_cbase = zhsc(i,j) - dscdepth(i,j)
          end if

          if ( k  <=  bl_levels-3 ) then
            ! need to test for K+1 to K+2 gradient in case profile is
            ! concave (would mess up the inversion diagnosis so best
            ! just to ignore)
            svl_lapse = max(zero,                                              &
                  ( svl(i,j,k+3) - svl(i,j,k+2) )*rdz(i,j,k+3) )
            if ( svl_lapse  >                                                  &
                  ( svl(i,j,k+2) - svl(i,j,k+1) )*rdz(i,j,k+2) )               &
                  svl_lapse = zero
          else
            svl_lapse = zero
          end if
          if ( k  >=  nbdsc(i,j)+2 ) then
              ! Use mean mixed layer gradient (if resolved) to allow
              ! for stablisation by gradient-adjustment
              ! Ignore level K in case inversion is dropping
            svl_lapse_base = ( svl(i,j,k-1)-svl(i,j,nbdsc(i,j)) )/             &
                          (z_tq(i,j,k-1)-z_tq(i,j,nbdsc(i,j)))
            svl_lapse_base = max( zero, svl_lapse_base )
          else
            svl_lapse_base = zero
          end if

          quad_a  = one_half*( svl_lapse - svl_lapse_base )
          quad_bm = svl(i,j,k+2) - svl(i,j,k)                                  &
               - svl_lapse * ( z_tq(i,j,k+2)-z_uv(i,j,k+2) )                   &
               - svl_lapse_base * ( z_uv(i,j,k+1)-z_tq(i,j,k) +                &
                                                      dzl(i,j,k+1) )
          quad_c  = dzl(i,j,k+1)*( svl(i,j,k+1) - svl(i,j,k) -                 &
               svl_lapse_base * (                                              &
                z_uv(i,j,k+1)-z_tq(i,j,k) + one_half*dzl(i,j,k+1) ) )

          if ( quad_bm  >   zero ) then
            if ( quad_c  <=  zero) then
                ! SVL extrapolated from K to K+1 is greater than
                ! the level K+1 value - inversion needs to rise
              dz_disc = dz_disc_min
            else if ( quad_bm*quad_bm  >=  4.0_r_bl*quad_a*quad_c ) then
                ! solve equation for DZ_DISC...
              if ( abs(quad_a) >= rbl_eps ) then
                  !   ...quadratic if QUAD_A ne 0
                dz_disc = ( quad_bm - sqrt( quad_bm*quad_bm                    &
                                         - 4.0_r_bl*quad_a*quad_c )            &
                                ) / (2.0_r_bl*quad_a)
              else
                  !   ...linear if QUAD_A = 0
                dz_disc = quad_c / quad_bm
              end if
            else
              dz_disc = 99999.9_r_bl  ! large dummy value
            end if

            if ( dz_disc  >   0.9_r_bl * dzl(i,j,k+1) ) then
              if ( ntdsc(i,j) == 2 ) then
                dz_disc = dzl(i,j,k+1)
              else
                ! ZHSC diagnosed very close to or below Z_UV(K+1):
                if ( svl(i,j,k)-svl(i,j,k-1)  >   zero) then
                  ! top of ML stably stratified so lower NTDSC but
                  ! set ZHSC only fractionally (DZ_DISC_MIN)
                  ! below the top of the inversion level.
                  ntdsc(i,j) = ntdsc(i,j) - 1
                  k=ntdsc(i,j)
                  dz_disc = dz_disc_min
                  dscdepth(i,j) = dscdepth(i,j) - dzl(i,j,k+1)
                  ! Note that all but DZ_DISC_MIN of this layer will
                  ! be added back on to DSCDEPTH a few lines below
                else
                  ! top of ML well-mixed so don't lower the inversion
                  ! level but set ZHSC just (DZ_DISC_MIN) above the
                  ! half-level to allow the inversion to subside if
                  ! necessary.
                  dz_disc = dzl(i,j,k+1) - dz_disc_min
                end if
              end if
            end if

          else  ! QUAD_BM le 0
            !.. ignoring lapse rates
            dsvl_top = svl(i,j,k+2) - svl(i,j,k)
            dz_disc = dzl(i,j,k+1) *                                           &
                            (svl(i,j,k+1)-svl(i,j,k)) / dsvl_top
          end if

          zhsc(i,j) = z_uv(i,j,k+2) - dz_disc

          if ( sc_diag_opt == sc_diag_cu_rh_max .or.                           &
               sc_diag_opt == sc_diag_all_rh_max ) then
            ! The initial zhsc can be between model-levels rather than exactly
            ! on a rho-level.  Assuming height of base of DSC layer stays the
            ! same, set new depth based on updated DSC top height:
            dscdepth(i,j) = zhsc(i,j) - z_cbase
          else
            ! The initial zhsc is always at rho-level ntdsc+1;
            ! increment the DSC depth consistent with this:
            dscdepth(i,j) = dscdepth(i,j) + zhsc(i,j) - z_uv(i,j,k+1)
          end if

          dsc_disc_inv(i,j) = 1  ! set flag to indicate disc inv found

          !-----------------------------------------------------------
          !..Calculate DSC inversion discontinuous jumps of SL and QW
          !-----------------------------------------------------------
                      ! Allow for lapse rate above inversion, if known
          dz_disc = z_tq(i,j,k+2) - zhsc(i,j)
          sl_lapse = zero
          qw_lapse = zero
          if ( k  <=  bl_levels-3 ) then
            sl_lapse = max( zero,                                              &
               ( sl(i,j,k+3) - sl(i,j,k+2) )*rdz(i,j,k+3) )
            qw_lapse = min( zero,                                              &
               ( qw(i,j,k+3) - qw(i,j,k+2) )*rdz(i,j,k+3) )

            if (l_wtrac) then      ! Store method
              if (qw_lapse >= 0.0) then
                qw_lapse_zero_dsc(i,j) = .true.
              end if
            end if

          end if

          !-----------------
          ! First SL jump
          !-----------------
          ! Only reduce 2 level jump by at most half
          dsl_dsc(i,j) = sl(i,j,k+2) - sl(i,j,k)
          dsl_dsc(i,j) = dsl_dsc(i,j) -                                        &
               min( one_half*dsl_dsc(i,j), sl_lapse*dz_disc )
          !-----------------
          ! Next QW jump
          !-----------------
          if ( qw(i,j,k+2) < qw(i,j,k+1) .and.                                 &
               qw(i,j,k+1) < qw(i,j,k) ) then
            ! QW monotonically decreasing across inversion
            ! Only allow for QW lapse rate if both it and the
            ! 2 grid-level jump are negative (expected sign)
            dqw_dsc(i,j) = qw(i,j,k+2) - qw(i,j,k)
            if (l_wtrac) dqw_dsc_meth(i,j) = 1

            if ( dqw_dsc(i,j) < zero ) then
              if (l_wtrac) then       ! Store method
                if ( one_half*dqw_dsc(i,j) > qw_lapse*dz_disc ) then
                  dqw_dsc_meth(i,j) = 2
                else
                  dqw_dsc_meth(i,j) = 3
                end if
              end if
              dqw_dsc(i,j) = dqw_dsc(i,j) -                                    &
                max( one_half*dqw_dsc(i,j), qw_lapse*dz_disc )
            end if
          else if ( qw(i,j,k+2) > qw(i,j,k+1) .and.                            &
                    qw(i,j,k+1) > qw(i,j,k) ) then
            ! QW monotonically increasing across inversion
            ! Suggests something unusual is going so not clear how
            ! to proceed, so currently leaving DQW as 2 level jump
            dqw_dsc(i,j) = qw(i,j,k+2) - qw(i,j,k)
            if (l_wtrac) dqw_dsc_meth(i,j) = 1
          end if

        end if  ! monotonic inversion
      end if ! test on K lt BL_LEVELS-2
    end if ! test on DSC
  end do
end do
!$OMP end do

!$OMP end PARALLEL

! Repeat necessary parts of last loop for water tracers
if (l_wtrac) then

  ! SML
  ! Create dummy field to provide SML equivalent of 'dsc' array which
  ! is needed in calc_dqw_inv_wtrac
  allocate(sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  sml(:,:) = .true.

  call calc_dqw_inv_wtrac(bl_levels, ntml_start, ntml, dqw_sml_meth,           &
                          sml_disc_inv, rdz, z_tq, zh, sml,                    &
                          qw_lapse_zero_sml, wtrac_bl, dqw_sml_wtrac)
  deallocate(sml)

  !DSC
  call calc_dqw_inv_wtrac(bl_levels, ntdsc_start, ntdsc, dqw_dsc_meth,         &
                          dsc_disc_inv, rdz, z_tq, zhsc, dsc,                  &
                          qw_lapse_zero_dsc, wtrac_bl, dqw_dsc_wtrac)
end if  ! l_wtrac

! Set flags used in the code to find max cloud-fraction at mixed-layer top...
if ( l_use_sml_dsc_fixes ) then
  ! When a cloudy inversion has almost risen to the next level up
  ! at z_rho(ntml+2), quite often the max cloud-fraction in the column can
  ! occur at theta-level ntml+1, since that level is almost fully within
  ! the mixed-layer at that point.
  ! But the original code didn't check for cloud at ntml+1 when setting
  ! the mixed-layer-top cloud-fraction, which is used to estimate
  ! various properties of the cloud-top-driven mixing.
  ! If fixes are switched on, check for cloud at ntml+1 (and ntdsc+1).
  l_check_ntp1 = .true.
  p1 = 1
else
  ! Otherwise, the default code only checks ntml and ntml-1
  l_check_ntp1 = .false.
  p1 = 0
end if

!-----------------------------------------------------------------------
! 4.  Calculate the within-layer vertical gradients of cloud liquid
!      and frozen water
! 4.1 Calculate gradient adjustment terms
!-----------------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP  private (i, j, jj, i_wt, k, kl, km, kp, kp2, kmax, wstar3, c_ws, w_m,   &
!$OMP  pr_neut, w_h, k_cff, virt_factor, z_cbase , zdsc_cbase, dsl_ga,         &
!$OMP  dqw_ga, cfl_ml, cff_ml, dqw, dsl, dqcl, dqcf, db_disc, cu_depth_fac,    &
!$OMP  k_rad_lim, z_rad_lim ,dfsw_inv, dflw_inv, dfsw_top, dsldz, cf_for_wb,   &
!$OMP  grad_t_adj_inv_rdz, grad_q_adj_inv_rdz)

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    grad_t_adj(i,j) = zero
    grad_q_adj(i,j) = zero
    if ( unstable(i,j) ) then
        ! Here this is an estimate of the gradient adjustment applied
        ! the previous timestep (assumes T1_SD has not changed much,
        ! which in turn assumes the surface fluxes have not)
      if (flux_grad  ==  Locketal2000) then
        grad_t_adj(i,j) = min( max_t_grad ,                                    &
                         a_grad_adj * t1_sd(i,j) / zh_prev(i,j) )
        grad_q_adj(i,j) = zero
      else if (flux_grad  ==  HoltBov1993) then
          ! Use constants from Holtslag and Boville (1993)
          ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
          ! Neut limit GAMMA_TH = 7.2*wstar*FTL1/(ustar^2*zh)
        wstar3 = fb_surf(i,j) * zh_prev(i,j)
        c_ws = 0.6_r_bl
        w_m =( v_s(i,j)**3 + c_ws*wstar3 )**one_third

        grad_t_adj(i,j) = a_ga_hb93*(wstar3**one_third)*ftl(i,j,1)             &
                          / ( rhostar_gb(i,j)*w_m*w_m*zh_prev(i,j) )
        ! GRAD_Q_ADJ(I,j) = A_GA_HB93*(WSTAR3**one_third)*FQW(I,j,1)
        !                  / ( RHOSTAR_GB(I,j)*W_M*W_M*ZH_PREV(I,j) )
        grad_q_adj(i,j) = zero
      else if (flux_grad  ==  LockWhelan2006) then
          ! Use constants from LockWhelan2006
          ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
          ! Neut limit GAMMA_TH = 7.5*FTL1/(ustar*zh)
        wstar3  = fb_surf(i,j) * zh_prev(i,j)
        c_ws    = 0.42_r_bl   !  = 0.75^3
        pr_neut = 0.75_r_bl
        w_h = ( ( v_s(i,j)**3+c_ws*wstar3 )**one_third )/ pr_neut

        grad_t_adj(i,j) = a_ga_lw06 * ftl(i,j,1)                               &
                           / ( rhostar_gb(i,j)*w_h*zh_prev(i,j) )
        grad_q_adj(i,j) = a_ga_lw06 * fqw(i,j,1)                               &
                           / ( rhostar_gb(i,j)*w_h*zh_prev(i,j) )
      end if
    end if  ! test on UNSTABLE

  end do
end do
!$OMP end do

! Water tracers assume flux_grad  ==  Locketal2000
if (l_wtrac) then
  do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        wtrac_bl(i_wt)%grad_q_adj(i,j) = zero
      end do
    end do
!$OMP end do
  end do
end if

!-----------------------------------------------------------------------
! 4.2 Calculate the within-layer vertical gradients of cloud liquid
!     and frozen water for the current layer
!-----------------------------------------------------------------------

! In the following loop the computation is divided into three sub-loops
! on the 'i' index. This is to help the compiler vectorising the
! middle loop, and thus producing higher performance, since it does not
! contain any conditional statement. The compiler directive NOFUSION
! is used to stop the compiler fusing the sub-loops, while the
! directive VECTOR ALWAYS indicates that the sub-loop can be
! vectorised. These directives work with both Cray and Intel
! compilers.

!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels
  do j = pdims%j_start, pdims%j_end

    do i = pdims%i_start, pdims%i_end
      if (k  <=  ntml_prev(i,j)) then
        dsldz(i) = -grcp + grad_t_adj(i,j)
      else
        dsldz(i) = -grcp
      end if
    end do

!DIR$ NOFUSION
!DIR$ VECTOR ALWAYS
    do i = pdims%i_start, pdims%i_end
      virt_factor = one + c_virtual*q(i,j,k) - qcl(i,j,k) -                    &
                          qcf(i,j,k)

      dqcldz(i,j,k) = -( dsldz(i)*dqsdt(i,j,k)                                 &
                     + g*qs(i,j,k)/(r*t(i,j,k)*virt_factor) )                  &
                      / ( one + lcrcp*dqsdt(i,j,k) )
      dqcfdz(i,j,k) = -( dsldz(i)*dqsdt(i,j,k)                                 &
                     + g*qs(i,j,k)/(r*t(i,j,k)*virt_factor) ) * fgf            &
                      / ( one + lsrcp*dqsdt(i,j,k) )
    end do

        ! limit calculation to greater than a small cloud fraction
!DIR$ NOFUSION
    do i = pdims%i_start, pdims%i_end
      if ( qcl(i,j,k) + qcf(i,j,k)  >   zero                                   &
           .and. cf(i,j,k)  >   1.0e-3_r_bl ) then
        cfl(i,j,k) = cf(i,j,k) * qcl(i,j,k) /                                  &
                     ( qcl(i,j,k) + qcf(i,j,k) )
        cff(i,j,k) = cf(i,j,k) * qcf(i,j,k) /                                  &
                     ( qcl(i,j,k) + qcf(i,j,k) )
      else
        cfl(i,j,k) = zero
        cff(i,j,k) = zero
      end if

    end do
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
! 5.  Calculate uniform mixed-layer cloud fraction and thence
!        estimate Sc layer cloud depth (not cloud fraction weighted).
!        (If DSC=.false. then NTDSC=0 and ZC_DSC remains equal to 0.)
!-----------------------------------------------------------------------
! First the SML
!---------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    cloud_base(i,j)= .false.
    zc(i,j)        = zero
    k_cbase(i,j)   = 0
    z_cf_base(i,j) = zh(i,j)
    z_ctop(i,j)    = zh(i,j)
      ! Use a single CF for whole mixed-layer (more realistic).
      ! Include NTML+1 if a subgrid inversion has been diagnosed
    if ( coupled(i,j) .or. cumulus(i,j) .or. ntml(i,j) == 1 ) then
      cf_sml(i,j)=zero
    else
      k = ntml(i,j)
      cf_sml(i,j) = max( cf(i,j,k), cf(i,j,k-1) )
      if ( l_check_ntp1 ) cf_sml(i,j) = max( cf_sml(i,j), cf(i,j,k+1) )
    end if
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! First find cloud-base as seen by the cloud scheme
! [K_LEVEL=first level below NTML with CF<SC_CFTOL and
!  height Z_CF_BASE=half-level above]
! to use as first guess or lower limit
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    k_level(i,j) = ntml(i,j)
    if ( cf_sml(i,j)  >   sc_cftol ) then
      if ( .not. l_check_ntp1 )  k_level(i,j) = ntml(i,j)-1
      do while ( cf(i,j,max(k_level(i,j),1))  >   sc_cftol                     &
                 .and. k_level(i,j)  >=  2 )
        k_level(i,j) = k_level(i,j) - 1
      end do
    end if
  end do
end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( cf_sml(i,j)  >   sc_cftol ) then
      if ( k_level(i,j)  ==  1 .and.                                           &
           cf(i,j,max(k_level(i,j),1))  >   sc_cftol) then
        z_cf_base(i,j) = zero
      else
        z_cf_base(i,j) = z_uv(i,j,k_level(i,j)+1)
      end if
      zc(i,j) = z_ctop(i,j) - z_cf_base(i,j)
    end if
  end do
end do
!$OMP end do NOWAIT

!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    do k = min(bl_levels, ntml(i,j)+p1), 1, -1
      if ( .not. cloud_base(i,j) .and.                                         &
           cf_sml(i,j)  >   sc_cftol ) then
             ! within cloudy boundary layer
        if ( k  ==  1) then
          cloud_base(i,j) = .true.
        else
          if ( cf(i,j,k-1)  <   cf(i,j,k) ) cloud_base(i,j) = .true.
        end if
        k_cbase(i,j) = k
      end if
    end do ! K
  end do ! I
end do ! j
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

        !--------------------------------------------------
        ! Use adiabatic qcl gradient to estimate cloud-base
        ! from in-cloud qcl in level K_CBASE
        ! If k_cbase = 0 then it hasn't been initialised
        !--------------------------------------------------

    if ( cloud_base(i,j) .and. k_cbase(i,j) /= 0 ) then
      z_cbase = z_tq(i,j,k_cbase(i,j)) -                                       &
                qcl(i,j,k_cbase(i,j)) /                                        &
                ( cf(i,j,k_cbase(i,j))*dqcldz(i,j,k_cbase(i,j)) )
      if ( dqcfdz(i,j,k_cbase(i,j))  >   zero ) then
        z_cbase = min( z_cbase, z_tq(i,j,k_cbase(i,j)) -                       &
               qcf(i,j,k_cbase(i,j)) /                                         &
               ( cf(i,j,k_cbase(i,j))*dqcfdz(i,j,k_cbase(i,j)) )               &
                     )
      else
            !---------------------------------------------------------
            ! No adiabatic QCF gradient so find lowest level, K_CFF,
            ! with CFF>SC_CFTOL and assume cloud-base within that leve
            !---------------------------------------------------------
        !Initialise K_CFF = lowest level with ice cloud
        k_cff = k_cbase(i,j)
        if (k_cff > 1) then
          do while ( cff(i,j,k_cff)  >   sc_cftol                              &
            .and. k_cff  >   1 )
            k_cff = k_cff - 1
          end do
        end if
        if ( cff(i,j,k_cff)  <=  sc_cftol .and.                                &
                      k_cff  <   k_cbase(i,j) )                                &
             k_cff = k_cff + 1
                 ! will want to raise K_CFF back up one level unless
                 ! level 1 is cloudy or no sig frozen cloud at all
        z_cbase = min( z_cbase, z_top(i,j,k_cff) -                             &
                  dzl(i,j,k_cff)                                               &
                * cff(i,j,k_cff)/cf(i,j,k_cff) )
      end if
          !------------------------------------------------------
          ! use cloud-base as seen by cloud scheme as lower limit
          ! and base of level NTML+1 as upper limit
          !------------------------------------------------------
      z_cbase = min( z_uv(i,j,ntml(i,j)+1),                                    &
                     max( z_cf_base(i,j), z_cbase) )

      zc(i,j) = z_ctop(i,j) - z_cbase
    end if

  end do
end do
!$OMP end do NOWAIT


!-----------------------------------------------------------------------
! Second DSC layer
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    cloud_base(i,j) = .false.
    zc_dsc(i,j) = zero
    k_cbase(i,j) = 0
    z_cf_base(i,j) = zhsc(i,j)
    z_ctop(i,j)    = zhsc(i,j)

    if ( dsc(i,j) ) then
      k = ntdsc(i,j)
      cf_dsc(i,j) = max( cf(i,j,k), cf(i,j,k-1) )
      if ( l_check_ntp1 ) cf_dsc(i,j) = max( cf_dsc(i,j), cf(i,j,k+1) )
    else
      cf_dsc(i,j) = zero
    end if
  end do
end do
!$OMP end do NOWAIT

!-------------------------------------------------------------
! Find cloud-base as seen by cloud scheme, Z_CF_BASE,
! to use as first guess or lower limit and find cloud top.
!-------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    k_level(i,j) = ntdsc(i,j)
    if ( cf_dsc(i,j)  >   sc_cftol ) then
        ! assume level NTDSC is cloudy so start from NTDSC-1
      if ( .not. l_check_ntp1 )  k_level(i,j) = max( 2, ntdsc(i,j) - 1 )
      do while ( cf(i,j,k_level(i,j))  >   sc_cftol                            &
               .and. k_level(i,j)  >=  2 )
        k_level(i,j) = k_level(i,j) - 1
      end do
    end if
  end do
end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( cf_dsc(i,j)  >   sc_cftol ) then
      if ( k_level(i,j)  ==  1 .and.                                           &
           cf(i,j,max(k_level(i,j),1))  >   sc_cftol) then
        z_cf_base(i,j) = zero
      else
        z_cf_base(i,j) = z_uv(i,j,k_level(i,j)+1)
      end if
      zc_dsc(i,j) = z_ctop(i,j) - z_cf_base(i,j)   ! first guess
    end if
  end do
end do
!$OMP end do NOWAIT

!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    do k = min(bl_levels,ntdsc(i,j)+p1), 1, -1
      if ( .not. cloud_base(i,j) .and.                                         &
           cf_dsc(i,j)  >   sc_cftol ) then
             ! within cloudy boundary layer
        if ( k  ==  1) then
          cloud_base(i,j) = .true.
        else
          if ( cf(i,j,k-1)  <   cf(i,j,k) ) cloud_base(i,j) = .true.
        end if
        k_cbase(i,j) = k
      end if
    end do ! K
  end do ! I
end do ! J
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

      !--------------------------------------------------
      ! use adiabatic qcl gradient to estimate cloud-base
      ! from in-cloud qcl in level K_CBASE
      !--------------------------------------------------
    if ( cloud_base(i,j) .and. k_cbase(i,j)  /= 0 ) then
      z_cbase = z_tq(i,j,k_cbase(i,j)) -                                       &
                qcl(i,j,k_cbase(i,j)) /                                        &
                ( cf(i,j,k_cbase(i,j))*dqcldz(i,j,k_cbase(i,j)) )
      if ( dqcfdz(i,j,k_cbase(i,j))  >   zero ) then
        z_cbase = min( z_cbase, z_tq(i,j,k_cbase(i,j)) -                       &
              qcf(i,j,k_cbase(i,j)) /                                          &
              ( cf(i,j,k_cbase(i,j))*dqcfdz(i,j,k_cbase(i,j)) )                &
                     )
      else
        ! Initialise K_CFF
        k_cff = k_cbase(i,j)
        if (k_cff > 1) then
          do while ( cff(i,j,k_cff)  >   sc_cftol                              &
            .and. k_cff  >   1)
            k_cff = k_cff - 1
          end do
        end if
        !----------------------------------------------------------
          ! No adiabatic QCF gradient so find lowest level, K_CFF,
          ! with CFF>SC_CFTOL and assume cloud-base within that level
          !----------------------------------------------------------
        if ( cff(i,j,k_cff)  <=  sc_cftol .and.                                &
                      k_cff  <   k_cbase(i,j) )                                &
             k_cff = k_cff + 1
               ! will want to raise K_CFF back up one level unless
               ! level 1 is cloudy or no sig frozen cloud at all
        z_cbase = min( z_cbase, z_top(i,j,k_cff) -                             &
                  dzl(i,j,k_cff)                                               &
                 * cff(i,j,k_cff)/cf(i,j,k_cff) )
      end if
        !------------------------------------------------------
        ! use cloud-base as seen by cloud scheme as lower limit
        ! and base of level NTDSC+1 as upper limit
        !------------------------------------------------------
      z_cbase = min( z_uv(i,j,ntdsc(i,j)+1),                                   &
                     max( z_cf_base(i,j) , z_cbase) )

      zc_dsc(i,j) = z_ctop(i,j) - z_cbase
    end if

    !-----------------------------------------------------------------
    !  Layer cloud depth cannot be > the layer depth itself.
    !-----------------------------------------------------------------

    zc_dsc(i,j) = min( zc_dsc(i,j), dscdepth(i,j) )

  end do !I
end do !J
!$OMP end do

!-----------------------------------------------------------------------
! 6. Calculate buoyancy flux factor used in the diagnosis of decoupling
!-----------------------------------------------------------------------

! In the following loop the computation is divided into two sub-loops
! on the 'i' index. This is to help the compiler vectorising the
! last loop, and thus producing higher performance, since it does not
! contain any conditional statement. The compiler directive NOFUSION
! is used to stop the compiler fusing the sub-loops, while the
! directive VECTOR ALWAYS indicates that the sub-loop can be
! vectorised. These directives work with both Cray and Intel
! compilers.

!$OMP do SCHEDULE(DYNAMIC)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end

    ! This is to help vectorization
    do i = pdims%i_start, pdims%i_end
      if ( k <= ntml_prev(i,j) .or. l_converge_ga ) then
        ! If using more accurate treatment of gradient adjustment in the
        ! buoyancy-flux integration, it needs to be set on all levels.
        grad_t_adj_inv_rdz(i) = grad_t_adj(i,j)/rdz(i,j,k)
        grad_q_adj_inv_rdz(i) = grad_q_adj(i,j)/rdz(i,j,k)
      else
        grad_t_adj_inv_rdz(i) = zero
        grad_q_adj_inv_rdz(i) = zero
      end if

      !----------------------------------------------------------
      ! CF_FOR_WB is uniform `bl' CF for use within cloud layers
      !----------------------------------------------------------
      cf_for_wb(i) = zero
      z_cbase = zh(i,j)-zc(i,j)
      zdsc_cbase = zhsc(i,j)-zc_dsc(i,j)
      if ( z_tq(i,j,k)  <=  zh(i,j) .and.                                      &
           z_tq(i,j,k)  >=  z_cbase) cf_for_wb(i) = cf_sml(i,j)
      if ( z_tq(i,j,k)  <=  zhsc(i,j) .and.                                    &
           z_tq(i,j,k)  >=  zdsc_cbase) cf_for_wb(i) = cf_dsc(i,j)
    end do

!DIR$ NOFUSION
!DIR$ VECTOR ALWAYS
    do i = pdims%i_start, pdims%i_end
      dqw = qw(i,j,k) - qw(i,j,k-1)
      dsl = sl(i,j,k) - sl(i,j,k-1)
      dsl_ga = dsl - grad_t_adj_inv_rdz(i)
      dqw_ga = dqw - grad_q_adj_inv_rdz(i)

        !----------------------------------------------------------
        ! WB = -K_SURF*(DB/DZ - gamma_buoy) - K_TOP*DB/DZ
        ! This is integrated in EXCF_NL, iterating the K profiles.
        ! Here the relevant integrated DB/DZ factors are calculated
        !----------------------------------------------------------
      db_ga_dry(i,j,k) = - g *                                                 &
               ( btm(i,j,k-1)*dsl_ga + bqm(i,j,k-1)*dqw_ga )
      db_noga_dry(i,j,k)  = - g *                                              &
               ( btm(i,j,k-1)*dsl + bqm(i,j,k-1)*dqw )
      db_ga_cld(i,j,k) = - g *                                                 &
               ( btm_cld(i,j,k-1)*dsl_ga + bqm_cld(i,j,k-1)*dqw_ga )
      db_noga_cld(i,j,k)  = - g *                                              &
               ( btm_cld(i,j,k-1)*dsl + bqm_cld(i,j,k-1)*dqw )
        !-------------------------------------------------------
        ! Weight cloud layer factors with cloud fraction
        !-------------------------------------------------------
      db_ga_cld(i,j,k) = db_ga_dry(i,j,k)*(one-cf_for_wb(i)) +                 &
                         db_ga_cld(i,j,k)*cf_for_wb(i)
      db_noga_cld(i,j,k)  = db_noga_dry(i,j,k)*(one-cf_for_wb(i)) +            &
                            db_noga_cld(i,j,k)*cf_for_wb(i)
    end do
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 7. Calculate inputs for the top of b.l. entrainment parametrization
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    zeta_r_dsc(i,j) = zero
    chi_s_dsct(i,j) = zero
    d_siems_dsc(i,j) = zero
    br_fback_dsc(i,j)= zero
    cld_factor_dsc(i,j) = zero
    bt_dsct(i,j) = zero
    btt_dsct(i,j) = zero
    db_dsct(i,j) = zero
    db_dsct_cld(i,j) = zero
    chi_s_top(i,j) = zero
    d_siems(i,j) = zero
    br_fback(i,j)= zero
    cld_factor(i,j) = zero
    bt_top(i,j) = zero
    btt_top(i,j) = zero
    btc_top(i,j) = zero
    db_top(i,j) = zero
    db_top_cld(i,j) = zero    ! default required if COUPLED
    z_cld(i,j) = zero
    z_cld_dsc(i,j) = zero
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
! 7.1 Calculate surface buoyancy flux
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

      ! use mixed-layer average of buoyancy parameters
    bflux_surf(i,j) = one_half * g * (                                         &
         (btm(i,j,1)+btm(i,j,ntml(i,j)))*ftl(i,j,1) +                          &
         (bqm(i,j,1)+bqm(i,j,ntml(i,j)))*fqw(i,j,1) )

    if ( bflux_surf(i,j)  >   zero ) then
      bflux_surf_sat(i,j) = one_half * g * (                                   &
         (btm_cld(i,j,1)+btm_cld(i,j,ntml(i,j)))*ftl(i,j,1) +                  &
         (bqm_cld(i,j,1)+bqm_cld(i,j,ntml(i,j)))*fqw(i,j,1) )
      if ( coupled(i,j) ) bflux_surf_sat(i,j) = one_half * g * (               &
         (btm_cld(i,j,1)+btm_cld(i,j,ntdsc(i,j)))*ftl(i,j,1) +                 &
         (bqm_cld(i,j,1)+bqm_cld(i,j,ntdsc(i,j)))*fqw(i,j,1) )
    else
      bflux_surf_sat(i,j) = zero
    end if

  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 7.2 Calculation of cloud fraction weighted thickness of
!     cloud in the SML and DSC layer (or to the surface if COUPLED)
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(DYNAMIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do k = 1, bl_levels
    do j = jj, min((jj+bl_segment_size)-1,pdims%j_end)
      do i = pdims%i_start, pdims%i_end
        if ( k  <=  ntml(i,j)+1 ) then
          z_cld(i,j) = z_cld(i,j) +                                            &
                     cf(i,j,k) * one_half * dzl(i,j,k) +                       &
                      min( cfl(i,j,k) * one_half * dzl(i,j,k) ,                &
                              qcl(i,j,k) / dqcldz(i,j,k) )
          if ( dqcfdz(i,j,k)  >   zero) then
            z_cld(i,j) = z_cld(i,j) +                                          &
                      min( cff(i,j,k) * one_half * dzl(i,j,k) ,                &
                              qcf(i,j,k) / dqcfdz(i,j,k) )
          else
            z_cld(i,j) = z_cld(i,j) + cff(i,j,k) * one_half * dzl(i,j,k)
          end if
        end if

        if ( dsc(i,j) .and. k <= ntdsc(i,j)+1 .and.                            &
             ( coupled(i,j) .or.                                               &
                   z_top(i,j,k) >= zhsc(i,j)-zc_dsc(i,j) ) ) then
          z_cld_dsc(i,j) = z_cld_dsc(i,j) +                                    &
                     cf(i,j,k) * one_half * dzl(i,j,k) +                       &
                      min( cfl(i,j,k) * one_half * dzl(i,j,k) ,                &
                              qcl(i,j,k) / dqcldz(i,j,k) )
          if ( dqcfdz(i,j,k)  >   zero) then
            z_cld_dsc(i,j) = z_cld_dsc(i,j) +                                  &
                      min( cff(i,j,k) * one_half * dzl(i,j,k) ,                &
                              qcf(i,j,k) / dqcfdz(i,j,k) )
          else
            z_cld_dsc(i,j) = z_cld_dsc(i,j) +                                  &
                                    cff(i,j,k) * one_half * dzl(i,j,k)
          end if
        end if
      end do
    end do
  end do
end do ! jj
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 7.3 Calculate the buoyancy jumps across the inversions
!----------------------------------------------------------------------
!..Where appropriate, replace grid-level jumps with jumps calculated
!..assuming a discontinuous subgrid inversion structure.
!----------------------------------------------------------------------

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
      !--------------------------
      ! First the SML
      !--------------------------
    qcl_ic_top(i,j) = zero
    qcf_ic_top(i,j) = zero
    cfl_ml = zero
    cff_ml = zero

    km = ntml(i,j)
    if ( sml_disc_inv(i,j) == 1 ) then
        !-----------------------------------------------------
        ! Extrapolate water contents from level with max CF,
        ! out of NTML and NTML-1 (ie near top of SML),
        ! to the top of the mixed layer
        !-----------------------------------------------------
      if (cf_sml(i,j)  >   zero) then
        kmax = km
        if (cf(i,j,km-1) > cf(i,j,km)) kmax = km-1
        if ( l_check_ntp1 .and. cf(i,j,km+1) > cf(i,j,kmax) ) kmax = km+1

        cfl_ml = cf_sml(i,j)*cfl(i,j,kmax)                                     &
                               /(cfl(i,j,kmax)+cff(i,j,kmax)+rbl_eps)
        cff_ml = cf_sml(i,j)*cff(i,j,kmax)                                     &
                               /(cfl(i,j,kmax)+cff(i,j,kmax)+rbl_eps)

        if (cfl_ml > 0.01_r_bl) qcl_ic_top(i,j) = qcl(i,j,kmax)/cfl_ml         &
                         + ( zh(i,j)-z_tq(i,j,km) )*dqcldz(i,j,km)
        if (cff_ml > 0.01_r_bl) qcf_ic_top(i,j) = qcf(i,j,kmax)/cff_ml         &
                         + ( zh(i,j)-z_tq(i,j,km) )*dqcfdz(i,j,km)
      end if

      dqw = dqw_sml(i,j)
      dsl = dsl_sml(i,j)
        ! ignore any cloud above the inversion
      dqcl = - cfl_ml*qcl_ic_top(i,j)
      dqcf = - cff_ml*qcf_ic_top(i,j)

      db_disc = g * ( btm(i,j,km)*dsl + bqm(i,j,km)*dqw +                      &
               (lcrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcl +                 &
               (lsrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcf   )

      if ( db_disc > 0.03_r_bl  ) then
          ! Diagnosed inversion statically stable and at least ~1K
        db_top(i,j) = db_disc
      else
          ! Diagnosed inversion statically UNstable
          ! Reset flag to use entrainment K (rather than fluxes)
        sml_disc_inv(i,j) = 0
        zh(i,j) = z_uv(i,j,ntml(i,j)+1)
      end if
    end if  ! disc inversion diagnosed

    if ( sml_disc_inv(i,j) == 0 ) then
         ! Calculate using simple grid-level differences
      kp = km+1
      dqw = qw(i,j,kp) - qw(i,j,km)
      dsl = sl(i,j,kp) - sl(i,j,km)
      qcl_ic_top(i,j) = qcl(i,j,km)
      qcf_ic_top(i,j) = qcf(i,j,km)
      dqcl = qcl(i,j,kp) - qcl_ic_top(i,j)
      dqcf = qcf(i,j,kp) - qcf_ic_top(i,j)
      db_top(i,j) = g * ( btm(i,j,km)*dsl + bqm(i,j,km)*dqw +                  &
                (lcrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcl +                &
                (lsrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcf )
    end if  ! no disc inversion diagnosed

    db_top_cld(i,j) = g * ( btm_cld(i,j,km)*dsl                                &
                          + bqm_cld(i,j,km)*dqw )
    denom = a_qsm(i,j,km)*dqw - a_dqsdtm(i,j,km)*dsl
    if (abs(denom) > rbl_eps) then
      chi_s_top(i,j) = -qcl_ic_top(i,j) / denom
      chi_s_top(i,j) = max( zero, min( chi_s_top(i,j), one) )
    end if

    if ( db_top(i,j)  <   0.003_r_bl ) then
        ! Diagnosed inversion statically unstable:
        ! ensure DB>0 so that entrainment is non-zero and
        ! instability can be removed.
      db_top(i,j) = 0.003_r_bl
      db_top_cld(i,j) = zero  ! set buoyancy reversal
      chi_s_top(i,j) = zero   ! term to zero
    end if

  end do
end do
!$OMP end do

        !--------------------------
        ! Then the DSC layer
        !--------------------------

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( dsc(i,j) ) then

      qcl_ic_top(i,j) = zero
      qcf_ic_top(i,j) = zero
      cfl_ml = zero
      cff_ml = zero

      km = ntdsc(i,j)
      if ( dsc_disc_inv(i,j) == 1 ) then
         !-----------------------------------------------------
         ! Extrapolate water contents from level with max CF,
         ! out of NTDSC and NTDSC-1 (ie near top of DSC),
         ! to the top of the mixed layer
         !-----------------------------------------------------
        if (cf_dsc(i,j) > zero) then
          kmax = km
          if (cf(i,j,km-1) > cf(i,j,km)) kmax = km-1
          if ( l_check_ntp1 .and. cf(i,j,km+1) > cf(i,j,kmax) ) kmax = km+1

          cfl_ml = cf_dsc(i,j)*cfl(i,j,kmax)                                   &
                                /(cfl(i,j,kmax)+cff(i,j,kmax)+rbl_eps)
          cff_ml = cf_dsc(i,j)*cff(i,j,kmax)                                   &
                                /(cfl(i,j,kmax)+cff(i,j,kmax)+rbl_eps)
          if (cfl_ml > 0.01_r_bl) qcl_ic_top(i,j) = qcl(i,j,kmax)/cfl_ml       &
                        + ( zhsc(i,j)-z_tq(i,j,km) )*dqcldz(i,j,km)
          if (cff_ml > 0.01_r_bl) qcf_ic_top(i,j) = qcf(i,j,kmax)/cff_ml       &
                        + ( zhsc(i,j)-z_tq(i,j,km) )*dqcfdz(i,j,km)

        end if

        dqw = dqw_dsc(i,j)
        dsl = dsl_dsc(i,j)
         ! ignore any cloud above the inversion
        dqcl = - cfl_ml*qcl_ic_top(i,j)
        dqcf = - cff_ml*qcf_ic_top(i,j)

        db_disc = g * ( btm(i,j,km)*dsl + bqm(i,j,km)*dqw +                    &
                (lcrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcl +                &
                (lsrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcf   )

        if ( db_disc > 0.03_r_bl  ) then
            ! Diagnosed inversion statically stable
          db_dsct(i,j) = db_disc
        else
            ! Diagnosed inversion statically UNstable
            ! Reset flag to use entrainment K (rather than fluxes)
          dsc_disc_inv(i,j) = 0
          zhsc(i,j) = z_uv(i,j,ntdsc(i,j)+1)
        end if
      end if  ! disc inversion diagnosed

      if ( dsc_disc_inv(i,j) == 0 ) then
         ! Calculate using simple grid-level differences
        kp = km+1
        dqw = qw(i,j,kp) - qw(i,j,km)
        dsl = sl(i,j,kp) - sl(i,j,km)
        qcl_ic_top(i,j) = qcl(i,j,km)
        qcf_ic_top(i,j) = qcf(i,j,km)
        dqcl = qcl(i,j,kp) - qcl_ic_top(i,j)
        dqcf = qcf(i,j,kp) - qcf_ic_top(i,j)
        db_dsct(i,j) = g * ( btm(i,j,km)*dsl + bqm(i,j,km)*dqw +               &
                  (lcrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcl +              &
                  (lsrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcf )
      end if  ! no disc inversion diagnosed

      db_dsct_cld(i,j) = g * ( btm_cld(i,j,km)*dsl                             &
                             + bqm_cld(i,j,km)*dqw )
      denom = a_qsm(i,j,km)*dqw - a_dqsdtm(i,j,km)*dsl
      if (abs(denom) > rbl_eps) then
        chi_s_dsct(i,j) = -qcl_ic_top(i,j) / denom
        chi_s_dsct(i,j) = max( zero, min( chi_s_dsct(i,j), one) )
      end if

      if ( db_dsct(i,j) < 0.003_r_bl ) then
         ! Diagnosed inversion statically unstable:
         ! ensure DB>0 so that entrainment is non-zero and
         ! instability can be removed.
        db_dsct(i,j) = 0.003_r_bl
        db_dsct_cld(i,j) = zero  ! set buoyancy reversal
        chi_s_dsct(i,j) = zero   ! term to zero
      end if
    end if  ! test on DSC

    !-----------------------------------------------------------------------
    ! 7.3 Calculation of other SML and DSC inputs to entr param.
    !     If COUPLED then SML are not used as no "entrainment" is then
    !     applied at ZH.
    !-----------------------------------------------------------------------

            !------------------------------------------------------
            ! Calculation of SML inputs.
            !------------------------------------------------------
    k = ntml(i,j)
    kp2=min(k+1+sml_disc_inv(i,j),bl_levels)
    cld_factor(i,j) = max( zero , cf_sml(i,j)-cf(i,j,kp2) )
    bt_top(i,j)  = g * btm(i,j,k)
    btt_top(i,j) = g * btm_cld(i,j,k)
    btc_top(i,j) = btt_top(i,j)
        !---------------------------------------------------
        ! Calculation of DSC inputs
        !---------------------------------------------------
    if (dsc(i,j)) then
      k = ntdsc(i,j)
      kp2=min(k+1+dsc_disc_inv(i,j),bl_levels)
      cld_factor_dsc(i,j) = max( zero , cf_dsc(i,j)-cf(i,j,kp2) )
      bt_dsct(i,j)  = g * btm(i,j,k)
      btt_dsct(i,j) = g * btm_cld(i,j,k)
    end if
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
! 7.4 Next those terms which depend on the presence of buoyancy reversal
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    z_cld(i,j) = min( z_cld(i,j), zh(i,j) )
    z_cld_dsc(i,j) = min( z_cld_dsc(i,j), zhsc(i,j) )
      !---------------------------------------------------------------
      ! First the surface mixed layer.
      !---------------------------------------------------------------
    if ( coupled(i,j) ) then
      zeta_s(i,j) = one - z_cld_dsc(i,j) / zhsc(i,j)
      zeta_r(i,j) = one - zc_dsc(i,j) / zhsc(i,j)
    else
      zeta_s(i,j) = one - z_cld(i,j) / zh(i,j)
      zeta_r(i,j) = one - zc(i,j) / zh(i,j)
    end if

    if (db_top_cld(i,j)  >=  zero) then
        !--------------------------------------------------
        ! i.e. no buoyancy reversal (or default if COUPLED)
        !--------------------------------------------------
      db_top_cld(i,j) = zero
      d_siems(i,j) = zero
      br_fback(i,j)= zero
    else
        !----------------------------
        ! if (DB_TOP_CLD(I,j)  <   0.0)
        ! i.e. buoyancy reversal
        !----------------------------
      db_top_cld(i,j) = -db_top_cld(i,j) * cld_factor(i,j)
      d_siems(i,j) = max( zero,                                                &
           chi_s_top(i,j) * db_top_cld(i,j) / (db_top(i,j)+rbl_eps) )
        ! Linear feedback dependence for D<0.1
      br_fback(i,j)= min( one, 10.0_r_bl*d_siems(i,j) )
      zeta_r(i,j)  = zeta_r(i,j) + (one-zeta_r(i,j))*br_fback(i,j)
    end if
      !---------------------------------------------------------------
      ! Now the decoupled Sc layer (DSC).
      !---------------------------------------------------------------
    if (dsc(i,j)) then
      if ( coupled(i,j) ) then
        zeta_r_dsc(i,j) = one - zc_dsc(i,j) / zhsc(i,j)
      else
        zeta_r_dsc(i,j) = one - zc_dsc(i,j) / dscdepth(i,j)
      end if

      if (db_dsct_cld(i,j)  >=  zero) then
          !----------------------------
          ! i.e. no buoyancy reversal
          !----------------------------
        db_dsct_cld(i,j) = zero
        d_siems_dsc(i,j) = zero
        br_fback_dsc(i,j)= zero
      else
          !----------------------------
          ! if (DB_DSCT_CLD(I,j)  <   0.0)
          ! i.e. buoyancy reversal
          !----------------------------
        db_dsct_cld(i,j) = -db_dsct_cld(i,j) * cld_factor_dsc(i,j)
        d_siems_dsc(i,j) = max( zero, chi_s_dsct(i,j)                          &
                        * db_dsct_cld(i,j) / (db_dsct(i,j)+rbl_eps) )
          ! Linear feedback dependence for D<0.1
        br_fback_dsc(i,j) = min( one, 10.0_r_bl*d_siems_dsc(i,j) )

        if ( entr_enhance_by_cu == Buoyrev_feedback                            &
             .and. cumulus(i,j)                                                &
             .and. d_siems_dsc(i,j) < 0.1_r_bl                                 &
             .and. d_siems_dsc(i,j) > rbl_eps ) then
            ! Assume mixing from cumulus can enhance the
            ! buoyancy reversal feedback in regime 0<D<0.1.
            ! Make enhancement dependent on Cu depth:
            !       Cu_depth_fac->0 below 400m, 1 above 1000m
          cu_depth_fac = one_half*( one+                                       &
                    tanh( ((zhpar(i,j)-zh(i,j))-700.0_r_bl)/100.0_r_bl) )
            ! BR_FBACK = unchanged for Cu<400m, ->1 for Cu>1000.
          br_fback_dsc(i,j) = cu_depth_fac +                                   &
                             (one-cu_depth_fac)*br_fback_dsc(i,j)
        end if

        zeta_r_dsc(i,j) = zeta_r_dsc(i,j) +                                    &
                         (one-zeta_r_dsc(i,j))*br_fback_dsc(i,j)

      end if
    end if
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
! 8. Calculate the radiative flux change across cloud top for mixed-
!    layer to ZH.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    k_cloud_top(i,j) = 0
    df_top_over_cp(i,j) = zero
    df_inv_sml(i,j) = zero
  end do
end do
!$OMP end do

if (l_new_kcloudtop) then
  !---------------------------------------------------------------------
  ! improved method of finding the k_cloud_dsct, the top of the mixed
  ! layer as seen by radiation
  !---------------------------------------------------------------------
  ! Find k_cloud_top as the equivalent to ntml (the top level of
  ! the SML) seen by radiation, which we take as the level two
  ! levels below the lowest level with free tropospheric radiative
  ! cooling.
  ! First find the level with maximum LW cooling, below z_rad_lim and
  ! in the upper half of the BL (ie, restrict search to `close' to ZH)
  !---------------------------------------------------------------------
!$OMP do SCHEDULE(DYNAMIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      k_rad_lim = ntml(i,j)+1
      z_rad_lim = max( z_tq(i,j,k_rad_lim)+0.1_r_bl, 1.2_r_bl*zh(i,j) )
      k = 1
      do while (z_tq(i,j,k) <  z_rad_lim .and. k < bl_levels)
        if (dflw_over_cp(i,j,k) > df_top_over_cp(i,j)                          &
            .and. z_tq(i,j,k) >  one_half*zh(i,j) ) then
          k_cloud_top(i,j) = k
          df_top_over_cp(i,j) = dflw_over_cp(i,j,k)
        end if
        k = k + 1
      end do ! k
      !-----------------------------------------------------------------
      ! If DF(K_CLOUD_TOP+1) is less than double DF(K_CLOUD_TOP+2) we
      ! assume DF(K_CLOUD_TOP+1) is actually typical of the free-trop and
      ! that the current K_CLOUD_TOP must be the inversion grid-level.
      ! Hence we lower K_CLOUD_TOP by one (it should mark the top of the
      ! mixed layer and cloud-top radiative cooling within the invesion
      ! grid-level will be included as DF_INV_SML below)
      !-----------------------------------------------------------------
      k = k_cloud_top(i,j)
      if ( k > 1 .and. k < bl_levels -1 ) then
        if (dflw_over_cp(i,j,k+1) < 1.5_r_bl*dflw_over_cp(i,j,k+2))            &
          k_cloud_top(i,j) = k-1
      end if
    end do ! i
    !-----------------------------------------------------------------
    !  Find bottom grid-level (K_LEVEL) for cloud-top radiative fux
    !  divergence: higher of base of LW radiatively cooled layer,
    !  0.5*ZH, since cooling must be in upper part of layer
    !-----------------------------------------------------------------
    do i = pdims%i_start, pdims%i_end
      k_level(i,j) = k_cloud_top(i,j)
      if ( k_cloud_top(i,j)  >   1 ) then
        k_rad_lim = 1
        k=k_cloud_top(i,j)-1
        kl=max(1,k)  ! only to avoid out-of-bounds compiler warning
        do while ( k  >   k_rad_lim                                            &
                  .and. dflw_over_cp(i,j,kl)  >   zero                         &
                  .and. z_tq(i,j,kl)  >   one_half*zh(i,j) )
          k_level(i,j) = k
          k = k-1
          kl=max(1,k)

        end do ! k
      end if
    end do ! i
  end do ! j
!$OMP end do

else

  ! original method of finding the k_cloud_dsct, the top of the mixed layer
  ! as seen by radiation, found to be resolution dependent

!$OMP do SCHEDULE(DYNAMIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ! restrict search to `close' to ZH
      k_rad_lim = ntml(i,j)+1
      do k = 1, min(bl_levels,k_rad_lim)
          !-------------------------------------------------------------
          ! Find the layer below K_RAD_LIM with the greatest LW
          ! radiative flux jump in the upper half of the BL
          ! and assume that this is the top of the SML.
          !-------------------------------------------------------------
        if (dflw_over_cp(i,j,k)  >   df_top_over_cp(i,j)                       &
                    .and. z_tq(i,j,k)  >   one_half*zh(i,j) ) then
          k_cloud_top(i,j) = k
          if ( k >  1 ) then
              ! Set K_CLOUD_TOP to the level below if its DF is
              ! greater than half the maximum.  DF in level
              ! K_CLOUD_TOP+1 is then included as DF_INV_SML below.
            if (dflw_over_cp(i,j,k-1)  >   one_half*dflw_over_cp(i,j,k))       &
              k_cloud_top(i,j) = k-1
          end if
          df_top_over_cp(i,j) = dflw_over_cp(i,j,k)
        end if

      end do ! k
    end do ! i

        !-----------------------------------------------------------------
        !  Find bottom grid-level (K_LEVEL) for cloud-top radiative fux
        !  divergence: higher of base of LW radiatively cooled layer,
        !  0.5*ZH, since cooling must be in upper part of layer
        !-----------------------------------------------------------------
    do i = pdims%i_start, pdims%i_end
      k_level(i,j) = k_cloud_top(i,j)
      if ( k_cloud_top(i,j)  >   1 ) then
        k_rad_lim = 1
        k=k_cloud_top(i,j)-1
        kl=max(1,k)  ! only to avoid out-of-bounds compiler warning
        do while ( k  >   k_rad_lim                                            &
                  .and. dflw_over_cp(i,j,kl)  >   zero                         &
                  .and. z_tq(i,j,kl)  >   one_half*zh(i,j) )
          k_level(i,j) = k
          k = k-1
          kl=max(1,k)

        end do
      end if
    end do ! i
  end do ! j
!$OMP end do

end if  ! test on l_new_kcloudtop

      !-----------------------------------------------------------------
      ! Calculate LW and SW flux divergences and combine into
      ! cloud-top turbulence forcing.
      ! Need to account for radiative divergences in cloud in inversion
      ! grid-level, DF_INV_SML. Assume DF_OVER_CP(K_cloud_top+2) is
      ! representative of clear-air rad divergence and so subtract this
      ! `clear-air' part from the grid-level divergence.
      !-----------------------------------------------------------------

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    if ( k_cloud_top(i,j)  >   0 ) then
      dflw_inv = zero
      dfsw_inv = zero
      if ( k_cloud_top(i,j) <  bl_levels ) then
        k = k_cloud_top(i,j)+1
        if ( k  <   bl_levels ) then
          dflw_inv = dflw_over_cp(i,j,k)                                       &
                     - dflw_over_cp(i,j,k+1)                                   &
                            * dzl(i,j,k)/dzl(i,j,k+1)
          dfsw_inv = dfsw_over_cp(i,j,k)                                       &
                     - dfsw_over_cp(i,j,k+1)                                   &
                            * dzl(i,j,k)/dzl(i,j,k+1)
        else
          dflw_inv = dflw_over_cp(i,j,k)
          dfsw_inv = dfsw_over_cp(i,j,k)
        end if
        dflw_inv = max( dflw_inv, zero )
        dfsw_inv = min( dfsw_inv, zero )
      end if
      df_inv_sml(i,j) = dflw_inv + dfsw_inv

      df_top_over_cp(i,j) = frad_lw(i,j,k_cloud_top(i,j)+1)                    &
                           - frad_lw(i,j,k_level(i,j))                         &
                           + dflw_inv

      dfsw_top = frad_sw(i,j,k_cloud_top(i,j)+1)                               &
               - frad_sw(i,j,k_level(i,j))                                     &
               + dfsw_inv

        !-----------------------------------------------------------
        ! Combine SW and LW cloud-top divergences into a net
        ! divergence by estimating SW flux divergence at a given
        ! LW divergence = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
        ! Empirically (from LEM data) a reasonable fit is found
        ! with A small and (1-exp{-A*kappa_sw/kappa_lw}) = dfsw_frac
        !-----------------------------------------------------------
      df_top_over_cp(i,j) = max( zero,                                         &
                    df_top_over_cp(i,j) + dfsw_frac * dfsw_top )
    end if
  end do
end do
!$OMP end do

! ------------------------------------------------------------------
! 9. Calculate the non-turbulent fluxes at the layer boundaries.
!  - the radiative flux at the inversion allows for an estimate
!    of the FA flux divergence within the inversion grid-level
!  - because the radiative time-step is usually longer the radiative
!    cloud-top grid-level (K_CLOUD_TOP) is allowed to differ from
!    the actual one (NTML)
!  - the subsidence flux at the inversion is taken from the
!    flux grid-level below it (assumes the divergence across
!    the inversion is physically above the BL)
!  - the microphysical flux at the inversion is taken from the
!    flux grid-level just above it (assumes the divergence across
!    the inversion grid-level is physically within the BL)
!  - if no rad cooling was identified in layer, need to set
!    K_CLOUD_TOP and K_CLOUD_DSCT to top level in mixed layer
! ------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( k_cloud_top(i,j)  ==  0 ) k_cloud_top(i,j) = ntml(i,j)

    ft_nt_zh(i,j)   = frad(i,j,k_cloud_top(i,j)+1)                             &
                         + df_inv_sml(i,j)
    ft_nt_zh(i,j)   = ft_nt_zh(i,j)   + fmic(i,j,ntml(i,j)+2,1)                &
                                      + fsubs(i,j,ntml(i,j),1)
    fq_nt_zh(i,j)   = fmic(i,j,ntml(i,j)+2,2)                                  &
                    + fsubs(i,j,ntml(i,j),2)

    ft_nt_zhsc(i,j) = zero
    ft_nt_zhsc(i,j) = zero
    fq_nt_zhsc(i,j) = zero
    if ( dsc(i,j) ) then
      if ( k_cloud_dsct(i,j)  ==  0 ) k_cloud_dsct(i,j) = ntdsc(i,j)
      ft_nt_zhsc(i,j) = frad(i,j,k_cloud_dsct(i,j)+1)                          &
                           + df_inv_dsc(i,j)
      ft_nt_zhsc(i,j) = ft_nt_zhsc(i,j) + fmic(i,j,ntdsc(i,j)+2,1)             &
                                        + fsubs(i,j,ntdsc(i,j),1)
      fq_nt_zhsc(i,j) = fmic(i,j,ntdsc(i,j)+2,2)                               &
                      + fsubs(i,j,ntdsc(i,j),2)
    end if
  end do
end do
!$OMP end do NOWAIT

! Repeat for water tracers
if (l_wtrac) then
  do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        fq_nt_zh_wtrac(i,j,i_wt) = fmic_wtrac(i,j,ntml(i,j)+2,i_wt)            &
                                   + fsubs_wtrac(i,j,ntml(i,j),i_wt)

        fq_nt_zhsc_wtrac(i,j,i_wt) = zero
        if ( dsc(i,j) ) then
          fq_nt_zhsc_wtrac(i,j,i_wt) = fmic_wtrac(i,j,ntdsc(i,j)+2,i_wt)       &
                                       + fsubs_wtrac(i,j,ntdsc(i,j),i_wt)
        end if
      end do
    end do
!$OMP end do NOWAIT
  end do
end if

!-----------------------------------------------------------------------
! 10.  Subroutine EXCF_NL.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ntml_save(i,j) = ntml(i,j)  ! needed to identify changes
    dsc_save(i,j)  = dsc(i,j)   !      in excf_nl
    dsc_removed(i,j) = 0
  end do
end do
!$OMP end do NOWAIT

!$OMP end PARALLEL

call excf_nl_9c (                                                              &
! in levels/switches
   bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                                     &
! in fields
   rdz,z_uv,z_tq,rho_mix,rho_wet_tq,rhostar_gb,v_s,fb_surf,zhpar,zh_prev,z_lcl,&
   btm,bqm,btm_cld,bqm_cld,cf_sml,cf_dsc,                                      &
   bflux_surf,bflux_surf_sat,zeta_s,svl_diff_frac,                             &
   df_top_over_cp,zeta_r,bt_top,btt_top,btc_top,                               &
   db_top,db_top_cld,chi_s_top,br_fback,                                       &
   df_dsct_over_cp,zeta_r_dsc,bt_dsct,btt_dsct,                                &
   db_dsct,db_dsct_cld,chi_s_dsct,br_fback_dsc,                                &
   db_ga_dry,db_noga_dry,db_ga_cld,db_noga_cld,                                &
   dsl_sml,dqw_sml,dsl_dsc,dqw_dsc,ft_nt,fq_nt,sl,qw,                          &
! INOUT fields
   dsc_removed,dsc,cumulus,coupled,ntml,zh,zhsc,dscdepth,ntdsc,zc,zc_dsc,      &
   ft_nt_zh, ft_nt_zhsc, fq_nt_zh, fq_nt_zhsc, dzh,                            &
   fq_nt_zh_wtrac, fq_nt_zhsc_wtrac,                                           &
! out fields
   rhokm, rhokh, rhokm_top, rhokh_top,                                         &
   rhokh_top_ent, rhokh_dsct_ent, rhokh_surf_ent,                              &
   rhof2,rhofsc,f_ngstress,tke_nl,zdsc_base,nbdsc                              &
  )

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP  private (i, j, i_wt, k, kl, kp, c_ws, c_tke, w_m, tothf_efl, totqf_efl, &
!$OMP  ml_tend, fa_tend, inv_tend, Prandtl, svl_lapse_rho,                     &
!$OMP  recip_svl_lapse, rhok_inv, svl_lapse, svl_target, svl_flux)

!-----------------------------------------------------------------------
!-adjust SML/DSC properties depending on diagnoses in EXCF_NL
! Note that the non-turbulent fluxes at inversions will have been
! swapped in EXCF_NL (ie. FT/Q_NT_ZH/ZHSC)
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( dsc(i,j) .and. .not. dsc_save(i,j) ) then
      !..decoupling diagnosed in EXCF_NL - change parameters around
      dsl_dsc(i,j) = dsl_sml(i,j)
      dqw_dsc(i,j) = dqw_sml(i,j)
      db_dsct(i,j) = db_top(i,j)  ! copy diagnostics across
      zc_dsc(i,j)  = zc(i,j)      !
      dsc_disc_inv(i,j) = sml_disc_inv(i,j)
      sml_disc_inv(i,j) = 0
      dsl_sml(i,j) = zero
      dqw_sml(i,j) = zero
      df_inv_dsc(i,j) = df_inv_sml(i,j)
      df_inv_sml(i,j) = zero
      res_inv(i,j) = 0  ! not clear what to do at dsc top (dzh diagnosed
                        ! assuming well-mixed) so do nothing for now!
      dzh(i,j) = rmdi   ! don't want to associate dzh with new zh
    end if
  end do
end do
!$OMP end do

! Repeat for water tracers
if (l_wtrac) then
  do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( dsc(i,j) .and. .not. dsc_save(i,j) ) then
          !..decoupling diagnosed in EXCF_NL - change parameters around
          dqw_dsc_wtrac(i,j,i_wt) = dqw_sml_wtrac(i,j,i_wt)
          dqw_sml_wtrac(i,j,i_wt) = zero
        end if
      end do
    end do
!$OMP end do
  end do
end if

if ( l_use_sml_dsc_fixes ) then

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( dsc_removed(i,j) == 1 ) then
        ! decoupled layer removed in EXCF_NL because it had no
        ! turbulence forcing
        dsl_dsc(i,j) = zero
        dqw_dsc(i,j) = zero
        dsc_disc_inv(i,j) = 0
        df_inv_dsc(i,j) = zero
      else if ( dsc_removed(i,j) == 2 ) then
        ! decoupled layer recoupled with surface layer
        dsl_sml(i,j) = dsl_dsc(i,j)
        dqw_sml(i,j) = dqw_dsc(i,j)
        dsl_dsc(i,j) = zero
        dqw_dsc(i,j) = zero
        sml_disc_inv(i,j) = dsc_disc_inv(i,j)
        dsc_disc_inv(i,j) = 0
        df_inv_sml(i,j) = df_inv_dsc(i,j)
        df_inv_dsc(i,j) = zero
        res_inv(i,j) = 0  ! dzh was diagnosed for a capping inversion
        dzh(i,j) = rmdi   ! subsumed into mixed layer so remove its flags
      end if
    end do
  end do
!$OMP end do

  ! Repeat for water tracers
  if (l_wtrac) then
    do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          if ( dsc_removed(i,j) == 1 ) then
            ! decoupled layer removed in EXCF_NL because it had no
            ! turbulence forcing
            dqw_dsc_wtrac(i,j,i_wt) = zero
          else if ( dsc_removed(i,j) == 2 ) then
            ! decoupled layer recoupled with surface layer
            dqw_sml_wtrac(i,j,i_wt) = dqw_dsc_wtrac(i,j,i_wt)
            dqw_dsc_wtrac(i,j,i_wt) = zero
          end if
        end do
      end do
!$OMP end do
    end do
  end if

else ! not l_use_sml_dsc_fixes

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( .not. dsc(i,j) .and. dsc_save(i,j) ) then
        !..decoupled layer removed in EXCF_NL; either...
        if ( ntml_save(i,j)  ==  ntml(i,j) ) then
          !...had no turbulence forcing
          dsl_dsc(i,j) = zero
          dqw_dsc(i,j) = zero
          dsc_disc_inv(i,j) = 0
          df_inv_dsc(i,j) = zero
        else
          !...recoupled with surface layer
          dsl_sml(i,j) = dsl_dsc(i,j)
          dqw_sml(i,j) = dqw_dsc(i,j)
          dsl_dsc(i,j) = zero
          dqw_dsc(i,j) = zero
          sml_disc_inv(i,j) = dsc_disc_inv(i,j)
          dsc_disc_inv(i,j) = 0
          df_inv_sml(i,j) = df_inv_dsc(i,j)
          df_inv_dsc(i,j) = zero
          res_inv(i,j) = 0  ! dzh was diagnosed for a capping inversion
          dzh(i,j) = rmdi   ! subsumed into mixed layer so remove its flags
        end if
      end if
    end do
  end do
!$OMP end do

  ! Repeat for water tracers
  if (l_wtrac) then
    do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          if ( .not. dsc(i,j) .and. dsc_save(i,j) ) then
            !..decoupled layer removed in EXCF_NL; either...
            if ( ntml_save(i,j)  ==  ntml(i,j) ) then
              !...had no turbulence forcing
              dqw_dsc_wtrac(i,j,i_wt) = zero
            else
              !...recoupled with surface layer
              dqw_sml_wtrac(i,j,i_wt) = dqw_dsc_wtrac(i,j,i_wt)
              dqw_dsc_wtrac(i,j,i_wt) = zero
            end if
          end if
        end do
      end do
!$OMP end do
    end do
  end if

end if ! not l_use_sml_dsc_fixes

!-----------------------------------------------------------------------
! 11.  Calculate "explicit" entrainment fluxes of SL and QW.
!-----------------------------------------------------------------------
! Calculate the non-turbulent fluxes at the DSC base
! ------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ft_nt_dscb(i,j) = ft_nt(i,j,1)
    if ( nbdsc(i,j)  >   1 ) then
      k = nbdsc(i,j)  ! NBDSC marks the lowest flux-level
                         !    within the DSC layer
                         ! Interpolate non-turb flux to base
                         !    of DSC layer:
      ft_nt_dscb(i,j) = ft_nt(i,j,k-1) +                                       &
                (ft_nt(i,j,k)-ft_nt(i,j,k-1))                                  &
               *(zdsc_base(i,j)-z_uv(i,j,k-1))/dzl(i,j,k-1)
    end if
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
!..Specify entrainment fluxes at NTML+1 and NTDSC+1 directly through FTL
!..and FQW (and set the entrainment RHOKH to zero).
!..The turbulent flux at the subgrid ZH is given by the entrainment
!..parametrization ( = -w_e*'jump across inversion').  Together with
!..the non-turbulent flux profile (rad+microphys+subs), this gives the
!..total flux at the subgrid ZH.  The linear total flux profile is then
!..interpolated onto the half-level below ZH (the entrainment flux
!..grid-level).  The total flux divergence across the inversion grid
!..level is then checked for consistency with the entrainment/subsidence
!..balance, eg. a falling inversion should warm the inversion grid-level
!..Finally, the non-turbulent flux is sutracted off to give the required
!..turbulent component at the entrainment flux grid-level.
!..For momentum, given the horizontal interpolation required, together
!..with the lack of accuracy in assuming a discontinuous inversion,
!..entrainment continues to be specified using the specified RHOKM.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!..First the surface-based mixed layer (if entraining and a
!..discontinuous inversion structure was diagnosed -
!..ie. the inversion is well-defined)
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    zh_np1(i,j)   = zero
    t_frac(i,j)   = zero
    zrzi(i,j)     = zero
    we_rho(i,j)   = zero
    tothf_zh(i,j) = zero
    totqf_zh(i,j) = zero
    k=ntml(i,j)+1

      ! Only RHOKH_ENT is passed out of EXCFNL so recalculate WE:
    we_parm(i,j) = rdz(i,j,k)*                                                 &
                       ( rhokh_top_ent(i,j)+rhokh_surf_ent(i,j) )              &
                                                  / rho_mix(i,j,k)

    if ( sml_disc_inv(i,j)  ==  1 .and. .not. coupled(i,j) .and.               &
         (rhokh_top_ent(i,j)+rhokh_surf_ent(i,j))  >   zero ) then

      !-----------------------------------------------------------------
      !..Calculate ZH at end of timestep, ZH_NP1
      !-----------------------------------------------------------------
      !..linearly interpolate vertical velocity to ZH
      if ( zh(i,j)  >=  z_tq(i,j,k) ) then
        w_ls(i,j) = w(i,j,k) + ( w(i,j,k+1) - w(i,j,k) )                       &
                    * (zh(i,j)-z_tq(i,j,k)) * rdz(i,j,k+1)
      else
        w_ls(i,j) = w(i,j,k) + ( w(i,j,k) - w(i,j,k-1) )                       &
                    * (zh(i,j)-z_tq(i,j,k)) * rdz(i,j,k)
      end if
      w_ls(i,j) = min ( w_ls(i,j), zero )
        ! only interested in subsidence

      zh_np1(i,j) = zh(i,j) +                                                  &
                      timestep * ( we_parm(i,j) + w_ls(i,j) )
      zh_np1(i,j) = max( zh_np1(i,j), z_uv(i,j,k-1) )
      if ( zh_np1(i,j)  >   z_top(i,j,k+1) ) then
          ! limit ZH and W_e (and therefore the entraiment fluxes)
          ! because the inversion cannot rise more than one level
          ! in a timestep.
        zh_np1(i,j) = z_top(i,j,k+1)
        we_parm(i,j) =                                                         &
                (z_top(i,j,k+1) - zh(i,j))/timestep - w_ls(i,j)
      end if
      !-----------------------------------------------------------------
      !..Decide on which grid-level to apply entrainment flux
      !-----------------------------------------------------------------
      if ( zh_np1(i,j)  >   z_uv(i,j,ntml(i,j)+2) ) then
          ! ZH risen above level K+1 so specify appropriate flux
          ! at this level and raise NTML by one (this means
          ! gradient-adjustment is also applied at half-level
          ! old_NTML+1).  Note KH profiles should already be
          ! calculated at level NTML+1 because ZH is above this level.
        ntml(i,j) = ntml(i,j) + 1
        k=ntml(i,j)+1
        sml_disc_inv(i,j)=2

          ! T_FRAC is fraction of timestep inversion is above
          ! the entrainment flux grid-level (at Z_UV(K))
        t_frac(i,j) = (zh_np1(i,j)-z_uv(i,j,k)) /                              &
                      (zh_np1(i,j)-zh(i,j))
          ! ZH_FRAC is the timestep-average fraction of mixed layer
          ! air in the inversion grid-level, level NTML+1
        zh_frac(i,j) = one_half*t_frac(i,j)*(zh_np1(i,j)-z_uv(i,j,k) )         &
                       / dzl(i,j,k)

      else if ( zh_np1(i,j)  >=  z_uv(i,j,ntml(i,j)+1) ) then
          ! ZH always between half-levels NTML+1 and NTML+2

        t_frac(i,j) = one
        zh_frac(i,j) = ( one_half*(zh(i,j)+zh_np1(i,j)) - z_uv(i,j,k) )        &
                       / dzl(i,j,k)

      else
          ! ZH falls below half-level NTML+1
          ! Keep implicit (diffusive) entrainment but apply
          ! at the level below
        if (ntml(i,j)  >=  2) then     ! ftl(k=1) is surface flux
          ntml(i,j) = ntml(i,j) - 1
          k=ntml(i,j)+1
          rhokh_top(i,j,k+1) = zero   ! also need to remove diffusion
          rhokh(i,j,k+1)     = zero   ! at old entrainment grid-level
        end if
        t_frac(i,j) = zero
        zh_frac(i,j) = zero
        sml_disc_inv(i,j) = 0

      end if  ! test on where to apply entrainment flux

      we_rho(i,j) = rho_mix(i,j,k) * we_parm(i,j)
      zrzi(i,j)   = z_uv(i,j,k)*2.0_r_bl/(zh(i,j)+zh_np1(i,j))

    end if   ! test on SML_DISC_INV, etc
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
!..Linearly interpolate between the known total (turb+rad+subs+micro)
!..flux at the surface and the parametrized flux at the inversion
!-----------------------------------------------------------------------
if (flux_grad  ==  LockWhelan2006) then
  c_ws = 0.42_r_bl
else if (flux_grad  ==  HoltBov1993) then
  c_ws = 0.6_r_bl
else
  c_ws = 0.25_r_bl
end if
c_tke = 1.33_r_bl/(vkman*c_ws**two_thirds)

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

      ! Entrainment flux applied to level NTML+1 which is the
      ! flux-level above the top of the SML
    k=ntml(i,j)+1

    if ( t_frac(i,j)  >   zero ) then

      rhokh_top(i,j,k) = zero   ! apply entrainment explicitly
      rhokh(i,j,k) = zero       !      "

      tothf_zh(i,j) = - we_rho(i,j)*dsl_sml(i,j) + ft_nt_zh(i,j)
        ! Linearly interpolate to entrainment flux grid-level
      tothf_efl = ft_nt(i,j,1) + ftl(i,j,1) +                                  &
                 ( tothf_zh(i,j)-ft_nt(i,j,1)-ftl(i,j,1) )*zrzi(i,j)
        ! Ensure total heat flux gradient in inversion grid-level is
        ! consistent with inversion rising (ie. implies cooling in
        ! level K relative to the mixed layer) or falling
        ! (implies warming)

      ml_tend = -( tothf_zh(i,j)-ft_nt(i,j,1)-ftl(i,j,1) ) / zh(i,j)
      fa_tend = zero
      if ( k+1  <=  bl_levels )                                                &
          fa_tend = - ( ft_nt(i,j,k+2) - ft_nt(i,j,k+1) )                      &
                      / dzl(i,j,k+1)
      inv_tend =       zh_frac(i,j) * ml_tend                                  &
               + (one-zh_frac(i,j)) * fa_tend
      if (we_parm(i,j)+w_ls(i,j)  >=  zero) then
          ! Inversion moving up so inversion level should cool
          ! Ensure it does cool relative to ML
        tothf_efl = min( tothf_efl,                                            &
                         ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
          ! Ensure inversion level won't end up colder than
          ! NTML by end of timestep.
          ! Set INV_TEND to max allowable cooling rate, also
          ! allowing for change in ML_TEND arising from this change
          ! to TOTHF_EFL:
        inv_tend = (sl(i,j,k-1)-sl(i,j,k))/timestep                            &
                    + (ft_nt(i,j,1)+ftl(i,j,1))/z_uv(i,j,k)
        tothf_efl = max( tothf_efl,                                            &
                         (ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                  &
                              /(one+ dzl(i,j,k)/z_uv(i,j,k)) )
      else  ! WE_PARM+W_LS < 0
          ! Ensure inversion level does warm relative to ML
        tothf_efl = max( tothf_efl,                                            &
                         ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
      end if
        ! Turbulent entrainment flux is then the residual of the total
        ! flux and the net flux from other processes
      ftl(i,j,k) =  t_frac(i,j) * ( tothf_efl - ft_nt(i,j,k) )
    else   ! not specifying entrainment flux but KH
        ! Include entrainment KH in K-profiles, if greater
        ! (for COUPLED layers these will be zero)
      rhokh_top(i,j,k) = max( rhokh_top(i,j,k), rhokh_top_ent(i,j) )
      rhokh(i,j,k)     = max( rhokh(i,j,k), rhokh_surf_ent(i,j) )

      if (res_inv(i,j) == 1) then
        Prandtl = min( rhokm(i,j,k)/(rbl_eps+rhokh_surf_ent(i,j)),             &
                       pr_max )
        if (BL_diag%l_tke) then
          ! need velocity scale for TKE diagnostic
          w_m = ( v_s(i,j)*v_s(i,j)*v_s(i,j) +                                 &
                    c_ws * zh(i,j) * fb_surf(i,j) ) ** one_third
        end if

        if (bl_res_inv == cosine_inv_flux) then
          svl_lapse_rho = (svl(i,j,k)-svl(i,j,k-1)) /                          &
                          ( (z_tq(i,j,k)-z_tq(i,j,k-1))*rho_mix(i,j,k) )
          kl=k+1
          do while ( z_uv(i,j,kl) < zh(i,j)+dzh(i,j) .and.                     &
                     kl <= bl_levels )
            recip_svl_lapse = (z_tq(i,j,kl)-z_tq(i,j,kl-1))/                   &
                              max( 0.01_r_bl, svl(i,j,kl)-svl(i,j,kl-1) )
            rhok_inv = rhokh_surf_ent(i,j) * svl_lapse_rho *                   &
                   rho_mix(i,j,kl) * recip_svl_lapse *                         &
                   cos(one_half*pi*(z_uv(i,j,kl)-zh(i,j))/dzh(i,j))
            rhok_inv = min( rhok_inv, 1000.0_r_bl )
            rhokh(i,j,kl) = max( rhokh(i,j,kl), rhok_inv )
            ! rescale for KM on staggered grid
            rhok_inv =  Prandtl * rhok_inv                                     &
                       * rdz(i,j,kl) * (z_uv(i,j,kl)-z_uv(i,j,kl-1))           &
                       * rho_wet_tq(i,j,kl-1) / rho_mix(i,j,kl)
            rhokm(i,j,kl) = max( rhokm(i,j,kl), rhok_inv )
            if (BL_diag%l_tke) then
              ! save Km/timescale for TKE diag, completed in bdy_expl2
              tke_nl(i,j,kl) = max( tke_nl(i,j,kl), rhok_inv*c_tke*w_m/zh(i,j))
            end if
            kl=kl+1
          end do
        else if (bl_res_inv == target_inv_profile) then
          svl_lapse = (svl(i,j,k)-svl(i,j,k-1)) /                              &
                      ( (z_tq(i,j,k)-z_tq(i,j,k-1)) )
          kp=k+1  ! kp marks the lowest level above the inversion
          do while ( z_uv(i,j,kp) < zh(i,j)+dzh(i,j) .and.                     &
                     kp <= bl_levels )
            kp=kp+1
          end do
          svl_flux(k) = - rhokh_surf_ent(i,j) * svl_lapse
          kl=k+1
          do while ( z_uv(i,j,kl) < zh(i,j)+dzh(i,j) .and.                     &
                     kl <= bl_levels )
            ! assume a linear target svl profile within inversion
            svl_target = svl(i,j,k-1) + (svl(i,j,kp)-svl(i,j,k-1)) *           &
                                            (z_uv(i,j,kl)-zh(i,j)) / dzh(i,j)
            rho_dz = rho_mix_tq(i,j,kl) * dzl(i,j,kl)
            svl_flux(kl) = svl_flux(kl-1) -                                    &
                                (svl_target-svl(i,j,kl))*rho_dz/timestep
            kl=kl+1
          end do
          ! linearly extrapolate flux to inversion top
          svl_flux(kp)=svl_flux(kp-1) + (svl_flux(kp-1)-svl_flux(kp-2))*       &
                                (zh(i,j)+dzh(i,j)-z_uv(i,j,kp-1))*rdz(i,j,kp-1)
          kl=k+1
          do while ( z_uv(i,j,kl) < zh(i,j)+dzh(i,j) .and.                     &
                     kl <= bl_levels )
            ! rescale svl_flux so as to have zero flux at the inversion top
            ! ie so svl_flux(kp)=0
            svl_flux(kl) = svl_flux(k)*( one -                                 &
                                      (svl_flux(kl)-svl_flux(k))/              &
                                      (svl_flux(kp)-svl_flux(k)) )
            recip_svl_lapse = (z_tq(i,j,kl)-z_tq(i,j,kl-1))/                   &
                              max( 0.01_r_bl, svl(i,j,kl)-svl(i,j,kl-1) )
            rhok_inv = - svl_flux(kl) * recip_svl_lapse

            rhok_inv = min( rhok_inv, 1000.0_r_bl )
            rhokh(i,j,kl) = max( rhokh(i,j,kl), rhok_inv )
            ! rescale for KM on staggered grid
            rhok_inv =  Prandtl * rhok_inv                                     &
                       * rdz(i,j,kl) * (z_uv(i,j,kl)-z_uv(i,j,kl-1))           &
                       * rho_wet_tq(i,j,kl-1) / rho_mix(i,j,kl)
            rhokm(i,j,kl) = max( rhokm(i,j,kl), rhok_inv )
            if (BL_diag%l_tke) then
              ! save Km/timescale for TKE diag, completed in bdy_expl2
              tke_nl(i,j,kl) = max( tke_nl(i,j,kl), rhok_inv*c_tke*w_m/zh(i,j))
            end if
            kl=kl+1
          end do
        end if  ! bl_res_inv option
      end if  ! res_inv
    end if  ! test on T_FRAC gt 0

  end do
end do
!$OMP end do NOWAIT

!-------------------------------------------------
!..Second the decoupled mixed layer, if entraining
!-------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    zhsc_np1(i,j)   = zero
    t_frac_dsc(i,j) = zero
    zrzi_dsc(i,j)   = zero
    we_rho_dsc(i,j) = zero
    tothf_zhsc(i,j) = zero
    totqf_zhsc(i,j) = zero

    k=ntdsc(i,j)+1
    we_dsc_parm(i,j) = rdz(i,j,k)*rhokh_dsct_ent(i,j)                          &
                                               / rho_mix(i,j,k)

    if ( dsc_disc_inv(i,j)  ==  1                                              &
                  .and. rhokh_dsct_ent(i,j)  >   zero ) then

      !-----------------------------------------------------------------
      !..Calculate ZHSC at end of timestep, ZHSC_NP1
      !-----------------------------------------------------------------
      !..interpolate vertical velocity to ZH
      if ( zhsc(i,j)  >=  z_tq(i,j,k) ) then
        w_ls_dsc(i,j) = w(i,j,k) + ( w(i,j,k+1) - w(i,j,k) ) *                 &
                         (zhsc(i,j)-z_tq(i,j,k)) * rdz(i,j,k+1)
      else
        w_ls_dsc(i,j) = w(i,j,k) + ( w(i,j,k) - w(i,j,k-1) ) *                 &
                         (zhsc(i,j)-z_tq(i,j,k)) * rdz(i,j,k)
      end if
      w_ls_dsc(i,j) = min ( w_ls_dsc(i,j), zero )
        ! only interested in subsidence

      zhsc_np1(i,j) = zhsc(i,j) +                                              &
            timestep * ( we_dsc_parm(i,j) + w_ls_dsc(i,j) )
      zhsc_np1(i,j) = max( zhsc_np1(i,j), z_uv(i,j,k-1) )
      if ( zhsc_np1(i,j)  >   z_top(i,j,k+1) ) then
          ! limit ZHSC and W_e (and therefore the entrainment fluxes)
          ! because the inversion cannot rise more than one level
          ! in a timestep.
        zhsc_np1(i,j) = z_top(i,j,k+1)
        we_dsc_parm(i,j) =                                                     &
           (z_top(i,j,k+1) - zhsc(i,j))/timestep - w_ls_dsc(i,j)
      end if
      !-----------------------------------------------------------------
      !..Decide on which grid-level to apply entrainment flux
      !-----------------------------------------------------------------
      if ( zhsc_np1(i,j)  >   z_uv(i,j,ntdsc(i,j)+2) ) then
          ! ZHSC risen above level K+1 so specify appropriate
          ! flux at this level and raise NTDSC by one

        ntdsc(i,j) = ntdsc(i,j) + 1
        k = ntdsc(i,j)+1
        dsc_disc_inv(i,j) = 2
        t_frac_dsc(i,j) = (zhsc_np1(i,j)-z_uv(i,j,k)) /                        &
                          (zhsc_np1(i,j)-zhsc(i,j))

        zhsc_frac(i,j) = one_half*t_frac_dsc(i,j)*                             &
                         ( zhsc_np1(i,j)-z_uv(i,j,k) )/ dzl(i,j,k)

      else if ( zhsc_np1(i,j)  >   z_uv(i,j,ntdsc(i,j)+1) ) then
          ! ZHSC always between half-levels NTDSC+1 and NTDSC+2

        t_frac_dsc(i,j) = one
        zhsc_frac(i,j) = ( one_half*(zhsc(i,j)+zhsc_np1(i,j))                  &
                                       - z_uv(i,j,k) )/ dzl(i,j,k)

      else
          ! ZHSC falls below half-level NTDSC+1
          ! Keep implicit (diffusive) entrainment but apply
          ! at the level below
        ntdsc(i,j) = ntdsc(i,j) - 1  ! could reduce NTDSC to 1
        k = ntdsc(i,j)+1
        rhokh_top(i,j,k+1) = zero
        rhokh(i,j,k+1)     = zero

        t_frac_dsc(i,j)   = zero
        zhsc_frac(i,j)    = zero
        dsc_disc_inv(i,j) = 0

      end if  ! test on where to apply entrainment flux

      we_rho_dsc(i,j) = rho_mix(i,j,k) * we_dsc_parm(i,j)
        ! for z'/z_i' assume height of DSC base is fixed in time
      zrzi_dsc(i,j) =( z_uv(i,j,k)-(zhsc(i,j)-dscdepth(i,j)) )                 &
                    /( dscdepth(i,j)+one_half*(zhsc_np1(i,j)-zhsc(i,j)) )

    end if   ! test on DSC_DISC_INV, etc
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
!..Linearly interpolate between the known total (turb+rad+subs+micro)
!..flux at the DSC base and the parametrized flux at the inversion
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

      ! Entrainment flux applied to level NTDSC+1 which is the
      ! flux-level above the top of the DSC layer
    k=ntdsc(i,j)+1

    if ( t_frac_dsc(i,j)  >   zero ) then

      rhokh_top(i,j,k) = zero   ! apply entrainment explicitly
      rhokh(i,j,k)     = zero   !      "

      tothf_zhsc(i,j) = - we_rho_dsc(i,j)*dsl_dsc(i,j)                         &
                              + ft_nt_zhsc(i,j)
      tothf_efl = ft_nt_dscb(i,j) +                                            &
                  ( tothf_zhsc(i,j)-ft_nt_dscb(i,j) )*zrzi_dsc(i,j)
        ! Ensure total heat flux gradient in inversion grid-level is
        ! consistent with inversion rising (implies cooling in
        ! level K, relative to mixed layer) or falling
        ! (implies warming)
      ml_tend = - ( tothf_zhsc(i,j)-ft_nt_dscb(i,j) )/ dscdepth(i,j)
      fa_tend = zero
      if ( k+1  <=  bl_levels )                                                &
          fa_tend = - ( ft_nt(i,j,k+2) - ft_nt(i,j,k+1) )                      &
                      / dzl(i,j,k+1)
      inv_tend =       zhsc_frac(i,j) * ml_tend                                &
               + (one-zhsc_frac(i,j)) * fa_tend

      if (we_dsc_parm(i,j)+w_ls_dsc(i,j)  >=  zero) then
          ! Inversion moving up so inversion level should cool
          ! Ensure it does cool relative to ML
        tothf_efl = min( tothf_efl,                                            &
                         ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
          ! Ensure inversion level won't end up colder than
          ! NTDSC by end of timestep.
        inv_tend = (sl(i,j,k-1)-sl(i,j,k))/timestep                            &
                   + ft_nt_dscb(i,j)/dscdepth(i,j)
        tothf_efl = max( tothf_efl,                                            &
                      (ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                     &
                       /(one+ dzl(i,j,k)/dscdepth(i,j))   )
      else   ! WE_DSC_PARM+W_LS_DSC < 0
        tothf_efl = max( tothf_efl,                                            &
                         ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
      end if
        ! Turbulent entrainment flux is then the residual of the total
        ! flux and the net flux from other processes
      ftl(i,j,k) = t_frac_dsc(i,j) * ( tothf_efl - ft_nt(i,j,k) )

    else if ( dsc(i,j) ) then

        ! Not specifying entrainment flux but KH
        ! Include entrainment KH in K-profile, if greater
      rhokh_top(i,j,k) = max( rhokh_top(i,j,k),rhokh_dsct_ent(i,j) )

    end if  ! if not DSC

  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! Specify QW entrainment fluxes
!-----------------------------------------------------------------------
! Calculate the non-turbulent fluxes at the layer boundaries
!  - the subsidence flux at the inversion is taken from the
!    flux grid-level below it (assumes the divergence across
!    the inversion is physically above the BL)
!  - the microphysical flux at the inversion is taken from the
!    flux grid-level just above it (assumes the divergence across
!    the inversion grid-level is physically within the BL)
! ------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    fq_nt_dscb(i,j) = fq_nt(i,j,1)
    if ( nbdsc(i,j)  >   1 ) then
      k = nbdsc(i,j)  ! NBDSC marks the lowest flux-level
                         !    within the DSC layer
                         ! Interpolate non-turb flux to base
                         !    of DSC layer:
      fq_nt_dscb(i,j) = fq_nt(i,j,k-1) +                                       &
                (fq_nt(i,j,k)-fq_nt(i,j,k-1))                                  &
               *(zdsc_base(i,j)-z_uv(i,j,k-1))/dzl(i,j,k-1)
    end if
  end do
end do
!$OMP end do NOWAIT

! Repeat for water tracers
if (l_wtrac) then
  do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        wtrac_bl(i_wt)%fq_nt_dscb(i,j) = wtrac_bl(i_wt)%fq_nt(i,j,1)
        if ( nbdsc(i,j)  >   1 ) then
          k = nbdsc(i,j)  ! NBDSC marks the lowest flux-level
                          !    within the DSC layer
                          ! Interpolate non-turb flux to base
                          !    of DSC layer:
          wtrac_bl(i_wt)%fq_nt_dscb(i,j) =  wtrac_bl(i_wt)%fq_nt(i,j,k-1) +    &
                ( wtrac_bl(i_wt)%fq_nt(i,j,k)- wtrac_bl(i_wt)%fq_nt(i,j,k-1))  &
                  *(zdsc_base(i,j)-z_uv(i,j,k-1))/dzl(i,j,k-1)
        end if
      end do
    end do
!$OMP end do
  end do
end if   ! l_wtrac

!-----------------------------------------------------------------------
! Calculate grid-level QW fluxes at inversion, ensuring the turbulent,
! microphysical and subsidence fluxes are correctly coupled.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    moisten(i,j) = .false.
    if (l_wtrac) totqf_efl_meth1(i,j) = 0
    if (l_wtrac) totqf_efl_meth2(i,j) = 0
    k = ntml(i,j)+1

    if ( t_frac(i,j)  >   zero ) then
        ! Calculate total (turb+micro+subs) QW flux at subgrid
        ! inversion height
      totqf_zh(i,j) = - we_rho(i,j)*dqw_sml(i,j) + fq_nt_zh(i,j)
        ! Interpolate to entrainment flux-level below
      totqf_efl = fq_nt(i,j,1) + fqw(i,j,1) +  zrzi(i,j) *                     &
                       ( totqf_zh(i,j) - fq_nt(i,j,1) - fqw(i,j,1) )
        ! Need to ensure the total QW flux gradient in inversion
        ! grid-level is consistent with inversion rising or falling.
        ! If QW(K) is drier than mixed layer then inversion rising
        ! implies moistening in level K relative to mixed layer
        ! while falling would imply relative drying of level K.
        ! If QW(K) is moister than ML then want opposite tendencies.
      ml_tend = - ( totqf_zh(i,j)-fq_nt(i,j,1)-fqw(i,j,1) ) /zh(i,j)
      fa_tend = zero
      if ( k+1  <=  bl_levels )                                                &
        fa_tend = - ( fq_nt(i,j,k+2)-fq_nt(i,j,k+1) )                          &
                    / dzl(i,j,k+1)
      inv_tend =       zh_frac(i,j) * ml_tend                                  &
               + (one-zh_frac(i,j)) * fa_tend

      if (we_parm(i,j)+w_ls(i,j) >=  zero) then
          ! inversion moving up so inversion will moisten/dry
          ! depending on relative QW in level below
        moisten(i,j) = ( qw(i,j,k) <= qw(i,j,k-1) )
      else
          ! inversion moving down so inversion will moisten/dry
          ! depending on relative QW in level above
        moisten(i,j) = ( qw(i,j,k) <= qw(i,j,k+1) )
      end if

      if ( moisten(i,j) ) then
          ! Ensure inversion level does moisten relative to ML

        if (l_wtrac .and. totqf_efl < (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k)) )   &
           totqf_efl_meth1(i,j) = 1         ! Store method

        totqf_efl = max( totqf_efl,                                            &
                        fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
        if (we_parm(i,j)+w_ls(i,j)  >=  zero) then
            ! Ensure inversion level won't end up more moist than
            ! NTML by end of timestep.
            ! Set INV_TEND to max allowable moistening rate, also
            ! allowing for change in ML_TEND arising from this change
            ! to TOTQF_EFL:
          inv_tend = (qw(i,j,k-1)-qw(i,j,k))/timestep                          &
                      + (fq_nt(i,j,1)+fqw(i,j,1))/z_uv(i,j,k)

          if (l_wtrac .and. totqf_efl >                                        &
                        ( (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                 &
                         /(one+ dzl(i,j,k)/z_uv(i,j,k)) ) )                    &
             totqf_efl_meth2(i,j) = 1       ! Store method

          totqf_efl = min( totqf_efl,                                          &
                        (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                   &
                         /(one+ dzl(i,j,k)/z_uv(i,j,k))   )
        end if
      else
        if (l_wtrac .and. totqf_efl > (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k)) )   &
          totqf_efl_meth1(i,j) = 1          ! Store method

        totqf_efl = min( totqf_efl,                                            &
                        fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
        if (we_parm(i,j)+w_ls(i,j)  >=  zero) then
            ! Ensure inversion level won't end up drier than
            ! NTML by end of timestep.
            ! Set INV_TEND to max allowable drying rate:
          inv_tend = (qw(i,j,k-1)-qw(i,j,k))/timestep                          &
                      + (fq_nt(i,j,1)+fqw(i,j,1))/z_uv(i,j,k)

          if (l_wtrac .and. totqf_efl <                                        &
                      ( (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                   &
                         /(one+ dzl(i,j,k)/z_uv(i,j,k)) )  )                   &
            totqf_efl_meth2(i,j) = 1          ! Store method

          totqf_efl = max( totqf_efl,                                          &
                        (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                   &
                         /(one+ dzl(i,j,k)/z_uv(i,j,k))   )
        end if
      end if
      fqw(i,j,k) = t_frac(i,j) *                                               &
                       ( totqf_efl - fq_nt(i,j,k) )
    end if

  end do
end do
!$OMP end do
!$OMP end PARALLEL

! Repeat the last block of code for water tracers
if (l_wtrac) then
  ! Set up temporary field for z_uv(:,:,ntml(i,j)+1) which is required in
  !  call to calc_fqw_inv_wtrac
  allocate(z_uv_ntmlp1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))

!$OMP  PARALLEL  do SCHEDULE(STATIC) DEFAULT(SHARED) private(i,j,k)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      k = ntml(i,j) + 1
      z_uv_ntmlp1(i,j) = z_uv(i,j,ntml(i,j)+1)
    end do
  end do
!$OMP end PARALLEL do

  call calc_fqw_inv_wtrac(bl_levels, ntml, totqf_efl_meth1,                    &
                          totqf_efl_meth2, t_frac, zh, zh_frac,                &
                          zrzi, z_uv_ntmlp1, dzl,                              &
                          we_rho, w_ls, we_parm,                               &
                          dqw_sml_wtrac, fq_nt_zh_wtrac, moisten,              &
                          'SML', wtrac_bl)
  deallocate(z_uv_ntmlp1)

end if   !l_wtrac

!-----------------------------------------------------------------------
! Now decoupled layer
!-----------------------------------------------------------------------

!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(SHARED)                            &
!$OMP  private ( i, j, k, totqf_efl,  ml_tend, fa_tend, inv_tend)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    moisten(i,j) = .false.
    if (l_wtrac) totqf_efl_meth1(i,j) = 0
    if (l_wtrac) totqf_efl_meth2(i,j) = 0
    if ( t_frac_dsc(i,j)  >   zero ) then

      k = ntdsc(i,j)+1

        ! Calculate total (turb+micro) QW flux at subgrid inversion
      totqf_zhsc(i,j) = - we_rho_dsc(i,j)*dqw_dsc(i,j)                         &
                          + fq_nt_zhsc(i,j)
        ! Interpolate to entrainment flux-level
      totqf_efl = fq_nt_dscb(i,j) +                                            &
                ( totqf_zhsc(i,j) - fq_nt_dscb(i,j) )*zrzi_dsc(i,j)

      ml_tend = - ( totqf_zhsc(i,j)-fq_nt_dscb(i,j) )/dscdepth(i,j)
      fa_tend = zero
      if ( k+1  <=  bl_levels )                                                &
         fa_tend = - ( fq_nt(i,j,k+2)-fq_nt(i,j,k+1) )                         &
                     / dzl(i,j,k+1)
      inv_tend =       zhsc_frac(i,j) * ml_tend                                &
               + (one-zhsc_frac(i,j)) * fa_tend

      if (we_dsc_parm(i,j)+w_ls_dsc(i,j) >=  zero) then
          ! inversion moving up so inversion will moisten/dry
          ! depending on relative QW in level below
        moisten(i,j) = ( qw(i,j,k) <= qw(i,j,k-1) )
      else
          ! inversion moving down so inversion will moisten/dry
          ! depending on relative QW in level above
        moisten(i,j) = ( qw(i,j,k) <= qw(i,j,k+1) )
      end if

      if ( moisten(i,j) ) then

        if (l_wtrac .and. (totqf_efl < fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k)) )   &
           totqf_efl_meth1(i,j) = 1             ! Store method

        totqf_efl = max( totqf_efl,                                            &
                          fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
        if (we_dsc_parm(i,j)+w_ls_dsc(i,j)  >=  zero) then
            ! Ensure inversion level won't end up more moist than
            ! NTDSC by end of timestep.
          inv_tend = (qw(i,j,k-1)-qw(i,j,k))/timestep                          &
                     + fq_nt_dscb(i,j)/dscdepth(i,j)

          if (l_wtrac .and. totqf_efl >                                        &
                      (  (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                  &
                         /(one+ dzl(i,j,k)/dscdepth(i,j)) ) )                  &
             totqf_efl_meth2(i,j) = 1            ! Store method

          totqf_efl = min( totqf_efl,                                          &
                        (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                   &
                         /(one+ dzl(i,j,k)/dscdepth(i,j))   )
        end if
      else
        if (l_wtrac .and. totqf_efl > (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k)) )   &
           totqf_efl_meth1(i,j) = 1               ! Store method

        totqf_efl = min( totqf_efl,                                            &
                          fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
        if (we_dsc_parm(i,j)+w_ls_dsc(i,j)  >=  zero) then
            ! Ensure inversion level won't end up drier than
            ! NTDSC by end of timestep.
            ! Set INV_TEND to max allowable drying rate:
          inv_tend = (qw(i,j,k-1)-qw(i,j,k))/timestep                          &
                     + fq_nt_dscb(i,j)/dscdepth(i,j)

          if (l_wtrac .and. totqf_efl <                                        &
                      ( (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                   &
                         /(one+ dzl(i,j,k)/dscdepth(i,j)))   )                 &
            totqf_efl_meth2(i,j) = 1             ! Store method

          totqf_efl = max( totqf_efl,                                          &
                        (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))                   &
                         /(one+ dzl(i,j,k)/dscdepth(i,j))   )
        end if
      end if
      fqw(i,j,k) = t_frac_dsc(i,j) * ( totqf_efl - fq_nt(i,j,k) )

    end if
  end do
end do
!$OMP end PARALLEL do

! Repeat last block of code for water tracers
if (l_wtrac) then

  call calc_fqw_inv_wtrac(bl_levels, ntdsc, totqf_efl_meth1,                   &
                          totqf_efl_meth2, t_frac_dsc, dscdepth, zhsc_frac,    &
                          zrzi_dsc, dscdepth, dzl,                             &
                          we_rho_dsc, w_ls_dsc, we_dsc_parm,                   &
                          dqw_dsc_wtrac, fq_nt_zhsc_wtrac, moisten,            &
                          'DSC', wtrac_bl)
end if  ! l_wtrac

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP  private (i, j, k, kp, w_var_inv, weight, tke_nl_rh, delta_tke,          &
!$OMP  w_s_ent, w_s_cubed, w_m, wstar3, w_h)

!-----------------------------------------------------------------------
! Estimate turbulent w-variance scale at discontinuous inversions
!-----------------------------------------------------------------------

if (BL_diag%l_tke) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

      ! We expect the turbulent temperature perturbation scale close
      ! to the inversion, Tl' ~ ftl / w_scale, to be at most 1/2 dsl_inv,
      ! where ftl is the entrainment heat flux, and w_scale = sqrt(w_var)
      ! => w_scale >= 2 ftl / dsl_inv
      ! and then square that to get w_var.

      ! Want to set non-local TKE to max of existing value and this
      ! w_scale^2, but this is complicated by the fact that
      ! tke_nl is on theta-levels, while w_scale is on rho-levels.
      ! Compare w_scale^2 with the existing non-local TKE
      ! interpolated to the same rho-level; increase TKE on the neighbouring
      ! theta-levels if needed to increase interpolated rho-level value
      ! to the desired inversion value.

      ! Note: heights z_tq are on actual theta-levels (surface at k=0),
      ! but tke_nl is offset (surface at k=1)

      ! Do the above for discontinuous inversion at top of surface mixed-layer
      if ( sml_disc_inv(i,j) >= 1 ) then
        k = ntml(i,j) + 1
        kp = min( k+1, bl_levels )

        w_var_inv = 2.0_r_bl * (ftl(i,j,k)/rho_mix(i,j,k)) / dsl_sml(i,j)
        w_var_inv = w_var_inv * w_var_inv * rho_mix(i,j,k)

        weight = ( z_uv(i,j,k) - z_tq(i,j,k-1) )                               &
               / ( z_tq(i,j,k) - z_tq(i,j,k-1) )
        tke_nl_rh = (one-weight) * tke_nl(i,j,k)                               &
                  +      weight  * tke_nl(i,j,kp)

        delta_tke = w_var_inv - tke_nl_rh
        if ( delta_tke > zero ) then
          tke_nl(i,j,k)  = tke_nl(i,j,k)  + delta_tke
          tke_nl(i,j,kp) = tke_nl(i,j,kp) + delta_tke
        end if

      end if

      ! Repeat for discontinuous inversion at top of decoupled Sc layer
      if ( dsc_disc_inv(i,j) >= 1 ) then
        k = ntdsc(i,j) + 1
        kp = min( k+1, bl_levels )

        w_var_inv = 2.0_r_bl * (ftl(i,j,k)/rho_mix(i,j,k)) / dsl_dsc(i,j)
        w_var_inv = w_var_inv * w_var_inv * rho_mix(i,j,k)

        weight = ( z_uv(i,j,k) - z_tq(i,j,k-1) )                               &
               / ( z_tq(i,j,k) - z_tq(i,j,k-1) )
        tke_nl_rh = (one-weight) * tke_nl(i,j,k)                               &
                  +      weight  * tke_nl(i,j,kp)

        delta_tke = w_var_inv - tke_nl_rh

        if ( delta_tke > zero ) then
          tke_nl(i,j,k)  = tke_nl(i,j,k)  + delta_tke
          tke_nl(i,j,kp) = tke_nl(i,j,kp) + delta_tke
        end if

      end if

    end do
  end do
!$OMP end do NOWAIT
end if  ! (BL_diag%l_tke)

!-----------------------------------------------------------------------
! Calculate effective entrainment (ie. reduced to allow for subsidence
! increments in the ML) for use in tracer mixing.  Take theta_l as a
! representative scalar field since jump should always be the same sign
! and code therefore simpler.
! In this version the inversion fluxes are only implemented at one
! grid-level so only one element of these 3D arrays is used.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    we_lim(i,j,1)    = zero
    t_frac_tr(i,j,1) = zero
    zrzi_tr(i,j,1)   = zero
    we_lim_dsc(i,j,1)    = zero
    t_frac_dsc_tr(i,j,1) = zero
    zrzi_dsc_tr(i,j,1)   = zero
    we_lim(i,j,3)    = zero
    t_frac_tr(i,j,3) = zero
    zrzi_tr(i,j,3)   = zero
    we_lim_dsc(i,j,3)    = zero
    t_frac_dsc_tr(i,j,3) = zero
    zrzi_dsc_tr(i,j,3)   = zero
  end do ! i
end do ! j
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    kent(i,j)   = ntml(i,j)+1
    t_frac_tr(i,j,2) = t_frac(i,j)
    zrzi_tr(i,j,2) = zrzi(i,j)
    if ( t_frac(i,j)  >   zero ) then
      w_s_ent = zero
      k = ntml(i,j)
      if ( abs( dsl_sml(i,j) )  >=  rbl_eps ) w_s_ent =                        &
          min( zero, -sls_inc(i,j,k) * dzl(i,j,k) /dsl_sml(i,j) )
        ! Only allow w_e to be reduced to zero!
      we_lim(i,j,2) = rho_mix(i,j,k+1) *                                       &
                        max( zero, we_parm(i,j) + w_s_ent )
    else
      we_lim(i,j,2) = zero
    end if
    kent_dsc(i,j)   = ntdsc(i,j)+1
    t_frac_dsc_tr(i,j,2) = t_frac_dsc(i,j)
    zrzi_dsc_tr(i,j,2) = zrzi_dsc(i,j)
    if ( t_frac_dsc(i,j)  >   zero ) then
      w_s_ent = zero
      k = ntdsc(i,j)
      if ( abs( dsl_dsc(i,j) )  >= rbl_eps ) w_s_ent =                         &
          min( zero, -sls_inc(i,j,k) * dzl(i,j,k) /dsl_dsc(i,j) )
        ! Only allow w_e to be reduced to zero!
      we_lim_dsc(i,j,2) = rho_mix(i,j,k) *                                     &
                        max( zero, we_dsc_parm(i,j) + w_s_ent )
    else
      we_lim_dsc(i,j,2) = zero
    end if
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 12. Update standard deviations and gradient adjustment to use this
!     timestep's ZH (code from SF_EXCH)
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( unstable(i,j) ) then
      if (flux_grad  ==  Locketal2000) then
        w_s_cubed = 0.25_r_bl * zh(i,j) * fb_surf(i,j)
        if (w_s_cubed  >   zero) then
          w_m  =                                                               &
         ( w_s_cubed + v_s(i,j) * v_s(i,j) * v_s(i,j) ) ** one_third
          t1_sd(i,j) = 1.93_r_bl * ftl(i,j,1) / (rhostar_gb(i,j) * w_m)
          q1_sd(i,j) = 1.93_r_bl * fqw(i,j,1) / (rhostar_gb(i,j) * w_m)
          tv1_sd(i,j) = t(i,j,1) *                                             &
            ( one + c_virtual*q(i,j,1) - qcl(i,j,1) - qcf(i,j,1) ) *           &
            ( bt(i,j,1)*t1_sd(i,j) + bq(i,j,1)*q1_sd(i,j) )
          t1_sd(i,j) = max ( zero , t1_sd(i,j) )
          q1_sd(i,j) = max ( zero , q1_sd(i,j) )
          if (tv1_sd(i,j)  <=  zero) then
            tv1_sd(i,j) = zero
            t1_sd(i,j) = zero
            q1_sd(i,j) = zero
          end if
        end if
        grad_t_adj(i,j) = min( max_t_grad ,                                    &
                         a_grad_adj * t1_sd(i,j) / zh(i,j) )
        grad_q_adj(i,j) = zero
      else if (flux_grad  ==  HoltBov1993) then
          ! Use constants from Holtslag and Boville (1993)
          ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
          ! Neut limit GAMMA_TH = 7.2*wstar*FTL1/(ustar^2*zh)
        wstar3 = fb_surf(i,j) * zh(i,j)
        w_m =( v_s(i,j)**3 + 0.6_r_bl*wstar3 )**one_third

        grad_t_adj(i,j) = a_ga_hb93*(wstar3**one_third)*ftl(i,j,1)             &
                          / ( rhostar_gb(i,j)*w_m*w_m*zh(i,j) )
        ! GRAD_Q_ADJ(I,j) = A_GA_HB93*(WSTAR3**one_third)*FQW(I,j,1)
        !                  / ( RHOSTAR_GB(I,j)*W_M*W_M*ZH(I,j) )
        ! Set q term to zero for same empirical reasons as Lock et al
        grad_q_adj(i,j) = zero
      else if (flux_grad  ==  LockWhelan2006) then
          ! Use constants LockWhelan2006
          ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
          ! Neut limit GAMMA_TH = 7.5*FTL1/(ustar*zh)
        wstar3 = fb_surf(i,j) * zh(i,j)
        w_h =( ((4.0_r_bl/3.0_r_bl)*v_s(i,j))**3 + wstar3 )**one_third

        grad_t_adj(i,j) = a_ga_lw06 * ftl(i,j,1)                               &
                           / ( rhostar_gb(i,j)*w_h*zh(i,j) )
        grad_q_adj(i,j) = a_ga_lw06 * fqw(i,j,1)                               &
                           / ( rhostar_gb(i,j)*w_h*zh(i,j) )
      end if
    end if  ! test on UNSTABLE
  end do
end do
!$OMP end do NOWAIT

! (Note, water tracers assume flux_grad  =  Locketal2000 so no need to
! update wtrac_bl%grad_q_adj as it is always zero)

!-----------------------------------------------------------------------
!- Save diagnostics
!-----------------------------------------------------------------------

if (BL_diag%l_dzh) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ! fill unset values (rmdi<0) with zero
      ! and in cumulus set dzh=zh-z_lcl
      if (kprof_cu >= on .and. cumulus(i,j)) then
        BL_diag%dzh(i,j)= max( zero, zh(i,j)-z_lcl(i,j) )
      else
        BL_diag%dzh(i,j)= max( zero, dzh(i,j) )
      end if
    end do
  end do
!$OMP end do NOWAIT
end if
if (BL_diag%l_dscbase) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( dsc(i,j) ) then
        BL_diag%dscbase(i,j)= zhsc(i,j)-dscdepth(i,j)
      else
        BL_diag%dscbase(i,j)= rmdi
      end if
    end do
  end do
!$OMP end do NOWAIT
end if
if (BL_diag%l_cldbase) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( dsc(i,j) ) then
        BL_diag%cldbase(i,j)= zhsc(i,j)-zc_dsc(i,j)
      else
        BL_diag%cldbase(i,j)= zh(i,j)-zc(i,j)
      end if
    end do
  end do
!$OMP end do NOWAIT
end if
if (BL_diag%l_weparm_dsc) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( dsc(i,j) ) then
        BL_diag%weparm_dsc(i,j)= we_dsc_parm(i,j)
      else
        BL_diag%weparm_dsc(i,j)= we_parm(i,j)
      end if
    end do
  end do
!$OMP end do NOWAIT
end if
if (BL_diag%l_weparm) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      BL_diag%weparm(i,j)= we_parm(i,j)
    end do
  end do
!$OMP end do NOWAIT
end if

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------

!End of OpenMP parallel region
!$OMP end PARALLEL

! Deallocate water tracer working arrays

deallocate(fq_nt_zhsc_wtrac)
deallocate(fq_nt_zh_wtrac)
if (l_wtrac) then
  deallocate(qw_lapse_zero_dsc)
  deallocate(qw_lapse_zero_sml)
  deallocate(totqf_efl_meth2)
  deallocate(totqf_efl_meth1)
  deallocate(dqw_dsc_meth)
  deallocate(dqw_sml_meth)
  deallocate(dqw_sml_wtrac)
  deallocate(dqw_dsc_wtrac)
  deallocate(qls_inc_wtrac)
  deallocate(dfsubs_wtrac)
  deallocate(dfmic_wtrac)
  deallocate(fmic_wtrac)
  deallocate(fsubs_wtrac)
  deallocate(ntdsc_start)
  deallocate(ntml_start)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine kmkhz_9c
end module kmkhz_9c_mod
