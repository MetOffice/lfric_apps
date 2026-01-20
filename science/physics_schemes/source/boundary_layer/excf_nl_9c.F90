! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  subroutine EXCF_NL ----------------------------------------------

!  Purpose: To calculate non-local exchange coefficients,
!           entrainment parametrization and non-gradient flux terms.

!  Programming standard: UMDP3

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module excf_nl_9c_mod

use um_types, only: rbl_eps, r_bl, real_64

implicit none

character(len=*), parameter, private :: ModuleName = 'EXCF_NL_9C_MOD'

contains

subroutine excf_nl_9c (                                                        &
! in levels/switches
 bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                                       &
! in fields
 rdz,z_uv,z_tq,rho_mix,rho_wet_tq,rhostar_gb,v_s,fb_surf,zhpar,zh_prev,z_lcl,  &
 btm,bqm,btm_cld,bqm_cld,cf_sml,cf_dsc,                                        &
 bflux_surf,bflux_surf_sat,zeta_s,svl_diff_frac,                               &
 df_top_over_cp,zeta_r,bt_top,btt_top,btc_top,                                 &
 db_top,db_top_cld,chi_s_top,br_fback,                                         &
 df_dsct_over_cp,zeta_r_dsc,bt_dsct,btt_dsct,                                  &
 db_dsct,db_dsct_cld,chi_s_dsct,br_fback_dsc,                                  &
 db_ga_dry,db_noga_dry,db_ga_cld,db_noga_cld,                                  &
 dsl_sml,dqw_sml,dsl_dsc,dqw_dsc,ft_nt, fq_nt, sl, qw,                         &
! INOUT fields
 dsc_removed,dsc,cumulus,coupled,ntml,zh,zhsc,dscdepth,ntdsc,zc,zc_dsc,        &
 ft_nt_zh, ft_nt_zhsc, fq_nt_zh, fq_nt_zhsc, dzh,                              &
 fq_nt_zh_wtrac, fq_nt_zhsc_wtrac,                                             &
! out fields
 rhokm, rhokh, rhokm_top, rhokh_top,                                           &
 rhokh_top_ent, rhokh_dsct_ent, rhokh_surf_ent,                                &
 rhof2,rhofsc,f_ngstress,tke_nl,zdsc_base,nbdsc                                &
)
use atm_fields_bounds_mod, only: pdims, pdims_s, ScmRowLen, ScmRow
use bl_diags_mod, only: strnewbldiag
use tuning_segments_mod, only: bl_segment_size
use bl_option_mod, only:                                                       &
    off, on, dec_thres_cloud, dec_thres_cu, ng_stress,                         &
    BrownGrant97, BrownGrant97_limited, BrownGrant97_original,                 &
    flux_grad, Locketal2000, HoltBov1993, LockWhelan2006, entr_smooth_dec,     &
    entr_taper_zh, kprof_cu, klcl_entr, buoy_integ, buoy_integ_low,            &
    max_cu_depth, bl_res_inv, a_ent_shr_nml,  a_ent_2, one_third, two_thirds,  &
    l_reset_dec_thres,                                                         &
    l_use_var_fixes, dzrad_disc_opt, dzrad_ntm1, dzrad_1p5dz,                  &
    num_sweeps_bflux, l_converge_ga, l_use_sml_dsc_fixes, zero, one, one_half
use model_domain_mod, only: model_type, mt_single_column
use missing_data_mod, only: rmdi
use planet_constants_mod, only: vkman => vkman_bl, g => g_bl
use s_scmop_mod,   only: default_streams, t_avg, d_bl, scmdiag_bl
use stochastic_physics_run_mod, only: l_rp2 ,  i_rp_scheme, i_rp2b,            &
                                      a_ent_1_rp, g1_rp,                       &
                                      a_ent_shr_rp, a_ent_shr_rp_max,          &
                                      rp_idx, rp_max_idx, rp_min_idx
use excfnl_cci_mod, only: excfnl_cci
use excfnl_compin_mod, only: excfnl_compin
use free_tracers_inputs_mod, only: l_wtrac, n_wtrac

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim


implicit none

! in fields
integer, intent(in) ::                                                         &
 bl_levels
                 ! in maximum number of boundary layer levels

! Additional variables for SCM diagnostics which are dummy in full UM
integer, intent(in) ::                                                         &
  nSCMDpkgs             ! No of SCM diagnostics packages

logical, intent(in) ::                                                         &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

!     Declaration of new BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

real(kind=r_bl), intent(in) ::                                                 &
 rdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                      &
                bl_levels),                                                    &
                              ! in Reciprocal of distance between
                              !    T,q-levels (m^-1). 1/RDZ(,K) is
                              !    the vertical distance from level
                              !    K-1 to level K, except that for
                              !    K=1 it is just the height of the
                              !    lowest atmospheric level.
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                     &
                bl_levels),                                                    &
                              ! in For a vertically staggered grid
                              !    with a u,v-level first above the
                              !    surface, Z_UV(*,K) is the height
                              !    of the k-th u,v-level (half level
                              !    k-1/2) above the surface
 z_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                     &
                bl_levels),                                                    &
                              ! in For a vertically staggered grid
                              !    with a u,v-level first above the
                              !    surface, Z_TQ(*,K) is the height
                              !    of the k-th T,q-level (full level
                              !    k) above the surface
 rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
                bl_levels),                                                    &
                              ! in For a vertically staggered grid
                              !    with a u,v-level first above the
                              !    surface, RHO_MIX(*,K) is the
                              !    density at the k-th u,v-level
                              !    above the surface
 rho_wet_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,               &
                bl_levels),                                                    &
                              ! in For a vertically staggered grid
                              !    with a u,v-level first above the
                              !    surface, RHO_WET_TQ(*,K) is the
                              !    density of the k-th T,q-level
                              !    above the surface
 bqm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                      &
                bl_levels),                                                    &
                              ! in Buoyancy parameters for clear and
 btm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                      &
                bl_levels),                                                    &
                              !    cloudy air on half levels
 bqm_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
                bl_levels),                                                    &
                              !    (*,K) elements are k+1/2 values
 btm_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
                bl_levels)                  !

real(kind=r_bl), intent(in) ::                                                 &
 rhostar_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                              ! in Surface density (kg/m3)
 v_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                     &
                              ! in Surface friction velocity (m/s).
 fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in Buoyancy flux at the surface over
                              !    density (m^2/s^3).
 zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                              ! in Height of top of NTPAR
                              !    NOTE: CAN BE ABOVE BL_LEVELS-1
 zh_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in previous timestep surface mixed-layer height
 z_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                              ! in Height of lifting condensation
                              !    level.
 bflux_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                              ! in Surface buoyancy flux (kg/m/s^3).
 bflux_surf_sat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! in Saturated-air surface buoyancy
                              !    flux.
 db_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                              ! in Buoyancy jump across the top of
                              !    the SML (m/s^2).
 df_top_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! in Radiative flux change at cloud top
                              !    divided by c_P (K.kg/m^2/s).
 bt_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                              ! in Buoyancy parameter at the top of
                              !    the b.l. (m/s^2/K).
 btt_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in In-cloud buoyancy parameter at
                              !    the top of the b.l. (m/s^2/K).
 btc_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in Cloud fraction weighted buoyancy
                              !    parameter at the top of the b.l.
 db_top_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                              ! in In-cloud buoyancy jump at the
                              !    top of the b.l. (m/s^2).
 chi_s_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                              ! in Mixing fraction of just saturated
                              !    mixture at top of the b.l.
 zeta_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                              ! in Non-cloudy fraction of mixing
                              !    layer for surface forced
                              !    entrainment term.
 zeta_r(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                              ! in Non-cloudy fraction of mixing
                              !    layer for cloud top radiative
                              !    cooling entrainment term.
 db_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in Buoyancy jump across the top of
                              !    the DSC layer (m/s^2).
 df_dsct_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                              ! in Radiative flux change at DSC top
                              !    divided by c_P (K.kg/m^2/s).
 bt_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in Buoyancy parameters at the top of
 btt_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                              !    the DSC layer (m/s^2/K)
 db_dsct_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                              ! in In-cloud buoyancy jump at the
                              !    top of the DSC (m/s^2).
 chi_s_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                              ! in Mixing fraction of just saturated
                              !    mixture at top of the DSC
 zeta_r_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                              ! in Non-cloudy fraction of DSC
                              !    for cloud top radiative
                              !    cooling entrainment term.
 br_fback(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
 br_fback_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                              ! in Weights for degree of buoyancy
                              !    reversal feedback
 dqw_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in QW change across SML disc inv
 dsl_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in SL change across SML disc inv
 dqw_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in QW change across DSC disc inv
 dsl_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! in SL change across DSC disc inv
 cf_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                              ! in cloud fraction of SML
 cf_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! in cloud fraction of DSC layer

real(kind=r_bl) , intent(in out) ::                                            &
 svl_diff_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! in Fractional svl decoupling difference

real(kind=r_bl) , intent(in) ::                                                &
 db_ga_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
                2:bl_levels),                                                  &
                              ! in Cloudy and cloud-free buoyancy
 db_noga_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
                2:bl_levels),                                                  &
                              !    jumps for flux integral
 db_ga_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
                2:bl_levels),                                                  &
                              !    calculation (m/s2):
 db_noga_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
                2:bl_levels)

real(kind=r_bl), intent(in) ::                                                 &
 ft_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
                bl_levels+1),                                                  &
                              ! in Non-turbulent heat (rho*Km/s) and
 fq_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
                bl_levels+1),                                                  &
                              !      moisture (rho*m/s) fluxes
 sl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),            &
                              ! TL + G*Z/CP (K)
 qw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels)
                              ! Total water content (kg per kg air).

! INOUT fields
logical, intent(in out) ::                                                     &
 coupled(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! INOUT Flag to indicate Sc layer
                              !       weakly decoupled
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! INOUT Flag for cumulus
 dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! INOUT Flag set if decoupled stratocu
                              !       layer found.

integer, intent(in out) ::                                                     &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                              ! INOUT  Number of turbulently mixed
                              !        layers.
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                              ! INOUT  Top level of any decoupled
                              !        turbulently mixed Sc layer
 dsc_removed(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! INOUT  Flag to indicate why dsc layer removed

real(kind=r_bl), intent(in out) ::                                             &
 zc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                      &
                              ! INOUT Cloud depth (not cloud fraction
                              !    weighted) (m).
 zc_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                              ! INOUT Cloud depth (not cloud fraction
                              !    weighted) (m).
 zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                              ! INOUT Cloud-layer height (m)
 zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                      &
                              ! INOUT Boundary layer height (m)
 dzh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                     &
                              ! INOUT Inversion thickness (m)
 dscdepth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                              ! INOUT Decoupled cloud-layer depth (m)
 ft_nt_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                              ! INOUT Non-turbulent heat (rho*Km/s)
 ft_nt_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                              !        and moisture (rho*m/s) fluxes
 fq_nt_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                              !        evaluated at the SML and DSC
 fq_nt_zhsc (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                              !        inversions
 fq_nt_zh_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,n_wtrac),  &
                              !        Water tracer fq_nt_zh
 fq_nt_zhsc_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,n_wtrac)
                              !        Water tracer fq_nt_zhsc

! out fields
integer, intent(out) ::                                                        &
 nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
       ! out Bottom level of any decoupled turbulently mixed Sc layer.

real(kind=r_bl), intent(out) ::                                                &
 rhokm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
                2:bl_levels),                                                  &
                              ! out Layer k-1 - to - layer k
                              !     turbulent mixing coefficient
                              !     for momentum (kg/m/s).
 rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
                2:bl_levels),                                                  &
                              ! out Layer k-1 - to - layer k
                              !     turbulent mixing coefficient
                              !     for heat and moisture (kg/m/s).
 rhokm_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
                2:bl_levels),                                                  &
                              ! out exchange coefficient for
                              !     momentum due to top-down mixing
 rhokh_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
                2:bl_levels),                                                  &
                              ! out exchange coefficient for
                              !     heat and moisture due to top-down
                              !     mixing
 rhof2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
                2:bl_levels),                                                  &
 rhofsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
                2:bl_levels),                                                  &
                              ! out f2 and fsc term shape profiles
                              !     multiplied by rho
 tke_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,2:bl_levels)
                              ! out Non-local TKE diag, currently times rho

real(kind=r_bl), intent(out) ::                                                &
 f_ngstress(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,       &
            2:bl_levels)      ! out dimensionless function for
                              !     non-gradient stresses

real(kind=r_bl), intent(out) ::                                                &
 zdsc_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                              ! out Height of base of K_top in DSC
 rhokh_surf_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! out SML surface-driven entrainment KH
 rhokh_top_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! out SML top-driven entrainment KH
 rhokh_dsct_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! out DSC top-driven entrainment KH
!  ---------------------------------------------------------------------
!    Local and other symbolic constants :-

character(len=*), parameter ::  RoutineName = 'EXCF_NL_9C'

real(kind=r_bl) :: a_ent_1,c_t,a_ent_shr,dec_thres_clear
real(kind=r_bl) :: g1
integer :: n_steps
parameter (                                                                    &
 c_t=one,                                                                      &
                              ! Parameter in Zilitinkevich term.
 dec_thres_clear=one,                                                          &
                              ! Decoupling threshold for cloud-free
                              ! boundary layers (larger makes
                              ! decoupling less likely)
 n_steps=3                                                                     &
                              ! Number of steps through the mixed
                              ! layer per sweep
)

real(kind=r_bl) :: s_m,a_ngs
parameter (                                                                    &
 s_m   = one,                                                                  &
                              ! empirical parameters in
 a_ngs = 2.7_r_bl                                                              &
                              ! non-gradient stresses
)

!  Define local storage.
!  (a) Workspace.

integer ::                                                                     &
 ksurf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
           ! First Theta-level above surface layer well-mixed SC layer

integer, allocatable :: iset_wtrac(:,:)    ! Used to set water tracer fields

logical ::                                                                     &
 ng_stress_calculate           ! Flag to do ng stress calculation

logical ::                                                                     &
 scbase(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                              ! Flag to signal base of CML reached
 test_well_mixed(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                              ! Flag to test wb integration
                              ! for a well-mixed layer
 ksurf_iterate(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! Flag to perform iteration to
                              ! find top of Ksurf
 ktop_iterate(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! Flag to perform iteration to
                              ! find base of Ktop
real(kind=r_bl) ::                                                             &
 kh_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,2:bl_levels)
                              ! Shape factor for non-local
                              ! turbulent mixing coefficient

real(kind=r_bl) ::                                                             &
        wbmix( pdims%i_start:pdims%i_start+scmrowlen-1                         &
             , pdims%j_start:pdims%j_start+scmrow-1                            &
             , pdims%k_start:pdims%k_end )
                              ! WB*DZ if were diag as mixed
real(kind=r_bl) ::                                                             &
        wbend( pdims%i_start:pdims%i_start+scmrowlen-1                         &
             , pdims%j_start:pdims%j_start+scmrow-1                            &
             , pdims%k_start:pdims%k_end )
                              ! WB*DZ after dec diag
real(kind=r_bl) ::                                                             &
        wbend_sml( pdims%i_start:pdims%i_start+scmrowlen-1                     &
             , pdims%j_start:pdims%j_start+scmrow-1                            &
             , pdims%k_start:pdims%k_end )
                              ! WB*DZ for sml if ksurf_iterate

real(kind=r_bl)::                                                              &
       w_m_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                    ! Turbulent velocity scale for momentum
                    !   evaluated at the top of the b.l.
       w_h_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                    ! Turbulent velocity scale for scalars
                    !   evaluated at the top of the b.l.
       prandtl_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                    ! Turbulent Prandtl number
                    !   evaluated at the top of the b.l.
       rhokm_surf_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                    ! SML surface-driven entrainment Km
       rhokm_top_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                    ! SML top-driven entrainment Km
       rhokm_dsct_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                    ! DSC top-driven entrainment Km
       kh_top_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                    ! Factors to ensure K_H and K_M profiles are
       km_top_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                    !   continuous at ZH and ZHSC
       kh_sct_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                    !               "
       km_sct_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                    !               "
       kh_dsct_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                    !               "
       km_dsct_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                    !               "
       v_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                    ! Velocity scale for top-down convection
       v_top_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &

       v_sum(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                    ! total velocity scale
       v_sum_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &

       zsml_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                    ! Height of top of surf-driven K in SML
       zsml_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                    ! Height of base of top-driven K in SML
       rhokh_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                    ! rho*KH at the LCL in cumulus
       cu_depth_scale(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                    ! Depth scale for decay of K profile
                    ! above LCL in cumulus (in m)
       scdepth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                    ! Depth of top-driven mixing in SML
       z_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                    ! inversion height (top of K profiles)
       wb_surf_int(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                    ! Estimate of wb integrated over surface layer
       wb_dzrad_int(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                    ! Estimate of wb integrated over cloud-top region
       dzrad(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                    ! Depth of cloud-top (radiatively cooled) region
       v_ktop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                    ! velocity scale for K_top profile
       v_ksum(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                    ! total velocity scale
       z_cbase(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                    ! cloud base height
       wb_ratio(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                    ! WBN_INT/WBP_INT
       dec_thres(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                    ! Local decoupling threshold
       wbp_int(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                    ! Positive part of buoyancy flux integral
       wbn_int(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                    ! Negative part of buoyancy flux integral
       zinv_pr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                    ! Height of layer top above surface
       khtop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                    ! temporary KH_top in wb integration
       khsurf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                    ! temporary KH_surf in wb integration
       zwb0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                    ! height at which wb assumed to go to zero
       z_top_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                    ! upper height limit on K profile
       z_bot_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                    ! lower height limit on K profile
       z_inc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                    ! Step size (m)
       rho_we(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                    ! rho*param.d entrainment rates...
       rho_we_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                    !  ...for surf and DSC layers (kg/m2/s)
       rho_we_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &

       cf_ml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                    ! Mixed layer cloud fraction
       df_ctop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                    ! Cloud-top radiative flux divergence
       dqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                    ! QW jump across inversion
       dsl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                    ! SL jump across inversion
       ft_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                    ! Non-turbulent heat flux at DSC base
       fq_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                    ! Non-turbulent moisture flux at DSC base
       tothf_zi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                    ! Total heat and moisture fluxes at
       totqf_zi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                    !   the inversion height, Zi

! Double precision friction velocity to avoid calcualtions going out of
! Single precision range
real(kind=real_64) ::                                                          &
     v_s_dbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Gradient-adjusted buoyancy intervals corrected for change in zh this timestep
real(kind=r_bl) ::  db_ga_dry_n                                                &
                           ( pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,                        &
                             2:bl_levels )
real(kind=r_bl) ::  db_ga_cld_n                                                &
                           ( pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end,                        &
                             2:bl_levels )
! Correction factor used to calculate the above
real(kind=r_bl) ::  ga_fac                                                     &
                           ( pdims%i_start:pdims%i_end,                        &
                             pdims%j_start:pdims%j_end )

!  (b) Scalars.

real(kind=r_bl) ::                                                             &
 Prandtl,                                                                      &
                  ! Turbulent Prandtl number.
 pr_neut,                                                                      &
                  ! Neutral limit for Prandtl number
 pr_conv,                                                                      &
                  ! Convective limit for Prandtl number
 zk_uv,                                                                        &
                  ! Height above surface of u,v-level.
 zk_tq,                                                                        &
                  ! Height above surface of T,q-level.
 wstar3,                                                                       &
                  ! Cube of free-convective velocity scale
 c_ws,                                                                         &
                  ! Empirical constant multiplying Wstar
 c_tke,                                                                        &
                  ! Empirical constant in tke diagnostic
 w_s_cubed_uv,                                                                 &
                  ! WSTAR for u,v-level
 w_s_cubed_tq,                                                                 &
                  !   and T,q-level
 w_m_uv,                                                                       &
                  ! Turbulent velocity scale for momentum: u,v-level
 w_m_tq,                                                                       &
                  !   and T,q-level
 w_h_uv,                                                                       &
                  ! Turbulent velocity scale for scalars: u,v-level
 w_h_tq,                                                                       &
                  !   and T,q-level
 w_m_hb_3,                                                                     &
                  ! Cube of W_M, as given by Holtslag and Boville, 93
 w_m_neut,                                                                     &
                  ! Neutral limit for W_M
 sf_term,                                                                      &
                  ! Surface flux term for entrainment parametrization.
 sf_shear_term,                                                                &
                  ! Surface shear term for entrainment paramn.
 ir_term,                                                                      &
                  ! Indirect radiative term for entrainment paramn.
 dr_term,                                                                      &
                  ! Direct radiative term for entrainment paramn.
 evap_term,                                                                    &
                  ! Evaporative term in entrainment parametrization.
 zil_corr,                                                                     &
                  ! Zilitinkevich correction term in entrn. paramn.
 zeta_s_fac,                                                                   &
                  ! Factor involving ZETA_S.
 zeta_r_sq,                                                                    &
                  ! ZETA_R squared.
 zr,                                                                           &
                  ! Ratio ZC/ZH.
 z_pr,                                                                         &
                  ! Height above surface layer
 zh_pr,                                                                        &
                  ! Height of layer top above surface
 z_ratio,                                                                      &
                  ! Ratio of heights
 rhokh_ent,                                                                    &
                  ! entrainment eddy viscosity
 rhokm_dsct,                                                                   &
                  ! top-down Km in decoupled layer
 frac_top,                                                                     &
                  ! Fraction of turbulent mixing driven from the top
 factor,                                                                       &
                  ! Temporary scalar
 alpha_t,                                                                      &
                  ! Parametrized fraction of cloud-top
                  ! radiative cooling within the inversion
 dz_inv,                                                                       &
                  ! Parametrizzed inversion thickness (m)
 l_rad,                                                                        &
                  ! Estimate of e-folding radiative flux
                  ! decay depth (assumed >= 25m)
 wb_cld,                                                                       &
                   ! Cloud layer buoyancy flux
 wb_scld,                                                                      &
                   ! Sub-cloud layer buoyancy flux
 cld_frac,                                                                     &
                   ! Vertical fraction of layer containing cloud
 zb_ktop,                                                                      &
                   ! height of base of K_top profile
 db_ratio,                                                                     &
                   ! Temporary in ZWB0 calculation
 gamma_wbs,                                                                    &
                   ! Surface layer wb gradient
 wsl_dzrad_int,                                                                &
                   ! Estimate of wsl and wqw integrated over
 wqw_dzrad_int,                                                                &
                   !   the cloud-top region
 wslng,                                                                        &
                   ! Non-gradient part of SL flux
 wqwng,                                                                        &
                   ! Non-gradient part of QW flux
 f2, fsc,                                                                      &
                   ! Shape functions for non-gradient fluxes
 lcl_fac,                                                                      &
                   ! fraction of LCL height to use as minimum mixing depth
 interp
                   ! Weight for vertical interpolation between model-levels

integer ::                                                                     &
 i,j,                                                                          &
                  ! Loop counter (horizontal field index).
 k,                                                                            &
                  ! Loop counter (vertical level index).
 n_sweep,                                                                      &
                  ! sweep counter
 ns,                                                                           &
                  ! step counter
 i_wt,                                                                         &
                  ! water tracer counter
 k_inv
                  ! Theta-level below the inversion

! 2D arrays for optimisation

integer :: ntop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                   ! top level of surf-driven K profile
           ntml_new(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                   ! temporary in NTML calculation
           kwb0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! level at which wb assumed to go to zero

integer :: up(pdims%i_end*pdims%j_end)
                   ! indicator of upward/downward sweep

logical :: status_ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Array introduced to calculate kwb0
logical :: kstatus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Flag for whether to apply surface-driven entrainment flux
logical :: l_apply_surf_ent

! Variables for vector compression

integer :: ic
integer :: c_len
logical :: to_do(pdims%i_end*pdims%j_end)
integer :: ind_todo(pdims%i_end*pdims%j_end)
integer :: c_len_i
logical :: todo_inner(pdims%i_end*pdims%j_end)
integer :: ind_todo_i(pdims%i_end*pdims%j_end)

integer :: i1, j1, l

integer            :: jj          !Cache blocking - loop index

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up water tracer field
if (l_wtrac) then
  allocate(iset_wtrac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
end if

!-----------------------------------------------------------------------
!     Set values of A_ENT_1, G1 and A_ENT_SHR

if (l_rp2) then
  a_ent_1 = a_ent_1_rp(rp_idx)       ! Entrainment parameter
  g1 = g1_rp(rp_idx)                 ! Velocity scale parameter
else
  a_ent_1 = 0.23_r_bl             ! Entrainment parameter
  g1 = 0.85_r_bl                  ! Velocity scale parameter
end if

if ( l_rp2 .and. i_rp_scheme == i_rp2b ) then
  !
  ! Alternative calculation of a_ent_shr to extend the range
  ! over which a_ent_shr varies with a_ent_1_rp.
  !
  if ( a_ent_1_rp(rp_min_idx) < 0.23_r_bl                                      &
         .and. a_ent_1_rp(rp_max_idx) > 0.23_r_bl ) then
    a_ent_shr =                                                                &
     ( ( a_ent_1 - a_ent_1_rp(rp_min_idx) ) /                                  &
       ( a_ent_1_rp(rp_max_idx)-0.23_r_bl ) ) *                                &
     ( ( a_ent_shr_rp_max*                                                     &
         ( a_ent_1 - 0.23_r_bl ) /                                             &
         ( a_ent_1_rp(rp_max_idx) - a_ent_1_rp(rp_min_idx) ) )                 &
     - ( a_ent_shr_rp*                                                         &
         ( a_ent_1 - a_ent_1_rp(rp_max_idx) ) /                                &
         ( 0.23_r_bl - a_ent_1_rp(rp_min_idx) ) ) )
  else
    ! Case a_ent_1_rp(rp_min_idx) = a_ent_1 = a_ent_1_rp(rp_max_idx)
    a_ent_shr = a_ent_shr_rp * a_ent_1 / 0.23_r_bl
  end if
else
  a_ent_shr = a_ent_shr_nml * a_ent_1 / 0.23_r_bl   ! Entrainment parameter.
end if

!-----------------------------------------------------------------------
! Index to subroutine EXCFNL9C

! 0. Calculate top-of-b.l. velocity scales and Prandtl number.
! 1. Calculate the top-of-b.l. entrainment parametrization
! 2. Estimate the depths of top-down and surface-up mixing.
!   2.1 First test for well-mixed boundary layer
!   2.2 Iterate to find top of surface-driven mixing, ZSML_TOP,
!   2.3 Iterate to find the base of the top-driven K profile, ZDSC_BASE.
! 3. Calculate factors required to ensure that the K profiles are
!    continuous at the inversion
! 4. Calculate height dependent turbulent transport coefficients
!    within the mixing layers.

!-----------------------------------------------------------------------
! 0.  Calculate top-of-b.l. velocity scales and Prandtl number.
!-----------------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP  private(i, j, k, jj, i_wt, c_ws, wstar3, pr_neut, pr_conv, w_m_neut,    &
!$OMP  zeta_s_fac, sf_term, sf_shear_term, zeta_r_sq, ir_term, zr,             &
!$OMP  evap_term, dz_inv, l_rad, alpha_t, dr_term, zil_corr,                   &
!$OMP  rhokh_ent, frac_top, zh_pr, factor, rhokm_dsct,                         &
!$OMP  wsl_dzrad_int, wqw_dzrad_int, db_ratio, zb_ktop, f2, fsc,               &
!$OMP  z_ratio, z_pr, wslng, wqwng, wb_scld, wb_cld, cld_frac, l,              &
!$OMP  j1, i1, ic,  w_m_hb_3, zk_uv, zk_tq, Prandtl, w_h_uv, w_h_tq,           &
!$OMP  w_m_uv, w_m_tq, w_s_cubed_tq, w_s_cubed_uv, gamma_wbs, c_tke,           &
!$OMP  n_sweep, ns, lcl_fac, ng_stress_calculate, l_apply_surf_ent, interp,    &
!$OMP  k_inv)

!cdir collapse
!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    rhokh_surf_ent(i,j) = zero
    rhokh_top_ent(i,j)  = zero
    rhokh_dsct_ent(i,j) = zero
    rhokm_surf_ent(i,j) = zero
    rhokm_top_ent(i,j)  = zero
    rhokm_dsct_ent(i,j) = zero
    v_top(i,j) = zero
    v_sum(i,j) = zero
    v_top_dsc(i,j) = zero
    v_sum_dsc(i,j) = zero
    rho_we(i,j)     = zero
    rho_we_sml(i,j) = zero
    rho_we_dsc(i,j) = zero

    if (fb_surf(i,j)  >=  zero) then

      ! The calculations below involve v_s**4. In very rare circumstances
      ! v_s can be order(10^-11), therefore v_s**4 is order(10^-44)
      ! which is outside the range of single precision calculations
      ! Hence we create a double precision version to enforce the correct
      ! calculation
      v_s_dbl(i,j) = v_s(i,j)

        ! Free-convective velocity scale cubed

      if (coupled(i,j)) then
        if ( entr_smooth_dec == entr_taper_zh ) then
          wstar3 = (      svl_diff_frac(i,j)  * zhsc(i,j)                      &
                   + (one-svl_diff_frac(i,j)) * zh(i,j)   )* fb_surf(i,j)
        else
          wstar3 = zhsc(i,j) * fb_surf(i,j)
        end if
      else
        wstar3 =   zh(i,j) * fb_surf(i,j)
      end if

      if (flux_grad  ==  Locketal2000) then

          ! Turbulent velocity scale for momentum

        c_ws = 0.25_r_bl
        w_m_top(i,j) = (v_s(i,j)*v_s(i,j)*v_s(i,j) +                           &
                        c_ws*wstar3)**one_third

          ! Turbulent Prandtl number and velocity scale for scalars
          ! gives 0.375<Pr<0.75 for convective to neutral conditions
        prandtl_top(i,j) = 0.75_r_bl *                                         &
                     ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +   &
                         (one/25.0_r_bl)*wstar3*w_m_top(i,j) ) /               &
                     ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +   &
                         (2.0_r_bl/25.0_r_bl)*wstar3*w_m_top(i,j) )
        w_h_top(i,j) = w_m_top(i,j) / prandtl_top(i,j)
      else if (flux_grad  ==  HoltBov1993) then
        c_ws = 0.6_r_bl
        w_m_top(i,j) = (v_s(i,j)*v_s(i,j)*v_s(i,j) +                           &
                        c_ws*wstar3)**one_third
          ! Using Lock et al interpolation but with
          ! HB93 range of 0.6<Pr<1
        pr_neut = one
        pr_conv = 0.6_r_bl
        prandtl_top(i,j) = pr_neut *                                           &
              ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +          &
                                (one/25.0_r_bl)*wstar3*w_m_top(i,j) ) /        &
              ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +          &
                  (pr_neut/(25.0_r_bl*pr_conv))*wstar3*w_m_top(i,j) )
        w_h_top(i,j) = w_m_top(i,j) / prandtl_top(i,j)
      else if (flux_grad  ==  LockWhelan2006) then
        pr_neut = 0.75_r_bl
        pr_conv = 0.6_r_bl
        c_ws    = 0.42_r_bl   !  ~ 0.75^3
          ! Slightly contrived notation since really we know W_H_TOP
          ! and Prandtl range but this makes similarity to
          ! Lock et al clearer (possibly!)
        w_m_neut = ( v_s(i,j)*v_s(i,j)*v_s(i,j) +                              &
                     c_ws*wstar3 )**one_third
        w_h_top(i,j) = w_m_neut / pr_neut

        prandtl_top(i,j) = pr_neut *                                           &
              ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +          &
                  (one/25.0_r_bl)*wstar3*w_m_neut ) /                          &
              ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +          &
                  (pr_neut/(25.0_r_bl*pr_conv))*wstar3*w_m_neut )

        w_m_top(i,j) = w_h_top(i,j) * prandtl_top(i,j)
      end if

    else
      w_m_top(i,j) = v_s(i,j)
      prandtl_top(i,j) = 0.75_r_bl
      w_h_top(i,j) = w_m_top(i,j) / prandtl_top(i,j)
    end if
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
! 1. Calculate the top-of-b.l. entrainment parametrization
!-----------------------------------------------------------------------
! 1.1 Initialise 3D arrays
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do k = 2, bl_levels
  !cdir collapse
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      rhokh(i,j,k) = zero
      rhokm(i,j,k) = zero
      rhokh_top(i,j,k) = zero
      rhokm_top(i,j,k) = zero
      rhof2(i,j,k)  = zero
      rhofsc(i,j,k) = zero
      f_ngstress(i,j,k) = zero
      kh_surf(i,j,k) = zero
      tke_nl(i,j,k) = zero
    end do
  end do
end do
!$OMP end do

!cdir collapse
!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    !-----------------------------------------------------------------------
    ! 1.2 Calculate top-of-b.l. entrainment mixing coefficients
    !     and store b.l. top quantities for later use.
    !-----------------------------------------------------------------------
    !      FIRST the top of the SML (if not coupled)
    !-----------------------------------------------
    k = ntml(i,j)+1
    if ( .not. coupled(i,j) .and. fb_surf(i,j)  >=  zero ) then
          !-----------------------------------------------------------
          ! Calculate the surface buoyancy flux term
          !-----------------------------------------------------------
      zeta_s_fac = (one - zeta_s(i,j)) * (one - zeta_s(i,j))
      sf_term = a_ent_1 * max ( zero ,                                         &
                          ( (one - zeta_s_fac) * bflux_surf(i,j)               &
                            + zeta_s_fac * bflux_surf_sat(i,j) ) )
          !-----------------------------------------------------------
          ! Calculate the surface shear term
          !-----------------------------------------------------------
      sf_shear_term =  a_ent_shr * v_s(i,j) * v_s(i,j) * v_s(i,j)              &
                      * rho_mix(i,j,k)  / zh(i,j)
          !-----------------------------------------------------------
          ! Calculate the indirect radiative term
          !-----------------------------------------------------------
      zeta_r_sq = zeta_r(i,j)*zeta_r(i,j)
      ir_term = ( bt_top(i,j)*zeta_r_sq +                                      &
                  btt_top(i,j)*(one-zeta_r_sq) )                               &
                * a_ent_1 * df_top_over_cp(i,j)
          !-----------------------------------------------------------
          ! Calculate the evaporative term
          !-----------------------------------------------------------
      if ( db_top(i,j)  >   zero) then
        zr = sqrt( zc(i,j) / zh(i,j) )
        evap_term = a_ent_2 * rho_mix(i,j,k)                                   &
                  * chi_s_top(i,j) * chi_s_top(i,j)                            &
                  * zr * zr * zr * db_top_cld(i,j)                             &
                  * sqrt( zh(i,j) * db_top(i,j) )
      else
        evap_term = zero
      end if
          !-----------------------------------------------------------
          ! Combine forcing terms to calculate the representative
          ! velocity scales
          !-----------------------------------------------------------
      v_sum(i,j) = ( (sf_term + sf_shear_term +                                &
                      ir_term + evap_term)                                     &
                   * zh(i,j) /(a_ent_1*rho_mix(i,j,k)) )**one_third
      v_top(i,j) = ( (ir_term+evap_term) * zh(i,j)                             &
                           / (a_ent_1*rho_mix(i,j,k)) )**one_third
          !-----------------------------------------------------------
          ! Calculate the direct radiative term
          !  can only calculate for DB_TOP > 0
          !-----------------------------------------------------------
      if ( db_top(i,j)  >   zero) then
        dz_inv  = min( v_sum(i,j)*v_sum(i,j) / db_top(i,j) ,100.0_r_bl )
        l_rad   = 15.0_r_bl * max( one, 200.0_r_bl/(zc(i,j)+rbl_eps) )
        alpha_t = one - exp(-one_half*dz_inv/l_rad)
           ! Make enhancement due to buoyancy reversal feedback
        alpha_t = alpha_t + br_fback(i,j)*(one-alpha_t)
        dr_term = btc_top(i,j) * alpha_t * df_top_over_cp(i,j)
           !----------------------------------------------------------
           ! Combine terms to calculate the entrainment
           ! mixing coefficients
           !----------------------------------------------------------
        zil_corr = c_t * ( (sf_term + sf_shear_term +                          &
                            ir_term + evap_term) /                             &
                    (rho_mix(i,j,k) * sqrt(zh(i,j))) )**two_thirds

        rho_we_sml(i,j) = (sf_term + sf_shear_term                             &
                    +      ir_term + evap_term + dr_term)                      &
                       / ( db_top(i,j) + zil_corr )

        rhokh_ent = rho_we_sml(i,j)/ rdz(i,j,k)

        frac_top = v_top(i,j) / ( v_top(i,j)+w_h_top(i,j)+rbl_eps )
        if ( l_use_sml_dsc_fixes ) then
          ! If bug-fix is on, leave surface-driven entrainment flux set to
          ! zero in cumulus layers if we are using the buoyancy-flux
          ! integration method to calculate the non-local mixing profile
          ! (code originally converged on a consistent negative buoyancy
          !  flux at the top and then added the entrainment flux onto
          !  this, which is double-counting).
          l_apply_surf_ent = ( .not. ( ( kprof_cu==buoy_integ .or.             &
                                         kprof_cu==buoy_integ_low ) .and.      &
                                       cumulus(i,j) ) )
        else
          ! If bug-fix is off, always add on surface-driven entrainment flux
          l_apply_surf_ent = .true.
        end if
        if ( l_apply_surf_ent ) then
          ! max to avoid rounding errors giving small negative numbers
          rhokh_surf_ent(i,j) = max(rhokh_ent * ( one - frac_top ), zero)
          rhokm_surf_ent(i,j) = prandtl_top(i,j) * rhokh_surf_ent(i,j)         &
                              * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))       &
                              * rho_wet_tq(i,j,k-1) / rho_mix(i,j,k)
        end if
        rhokh_top_ent(i,j) = rhokh_ent * frac_top
        rhokm_top_ent(i,j) = 0.75_r_bl * rhokh_top_ent(i,j)                    &
                             * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))        &
                             * rho_wet_tq(i,j,k-1) / rho_mix(i,j,k)
        if (.not. l_use_var_fixes) then
          rhokm(i,j,k) = rhokm_surf_ent(i,j)
          rhokm_top(i,j,k) = rhokm_top_ent(i,j)
        end if

      end if    ! test on DB_TOP gt 0
    end if
    !----------------------------------------------------------------
    !      then the top of the DSC (if coupled use ZHSC length-scale)
    !----------------------------------------------------------------
    if ( ntdsc(i,j)  >   0 ) then
      k = ntdsc(i,j)+1
      if (coupled(i,j)) then
            !--------------------------------------------------------
            ! Calculate the surface buoyancy flux term
            !--------------------------------------------------------
        zeta_s_fac = (one - zeta_s(i,j)) * (one - zeta_s(i,j))
        sf_term = a_ent_1 * max ( zero ,                                       &
                          ( (one - zeta_s_fac) * bflux_surf(i,j)               &
                            + zeta_s_fac * bflux_surf_sat(i,j) ) )
            !--------------------------------------------------------
            ! Calculate the surface shear term
            !--------------------------------------------------------
        if (fb_surf(i,j)  >=  zero) then
          sf_shear_term = a_ent_shr * v_s(i,j)*v_s(i,j)*v_s(i,j)               &
                          * rho_mix(i,j,k)  / zhsc(i,j)
        else
          sf_shear_term = zero
        end if
        if ( entr_smooth_dec == on .or. entr_smooth_dec == entr_taper_zh ) then
          ! taper surface terms to zero depending on svl_diff_frac
          sf_term =       sf_term       * svl_diff_frac(i,j)
          sf_shear_term = sf_shear_term * svl_diff_frac(i,j)
        end if
      else
        sf_term = zero
        sf_shear_term = zero
      end if
        !-----------------------------------------------------------
        ! Calculate the indirect radiative term
        !-----------------------------------------------------------
      zeta_r_sq = zeta_r_dsc(i,j)*zeta_r_dsc(i,j)
      ir_term = ( bt_dsct(i,j)*zeta_r_sq +                                     &
                    btt_dsct(i,j)*(one-zeta_r_sq) )                            &
                    * a_ent_1 * df_dsct_over_cp(i,j)
        !-----------------------------------------------------------
        ! Calculate the evaporative term
        !-----------------------------------------------------------
      if (db_dsct(i,j)  >   zero) then
        zr = sqrt( zc_dsc(i,j) / dscdepth(i,j) )
        evap_term = a_ent_2 * rho_mix(i,j,k)                                   &
                  * chi_s_dsct(i,j) * chi_s_dsct(i,j)                          &
                  * zr * zr * zr * db_dsct_cld(i,j)                            &
                  * sqrt( dscdepth(i,j) * db_dsct(i,j) )
      else
        evap_term = zero
      end if
        !-----------------------------------------------------------
        ! Combine forcing terms to calculate the representative
        ! velocity scales
        !-----------------------------------------------------------
      v_sum_dsc(i,j) = ( (sf_term + sf_shear_term +                            &
                            ir_term + evap_term)                               &
                * dscdepth(i,j) / (a_ent_1*rho_mix(i,j,k)) )**one_third
      v_top_dsc(i,j) =( (ir_term + evap_term) * dscdepth(i,j) /                &
                           (a_ent_1*rho_mix(i,j,k)) )**one_third
        !-----------------------------------------------------------
        ! Calculate the direct radiative term
        !-----------------------------------------------------------
      if (db_dsct(i,j)  >   zero) then
        dz_inv  = min( v_sum_dsc(i,j)*v_sum_dsc(i,j)/db_dsct(i,j),             &
                       100.0_r_bl )
        l_rad   = 15.0_r_bl * max( one, 200.0_r_bl/(zc_dsc(i,j)+one) )
        alpha_t = one - exp(-one_half*dz_inv/l_rad)
           ! Make enhancement due to buoyancy reversal feedback
        alpha_t = alpha_t + br_fback_dsc(i,j)*(one-alpha_t)
        dr_term = btc_top(i,j) * alpha_t * df_dsct_over_cp(i,j)
           !----------------------------------------------------------
           ! Finally combine terms to calculate the entrainment
           ! rate and mixing coefficients
           !----------------------------------------------------------
        zil_corr = c_t * ( (sf_term + sf_shear_term +                          &
                            ir_term + evap_term) /                             &
                 (rho_mix(i,j,k) * sqrt(dscdepth(i,j))) )**two_thirds
        rho_we_dsc(i,j) = ( sf_term + sf_shear_term                            &
                +           ir_term + evap_term + dr_term )                    &
                           / ( db_dsct(i,j) + zil_corr )
        rhokh_dsct_ent(i,j) = rho_we_dsc(i,j)/ rdz(i,j,k)
        rhokm_dsct_ent(i,j) = 0.75_r_bl * rhokh_dsct_ent(i,j)                  &
                              * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))       &
                              * rho_wet_tq(i,j,k-1) / rho_mix(i,j,k)
        if (.not. l_use_var_fixes) rhokm_top(i,j,k) = rhokm_dsct_ent(i,j)
      end if   ! test on DB_DSCT gt 0
    end if
  end do
end do
!$OMP end do

!  If there is no turbulence generation in DSC layer, ignore it.

!cdir collapse
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( dsc(i,j) .and. v_top_dsc(i,j)  <=  zero ) then
      dsc_removed(i,j) = 1
      dsc(i,j) = .false.
      ntdsc(i,j) = 0
      zhsc(i,j) = zero
      zc_dsc(i,j) = zero
      dscdepth(i,j) = zero
      coupled(i,j) = .false.
    end if
  end do
end do
!$OMP end do

! ----------------------------------------------------------------------
! 2.0_r_bl Estimate the depths of top-down and surface-up mixing.
!     These amount to diagnoses of recoupling and decoupling.
!     The K_H profiles are applied over layers such that the ratio
!        WBN_INT/WBP_INT = DEC_THRES (parameter),
!     where WBN_INT and WBP_INT are the magnitudes of the integrals of
!     the negative and positive parts, respectively, of the resulting
!     buoyancy flux profile (given by - KH * DB_FOR_FLUX).
! ----------------------------------------------------------------------
! 2.1 First test for well-mixed boundary layer
!     (ie. both KH profiles extending from cloud-top to the surface).
!     If the parcel ascent diagnosed:
!        DSC    - test for well-mixed up to ZHSC = recoupling
!        no DSC - test for well-mixed up to ZH   = decoupling
! -----------------------------------------------------------

if (model_type == mt_single_column) then
!$OMP do SCHEDULE(STATIC)
  ! Need to use pdims%k_end(=bl_levels-1) here because bl_levels
  ! in this routine is actually nl_bl_levels (< bl_levels)
  do k = pdims%k_start, pdims%k_end
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        wbmix(i,j,k) = zero  ! WB if were diag as well-mixed
        wbend(i,j,k) = zero  ! WB after dec diag
        wbend_sml(i,j,k) = zero ! WB after ksurf_iterate
      end do
    end do
  end do
!$OMP end do
end if ! model_type


! Default settings

!cdir collapse
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    zsml_top(i,j)  = zh(i,j)
    zsml_base(i,j) = 0.1_r_bl * zh(i,j)
    zdsc_base(i,j) = 0.1_r_bl * zhsc(i,j)
    z_top_lim(i,j) = zero ! initialising
    z_inv(i,j) = zero     ! inversion height (top of K profiles)
    zwb0(i,j)  = zero     ! height at which WB goes to zero
    wbp_int(i,j) = zero
    wbn_int(i,j) = zero
    wb_surf_int(i,j) = zero
    wb_dzrad_int(i,j) = zero
    dzrad(i,j) = 100.0_r_bl
    tothf_zi(i,j) = zero
    totqf_zi(i,j) = zero
    kstatus(i,j)= .true.
    kwb0(i,j)  = 2
    ntop(i,j)  = -1
    ksurf(i,j) = 1
    dec_thres(i,j) = dec_thres_cloud ! use cloudy by default
  end do
end do
!$OMP end do

! Find KSURF, the first theta-level above the surface layer

!$OMP do SCHEDULE(STATIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do k = 2, bl_levels
    !cdir collapse
    do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
      do i = pdims%i_start, pdims%i_end
        if ( z_tq(i,j,k-1)  <   0.1_r_bl*zh(i,j) ) ksurf(i,j) = k
      end do
    end do
  end do
end do
!$OMP end do

! Set flags for iterating wb integral to calculate depth of mixing,
! one each for KSURF and K_TOP.  Note these will be updated depending on
! what happpens on testing for a well-mixed layer in section 2.2.

!cdir collapse
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    test_well_mixed(i,j) = .false.
    ksurf_iterate(i,j)= .false.
    ktop_iterate(i,j) = .false.
    if ( ntdsc(i,j)  >   2 ) then
      ktop_iterate(i,j) = .true.
    end if
  end do ! I
end do ! J
!$OMP end do

!cdir collapse
!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( bflux_surf(i,j)  >   zero) then
        ! can only be coupled to an unstable SML
      if ( ntdsc(i,j)  >   2 ) then
          ! well-resolved DSC layer
          ! ...test for recoupling with SML
        test_well_mixed(i,j) = .true.
        z_inv(i,j)  = zhsc(i,j)
        z_cbase(i,j)= z_inv(i,j) - zc_dsc(i,j)
        cf_ml(i,j)  = cf_dsc(i,j)
        v_ktop(i,j) = v_top_dsc(i,j)
        v_ksum(i,j) = v_sum_dsc(i,j)
        df_ctop(i,j)= df_dsct_over_cp(i,j)
        rho_we(i,j) = rho_we_dsc(i,j)
        dsl(i,j)    = dsl_dsc(i,j)
        dqw(i,j)    = dqw_dsc(i,j)
        tothf_zi(i,j)= - rho_we(i,j)*dsl(i,j) + ft_nt_zhsc(i,j)
        totqf_zi(i,j)= - rho_we(i,j)*dqw(i,j) + fq_nt_zhsc(i,j)
        ntop(i,j)   = ntdsc(i,j)
          ! assuming wb goes to zero by the lowest of ZH
          ! or cloud-base, but above surface layer
      else if ( .not. dsc(i,j) .and. .not. cumulus(i,j) .and.                  &
                ntml(i,j)  >   2) then
          ! well-resolved SML
          ! ...test for decoupling
          ! Note: code can only deal with one DSC layer at a time so
          ! can't decouple SML if a DSC layer already exists.
          ! ---------------------------------------------------------
          ! If the BL layer is cloud-free then use a less restrictive
          ! threshold - ideally, the parcel ascent would have
          ! found the correct BL top in this case but this test is
          ! kept to keep negative buoyancy fluxes under control
          ! (ie. DEC_THRES_CLEAR=1 ensures wbn_int < |wbp_int|)
        if (abs(zc(i,j)) < rbl_eps) dec_thres(i,j) = dec_thres_clear
          ! ---------------------------------------------------------
        test_well_mixed(i,j) = .true.
        z_inv(i,j)  = zh(i,j)
        z_cbase(i,j)= z_inv(i,j) - zc(i,j)
        cf_ml(i,j)  = cf_sml(i,j)
        v_ktop(i,j) = v_top(i,j)
        v_ksum(i,j) = v_sum(i,j)
        df_ctop(i,j)= df_top_over_cp(i,j)
        rho_we(i,j) = rho_we_sml(i,j)
        dsl(i,j)    = dsl_sml(i,j)
        dqw(i,j)    = dqw_sml(i,j)
        tothf_zi(i,j)= - rho_we(i,j)*dsl(i,j) + ft_nt_zh(i,j)
        totqf_zi(i,j)= - rho_we(i,j)*dqw(i,j) + fq_nt_zh(i,j)
        ntop(i,j)   = ntml(i,j)
      end if
    end if
  end do ! I
end do ! J
!$OMP end do

if ( l_converge_ga ) then
  ! Set gradient adjustment terms consistent with zh = z_inv
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( test_well_mixed(i,j) ) then
        ! The gradient adjustment scales with 1/zh.  The terms were calculated
        ! in kmkhz_9c using zh_prev (mixed-layer height from previous timestep),
        ! but here we are setting the mixed-layer height to z_inv,
        ! so need to scale the gradient adjustment by zh_prev/z_inv
        ga_fac(i,j) = zh_prev(i,j) / z_inv(i,j)
      else
        ! Set to 1 where not used.
        ga_fac(i,j) = one
      end if
    end do
  end do
!$OMP end do
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        ! Updated gradient-adjusted db is non-adjusted value + ga_fac times
        ! the original adjustment.
        db_ga_dry_n(i,j,k) = db_noga_dry(i,j,k)                                &
          + ( db_ga_dry(i,j,k) - db_noga_dry(i,j,k) ) * ga_fac(i,j)
        db_ga_cld_n(i,j,k) = db_noga_cld(i,j,k)                                &
          + ( db_ga_cld(i,j,k) - db_noga_cld(i,j,k) ) * ga_fac(i,j)
      end do
    end do
  end do
!$OMP end do
else  ! ( .not. l_converge_ga )
  ! Use original gradient adjustment set consistent with zh_prev
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        db_ga_dry_n(i,j,k) = db_ga_dry(i,j,k)
        db_ga_cld_n(i,j,k) = db_ga_cld(i,j,k)
      end do
    end do
  end do
!$OMP end do
end if  ! ( l_converge_ga )

! Find kwb0, level with lowest positive cloud-free buoyancy gradient

!$OMP do SCHEDULE(STATIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do k = 2, bl_levels
    do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
      do i = pdims%i_start, pdims%i_end
        if (kstatus(i,j)) then
          if ( (db_ga_dry_n(i,j,k) <=  zero) .or.                              &
               (k >= ntml(i,j)) ) then
            kstatus(i,j)=.false.
            kwb0(i,j)=k
          end if
        end if
      end do
    end do
  end do
end do
!$OMP end do

! ----------------------------------------------------------------------
! 2.1.1 Estimate wb integral over radiatively cooled cloud-top region,
!       from Z_INV to Z_INV-DZRAD.
!       DZRAD taken to be constant (100m) for simplicity, but also
!       integration depth taken to extend down to at least
!       Z_TQ(NTML-1) but without going below cloud-base.
!       This code also invoked for cloud-free cases
!           - does this matter?
!           - only ignoring top grid level of integration which
!             shouldn't be important/relevant for cloud-free layers.
! ----------------------------------------------------------------------
!cdir collapse
!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( test_well_mixed(i,j) ) then
      dzrad(i,j) = 100.0_r_bl

      if ( dzrad_disc_opt == dzrad_ntm1 ) then
        ! Original discretisation of cloud-top radiatively-cooled layer;
        ! set the base of the layer at theta-level ntdsc-1
        ! (which is between 1.5 and 2.5 model-levels below cloud-top,
        !  depending on where zhsc is between z_uv(ntdsc+1) and z_uv(ntdsc+2)),
        ! then lower it by extra whole model-levels if needed to increase the
        ! depth to 100m.

        ntop(i,j) = ntop(i,j) - 1
        if ( ntop(i,j)  >   2 ) then
          inner_loop1: do while ( z_tq(i,j,ntop(i,j))  >  z_inv(i,j)-dzrad(i,j)&
               .and. z_tq(i,j,ntop(i,j)-1)  >   z_inv(i,j)-z_cbase(i,j))
            ntop(i,j) = ntop(i,j) - 1
            if (ntop(i,j) == 2 ) then
              exit inner_loop1
            end if
          end do inner_loop1
        end if
        dzrad(i,j) = z_inv(i,j) - z_tq(i,j,ntop(i,j))

      else if ( dzrad_disc_opt == dzrad_1p5dz ) then
        ! Alternative method;
        ! set the base of the layer 1.5 model-levels below cloud-top so that
        ! it always moves smoothly with zhsc instead of jumping suddenly
        ! (1.5 levels is the minimum distance below zhsc which will always
        !  keep the base of the layer below theta-level ntdsc, which is a
        !  requirement of the discretisation).
        ! The depth is also forced to be at least 100m.

        ! Find smoothly-varying height 1.5 model-levels below z_inv
        k_inv = ntop(i,j)
        if ( z_inv(i,j) < z_tq(i,j,k_inv+1) ) then
          interp = ( z_inv(i,j)        - z_uv(i,j,k_inv+1) )                   &
                 / ( z_tq(i,j,k_inv+1) - z_uv(i,j,k_inv+1) )
          z_pr = (one-interp) * z_tq(i,j,k_inv-1)                              &
               +      interp  * z_uv(i,j,k_inv)
        else
          interp = ( z_inv(i,j)        - z_tq(i,j,k_inv+1) )                   &
                 / ( z_uv(i,j,k_inv+2) - z_tq(i,j,k_inv+1) )
          z_pr = (one-interp) * z_uv(i,j,k_inv)                                &
               +      interp  * z_tq(i,j,k_inv)
        end if
        ! Go down to 100m below z_inv (or cloud-base) if that is lower,
        ! but don't allow lower than theta-level 1
        z_pr = min( z_pr, max( z_inv(i,j) - dzrad(i,j), z_cbase(i,j) ) )
        z_pr = max( z_pr, z_tq(i,j,1) )
        ! Set ntop to the rho-level straddling the base of the layer
        do while ( z_tq(i,j,ntop(i,j)-1) > z_pr )
          ntop(i,j) = ntop(i,j) - 1
        end do
        ! Set new cloud-top radiatively-cooled layer-depth
        dzrad(i,j) = z_inv(i,j) - z_pr

        ! Fraction of rho-level ntop below the dzrad layer
        interp = ( z_pr                - z_tq(i,j,ntop(i,j)-1) )               &
               / ( z_tq(i,j,ntop(i,j)) - z_tq(i,j,ntop(i,j)-1) )
        ! Add on SL and qw differences between the discontinuous inversion base
        ! and the base of the dzrad layer
        dsl(i,j) = dsl(i,j) + sl(i,j,k_inv) - (one-interp)*sl(i,j,ntop(i,j)-1) &
                                            -      interp *sl(i,j,ntop(i,j))
        dqw(i,j) = dqw(i,j) + qw(i,j,k_inv) - (one-interp)*qw(i,j,ntop(i,j)-1) &
                                            -      interp *qw(i,j,ntop(i,j))

      end if  ! ( dzrad_disc_opt )

      wsl_dzrad_int = dzrad(i,j) *                                             &
                      ( 0.66_r_bl*df_ctop(i,j) - rho_we(i,j)*dsl(i,j) )
      wqw_dzrad_int = - dzrad(i,j) * rho_we(i,j) * dqw(i,j)

      wb_dzrad_int(i,j) = wsl_dzrad_int * (                                    &
                            (one-cf_ml(i,j))*btm(i,j,ntop(i,j)+1) +            &
                              cf_ml(i,j)*btm_cld(i,j,ntop(i,j)+1) )            &
                        + wqw_dzrad_int * (                                    &
                            (one-cf_ml(i,j))*bqm(i,j,ntop(i,j)+1) +            &
                              cf_ml(i,j)*bqm_cld(i,j,ntop(i,j)+1) )
      wb_dzrad_int(i,j) = g * wb_dzrad_int(i,j)
      wb_dzrad_int(i,j) = max( zero, wb_dzrad_int(i,j) )

        ! Include WB_DZRAD_INT in WBP_INT as it set to be >0
      wbp_int(i,j) = wbp_int(i,j) + wb_dzrad_int(i,j)

    else

      wb_dzrad_int(i,j) = -one  ! To identify not calculated

    end if
  end do ! I
end do ! J
!$OMP end do

! For WB diagnostics, convert integrated WB to uniform profile
if (model_type == mt_single_column) then

!$OMP do SCHEDULE(STATIC)
  do k=2, bl_levels
    do j=pdims%j_start, pdims%j_end
      do i=pdims%i_start, pdims%i_end
        if ( test_well_mixed(i,j) .and. k >=  ntop(i,j)+1 ) then
          if ( k <= ntdsc(i,j)+1 ) then
            wbend(i,j,k) = wb_dzrad_int(i,j)/dzrad(i,j)
            wbmix(i,j,k) = wb_dzrad_int(i,j)/dzrad(i,j)
          else if ( k <= ntml(i,j)+1 ) then
            wbend(i,j,k) = wb_dzrad_int(i,j)/dzrad(i,j)
            wbmix(i,j,k) = wb_dzrad_int(i,j)/dzrad(i,j)
          end if
        end if
      end do ! i
    end do ! j
  end do ! k
!$OMP end do

end if ! model_type


! ----------------------------------------------------------------------
! 2.1.2 Estimate wb integral over surface layer
!       (and up to next theta-level, namely Z_TQ(KSURF) )
!       assuming a linear profile going to zero at ZWB0
! ----------------------------------------------------------------------
!cdir collapse

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    if ( test_well_mixed(i,j) ) then
      if ( kwb0(i,j)  ==  ntml(i,j) ) then
        zwb0(i,j) = zh(i,j)
      else if ( kwb0(i,j)  ==  2 ) then
        zwb0(i,j) = z_uv(i,j,2)
      else
        k=kwb0(i,j)
          ! now DB_GA_DRY(K) le 0 and DB_GA_DRY(K-1) gt 0
          ! so interpolate:
        db_ratio = db_ga_dry_n(i,j,k-1)                                        &
                / ( db_ga_dry_n(i,j,k-1) - db_ga_dry_n(i,j,k) )
        db_ratio = max( zero, db_ratio )  ! trap for rounding error
        zwb0(i,j)=z_uv(i,j,k-1) +                                              &
                  db_ratio * (z_uv(i,j,k)-z_uv(i,j,k-1))
      end if
      wb_surf_int(i,j) = bflux_surf(i,j) * z_tq(i,j,ksurf(i,j)) *              &
                    ( one - z_tq(i,j,ksurf(i,j))/(2.0_r_bl*zwb0(i,j)))
      wb_surf_int(i,j) = max( rbl_eps, wb_surf_int(i,j) )
    else
        ! only include surface layer contribution for unstable mixing
      wb_surf_int(i,j) = rbl_eps
    end if

    wbp_int(i,j) = wbp_int(i,j) + wb_surf_int(i,j) ! must be >0

  end do ! I
end do ! J
!$OMP end do


if (model_type == mt_single_column) then
!$OMP do SCHEDULE(STATIC)
  ! Save surface and bl-top layer integral for diagnostics
  do j=pdims%j_start, pdims%j_end
    do i=pdims%i_start, pdims%i_end
      wbmix(i,j,ksurf(i,j)) = wb_surf_int(i,j)
      wbend(i,j,ksurf(i,j)) = wb_surf_int(i,j)
    end do ! i
  end do ! j
!$OMP end do
end if ! model_type


! ----------------------------------------------------------------------
! 2.1.3 Loop over well-mixed boundary layer integrating WB
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(DYNAMIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do k = 2, bl_levels-1
    do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
      do i = pdims%i_start, pdims%i_end

        if ( test_well_mixed(i,j) ) then
          ! ----------------------------------------------
          ! worth testing layer as well-mixed to cloud-top
          ! ----------------------------------------------
          zb_ktop = 0.1_r_bl*z_inv(i,j)
          zinv_pr(i,j) = z_inv(i,j) - zb_ktop
          ! DB(K)is the K to K-1 difference and already
          ! integrated up to KSURF, so start this loop at KSURF+1
          if ( (k >= ksurf(i,j)+1) .and. (k <= ntop(i,j)) ) then

            khtop(i,j) = zero
            khsurf(i,j)= zero
            f2         = zero
            fsc        = zero

            z_pr = z_uv(i,j,k) - zb_ktop
            if (z_pr  >   zero .and. z_pr  <   zinv_pr(i,j)) then
              z_ratio = z_pr/zinv_pr(i,j)

              if (flux_grad  ==  LockWhelan2006) then
                khtop(i,j) = 3.6_r_bl * vkman * rho_mix(i,j,k) *v_ktop(i,j)    &
                                 * zinv_pr(i,j) * (z_ratio**3)                 &
                                 * ( (one-z_ratio)*(one-z_ratio) )
                f2 = rho_mix(i,j,k) * one_half * z_ratio                       &
                                         * 2.0_r_bl**( z_ratio**4 )
                if ( v_ksum(i,j)  >   zero ) then
                  fsc = rho_mix(i,j,k) * 3.5_r_bl*(v_ktop(i,j)/v_ksum(i,j))    &
                        * (z_ratio**3) * (one-z_ratio)
                end if
              else
                ! max to avoid rounding errors giving small negative numbers
                khtop(i,j) = g1 * vkman * rho_mix(i,j,k) * v_ktop(i,j)         &
                                * ((max(one - z_ratio,zero))**0.8_r_bl)        &
                                * z_pr * z_ratio
              end if
            end if

            z_pr = z_uv(i,j,k)
            if ( z_pr  <   z_inv(i,j)) then
              !--------------------------------
              ! include surface-driven profile
              !--------------------------------
              khsurf(i,j) = vkman * rho_mix(i,j,k) *                           &
                       w_h_top(i,j)*z_pr*( one - z_pr/z_inv(i,j) )             &
                                        *( one - z_pr/z_inv(i,j) )
            end if

            if (flux_grad  ==  LockWhelan2006) then
              wslng = (f2+fsc)*tothf_zi(i,j) - ft_nt(i,j,k)
              wqwng = (f2+fsc)*totqf_zi(i,j) - fq_nt(i,j,k)
            end if

            if ( z_tq(i,j,k)  <=  z_cbase(i,j) ) then
              ! Completely below cloud-base so use cloud-free formula
              wb_scld = khsurf(i,j) * db_ga_dry_n(i,j,k) +                     &
                        khtop(i,j) * db_noga_dry(i,j,k)
              if (flux_grad  ==  LockWhelan2006) then
                wb_scld = wb_scld + ( g/rdz(i,j,k) ) *                         &
                  ( btm(i,j,k-1)*wslng + bqm(i,j,k-1)*wqwng )
              end if
              wb_cld  = zero
            else if (z_tq(i,j,k-1)  >=  z_cbase(i,j)) then
              ! Completely above cloud-base so use cloudy formula
              wb_cld = ( khsurf(i,j) * db_ga_cld_n(i,j,k) +                    &
                         khtop(i,j)  * db_noga_cld(i,j,k) )
              if (flux_grad  ==  LockWhelan2006) then
                wb_cld = wb_cld + ( g/rdz(i,j,k) ) * (                         &
                         ( btm(i,j,k-1)*(one-cf_ml(i,j)) +                     &
                           btm_cld(i,j,k-1)*cf_ml(i,j) )*wslng +               &
                         ( bqm(i,j,k-1)*(one-cf_ml(i,j)) +                     &
                           bqm_cld(i,j,k-1)*cf_ml(i,j) )*wqwng )
              end if
              wb_scld = zero
            else
              ! cloud-base within this integration range
              ! so treat cloud and sub-cloud layer wb separately
              wb_scld = khsurf(i,j) * db_ga_dry_n(i,j,k) +                     &
                        khtop(i,j) * db_noga_dry(i,j,k)
              wb_cld = ( khsurf(i,j) * db_ga_cld_n(i,j,k) +                    &
                         khtop(i,j)  * db_noga_cld(i,j,k) )
              if (flux_grad  ==  LockWhelan2006) then
                wb_scld = wb_scld + ( g/rdz(i,j,k) ) *                         &
                  ( btm(i,j,k-1)*wslng + bqm(i,j,k-1)*wqwng )
                wb_cld = wb_cld + ( g/rdz(i,j,k) ) * (                         &
                         ( btm(i,j,k-1)*(one-cf_ml(i,j)) +                     &
                           btm_cld(i,j,k-1)*cf_ml(i,j) )*wslng +               &
                         ( bqm(i,j,k-1)*(one-cf_ml(i,j)) +                     &
                           bqm_cld(i,j,k-1)*cf_ml(i,j) )*wqwng )
              end if
              cld_frac = (z_tq(i,j,k)-z_cbase(i,j))                            &
                        /(z_tq(i,j,k)-z_tq(i,j,k-1))
              wb_cld  = cld_frac * wb_cld
              wb_scld = (one-cld_frac) * wb_scld
            end if

            if ( dzrad_disc_opt == dzrad_1p5dz .and. k == ntop(i,j) ) then
              ! At top of Sc layer, only include the part of the integral
              ! that is below the radiatively-cooled cloud-top layer.
              interp = ( z_inv(i,j) - dzrad(i,j) - z_tq(i,j,k-1) )             &
                     / ( z_tq(i,j,k) - z_tq(i,j,k-1) )
              wb_cld  = wb_cld  * interp
              wb_scld = wb_scld * interp
            end if

            if (wb_cld  >=  zero) then
              wbp_int(i,j) = wbp_int(i,j) + wb_cld
            else
              wbn_int(i,j) = wbn_int(i,j) - wb_cld
            end if
            if (wb_scld  >=  zero) then
              wbp_int(i,j) = wbp_int(i,j) + wb_scld
            else
              wbn_int(i,j) = wbn_int(i,j) - wb_scld
            end if

            if (model_type == mt_single_column) then
              wbmix(i,j,k) = wb_cld+wb_scld
              wbend(i,j,k) = wb_cld+wb_scld
            end if

          end if ! K

        end if ! TEST_WELL_MIXED
      end do ! I
    end do ! J
  end do ! K
end do
!$OMP end do

! ----------------------------------------------------------------------
! 2.1.4 Test WB_Ratio to see if well-mixed layer allowed
!-----------------------------------------------------------------------
!cdir collapse

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if (l_wtrac) iset_wtrac(i,j) = 0
    if ( test_well_mixed(i,j) ) then

      wb_ratio(i,j) = wbn_int(i,j)/wbp_int(i,j)

      if ( wb_ratio(i,j)  <=  dec_thres(i,j) ) then
          ! No need to test depth of mixing any further as within
          ! well-mixed layer buoyancy flux integral criteria.
          ! SML will simply stay well-mixed (and so use defaults)
        ksurf_iterate(i,j)= .false.
        ktop_iterate(i,j) = .false.
        if ( dsc(i,j) ) then
          ! Recouple DSC with SML:
          dsc_removed(i,j) = 2 ! set indicator flag
          ! move surface driven entrainment
          ! RHOKH(z_i) = rho * w_e * DZL and w_e ~ 1/DB_TOP, so:
          if ( db_top(i,j) >  zero .and. db_dsct(i,j) >  0.01_r_bl ) then
                                        ! can't calc Zil. term
            rhokh_surf_ent(i,j) = rhokh_surf_ent(i,j) *                        &
                   ( rho_mix(i,j,ntdsc(i,j)+1) * db_top(i,j) *                 &
                     rdz(i,j,ntml(i,j)+1) ) /                                  &
                   ( rho_mix(i,j,ntml(i,j)+1) * db_dsct(i,j) *                 &
                                        rdz(i,j,ntdsc(i,j)+1) )
            rhokm_surf_ent(i,j) = rhokm_surf_ent(i,j) *                        &
                   ( rho_wet_tq(i,j,ntdsc(i,j)) * db_top(i,j) *                &
                     rdz(i,j,ntml(i,j)+1) ) /                                  &
                   ( rho_wet_tq(i,j,ntml(i,j)) * db_dsct(i,j) *                &
                                      rdz(i,j,ntdsc(i,j)+1) )
            if (.not. l_use_var_fixes)                                         &
                  rhokm(i,j,ntdsc(i,j)+1) = rhokm_surf_ent(i,j)
          end if
          ! redesignate top-driven entrainment at ZHSC
          ! (ignore that calculated at ZH)
          rhokh_top_ent(i,j) = rhokh_dsct_ent(i,j)
          rhokm_top_ent(i,j) = rhokm_dsct_ent(i,j)
          zh(i,j) = zhsc(i,j)
          ntml(i,j) = ntdsc(i,j)
          v_top(i,j) = v_top_dsc(i,j)
          zsml_base(i,j) = 0.1_r_bl * zh(i,j)
          zc(i,j) = zc_dsc(i,j)
          zhsc(i,j) = zero
          ntdsc(i,j) = 0
          v_top_dsc(i,j) = zero
          zdsc_base(i,j) = zero
          zc_dsc(i,j)    = zero
          ft_nt_zh(i,j)   = ft_nt_zhsc(i,j)
          ft_nt_zhsc(i,j) = zero
          fq_nt_zh(i,j)   = fq_nt_zhsc(i,j)
          fq_nt_zhsc(i,j) = zero
          if (l_wtrac) iset_wtrac(i,j) = 1
          dsc(i,j) = .false.
          dzh(i,j) = rmdi   ! SML inversion subsumed into new mixed layer
          cumulus(i,j) = .false.
          coupled(i,j) = .false.
        end if  ! recoupled DSC layer
      else   ! buoyancy flux threshold violated

        if ( l_use_sml_dsc_fixes ) then
          ! Under fix l_use_sml_dsc_fixes:
          ! If there is already a DSC layer diagnosed above the existing
          ! ntml, we have so far just tested for recoupling of that DSC
          ! with the SML; we have not tested whether the existing SML
          ! should split off another decoupled layer.
          ! In this case, we don't need to turn on ktop_iterate
          ! (already setup earlier), and we should only turn on ksurf_iterate
          ! if the SML is cumulus capped and using buoyancy flux integration
          ! at cumulus points in general.

          if ( dsc(i,j) ) then
            ! Existing DSC layer confirmed not to be merged with the SML.
            if ( cumulus(i,j) .and. ( kprof_cu == buoy_integ .or.              &
                                      kprof_cu == buoy_integ_low ) ) then
              ! Cumulus layer under the DSC; only turn on ksurf_iterate if
              ! using appropriate kprof_cu option.
              ksurf_iterate(i,j) = .true.
              dec_thres(i,j) = dec_thres_cu
            end if
          else
            ! No existing DSC so we're splitting a new DSC off the SML;
            ! turn on ktop_iterate to find the new DSC's base
            ktop_iterate(i,j)  = .true.
            if ( (.not. cumulus(i,j)) .or. kprof_cu == buoy_integ .or.         &
                                           kprof_cu == buoy_integ_low ) then
              ! Switch to using b-flux integration to find new lower SML-top,
              ! but not in Cu layers unless using appropriate kprof_cu option.
              ksurf_iterate(i,j) = .true.
              if ( cumulus(i,j) )  dec_thres(i,j) = dec_thres_cu
            end if
          end if

        else  ! ( .not. l_use_sml_dsc_fixes )
          ! If the bug-fix is off, we turn on ksurf_iterate at non-cumulus
          ! points even if dsc was already true
          ! (when all we've done is confirm the DSC shouldn't merge with the
          !  SML, so we had no reason to be meddling with the SML).
          ! This creates some inconsistencies down the line...

          !---------------------------------
          ! Extent of mixing must be reduced
          !---------------------------------
          if ( .not. cumulus(i,j) .or. kprof_cu == buoy_integ .or.             &
                                       kprof_cu == buoy_integ_low)             &
                 ksurf_iterate(i,j) = .true.
          if ( cumulus(i,j) .and. ( kprof_cu == buoy_integ .or.                &
                                    kprof_cu == buoy_integ_low) )              &
                 dec_thres(i,j) = dec_thres_cu
          ktop_iterate(i,j)  = .true.

        end if  ! ( l_use_sml_dsc_fixes )

        if ( .not. dsc(i,j) ) then
          if (l_use_var_fixes .and. v_top(i,j) <= zero) then
            ! If v_top=0 we have no mechanism to generate turbulence in a DSC
            ktop_iterate(i,j) = .false.
            z_top_lim(i,j)    = zh(i,j) ! needs to be set here as zhsc not used
          else
            ! Set up a `COUPLED' decoupled layer,
            !   implies no explicit `entrainment' at ZH.
            ! Note a new ZH (and thence NTML) will be calculated by
            ! wb integral iteration.
            dsc(i,j) = .true.
            coupled(i,j) = .true.
            ntdsc(i,j) = ntml(i,j)
            zhsc(i,j) = zh(i,j)
            zc_dsc(i,j) = zc(i,j)
            v_top_dsc(i,j) = v_top(i,j)
            if ( l_use_sml_dsc_fixes ) then
              ! Need to ensure SML-top entrainment flux is turned off when
              ! using buoyancy-flux integration to find a new zh
              v_top(i,j) = zero
              ! Set coupling strength to maximal when just diagnosed decoupled
              svl_diff_frac(i,j) = one
            end if
            v_sum_dsc(i,j) = v_sum(i,j)
            ft_nt_zhsc(i,j) = ft_nt_zh(i,j)
            fq_nt_zhsc(i,j) = fq_nt_zh(i,j)
            if (l_wtrac) iset_wtrac(i,j) = 2
            ! put all entrainment into RHOKH_TOP
            rhokh_dsct_ent(i,j) = rhokh_top_ent(i,j)                           &
                                + rhokh_surf_ent(i,j)
            rhokh_top_ent(i,j) = zero
            rhokh_surf_ent(i,j) = zero
            rhokm_dsct_ent(i,j) = rhokm_top_ent(i,j)                           &
                                + rhokm_surf_ent(i,j)
            rhokm_top_ent(i,j) = zero
            rhokm_surf_ent(i,j) = zero
            if (.not. l_use_var_fixes) then
              rhokm_top(i,j,ntml(i,j)+1) = rhokm_top(i,j,ntml(i,j)+1)          &
                                         + rhokm(i,j,ntml(i,j)+1)
              rhokm(i,j,ntml(i,j)+1) = zero
            end if
          end if ! test on l_use_var_fixes and v_top <= zero
        end if ! test no not.dsc
      end if   ! test on WB_RATIO le DEC_THRES
    end if   ! testing for well-mixed layer (TEST_WELL_MIXED)

  end do ! I
end do ! J
!$OMP end do

! Set water tracer fields in the same manner as water in last loop
if (l_wtrac) then
  do i_wt = 1, n_wtrac
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if (iset_wtrac(i,j) == 1) then
          fq_nt_zh_wtrac(i,j,i_wt)   = fq_nt_zhsc_wtrac(i,j,i_wt)
          fq_nt_zhsc_wtrac(i,j,i_wt) = zero
        else if (iset_wtrac(i,j) == 2) then
          fq_nt_zhsc_wtrac(i,j,i_wt) = fq_nt_zh_wtrac(i,j,i_wt)
        end if
      end do
    end do
!$OMP end do
  end do

end if    ! l_wtrac

! ----------------------------------------------------------------------
! 2.2 Start iteration to find top of surface-driven mixing, ZSML_TOP,
!     within predetermined maximum and minimum height limits.
!     The solution is the height that gives WB_RATIO = DEC_THRES.
!     Procedure used makes 3 sweeps (up, down and up again), using
!     progressively smaller increments (Z_INC), each time stopping when
!     the buoyancy flux threshold or the height limits are reached.
!--------------------------------------------------------------------
!     If boundary layer is stable then ignore surface driven mixing.
!--------------------------------------------------------------------
!cdir collapse

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    if ( ksurf_iterate(i,j) ) then
      !-----------------------------------------------------------
      ! Mixing must extend just above surface layer
      ! (not clear precisely how to define this here: for now use
      !  KSURF calculated from ZH)
      !-----------------------------------------------------------
      z_bot_lim(i,j)=z_uv(i,j,ksurf(i,j)+1)                                    &
             + 0.1_r_bl * (z_uv(i,j,ksurf(i,j)+2)-z_uv(i,j,ksurf(i,j)+1))
      if ( l_converge_ga ) then
        ! Allow K-surf to go up to the Sc-top (increasing zsml_top within
        ! the cloud-top layer can still make the buoyancy fluxes below more
        ! negative because this reduces the gradient adjustment)
        z_top_lim(i,j)=max( z_top_lim(i,j),                                    &
                            max( z_bot_lim(i,j), zhsc(i,j) ) )
        ! ksurf_iterate currently only true if test_well_mixed was true and
        ! decoupling was diagnosed
        if (cumulus(i,j)) then
          ! ignore cloudy component of wb for cumulus sml iteration as the
          ! cloudy flux is carried by the convection scheme
          z_cbase(i,j) = z_top_lim(i,j)
        else
          ! sml has decoupled from the cloud layer above.  For consistency with
          ! the test_well_mixed calculation, continue to use saturated wb above
          ! the now decoupled layer's cloud base.
          z_cbase(i,j) = zhsc(i,j) - zc_dsc(i,j)
        end if
      else  ! ( .not. l_converge_ga )
        ! limit K-surf to below cloud-top radiatively cooled layer
        ! or original zh (held in z_top_lim) if no top-driven turbulence
        z_top_lim(i,j)=max( z_top_lim(i,j),                                    &
                            max( z_bot_lim(i,j), zhsc(i,j) - dzrad(i,j) ) )
        ! Always use saturated wb above decoupled-layer cloud-base
        z_cbase(i,j) = zhsc(i,j) - zc_dsc(i,j)
      end if  ! ( l_converge_ga )
      !-----------------------------------------------------
      ! Initial increment to ZSML_TOP found by dividing
      ! up depth of layer within which it is allowed:
      ! Start with ZSML_TOP at lower limit and work upwards
      !-----------------------------------------------------
      z_inc(i,j)=(z_top_lim(i,j)-z_bot_lim(i,j)) / real(n_steps, r_bl)
      zsml_top(i,j) = z_bot_lim(i,j)

      wb_ratio(i,j) = dec_thres(i,j) - one ! to be < DEC_THRES

    end if ! KSURF_ITERATE

  end do
end do
!$OMP end do

! ----------------------------------------------------------------------
! 2.2.1 Include ksurf iteration for SML under cumulus too
!-----------------------------------------------------------------------
if (kprof_cu == buoy_integ .or. kprof_cu == buoy_integ_low) then

  lcl_fac = one
  if (kprof_cu == buoy_integ_low) lcl_fac = one_half

!$OMP do SCHEDULE(DYNAMIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( cumulus(i,j) .and. .not. dsc(i,j) .and.                             &
                              .not. ksurf_iterate(i,j)) then
        ! pure cumulus layer and necessary parameters not already set up
        ksurf_iterate(i,j)  = .true.
        dec_thres(i,j) = dec_thres_cu ! use cu threshold
        ! limit top of K_surf to be between...
        z_bot_lim(i,j) = lcl_fac*z_lcl(i,j)  ! ...a fraction of the LCL and...
        z_top_lim(i,j) = z_lcl(i,j) + 1000.0_r_bl ! ...1km above the LCL
        ntop(i,j)=ntml(i,j)
        do while ( z_uv(i,j,ntop(i,j)+1) <= z_top_lim(i,j) .and.               &
                   ntop(i,j)+1 < bl_levels-2 )
          ntop(i,j) = ntop(i,j) + 1  ! z_uv(ntop+1) > z_top_lim
        end do
        z_top_lim(i,j)=z_uv(i,j,ntop(i,j)) ! highest z_uv below z_top_lim
        z_cbase(i,j) = z_top_lim(i,j)
        !-----------------------------------------------------
        ! Initial increment to ZSML_TOP found by dividing
        ! up depth of layer within which it is allowed:
        ! Start with ZSML_TOP at lower limit and work upwards
        !-----------------------------------------------------
        z_inc(i,j)=(z_top_lim(i,j)-z_bot_lim(i,j))                             &
                     / real(n_steps, r_bl)
        zsml_top(i,j) = z_bot_lim(i,j)
        wb_ratio(i,j) = dec_thres(i,j) - one ! to be < DEC_THRES
      end if
    end do ! I
  end do ! J
!$OMP end do

end if  ! test in kprof_cu

! ----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    l = i - pdims%i_start + 1 + pdims%i_len * (j - pdims%j_start)
    ind_todo(l) = l
    up(l)       = 1
    if (ksurf_iterate(i,j)) then
      to_do(l) = .true.
    else
      to_do(l) = .false.
    end if
  end do
end do
!$OMP end do

!$OMP MASTER
c_len=pdims%i_len*pdims%j_len
!$OMP end MASTER
!$OMP BARRIER

do n_sweep = 1, num_sweeps_bflux
!$OMP BARRIER

!$OMP MASTER

          ! Compress to_do and ind_todo (will have new length c_len)
  call excfnl_cci(c_len, to_do, ind_todo)

      ! Restart inner interation with the points of outer
  c_len_i = c_len
  todo_inner(1:c_len_i) = to_do(1:c_len_i)
  ind_todo_i(1:c_len_i) = ind_todo(1:c_len_i)

!$OMP end MASTER
!$OMP BARRIER

  do ns = 1, n_steps
!$OMP BARRIER

!$OMP MASTER

              ! Calculate active elements and compress
    call excfnl_compin(up, wb_ratio, dec_thres, 1,                             &
                       c_len_i, ind_todo_i, todo_inner)

!$OMP end MASTER
!$OMP BARRIER

    !cdir nodep
!$OMP do SCHEDULE(STATIC)
    do ic = 1, c_len_i
      j1=(ind_todo_i(ic)-1)/pdims%i_end+1
      i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

      zsml_top(i1,j1)=zsml_top(i1,j1)+z_inc(i1,j1)

      ! The gradient adjustment scales with 1/zh.  The terms were calculated
      ! in kmkhz_9c using zh_prev (mixed-layer height from previous timestep),
      ! but here we are setting the mixed-layer height to a new zsml_top,
      ! so need to scale the gradient adjustment by zh_prev/zsml_top
      ga_fac(i1,j1) = zh_prev(i1,j1) / zsml_top(i1,j1)

          ! assume wb goes to zero at ZSML_TOP
      wb_surf_int(i1,j1) =                                                     &
           bflux_surf(i1,j1) * z_tq(i1,j1,ksurf(i1,j1)) *                      &
           ( one - z_tq(i1,j1,ksurf(i1,j1))/                                   &
                                      (2.0_r_bl*zsml_top(i1,j1)) )
      wb_surf_int(i1,j1) = max(rbl_eps,wb_surf_int(i1,j1))

          ! Note: WB_DZRAD_INT not included as K_SURF restricted
          !       to below zi-dzrad
      wbp_int(i1,j1) = wb_surf_int(i1,j1)  ! must be > 0
      wbn_int(i1,j1) = zero

      z_inv(i1,j1) = zsml_top(i1,j1)

      if (model_type == mt_single_column) then
        wbend_sml(i1,j1,ksurf(i1,j1)) = wb_surf_int(i1,j1)
      end if

    end do ! ic c_len_i
!$OMP end do

    !..Integrate buoyancy flux profile given this ZSML_TOP

!$OMP do SCHEDULE(STATIC)
    do jj = 1, c_len_i, bl_segment_size
      do k = 2, bl_levels-1
        !cdir nodep
        do ic = jj, min(jj+bl_segment_size-1, c_len_i)
          j1=(ind_todo_i(ic)-1)/pdims%i_end+1
          i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

          if ( k  >=  ksurf(i1,j1)+1 .and.                                     &
               k  <=  ntop(i1,j1) ) then

            if ( l_converge_ga ) then
              ! Updated gradient-adjusted db is non-adjusted value + ga_fac
              ! times the original adjustment.
              db_ga_dry_n(i1,j1,k) = db_noga_dry(i1,j1,k)                      &
                + ( db_ga_dry(i1,j1,k) - db_noga_dry(i1,j1,k) ) * ga_fac(i1,j1)
              db_ga_cld_n(i1,j1,k) = db_noga_cld(i1,j1,k)                      &
                + ( db_ga_cld(i1,j1,k) - db_noga_cld(i1,j1,k) ) * ga_fac(i1,j1)
            end if

            z_pr = z_uv(i1,j1,k)
            if (z_pr  <   z_inv(i1,j1)) then
              if ( l_use_sml_dsc_fixes ) then
                ! Include factors vkman * rho_mix, to be consistent with
                ! the test_well_mixed code.
                kh_surf(i1,j1,k) = vkman * rho_mix(i1,j1,k) * w_h_top(i1,j1)   &
                                    *z_pr*(one-z_pr/z_inv(i1,j1) )             &
                                         *(one-z_pr/z_inv(i1,j1) )
              else
                ! Original version is missing the factors vkman * rho_mix,
                ! which makes kh_surf about a factor of 2 larger here than
                ! it is in the test_well_mixed code.
                ! This only matters when the kh_surf and DSC kh profiles
                ! overlap in the vertical, as what is stored in kh_surf
                ! gets added onto the DSC kh in the ktop iteration code
                ! (and DSC kh does include these factors).
                kh_surf(i1,j1,k) = w_h_top(i1,j1)                              &
                                    *z_pr*(one-z_pr/z_inv(i1,j1) )             &
                                         *(one-z_pr/z_inv(i1,j1) )
              end if
            else
              kh_surf(i1,j1,k) = zero
            end if
              !-----------------------------------------------------
              ! No F2 or FSc terms here because we're only
              ! considering effects driven from the surface
              !-----------------------------------------------------
            if (z_cbase(i1,j1)  >   z_tq(i1,j1,k)) then
                ! cloud-base above this range so use dry WB
              wb_scld= kh_surf(i1,j1,k) * db_ga_dry_n(i1,j1,k)
              wb_cld = zero
            else if (z_cbase(i1,j1)  <   z_tq(i1,j1,k-1)) then
                ! cloud-base below this range so use cloudy WB
              wb_cld = kh_surf(i1,j1,k) * db_ga_cld_n(i1,j1,k)
              wb_scld=zero
            else
                ! cloud-base within this integration range
                ! so treat cloud and sub-cloud layer wb separately
              cld_frac = (z_tq(i1,j1,k)-z_cbase(i1,j1))                        &
                        /(z_tq(i1,j1,k)-z_tq(i1,j1,k-1))
              wb_cld  = cld_frac                                               &
                         * kh_surf(i1,j1,k)*db_ga_cld_n(i1,j1,k)
              wb_scld = (one-cld_frac)                                         &
                         * kh_surf(i1,j1,k)*db_ga_dry_n(i1,j1,k)
            end if

            if (wb_cld  >=  zero) then
              wbp_int(i1,j1)= wbp_int(i1,j1) + wb_cld
            else
              wbn_int(i1,j1)= wbn_int(i1,j1) - wb_cld
            end if
            if (wb_scld  >=  zero) then
              wbp_int(i1,j1)= wbp_int(i1,j1) + wb_scld
            else
              wbn_int(i1,j1) = wbn_int(i1,j1)- wb_scld
            end if

            if (model_type == mt_single_column) then
              wbend_sml(i1,j1,k) = wb_cld + wb_scld
            end if

          end if ! K

        end do ! ic c_len_i
      end do ! K

      do ic = jj, min(jj+bl_segment_size-1, c_len_i)
        j1=(ind_todo_i(ic)-1)/pdims%i_end+1
        i1=ind_todo_i(ic)-(j1-1)*pdims%i_end
        wb_ratio(i1,j1) = wbn_int(i1,j1)/wbp_int(i1,j1)
      end do ! ic c_len_i
    end do
!$OMP end do

  end do  ! loop stepping up through ML (N_steps)

  !cdir nodep
!$OMP do SCHEDULE(STATIC)
  do ic = 1, c_len
    l=ind_todo(ic)
    j1=(l-1)/pdims%i_end+1
    i1=l-(j1-1)*pdims%i_end

    !..sub-divide current Z_INC into one more part than there will be steps
    !..as there is no need to calculate WB for ZSML_TOP at a current Z_INC
    z_inc(i1,j1)= z_inc(i1,j1)/real(n_steps+1, r_bl)

    if ((up(l) == 1 .and. wb_ratio(i1,j1) >= dec_thres(i1,j1)) .or.            &
               ! hit thres while working up
        (up(l) == 0 .and. wb_ratio(i1,j1) <= dec_thres(i1,j1))) then
               ! hit thres while working down
      up(l) = 1-up(l)   ! change direction of sweep
      z_inc(i1,j1)= - z_inc(i1,j1)
    else if (zsml_top(i1,j1) >= z_top_lim(i1,j1)-one) then
               ! hit upper height limit (give-or-take 1m) without
               ! reaching threshold
      to_do(ic)=.false.
      zsml_top(i1,j1) = z_top_lim(i1,j1)
    else if (zsml_top(i1,j1)  <=  z_bot_lim(i1,j1)+ one) then
               ! hit lower height limit (give-or-take 1m) without
               ! reaching threshold
      to_do(ic)=.false.
      zsml_top(i1,j1) = z_bot_lim(i1,j1)
    end if

    !..Note that if the threshold has not been passed then the next sweep
    !..continues in the same direction (but with reduced increment).

  end do ! c_len
!$OMP end do

end do ! n_sweep

!$OMP do SCHEDULE(DYNAMIC)
do jj = pdims%j_start, pdims%j_end,bl_segment_size
  do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
    do i = pdims%i_start, pdims%i_end
      ntml_new(i,j) = 2
      status_ntml(i,j)=.true.
    end do
  end do

  do k = 2, bl_levels-2
    !cdir collapse
    do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
      do i = pdims%i_start, pdims%i_end
        if ( ksurf_iterate(i,j) .and. status_ntml(i,j) ) then
          ! -------------
          ! find new NTML
          ! -------------
          if (z_uv(i,j,k+1)  <  zsml_top(i,j)) then
            ntml_new(i,j) = k+1
          else
            status_ntml(i,j)=.false.
          end if
          ! --------------------------------------------------------
          ! Rounding error previously found to give
          !      ZSML_TOP > Z_TOP_LIM = ZHSC
          ! Test on ZSML_TOP hitting thresholds consequently changed
          ! but also include the following failsafe tests here.
          ! --------------------------------------------------------
          if (dsc(i,j)) then
            ntml(i,j) = min( ntdsc(i,j), ntml_new(i,j)-1 )
            zh(i,j)   = min(  zhsc(i,j), zsml_top(i,j) )
          else if (.not. status_ntml(i,j)) then
            if (kprof_cu == buoy_integ_low) then
              ! just use zh from buoyancy flux integration
              ntml(i,j) = ntml_new(i,j)-1
              zh(i,j)   = zsml_top(i,j)
            else
              ! for kprof_cu eq buoy_integ use max of ksurf top calculated
              ! from buoyancy flux integration and original LCL based values
              ntml(i,j) =  max( ntml(i,j), ntml_new(i,j)-1 )
              if (l_use_var_fixes) then
                zh(i,j) = max( zh(i,j), zsml_top(i,j) )
              else
                zh(i,j) = zsml_top(i,j)
              end if
            end if
          end if

        end if  ! KSURF_ITERATE true

      end do
    end do
  end do
end do
!$OMP end do

! ----------------------------------------------------------------------
! 2.3 Now repeat the above procedure to find the base of the
!     top-driven K profile, ZDSC_BASE, in a decoupled cloud layer
! ----------------------------------------------------------------------
!cdir collapse
!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end

    if ( ktop_iterate(i,j) ) then

        ! Lower limit on base of DSC layer
      z_bot_lim(i,j) = 0.1_r_bl * zh(i,j)
        ! Upper limit on base of DSC layer
      z_top_lim(i,j) = zhsc(i,j) - dzrad(i,j)
        ! If cumulus limit base of top-driven mixing to above ZH
        ! (not applying this limit if converging the gradient adjustment,
        !  as in this case we need to allow kh_surf to go as high as it wants
        !  so the DSC kh needs to be able to overlap with it).
      if ( cumulus(i,j) .and. (.not. l_converge_ga) ) then
        z_bot_lim(i,j) = min( zh(i,j), z_top_lim(i,j) )
      end if

      z_cbase(i,j) = zhsc(i,j) - zc_dsc(i,j)

      !..Divide up depth of layer within which ZDSC_BASE is allowed
      z_inc(i,j) = (z_top_lim(i,j)-z_bot_lim(i,j)) / real(n_steps, r_bl)
      zdsc_base(i,j) = z_bot_lim(i,j)
                         ! will start at Z_BOT_LIM+Z_INC

      if (l_reset_dec_thres) then
        dec_thres(i,j) = dec_thres_cloud  ! reset after potential use
                                          ! of dec_thres_cu or
                                          ! dec_thres_clear above
      end if
      wb_ratio(i,j) = dec_thres(i,j) + one ! to be > DEC_THRES

    end if ! KTOP_ITERATE

  end do

  do i = pdims%i_start, pdims%i_end
    if ( ktop_iterate(i,j) .and. wb_dzrad_int(i,j)  <   zero ) then
        !-------------------------------------------------------------
        ! Estimation of wb integral over radiatively cooled cloud-top
        ! region not yet performed (ie. DSC over stable surface) so
        ! do it now.
        !-------------------------------------------------------------
      z_inv(i,j)  = zhsc(i,j)
      z_cbase(i,j)= z_inv(i,j) - zc_dsc(i,j)
      cf_ml(i,j)  = cf_dsc(i,j)
      df_ctop(i,j)= df_dsct_over_cp(i,j)
      rho_we(i,j) = rho_we_dsc(i,j)
      dsl(i,j)    = dsl_dsc(i,j)
      dqw(i,j)    = dqw_dsc(i,j)
      tothf_zi(i,j)= - rho_we(i,j)*dsl(i,j) + ft_nt_zhsc(i,j)
      totqf_zi(i,j)= - rho_we(i,j)*dqw(i,j) + fq_nt_zhsc(i,j)
      ntop(i,j)   = ntdsc(i,j)
      dzrad(i,j)  = 100.0_r_bl
      if ( dzrad_disc_opt == dzrad_ntm1 ) then
        ! Original discretisation of cloud-top radiatively-cooled layer;
        ! set the base of the layer at theta-level ntdsc-1
        ! (which is between 1.5 and 2.5 model-levels below cloud-top,
        !  depending on where zhsc is between z_uv(ntdsc+1) and z_uv(ntdsc+2)),
        ! then lower it by extra whole model-levels if needed to increase the
        ! depth to 100m.

        ntop(i,j) = ntop(i,j) - 1
        if ( ntop(i,j)  >   1 ) then
          inner_loop2: do while ( z_tq(i,j,ntop(i,j))  >  z_inv(i,j)-dzrad(i,j)&
               .and. z_tq(i,j,ntop(i,j)-1)  >   z_inv(i,j)-z_cbase(i,j))
            ntop(i,j) = ntop(i,j) - 1
            if (ntop(i,j) == 1 ) then
              exit inner_loop2
            end if
          end do inner_loop2
        end if
        dzrad(i,j) = z_inv(i,j) - z_tq(i,j,ntop(i,j))

      else if ( dzrad_disc_opt == dzrad_1p5dz ) then
        ! Alternative method;
        ! set the base of the layer 1.5 model-levels below cloud-top so that
        ! it always moves smoothly with zhsc instead of jumping suddenly
        ! (1.5 levels is the minimum distance below zhsc which will always
        !  keep the base of the layer below theta-level ntdsc, which is a
        !  requirement of the discretisation).
        ! The depth is also forced to be at least 100m.

        ! Find smoothly-varying height 1.5 model-levels below z_inv
        k_inv = ntop(i,j)
        if ( z_inv(i,j) < z_tq(i,j,k_inv+1) ) then
          interp = ( z_inv(i,j)        - z_uv(i,j,k_inv+1) )                   &
                 / ( z_tq(i,j,k_inv+1) - z_uv(i,j,k_inv+1) )
          z_pr = (one-interp) * z_tq(i,j,k_inv-1)                              &
               +      interp  * z_uv(i,j,k_inv)
        else
          interp = ( z_inv(i,j)        - z_tq(i,j,k_inv+1) )                   &
                 / ( z_uv(i,j,k_inv+2) - z_tq(i,j,k_inv+1) )
          z_pr = (one-interp) * z_uv(i,j,k_inv)                                &
               +      interp  * z_tq(i,j,k_inv)
        end if
        ! Go down to 100m below z_inv (or cloud-base) if that is lower,
        ! but don't allow lower than theta-level 1
        z_pr = min( z_pr, max( z_inv(i,j) - dzrad(i,j), z_cbase(i,j) ) )
        z_pr = max( z_pr, z_tq(i,j,1) )
        ! Set ntop to the rho-level straddling the base of the layer
        do while ( z_tq(i,j,ntop(i,j)-1) > z_pr )
          ntop(i,j) = ntop(i,j) - 1
        end do
        ! Set new cloud-top radiatively-cooled layer-depth
        dzrad(i,j) = z_inv(i,j) - z_pr

        ! Fraction of rho-level ntop below the dzrad layer
        interp = ( z_pr                - z_tq(i,j,ntop(i,j)-1) )               &
               / ( z_tq(i,j,ntop(i,j)) - z_tq(i,j,ntop(i,j)-1) )
        ! Add on SL and qw differences between the discontinuous inversion base
        ! and the base of the dzrad layer
        dsl(i,j) = dsl(i,j) + sl(i,j,k_inv) - (one-interp)*sl(i,j,ntop(i,j)-1) &
                                            -      interp *sl(i,j,ntop(i,j))
        dqw(i,j) = dqw(i,j) + qw(i,j,k_inv) - (one-interp)*qw(i,j,ntop(i,j)-1) &
                                            -      interp *qw(i,j,ntop(i,j))

      end if  ! ( dzrad_disc_opt )

      wsl_dzrad_int = dzrad(i,j) *                                             &
                      ( 0.66_r_bl*df_ctop(i,j) - rho_we(i,j)*dsl(i,j) )
      wqw_dzrad_int = - dzrad(i,j) * rho_we(i,j) * dqw(i,j)

      wb_dzrad_int(i,j) = wsl_dzrad_int * (                                    &
                            (one-cf_ml(i,j))*btm(i,j,ntop(i,j)+1) +            &
                              cf_ml(i,j)*btm_cld(i,j,ntop(i,j)+1) )            &
                        + wqw_dzrad_int * (                                    &
                            (one-cf_ml(i,j))*bqm(i,j,ntop(i,j)+1) +            &
                              cf_ml(i,j)*bqm_cld(i,j,ntop(i,j)+1) )
      wb_dzrad_int(i,j) = g * wb_dzrad_int(i,j)
      wb_dzrad_int(i,j) = max( zero, wb_dzrad_int(i,j) )

      if (model_type == mt_single_column) then
        do k = ntop(i,j)+1, ntdsc(i,j)+1
          wbend(i,j,k) = wb_dzrad_int(i,j)/dzrad(i,j)
        end do
      end if

    end if

    wb_dzrad_int(i,j) = max( rbl_eps, wb_dzrad_int(i,j) )

  end do ! I
end do ! J
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    l = i - pdims%i_start + 1 + pdims%i_len * (j - pdims%j_start)
    ind_todo(l) = l
    up(l)       = 1
    if (ktop_iterate(i,j)) then
      to_do(l) = .true.
    else
      to_do(l) = .false.
    end if
  end do
end do
!$OMP end do

!$OMP MASTER
c_len=pdims%i_len*pdims%j_len
!$OMP end MASTER
!$OMP BARRIER

do n_sweep = 1, num_sweeps_bflux
!$OMP BARRIER

!$OMP MASTER

          ! Compress to_do and ind_todo (will have new length c_len)
  call excfnl_cci(c_len, to_do, ind_todo)

      ! Restart inner interation with the points of outer
  c_len_i = c_len
  todo_inner(1:c_len_i) = to_do(1:c_len_i)
  ind_todo_i(1:c_len_i) = ind_todo(1:c_len_i)

!$OMP end MASTER
!$OMP BARRIER

  do ns = 1, n_steps
!$OMP BARRIER

!$OMP MASTER

              ! Calculate active elements and compress
    call excfnl_compin(up, wb_ratio, dec_thres, 2,                             &
                       c_len_i, ind_todo_i, todo_inner)
!$OMP end MASTER
!$OMP BARRIER

    !cdir nodep
!$OMP do SCHEDULE(STATIC)
    do ic = 1, c_len_i
      j1=(ind_todo_i(ic)-1)/pdims%i_end+1
      i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

      zdsc_base(i1,j1) = zdsc_base(i1,j1)+z_inc(i1,j1)
      scbase(i1,j1)    = .true.
          ! Flag for NBDSC found (only needed for LockWhelan2006)
      if (flux_grad  ==  LockWhelan2006) scbase(i1,j1) = .false.
      ft_nt_dscb(i1,j1)= zero
      fq_nt_dscb(i1,j1)= zero

      wbn_int(i1,j1) = zero
      wbp_int(i1,j1) = wb_dzrad_int(i1,j1)

      if (model_type == mt_single_column) then
        wbend(i1,j1,ksurf(i1,j1)) = zero
      end if

      if ( ksurf_iterate(i1,j1) .and.                                          &
           zdsc_base(i1,j1)  <   zsml_top(i1,j1) ) then
            ! only include surface flux if K_SURF is included
            ! in the wb calculation and K profiles overlap
        wbp_int(i1,j1) = wbp_int(i1,j1) + wb_surf_int(i1,j1)
        wbn_int(i1,j1) = zero

        if (model_type == mt_single_column) then
          wbend(i1,j1,ksurf(i1,j1)) = wb_surf_int(i1,j1)
        end if

      end if
      zinv_pr(i1,j1) = zhsc(i1,j1)-zdsc_base(i1,j1)

    end do ! ic c_len_i
!$OMP end do

    !..Integrate buoyancy flux profile given this ZDSC_BASE

!$OMP do SCHEDULE(STATIC)
    do jj = 1, c_len_i, bl_segment_size
      do k = 2, bl_levels-1
        !cdir nodep
        do ic = jj, min(jj+bl_segment_size-1, c_len_i)
          j1=(ind_todo_i(ic)-1)/pdims%i_end+1
          i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

          if ((k >= ksurf(i1,j1)+1) .and. (k <= ntop(i1,j1))) then

            khtop(i1,j1) = zero
            f2           = zero
            fsc          = zero
            if (.not. scbase(i1,j1) ) then
              ft_nt_dscb(i1,j1) = ft_nt(i1,j1,k)
              fq_nt_dscb(i1,j1) = fq_nt(i1,j1,k)
            end if
            z_pr = z_uv(i1,j1,k) - zdsc_base(i1,j1)

            if (z_pr >   zero .and. z_pr <  zinv_pr(i1,j1)) then

              if (.not. scbase(i1,j1) ) then
                scbase(i1,j1) = .true.
                z_ratio = (zdsc_base(i1,j1)-z_uv(i1,j1,k-1))                   &
                         /(z_uv(i1,j1,k)-z_uv(i1,j1,k-1))
                ft_nt_dscb(i1,j1) = ft_nt(i1,j1,k-1) +                         &
                      (ft_nt(i1,j1,k)-ft_nt(i1,j1,k-1))*z_ratio
                fq_nt_dscb(i1,j1) = fq_nt(i1,j1,k-1) +                         &
                      (fq_nt(i1,j1,k)-fq_nt(i1,j1,k-1))*z_ratio
              end if

              z_ratio = z_pr/zinv_pr(i1,j1)

              if (flux_grad  ==  LockWhelan2006) then
                khtop(i1,j1) = 3.6_r_bl * vkman * rho_mix(i1,j1,k)             &
                             * v_top_dsc(i1,j1) * zinv_pr(i1,j1)               &
                  * (z_ratio**3) * ( (one-z_ratio)*(one-z_ratio) )
                f2 = rho_mix(i1,j1,k) * one_half * z_ratio                     &
                                     * 2.0_r_bl**( z_ratio**4 )
                if ( v_sum_dsc(i1,j1)  >   zero ) then
                  fsc = 3.5_r_bl * rho_mix(i1,j1,k)                            &
                           * (v_top_dsc(i1,j1)/v_sum_dsc(i1,j1))               &
                           * (z_ratio**3) * (one-z_ratio)
                end if
              else
                ! max to avoid rounding errors giving small negative numbers
                khtop(i1,j1) = g1 * vkman * rho_mix(i1,j1,k)                   &
                    * v_top_dsc(i1,j1) * ((max(one - z_ratio,zero))**0.8_r_bl) &
                       * z_pr * z_ratio
              end if

            end if ! 0 < z_pr < zinv_pr

            khsurf(i1,j1) = zero
            if ( zdsc_base(i1,j1)  <   zsml_top(i1,j1) ) then
                ! only include K_surf if profiles overlap
                ! otherwise layers are independent
              khsurf(i1,j1) = kh_surf(i1,j1,k)
            end if

            if (flux_grad  ==  LockWhelan2006) then
              wslng = (f2+fsc)*(tothf_zi(i1,j1)-ft_nt_dscb(i1,j1))             &
                          - ( ft_nt(i1,j1,k)-ft_nt_dscb(i1,j1) )
              wqwng = (f2+fsc)*(totqf_zi(i1,j1)-fq_nt_dscb(i1,j1))             &
                          - ( fq_nt(i1,j1,k)-fq_nt_dscb(i1,j1) )
            end if

            if ( z_tq(i1,j1,k)  <=  z_cbase(i1,j1) ) then
                ! Completely below cloud-base so use cloud-free form
              wb_scld = khsurf(i1,j1)* db_ga_dry_n(i1,j1,k) +                  &
                        khtop(i1,j1) * db_noga_dry(i1,j1,k)
              if (flux_grad  ==  LockWhelan2006) then
                wb_scld = wb_scld + ( g/rdz(i1,j1,k) ) *                       &
                   ( btm(i1,j1,k-1)*wslng + bqm(i1,j1,k-1)*wqwng )
              end if
              wb_cld  = zero
            else if (z_tq(i1,j1,k-1)  >=  z_cbase(i1,j1)) then
                ! Completely above cloud-base so use cloudy formula
              wb_cld = ( khsurf(i1,j1) * db_ga_cld_n(i1,j1,k) +                &
                         khtop(i1,j1)  * db_noga_cld(i1,j1,k) )
              if (flux_grad  ==  LockWhelan2006) then
                wb_cld = wb_cld + ( g/rdz(i1,j1,k) ) * (                       &
                     ( btm(i1,j1,k-1)*(one-cf_ml(i1,j1)) +                     &
                       btm_cld(i1,j1,k-1)*cf_ml(i1,j1) )*wslng +               &
                     ( bqm(i1,j1,k-1)*(one-cf_ml(i1,j1)) +                     &
                       bqm_cld(i1,j1,k-1)*cf_ml(i1,j1) )*wqwng )
              end if
              wb_scld = zero
            else
                ! cloud-base within this integration range
                ! so treat cloud and sub-cloud layer wb separately
              wb_scld = khsurf(i1,j1) * db_ga_dry_n(i1,j1,k) +                 &
                        khtop(i1,j1) * db_noga_dry(i1,j1,k)
              wb_cld = ( khsurf(i1,j1) * db_ga_cld_n(i1,j1,k) +                &
                         khtop(i1,j1)  * db_noga_cld(i1,j1,k) )
              if (flux_grad  ==  LockWhelan2006) then
                wb_scld = wb_scld + ( g/rdz(i1,j1,k) ) *                       &
                   ( btm(i1,j1,k-1)*wslng + bqm(i1,j1,k-1)*wqwng )
                wb_cld = wb_cld + ( g/rdz(i1,j1,k) ) * (                       &
                     ( btm(i1,j1,k-1)*(one-cf_ml(i1,j1)) +                     &
                       btm_cld(i1,j1,k-1)*cf_ml(i1,j1) )*wslng +               &
                     ( bqm(i1,j1,k-1)*(one-cf_ml(i1,j1)) +                     &
                       bqm_cld(i1,j1,k-1)*cf_ml(i1,j1) )*wqwng )
              end if
              cld_frac = (z_tq(i1,j1,k)-z_cbase(i1,j1))                        &
                        /(z_tq(i1,j1,k)-z_tq(i1,j1,k-1))
              wb_cld  = cld_frac * wb_cld
              wb_scld = (one-cld_frac) * wb_scld
            end if

            if ( dzrad_disc_opt == dzrad_1p5dz .and. k == ntop(i1,j1) ) then
              ! At top of Sc layer, only include the part of the integral
              ! that is below the radiatively-cooled cloud-top layer.
              interp = ( zhsc(i1,j1) - dzrad(i1,j1) - z_tq(i1,j1,k-1) )        &
                     / ( z_tq(i1,j1,k) - z_tq(i1,j1,k-1) )
              wb_cld  = wb_cld  * interp
              wb_scld = wb_scld * interp
            end if

            if (wb_cld  >=  zero) then
              wbp_int(i1,j1) = wbp_int(i1,j1)+wb_cld
            else
              wbn_int(i1,j1) = wbn_int(i1,j1)-wb_cld
            end if
            if (wb_scld  >=  zero) then
              wbp_int(i1,j1) = wbp_int(i1,j1)+wb_scld
            else
              wbn_int(i1,j1) = wbn_int(i1,j1)-wb_scld
            end if

            if (model_type == mt_single_column) then
              wbend(i1,j1,k) = wb_cld + wb_scld
            end if

          end if ! K
        end do ! ic c_len_i
      end do ! K
    end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
    do ic = 1, c_len_i
      j1=(ind_todo_i(ic)-1)/pdims%i_end+1
      i1=ind_todo_i(ic)-(j1-1)*pdims%i_end
      wb_ratio(i1,j1)=wbn_int(i1,j1)/wbp_int(i1,j1)
    end do ! ic c_len_i
!$OMP end do

  end do  ! loop stepping up through ML

  !cdir nodep
!$OMP do SCHEDULE(STATIC)
  do ic = 1, c_len
    l=ind_todo(ic)
    j1=(l-1)/pdims%i_end+1
    i1=l-(j1-1)*pdims%i_end

    !..sub-divide current Z_INC into one more part than there will be steps
    !..as there is no need to recalculate WB at a current Z_INC
    z_inc(i1,j1)= z_inc(i1,j1)/real(n_steps+1, r_bl)

    if (                                                                       &
       (up(l) == 1 .and. wb_ratio(i1,j1) <= dec_thres(i1,j1)) .or.             &
              ! hit thres while working up
       (up(l) == 0 .and. wb_ratio(i1,j1) >= dec_thres(i1,j1))) then
              ! hit thres while working down
      up(l) = 1-up(l)   ! change direction of sweep
      z_inc(i1,j1)=- z_inc(i1,j1)
    else
      if ( l_use_sml_dsc_fixes ) then
        ! Set base height to limit when within 1m of limit
        ! (consistent with what is done for zsml_top)
        if ( zdsc_base(i1,j1) >= z_top_lim(i1,j1)-one ) then
               ! hit upper height limit (give-or-take 1m) without
               ! reaching threshold
          to_do(ic)=.false.
          zdsc_base(i1,j1) = z_top_lim(i1,j1)
        else if (zdsc_base(i1,j1) <=  z_bot_lim(i1,j1)+one ) then
               ! hit lower height limit (give-or-take 1m) without
               ! reaching threshold
          to_do(ic)=.false.
          zdsc_base(i1,j1) = z_bot_lim(i1,j1)
        end if
      else  ! ( .not. l_use_sml_dsc_fixes )
        ! Original code leaves base height close to but not at limit
        if (  zdsc_base(i1,j1) >= z_top_lim(i1,j1)-one .or.                    &
              zdsc_base(i1,j1) <=  z_bot_lim(i1,j1)+one ) then
          ! hit height limits (give-or-take 1m) without
          ! reaching threshold
          to_do(ic)=.false.
        end if
      end if  ! ( l_use_sml_dsc_fixes )
    end if

    !..Note that if the threshold has not been passed then the next sweep
    !..continues in the same direction (but with reduced increment).

  end do ! c_len
!$OMP end do

end do  ! loop over sweeps

! Convert integrated WB to profiles of WB itself for diagnostics
if (model_type == mt_single_column) then

!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels-1
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        ! Include sml
        if ( .not. zdsc_base(i,j) < zsml_top(i,j) ) then
          ! Note: wbend already includes the SML buoyancy flux
          ! in the case where the SML and DSC mixing profiles overlap
          ! (since the DSC b-flux integration sums both the SML and DSC
          !  mixing contributions and stores this in wbend).
          ! So only add the SML contribution back on if they do not overlap.
          wbend(i,j,k) = wbend(i,j,k) + wbend_sml(i,j,k)
        end if
      end do
    end do
  end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  ! Note parallelised over j as k isn't independent
  do j = pdims%j_start, pdims%j_end
    do k = 1, bl_levels-1
      do i = pdims%i_start, pdims%i_end
        ! convert to m2/s-3
        if ( k <= ntop(i,j) ) then
          if ( k <= ksurf(i,j) ) then
            gamma_wbs = ( (wbmix(i,j,ksurf(i,j))/z_tq(i,j,ksurf(i,j)))         &
                      - bflux_surf(i,j)  )*2.0_r_bl/z_tq(i,j,ksurf(i,j))
            wbmix(i,j,k) = bflux_surf(i,j) + gamma_wbs*z_uv(i,j,k)

            gamma_wbs = ( (wbend(i,j,ksurf(i,j))/z_tq(i,j,ksurf(i,j)))         &
                      - bflux_surf(i,j)  )*2.0_r_bl/z_tq(i,j,ksurf(i,j))
            wbend(i,j,k) =  bflux_surf(i,j) + gamma_wbs*z_uv(i,j,k)
          else
            wbmix(i,j,k)=wbmix(i,j,k)*rdz(i,j,k)
            wbend(i,j,k)=wbend(i,j,k)*rdz(i,j,k)
          end if
        end if
      end do
    end do
  end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        ! convert to m2/s-3
        if ( k >=  ntop(i,j)+1 .and. k <= ntdsc(i,j)+1 ) then
          wbend(i,j,k) = wb_dzrad_int(i,j)/dzrad(i,j)
        end if
        if ( k >=  ntop(i,j)+1 .and. k <= ntml(i,j)+1 ) then
          wbmix(i,j,k) = wb_dzrad_int(i,j)/dzrad(i,j)
        end if
      end do
    end do
  end do
!$OMP end do
end if ! model_type

! ----------------------------------------------------------------------
! 2.4 Set depth of cloud-top driven mixing in SML when there is a DSC
!     layer above (eg. fog under Sc) to be the SML layer depth
!     and, if option selected, the top of any K profile above the LCL
!     in cumulus to be the lower of the parcel top, the DSC base and
!     500m above the LCL but at least 1.1*z_lcl
! ----------------------------------------------------------------------
if ( entr_smooth_dec == on .or. entr_smooth_dec == entr_taper_zh ) then

  !cdir collapse
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( cumulus(i,j) .or.                                                   &
           ( dsc(i,j) .and. zdsc_base(i,j) < zh(i,j) ) ) then
        ! ignore SML `cloud-top' driven mixing
        zsml_base(i,j) = zh(i,j)
        v_top(i,j)     = zero
      else
        zsml_base(i,j) = 0.1_r_bl*zh(i,j)
      end if
    end do  ! loop over j
  end do  ! loop over I
!$OMP end do

else ! entr_smooth_dec off

  !cdir collapse
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( cumulus(i,j) .or. coupled(i,j) ) then
        ! ignore SML `cloud-top' driven mixing
        zsml_base(i,j) = zh(i,j)
        v_top(i,j)     = zero
      else
        zsml_base(i,j) = 0.1_r_bl*zh(i,j)
      end if
    end do  ! loop over j
  end do  ! loop over I
!$OMP end do

end if  ! test on entr_smooth_dec

!cdir collapse
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    zsml_top(i,j) = zh(i,j)
    if (bl_res_inv /= off .and. dzh(i,j) > zero)                               &
                                zsml_top(i,j) = zh(i,j)+dzh(i,j)
  end do  ! loop over j
end do  ! loop over I
!$OMP end do

if ( kprof_cu == klcl_entr ) then
  !cdir collapse
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( cumulus(i,j) ) then
        zsml_top(i,j) = zh(i,j)+max_cu_depth
        ! Keep top of Kprof for cu below DSC base:
        if (dsc(i,j)) then
          zsml_top(i,j) = min( zsml_top(i,j), zdsc_base(i,j) )
        end if
        ! ...but at least 1.1*z_lcl:
        zsml_top(i,j) = max( 1.1_r_bl*z_lcl(i,j), zsml_top(i,j) )
        ! ...but no higher than parcel top:
        zsml_top(i,j) = min( zhpar(i,j), zsml_top(i,j) )

        cu_depth_scale(i,j) = (zsml_top(i,j)-zh(i,j))/3.0_r_bl

        rhokh_lcl(i,j) = zero
        ! Use the BL entrainment parametrization as calculated above
        rhokh_lcl(i,j) = min( rhokh_surf_ent(i,j), 5.0_r_bl)
      end if
    end do  ! loop over j
  end do  ! loop over I
!$OMP end do
end if  ! test on kprof_cu
!-----------------------------------------------------------------------
! 3.  Calculate factors required to ensure that the non-local turbulent
!     mixing coefficient profiles are continuous as the entrainment
!     level is approached.
!-----------------------------------------------------------------------
!cdir collapse

!$OMP do SCHEDULE(DYNAMIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    k=ntml(i,j)+1
    kh_top_factor(i,j) = max( 0.7_r_bl , one - sqrt(                           &
             rhokh_surf_ent(i,j) /                                             &
                   ( rho_mix(i,j,k)*w_h_top(i,j)*vkman*zh(i,j) ) ) )
    if (l_use_var_fixes) then
      km_top_factor(i,j) = max( 0.7_r_bl, one-sqrt( rhokm_surf_ent(i,j) /      &
               ( rho_wet_tq(i,j,k-1)*w_m_top(i,j)*vkman*zh(i,j) ) ) )
    else
      km_top_factor(i,j) = max( 0.7_r_bl , one - sqrt( rhokm(i,j,k) /          &
               ( rho_wet_tq(i,j,k-1)*w_m_top(i,j)*vkman*zh(i,j) ) ) )
    end if
    scdepth(i,j) = zh(i,j) - zsml_base(i,j)
    factor = g1 * rho_mix(i,j,k) * v_top(i,j) *vkman *scdepth(i,j)
    if ( factor  >   zero) then
      kh_sct_factor(i,j) = one - ( rhokh_top_ent(i,j) / factor )**1.25_r_bl
                                                        ! 1.25=1/0.8
    else
      kh_sct_factor(i,j) = one
    end if
    factor = g1 * rho_wet_tq(i,j,k-1) * v_top(i,j) *                           &
                    vkman * scdepth(i,j) * 0.75_r_bl
    if ( factor  >   zero) then
      if (l_use_var_fixes) then
        km_sct_factor(i,j) = one -                                             &
             ( rhokm_top_ent(i,j) / factor )**1.25_r_bl
                                                        ! 1.25=1/0.8
      else
        km_sct_factor(i,j) = one -                                             &
             ( rhokm_top(i,j,k) / factor )**1.25_r_bl
      end if
    else
      km_sct_factor(i,j) = one
    end if

    if (ntdsc(i,j)  >   0) then
      !-------------------------------------------------------------
      ! Set up factors to ensure K profile continuity at ZHSC;
      ! no need to limit size of factor as precise shape of top-down
      ! mixing profile not important.
      ! Only calculate _DSCT_FACTORs when a decoupled stratocumulus
      ! layer exists, i.e. NTDSC > 0.
      !-------------------------------------------------------------
      k=ntdsc(i,j)+1
      dscdepth(i,j) = zhsc(i,j) - zdsc_base(i,j)
      factor = g1*rho_mix(i,j,k)*v_top_dsc(i,j)*vkman*dscdepth(i,j)
      if ( factor  >   zero) then
        kh_dsct_factor(i,j) = one -                                            &
                            ( rhokh_dsct_ent(i,j) / factor )**1.25_r_bl
                                                        ! 1.25=1/0.8
      else
        kh_dsct_factor(i,j) = one
      end if

      factor = 0.75_r_bl * g1 * rho_wet_tq(i,j,k-1) * v_top_dsc(i,j) *         &
                           vkman * dscdepth(i,j)
      if ( factor  >   zero) then
        if (l_use_var_fixes) then
          km_dsct_factor(i,j) = one -                                          &
               ( rhokm_dsct_ent(i,j) / factor )**1.25_r_bl
                                                        ! 1.25=1/0.8
        else
          km_dsct_factor(i,j) = one -                                          &
               ( rhokm_top(i,j,k) / factor )**1.25_r_bl
        end if
      else
        km_dsct_factor(i,j) = one
      end if
    end if
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
! 4.  Calculate height dependent turbulent
!     transport coefficients within the mixing layer.
!-----------------------------------------------------------------------

! Reset identifiers of base of decoupled layer mixing

!$OMP MASTER

!cdir collapse
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    scbase(i,j) = .false.
    nbdsc(i,j)  = 0
  end do
end do

!$OMP end MASTER
!$OMP BARRIER

!-------------------------------------------------------------
! Calculate RHOK(H/M)_TOP, top-down turbulent mixing profiles
! for the surface mixed layer.
! This is a variation on an up-side-down version of the cubic
! surface-forced profiles below.  Implement between at least
! the top of the `surface layer' (at Z=0.1*ZH) and ZH.
! Note this may well include NTML+1: entrainment fluxes will
! be dealt with in KMKHZ.
!-------------------------------------------------------------
c_tke = 1.33_r_bl/(vkman*g1)

!$OMP do SCHEDULE(STATIC)
do jj = pdims%j_start, pdims%j_end, bl_segment_size
  do k = 2, bl_levels
    !cdir collapse
    do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
      do i = pdims%i_start, pdims%i_end

        !           Calculate the height of u,v-level above the surface
        ! *APL: z0m removed from z in K(z)
        zk_uv = z_uv(i,j,k)

        !           Calculate the height of T,q-level above the surface

        zk_tq = z_tq(i,j,k-1)

        if ( zk_uv  <   zh(i,j) .and.                                          &
             zk_uv  >   zsml_base(i,j) ) then
          z_pr  = zk_uv - zsml_base(i,j)
          zh_pr = zh(i,j) - zsml_base(i,j)
          z_ratio = z_pr/zh_pr
          if (flux_grad  ==  LockWhelan2006) then

            rhokh_top(i,j,k) = 3.6_r_bl*vkman * rho_mix(i,j,k) * v_top(i,j)    &
                               * zh_pr * (z_ratio**3)                          &
                               * (( one - z_ratio )**2)

            if ( .not. coupled(i,j) ) then
              rhof2(i,j,k)  = rho_mix(i,j,k) * one_half * z_ratio              &
                                          * 2.0_r_bl**( z_ratio**4 )
              if ( v_sum(i,j)  >   zero ) then
                rhofsc(i,j,k) = 3.5_r_bl * rho_mix(i,j,k)                      &
                                * (v_top(i,j)/v_sum(i,j))                      &
                                * (z_ratio**3) * (one-z_ratio)
              end if
            end if

          else  ! Not LockWhelan2006
            ! max to avoid rounding errors giving small negative numbers
            rhokh_top(i,j,k) = rho_mix(i,j,k) * v_top(i,j) * g1 *              &
              vkman * ((max(one - kh_sct_factor(i,j)*z_ratio,zero))**0.8_r_bl )&
                                             * z_pr * z_ratio
          end if

        end if
        !-------------------------------------------------------
        !   For LockWhelan2006, KM_TOP could be changed to match
        !   the shape of KH_TOP.  This has not been done on the
        !   grounds that the change in shape arises with the
        !   inclusion of the other non-gradient terms.
        !-------------------------------------------------------
        if ( zk_tq  <   zh(i,j) .and.                                          &
             zk_tq  >   zsml_base(i,j) ) then
          z_pr = zk_tq - zsml_base(i,j)
          zh_pr = zh(i,j) - zsml_base(i,j)
          rhokm_top(i,j,k) = 0.75_r_bl * rho_wet_tq(i,j,k-1) * v_top(i,j) *    &
                g1 * vkman *                                                   &
                ( (max(one - km_sct_factor(i,j)*z_pr/zh_pr, zero))**0.8_r_bl ) &
                                           * z_pr * z_pr / zh_pr
                                                    ! PRANDTL=0.75
          if (BL_diag%l_tke) then
            ! save Km/timescale for TKE diag, completed in bdy_expl2
            tke_nl(i,j,k)=rhokm_top(i,j,k)*c_tke*v_top(i,j)/zh(i,j)
          end if
        end if
        !-------------------------------------------------------------
        ! Add contribution to top-down mixing coefficient
        ! profiles for decoupled stratocumulus layers when
        ! one exists
        !-------------------------------------------------------------
        if ( zk_uv  <   zhsc(i,j) .and.                                        &
                zk_uv  >   zdsc_base(i,j) ) then
          if (.not. scbase(i,j) ) then
            scbase(i,j) = .true.
            ! identifies lowest layer below which there is mixing
            nbdsc(i,j) = k
          end if
          !-----------------------------------------------------------
          ! Calculate RHOK(H/M)_TOP, top-down turbulent mixing
          ! profiles and add to any generated in the surface mixing
          ! layer.
          ! This is a variation on an up-side-down version of the
          ! cubic surface-forced profiles above.  Implement between
          ! at least the top of the `surface layer' (at Z=0.1*ZH) and
          ! ZHSC.
          !-----------------------------------------------------------
          z_pr = zk_uv - zdsc_base(i,j)
          zh_pr = zhsc(i,j) - zdsc_base(i,j)
          z_ratio = z_pr/zh_pr

          if (flux_grad  ==  LockWhelan2006) then

            rhokh_top(i,j,k) = 3.6_r_bl*vkman * rho_mix(i,j,k)                 &
                             * v_top_dsc(i,j) * zh_pr * (z_ratio**3)           &
                             * (( one - z_ratio )**2)

            rhof2(i,j,k)  = rho_mix(i,j,k) * one_half * z_ratio                &
                                          * 2.0_r_bl**( z_ratio**4 )
            if ( v_sum_dsc(i,j)  >   zero ) then
              rhofsc(i,j,k) = 3.5_r_bl * rho_mix(i,j,k)                        &
                                * (v_top_dsc(i,j)/v_sum_dsc(i,j))              &
                                * (z_ratio**3) * (one-z_ratio)
            end if

          else  ! Not LockWhelan2006
            ! max to avoid rounding errors giving small negative numbers
            rhokh_top(i,j,k) = rhokh_top(i,j,k) +                              &
               rho_mix(i,j,k)*v_top_dsc(i,j)*g1*vkman*                         &
                  ( (max(one - kh_dsct_factor(i,j)*z_ratio,zero))**0.8_r_bl )  &
                                             * z_pr * z_ratio
          end if
        end if
        !-------------------------------------------------------------
        ! Now momentum
        !-------------------------------------------------------------
        if ( zk_tq  <   zhsc(i,j) .and.                                        &
             zk_tq  >   zdsc_base(i,j) ) then
          z_pr = zk_tq - zdsc_base(i,j)
          zh_pr = zhsc(i,j) - zdsc_base(i,j)
          ! max to avoid rounding errors giving small negative numbers
          rhokm_dsct =                                                         &
             0.75_r_bl*rho_wet_tq(i,j,k-1)*v_top_dsc(i,j)*g1*vkman*            &
                ( (max(one - km_dsct_factor(i,j)*z_pr/zh_pr,zero))**0.8_r_bl ) &
                                        * z_pr * z_pr / zh_pr
          if (BL_diag%l_tke) then
            ! save Km/timescale for TKE diag, completed in bdy_expl2
            tke_nl(i,j,k) = tke_nl(i,j,k) +                                    &
                            rhokm_dsct*c_tke*v_top_dsc(i,j)/dscdepth(i,j)
          end if
          rhokm_top(i,j,k) = rhokm_top(i,j,k) + rhokm_dsct
        end if

      end do
    end do
  end do
end do
!$OMP end do

!----------------------------------------------------
! Now K_SURF profiles
!----------------------------------------------------
if (flux_grad  ==  LockWhelan2006) then
  !----------------------------------------------------
  ! Lock and Whelan formulation
  !----------------------------------------------------

  c_ws = 0.42_r_bl     ! ~ PR_NEUT^3 by design
  pr_neut = 0.75_r_bl
  pr_conv = 0.6_r_bl
  c_tke = 1.33_r_bl/(vkman*c_ws**two_thirds)

!$OMP do SCHEDULE(STATIC)
  do jj = pdims%j_start, pdims%j_end, bl_segment_size
    do k = 2, bl_levels
      !cdir collapse
      do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
        do i = pdims%i_start, pdims%i_end

          zk_uv = z_uv(i,j,k)
          zk_tq = z_tq(i,j,k-1)

          if (fb_surf(i,j)  >=  zero) then

            ! Calculate the free-convective scaling velocity at z(k)

            if (coupled(i,j)) then  !  coupled
              if ( entr_smooth_dec == entr_taper_zh ) then
                wstar3 = (      svl_diff_frac(i,j)  * zhsc(i,j)                &
                         + (one-svl_diff_frac(i,j)) * zh(i,j)   ) *fb_surf(i,j)
              else
                wstar3 = zhsc(i,j) * fb_surf(i,j)
              end if
            else
              wstar3 = zh(i,j) * fb_surf(i,j)
            end if

            if (zk_uv  <=  0.1_r_bl*zh(i,j)) then
              !             Surface layer calculation
              w_s_cubed_uv = 10.0_r_bl*c_ws * zk_uv * fb_surf(i,j)
            else
              !             Outer layer calculation
              w_s_cubed_uv = c_ws * wstar3
            end if

            if (zk_tq  <=  0.1_r_bl*zh(i,j)) then
              !             Surface layer calculation
              w_s_cubed_tq = 10.0_r_bl*c_ws * zk_tq * fb_surf(i,j)
            else
              !             Outer layer calculation
              w_s_cubed_tq = c_ws * wstar3
            end if

            !           Turbulent velocity scale for scalars

            w_m_neut = ( v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_uv )           &
                                   **one_third
            w_h_uv = w_m_neut/pr_neut

            ! Also calc on TQ levels for W_M_TQ
            w_m_neut = ( v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_tq )           &
                                   **one_third
            w_h_tq = w_m_neut/pr_neut

            ! The calculations below involve v_s**4. In very rare circumstances
            ! v_s can be order(10^-11), therefore v_s**4 is order(10^-44)
            ! which is outside the range of single precision calculations
            ! Hence we create a double precision version to enforce the correct
            ! calculation
            v_s_dbl(i,j) = v_s(i,j)

            ! Turbulent Prandtl number and velocity scale for scalars

            Prandtl = pr_neut*                                                 &
            ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +            &
             (one/(c_ws*25.0_r_bl))*w_s_cubed_tq*w_m_neut ) /                  &
            ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +            &
             (one/(c_ws*25.0_r_bl))*(pr_neut/pr_conv)*w_s_cubed_tq*            &
             w_m_neut )

            w_m_tq = Prandtl * w_h_tq

            if ( zk_uv  <   zh(i,j) ) then
              !---------------------------------------------------------
              ! Calculate RHOKH(w_h,z/z_h)
              !---------------------------------------------------------

              rhokh(i,j,k) = rho_mix(i,j,k) * w_h_uv * vkman * zk_uv *         &
                                    ( one - ( zk_uv / zh(i,j) ) ) *            &
                                    ( one - ( zk_uv / zh(i,j) ) )

            end if
            if ( zk_tq  <   zh(i,j) ) then
              !---------------------------------------------------------
              ! Calculate RHOKM(w_m,z/z_h)
              !---------------------------------------------------------

              rhokm(i,j,k) = rho_wet_tq(i,j,k-1)*w_m_tq*vkman*zk_tq *          &
                                    ( one - ( zk_tq / zh(i,j) ) ) *            &
                                    ( one - ( zk_tq / zh(i,j) ) )

              if (BL_diag%l_tke) then
                ! save Km/timescale for TKE diag, completed in bdy_expl2
                tke_nl(i,j,k) = tke_nl(i,j,k) +                                &
                                   rhokm(i,j,k)*c_tke*w_m_tq/zh(i,j)
              end if
            end if
          end if
        end do
      end do
    end do
  end do
!$OMP end do

else
  !----------------------------------------------------
  ! Lock et al and Holtstalg and Boville formulations
  !----------------------------------------------------

        ! Default to Lock et al
  c_ws = 0.25_r_bl
  pr_neut = 0.75_r_bl
  pr_conv = 0.375_r_bl
  if (flux_grad  ==  HoltBov1993) then
    c_ws = 0.6_r_bl
    pr_neut = one
    pr_conv = 0.6_r_bl
  end if
  c_tke = 1.33_r_bl/(vkman*c_ws**two_thirds)

!$OMP do SCHEDULE(STATIC)
  do jj = pdims%j_start, pdims%j_end, bl_segment_size
    do k = 2, bl_levels
      !cdir collapse
      do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
        do i = pdims%i_start, pdims%i_end

          !         Calculate the height of u,v-level above the surface
          ! *APL: z0m removed from z in K(z)
          zk_uv = z_uv(i,j,k)

          !         Calculate the height of T,q-level above the surface

          zk_tq = z_tq(i,j,k-1)

          if (fb_surf(i,j)  >=  zero) then

            ! Calculate the free-convective scaling velocity at z(k)

            if (zk_uv  <=  0.1_r_bl*zh(i,j)) then

              !             Surface layer calculation

              w_s_cubed_uv = 10.0_r_bl*c_ws * zk_uv * fb_surf(i,j)
            else

              !             Outer layer calculation

              if (coupled(i,j)) then  !  coupled and cloudy
                if ( entr_smooth_dec == entr_taper_zh ) then
                  w_s_cubed_uv = c_ws                                          &
                     * (      svl_diff_frac(i,j)  * zhsc(i,j)                  &
                       + (one-svl_diff_frac(i,j)) * zh(i,j)   ) * fb_surf(i,j)
                else
                  w_s_cubed_uv = c_ws * zhsc(i,j) * fb_surf(i,j)
                end if
              else
                w_s_cubed_uv = c_ws * zh(i,j) * fb_surf(i,j)
              end if
            end if

            if (zk_tq  <=  0.1_r_bl*zh(i,j)) then

              !             Surface layer calculation

              w_s_cubed_tq = 10.0_r_bl*c_ws * zk_tq * fb_surf(i,j)
            else

              !             Outer layer calculation

              if (coupled(i,j)) then  !  coupled and cloudy
                if ( entr_smooth_dec == entr_taper_zh ) then
                  w_s_cubed_tq = c_ws                                          &
                     * (      svl_diff_frac(i,j)  * zhsc(i,j)                  &
                       + (one-svl_diff_frac(i,j)) * zh(i,j)   ) * fb_surf(i,j)
                else
                  w_s_cubed_tq = c_ws * zhsc(i,j) * fb_surf(i,j)
                end if
              else
                w_s_cubed_tq = c_ws * zh(i,j) * fb_surf(i,j)
              end if
            end if

            !           Turbulent velocity scale for momentum

            w_m_uv = (v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_uv)               &
                                   **one_third

            w_m_tq = (v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_tq)               &
                                   **one_third

            ! The calculations below involve v_s**4. In very rare circumstances
            ! v_s can be order(10^-11), therefore v_s**4 is order(10^-44)
            ! which is outside the range of single precision calculations
            ! Hence we create a double precision version to enforce the correct
            ! calculation
            v_s_dbl(i,j) = v_s(i,j)

            !           Turbulent Prandtl number and velocity scale for scalars

            Prandtl = pr_neut*                                                 &
               ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +         &
               (one/(c_ws*25.0_r_bl))*w_s_cubed_uv*w_m_uv ) /                  &
               ( v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j)*v_s_dbl(i,j) +         &
               (one/(c_ws*25.0_r_bl))*(pr_neut/pr_conv)*w_s_cubed_uv*          &
               w_m_uv )
            w_h_uv = w_m_uv / Prandtl

            if ( zk_uv  <   zh(i,j) ) then
              !---------------------------------------------------------
              ! Calculate RHOKH(w_h,z/z_h)
              !---------------------------------------------------------

              rhokh(i,j,k) = rho_mix(i,j,k) * w_h_uv * vkman * zk_uv *         &
                  ( one - kh_top_factor(i,j) * ( zk_uv / zh(i,j) ) ) *         &
                  ( one - kh_top_factor(i,j) * ( zk_uv / zh(i,j) ) )
            else if ( kprof_cu == klcl_entr .and. cumulus(i,j) ) then
              if ( zk_uv < zsml_top(i,j) ) then
                ! Exponential decay from ZH but tends to zero
                !  at zsml_top
                rhokh(i,j,k) = rhokh_lcl(i,j) *                                &
                       exp(-(zk_uv-zh(i,j))/cu_depth_scale(i,j)) *             &
                       (one-(zk_uv-zh(i,j))/(zsml_top(i,j)-zh(i,j)))
              end if
            end if
            if ( zk_tq  <   zh(i,j) ) then
              !---------------------------------------------------------
              ! Calculate RHOKM(w_m,z/z_h)
              !---------------------------------------------------------

              rhokm(i,j,k) = rho_wet_tq(i,j,k-1)*w_m_tq*vkman* zk_tq *         &
                  ( one - km_top_factor(i,j) * ( zk_tq / zh(i,j) ) ) *         &
                  ( one - km_top_factor(i,j) * ( zk_tq / zh(i,j) ) )

              if (BL_diag%l_tke) then
                ! save Km/timescale for TKE diag, completed in bdy_expl2
                tke_nl(i,j,k) = tke_nl(i,j,k) +                                &
                                   rhokm(i,j,k)*c_tke*w_m_tq/zh(i,j)
              end if
            else if ( kprof_cu == klcl_entr .and. cumulus(i,j) ) then
              if ( zk_tq < zsml_top(i,j) ) then
                ! Exponential decay from ZH but tends to zero
                !  at zsml_top
                rhokm(i,j,k) = prandtl_top(i,j) * rhokh_lcl(i,j) *             &
                       exp(-(zk_tq-zh(i,j))/cu_depth_scale(i,j)) *             &
                       (one-(zk_tq-zh(i,j))/(zsml_top(i,j)-zh(i,j)))
                if (BL_diag%l_tke) then
                  ! save Km/timescale for TKE diag, completed in bdy_expl2
                  tke_nl(i,j,k) = tke_nl(i,j,k) +                              &
                                     rhokm(i,j,k)*c_tke*w_m_tq/zh(i,j)
                end if
              end if
            end if

          end if
        end do
      end do
    end do
  end do
!$OMP end do

end if  ! Test on Flux_grad

if ( ng_stress  ==  BrownGrant97 .or.                                          &
     ng_stress  ==  BrownGrant97_limited .or.                                  &
     ng_stress  ==  BrownGrant97_original ) then

!$OMP do SCHEDULE(STATIC)
  do jj = pdims%j_start, pdims%j_end, bl_segment_size
    do k = 2, bl_levels
      !cdir collapse
      do j = jj, min(jj+bl_segment_size-1, pdims%j_end)
        do i = pdims%i_start, pdims%i_end
          zk_tq = z_tq(i,j,k-1)   ! stresses are calc on theta-levs
          if ( fb_surf(i,j) > zero .and. zk_tq < zh(i,j) ) then
            !---------------------------------------------------------
            ! Calculate non-gradient stress function
            ! (Brown and Grant 1997)
            ! Shape function chosen such that non-gradient stress
            ! either goes to zero at 0.1*ZH and ZH (ng_stress==BrownGrant97
            ! and ng_stress==BrownGrant97_limited), or goes to zero at
            ! at the ground surface and ZH (ng_stress==BrownGrant97_original)
            !---------------------------------------------------------
            ng_stress_calculate = .false.

            if (ng_stress == BrownGrant97 .or.                                 &
                ng_stress == BrownGrant97_limited) then
              if ( zk_tq  >   0.1_r_bl*zh(i,j) ) then
                z_pr = zk_tq - 0.1_r_bl*zh(i,j)
                zh_pr = 0.9_r_bl*zh(i,j)
                ng_stress_calculate = .true.
              end if
            else if (ng_stress==BrownGrant97_original) then
              z_pr = zk_tq
              zh_pr = zh(i,j)
              ng_stress_calculate = .true.
            end if

            if ( ng_stress_calculate ) then
              ! Outer layer calculation
              if (coupled(i,j)) then  !  coupled and cloudy
                if ( entr_smooth_dec == entr_taper_zh ) then
                  wstar3 = (      svl_diff_frac(i,j)  * zhsc(i,j)              &
                           + (one-svl_diff_frac(i,j)) * zh(i,j) ) *fb_surf(i,j)
                else
                  wstar3 = zhsc(i,j) * fb_surf(i,j)
                end if
              else
                wstar3 =   zh(i,j) * fb_surf(i,j)
              end if

              ! Use the Holtslag and Boville velocity scale for
              ! non-gradient stress stability dependence, as in BG97
              w_m_hb_3 = v_s(i,j)*v_s(i,j)*v_s(i,j) + 0.6_r_bl*wstar3
              f_ngstress(i,j,k) =(rho_wet_tq(i,j,k-1)/rhostar_gb(i,j))         &
                * s_m * ( a_ngs * wstar3 / w_m_hb_3 )                          &
                   * ( z_pr / zh_pr ) * ( one -  ( z_pr / zh_pr ) ) *          &
                                        ( one -  ( z_pr / zh_pr ) )
            end if  ! ng_stress_calculate
          end if  ! fb_surf>0 and z<zh
        end do
      end do
    end do
  end do
!$OMP end do
end if

!$OMP end PARALLEL

if (l_wtrac) deallocate(iset_wtrac)

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
!------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine excf_nl_9c
end module excf_nl_9c_mod
