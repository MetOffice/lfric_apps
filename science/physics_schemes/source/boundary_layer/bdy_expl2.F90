! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Calculate the explicit turbulent fluxes of heat, moisture
!           and momentum between atmospheric levels
!           within the boundary layer, and/or the effects of these
!           fluxes on the primary model variables.

!  Programming standard : UMDP 3

!  Documentation: UMDP 24.

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module bdy_expl2_mod

use fm_drag_mod, only: fm_drag

use tuning_segments_mod, only: bl_segment_size

use um_types, only: r_bl

use timer_mod, only: timer

implicit none

!Automatic segment size tuning

character(len=*), parameter, private :: ModuleName = 'BDY_EXPL2_MOD'
contains

subroutine bdy_expl2 (                                                         &
! in values defining vertical grid of model atmosphere :
 bl_levels,p_theta_levels,land_pts,land_index,                                 &
 r_theta_levels, r_rho_levels, cos_theta_latitude,                             &
! in U, V and W momentum fields.
 u_p,v_p,u_0_px,v_0_px,                                                        &
! in from other part of explicit boundary layer code
 rho_mix,rho_wet_tq,rho_mix_tq,dzl_charney,rdz,rdz_charney_grid,               &
 z_tq,z_uv,rhostar,bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt,         &
 recip_l_mo_sea,flandg,rib_gb, sil_orog_land, z0m_eff_gb,                      &
! in cloud/moisture data :
 cf_bulk,q,qcf,qcl,t,qw,tl,                                                    &
! in everything not covered so far :
 rad_hr,micro_tends,fb_surf,u_s,pstar,tstar,h_blend_orog,                      &
 zh_prev,zhpar,z_lcl,ho2r2_orog,sd_orog,wtrac_as,                              &
! 2 in 3 INOUT for Smagorinsky
 delta_smag, max_diff, rneutml_sq, visc_m, visc_h,                             &
! SCM Diagnostics (dummy values in full UM) & stash diagnostics
 nSCMDpkgs,L_SCMDiags,BL_diag,                                                 &
! INOUT variables
 zh,dzh,ntml,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,w,etadot,             &
 t1_sd,q1_sd,wtrac_bl,                                                         &
! out new variables for message passing
 tau_fd_x, tau_fd_y, f_ngstress,                                               &
! out Diagnostic not requiring STASH flags :
 zhnl,shallowc,cu_over_orog,                                                   &
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,bl_type_7,        &
! Out data for turbulent generation of mixed-phase cloud:
 bl_w_var,                                                                     &
! out data required for tracer mixing :
 kent, we_lim, t_frac, zrzi, kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,       &
! out data required elsewhere in UM system :
 zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0,rhcpt,                                   &
 tgrad_bm, mix_len_bm                                                          &
 )

use atm_fields_bounds_mod, only: pdims, tdims, wdims, tdims_l, pdims_s


use bl_diags_mod, only: strnewbldiag
use bl_option_mod, only:                                                       &
    off, max_t_grad, a_grad_adj, sg_orog_mixing, l_use_surf_in_ri,             &
    h_scale, t_drain, idyndiag, DynDiag_ZL, l_noice_in_turb,                   &
    DynDiag_ZL_corrn, DynDiag_ZL_CuOnly, DynDiag_Ribased,                      &
    RiCrit_sharp, zhloc_depth_fac, non_local_bl, on, l_full_lambdas,           &
    nl_bl_levels, local_fa, free_trop_layers, to_sharp_across_1km,             &
    sbl_op, equilibrium_sbl, one_third, two_thirds, blending_option,           &
    blend_except_cu, blend_cth_shcu_only, sg_shear, sg_shear_enh_lambda,       &
    max_tke, var_diags_opt, original_vars, split_tke_and_inv, tke_diag_fac,    &
    i_interp_local, i_interp_local_gradients, i_interp_local_cf_dbdz,          &
    shallow_cu_maxtop, sc_cftol, near_neut_z_on_l, zero, one, one_half
use cloud_inputs_mod, only: i_rhcpt, forced_cu, i_cld_vn, i_pc2_init_method,   &
    ez_max_bm, i_bm_ez_opt, i_bm_ez_subcrit, i_bm_ez_entpar
use cderived_mod, only: delta_lambda, delta_phi
use cv_run_mod, only: l_param_conv, l_wvar_for_conv
use gen_phys_inputs_mod, only: l_mr_physics
use jules_surface_mod, only: formdrag, explicit_stress
use missing_data_mod, only: rmdi
use mixing_config_mod, only: smag_l_calc, smag_l_calc_use_geo
use model_domain_mod, only: model_type, mt_single_column
use mphys_inputs_mod, only: l_subgrid_qcl_mp
use pc2_constants_mod, only: rhcpt_tke_based, i_cld_bimodal, i_cld_pc2,        &
                             pc2init_bimodal, bm_negative_init, bm_tiny
use planet_constants_mod, only: vkman => vkman_bl, grcp => grcp_bl,            &
     pref => pref_bl, kappa => kappa_bl, g => g_bl
use s_scmop_mod,   only: default_streams,                                      &
                         t_inst, t_avg, d_bl, d_sl, d_point, scmdiag_bl
use science_fixes_mod, only: l_fix_dyndiag, l_fix_zh
use stochastic_physics_run_mod, only: l_rp2, par_mezcla_rp, rp_idx,            &
    i_rp_scheme, i_rp2b, cs_rp, zhloc_depth_fac_rp
use turb_diff_mod, only:                                                       &
    l_subfilter_vert, l_subfilter_horiz, mix_factor,                           &
    l_blend_isotropic, turb_startlev_vert, turb_endlev_vert

use qsat_mod, only: qsat, qsat_mix, qsat_wat, qsat_wat_mix

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

!$ use omp_lib, only: omp_get_max_threads

use btq_int_mod, only: btq_int
use ex_coef_mod, only: ex_coef
use ex_flux_tq_mod, only: ex_flux_tq
use kmkh_mod, only: kmkh
use kmkhz_9c_mod, only: kmkhz_9c

use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use wtrac_bl_mod,            only: bl_wtrac_type
use wtrac_atm_step_mod,      only: atm_step_wtrac_type

implicit none

!  Inputs :-
integer, intent(in) ::                                                         &
 land_pts,                                                                     &
                             ! No.of land points in whole grid.
 bl_levels
                             ! in Max. no. of "boundary" levels

!     Declaration of new BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

real(kind=r_bl), intent(in) ::                                                 &
  p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
                 0:bl_levels+1),                                               &
  r_theta_levels(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,  &
                 0:bl_levels),                                                 &
  r_rho_levels(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,    &
               bl_levels),                                                     &
                                 ! in height of rho and theta levels
 cos_theta_latitude(pdims_s%i_start:pdims_s%i_end,                             &
                    pdims_s%j_start:pdims_s%j_end),                            &
                                 ! in cosine of the latitude
 rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels+1),                                                         &
                                 ! in density on UV (ie. rho) levels;
                                 !    used in RHOKH so dry density if
                                 !    L_mr_physics is true
 rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                                 ! in density on TQ (ie. theta) levels;
                                 !    used in RHOKM so wet density
 rho_mix_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                                 ! in density on TQ (ie. theta) levels;
                                 !    used in non-turb flux integration
                                 !    so dry density if L_mr_physics is true
 dzl_charney(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             bl_levels),                                                       &
                                 ! in DZL(,K) is depth in m of theta
                                 !    level K, i.e. distance from
                                 !    boundary K-1/2 to boundary K+1/2
 rdz( pdims_s%i_start:pdims_s%i_end,                                           &
      pdims_s%j_start:pdims_s%j_end, bl_levels ),                              &
                                 ! in RDZ(,1) is the reciprocal of
                                 !    the height of level 1, i.e. of
                                 !    the middle of layer 1.  For
                                 !    K > 1, RDZ(,K) is the
                                 !    reciprocal of the vertical
                                 !    distance from level K-1 to
                                 !    level K.
 rdz_charney_grid(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  bl_levels),                                                  &
                                 ! in RDZ(,1) is the reciprocal of
                                 !       the height of level 1,
                                 !       i.e. of the middle of layer 1
                                 !       For K > 1, RDZ(,K) is the
                                 !       reciprocal of the vertical
                                 !       distance from level K-1 to
                                 !       level K.
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                                 ! in Z_tq(*,K) is height of full
                                 !    level k.
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1),        &
                                  ! out Z_uv(*,K) is height of half
                                  ! level k-1/2.
 rhostar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                                 ! in Surface air density
 u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),           &
                                 ! in U on P-grid.
 v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),           &
                                 ! in V on P-grid.
 bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                                 ! in A buoyancy parameter for clear
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                                 ! in A buoyancy parameter for clear
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bt_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                                 ! in A buoyancy parameter for cloudy
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bq_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                                 ! in A buoyancy parameter for cloudy
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                                 ! in A grid-box mean buoyancy param
                                 ! on p,T,q-levels (full levels).
 bq_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                                 ! in A grid-box mean buoyancy param
                                 ! on p,T,q-levels (full levels).
 a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                                 ! in Saturated lapse rate factor
                                 !    on p,T,q-levels (full levels).
 a_dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                                 ! in Saturated lapse rate factor
                                 !    on p,T,q-levels (full levels).
 dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                 ! in Derivative of q_SAT w.r.t. T

real(kind=r_bl), intent(in) ::                                                 &
 recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                 ! in Reciprocal of the surface
                                 !    Obukhov length over sea (m^-1).
 flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
                                 ! in Land fraction on all tiles

! (f) Atmospheric + any other data not covered so far, incl control.
real(kind=r_bl), intent(in) ::                                                 &
 rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2,bl_levels),                                                          &
                                  ! in (LW,SW) rad heating rate (K/s)
  micro_tends(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
              2, bl_levels),                                                   &
                         ! Tendencies from microphys within BL levels
                         ! (TL, K/s; QW, kg/kg/s)
 fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
                                  ! in Surface flux buoyancy over
                                  ! density (m^2/s^3)

 u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                     &
                                  ! in Surface friction velocity
                                  !    (m/s)
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                                  ! in Surface pressure (Pascals).
 tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                                  ! in Surface temperature (K).
 h_blend_orog(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                  ! in Blending height used as part
                                  ! of effective roughness scheme
 zh_prev(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
                                  ! in boundary layer height from
                                  !    previous timestep
 rib_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
                               ! in  Bulk Richardson number for lowest
                               ! layer
 sil_orog_land(land_pts),                                                      &
                               ! in Silhouette area of unresolved
                               ! orography per unit horizontal area
 zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                               ! in Height of top of initial
                               !     parcel ascent
 z_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! in Height of LCL

! Water tracer structure containing 'micro_tends' fields
type(atm_step_wtrac_type), intent(in) :: wtrac_as(n_wtrac)

! Additional variables for SCM diagnostics which are dummy in full UM
integer, intent(in) ::                                                         &
 nSCMDpkgs             ! No of SCM diagnostics packages

logical, intent(in) ::                                                         &
 L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

real(kind=r_bl), intent(in) ::                                                 &
 u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),          &
                                 ! in W'ly component of surface
!                                       current (m/s). P grid
   v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
                                   ! in S'ly component of surface
!                                       current (m/s). P grid
   ho2r2_orog(land_pts),                                                       &
                                   ! in peak to trough height of
!                                       unresolved orography
!                                       on land points only (m)
   sd_orog(land_pts),                                                          &
                                   ! in Standard Deviation of unresolved
!                                       orography on land points only (m)
   z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! in Effective grid-box roughness
!                                 length for momentum

integer, intent(in) ::                                                         &
 land_index(land_pts)        ! in LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.
! (e) Cloud data.
real(kind=r_bl), intent(in) ::                                                 &
 cf_bulk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                                        ! in Cloud fraction (decimal).
 qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                                   ! in Cloud ice (kg per kg air)
 qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                                   ! in Cloud liquid water
 q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,                &
   tdims_l%k_start:bl_levels),                                                 &
                                   ! in specific humidity
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),             &
                                   ! in temperature
 qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),           &
                                 ! in Total water content
 tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),           &
                                 ! in Ice/liquid water temperature
 delta_smag(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                                 ! in delta_x used by Smagorinsky
 max_diff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                 ! in maximum diffusion coefficient

! INOUT variables
real(kind=r_bl), intent(in out) ::                                             &
 zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                      &
                                 ! INOUT Height above surface of top
                                 !       of boundary layer (metres).
 dzh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                     &
                                 ! INOUT inversion thickness (m)
 fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                                 ! INOUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                                   ! INOUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),       &
                                   ! INOUT Exchange coeffs for moisture.
    w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,0:bl_levels),        &
    etadot(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,                &
           0:bl_levels),                                                       &
   t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
                                    ! INOUT Standard deviation of
                                    ! turbulent fluctuations of layer 1
                                    ! temperature; for use in
                                    ! initiating convection.
   q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                    ! INOUT Standard deviation of turbulent
                                    !    fluctuations of layer 1
                                    !    humidity; for use in initiating
                                    !    convection.

real(kind=r_bl), intent(in out) ::                                             &
 rhokm(pdims_s%i_start:pdims_s%i_end,                                          &
       pdims_s%j_start:pdims_s%j_end ,bl_levels),                              &
!            Exchange coefficients for momentum on P-grid
 rneutml_sq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
! Square of the neutral mixing length for Smagorinsky
 visc_m(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,bl_levels),&
! Smagorinsky diffusion coefficient for momentum
 visc_h(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,bl_levels)
! Smagorinsky diffusion coefficient for heat

logical, intent(in out) ::                                                     &
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                                 ! INOUT Logical switch for trade Cu
 l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! INOUT Flag to indicate shallow
                                 !     convection

integer, intent(in out) ::                                                     &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                               ! INOUT Number of model layers in the
                               !    turbulently mixed layer
 ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! INOUT Top level of initial parcel
                               !  ascent. Used in convection scheme.

! Water tracer structure containing boundary layer fields
type(bl_wtrac_type), intent(in out) :: wtrac_bl(n_wtrac)

!  Outputs :-
!  (a) Calculated anyway (use STASH space from higher level) :-
real(kind=r_bl), intent(out) ::                                                &
 f_ngstress(pdims_s%i_start:pdims_s%i_end,                                     &
            pdims_s%j_start:pdims_s%j_end,2:bl_levels),                        &
 tau_fd_x(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,         &
          bl_levels),                                                          &
 tau_fd_y(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,         &
          bl_levels),                                                          &
  bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! out Indicator set to 1.0 if stable
                                 !     b.l. diagnosed, 0.0 otherwise.
  bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! out Indicator set to 1.0 if Sc over
                                 !     stable surface layer diagnosed,
                                 !     0.0 otherwise.
  bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! out Indicator set to 1.0 if well
                                 !     mixed b.l. diagnosed,
                                 !     0.0 otherwise.
  bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! out Indicator set to 1.0 if
                                 !     decoupled Sc layer (not over
                                 !     cumulus) diagnosed,
                                 !     0.0 otherwise.
  bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! out Indicator set to 1.0 if
                                 !     decoupled Sc layer over cumulus
                                 !     diagnosed, 0.0 otherwise.
  bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! out Indicator set to 1.0 if a
                                 !     cumulus capped b.l. diagnosed,
                                 !     0.0 otherwise.
  bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! out Indicator set to 1.0 if a
                                 !     Shear-dominated unstable b.l.
                                 !     diagnosed, 0.0 otherwise.

real(kind=r_bl), intent(out) ::                                                &
                     bl_w_var(   tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             2 : tdims%k_end+1 )

real(kind=r_bl), intent(out) ::                                                &
  zhnl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                                 ! out non-local PBL depth
  wstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                                 ! out Convective velocity scale (m/s)
  wthvs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                                 ! out surface flux of thv (Km/s)
  shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                                 ! out Shallow Cu diagnostic
                                 !   Indicator set to 1.0 if shallow,
                                 !   0.0 if not shallow or not cumulus
  cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! out Indicator for cumulus
                                 !     over steep orography
                                 !   Indicator set to 1.0 if true,
                                 !   0.0 if false. Exclusive.
  we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),               &
                                  ! out rho*entrainment rate implied b
                                  !     placing of subsidence
  zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),                 &
                                  ! out (z-z_base)/(z_i-z_base)
  t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),               &
                                  ! out a fraction of the timestep
  we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),           &
                                  ! out rho*entrainment rate implied b
                                  !     placing of subsidence
  zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),             &
                                  ! out (z-z_base)/(z_i-z_base)
  t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),           &
                                  ! out a fraction of the timestep
  zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! out Top of decoupled layer

integer, intent(out) ::                                                        &
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                                 ! out Top level for turb mixing in
                                 !     any decoupled Sc layer
 nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                                 ! out Bottom level of any decoupled
                                 !     turbulently-mixed Sc layer.
 kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                                 ! out grid-level of SML inversion
 kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! out grid-level of DSC inversion

!-2 Genuinely output, needed by other atmospheric routines :-
real(kind=r_bl), intent(out) ::                                                &
  uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                           ! out U-component of surface wind stress
                           !     on P-grid
  vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                           ! out V-component of surface wind stress
                           !     on P-grid
real(kind=r_bl), intent(out) ::                                                &
                     rhcpt(tdims%i_start:tdims%i_end,                          &
                           tdims%j_start:tdims%j_end,1:tdims%k_end)
real(kind=r_bl), intent(out) :: tgrad_bm(tdims%i_start:tdims%i_end,            &
                                        tdims%j_start:tdims%j_end,1:tdims%k_end)
!       Gradient of liquid potential temperature interpolated to theta levels  &
!       for bimodal cloud scheme
real(kind=r_bl), intent(out) :: mix_len_bm                                     &
                                       (tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end)
!       Turbulent mixing length for bimodal cloud scheme

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

character(len=*), parameter ::  RoutineName = 'BDY_EXPL2'


! Parameters also passed to EX_COEF
! Layer interface K_LOG_LAYR-1/2 is the highest which requires log
! profile correction factors to the vertical finite differences.
! The value should be reassessed if the vertical resolution is changed.
! We could set K_LOG_LAYR = BL_LEVELS and thus apply the correction
! factors for all the interfaces treated by the boundary layer scheme;
! this would be desirable theoretically but expensive computationally
! because of the use of the log function.
integer ::    k_log_layr
parameter (k_log_layr=2)
!-----------------------------------------------------------------------
!  Workspace :-
real(kind=r_bl) ::                                                             &
 recip_time_sbl, recip_time_cbl ! inverse timescales for TKE diagnostic
real(kind=r_bl) ::                                                             &
 a_dqsdtm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          bl_levels),                                                          &
                              ! Saturated lapse rate factor
                              ! on intermediate levels (half levels).
 a_qsm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                              ! Saturated lapse rate factor
                              ! on intermediate levels (half levels).
 bqm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                              ! A buoyancy parameter for clear air
                              ! on intermediate levels (half levels).
 bqm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                              ! A buoyancy parameter for cloudy air
                              ! on intermediate levels (half levels).
 btm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                              ! A buoyancy parameter for clear air
                              ! on intermediate levels (half levels).
 btm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                              ! A buoyancy parameter for cloudy air
                              ! on intermediate levels (half levels).
 dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),        &
                              ! Buoyancy gradient across layer
                              !  interface.
 dbdz_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                              ! Buoyancy gradient across layer
                              !  interface, inc gradient adjustment
 dvdzm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
       2:bl_levels),                                                           &
                              ! Modulus of wind shear.
 rmlmax2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, bl_levels),      &
                              ! Square of asymptotic mixing length
                              ! for Smagorinsky scheme
 ri(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),          &
                              ! Local Richardson number.
 ri_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       2:bl_levels),                                                           &
                              ! Local Richardson number, inc grad adj
 grad_q_adj(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                              ! Humidity gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
   grad_t_adj(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                ! Temperature gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
   rhokhz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                 &
          2:bl_levels),                                                        &
                                ! Non-local turbulent mixing
!                                 coefficient for heat and moisture.
   rhokh_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
             2:bl_levels),                                                     &
                                ! Non-local turbulent mixing coefficient
                                ! for top-down mixing of heat and
                                ! moisture.
   rhokh_th(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,               &
            bl_levels),                                                        &
                                ! local scheme rhokh on th-levels,
                                ! index k held on th-level(k-1),
                                ! same as rhokm
   rhokmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          2:bl_levels),                                                        &
                                ! Non-local turbulent mixing
!                                 coefficient for momentum.
   tke_loc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels),                                                       &
                                ! Ri-based scheme diagnosed TKE
   tke_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                 &
             2:bl_levels),                                                     &
                                ! Non-local scheme TKE diag times rho
   rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             2:bl_levels),                                                     &
                                ! Non-local turbulent mixing coefficient
                                ! for top-down mixing of momentum.
   weight_1dbl(pdims%i_start:pdims%i_end,                                      &
               pdims%j_start:pdims%j_end ,bl_levels),                          &
                                ! Weighting applied to 1D BL scheme
                                ! to blend with Smagorinsky scheme,
                                ! index k held on theta level (k-1)
   weight_1dbl_rho(pdims%i_start:pdims%i_end,                                  &
                   pdims%j_start:pdims%j_end,bl_levels),                       &
                                ! weight_1dbl interpolated to rho levels
   elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),       &
                                ! Mixing length for momentum
   elh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),       &
   elh_rho(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           2:bl_levels),                                                       &
                                ! Mixing length for heat (m),
                                ! held on theta and rho levels, resp.
   fm_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),       &
                                ! stability function for momentum transport
                                ! level 1 value is dummy
   fh_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),       &
                                ! stability function for heat and moisture.
                                ! level 1 value is dummy
   sigma_h(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
                                ! Standard deviation of subgrid
                                ! orography for sg mixing options (m)
   dbdz_rh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,1:bl_levels),   &
                                ! Grid-mean static stability on rho-levels
   dbdz_ga_rh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
              1:bl_levels),                                                    &
                                ! Gradient-adjusted dbdz on rho-levels
   supersat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            1:bl_levels),                                                      &
                                ! Supersaturation
                                ! (qw - qsat(Tl))/(1 + Lc/cp dqsat/dT)
                                ! gradient of this is used for estimating
                                ! saturated fraction on rho-levels
   qs_tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,1:bl_levels)
                                ! qsat(Tl) on current level
                                ! Vectorised over levels (k) for open-mp ease

real(kind=r_bl) ::                                                             &
   frac_sat, frac_dry, frac_edg, frac_lev, qc_tot, bt_rh, bq_rh
   ! Temporary variables used to compute buoyancy coefficients on rho-levels

real(kind=r_bl), allocatable :: visc_h_rho (:,:,:)
                                ! visc_h on rho levels

    ! Terms for non-gradient flux parametrization
    !  (=0 unless using 9C code with FLUX_GRAD=LockWhelan2006)
real(kind=r_bl) ::                                                             &
  ft_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels+1),                                                          &
                              ! Non-turbulent heat and moisture flux
  fq_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels+1)          !  (on rho levels, surface flux(K=1)=0)
real(kind=r_bl) ::                                                             &
  rhof2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        2:bl_levels),                                                          &
                              ! f2 and fsc term shape profiles
  rhofsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         2:bl_levels)

real(kind=r_bl) ::                                                             &
  tothf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                              ! Total heat fluxes at inversions
  tothf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &

  totqf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                              ! Total moisture fluxes at inversions
  totqf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &

  ft_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                              ! Non-turbulent heat and moisture flux
  fq_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              !    at the base of the DSC layer.

real(kind=r_bl) ::                                                             &
zh_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! Height above surface of top of
                              !  boundary layer (metres) as
                              !  determined from the local
                              !  Richardson number profile.
dsldz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                     &
      bl_levels),                                                              &
                              ! TL+gz/cp gradient between
                              ! levels K and K-1
dsldz_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         bl_levels),                                                           &
                              ! TL+gz/cp gradient between
                              ! levels K and K-1, inc gradient adjust
dqwdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                              ! QW gradient between levels K and K-1


integer ::                                                                     &
 ntml_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                                 ! Number of model layers in the
                                 ! turbulently mixed layer as
                                 ! determined from the local
                                 ! Richardson number profile.
 ntml_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                                 ! Number of model layers in the
                                 ! turbulently mixed layer as
                                 ! determined from the parcel ascent.
 ntml_save(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                                 ! saved copy of ntml on entry
 sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                 ! Flags for whether discontinuous
 dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                 ! inversions are diagnosed
 kplume(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! Start grid-level for surface-driven plume

logical ::                                                                     &
 unstable(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                               ! Logical switch for unstable
                               !    surface layer.
 dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                     &
                               ! Flag set if decoupled
                               ! stratocumulus layer found
 coupled(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                               ! Flag to indicate Sc layer weakly
                               ! coupled to surface (ie weakly
                               ! decoupled)
 dynamic_bl_diag(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                               ! Flag to indicate the dynamic
                               ! diagnosis (iDynDiag) has
                               ! determined the BL type
 topbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                               ! Flag for having reached
                               ! the top of the turbulently mixed
                               ! layer.
 l_shallow_cth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                               ! Flag to indicate shallow convection based on
                               ! cf_bulk cloud top height
 cloud_base_found(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                               ! Flag for having reached cloud base
 cloud_top_found(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! Flag for having reached cloud top

!  Local scalars :-
real(kind=r_bl) ::                                                             &
 dzu,                                                                          &
            ! Westerly wind shear between levels K+1 and K.
 dzv,                                                                          &
            ! Southerly wind shear between levels K+1 and K.
 lambda_min,                                                                   &
            ! Min value of length scale LAMBDA.
 lambdah,                                                                      &
            ! Asymptotic mixing length for turbulent transport
            ! of heat/moisture.
 vkz,                                                                          &
            ! Temporary in calculation of ELH.
 f_log,                                                                        &
            ! Temporary in calculation of logarithmic correction
 zmaxb_for_dsc,                                                                &
 zmaxt_for_dsc
            ! Max heights to look for DSC cloud base and top

real(kind=r_bl) ::                                                             &
  weight1,                                                                     &
  weight2,                                                                     &
  weight3,                                                                     &
  z_scale,                                                                     &
             ! scaling with height
  zpr,                                                                         &
             ! z/sigma_h
  slope,                                                                       &
             ! subgrid orographic slope
  dsldzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
             ! TL+gz/cp gradient interpolated to Z_TQ
  dsldzm_ga,                                                                   &
             ! dsldzm with gradient adjustment terms
  dqwdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
             ! QW gradient interpolated to Z_TQ
  qssurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
             ! qsat of surface

! variables for rhcrit parametrization
real(kind=r_bl) ::                                                             &
 b2, sh, exner, root6, delta_x, var_fac, sl_var, qw_var, sl_qw,                &
 sgm(tdims%i_start:tdims%i_end),                                               &
 qsw_arr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,1:bl_levels-1),   &
 max_rhcpt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
 min_rhcpt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

integer ::                                                                     &
 i,j,                                                                          &
                ! LOCAL Loop counter (horizontal field index).
 k,km,kp,ient,                                                                 &
                ! LOCAL Loop counter (vertical level index).
 kmax,                                                                         &
                ! level of max rhokm,
 l,                                                                            &
                ! LOCAL Loop counter for land points
 ntop,                                                                         &
                ! NTPAR restricted below BL_LEVELS
 i_wt           ! Water tracer loop counter


integer ::                                                                     &
 pdims_omp_block, pdims_seg_block,                                                         &
                  ! for pdims open-mp blocking
 tdims_omp_block, tdims_seg_block,                                                         &
                  ! for tdims open mp blocking
 ii,                                                                           &
                  ! for indexing over open-mp block
 max_threads,                                                                  &
                  ! Hold the number of threads for open-mp
 seg_slice_start, seg_slice_end,                                               &
                  ! Segmentaion of calls into qstat, start and end
 bl_segment_range ! Range of segmentaion, derived from start and end

real(kind=r_bl), parameter :: max_abs_obkhov = 1.0e6_r_bl
                 ! Maximum permitted magnitude of the Obukhov
                 ! length (m).
real(kind=r_bl), parameter :: small_sh = 0.01_r_bl
                 ! Minimum  value of sh used in
                 ! original_vars variance diagnostics
real(kind=r_bl), parameter :: small_tke = 1.0e-6_r_bl
                 ! Minimum required value of TKE before
                 ! split_tke_and_inv variance diagnostics are calculated
real(kind=r_bl), parameter :: max_ri = 0.01_r_bl*sqrt(huge(one))
                 ! Maximum (absolute) Richardson number which ensures that
                 ! the stability functions (~ri^2) remain real-valued at
                 ! the given model precision

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
call timer('bdy_expl2')

! Set up blocking for OMP sections
max_threads = 1
!$ max_threads = omp_get_max_threads()

!-----------------------------------------------------------------------
! 1) Set various diagnostics and switches
!-----------------------------------------------------------------------
! checking of nl_bl_levels has been moved to readsize/scm_shell
! so that it is only executed at initialisation

!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( pdims, dynamic_bl_diag, ntml_save, ntml, bl_levels,             &
!$OMP          weight_1dbl, weight_1dbl_rho, BL_diag, u_s, fb_surf )           &
!$OMP  private( i, j, k )
!$OMP  do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    dynamic_bl_diag(i,j) = .false.
    ntml_save(i,j) = ntml(i,j)
  end do
end do
!$OMP end do NOWAIT

!------------------------------------------------------------------
!  Initialize weighting applied to 1d BL scheme
!  (used to blend with 3D Smagorinsky scheme)
!------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      weight_1dbl(i,j,k) = one
      weight_1dbl_rho(i,j,k) = one
    end do
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! Set surface scaling diagnostics
!-----------------------------------------------------------------------
 ! Obukhov length
if (BL_diag%l_oblen) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      !       Limit the magnitude of the Obukhov length to avoid
      !       problems with packing.
      BL_diag%oblen(i,j)= u_s(i,j)*u_s(i,j)*u_s(i,j)
      if ( BL_diag%oblen(i,j) <                                                &
           max_abs_obkhov*abs(vkman*fb_surf(i,j)) ) then
        BL_diag%oblen(i,j)=-BL_diag%oblen(i,j)/(vkman*fb_surf(i,j))
      else
        BL_diag%oblen(i,j)=-sign(max_abs_obkhov, fb_surf(i,j))
      end if
    end do
  end do
!$OMP end do NOWAIT
end if

! Ustar
if (BL_diag%l_ustar) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      BL_diag%ustar(i,j)=u_s(i,j)
    end do
  end do
!$OMP end do NOWAIT
end if

! Surface buoyancy flux
if (BL_diag%l_wbsurf) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      BL_diag%wbsurf(i,j)=fb_surf(i,j)
    end do
  end do
!$OMP end do
end if
!$OMP end PARALLEL
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.  Interpolate BT and BQ to half levels and calculate Ri
!-----------------------------------------------------------------------
call btq_int (                                                                 &
! in levels
   bl_levels,                                                                  &
! in fields
   z_tq,z_uv,bq,bt,bq_cld,bt_cld,a_qs,a_dqsdt,                                 &
! out fields
   bqm,btm,bqm_cld,btm_cld,a_qsm,a_dqsdtm                                      &
    )
!-----------------------------------------------------------------------
! Calculate lapse rates
!-----------------------------------------------------------------------
! Sub divide pdims domain len over the number of threads
pdims_omp_block = ceiling(real(pdims%i_len)/max_threads)
! Pick the smallest option. 
! Work should be split accross the threads, this provides a cap.
! In the occurence where the pdims domain len is very small, it will be seleted.
! Otherwise bl_segment_size will be selected.
pdims_seg_block = min(bl_segment_size, pdims_omp_block, pdims%i_len)
!$OMP  PARALLEL DEFAULT(SHARED) private(ii, i, j, k, weight1,  weight2,        &
!$OMP  weight3, zpr, dzv, dzu, l, slope, dsldzm_ga,                            &
!$OMP  qs_tl, frac_sat, frac_dry, frac_edg, frac_lev, qc_tot, bt_rh, bq_rh,    &
!$OMP  seg_slice_start, seg_slice_end, bl_segment_range)
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    grad_t_adj(i,j) = min( max_t_grad,                                         &
                           a_grad_adj * t1_sd(i,j) / zh_prev(i,j) )
  end do
end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      dsldz(i,j,k)    = ( tl(i,j,k) - tl(i,j,k-1) )                            &
                              *rdz_charney_grid(i,j,k) + grcp
      dsldz_ga(i,j,k) = dsldz(i,j,k)
      if ( z_tq(i,j,k) <= zh_prev(i,j) ) then
        dsldz_ga(i,j,k) = dsldz_ga(i,j,k) - grad_t_adj(i,j)
      end if
      dqwdz(i,j,k)    = ( qw(i,j,k) - qw(i,j,k-1) )                            &
                             * rdz_charney_grid(i,j,k)
    end do
  end do
end do
!$OMP end do

! We need a grid-box mean surface-to-level-1 buoyancy gradient in order
! to construct Ri on level 1.  Either the gradient from level 1 to 2
! can be extrapolated, or surface properties can be used.  Over a
! heterogeneous land surface this is poorly defined and we can't use Rib
! from the surface scheme as vertically averaging Ri is numerically
! unstable.  So, over land, only the average temperature gradient is used
if ( .not. l_use_surf_in_ri .and. (var_diags_opt > original_vars .or.          &
                             i_interp_local == i_interp_local_cf_dbdz ) ) then
  ! if not using surface variables in Ri (l_use_surf_in_ri=false) we
  ! extrapolate dbdz itself from level 2, but the sl and qw gradients are
  ! used in the newer variance calculations and with i_interp_local_cf_dbdz
  ! so extrapolate them here
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      dsldz(i,j,1)    = dsldz(i,j,2)
      dsldz_ga(i,j,1) = dsldz_ga(i,j,2)
      dqwdz(i,j,1)    = dqwdz(i,j,2)
    end do
  end do
!$OMP end do

else ! l_use_surf_in_ri = true
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do ii = pdims%i_start, pdims%i_end, pdims_seg_block
      seg_slice_start  = ii
      seg_slice_end = MIN(((ii+pdims_seg_block)-1), pdims%i_end)
      bl_segment_range = (seg_slice_end - seg_slice_start) + 1

      if (l_noice_in_turb) then
        ! use qsat_wat
        if ( l_mr_physics ) then
          call qsat_wat_mix(qssurf(seg_slice_start:seg_slice_end,j), &
                            tstar(seg_slice_start:seg_slice_end,j),  &
                            pstar(seg_slice_start:seg_slice_end,j),  &
                            bl_segment_range) ! No halos
        else
          call qsat_wat(qssurf(seg_slice_start:seg_slice_end,j),  &
                        tstar(seg_slice_start:seg_slice_end,j),   &              
                        pstar(seg_slice_start:seg_slice_end,j),   &
                        bl_segment_range) ! No halos
        end if
      else ! l_noice_in_turb
        ! use qsat
        if ( l_mr_physics ) then
          call qsat_mix(qssurf(seg_slice_start:seg_slice_end,j),  &
                        tstar(seg_slice_start:seg_slice_end,j),   &
                        pstar(seg_slice_start:seg_slice_end,j),   &
                        bl_segment_range) ! No halos
        else
          call qsat(qssurf(seg_slice_start:seg_slice_end,j),  &
                    tstar(seg_slice_start:seg_slice_end,j),   &
                    pstar(seg_slice_start:seg_slice_end,j),   &
                    bl_segment_range) ! No halos
        end if
      end if ! l_noice_in_turb
    end do ! ii
  end do ! j
!$OMP end do NOWAIT
  k=1
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      dsldz(i,j,k)    = ( tl(i,j,k) - tstar(i,j) )                             &
                              *rdz_charney_grid(i,j,k) + grcp
      dsldz_ga(i,j,k) = dsldz(i,j,k) ! no GA below level 1
      if ( flandg(i,j) < 0.2_r_bl) then
        dqwdz(i,j,k)  = ( qw(i,j,k) - qssurf(i,j) )                            &
                             * rdz_charney_grid(i,j,k)
      else
        dqwdz(i,j,k)  = dqwdz(i,j,2) ! extrapolate qw if mainly land
      end if
    end do
  end do
!$OMP end do

end if

! Output gradient adjusted SL

! Local Ri-based calculation of RHOKM and RHOKH:

select case(i_interp_local)
case (i_interp_local_gradients)
  ! Interpolate gradients to theta-levels, then
  ! calculate `buoyancy' gradient, DBDZ, on theta-levels
  ! NOTE: DBDZ(K) is on theta-level K-1

!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        weight1 = r_rho_levels(i,j,k) -                                        &
                  r_rho_levels(i,j,k-1)
        weight2 = r_theta_levels(i,j,k-1)-                                     &
                  r_rho_levels(i,j,k-1)
        weight3 = r_rho_levels(i,j,k) -                                        &
                  r_theta_levels(i,j,k-1)
        dsldzm(i,j,k) = weight2 * dsldz(i,j,k)                                 &
               + weight3 * dsldz(i,j,k-1)
        dqwdzm(i,j,k) = weight2 * dqwdz(i,j,k)                                 &
               + weight3 * dqwdz(i,j,k-1)
        dbdz(i,j,k) = g*( bt_gb(i,j,k-1)*dsldzm(i,j,k) +                       &
                          bq_gb(i,j,k-1)*dqwdzm(i,j,k) )/weight1
        !       ! Now with gradient adjustment
        dsldzm_ga = weight2 * dsldz_ga(i,j,k)                                  &
               + weight3 * dsldz_ga(i,j,k-1)
        dbdz_ga(i,j,k) = g*( bt_gb(i,j,k-1)*dsldzm_ga +                        &
                             bq_gb(i,j,k-1)*dqwdzm(i,j,k) )/weight1
        ! Complete calculation of interpolated gradients
        dsldzm(i,j,k) = dsldzm(i,j,k) / weight1
        dqwdzm(i,j,k) = dqwdzm(i,j,k) / weight1
      end do
    end do
  end do
!$OMP end do

  ! Overwrite element 2 with level 1 to 2 gradients
  ! (ie instead of using centred value that uses surface parameters)
  if ( .not. l_use_surf_in_ri ) then
    k = 2
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        dbdz(i,j,k) = g*( bt_gb(i,j,k-1)*dsldz(i,j,k) +                        &
                          bq_gb(i,j,k-1)*dqwdz(i,j,k) )
        dbdz_ga(i,j,k) = g*( bt_gb(i,j,k-1)*dsldz_ga(i,j,k) +                  &
                             bq_gb(i,j,k-1)*dqwdz(i,j,k) )
      end do
    end do
!$OMP end do
  end if

case (i_interp_local_cf_dbdz)
  ! Alternative ordering of calculation of dbdz;
  ! Compute dbdz on rho-levels, to avoid having to use interpolated
  ! sl,qw gradients, and then interpolate dbdz onto theta-levels.
  ! Requires calculation of buoyancy coefficients on rho-levels,
  ! with some complexity around how to estimate the cloud-fraction
  ! on rho-levels (depends on gradient of total-water supersaturation
  ! qw - qsat(Tl) between the neighbouring theta-levels).
  ! Crucially this avoids a feature of the previous option, where
  ! applying the cloudy buoyancy coefficients at a cloud top level k
  ! to a strong qw gradient (between k-1 and k+1) can yield an unstable
  ! dbdz despite strong static stability.

  ! Calculate total-water supersaturation qw - qsat(Tl), on theta-levels
!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do ii = tdims%i_start, tdims%i_end, bl_segment_size
      seg_slice_start  = ii
      seg_slice_end = MIN(((ii+bl_segment_size)-1), tdims%i_end)
      bl_segment_range = (seg_slice_end - seg_slice_start) + 1

      ! Calculate qsat(Tl)...
      if ( l_mr_physics ) then
        call qsat_mix(qs_tl(seg_slice_start:seg_slice_end,:,k),                &
                      tl(seg_slice_start:seg_slice_end,:,k),                   &
                      p_theta_levels(seg_slice_start:seg_slice_end,:,k),       &
                      bl_segment_range,                                        &
                      tdims%j_len )
      else
        call qsat(qs_tl(seg_slice_start:seg_slice_end,:,k),                    &
                  tl(seg_slice_start:seg_slice_end,:,k),                       &
                  p_theta_levels(seg_slice_start:seg_slice_end,:,k),           &
                  bl_segment_range,                                            &
                  tdims%j_len )
      end if
    end do ! ii
    ! ...then subtract from qw to get supersaturation, and multiply by
    !    1/(1 + Lc/cp dqsat/dT)
    !    in order to compare with values of qcl+qcf.
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        supersat(i,j,k) = a_qs(i,j,k) * ( qw(i,j,k) - qs_tl(i,j,k) )
      end do !j
    end do !i
  end do !k
!$OMP end do

  ! Calculate dbdz on rho-levels, so between th-levels k-1 and k
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

        ! Interpolate the buoyancy coefficients onto rho-levels...
        ! Split the interpolation into weighted contributions from:
        ! a) Area that is cloudy at k-1 and at k
        !    (just use saturated coefficients).
        ! b) Area that is cloudy at one height but not at the other
        !    (interpolate qw-qsat to find fraction of level that is saturated).
        ! c) Area that is cloud-free at k-1 and k
        !    (just use cloud-free coefficients).

        ! Area (a) is min of cloud-fractions from both levels (max overlap)
        frac_sat = min( cf_bulk(i,j,k-1), cf_bulk(i,j,k) )
        ! Area (c) is 1 - max of both levels
        frac_dry = one - max( cf_bulk(i,j,k-1), cf_bulk(i,j,k) )
        ! Area (b) (cloud-edge between k-1 and k) is the remainder
        frac_edg = one - (frac_sat + frac_dry)

        ! Find volume fraction of the model-level that is saturated within
        ! the "edge" area
        if ( cf_bulk(i,j,k-1) > cf_bulk(i,j,k) ) then
          ! There is less cloud area at k than k-1
          qc_tot = (qcl(i,j,k-1)+qcf(i,j,k-1))/cf_bulk(i,j,k-1)
          if ( supersat(i,j,k-1) - supersat(i,j,k) > qc_tot ) then
            frac_lev = qc_tot / ( supersat(i,j,k-1) - supersat(i,j,k) )
          else
            frac_lev = one
          end if
        else if ( cf_bulk(i,j,k) > cf_bulk(i,j,k-1) ) then
          ! There is more cloud area at k than k-1
          qc_tot = (qcl(i,j,k)+qcf(i,j,k))/cf_bulk(i,j,k)
          if ( supersat(i,j,k) - supersat(i,j,k-1) > qc_tot ) then
            frac_lev = qc_tot / ( supersat(i,j,k) - supersat(i,j,k-1) )
          else
            frac_lev = one
          end if
        else
          ! CF equal on both levels (edge contribution will be zero)
          frac_lev = zero
        end if

        ! Find saturated volume-fraction of the rho-level
        weight2 = frac_sat + frac_lev * frac_edg

        ! Find vertical interpolation weight
        weight1 = ( z_uv(i,j,k) - z_tq(i,j,k-1) )                              &
                / ( z_tq(i,j,k) - z_tq(i,j,k-1) )

        ! Interpolate to find combined volume-average buoyancy coefficients
        bt_rh =      weight2  * ( (one-weight1) * bt_cld(i,j,k-1)              &
                                      + weight1 * bt_cld(i,j,k) )              &
              + (one-weight2) * ( (one-weight1) * bt(i,j,k-1)                  &
                                      + weight1 * bt(i,j,k) )
        bq_rh =      weight2  * ( (one-weight1) * bq_cld(i,j,k-1)              &
                                      + weight1 * bq_cld(i,j,k) )              &
              + (one-weight2) * ( (one-weight1) * bq(i,j,k-1)                  &
                                      + weight1 * bq(i,j,k) )

        ! Compute dbdz
        dbdz_rh(i,j,k) = g * ( bt_rh * dsldz(i,j,k)                            &
                             + bq_rh * dqwdz(i,j,k) )
        ! Also compute gradient-adjusted version
        dbdz_ga_rh(i,j,k) = g * ( bt_rh * dsldz_ga(i,j,k)                      &
                                + bq_rh * dqwdz(i,j,k) )

      end do
    end do
  end do
!$OMP end do NOWAIT
  k = 1
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! At surface, just use buoyancy coefficients from the first theta-level
      dbdz_rh(i,j,k) = g * ( bt_gb(i,j,k) * dsldz(i,j,k)                       &
                           + bq_gb(i,j,k) * dqwdz(i,j,k) )
      dbdz_ga_rh(i,j,k) = g * ( bt_gb(i,j,k) * dsldz_ga(i,j,k)                 &
                              + bq_gb(i,j,k) * dqwdz(i,j,k) )
    end do
  end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

        ! Interpolate dbdz to theta-levels
        ! Note: dbdz(k) is defined on theta-level k-1
        weight1 = ( z_tq(i,j,k-1) - z_uv(i,j,k-1) )                            &
                / ( z_uv(i,j,k)   - z_uv(i,j,k-1) )
        dbdz(i,j,k) = (one-weight1) * dbdz_rh(i,j,k-1)                         &
                    +      weight1  * dbdz_rh(i,j,k)
        dbdz_ga(i,j,k) = (one-weight1) * dbdz_ga_rh(i,j,k-1)                   &
                       +      weight1  * dbdz_ga_rh(i,j,k)

        ! Also need to interpolate dsldz, dqwdz onto theta-levels for use
        ! in the RHcrit calculation
        dsldzm(i,j,k) = (one-weight1) * dsldz(i,j,k-1)                         &
                      +      weight1  * dsldz(i,j,k)
        dqwdzm(i,j,k) = (one-weight1) * dqwdz(i,j,k-1)                         &
                      +      weight1  * dqwdz(i,j,k)

      end do
    end do
  end do
!$OMP end do

end select  ! (i_interp_local)


!--------------------------------------------------
! Calculate modulus of shear on theta-levels
! dvdzm(k) is on theta-level(k-1)
!--------------------------------------------------
if (.not. l_subfilter_vert) then

!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        dzu = u_p(i,j,k) - u_p(i,j,k-1)
        dzv = v_p(i,j,k) - v_p(i,j,k-1)
        dvdzm(i,j,k) = max( 1.0e-12_r_bl,                                      &
                            sqrt(dzu*dzu + dzv*dzv) * rdz(i,j,k) )
      end do
    end do
  end do
!$OMP end do

else

  if ( model_type == mt_single_column ) then
    ! visc_m, visc_h need to be the 1D shear(k) on theta-level(k)

!$OMP do SCHEDULE(STATIC)
    do k = 2, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          dzu = u_p(i,j,k) - u_p(i,j,k-1)
          dzv = v_p(i,j,k) - v_p(i,j,k-1)
          dvdzm(i,j,k) = max( 1.0e-12_r_bl,                                    &
                              sqrt(dzu*dzu + dzv*dzv) * rdz(i,j,k) )
          visc_m(i,j,k-1) = dvdzm(i,j,k)
          visc_h(i,j,k-1) = dvdzm(i,j,k)
        end do
      end do
    end do
!$OMP end do

  else
    ! On entry, visc_m is 3D shear(k) on theta-level(k)

!$OMP do SCHEDULE(STATIC)
    do k = 2, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          dvdzm(i,j,k) = max( 1.0e-12_r_bl, visc_m(i,j,k-1) )
        end do
      end do
    end do
!$OMP end do
  end if  ! test on SCM

end if  ! test on l_subfilter_vert

if (l_subfilter_horiz .or. l_subfilter_vert) then
  if (smag_l_calc == smag_l_calc_use_geo) then
    ! use 3d grid geometric mean

!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          rmlmax2(i,j,k) = ( delta_smag(i,j) * delta_smag(i,j) *               &
                               dzl_charney(i,j,k) )**one_third
        end do
      end do
    end do
!$OMP end do

  else
    ! only use horizontal grid
!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          rmlmax2(i,j,k) = delta_smag(i,j)
        end do
      end do
    end do
!$OMP end do

 end if

  if ( l_rp2 .and. i_rp_scheme == i_rp2b ) then

!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          rmlmax2(i,j,k) = ( cs_rp(rp_idx) * rmlmax2(i,j,k) )**2
        end do
      end do
    end do
!$OMP end do

  else

!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          rmlmax2(i,j,k) = ( mix_factor * rmlmax2(i,j,k) )**2
        end do
      end do
    end do
!$OMP end do

  end if

!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        rneutml_sq(i,j,k) = one / (                                            &
                 one/( vkman*(z_tq(i,j,k) + z0m_eff_gb(i,j)) )**2              &
               + one/rmlmax2(i,j,k) )
      end do
    end do
  end do
!$OMP end do

end if
!-----------------------------------------------------------------------
! 2.1 Orographic enhancement of subgrid mixing
!-----------------------------------------------------------------------
!  Set-up 2D array for standard deviation of subgrid orography.
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    sigma_h(i,j) = zero
  end do
end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do l = 1, land_pts
  j=(land_index(l)-1)/pdims%i_end + 1
  i=land_index(l) - (j-1)*pdims%i_end
  sigma_h(i,j) = min( sd_orog(l), 300.0_r_bl )
end do
!$OMP end do
!-----------------------------------------------------------------------
!  Enhance resolved shear through unresolved subgrid drainage flows.
!-----------------------------------------------------------------------
if (sg_orog_mixing == sg_shear .or.                                            &
    sg_orog_mixing == sg_shear_enh_lambda) then

!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end

        if (sigma_h(i,j) > one ) then
          zpr = z_tq(i,j,k-1)/sigma_h(i,j)
          ! Height dependence, to reduce effect to zero with height
          !   gives z_scale~[1,0.95,0.5,0] at zpr=[0,0.6,1,1.7]
          weight1 = one_half*( one - tanh(4.0_r_bl*(zpr-one) ) )

          ! Take slope ~ sd/h_scale for small sd;
          !            tends to 0.2 for large sd
          slope = one / sqrt( 25.0_r_bl + (h_scale/sigma_h(i,j))**2 )

          dvdzm(i,j,k) = max ( dvdzm(i,j,k),                                   &
                               weight1*slope*t_drain*dbdz(i,j,k) )

          if (k==2 .and. BL_diag%l_dvdzm)                                      &
            BL_diag%dvdzm(i,j,1)=weight1*slope*t_drain*dbdz(i,j,k)

        end if
      end do
    end do
  end do
!$OMP end do
end if     ! sg_orog_mixing

!$OMP do SCHEDULE(STATIC)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ri(i,j,k)    = dbdz(i,j,k)    / ( dvdzm(i,j,k)*dvdzm(i,j,k) )
      ri(i,j,k)    = max(min(ri(i,j,k),max_ri),-max_ri)
      ri_ga(i,j,k) = dbdz_ga(i,j,k) / ( dvdzm(i,j,k)*dvdzm(i,j,k) )
    end do
  end do
end do
!$OMP end do

if (BL_diag%l_gradrich) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%gradrich(i,j,k)=ri(i,j,k)
      end do
    end do
  end do
!$OMP end do
end if

if (BL_diag%l_dbdz) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%dbdz(i,j,k)=dbdz(i,j,k)
      end do
    end do
  end do
!$OMP end do
end if

if (BL_diag%l_dvdzm) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%dvdzm(i,j,k)=dvdzm(i,j,k)
      end do
    end do
  end do
!$OMP end do
end if
!$OMP end PARALLEL
!-----------------------------------------------------------------------
! 3.  Orographic formdrag - distributed drag option
!-----------------------------------------------------------------------
if (formdrag ==  explicit_stress) then
  !------------------------------------------------------------------
  !      Calculate stress profiles
  !------------------------------------------------------------------
  call fm_drag (                                                               &
  ! in levels
        land_pts, land_index, bl_levels,                                       &
  ! in fields
        u_p, v_p, tl, qw, bt_gb, bq_gb, rho_wet_tq,                            &
        z_uv, z_tq, z0m_eff_gb, zh_prev, rib_gb, sil_orog_land,                &
  ! out fields
        tau_fd_x(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
                 1:bl_levels),                                                 &
        tau_fd_y(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
                 1:bl_levels)                                                  &
        )
  !------------------------------------------------------------------
  !      Orographic stress diagnostics
  !------------------------------------------------------------------
  if (BL_diag%l_ostressx) then
!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(k,j,i)                                                           &
!$OMP SHARED(bl_levels,tdims,BL_diag,tau_fd_x)
    do k = 1, bl_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          BL_diag%ostressx(i,j,k)=tau_fd_x(i,j,k)
        end do
      end do
    end do
!$OMP end PARALLEL do
  end if

  if (BL_diag%l_ostressy) then
!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(k,j,i)                                                           &
!$OMP SHARED(bl_levels,tdims,BL_diag,tau_fd_y)
    do k = 1, bl_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          BL_diag%ostressy(i,j,k)=tau_fd_y(i,j,k)
        end do
      end do
    end do
!$OMP end PARALLEL do
  end if

end if
!-----------------------------------------------------------------------
! 4. Apply dynamic diagnosis of shear-driven layers.
!-----------------------------------------------------------------------
!       In cases where the parcel ascent continues right through the
!       boundary layer, we diagnose cumulus only where the surface
!       buoyancy flux is sufficiently unstable, testing the ratio of
!       the depth of the inversion to the Obukhov length. A value of
!       1 -- 2 is reasonable for this test and 1.6 is selected, but
!       no great precision is attached to this value. Since this is
!       of importance mainly at sea-points, to avoid complications
!       with coastal tiling, the scheme operates only at points
!       where the land fraction is below 0.5.

!$OMP PARALLEL DEFAULT(SHARED) private(i, j, k, ii, ntop, z_scale)

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ! First override the provisional cumulus diagnosis if the
    ! actual surface buoyancy flux is stable
    if ( fb_surf(i,j) < zero ) then
      cumulus(i,j)   = .false.
      l_shallow(i,j) = .false.
      ntml(i,j)      = 1
      zh(i,j)        = z_uv(i,j,2)
      dzh(i,j)       = rmdi
    end if
  end do
end do
!$OMP end do


if (idyndiag == DynDiag_ZL) then

  ! Original version - causes spuriously deep boundary layers if
  ! BL_LEVELS is >> 3km

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

      if ( flandg(i,j) < one_half ) then
        ntop = min(ntpar(i,j),bl_levels-1)
        if ( -z_uv(i,j,ntop+1) * recip_l_mo_sea(i,j) < near_neut_z_on_l ) then
          cumulus(i,j)   = .false.
          l_shallow(i,j) = .false.
          ntml(i,j)      = ntop
          zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
        end if
      end if

    end do
  end do
!$OMP end do

else if (idyndiag == DynDiag_ZL_corrn) then

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

        !------------------------------------------------------------
        ! As original code above, except restrict depth scale to <3km
        ! and, if near-neutral, ignore the BL depth diagnosed by the
        ! adiabatic parcel (ie ZH, NTML) completely.
        !------------------------------------------------------------
      if ( flandg(i,j) < one_half ) then
        ntop = min(ntpar(i,j),bl_levels-1)
        z_scale = min( 3000.0_r_bl, z_uv(i,j,ntop+1) )
        if ( -z_scale*recip_l_mo_sea(i,j) < near_neut_z_on_l ) then
          cumulus(i,j)   = .false.
          l_shallow(i,j) = .false.
          ntml(i,j)      = 1
          zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
        end if
      end if

    end do
  end do
!$OMP end do

else if (idyndiag == DynDiag_ZL_CuOnly) then

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

        !------------------------------------------------------------
        ! As DynDiag_ZL_corrn but only affects cumulus points and
        ! points that are entirely sea.  Note that DynDiag_ZL_corrn
        ! has been found to switch off non-local mixing in
        ! stratocumulus, where surface fluxes are typically small
        !------------------------------------------------------------
      if ( cumulus(i,j) .and. flandg(i,j) < 0.01_r_bl ) then
        ntop = min(ntpar(i,j),bl_levels-1)
        z_scale = min( 3000.0_r_bl, z_uv(i,j,ntop+1) )
        if ( -z_scale*recip_l_mo_sea(i,j) < near_neut_z_on_l ) then
          ! - ZH/L indicates BL close to neutral
          dynamic_bl_diag(i,j) = .true.
          cumulus(i,j)   = .false.
          l_shallow(i,j) = .false.
          if (l_fix_dyndiag) then
            ntml(i,j) = 1
          else
            ! note that ntop can be very high making a spuriously
            ! deep well-mixed PBL - hence depricated
            ntml(i,j) = ntop
          end if
          zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
        end if
      end if

    end do
  end do
!$OMP end do

else if (idyndiag == DynDiag_Ribased ) then
  !------------------------------------------------------------
  ! As DynDiag_ZL_CuOnly but also allow ZH(Ri) to overrule the
  ! Cumulus diagnosis
  !------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      topbl(i,j)           = .false.
    end do
  end do
!$OMP end do
  !---------------------------------------------------------------
  !  Loop over levels to find Ri > RiCrit_sharp (=0.25) to find
  !  level to which Ri really is close to neutral
  !---------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
  do ii = pdims%i_start, pdims%i_end, pdims_omp_block
    do k = 2, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = ii, min(((ii+pdims_omp_block)-1),pdims%i_end)
          if ( .not. topbl(i,j) .and.                                          &
            (ri_ga(i,j,k) >  RiCrit_sharp .or. k > bl_levels-1) ) then
            topbl(i,j) = .true.
            zh_local(i,j) = z_uv(i,j,k)
          end if
        end do  ! Loop over points
      end do  ! Loop over points
    end do  ! Loop over levels
  end do
!$OMP end do
  !---------------------------------------------------------------
  !  Overrule Cumulus flag where close to neutral BL
  !---------------------------------------------------------------
  if (l_fix_dyndiag) then

    ! If the RP2b scheme is switched on,
    ! update value of zhloc_depth_fac
    if ( l_rp2 .and. i_rp_scheme == i_rp2b ) then
      zhloc_depth_fac = zhloc_depth_fac_rp(rp_idx)
    end if

!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end

        if ( cumulus(i,j) .and. flandg(i,j) < 0.01_r_bl ) then
          ntop = min(ntpar(i,j),bl_levels-1)
          z_scale = min( 3000.0_r_bl, z_uv(i,j,ntop+1) )
          if ( zh_local(i,j)                                                   &
                  > zh(i,j)+zhloc_depth_fac*(zhpar(i,j)-zh(i,j)) ) then
                ! ZH(Ri>RiCrit) more than zhloc_depth_fac up the
                ! cloud layer, indicating significant shear disruption
            dynamic_bl_diag(i,j) = .true.
            cumulus(i,j)   = .false.
            l_shallow(i,j) = .false.
            ntml(i,j)      = ntop
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          else if ( -z_scale*recip_l_mo_sea(i,j) < near_neut_z_on_l ) then
            ! - ZH/L indicates BL close to neutral so overrule
            ! cumulus diagnosis and leave mixing to the local
            ! Ri-based scheme, given no better information on where
            ! unstable bl top should be
            dynamic_bl_diag(i,j) = .true.
            cumulus(i,j)   = .false.
            l_shallow(i,j) = .false.
            ntml(i,j)      = 1
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          end if
        end if ! Cu over sea

      end do
    end do
!$OMP end do

  else  ! old version: not good to set ntml=ntop for small zh/L

!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end

        if ( cumulus(i,j) .and. flandg(i,j) < 0.01_r_bl ) then
          ntop = min(ntpar(i,j),bl_levels-1)
          z_scale = min( 3000.0_r_bl, z_uv(i,j,ntop+1) )
          if ( -z_scale*recip_l_mo_sea(i,j) < near_neut_z_on_l .or.            &
                ! - ZH/L indicates BL close to neutral
             zh_local(i,j) > zh(i,j)+zhloc_depth_fac*(z_scale-zh(i,j))         &
                ! ZH(Ri>RiCrit) more than zhloc_depth_fac up the
                ! cloud layer, indicating significant shear disruption
            ) then
            dynamic_bl_diag(i,j) = .true.
            cumulus(i,j)   = .false.
            l_shallow(i,j) = .false.
            ntml(i,j)      = ntop
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          end if
        end if ! Cu over sea

      end do
    end do
!$OMP end do

  end if  ! l_fix_dyndiag




end if  ! tests on iDynDiag

!$OMP end PARALLEL
!-----------------------------------------------------------------------
! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------
! 5.1  Calculate the non-local terms and diffusion coefficients
!-----------------------------------------------------------------------
! Set NTML_NL to NTML as passed in from initial diagnosis routine
!-----------------------------------------------------------------------

!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(j,i)                                                             &
!$OMP SHARED(pdims,ntml_nl,ntml)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ntml_nl(i,j) = ntml(i,j)
  end do
end do
!$OMP end PARALLEL do

if (nl_bl_levels < bl_levels) then
      ! Set to huge value to make if-test in KMKHZ redundent
  zmaxb_for_dsc = 1.0e10_r_bl
  zmaxt_for_dsc = zmaxb_for_dsc

!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(j,i)                                                             &
!$OMP SHARED(pdims,ntml_nl,nl_bl_levels,zh,z_uv)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( ntml_nl(i,j) > nl_bl_levels-1 ) then
        ntml_nl(i,j) = nl_bl_levels-1
        zh(i,j)      = z_uv(i,j,ntml_nl(i,j)+1)
      end if
    end do
  end do
!$OMP end PARALLEL do
else
  zmaxb_for_dsc = 2500.0_r_bl
  zmaxt_for_dsc = 3000.0_r_bl
end if
!-----------------------------------------------------------------------
! Initialise non-local K and fluxes to zero; necessary for levels
! above NL_BL_LEVELS
!-----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( bl_levels, pdims, ftl, fqw, rhokmz, rhokhz, rhokm_top,          &
!$OMP          rhokh_top, tke_nl, f_ngstress, rhof2, rhofsc, ft_nt, fq_nt,     &
!$OMP          tothf_zh, tothf_zhsc, totqf_zh, totqf_zhsc, ft_nt_dscb,         &
!$OMP          zhnl, zh, fq_nt_dscb, tke_loc ) private( i, j, k )
!$OMP  do SCHEDULE(STATIC)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ftl(i,j,k) = zero
      fqw(i,j,k) = zero
      rhokmz(i,j,k) = zero
      rhokhz(i,j,k) = zero
      rhokm_top(i,j,k) = zero
      rhokh_top(i,j,k) = zero
      tke_nl(i,j,k)    = zero
      tke_loc(i,j,k)   = zero
      f_ngstress(i,j,k) = zero
      rhof2(i,j,k)  = zero
      rhofsc(i,j,k) = zero
      ! Initialise Lock-Whelan non-gradient terms to zero
      ! Only calculated in KMKHZ9C
      ft_nt(i,j,k)  = zero
      fq_nt(i,j,k)  = zero
    end do
  end do
end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ft_nt(i,j,1)  = zero
    fq_nt(i,j,1)  = zero
    tothf_zh(i,j)   = zero
    tothf_zhsc(i,j) = zero
    totqf_zh(i,j)   = zero
    totqf_zhsc(i,j) = zero
    ft_nt_dscb(i,j) = zero
    fq_nt_dscb(i,j) = zero
    zhnl(i,j) = zh(i,j)  ! initialise non-local PBL depth
                         ! to that diagnosed in conv_diag
  end do
end do
!$OMP end do
!$OMP end PARALLEL

! Repeat relevant parts of last 2 loops for water tracers
if (l_wtrac) then
  do i_wt = 1, n_wtrac

!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( i_wt, bl_levels, pdims, wtrac_bl) private( i, j, k)
!$OMP  do SCHEDULE(STATIC)
    do k = 2, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          wtrac_bl(i_wt)%fqw(i,j,k)   = zero
          wtrac_bl(i_wt)%fq_nt(i,j,k) = zero
        end do
      end do
    end do
!$OMP end do

!$OMP  do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        wtrac_bl(i_wt)%fq_nt(i,j,1)    = zero
        wtrac_bl(i_wt)%fq_nt_dscb(i,j) = zero
        wtrac_bl(i_wt)%totqf_zh(i,j)   = zero
        wtrac_bl(i_wt)%totqf_zhsc(i,j) = zero
      end do
    end do
!$OMP end do
!$OMP end PARALLEL

  end do   ! i_wt
end if

! switching off non_local_bl based on vertical Smagorinsky has been moved
! to readlsta/readlsta_4a/scm_shell so that it is only executed once

if (non_local_bl == on) then


  call kmkhz_9c (                                                              &
    !     in levels/switches
             nl_bl_levels,BL_diag, nSCMDpkgs,L_SCMDiags,                       &
    !     in fields
             p_theta_levels,rho_wet_tq,rho_mix,rho_mix_tq,t,q,qcl,qcf,         &
             cf_bulk, qw,tl, dzl_charney,rdz_charney_grid,z_tq,z_uv,           &
             rad_hr,micro_tends,                                               &
             bt,bq,btm,bqm,dqsdt,btm_cld,bqm_cld,a_qs,a_qsm,a_dqsdtm,          &
             u_s,fb_surf,rhostar,ntpar,zh_prev,                                &
             zhpar,z_lcl,zmaxb_for_dsc,zmaxt_for_dsc,l_shallow,wtrac_as,       &
    !     INOUT fields
             ftl,fqw,zhnl,dzh,cumulus,ntml_nl,w,etadot,t1_sd,q1_sd,wtrac_bl,   &
    !     out fields
             rhokmz,rhokhz,rhokm_top,rhokh_top,zhsc,                           &
             unstable,dsc,coupled,sml_disc_inv,dsc_disc_inv,                   &
             ntdsc,nbdsc,f_ngstress,tke_nl,                                    &
             grad_t_adj, grad_q_adj,                                           &
             rhof2, rhofsc, ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb,              &
             tothf_zh, tothf_zhsc, totqf_zh, totqf_zhsc,                       &
             kent, we_lim, t_frac, zrzi,                                       &
             kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,                       &
             kplume                                                            &
             )


else   ! not NON_LOCAL_BL

       !-------------------------------------------------------------
       ! Set all variables from the non-local scheme to zero or "off"
       !  - reset all fluxes and K's arising from the non-local scheme
       !-------------------------------------------------------------

!$OMP PARALLEL DEFAULT(none) private(j,i,ient)                                 &
!$OMP SHARED(pdims,unstable,fb_surf,cumulus,l_shallow,sml_disc_inv,ntpar,      &
!$OMP        ntml_nl,zhnl,grad_t_adj,grad_q_adj,dsc,dsc_disc_inv,ntdsc,nbdsc,  &
!$OMP        zhsc,dzh,coupled,kent,kent_dsc,t_frac,zrzi,we_lim,t_frac_dsc,     &
!$OMP        zrzi_dsc,we_lim_dsc,kplume)
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ! surface mixed layer
      unstable(i,j) = (fb_surf(i,j) >  zero)
      cumulus(i,j) = .false.
      l_shallow(i,j) = .false.
      sml_disc_inv(i,j) = 0
      ntpar(i,j)   = 0
      ntml_nl(i,j) = -1    ! to ensure correct diagnostics
      zhnl(i,j)      = zero
      dzh(i,j)       = zero
      grad_t_adj(i,j) = zero
      grad_q_adj(i,j) = zero
      ! decoupled mixed layer
      dsc(i,j)     = .false.
      dsc_disc_inv(i,j) = 0
      ntdsc(i,j)   = 0
      nbdsc(i,j)   = 0
      zhsc(i,j)    = zero
      coupled(i,j) = .false.
      ! entrainment variables for non-local tracer mixing
      kent(i,j) = 2
      kent_dsc(i,j) = 2
      ! kplume for bi-modal cloud scheme
      kplume(i,j) = 1
    end do
  end do
!$OMP end do NOWAIT
  do ient = 1, 3
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        t_frac(i,j,ient) = zero
        zrzi(i,j,ient)   = zero
        we_lim(i,j,ient) = zero
        t_frac_dsc(i,j,ient) = zero
        zrzi_dsc(i,j,ient)   = zero
        we_lim_dsc(i,j,ient) = zero
      end do
    end do
!$OMP end do NOWAIT
  end do
!$OMP end PARALLEL
end if  ! test on NON_LOCAL_BL

! default values
!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(j,i)                                                             &
!$OMP SHARED(pdims,l_shallow_cth)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    l_shallow_cth(i,j) = .false.
  end do
end do
!$OMP end PARALLEL do

if (blending_option == blend_cth_shcu_only) then
  ! only going to use the parcel top as the length scale in blending
  ! if the convection is shallow, where we define shallow convection here
  ! as the resolved cloud top (from cf_bluk) being below shallow_cu_maxtop
!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(none)                                                            &
!$OMP private(i,j,k)                                                           &
!$OMP SHARED(pdims,bl_levels,z_uv,l_shallow_cth,cumulus,ntml,cf_bulk,          &
!$OMP        shallow_cu_maxtop,cloud_base_found,cloud_top_found)
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      cloud_base_found(i,j) = .false.
      cloud_top_found(i,j)  = .false.
    end do
  end do
!$OMP end do
  do k = 3, bl_levels-1
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( cumulus(i,j) .and. k > ntml(i,j)+1 .and.                          &
             .not. cloud_base_found(i,j) .and. cf_bulk(i,j,k) > sc_cftol ) then
             ! found a cumulus cloud base
          cloud_base_found(i,j) = .true.
        end if
        if ( cloud_base_found(i,j) .and. .not. cloud_top_found(i,j) .and.      &
          ! found cloud-base but not yet reached cloud-top
             cf_bulk(i,j,k+1) < sc_cftol ) then
          ! got to cloud-top
          cloud_top_found(i,j) = .true.
          l_shallow_cth(i,j)   = z_uv(i,j,k+1) < shallow_cu_maxtop
        end if
      end do
    end do
!$OMP end do
  end do
!$OMP end PARALLEL
end if  ! blending_option
!-----------------------------------------------------------------------
! 5.1  Call local coeff calculation for levels 2 to bl_levels
!-----------------------------------------------------------------------
call ex_coef (                                                                 &
! in levels/logicals
   bl_levels,k_log_layr,nSCMDpkgs,L_SCMDiags,BL_diag,                          &
! in fields
   sigma_h,flandg,dbdz,dvdzm,ri,rho_wet_tq,z_uv,z_tq,z0m_eff_gb,               &
   h_blend_orog,zhpar,ntpar,ntml_nl,ntdsc,nbdsc,u_p,v_p,u_s,fb_surf,           &
   qw,tl,l_shallow_cth,rmlmax2, rneutml_sq, delta_smag,                        &
! in/out fields
   cumulus,weight_1dbl,                                                        &
! out fields
   lambda_min,zh_local,ntml_local,elm,elh,elh_rho,rhokm,rhokh_th,              &
   fm_3d,fh_3d,tke_loc                                                         &
   )

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
! interpolate rhokh_th to rho levels 2 to bl_levels

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP  private(i, j, k, weight1, weight2, weight3, lambdah, z_scale,           &
!$OMP          vkz, f_log)
!$OMP  do SCHEDULE(STATIC)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      weight1 = r_theta_levels(i,j,k) -                                        &
              r_theta_levels(i,j, k-1)
      weight2 = r_theta_levels(i,j,k) -                                        &
              r_rho_levels(i,j,k)
      weight3 = r_rho_levels(i,j,k) -                                          &
              r_theta_levels(i,j,k-1)
      if ( k  ==  bl_levels ) then
              ! assume rhokh_th(BL_LEVELS+1) is zero
        rhokh(i,j,k) = ( weight2/weight1 ) * rhokh_th(i,j,k)
        if (blending_option /= off) weight_1dbl_rho(i,j,k) =                   &
                               (weight2/weight1) * weight_1dbl(i,j,k)
      else
        rhokh(i,j,k) = (weight3/weight1) * rhokh_th(i,j,k+1)                   &
                     + (weight2/weight1) * rhokh_th(i,j,k)
        if (blending_option /= off) weight_1dbl_rho(i,j,k) =                   &
                              (weight3/weight1)*weight_1dbl(i,j,k+1)           &
                            + (weight2/weight1)*weight_1dbl(i,j,k)
      end if

      if (local_fa == free_trop_layers) then
        ! elh already included in rhokh_th so no need to calculate
        ! here, but interpolate elh separately for diagnostic
        if (BL_diag%l_elh3d) then
          if ( k  ==  bl_levels ) then
            ! assume rhokh_th(BL_LEVELS+1) is zero
            elh_rho(i,j,k) = ( weight2/weight1 ) * elh(i,j,k)
          else
            elh_rho(i,j,k) =                                                   &
              weight3/weight1 *                                                &
                      elh(i,j,k+1)                                             &
             +weight2/weight1 *                                                &
                      elh(i,j,k)
          end if
          BL_diag%elh3d(i,j,k)=elh_rho(i,j,k)
        end if
      else
        if ((sbl_op/=equilibrium_sbl) .or. (fb_surf(i,j) >  zero)) then
          !------------------------------------------------------------
          !  Include mixing length, ELH, in RHOKH.
          !  Code moved from EX_COEF to avoid interpolation
          !------------------------------------------------------------
          if ( k >= ntml_local(i,j)+2 .and. l_full_lambdas .and.               &
               local_fa == to_sharp_across_1km ) then
            ! Assuming only LOCAL_FA = "to_sharp_across_1km" option
            ! will have L_FULL_LAMBDAS.
            ! If other LOCAL_FA options are coded here then
            ! changes must be included in section 2.1 of ex_coef
            if (l_rp2) then
              lambdah = max ( lambda_min ,                                     &
                              par_mezcla_rp(rp_idx)*zh_local(i,j) )
            else
              lambdah = max ( lambda_min , 0.15_r_bl*zh_local(i,j) )
            end if
            z_scale = 1000.0_r_bl
            weight1 = one_half*( one -                                         &
                        tanh(3.0_r_bl*((z_uv(i,j,k)/z_scale )-one) ) )
            lambdah = lambdah * weight1                                        &
                         + lambda_min*( one -  weight1)
            ! no need to do log profile correction as klog_layr eq 2
            vkz = vkman * ( z_uv(i,j,k) + z0m_eff_gb(i,j) )
            elh_rho(i,j,k) = vkz / (one + vkz/lambdah )
          end if
          ! Reinstate UKV drainage flow bug here, where lambdah was not
          ! enhanced as intended (and as was done in ex_coef)!
          if (sg_orog_mixing == sg_shear_enh_lambda) then
            if (l_rp2) then
              lambdah = max ( lambda_min ,                                     &
                              par_mezcla_rp(rp_idx)*zh_local(i,j) )
            else
              lambdah = max ( lambda_min , 0.15_r_bl*zh_local(i,j) )
            end if
            if (k >= ntml_local(i,j)+2 .and. .not. l_full_lambdas) then
              lambdah = lambda_min
            end if
            if (k <= k_log_layr) then
              vkz   = vkman * ( z_tq(i,j,k) - z_tq(i,j,k-1) )
              f_log = log( ( z_tq(i,j,k) + z0m_eff_gb(i,j)   ) /               &
                           ( z_tq(i,j,k-1) + z0m_eff_gb(i,j) ) )
              elh_rho(i,j,k) = vkz / ( f_log + vkz / lambdah )
            else
              vkz = vkman * ( z_uv(i,j,k) + z0m_eff_gb(i,j) )
              elh_rho(i,j,k) = vkz / (one + vkz/lambdah )
            end if
          end if
          ! End of UKV bug!

          if (BL_diag%l_elh3d) BL_diag%elh3d(i,j,k)=elh_rho(i,j,k)

          rhokh(i,j,k) = elh_rho(i,j,k) * rhokh(i,j,k)

        end if   ! test on sbl_op
      end if   ! test on local_fa = free_trop_layers

          ! Finally multiply RHOKH by dry density
      if (l_mr_physics)                                                        &
         rhokh(i,j,k) = rho_mix(i,j,k) * rhokh(i,j,k)

    end do
  end do
end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
      !----------------------------------------------------------
      ! Use local NTML if significantly higher (to allow for
      ! local averaging) than the non-local or if the non-local
      ! is on the ground (=1)
      !----------------------------------------------------------
    if ( .not. cumulus(i,j) .and.                                              &
              ( ntml_local(i,j)  >   ntml_nl(i,j)+1                            &
                .or. ntml_nl(i,j)  ==  1 )            ) then
      ntml(i,j) = ntml_local(i,j)
      sml_disc_inv(i,j) = 0   ! reset flag for subgrid inversion
    else
      ntml(i,j) = ntml_nl(i,j)
    end if
      !----------------------------------------------------------
      ! If local NTML is higher than NTDSC then ignore DSC layer
      ! for diagnostics but keep mixing associated with it
      !----------------------------------------------------------
    if ( ntml_local(i,j)  >   ntdsc(i,j)+1 ) then
      dsc_disc_inv(i,j) = 0
      ntdsc(i,j) = 0
      nbdsc(i,j) = 0
      zhsc(i,j)  = zero
      dsc(i,j)   = .false.
      coupled(i,j) = .false.
    end if
  end do
end do
!$OMP end do
!$OMP end PARALLEL

! Calculate max of two coeffs
call kmkh (                                                                    &
! in data
   bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                                     &
   ntml,cumulus,ntdsc,sml_disc_inv,dsc_disc_inv,                               &
   weight_1dbl, weight_1dbl_rho,                                               &
! INOUT data
   rhokm,rhokh,rhokmz,rhokhz,rhokm_top,rhokh_top,tke_loc                       &
   )

if (BL_diag%l_weight1d) then
!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(k,j,i)                                                           &
!$OMP SHARED(bl_levels,pdims,BL_diag,weight_1dbl)
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%weight1d(i,j,k)=weight_1dbl(i,j,k)
      end do
    end do
  end do
!$OMP end PARALLEL do
end if

if (BL_diag%l_rhokm) then
!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(k,j,i)                                                           &
!$OMP SHARED(bl_levels,pdims,BL_diag,rhokm)
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokm(i,j,k)=rhokm(i,j,k)
      end do
    end do
  end do
!$OMP end PARALLEL do
end if

if (BL_diag%l_rhokh) then
!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(k,j,i)                                                           &
!$OMP SHARED(bl_levels,pdims,BL_diag,rhokh)
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokh(i,j,k)=rhokh(i,j,k)
      end do
    end do
  end do
!$OMP end PARALLEL do
end if

! Calculation of TKE diagnostic.
! Stored on theta-levels with TKE(K) on theta-level(K-1),
! consistent with RHOKM(K), RI(K), etc.
! The K=1 value could be set to a diagnosed surface value (eg as a
! function of ustar, wstar) but is currently just set to zero
if (BL_diag%l_tke) then

  if (var_diags_opt == original_vars) then
    ! BL_diag%tke currently contains (the reciprocal of) the non-local
    ! (SML and DSC) mixed layer timescale, calculated in excf_nl

!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED(BL_diag, rho_wet_tq, rhokm, dbdz, bl_levels, pdims)              &
!$OMP  private(i, j, k, recip_time_sbl, recip_time_cbl)
    do k = 2, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end

          ! stable timescale
          recip_time_sbl = sqrt( max(dbdz(i,j,k),1.0e-5_r_bl) )/0.7_r_bl
          recip_time_cbl = BL_diag%tke(i,j,k)

          ! TKE diagnostic - taking 5 m2/s2 as a suitable maximum
          BL_diag%tke(i,j,k)=  min( max_tke,                                   &
                ( rhokm(i,j,k) / rho_wet_tq(i,j,k-1) )*(                       &
                   recip_time_cbl + recip_time_sbl ) )

        end do
      end do
    end do
!$OMP end PARALLEL do

  else if (var_diags_opt == split_tke_and_inv) then
    ! Combine the, separately calculated, local and non-local TKE diagnostics

!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED(BL_diag, tke_nl, tke_loc, rho_wet_tq, weight_1dbl,               &
!$OMP         tke_diag_fac, bl_levels, pdims)                                  &
!$OMP  private(i, j, k)
    do k = 2, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end

          ! TKE diagnostic
          ! Currently tke_nl is strictly rho*sigma_w^2 while tke_loc is TKE
          ! Assume isotropic turb so TKE_nl = 3/2 sigma_w^2
          tke_nl(i,j,k) = weight_1dbl(i,j,k) * 1.5_r_bl *                      &
                                            tke_nl(i,j,k)/rho_wet_tq(i,j,k-1)
          BL_diag%tke(i,j,k) = max( tke_loc(i,j,k), tke_nl(i,j,k) )

          ! Multiply by tuning factor
          BL_diag%tke(i,j,k) = tke_diag_fac * BL_diag%tke(i,j,k)

          ! Keep TKE below a sensible max value of max_tke
          BL_diag%tke(i,j,k) = min( max_tke, BL_diag%tke(i,j,k) )
          ! Applying this limit can occasionally cause the length-scale
          ! Km / sqrt(w_var) to become unrealistically large, since no
          ! equivalent limiting is done on Km.

        end do
      end do
    end do
!$OMP end PARALLEL do

  end if  ! var_diags_opt

  if ( i_bm_ez_opt == i_bm_ez_entpar ) then
    ! Calculate mixing-length to pass to bimodal cloud scheme,
    ! using Km and TKE.
!$OMP PARALLEL DEFAULT(none) private( i, j, k, weight1 )                       &
!$OMP SHARED( bl_diag, mix_len_bm, rhokm, rho_wet_tq, pdims, tdims, bl_levels, &
!$OMP         sml_disc_inv, ntml_nl, zhnl, dsc_disc_inv, ntdsc, zhsc, z_uv )
!$OMP do SCHEDULE(STATIC)
    do k = 1, bl_levels-1
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          ! TKE and rhokm are on boundary-layer theta-levels (surface at k=1),
          ! mix_len_bm and rho_wet_tq are on UM theta-levels (surface at k=0).
          if ( BL_diag%tke(i,j,k+1) > 0.0 ) then
            mix_len_bm(i,j,k) = max( ( rhokm(i,j,k+1) / rho_wet_tq(i,j,k) )    &
                                   / sqrt( two_thirds * BL_diag%tke(i,j,k+1) ),&
                                     bm_tiny )
          else
            mix_len_bm(i,j,k) = bm_tiny
          end if
        end do
      end do
    end do
!$OMP end do
    ! Recalculate at sub-grid inversions, for improved numerical behaviour.
    ! Note we should always have z_uv(ntml+1) < zh < z_uv(ntml+2).
    ! However there is some code towards the end of kmkhz_9c which can
    ! increment ntml/ntdsc due to the inversion rise or fall over the
    ! current timestep, but doesn't increment zhnl/zhsc to keep it consistent.
    ! We need to check for this to get sensible interpolation weights...
!$OMP do SCHEDULE(STATIC)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( sml_disc_inv(i,j) > 0 .and. ntml_nl(i,j) > 1 ) then
          k = ntml_nl(i,j) + 1
          if ( zhnl(i,j) < z_uv(i,j,k)   )  k = k - 1
          if ( zhnl(i,j) > z_uv(i,j,k+1) )  k = k + 1
          weight1 = ( zhnl(i,j) - z_uv(i,j,k) )/( z_uv(i,j,k+1) - z_uv(i,j,k) )
          mix_len_bm(i,j,k) = mix_len_bm(i,j,k)                                &
                            + weight1 * max( 0.5 * mix_len_bm(i,j,k-1)         &
                                                 - mix_len_bm(i,j,k), 0.0 )
        end if
        if ( dsc_disc_inv(i,j) > 0 .and. ntdsc(i,j) > 1 ) then
          k = ntdsc(i,j) + 1
          if ( zhsc(i,j) < z_uv(i,j,k)   )  k = k - 1
          if ( zhsc(i,j) > z_uv(i,j,k+1) )  k = k + 1
          weight1 = ( zhsc(i,j) - z_uv(i,j,k) )/( z_uv(i,j,k+1) - z_uv(i,j,k) )
          mix_len_bm(i,j,k) = mix_len_bm(i,j,k)                                &
                            + weight1 * max( 0.5 * mix_len_bm(i,j,k-1)         &
                                                 - mix_len_bm(i,j,k), 0.0 )
        end if
      end do
    end do
!$OMP end do NOWAIT
    ! Set values above bl_levels
!$OMP do SCHEDULE(STATIC)
    do k = bl_levels, tdims%k_end
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          mix_len_bm(i,j,k) = bm_tiny
        end do
      end do
    end do
!$OMP end do NOWAIT
!$OMP end PARALLEL
  end if  ! ( i_bm_ez_opt == i_bm_ez_entpar )

  if ( l_subgrid_qcl_mp .or. l_wvar_for_conv .or.                              &
       i_cld_vn == i_cld_bimodal .or.                                          &
       (i_cld_vn == i_cld_pc2 .and.                                            &
        i_pc2_init_method == pc2init_bimodal) ) then

    ! Set bl_w_var to mimimum value.
    ! Prevents any unset values in the prognostic
    ! or any very low values being passed through
    ! to the microphysics turbulence call

!$OMP PARALLEL DEFAULT(none) private(k,j,i)                                    &
!$OMP SHARED(tdims, bl_levels,pdims,BL_diag,bl_w_var,var_diags_opt)

!$OMP do SCHEDULE(STATIC)
    do k = 2, tdims%k_end+1
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          bl_w_var(i,j,k) = 1.0e-12_r_bl
        end do
      end do
    end do
!$OMP end do

    if (var_diags_opt == original_vars) then
      ! Original version: w_var = TKE
!$OMP do SCHEDULE(STATIC)
      do k = 2, bl_levels
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            if (BL_diag%tke(i,j,k) > 1.0e-12_r_bl) then
              bl_w_var(i,j,k) = BL_diag%tke(i,j,k)
            end if
          end do
        end do
      end do
!$OMP end do
    else if (var_diags_opt == split_tke_and_inv) then
      ! Improved version: w_var = 2/3 TKE
!$OMP do SCHEDULE(STATIC)
      do k = 2, bl_levels
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            if (BL_diag%tke(i,j,k) > 1.0e-12_r_bl) then
              bl_w_var(i,j,k) = two_thirds * BL_diag%tke(i,j,k)
            end if
          end do
        end do
      end do
!$OMP end do
    end if ! test on var_diags_opt
!$OMP end PARALLEL

  end if ! l_subgrid_qcl_mp .or. l_wvar_for_conv

  if (var_diags_opt == original_vars) then
    ! at this point, BL_diag%tke really contains 1.5*sigma_w^2. To make it look
    ! a bit more like TKE near the surface, we will keep it constant below
    ! the max of rhokm_surf
!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED( pdims, rhokmz, BL_diag ) private(i, j, k, kmax )
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        kmax = maxloc(rhokmz(i,j,:), dim=1)
        do k = 2, kmax
          BL_diag%tke(i,j,k) = BL_diag%tke(i,j,kmax)
        end do
      end do
    end do
!$OMP end PARALLEL do
  else if (var_diags_opt == split_tke_and_inv) then
    ! at this point, tke_nl really contained 1.5*sigma_w^2. To make it look
    ! a bit more like TKE near the surface, we will keep it constant below
    ! the max of rhokm_surf (here we find the first local maximum in case there
    ! is a larger rhokmz in a resolved inversion
!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED( pdims,  bl_levels, rhokmz, BL_diag, tke_loc, tke_nl )           &
!$OMP  private(i, j, k, kp)
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        kp=2
        do while ( rhokmz(i,j,kp+1) > rhokmz(i,j,kp) .and. kp < bl_levels )
          kp = kp +1
        end do
        do k = 2, kp-1
          BL_diag%tke(i,j,k) = max( tke_loc(i,j,k), tke_nl(i,j,kp) )
        end do
      end do
    end do
!$OMP end PARALLEL do
  end if ! test on var_diags_opt

end if    ! BL_diag%L_tke

if (blending_option /= off .and. l_blend_isotropic) then
  !   ! Blended diffusion coefficients now held in rhokm and rhokh
  !   ! so copy to visc_m,h for horizontal diffusion too.
  !   ! Need to interpolate rhokh back to theta levels for visc_h

!$OMP  PARALLEL DEFAULT(SHARED) private(i, j, k, weight1, weight2,             &
!$OMP  weight3)

!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels-1
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( blending_option == blend_except_cu .and.                          &
             cumulus(i,j) .and. ntdsc(i,j) == 0) then
          ! pure cumulus layer so revert to Smag scheme
          visc_m(i,j,k) = visc_m(i,j,k)                                        &
                           *rneutml_sq(i,j,k)*fm_3d(i,j,k+1)
          visc_h(i,j,k) = visc_h(i,j,k)                                        &
                           *rneutml_sq(i,j,k)*fh_3d(i,j,k+1)
          if (k == bl_levels-1) then
            ! need to do bl_levels too apparently
            ! (see non-blend code below)
            visc_m(i,j,k+1) = visc_m(i,j,k+1)*rneutml_sq(i,j,k+1)
            visc_h(i,j,k+1) = visc_h(i,j,k+1)*rneutml_sq(i,j,k+1)
          end if
        else
          visc_m(i,j,k) = rhokm(i,j,k+1)/rho_wet_tq(i,j,k)
        end if
      end do
    end do
  end do
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels-1
    if ( k > 1 ) then
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          if (.not. ( blending_option == blend_except_cu .and.                 &
               cumulus(i,j) .and. ntdsc(i,j) == 0)) then
            ! not a pure cumulus layer or blending_option ne
            ! blend_except_cu, so use standard blending
            weight1 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
            weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
            visc_h(i,j,k) = ( weight2*(rhokh(i,j,k+1)/rho_mix(i,j,k+1))        &
                            + weight3*(rhokh(i,j,k)  /rho_mix(i,j,k)  ))       &
                                           / weight1
          end if
        end do
      end do
    else
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          if (.not. ( blending_option == blend_except_cu .and.                 &
                 cumulus(i,j) .and. ntdsc(i,j) == 0)) then
            ! not a pure cumulus layer or blending_option ne
            ! blend_except_cu, so use standard blending
            visc_h(i,j,k) = rhokh_th(i,j,k+1)/rho_wet_tq(i,j,k)
          end if
        end do
      end do
    end if
  end do
!$OMP end do
!$OMP end PARALLEL

else if (l_subfilter_horiz .or. l_subfilter_vert) then

  ! visc_m,h on in are just S and visc_m,h(k) are co-located with w(k)

!$OMP PARALLEL do SCHEDULE(STATIC)                                             &
!$OMP DEFAULT(none)                                                            &
!$OMP private(k,j,i)                                                           &
!$OMP SHARED(bl_levels,pdims,visc_m,rneutml_sq,visc_h,fh_3d,fm_3d,max_diff)

  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        visc_m(i,j,k) = visc_m(i,j,k)*rneutml_sq(i,j,k)
        visc_h(i,j,k) = visc_h(i,j,k)*rneutml_sq(i,j,k)
      end do
    end do

    if ( k < bl_levels ) then
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          ! stability functions are indexed with Ri, fm(k) on w(k-1)
          visc_h(i,j,k) = visc_h(i,j,k)*fh_3d(i,j,k+1)
          visc_m(i,j,k) = visc_m(i,j,k)*fm_3d(i,j,k+1)
          ! APL why apply this cap here for implicit vertical diffusion?
          ! (also applied in atm_step_phys_init for horiz diffn,
          ! that actually needs it)?
          visc_h(i,j,k) = min(visc_h(i,j,k),max_diff(i,j))
          visc_m(i,j,k) = min(visc_m(i,j,k),max_diff(i,j))
        end do
      end do
    end if
  end do
!$OMP end PARALLEL do

  ! visc_m and visc _h are now lambda^2*S*FM and lambda^2*S*FH

  if (l_subfilter_vert .and. blending_option == off) then

    ! visc_h_rho(k) is held on rho(k), same as BL's rhokh
    allocate (visc_h_rho(pdims%i_start:pdims%i_end,                            &
                         pdims%j_start:pdims%j_end, bl_levels))

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP SHARED(bl_levels, pdims, r_theta_levels, r_rho_levels, visc_h_rho,       &
!$OMP        visc_h, turb_startlev_vert, turb_endlev_vert, rhokm, visc_m,      &
!$OMP        rho_wet_tq, rhokh, rho_mix)                                       &
!$OMP private(i, j, k, weight1, weight2, weight3)
    do k = 2, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          weight1 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
          weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
          weight3 = r_rho_levels(i,j,k)   - r_theta_levels(i,j,k-1)
          if ( k  ==  bl_levels ) then
            ! assume visc_h(bl_levels) is zero (Ri and thence f_h not defined)
            visc_h_rho(i,j,k) = ( weight2/weight1 ) * visc_h(i,j,k-1)
          else
            visc_h_rho(i,j,k) = ( weight3/weight1 ) * visc_h(i,j,k)            &
                              + ( weight2/weight1 ) * visc_h(i,j,k-1)
          end if
        end do
      end do

      ! Overwrite the diffusion coefficients from the BL scheme
      !(RHOKM and RHOKH) with those obtained from the Smagorinsky scheme

      if (k >= turb_startlev_vert .and.                                        &
          k <= turb_endlev_vert) then
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            rhokm(i,j,k) = visc_m(i,j,k-1)*rho_wet_tq(i,j,k-1)
            rhokh(i,j,k) = visc_h_rho(i,j,k)*rho_mix(i,j,k)
          end do
        end do
      else
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            rhokm(i,j,k) = zero
            rhokh(i,j,k) = zero
          end do
        end do
      end if
    end do
!$OMP end PARALLEL do
    deallocate (visc_h_rho)

  end if ! L_subfilter_vert but not blend
end if ! L_subfilter_horiz or L_subfilter_vert or isotropic blending_option

!-----------------------------------------------------------------------
! Diagnose boundary layer type.
!      Seven different types are considered:
!      1 - Stable b.l.
!      2 - Stratocumulus over a stable surface layer.
!      3 - Well mixed buoyancy-driven b.l. (possibly with stratocumulus)
!      4 - Decoupled stratocumulus (not over cumulus).
!      5 - Decoupled stratocumulus over cumulus.
!      6 - Cumulus capped b.l.
!      7 - Shear-dominated unstable b.l.
!-----------------------------------------------------------------------
!      First initialise the type variables and set the depth diagnostics

! Top of surface mixed layer (Ksurf profile)
if (BL_diag%l_smltop) then
  BL_diag%smltop=zhnl
end if
! Top of decoupled stratocu layer
if (BL_diag%l_dsctop) then
  BL_diag%dsctop=zhsc
end if
! Height of diagnosis parcel top
if (BL_diag%l_zhpar) then
  BL_diag%zhpar=zhpar
end if

!$OMP PARALLEL DEFAULT(SHARED) private(i, j, k)

if (BL_diag%l_zht) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      bl_diag%zht(i,j) = max( zhnl(i,j) , zhsc(i,j) )
      if ( ntml(i,j)  >   ntml_nl(i,j) ) then
        bl_diag%zht(i,j) = max( bl_diag%zht(i,j), zh_local(i,j) )
      end if
    end do
  end do
!$OMP end do NOWAIT
end if

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ! Max height of BL turbulent mixing
    ! PBL depth diagnostic: start with top of non-local BL
    zh(i,j) = zhnl(i,j)
    if ( ntml(i,j)  >   ntml_nl(i,j) ) then
      ! Higher local K allowed so reset ZH, ZHT diagnostics
      zh(i,j)  = max( zh(i,j) , zh_local(i,j) )
      if (.not. l_fix_zh) then
        if (forced_cu == off) then
          ! going to use zhnl in tr_mix and ex_flux_uv which
          ! was formerly spuriously set to zh
          zhnl(i,j)=zh(i,j)
        else
          ! spuriously overwriting zh with zhnl
          zh(i,j)=zhnl(i,j)
        end if
      end if
    end if
    bl_type_1(i,j) = zero
    bl_type_2(i,j) = zero
    bl_type_3(i,j) = zero
    bl_type_4(i,j) = zero
    bl_type_5(i,j) = zero
    bl_type_6(i,j) = zero
    bl_type_7(i,j) = zero
  end do
end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if (dynamic_bl_diag(i,j)) then
      !       ! shear-dominated, via iDynDiag option
      bl_type_7(i,j) = one
    else
      if (.not. unstable(i,j) .and. .not. dsc(i,j) .and.                       &
          .not. cumulus(i,j)) then
        !         ! Stable b.l.
        bl_type_1(i,j) = one
      else if (.not. unstable(i,j) .and. dsc(i,j) .and.                        &
               .not. cumulus(i,j)) then
        !         ! Stratocumulus over a stable surface layer
        bl_type_2(i,j) = one
      else if (unstable(i,j) .and. .not. cumulus(i,j) .and.                    &
              .not. dsc(i,j) ) then
        !         ! Well mixed b.l. (possibly with stratocumulus)
        if ( ntml(i,j)  >   ntml_nl(i,j) ) then
          ! shear-dominated - currently identified
          ! by local NTML overriding non-local
          bl_type_7(i,j) = one
        else
          ! buoyancy-dominated
          bl_type_3(i,j) = one
        end if
      else if (unstable(i,j) .and. dsc(i,j) .and.                              &
                                      .not. cumulus(i,j)) then
        !         ! Decoupled stratocumulus (not over cumulus)
        bl_type_4(i,j) = one
      else if (dsc(i,j) .and. cumulus(i,j)) then
        !         ! Decoupled stratocumulus over cumulus
        bl_type_5(i,j) = one
      else if (.not. dsc(i,j) .and. cumulus(i,j)) then
        !         ! Cumulus capped b.l.
        bl_type_6(i,j) = one
        ! Label this shallow regime as 2.0_r_bl here, to be able to identify it
        ! in diagnostics_bl, but the "cumulus" stash output will still be 1.0
        if (l_shallow_cth(i,j)) bl_type_6(i,j) = 2.0_r_bl
      end if
    end if

  end do
end do
!$OMP end do NOWAIT
!$OMP end PARALLEL

!-----------------------------------------------------------------------
! 5.5 Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------
call ex_flux_tq (                                                              &
! in levels etc
    bl_levels,nSCMDpkgs,L_SCMDiags, BL_diag,                                   &
! in fields
    tl,qw,rdz_charney_grid, rhokh, rhokhz,grad_t_adj,grad_q_adj,               &
    rhof2, rhofsc, ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb, tothf_zh,             &
    tothf_zhsc, totqf_zh, totqf_zhsc, weight_1dbl_rho,                         &
    ntml_nl, ntdsc, nbdsc,                                                     &
! INOUT fields
    ftl,fqw,wtrac_bl                                                           &
    )

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(k, j, i)              &
!$OMP SHARED(bl_levels, pdims, f_ngstress, weight_1dbl)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      f_ngstress(i,j,k) = weight_1dbl(i,j,k) * f_ngstress(i,j,k)
    end do
  end do
end do
!$OMP end PARALLEL do

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 5.6.1 Calculate explicit surface fluxes of U and V on
!       P-grid for convection scheme
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(none)  private(j,i,l)                                   &
!$OMP SHARED(pdims,uw0,rhokm,u_p,u_0_px,vw0,v_p,v_0_px,wstar,wthvs,            &
!$OMP        non_local_bl,cu_over_orog,cumulus,fb_surf,zh,g,bt,l_param_conv,   &
!$OMP        ntml,ntml_nl,ntml_local,formdrag,tau_fd_x,tau_fd_y,ntdsc,BL_diag, &
!$OMP        land_pts,land_index,ho2r2_orog,l_shallow,bl_type_5,bl_type_6,     &
!$OMP        ntml_save,shallowc)

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    uw0(i,j) = -rhokm(i,j,1) *                                                 &
                      ( u_p(i,j,1) - u_0_px(i,j) )
    vw0(i,j) = -rhokm(i,j,1) *                                                 &
                      ( v_p(i,j,1) - v_0_px(i,j) )
  end do
end do
!$OMP end do NOWAIT
if (formdrag ==  explicit_stress) then
!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      uw0(i,j) = uw0(i,j) - tau_fd_x(i,j,1)
      vw0(i,j) = vw0(i,j) - tau_fd_y(i,j,1)
    end do
  end do
!$OMP end do NOWAIT
end if

!------------------------------------------------------------------------
! 5.7 Reset NTML to max of non-local scheme's turbulently mixed layers
!     when not cumulus, to stop mid-level from operating withing the PBL,
!     and calculate quantities to pass to convection scheme.
!------------------------------------------------------------------------

if (non_local_bl == on) then

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( cumulus(i,j) ) then
        if (.not. l_param_conv) then
          ntml(i,j) = max( 2, ntml_nl(i,j) - 1 )
        end if
      else
        ntml(i,j) = max( ntml_nl(i,j) , ntdsc(i,j) )
      end if
    end do
  end do
!$OMP end do NOWAIT

else ! non-local scheme not being used so always use local scheme values

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ntml(i,j) = ntml_local(i,j)
    end do
  end do
!$OMP end do NOWAIT

end if

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    wstar(i,j) = zero
    wthvs(i,j) = zero
    cu_over_orog(i,j) = zero
    if ( cumulus(i,j) ) then
      if ( fb_surf(i,j)  >   zero ) then
        wstar(i,j) = ( zh(i,j)*fb_surf(i,j) )**one_third
        wthvs(i,j) = fb_surf(i,j) / ( g * bt(i,j,1) )
      end if
      wstar(i,j) = max( 0.1_r_bl, wstar(i,j) )
    end if
    ! Limit explicitly calculated surface stresses
    ! to a physically plausible level.
    if ( uw0(i,j)  >=  5.0_r_bl ) then
      uw0(i,j) =  5.0_r_bl
    else if ( uw0(i,j)  <=  -5.0_r_bl ) then
      uw0(i,j) = -5.0_r_bl
    end if
    if ( vw0(i,j)  >=  5.0_r_bl ) then
      vw0(i,j) =  5.0_r_bl
    else if ( vw0(i,j)  <=  -5.0_r_bl ) then
      vw0(i,j) = -5.0_r_bl
    end if
    if (BL_diag%l_wstar .and. (fb_surf(i,j) >zero)) then
      BL_diag%wstar(i,j)= (zh(i,j)*fb_surf(i,j))**one_third
    end if
  end do
end do
!$OMP end do

if (l_param_conv) then

  ! Check for CUMULUS having been diagnosed over steep orography.
  ! Reset to false but keep NTML at NLCL (though decrease by 2 so that
  ! coupling between BL and convection scheme can be maintained).
  ! Reset type diagnostics.

!$OMP do SCHEDULE(STATIC)
  do l = 1, land_pts
    j=(land_index(l)-1)/pdims%i_end + 1
    i=land_index(l) - (j-1)*pdims%i_end
    if (cumulus(i,j) .and. ho2r2_orog(l)  >   900.0_r_bl) then
      cumulus(i,j)      = .false.
      l_shallow(i,j)    = .false.
      bl_type_5(i,j)    = zero
      bl_type_6(i,j)    = zero
      cu_over_orog(i,j) = one
      if (ntml(i,j)  >=  3) ntml(i,j) = ntml(i,j) - 2
    end if
  end do
!$OMP end do

  ! Check that CUMULUS and L_SHALLOW are still consistent
  ! and reset ntml back to where it was originally diagnosed
  ! if cumulus is still true in order to trigger convection
  ! at the correct level

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( .not. cumulus(i,j) ) l_shallow(i,j) = .false.
      if ( cumulus(i,j) ) ntml(i,j) = ntml_save(i,j)
    end do
  end do
!$OMP end do NOWAIT

end if    ! (l_param_conv)
!-----------------------------------------------------------------------
!     Set shallow convection diagnostic: 1.0 if L_SHALLOW (and CUMULUS)
!                                        0.0 if .not. CUMULUS
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( cumulus(i,j) .and. l_shallow(i,j) ) then
      shallowc(i,j) = one
    else
      shallowc(i,j) = zero
    end if
  end do
end do
!$OMP end do
!$OMP end PARALLEL

! ------rhcrit parametrization---------------------
! rhcpt is on theta-level K
! a_qs = a in documentation on theta level K
! a_dqsdt = b/exner in documentation on theta level K
! rho_wet_tq = rho on theta level K
! qsw_arr = qsat(Tl) on theta level K-1

if (i_rhcpt == rhcpt_tke_based .or. BL_diag%l_slvar .or. BL_diag%l_qwvar       &
    .or. BL_diag%l_slqw) then
  b2 = 15.0_r_bl
  root6 = sqrt(6.0_r_bl)

  ! slvar,qwvar,slqw are on theta-level K
  ! tke is on theta-level K-1
  ! exner is on theta-level K-1

  ! For var_diags_opt = original_vars we use the following variables:
  ! dsldzm is on theta-level K-1
  ! dqwdzm is on theta-level K-1
  ! elm is on theta-level K-1
  ! rhokh is on rho-level K

  ! For var_diags_opt = split_tke_and_inv
  ! dsldz is on rho-level K
  ! dqwdz is on rho-level K
  ! ftl is on rho-level k
  ! rhokm is on theta-level K-1

  ! Sub divide tdims domain len over the number of threads.
  tdims_omp_block = ceiling(real(tdims%i_len)/max_threads)
  ! Pick the smallest option. 
  ! Work should be split accross the threads, this provides a cap.
  ! In the occurence where the tdims domain len is very small, it will be seleted.
  ! Otherwise bl_segment_size will be selected.
  tdims_seg_block = min(bl_segment_size, tdims_omp_block, tdims%i_len)
  ! level 2 needs special treatment because of the surface

!$OMP PARALLEL DEFAULT(SHARED)                                                 &
!$OMP private(k,j,i,ii,km,kp,sh,exner,sgm,qsw_arr,weight1,weight2,weight3,     &
!$OMP         var_fac,sl_var,qw_var,sl_qw,delta_x,                             &
!$OMP         seg_slice_start, seg_slice_end, bl_segment_range)
  k = 2
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do ii = tdims%i_start, tdims%i_end, tdims_seg_block
      seg_slice_start  = ii
      seg_slice_end = MIN(((ii+tdims_seg_block)-1), tdims%i_end)
      bl_segment_range = (seg_slice_end - seg_slice_start) + 1

      if ( l_mr_physics ) then
        call qsat_wat_mix(qsw_arr(seg_slice_start:seg_slice_end, j, k-1),       &
                          tl(seg_slice_start:seg_slice_end,j,k-1),              &
                          p_theta_levels(seg_slice_start:seg_slice_end,j,k-1),  &
                          bl_segment_range )
      else
        call qsat_wat(qsw_arr(seg_slice_start:seg_slice_end, j, k-1),            &
                      tl(seg_slice_start:seg_slice_end,j,k-1),                   &
                      p_theta_levels(seg_slice_start:seg_slice_end,j,k-1),       &
                      bl_segment_range )
      end if

    end do ! ii

    if (var_diags_opt == original_vars) then

      do i = tdims%i_start, tdims%i_end
        ! calculate sh, don't interpolate with the surface flux or divide by 0
        if ( BL_diag%tke(i,j,k) > 1.0e-10_r_bl ) then
          sh = rhokh(i,j,k)                                                    &
               / ( rho_wet_tq(i,j,k-1) * elm(i,j,k)                            &
               * sqrt( 2.0_r_bl * BL_diag%tke(i,j,k) ) )
        else
          sh = small_sh
        end if
        if ( sh < small_sh ) sh = small_sh
        ! calculate exner
        exner = ( p_theta_levels(i,j,k-1) / pref )**kappa
        ! calculate the variance, use gradient interpolated between levs 1 and 2
        sgm(i) = a_dqsdt(i,j,k-1)**2 * exner**2 * b2 * sh                      &
              * elm(i,j,k)**2 * dsldz(i,j,k)**2                                &
            + a_qs(i,j,k-1)**2 * b2 * sh                                       &
              * elm(i,j,k)**2 * dqwdz(i,j,k)**2                                &
            - 2.0_r_bl * a_qs(i,j,k-1) * a_dqsdt(i,j,k-1) * exner * b2         &
              * sh * elm(i,j,k)**2 * dsldz(i,j,k) * dqwdz(i,j,k)
        if (BL_diag%l_slvar) then
          BL_diag%slvar(i,j,k-1) = b2 * sh * elm(i,j,k)**2 * dsldz(i,j,k)**2
        end if
        if (BL_diag%l_qwvar) then
          BL_diag%qwvar(i,j,k-1) = b2 * sh * elm(i,j,k)**2 * dqwdz(i,j,k)**2
        end if
        if (BL_diag%l_slqw) then
          BL_diag%slqw(i,j,k-1) = b2 * sh * elm(i,j,k)**2 * dsldz(i,j,k)       &
                                     * dqwdz(i,j,k)
        end if
        ! do this for safety, not sure if it's really needed
        sgm(i) = sqrt ( max( sgm(i), zero ) )
      end do

    else ! var_diags_opt

      do i = tdims%i_start, tdims%i_end
        ! calculate the variance
        sl_var = zero
        qw_var = zero
        sl_qw = zero
        if ( BL_diag%tke(i,j,k) > small_tke ) then
          ! vertical interpolation weights
          weight1 = r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)
          weight2 = r_theta_levels(i,j,k-1) - r_rho_levels(i,j,k-1)
          weight3 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
          ! var_fac=b2*L/sqrt(TKE)
          var_fac = b2 * rhokm(i,j,k) / ( weight1 * BL_diag%tke(i,j,k) *       &
                                    rho_wet_tq(i,j,k-1) * rho_mix_tq(i,j,k-1) )
          ! Note that flux*gradient can be negative so the absolute values
          ! are used
          qw_var= abs( -var_fac*( weight2 * fqw(i,j,k)  * dqwdz(i,j,k) +       &
                                  weight3 * fqw(i,j,k-1)* dqwdz(i,j,k-1) ) )
          sl_var= abs( -var_fac*( weight2 * ftl(i,j,k)  * dsldz(i,j,k) +       &
                                  weight3 * ftl(i,j,k-1)* dsldz(i,j,k-1) ) )
          sl_qw = - one_half * var_fac*(                                       &
           weight2*( ftl(i,j,k)*dqwdz(i,j,k) + fqw(i,j,k)*dsldz(i,j,k) ) +     &
           weight3*( ftl(i,j,k-1)*dqwdz(i,j,k-1) + fqw(i,j,k-1)*dsldz(i,j,k-1))&
                                  )
          if (BL_diag%l_slvar) then
            BL_diag%slvar(i,j,k-1) = sl_var
          end if
          if (BL_diag%l_qwvar) then
            BL_diag%qwvar(i,j,k-1) = qw_var
          end if
          if (BL_diag%l_slqw) then
            BL_diag%slqw(i,j,k-1)  = sl_qw
          end if
        end if
        exner = ( p_theta_levels(i,j,k-1) / pref )**kappa
        sgm(i) = a_dqsdt(i,j,k-1)**2 * exner**2 * sl_var                       &
               + a_qs(i,j,k-1)**2 * qw_var                                     &
               - 2.0_r_bl * a_qs(i,j,k-1) * a_dqsdt(i,j,k-1) * exner * sl_qw
        !  do this for safety, not sure if it's really needed
        sgm(i) = sqrt ( max( sgm(i), zero ) )

      end do !i
    end if ! var_diags_opt

    if (i_rhcpt == rhcpt_tke_based) then
      do i = tdims%i_start, tdims%i_end
        ! calculate rhcrit, with appropriate limits
        ! calculate grid-box size, just take surface for simplicity
        delta_x = sqrt( r_theta_levels(i,j,k-1) * delta_lambda *               &
                      r_theta_levels(i,j,k-1) * delta_phi *                    &
                      cos_theta_latitude(i,j) )
        ! max limit, based on curve fitted to aircraft observations
        max_rhcpt(i,j) = min( 0.99_r_bl, 0.997_r_bl - 0.0078_r_bl *            &
                                  log( 0.001_r_bl * delta_x ) )
        ! min limit, based on curve fitted to aircraft observations
        min_rhcpt(i,j) = max( 0.6_r_bl, 0.846_r_bl - 0.065_r_bl *              &
                                 log( 0.001_r_bl * delta_x ) )
        ! full expression
        rhcpt(i,j,k-1) = min( max_rhcpt(i,j), max( min_rhcpt(i,j),             &
                      one - root6 * sgm(i) / (a_qs(i,j,k-1) * qsw_arr(i,j,k-1))))
      end do !i
    end if
  end do   !j
!$OMP end do

  ! remaining bl levels use proper interpolation

  do k = 3, bl_levels-1
!$OMP do SCHEDULE(STATIC)
    do j = tdims%j_start, tdims%j_end
      do ii = tdims%i_start, tdims%i_end, tdims_seg_block
        seg_slice_start  = ii
        seg_slice_end = MIN(((ii+tdims_seg_block)-1), tdims%i_end)
        bl_segment_range = (seg_slice_end - seg_slice_start) + 1

        if ( l_mr_physics ) then
          call qsat_wat_mix(qsw_arr(seg_slice_start:seg_slice_end, j, k-1),      &
                            tl(seg_slice_start:seg_slice_end,j,k-1),             &
                            p_theta_levels(seg_slice_start:seg_slice_end,j,k-1), &
                            bl_segment_range)
        else
          call qsat_wat(qsw_arr(seg_slice_start:seg_slice_end, j, k-1),          &
                        tl(seg_slice_start:seg_slice_end,j,k-1),                 &
                        p_theta_levels(seg_slice_start:seg_slice_end,j,k-1),     &
                        bl_segment_range)
        end if

      end do ! ii

      if (var_diags_opt == original_vars) then

        do i = tdims%i_start, tdims%i_end
          ! calculate sh, don't divide by 0
          if ( BL_diag%tke(i,j,k) > 1.0e-10_r_bl ) then
            weight1 = r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)
            weight2 = r_theta_levels(i,j,k-1)-r_rho_levels(i,j,k-1)
            weight3 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            sh = ( weight2 * rhokh(i,j,k) + weight3 * rhokh(i,j,k-1) )         &
                 / ( weight1 * rho_wet_tq(i,j,k-1) * elm(i,j,k)                &
                 * sqrt( 2.0_r_bl * BL_diag%tke(i,j,k) ) )
          else
            sh = small_sh
          end if
          if ( sh < small_sh ) sh = small_sh
          ! calculate exner
          exner = ( p_theta_levels(i,j,k-1) / pref )**kappa
          ! calculate the variance
          sgm(i) = a_dqsdt(i,j,k-1)**2 * exner**2 * b2 * sh                    &
                * elm(i,j,k)**2 * dsldzm(i,j,k)**2                             &
              + a_qs(i,j,k-1)**2 * b2 * sh                                     &
                * elm(i,j,k)**2 * dqwdzm(i,j,k)**2                             &
              - 2.0_r_bl * a_qs(i,j,k-1) * a_dqsdt(i,j,k-1) * exner * b2       &
                * sh * elm(i,j,k)**2 * dsldzm(i,j,k) * dqwdzm(i,j,k)
          if (BL_diag%l_slvar) then
            BL_diag%slvar(i,j,k-1) = b2 * sh * elm(i,j,k)**2 * dsldz(i,j,k)**2
          end if
          if (BL_diag%l_qwvar) then
            BL_diag%qwvar(i,j,k-1) = b2 * sh * elm(i,j,k)**2 * dqwdz(i,j,k)**2
          end if
          if (BL_diag%l_slqw) then
            BL_diag%slqw(i,j,k-1) = b2 * sh * elm(i,j,k)**2 * dsldzm(i,j,k)    &
                                       * dqwdz(i,j,k)
          end if
        end do !i
        if (i_rhcpt == rhcpt_tke_based) then
          do i = tdims%i_start, tdims%i_end
            ! do this for safety, not sure if it's really needed
            sgm(i) = sqrt ( max( sgm(i), zero ) )
            ! calculate rhcrit, with appropriate limits
            rhcpt(i,j,k-1) = min( max_rhcpt(i,j), max( min_rhcpt(i,j),         &
                      one - root6 * sgm(i) / (a_qs(i,j,k-1) * qsw_arr(i,j,k-1))))
          end do !i
        end if

      else ! var_diags_opt

        do i = tdims%i_start, tdims%i_end
          ! calculate the variance
          sl_var = zero
          qw_var = zero
          sl_qw  = zero
          if ( BL_diag%tke(i,j,k) > small_tke ) then
            ! vertical interpolation weights
            weight1 = r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)
            weight2 = r_theta_levels(i,j,k-1) - r_rho_levels(i,j,k-1)
            weight3 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            var_fac = b2 * rhokm(i,j,k) / ( weight1*BL_diag%tke(i,j,k) *       &
                                     rho_wet_tq(i,j,k-1) * rho_mix_tq(i,j,k-1))
            kp=k
            km=k-1
            ! Don't use level ntml (or ntdsc) if disc_inv=2 as this indicates
            ! the inversion has just risen and the gradients between ntml and
            ! ntml-1 are likely to give excessive variances
            if ( (kp == ntml(i,j)  .and. sml_disc_inv(i,j) == 2) .or.          &
                 (kp == ntdsc(i,j) .and. dsc_disc_inv(i,j) == 2) ) kp = km
            if ( (km == ntml(i,j)  .and. sml_disc_inv(i,j) == 2) .or.          &
                 (km == ntdsc(i,j) .and. dsc_disc_inv(i,j) == 2) ) km = kp
            ! Note that flux*gradient can be negative so the absolute values
            ! are used
            qw_var= abs( -var_fac*( weight2 * fqw(i,j,kp) * dqwdz(i,j,kp) +    &
                                    weight3 * fqw(i,j,km) * dqwdz(i,j,km) ) )
            sl_var= abs( -var_fac*( weight2 * ftl(i,j,kp) * dsldz(i,j,kp) +    &
                                    weight3 * ftl(i,j,km) * dsldz(i,j,km) ) )
            sl_qw = - one_half * var_fac*(                                     &
           weight2*( ftl(i,j,k)*dqwdz(i,j,k) + fqw(i,j,k)*dsldz(i,j,k) ) +     &
           weight3*( ftl(i,j,k-1)*dqwdz(i,j,k-1) + fqw(i,j,k-1)*dsldz(i,j,k-1))&
                                    )
            if (BL_diag%l_slvar) then
              BL_diag%slvar(i,j,k-1) = sl_var
            end if
            if (BL_diag%l_qwvar) then
              BL_diag%qwvar(i,j,k-1) = qw_var
            end if
            if (BL_diag%l_slqw) then
              BL_diag%slqw(i,j,k-1) = sl_qw
            end if
          end if
          exner = ( p_theta_levels(i,j,k-1) / pref )**kappa
          sgm(i) = a_dqsdt(i,j,k-1)**2 * exner**2 * sl_var                     &
                 + a_qs(i,j,k-1)**2 * qw_var                                   &
                 - 2.0_r_bl * a_qs(i,j,k-1) * a_dqsdt(i,j,k-1) * exner * sl_qw
        end do !i
        if (i_rhcpt == rhcpt_tke_based) then
          do i = tdims%i_start, tdims%i_end
            ! do this for safety, not sure if it's really needed
            sgm(i) = sqrt ( max( sgm(i), zero ) )
            ! calculate rhcrit, with appropriate limits
            rhcpt(i,j,k-1) = min( max_rhcpt(i,j), max( min_rhcpt(i,j),         &
                      one - root6 * sgm(i) / (a_qs(i,j,k-1) * qsw_arr(i,j,k-1))))
          end do !i
        end if

      end if ! var_diags_opt

    end do   !j
!$OMP end do NOWAIT
  end do     !k
!$OMP BARRIER

  if (i_rhcpt == rhcpt_tke_based) then
    ! just use rhcpt=0.8 above bl_levels
!$OMP do SCHEDULE(STATIC)
    do k = bl_levels-1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          rhcpt(i,j,k) = 0.8_r_bl
        end do
      end do
    end do
!$OMP end do
  end if
!$OMP end PARALLEL
end if !i_rhcpt or variance diagnostics

! ------ liquid potential temperature gradient for bi-modal cloud scheme -------

! Initialise values
if ( i_cld_vn == i_cld_bimodal .or.                                            &
    (i_cld_vn == i_cld_pc2 .and. i_pc2_init_method == pc2init_bimodal) ) then

!$OMP PARALLEL DEFAULT(none) private(k,j,i)                                    &
!$OMP SHARED(tdims,tgrad_bm,bl_levels,kplume,dsldzm,z_tq,zhnl,zh,zhsc,         &
!$OMP        i_bm_ez_opt,ez_max_bm)
!$OMP do SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        ! Initialise liquid potential temperature gradient
        tgrad_bm(i,j,k) = bm_negative_init
      end do
    end do
  end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels-1
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        ! Calculate liquid potential temperature gradient for bimodal cloud
        ! scheme entrainment zone identification. Only calculate this gradient
        ! if the level is within the top half of the boundary layer and resides
        ! below 3000 m (low cloud), or if the level is  within ez_max_bm above
        ! the top of the mixed layer (defined by zh or zhsc). Else keep the
        ! gradient at its initialised zero-value. This change avoids unrealistc
        ! free-tropospheric deep entrainment zones. Note that dsldzm is held at
        ! k+1
        if ( k-1 > kplume(i,j) .and. k > 2 .and.                               &
             z_tq(i,j,k-1) > one_half*zhnl(i,j) .and.                          &
             (i_bm_ez_opt==i_bm_ez_subcrit .or.                                &
             (z_tq(i,j,k-1) < 3.0e3_r_bl) .or.                                 &
             ((z_tq(i,j,k-1)-zh(i,j)) < ez_max_bm) .or.                        &
             ((z_tq(i,j,k-1)-zhsc(i,j)) < ez_max_bm))) then
          tgrad_bm(i,j,k-1) = dsldzm(i,j,k)
        end if
      end do
    end do
  end do
!$OMP end do
!$OMP end PARALLEL

end if ! i_cld_vn == i_cld_bimodal

!If autotuning is active, decide what to do with the
!trial segment size and report the current status.
call timer('bdy_expl2')
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bdy_expl2
end module bdy_expl2_mod
