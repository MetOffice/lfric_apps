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

!  Documentation: UMDP 25.

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module bdy_expl2_1a_mod

use UM_ParCore, only: parcore_mype => mype
use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'BDY_EXPL2_1A_MOD'
contains

subroutine bdy_expl2_1a (                                                      &
! IN values defining vertical grid of model atmosphere :
 bl_levels,p_theta_levels,land_pts,land_index,                                 &
! IN U, V and W momentum fields.
 u_p,v_p,u_0_px,v_0_px,                                                        &
! IN variables for TKE scheme
 pstar,p_rho_levs,                                                             &
! IN from other part of explicit boundary layer code
 rho_mix,rho_wet_tq,rdz,rdz_charney_grid,                                      &
 z_tq,z_uv,bt,bt_gb,bq_gb,                                                     &
 flandg,rib_gb, sil_orog_land, z0m_eff_gb,                                     &
! IN cloud/moisture data :
 q,qcf,qcl,t,qw,tl,                                                            &
! IN everything not covered so far :
 fb_surf,u_s,                                                                  &
 zh_prev,ho2r2_orog,sd_orog,                                                   &
! 2 IN for Smagorinsky
  delta_smag, shear,                                                           &
! stash diagnostics
 BL_diag,                                                                      &
! INOUT variables
 zh,ntml,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,                          &
! INOUT variables on TKE based turbulence schemes
 e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                                 &
! OUT new variables for message passing
 tau_fd_x, tau_fd_y, visc_m, visc_h, rhogamu, rhogamv,                         &
! OUT Diagnostic not requiring STASH flags :
 shallowc,cu_over_orog,                                                        &
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,bl_type_7,        &
! OUT data required for tracer mixing :
 kent, we_lim, t_frac, zrzi, kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,       &
! OUT data required elsewhere in UM system :
 zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0,leonard_kl_tke                           &
     )

use atm_fields_bounds_mod, only: pdims, tdims, tdims_l,                        &
    pdims_s
use bl_option_mod, only: t_drain, h_scale, sg_orog_mixing, local_fa,           &
      free_trop_layers, smooth_to_bdys, one_third, sg_shear,                   &
      sg_shear_enh_lambda
use bl_diags_mod, only: strnewbldiag
use cv_run_mod, only: l_param_conv
use gen_phys_inputs_mod, only: l_mr_physics
use jules_surface_mod, only: formdrag, explicit_stress
use mym_option_mod, only:                                                      &
   bdy_tke, deardorff, mymodel25, mymodel3, tke_levels,                        &
   l_local_above_tkelvs, l_3dtke
use planet_constants_mod, only: cp, g, vkman
use turb_diff_mod, only:                                                       &
    l_subfilter_vert, l_subfilter_horiz, mix_factor,                           &
    turb_startlev_vert, turb_endlev_vert
use umPrintMgr, only: umPrint, umMessage
use water_constants_mod, only: lc

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

use ddf_ctl_mod, only: ddf_ctl
use ex_coef_mod, only: ex_coef
use mym_ctl_mod, only: mym_ctl
use mym_ex_flux_tq_mod, only: mym_ex_flux_tq
use fm_drag_mod, only: fm_drag

implicit none

!  Inputs :-
integer, intent(in) ::                                                         &
 land_pts,                                                                     &
                             ! No.of land points in whole grid.
 bl_levels
                             ! IN Max. no. of "boundary" levels

!     Declaration of new BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

real(kind=r_bl), intent(in) ::                                          &
 p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
                0:bl_levels+1),                                                &
 rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels+1),                                                         &
                                 ! IN density on UV (ie. rho) levels;
                                 !    used in RHOKH so dry density if
                                 !    L_mr_physics is true
 rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                                 ! IN density on TQ (ie. theta) levels;
                                 !    used in RHOKM so wet density
 rdz( pdims_s%i_start:pdims_s%i_end,                                           &
      pdims_s%j_start:pdims_s%j_end, bl_levels ),                              &
                                 ! IN RDZ(,1) is the reciprocal of
                                 !    the height of level 1, i.e. of
                                 !    the middle of layer 1.  For
                                 !    k > 1, RDZ(,k) is the
                                 !    reciprocal of the vertical
                                 !    distance from level k-1 to
                                 !    level k.
 rdz_charney_grid(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                  bl_levels),                                                  &
                                 ! IN RDZ(,1) is the reciprocal of
                                 !       the height of level 1,
                                 !       i.e. of the middle of layer 1
                                 !       For k > 1, RDZ(,k) is the
                                 !       reciprocal of the vertical
                                 !       distance from level k-1 to
                                 !       level k.
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                                 ! IN Z_tq(*,k) is height of full
                                 !    level k.
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1),        &
                                  ! OUT Z_uv(*,k) is height of half
                                  ! level k-1/2.
 u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),           &
                                 ! IN U on P-grid.
 v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),           &
                                 ! IN V on P-grid.
 bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),            &
                                 ! IN A buoyancy parameter for clear
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                                 ! IN A grid-box mean buoyancy param
                                 ! on p,T,q-levels (full levels).
 bq_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                 ! IN A grid-box mean buoyancy param
                                 ! on p,T,q-levels (full levels).

real(kind=r_bl), intent(in) ::                                          &
 flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),          &
                                 ! IN Land fraction on all tiles
 p_rho_levs(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,       &
            pdims_s%k_start:bl_levels+1),                                      &
                              ! IN p_rho_levs(*,k) is pressure at half
                              ! level k-1/2.
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! IN Surface pressure (Pascals).

! (f) Atmospheric + any other data not covered so far, incl control.

real(kind=r_bl), intent(in) ::                                          &
 fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
                                  ! IN Surface flux buoyancy over
                                  ! density (m^2/s^3)

 u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                     &
                                  ! IN Surface friction velocity
                                  !    (m/s)
 zh_prev(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
                                  ! IN boundary layer height from
                                  !    previous timestep
 rib_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
                               ! IN  Bulk Richardson number for lowest
                               ! layer
 sil_orog_land(land_pts),                                                      &
                               ! IN Silhouette area of unresolved
                               ! orography per unit horizontal area
 delta_smag(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                                 ! IN delta_x used by Smagorinsky
 shear(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,bl_levels)
                                 ! IN 3D Wind shear parameter

real(kind=r_bl), intent(in) ::                                          &
 u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),          &
                                 ! IN W'ly component of surface
!                                       current (m/s). P grid
   v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
                                   ! IN S'ly component of surface
!                                       current (m/s). P grid
   ho2r2_orog(land_pts),                                                       &
                                   ! IN peak to trough height of
!                                       unresolved orography
!                                       on land points only (m)
   sd_orog(land_pts),                                                          &
                                   ! IN Standard Deviation of unresolved
!                                       orography on land points only (m)
   z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! IN Effective grid-box roughness
!                                 length for momentum

integer, intent(in) ::                                                         &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.
! (e) Cloud data.
real(kind=r_bl), intent(in) ::                                          &
 qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                                   ! IN Cloud ice (kg per kg air)
 qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,              &
     tdims_l%k_start:bl_levels),                                               &
                                   ! IN Cloud liquid water
 q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,                &
   tdims_l%k_start:bl_levels),                                                 &
                                   ! IN specific humidity
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),             &
                                   ! IN temperature
 qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),           &
                                 ! IN Total water content
 tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)
                                 ! IN Ice/liquid water temperature

! INOUT variables
real(kind=r_bl), intent(in out) ::                                      &
 zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                      &
                                 ! INOUT Height above surface of top
                                 !       of boundary layer (metres).
 fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),           &
                                 ! INOUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                                   ! INOUT FTL(,k) contains net turbulent
!                                     sensible heat flux into layer k
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                   ! INOUT Exchange coeffs for moisture.

real(kind=r_bl), intent(in out) ::                                      &
 rhokm(pdims_s%i_start:pdims_s%i_end,                                          &
       pdims_s%j_start:pdims_s%j_end ,bl_levels)
                              ! Exchange coefficients for momentum on P-grid
! INOUT but not used: variables used in the 1A version (TKE-based schemes)
real(kind=r_bl), intent(in out) ::                                      &
  e_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
      bl_levels),                                                              &
  tsq_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
      bl_levels),                                                              &
  qsq_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
      bl_levels),                                                              &
  cov_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
      bl_levels),                                                              &
  zhpar_shcu(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

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

!  Outputs :-
!  (a) Calculated anyway (use STASH space from higher level) :-
real(kind=r_bl), intent(out) ::                                         &
 visc_m(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,bl_levels),&
                              ! Diffusion coefficient for momentum
 visc_h(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,bl_levels),&
                              ! Diffusion coefficient for heat and moisture
 rhogamu(pdims_s%i_start:pdims_s%i_end,                                        &
         pdims_s%j_start:pdims_s%j_end,2:bl_levels),                           &
                  ! Counter gradient terms for u
                  ! defined at theta level k-1
 rhogamv(pdims_s%i_start:pdims_s%i_end,                                        &
         pdims_s%j_start:pdims_s%j_end,2:bl_levels),                           &
                  ! Counter gradient terms for v
                  ! defined at theta level k-1
 tau_fd_x(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,         &
          bl_levels),                                                          &
 tau_fd_y(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,         &
          bl_levels),                                                          &
  bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! OUT Indicator set to 1.0 if stable
                                 !     b.l. diagnosed, 0.0 otherwise.
  bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! OUT Indicator set to 1.0 if Sc over
                                 !     stable surface layer diagnosed,
                                 !     0.0 otherwise.
  bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! OUT Indicator set to 1.0 if well
                                 !     mixed b.l. diagnosed,
                                 !     0.0 otherwise.
  bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! OUT Indicator set to 1.0 if
                                 !     decoupled Sc layer (not over
                                 !     cumulus) diagnosed,
                                 !     0.0 otherwise.
  bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! OUT Indicator set to 1.0 if
                                 !     decoupled Sc layer over cumulus
                                 !     diagnosed, 0.0 otherwise.
  bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! OUT Indicator set to 1.0 if a
                                 !     cumulus capped b.l. diagnosed,
                                 !     0.0 otherwise.
  bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! OUT Indicator set to 1.0 if a
                                 !     Shear-dominated unstable b.l.
                                 !     diagnosed, 0.0 otherwise.

real(kind=r_bl), intent(out) ::                                         &
  wstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                                 ! OUT Convective velocity scale (m/s)
  wthvs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                                 ! OUT surface flux of thv (Km/s)
  shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                                 ! OUT Shallow Cu diagnostic
                                 !   Indicator set to 1.0 if shallow,
                                 !   0.0 if not shallow or not cumulus
  cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! OUT Indicator for cumulus
                                 !     over steep orography
                                 !   Indicator set to 1.0 if true,
                                 !   0.0 if false. Exclusive.
  we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),               &
                                  ! OUT rho*entrainment rate implied b
                                  !     placing of subsidence
  zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),                 &
                                  ! OUT (z-z_base)/(z_i-z_base)
  t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),               &
                                  ! OUT a fraction of the timestep
  we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),           &
                                  ! OUT rho*entrainment rate implied b
                                  !     placing of subsidence
  zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),             &
                                  ! OUT (z-z_base)/(z_i-z_base)
  t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),           &
                                  ! OUT a fraction of the timestep
  zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! OUT Top of decoupled layer

integer, intent(out) ::                                                        &
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                                 ! OUT Top level for turb mixing in
!                                           any decoupled Sc layer
   nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                                   ! OUT Bottom level of any decoupled
                                   !     turbulently-mixed Sc layer.
   kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                                    ! OUT grid-level of SML inversion
   kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! OUT grid-level of DSC inversion

!-2 Genuinely output, needed by other atmospheric routines :-
real(kind=r_bl), intent(out) ::                                         &
  uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                           ! OUT U-component of surface wind stress
                           !     on P-grid
  vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                           ! OUT V-component of surface wind stress
                           !     on P-grid
real(kind=r_bl), intent(out) ::                                                &
         leonard_kl_tke( tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                         1:bl_levels, 2 )
!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

! Parameters also passed to EX_COEF
! Layer interface K_LOG_LAYR-1/2 is the highest which requires log
! profile correction factors to the vertical finite differences.
! The value should be reassessed if the vertical resolution is changed.
! We could set K_LOG_LAYR = BL_LEVELS and thus apply the correction
! factors for all the interfaces treated by the boundary layer scheme;
! this would be desirable theoretically but expensive computationally
! because of the use of the log function.
integer, parameter ::    k_log_layr = 2
!-----------------------------------------------------------------------
!  Workspace :-
real(kind=r_bl) ::                                                      &
 dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                     &
      2:bl_levels),                                                            &
                              ! Buoyancy gradient across layer
                              !  interface.
 dvdzm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                    &
       2:bl_levels),                                                           &
                              ! Modulus of wind shear.
 rmlmax2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
                              ! Square of asymptotic mixing length
                              ! for Smagorinsky scheme
 ri(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),          &
                              ! Local Richardson number.
 rhokh_th_ri(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
          bl_levels),                                                          &
 rhokm_ri(pdims_s%i_start:pdims_s%i_end,                                       &
            pdims_s%j_start:pdims_s%j_end ,bl_levels),                         &
                              ! Exchange coefficients for momentum and
                              ! heat on theta-levels as calculated by
                              ! the local Ri-based scheme
 weight_1dbl(pdims%i_start:pdims%i_end,                                        &
             pdims%j_start:pdims%j_end ,bl_levels),                            &
                              ! Weighting applied to 1D BL scheme
                              ! to blend with Smagorinsky scheme,
                              ! index k held on theta level (k-1)
 weight_1dbl_rho(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end,bl_levels),                         &
                              ! weight_1dbl interpolated to rho levels
 elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),         &
                              ! Mixing length for momentum as
                              ! calculated by the Ri-based scheme
 elh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),         &
 elh_rho(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                              ! Mixing length for heat (m),
                              ! held on theta and rho levels, resp.
 tke_loc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
            2:bl_levels),                                                      &
                              ! Ri-based scheme diagnosed TKE
 fm_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                              ! stability function for momentum transport
                              ! level 1 value is dummy
 fh_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                              ! stability function for heat and moisture.
                              ! level 1 value is dummy
 sigma_h(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
                              ! Standard deviation of subgrid
                              ! orography (m) [= 2root2 * ho2r2_orog]
 p_half(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),        &
 rneutml_sq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                              ! Square of the neutral mixing length scale

real(kind=r_bl), allocatable :: visc_h_rho (:,:,:)
                                                       ! visc_h on rho levels

real(kind=r_bl) ::                                                      &
   zh_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                              ! Height above surface of top of
                              !  boundary layer (metres) as
                              !  determined from the local
                              !  Richardson number profile.
   zhnl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                  &
                              ! non-local PBL depth
   zdsc_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                              ! Height of base of K_top in DSC
   dtldz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
      2:bl_levels),                                                            &
                              ! TL+gz/cp gradient between
                              ! levels k and k-1
   dqwdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
      2:bl_levels),                                                            &
                              ! QW gradient between
   dtldzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          2:bl_levels),                                                        &
                  ! gradient of TL across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,k) represents the value on theta level k-1
   dqwdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          2:bl_levels),                                                        &
                  ! gradient of QW across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,k) represents the value on theta level k-1
   dudz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:bl_levels),                                                          &
                  ! Gradient of u at theta levels.
                  !(:,:,k) represents the value on theta level k-1
   dvdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:bl_levels)
                  ! Gradient of v at theta levels.
                  !(:,:,k) repserents the value on theta level k-1

integer ::                                                                     &
 ntml_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! Number of model layers in the
!                                    turbulently mixed layer as
!                                    determined from the local
!                                    Richardson number profile.
   ntml_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! Number of model layers in the
!                                    turbulently mixed layer as
!                                    determined from the parcel ascent.

logical ::                                                                     &
 unstable(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                               ! Logical switch for unstable
                               !    surface layer.
 dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! Flag set if decoupled
                               ! stratocumulus layer found

real(kind=r_bl) ::                                                      &
   rhogamt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels),                                                       &
                  ! Counter gradient terms for TL
                  ! defined at rho levels
   rhogamq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                &
           2:bl_levels)
                  ! Counter gradient terms for QW
                  ! defined at rho levels

real(kind=r_bl) ::                                                      &
 lambda_min
            ! Min value of length scale LAMBDA.

!  Local scalars :-
real(kind=r_bl) ::                                                      &
   weight1,                                                                    &
   weight2,                                                                    &
   weight3,                                                                    &
   r_weight1,                                                                  &
   zpr,                                                                        &
             ! z/sigma_h
   slope,                                                                      &
             ! subgrid orographic slope
   grcp      ! G/CP

integer  ::                                                                    &
   i,j,                                                                        &
                     ! LOCAL Loop counter (horizontal field index).
   k,ient,                                                                     &
                     ! LOCAL Loop counter (vertical level index).
   l
! LOCAL Loop counter for land points

real(kind=r_bl), parameter :: max_ri = 0.01_r_bl*sqrt(huge(1.0_r_bl))
                 ! Maximum (absolute) Richardson number which ensures that
                 ! the stability functions (~ri^2) remain real-valued at
                 ! the given model precision

character(len=*), parameter ::  RoutineName = 'BDY_EXPL2_1A'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Parameter check
! error checking here moved to readsize/scm_shell

! set pressure array.
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    p_half(i,j,1) = pstar(i,j)
  end do
end do
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      p_half(i,j,k) = p_rho_levs(i,j,k)
    end do
  end do
end do  ! end of loop over bl_levels

!-----------------------------------------------------------------------
if (formdrag ==  explicit_stress) then
  !------------------------------------------------------------------
  !      Calculate stress profiles
  !------------------------------------------------------------------
  call fm_drag (                                                               &
  ! IN levels
        land_pts, land_index, bl_levels,                                       &
  ! IN fields
        u_p, v_p, tl, qw, bt_gb, bq_gb, rho_wet_tq,                            &
        z_uv, z_tq, z0m_eff_gb, zh_prev, rib_gb, sil_orog_land,                &
  ! OUT fields
        tau_fd_x(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
                 1:bl_levels),                                                 &
        tau_fd_y(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
                 1:bl_levels)                                                  &
        )
  !------------------------------------------------------------------
  !      Orographic stress diagnostics
  !------------------------------------------------------------------
  if (BL_diag%l_ostressx) then
    do k = 1, bl_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          BL_diag%ostressx(i,j,k)=tau_fd_x(i,j,k)
        end do
      end do
    end do
  end if
  if (BL_diag%l_ostressy) then
    do k = 1, bl_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          BL_diag%ostressy(i,j,k)=tau_fd_y(i,j,k)
        end do
      end do
    end do
  end if

end if

!------------------------------------------------------------------
!  Initialize weighting applied to 1d BL scheme
!  (used to blend 1D with 3D Smagorinsky scheme)
!------------------------------------------------------------------
do k = 1, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      weight_1dbl(i,j,k) = 1.0
      weight_1dbl_rho(i,j,k) = 1.0     ! dummy here
    end do
  end do
end do
!------------------------------------------------------------------
!  Initialize fluxes
!------------------------------------------------------------------
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ftl(i,j,k) = 0.0
      fqw(i,j,k) = 0.0
    end do
  end do
end do
!-------------------------------------------------------------
! Set all variables from the non-local scheme to zero or "off"
!  - reset all fluxes and K's arising from the non-local scheme
!-------------------------------------------------------------
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    ntml_nl(i,j) = ntml(i,j)
        ! decoupled mixed layer
    dsc(i,j)     = .false.
    ntdsc(i,j)   = 0
    nbdsc(i,j)   = 0
    zhsc(i,j)    = 0.0
    zhnl(i,j)    = 0.0
    zdsc_base(i,j) = 0.0
        ! entrainment variables for non-local tracer mixing
    kent(i,j) = 2
    kent_dsc(i,j) = 2
    do ient = 1, 3
      t_frac(i,j,ient) = 0.0
      zrzi(i,j,ient)   = 0.0
      we_lim(i,j,ient) = 0.0
      t_frac_dsc(i,j,ient) = 0.0
      zrzi_dsc(i,j,ient)   = 0.0
      we_lim_dsc(i,j,ient) = 0.0
    end do
    unstable(i,j) = (fb_surf(i,j) >  0.0)
  end do
end do

! for compatibility to the original bdy_expl2
if (l_subfilter_vert) then
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      cumulus(i,j) = .false.
      l_shallow(i,j) = .false.
      ntpar(i,j)   = 0
      ntml_nl(i,j) = -1    ! to ensure correct diagnostics
    end do
  end do
end if
!-----------------------------------------------------------------------
! Calculate lapse rates
!-----------------------------------------------------------------------
grcp = g/cp
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      dtldz(i,j,k) = ( tl(i,j,k) - tl(i,j,k-1) )                               &
                                    * rdz_charney_grid(i,j,k) + grcp
      dqwdz(i,j,k) = ( qw(i,j,k) - qw(i,j,k-1) )                               &
                                    * rdz_charney_grid(i,j,k)
    end do
  end do
end do


! Calculate `buoyancy' gradient, DBDZ, on theta-levels
! NOTE: DBDZ(k) is on theta-level k-1
do k = 3, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      r_weight1 = 1.0 / (z_uv(i,j,k) - z_uv(i,j,k-1))
      weight2 = z_tq(i,j,k-1)- z_uv(i,j,k-1)
      weight3 = z_uv(i,j,k) - z_tq(i,j,k-1)
      dtldzm(i, j, k) = (weight2 * dtldz(i,j,k)                                &
                       + weight3 * dtldz(i,j,k-1)) * r_weight1
      dqwdzm(i, j, k) = (weight2 * dqwdz(i,j,k)                                &
                           + weight3 * dqwdz(i,j,k-1)) * r_weight1
      dbdz(i,j,k) = g*( bt_gb(i,j,k-1)*dtldzm(i, j, k) +                       &
                                  bq_gb(i,j,k-1)*dqwdzm(i, j, k))
    end do
  end do
end do

k = 2
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    dtldzm(i,j,k) = dtldz(i,j,k)
    dqwdzm(i,j,k) = dqwdz(i,j,k)
    dbdz(i,j,k) = g*( bt_gb(i,j,k-1)*dtldz(i,j,k) +                            &
                              bq_gb(i,j,k-1)*dqwdz(i,j,k) )
  end do
end do

!--------------------------------------------------
! Calculate modulus of shear on theta-levels
! dvdzm(k) is on theta-level(k-1)
!--------------------------------------------------
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      ! Calculation of dudz, dvdz is temporary until use of shear terms
      ! is sorted.
      dudz(i, j, k) = (u_p(i,j,k) - u_p(i,j,k-1)) * rdz(i, j, k)
      dvdz(i, j, k) = (v_p(i,j,k) - v_p(i,j,k-1)) * rdz(i, j, k)
    end do
  end do
end do

if ((.not. l_subfilter_vert) .and. (.not. l_3dtke) ) then

  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        dvdzm(i, j, k) = max ( 1.0e-12 ,                                       &
                     sqrt(dudz(i, j, k) ** 2 + dvdz(i, j, k) ** 2))
      end do
    end do
  end do

else

  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        dvdzm(i,j,k) = max( 1.0e-12 , shear(i,j,k-1) )
      end do
    end do
  end do

end if

if (l_subfilter_horiz .or. l_subfilter_vert .or. l_3dtke) then

  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      rmlmax2(i,j) = ( mix_factor * delta_smag(i,j) )**2
    end do
  end do

  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        rneutml_sq(i,j,k) = 1.0 / (                                            &
                 1.0/( vkman*(z_tq(i,j,k) + z0m_eff_gb(i,j)) )**2              &
               + 1.0/rmlmax2(i,j) )
      end do
    end do
  end do

end if
!-----------------------------------------------------------------------
! Orographic enhancement of subgrid mixing
!-----------------------------------------------------------------------
! Calculate 2D array for standard deviation of subgrid orography.
!-----------------------------------------------------------------------
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    sigma_h(i,j) = 0.0
  end do
end do
do l = 1, land_pts
  j=(land_index(l)-1)/pdims%i_end + 1
  i=land_index(l) - (j-1)*pdims%i_end
  sigma_h(i,j) =  min( sd_orog(l), 300.0 )
end do
!-----------------------------------------------------------------------
!  Enhance resolved shear through unresolved subgrid drainage flows.
!-----------------------------------------------------------------------
if (sg_orog_mixing == sg_shear .or.                                            &
    sg_orog_mixing == sg_shear_enh_lambda) then

  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end

        if (sigma_h(i,j) > 1.0 ) then
          zpr = z_tq(i,j,k-1)/sigma_h(i,j)
          ! Height dependence, to reduce effect to zero with height
          !   gives z_scale~[1,0.95,0.5,0] at zpr=[0,0.6,1,1.7]
          weight1 = 0.5*( 1.0 - tanh(4.0*(zpr-1.0) ) )

          ! Take slope ~ sd/h_scale for small sd;
          !            tends to 0.2 for large sd
          slope = 1.0 / sqrt( 25.0 + (h_scale/sigma_h(i,j))**2 )

          dvdzm(i,j,k) = max ( dvdzm(i,j,k),                                   &
                               weight1*slope*t_drain*dbdz(i,j,k) )

          if (k==2 .and. BL_diag%l_dvdzm)                                      &
            BL_diag%dvdzm(i,j,1)=weight1*slope*t_drain*dbdz(i,j,k)

        end if
      end do
    end do
  end do

end if      ! sg_orog_mixing


!------------------------------------------------------------------
!  call main subroutines
!------------------------------------------------------------------
if (bdy_tke == mymodel25 .or. bdy_tke == mymodel3) then
  call mym_ctl(                                                                &
  !in levels/switches
            bl_levels, bdy_tke, BL_diag,                                       &
  !in fields
            z_uv,z_tq, u_p, v_p, qw, tl, t, q, qcl, qcf, bq_gb, bt_gb,         &
            rho_mix, rho_wet_tq, fqw, ftl,                                     &
            dtldzm, dqwdzm, dudz, dvdz, dbdz, dvdzm, delta_smag,               &
            p_theta_levels, p_half, u_s, fb_surf, pstar,                       &
  ! inout
            e_trb, tsq_trb, qsq_trb, cov_trb, rhokm, rhokh, zhpar_shcu,        &
  ! out
            visc_m, visc_h, rhogamu, rhogamv, rhogamt, rhogamq, leonard_kl_tke)
else if (bdy_tke == deardorff) then
  call ddf_ctl(                                                                &
  ! IN levels/switches
          bl_levels, BL_diag,                                                  &
  ! IN fields
          z_uv,z_tq, u_p, v_p, qw, tl, t, q, qcl,                              &
          qcf, p_theta_levels, p_half,bq_gb, bt_gb, rho_mix, rho_wet_tq,       &
          dtldzm, dqwdzm, dudz, dvdz, dbdz, dvdzm, delta_smag,                 &
          u_s, fb_surf, pstar,                                                 &
  ! INOUT fields
          e_trb, rhokm, rhokh, zhpar_shcu)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        rhogamu(i, j, k) = 0.0
        rhogamv(i, j, k) = 0.0
        rhogamt(i, j, k) = 0.0
        rhogamq(i, j, k) = 0.0
      end do
    end do
  end do
end if

! RHOKM and RHOKH could be changed by the subgrid turbulence
! scheme, but BL_diag%rhokm, rhokh are the exchange coefficients
! by the TKE schemes, which is the same sense in bdy_expl2 for
! the UM BL scheme.

if (BL_diag%l_rhokm) then
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokm(i,j,k)=rhokm(i,j,k)
      end do
    end do
  end do
end if

if (BL_diag%l_rhokh) then
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokh(i,j,k)=rhokh(i,j,k)
      end do
    end do
  end do
end if

!-----------------------------------------------------------------------
! The purpose of this block is to calculate local mixing above tke_levels
! and the stability functions FM_3D and FM_3H with EX_COEF.
!-----------------------------------------------------------------------
if (l_subfilter_horiz .or. l_subfilter_vert .or.                               &
        (tke_levels < bl_levels .and. l_local_above_tkelvs)) then

  ! call local coeff calculation for levels 2 to bl_levels
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        ri(i, j, k) = dbdz(i, j, k)                                            &
                      / ( dvdzm(i, j, k) * dvdzm(i ,j, k) )
        ri(i, j, k) = max(min(ri(i,j,k),max_ri),-max_ri)
      end do
    end do
  end do

  if (BL_diag%l_gradrich) then
    do k = 2, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          BL_diag%gradrich(i,j,k)=ri(i,j,k)
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------
  ! call local coeff calculation for levels 2 to bl_levels
  !-----------------------------------------------------------------------
  call ex_coef (                                                               &
  ! IN levels/logicals
       bl_levels,k_log_layr,BL_diag,                                           &
  ! IN fields
      sigma_h,flandg,dvdzm,ri,rho_wet_tq,z_uv,z_tq,z0m_eff_gb,zhnl,zhpar_shcu, &
      zhsc,zdsc_base,ntpar,ntml_nl,ntdsc,nbdsc,l_shallow,rmlmax2,rneutml_sq,   &
      delta_smag,                                                              &
  ! IN/OUT fields
      cumulus,weight_1dbl,                                                     &
  ! OUT fields
      lambda_min,zh_local,ntml_local,elm,elh,elh_rho,rhokm_ri,                 &
      rhokh_th_ri,fm_3d,fh_3d,tke_loc                                          &
    )
  !------------------------------------------------------------------
  !  set diffusion coefs between tke_levels + 1 and bl_levels
  !  with ones by the local scheme (EX_COEF)
  !------------------------------------------------------------------
  if (tke_levels < bl_levels .and. l_local_above_tkelvs) then
    do k = tke_levels + 1, bl_levels
      do j = pdims%j_start, pdims%j_end
        do i = pdims%i_start, pdims%i_end
          rhokm(i, j, k) = rhokm_ri(i, j, k)

          weight1 = z_tq(i,j,k) - z_tq(i,j, k-1)
          weight2 = z_tq(i,j,k) - z_uv(i,j,k)
          weight3 = z_uv(i,j,k) - z_tq(i,j,k-1)
          if ( k  ==  bl_levels ) then
              ! assume RHOKH_uv(BL_LEVELS+1) is zero
            rhokh(i,j,k) = ( weight2/weight1 ) * rhokh_th_ri(i,j,k)
          else
            rhokh(i,j,k) =    weight3/weight1 *                                &
                                      rhokh_th_ri(i,j,k+1)                     &
                             +weight2/weight1 *                                &
                                      rhokh_th_ri(i,j,k)
          end if

          if ((local_fa /= free_trop_layers) .and. &
              (local_fa /= smooth_to_bdys)) then
            !--------------------------------------------------------
            !  Code moved from EX_COEF to avoid interpolation:
            !  Include mixing length, ELH, in RHOKH.
            !  Here only use free trop mixing length, lambda_min
            !--------------------------------------------------------
            rhokh(i,j,k) = lambda_min * rhokh(i,j,k)
          end if   ! test on local_fa NE free_trop_layers

          ! Finally multiply RHOKH by dry density
          if (l_mr_physics) rhokh(i,j,k) = rho_mix(i,j,k) * rhokh(i,j,k)

        end do
      end do
    end do
  end if

  if (l_subfilter_horiz .or. l_subfilter_vert) then

    ! visc_m and visc_h for levels below tke_levels are set in mym_ctl.

    if (l_3dtke .and.                                                          &
        (tke_levels < bl_levels .and. l_local_above_tkelvs)) then

      do k = tke_levels, bl_levels
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            visc_m(i,j,k) = shear(i,j,k)*rneutml_sq(i,j,k)
            visc_h(i,j,k) = shear(i,j,k)*rneutml_sq(i,j,k)
          end do
        end do
      end do

      do k = tke_levels, bl_levels-1
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            ! stability functions are indexed with Ri, fm(k) on w(k-1)
            visc_m(i,j,k) = visc_m(i,j,k)*fm_3d(i,j,k+1)
            visc_h(i,j,k) = visc_h(i,j,k)*fh_3d(i,j,k+1)
          end do
        end do
      end do

    else if (.not. l_3dtke) then

      ! visc_m,h on IN are just S and visc_m,h(k) are co-located with w(k)
      do k = 1, bl_levels
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            visc_m(i,j,k) = shear(i,j,k)*rneutml_sq(i,j,k)
            visc_h(i,j,k) = shear(i,j,k)*rneutml_sq(i,j,k)
          end do
        end do
      end do

      do k = 1, bl_levels-1
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            ! stability functions are indexed with Ri, fm(k) on w(k-1)
            visc_m(i,j,k) = visc_m(i,j,k)*fm_3d(i,j,k+1)
            visc_h(i,j,k) = visc_h(i,j,k)*fh_3d(i,j,k+1)
          end do
        end do
      end do

    end if
    ! visc_m and visc _h are now lambda^2*S*FM and lambda^2*S*FH

    if (l_subfilter_vert) then

      ! visc_h_rho(k) is held on rho(k), same as BL's rhokh
      allocate (visc_h_rho(pdims%i_start:pdims%i_end,                          &
                           pdims%j_start:pdims%j_end, bl_levels))

      do k = 2, bl_levels
        do j = pdims%j_start, pdims%j_end
          do i = pdims%i_start, pdims%i_end
            weight1 = z_tq(i,j,k) - z_tq(i,j, k-1)
            weight2 = z_tq(i,j,k) - z_uv(i,j,k)
            weight3 = z_uv(i,j,k) - z_tq(i,j,k-1)
            if ( k  ==  bl_levels ) then
              ! assume visc_h(bl_levels) is zero
              ! (Ri and thence f_h not defined)
              visc_h_rho(i,j,k) = (weight2/weight1) * visc_h(i,j,k-1)
            else
              visc_h_rho(i,j,k) = (weight3/weight1) * visc_h(i,j,k)            &
                                + (weight2/weight1) * visc_h(i,j,k-1)
            end if
          end do
        end do
      end do

      ! Overwrite the diffusion coefficients from the local BL scheme
      !(RHOKM and RHOKH) with those obtained from the Smagorinsky scheme.

      do k = 2, bl_levels
        if (k >= turb_startlev_vert .and.                                      &
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
              rhokm(i,j,k) = 0.0
              rhokh(i,j,k) = 0.0
            end do
          end do
        end if
      end do

      deallocate (visc_h_rho)

    end if ! L_subfilter_vert
  end if ! L_subfilter_horiz or L_subfilter_vert
end if ! Main if-test for calling Ri-based scheme

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

! Note that this part is exactly the same as the original bdy_expl2,
! but diagnosed BL types can be only 1, 3, and 6.
!-----------------------------------------------------------------------
!      First initialise the type variables and set the diagnostic ZHT.

if (BL_diag%l_zht) then
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      bl_diag%zht(i,j) = max( zh(i,j) , zhsc(i,j) )
    end do
  end do
end if
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    bl_type_1(i,j) = 0.0
    bl_type_2(i,j) = 0.0
    bl_type_3(i,j) = 0.0
    bl_type_4(i,j) = 0.0
    bl_type_5(i,j) = 0.0
    bl_type_6(i,j) = 0.0
    bl_type_7(i,j) = 0.0
  end do
end do
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if (.not. unstable(i,j) .and. .not. dsc(i,j) .and.                         &
               .not. cumulus(i,j)) then
      !         Stable b.l.
      bl_type_1(i,j) = 1.0
    else if (.not. unstable(i,j) .and. dsc(i,j) .and.                          &
                .not. cumulus(i,j)) then
      !         Stratocumulus over a stable surface layer
      bl_type_2(i,j) = 1.0
    else if (unstable(i,j) .and. .not. cumulus(i,j) .and.                      &
                .not. dsc(i,j) ) then
      !         Well mixed b.l. (possibly with stratocumulus)
      if ( ntml(i,j)  >   ntml_nl(i,j) ) then
          ! shear-dominated - currently identified
          ! by local NTML overriding non-local
        bl_type_7(i,j) = 1.0
      else
          ! buoyancy-dominated
        bl_type_3(i,j) = 1.0
      end if
    else if (unstable(i,j) .and. dsc(i,j) .and.                                &
                                            .not. cumulus(i,j)) then
      !         Decoupled stratocumulus (not over cumulus)
      bl_type_4(i,j) = 1.0
    else if (dsc(i,j) .and. cumulus(i,j)) then
      !         Decoupled stratocumulus over cumulus
      bl_type_5(i,j) = 1.0
    else if (.not. dsc(i,j) .and. cumulus(i,j)) then
      !         Cumulus capped b.l.
      bl_type_6(i,j) = 1.0
    end if
  end do
end do
!-----------------------------------------------------------------------
! Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------
write(umMessage,*)'Entering mym_ex_flux_tq.'
call umPrint(umMessage)

call mym_ex_flux_tq(                                                           &
      bl_levels,                                                               &
      tl, qw, rhokh, rhogamt, rhogamq, rdz_charney_grid,                       &
      ftl, fqw)


if (BL_diag%l_rhogamu) then
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhogamu(i, j, k) = rhogamu(i, j, k)
      end do
    end do
  end do
end if

if (BL_diag%l_rhogamv) then
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhogamv(i, j, k) = rhogamv(i, j, k)
      end do
    end do
  end do
end if

if (BL_diag%l_rhogamt) then
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhogamt(i, j, k) = - cp * rhogamt(i, j, k)
      end do
    end do
  end do
end if

if (BL_diag%l_rhogamq) then
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhogamq(i, j, k) = - lc * rhogamq(i, j, k)
      end do
    end do
  end do
end if

!-----------------------------------------------------------------------
! Calculate explicit surface fluxes of U and V on
! P-grid for convection scheme
!-----------------------------------------------------------------------
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    uw0(i,j) = -rhokm(i,j,1) *                                                 &
                             ( u_p(i,j,1) - u_0_px(i,j) )
    vw0(i,j) = -rhokm(i,j,1) *                                                 &
                             ( v_p(i,j,1) - v_0_px(i,j) )
  end do
end do
!-----------------------------------------------------------------------
! Set NTML to max number of turbulently mixed layers
! Calculate quantities to pass to convection scheme.
!-----------------------------------------------------------------------
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    wstar(i,j) = 0.0
    wthvs(i,j) = 0.0
    cu_over_orog(i,j) = 0.0
    if ( cumulus(i,j) ) then
      if ( fb_surf(i,j)  >   0.0 ) then
        wstar(i,j) = ( zh(i,j)*fb_surf(i,j) )**one_third
        wthvs(i,j) = fb_surf(i,j) / ( g * bt(i,j,1) )
      end if
      wstar(i,j) = max( 0.1, wstar(i,j) )
      if (.not. l_param_conv) then
        ntml(i,j) = max( 2, ntml_nl(i,j) - 1 )
      end if
    else
      ntml(i,j) = max( ntml_nl(i,j) , ntdsc(i,j) )
    end if
    ! Limit explicitly calculated surface stresses
    ! to a physically plausible level.
    if ( uw0(i,j)  >=  5.0 ) then
      uw0(i,j) =  5.0
    else if ( uw0(i,j)  <=  -5.0 ) then
      uw0(i,j) = -5.0
    end if
    if ( vw0(i,j)  >=  5.0 ) then
      vw0(i,j) =  5.0
    else if ( vw0(i,j)  <=  -5.0 ) then
      vw0(i,j) = -5.0
    end if
    if (BL_diag%l_wstar .and. (fb_surf(i,j) >0.0)) then
      BL_diag%wstar(i,j)= (zh(i,j)*fb_surf(i,j))**one_third
    end if
  end do
end do

if (l_param_conv) then

  ! Check for CUMULUS having been diagnosed over steep orography.
  ! Reset to false but keep NTML at NLCL (though decrease by 2 so that
  ! coupling between BL and convection scheme can be maintained).
  ! Reset type diagnostics.

  do l = 1, land_pts
    j=(land_index(l)-1)/pdims%i_end + 1
    i=land_index(l) - (j-1)*pdims%i_end
    if (cumulus(i,j) .and. ho2r2_orog(l)  >   900.0) then
      cumulus(i,j) = .false.
      l_shallow(i,j) = .false.
      bl_type_5(i,j) = 0.0
      bl_type_6(i,j) = 0.0
      cu_over_orog(i,j) = 1.0
      if (ntml(i,j)  >=  3) ntml(i,j) = ntml(i,j) - 2
    end if
  end do

  ! Check that CUMULUS and L_SHALLOW are still consistent

  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      if ( .not. cumulus(i,j) ) l_shallow(i,j) = .false.
    end do
  end do

end if    ! (l_param_conv)
!-----------------------------------------------------------------------
!     Set shallow convection diagnostic: 1.0 if L_SHALLOW (and CUMULUS)
!                                        0.0 if .NOT. CUMULUS
!-----------------------------------------------------------------------
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    if ( cumulus(i,j) .and. l_shallow(i,j) ) then
      shallowc(i,j) = 1.0
    else
      shallowc(i,j) = 0.0
    end if
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bdy_expl2_1a
end module bdy_expl2_1a_mod
