! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! To diagnose convective occurrence and type
!
module conv_diag_comp_6a_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'CONV_DIAG_COMP_6A_MOD'
contains

subroutine conv_diag_comp_6a(                                                  &

! in values defining field dimensions and subset to be processed :
       row_length, rows                                                        &

! in level and points info
      , bl_levels, nunstable                                                   &
      , index_i, index_j                                                       &

! in grid information
      , p, P_theta_lev, exner_rho                                              &
      , rho_only, rho_theta, z_full, z_half, r_theta_levels, r_rho_levels      &

! in Cloud data :
      , qcf, qcl, cloud_fraction                                               &

! in everything not covered so far :

      , pstar, q, theta, exner_theta_levels                                    &
      , frac_land                                                              &
      , w_copy, w_max, tv1_sd, bl_vscale2                                      &
      , conv_prog_precip                                                       &

! SCM Diagnostics (dummy values in full UM)
      , nSCMDpkgs, L_SCMDiags                                                  &

! INOUT data required elsewhere in UM system :

      , ntml, ntpar, nlcl                                                      &
      , cumulus, L_shallow, l_congestus, l_congestus2 , conv_type              &
      , zh,zhpar,dzh,qcl_inv_top,z_lcl,z_lcl_uv,delthvu,ql_ad                  &
      , cin, cape, entrain_coef                                                &
      , qsat_lcl                                                               &
       )

! Modules used
! Definitions of prognostic variable array sizes
use atm_fields_bounds_mod, only:                                               &
  pdims_s, tdims_s, tdims, tdims_l


use cv_diag_param_mod, only:                                                   &
  a_plume, b_plume, max_t_grad

use cv_param_mod, only:                                                        &
   max_diag_thpert, sh_grey_closure, refqsat, ae2, entcoef

use bl_option_mod, only: bl_res_inv, off, on, zero, one, one_half

use planet_constants_mod, only:                                                &
    c_virtual => c_virtual_bl, g => g_bl, planet_radius => planet_radius_bl,   &
    lsrcp => lsrcp_bl, lcrcp => lcrcp_bl, gamma_dry => grcp_bl
use water_constants_mod, only: tm => tm_bl

use cv_run_mod, only:                                                          &
    iconv_congestus, icvdiag, cldbase_opt_dp, w_cape_limit, iconv_deep,        &
    limit_pert_opt, cvdiag_inv, cvdiag_sh_wtest, iwtcs_diag1,                  &
    iwtcs_diag2, cldbase_opt_sh, ent_fac_dp,                                   &
    prog_ent_grad, prog_ent_int, prog_ent_max, prog_ent_min,                   &
    l_reset_neg_delthvu

use cloud_inputs_mod, only: forced_cu

use column_rh_mod, only:                                                       &
    calc_column_rh

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

! subroutines
use cumulus_test_mod, only: cumulus_test
use mean_w_layer_mod, only: mean_w_layer
use parcel_ascent_mod, only: parcel_ascent



use missing_data_mod, only: rmdi, imdi
use ereport_mod,      only: ereport
use model_domain_mod, only: model_type, mt_single_column
use s_scmop_mod,      only: default_streams,                                   &
                            t_inst, d_wet, d_point,                            &
                            scmdiag_bl, scmdiag_conv
use nlsizes_namelist_mod,   only: model_levels

use errormessagelength_mod, only: errormessagelength

use cv_parcel_neutral_dil_mod, only: cv_parcel_neutral_dil
use cv_parcel_neutral_inv_mod, only: cv_parcel_neutral_inv
use lift_cond_lev_mod,         only: lift_cond_lev
implicit none

! ------------------------------------------------------------------------------
! Description:
!   To diagnose convection occurrence and type on just unstable points.
!
!   Called by conv_diag
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

! (a) Defining horizontal grid and subset thereof to be processed.

integer, intent(in) ::                                                         &
  row_length           & ! Local number of points on a row
 ,rows                   ! Local number of rows in a theta field

! (b) Defining vertical grid of model atmosphere.

integer, intent(in) ::                                                         &
  bl_levels            & ! Max. no. of "boundary" levels allowed.
 ,nunstable              ! number of unstable points

integer, intent(in) ::                                                         &
  index_i(nunstable)   & ! column number of unstable points
 ,index_j(nunstable)     ! row number of unstable points

real(kind=r_bl), intent(in) ::                                                 &
  p(pdims_s%i_start:pdims_s%i_end,    & ! pressure  on rho levels (Pa)
    pdims_s%j_start:pdims_s%j_end,                                             &
    pdims_s%k_start:pdims_s%k_end)                                             &
 ,p_theta_lev(tdims%i_start:tdims%i_end,    & ! P on theta lev (Pa)
              tdims%j_start:tdims%j_end,                                       &
                          1:tdims%k_end)                                       &
 ,exner_rho(pdims_s%i_start:pdims_s%i_end,  & ! Exner on rho level
            pdims_s%j_start:pdims_s%j_end,  & !
            pdims_s%k_start:pdims_s%k_end)                                     &
 ,exner_theta_levels(tdims%i_start:tdims%i_end, & ! exner pressure theta lev
                     tdims%j_start:tdims%j_end, & !  (Pa)
                                 1:tdims%k_end)

real(kind=r_bl), intent(in) ::                                                 &
  rho_only(row_length,rows,model_levels)    & ! density (kg/m3)
 ,rho_theta(row_length,rows,model_levels-1) & ! density th lev (kg/m3)
 , r_theta_levels(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end, &
                   0:tdims%k_end) & ! dist of theta lev from Earth centre (m)
 , r_rho_levels(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
                 tdims%k_end)     & ! dist of rho lev from Earth centre (m)
 ,z_full(row_length,rows,model_levels)      & ! height th lev (m)
 ,z_half(row_length,rows,model_levels)        ! height rho lev (m)

! (c) Cloud data.
real(kind=r_bl), intent(in) ::                                                 &
  qcf(row_length,rows,model_levels)     & ! Cloud ice (kg/kg air)
 ,qcl(row_length,rows,model_levels)     & ! Cloud liquid water (kg/kg air)
 ,cloud_fraction(row_length, rows, model_levels)   ! Cloud fraction

! (d) Atmospheric + any other data not covered so far, incl control.
real(kind=r_bl), intent(in) ::                                                 &
  pstar(row_length, rows)                & ! Surface pressure (Pascals).
 ,w_copy(row_length,rows,0:model_levels) & ! vertical velocity (m/s)
 ,w_max(row_length,rows)                 & ! col max vertical velocity (m/s)
 ,q(row_length,rows,model_levels)        & ! water vapour (kg/kg)
 ,theta(tdims%i_start:tdims%i_end,              & ! Theta (Kelvin)
        tdims%j_start:tdims%j_end,                                             &
                    1:tdims%k_end)

real(kind=r_bl), intent(in) ::                                                 &
  frac_land(nunstable)       & ! fraction of land in gridbox
 ,tv1_sd(nunstable)          & ! Approx to standard dev of level of
                               ! 1 virtual temperature (K).
 ,bl_vscale2(nunstable)        ! Velocity scale squared for boundary
                               ! layer eddies (m/s)

! Surface precipitation based 3d convective prognostic in kg/m2/s
real(kind=r_bl), intent(in) ::                                                 &
                    conv_prog_precip( tdims_s%i_start:tdims_s%i_end,           &
                                      tdims_s%j_start:tdims_s%j_end,           &
                                      tdims_s%k_start:tdims_s%k_end )

! Additional variables for SCM diagnostics which are dummy in full UM
integer, intent(in) ::                                                         &
  nSCMDpkgs          ! No of SCM diagnostics packages

logical, intent(in) ::                                                         &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages


integer, intent(in out) ::                                                     &
  ntml(row_length,rows)          & ! Number of model levels in the
                                   ! turbulently mixed layer.
 ,ntpar(row_length,rows)         & ! Max levels for parcel ascent
 ,nlcl(row_length,rows)            ! No. of model layers below the
                                   ! lifting condensation level.

logical, intent(in out) ::                                                     &
  cumulus(row_length,rows)       & ! Logical indicator for convection
 ,l_shallow(row_length,rows)     & ! Logical indicator for shallow Cu
 ,l_congestus(row_length,rows)   & ! Logical indicator for congestus Cu
 ,l_congestus2(row_length,rows)    ! Logical ind 2 for congestus Cu

! Convective type array ::
integer, intent(out) ::                                                        &
  conv_type(row_length, rows)
                                   ! Integer index describing convective type:
                                   !    0=no convection
                                   !    1=non-precipitating shallow
                                   !    2=drizzling shallow
                                   !    3=warm congestus
                                   !    ...
                                   !    8=deep convection

real(kind=r_bl), intent(in out) ::                                             &
  zh(row_length,rows)              ! Height above surface of top of boundary
                                   ! layer (metres).

real(kind=r_bl), intent(out) ::                                                &
  zhpar(row_length,rows)         & ! Height of max parcel ascent (m)
 ,z_lcl(row_length,rows)         & ! Height of lifting condensation
                                   !  level (not a model level) (m)
 ,z_lcl_uv(row_length,rows)      & ! Height of lifting condensation
                                   ! level on nearest uv level (m)
 ,dzh(row_length,rows)           & ! Inversion thickness (m)
 ,qcl_inv_top(row_length,rows)   & ! Parcel water content at inversion top
 ,delthvu(row_length,rows)       & ! Integral of undilute parcel buoyancy
                                   ! over convective cloud layer
                                   ! (for convection scheme)
 ,ql_ad(row_length,rows)         & ! adiabatic liquid water content at
                                   ! inversion or cloud top (kg/kg)
 ,cape(row_length, rows)         & ! CAPE from undilute parcel ascent (m2/s2)
 ,cin(row_length, rows)          & ! CIN from undilute parcel ascent (m2/s2)
 ,entrain_coef(row_length,rows)  & ! Entrainment coefficient
 ,qsat_lcl(row_length,rows)        ! qsat at cloud base


!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
integer ::                                                                     &
  i,j        & ! LOCAL Loop counter (horizontal field index).
 ,ii         & ! Local compressed array counter.
 ,k          & ! LOCAL Loop counter (vertical level index).
 ,ntpar_l    & ! cloud top level
 ,mbl          ! Maximum number of model layers allowed in the
               ! mixing layer; set to bl_levels-1.

integer ::                                                                     &
  icode        ! error code

! Arrays holding various key model levels

integer ::                                                                     &
  kshmin(nunstable)        & ! Position of buoyancy minimum above
                             ! topbl (used for shallow Cu diag)
 ,k_plume(nunstable)       & ! start level for surface-driven plume
 ,k_max(nunstable)         & ! level of max parcel buoyancy
 ,k_max_dil(nunstable)     & ! level of max parcel buoyancy - dilute
 ,nlcl_min(nunstable)      & ! minimum allowed level of LCL
 ,k_neutral(nunstable)     & ! level of neutral parcel buoyancy
 ,k_neutral_dil(nunstable) & ! level of neutral parcel buoyancy - dilute
 ,k_inv(nunstable)         & ! level from inversion testing
 ,freeze_lev(nunstable)      ! freezing level


! Compressed arrays store only values for unstable columns
! names as original array plus _c

integer ::                                                                     &
  ntml_c(nunstable)                                                            &
 ,nlcl_c(nunstable)                                                            &
 ,nlfc_c(nunstable)                                                            &
 ,ntml_save(nunstable)                                                         &
 ,ntinv_c(nunstable)

real(kind=r_bl) ::                                                             &
  q_c(nunstable, model_levels)                                                 &
 ,qcl_c(nunstable, model_levels)                                               &
 ,qcf_c(nunstable, model_levels)                                               &
 ,z_full_c(nunstable, model_levels)                                            &
 ,z_half_c(nunstable, model_levels)                                            &
 ,exner_theta_levels_c(nunstable, model_levels)                                &
 ,exner_rho_c(nunstable, model_levels)                                         &
 ,P_theta_lev_c(nunstable, model_levels)                                       &
 ,P_c(nunstable, model_levels)                                                 &
 ,cloud_fraction_c(nunstable,model_levels)                                     &
 ,zh_c(nunstable)                                                              &
 ,zh_save(nunstable)                                                           &
 ,zh_itop_c(nunstable)                                                         &
 ,delthvu_c(nunstable)                                                         &
 ,cape_c(nunstable)                                                            &
 ,cin_c(nunstable)                                                             &
 ,ql_ad_c(nunstable)                                                           &
 ,pstar_c(nunstable)                                                           &
 ,qsat_lcl_c(nunstable)

real(kind=r_bl) ::                                                             &
  z_lcl_c(nunstable)     & ! LCL height rounded to nearest model level
 ,zlcl_c(nunstable)        ! Different variable exact height of LCL

real(kind=r_bl) :: column_rh(nunstable)
                           ! Column integrated relative humidity
real(kind=r_bl) :: column_rh_bl(nunstable)
                           ! Column integrated relative humidity in BL
real(kind=r_bl) :: column_q(nunstable)           ! Column integrated q
real(kind=r_bl) :: qsat_c(nunstable, model_levels) ! Qsat

logical ::                                                                     &
  l_dilute                 ! if true also do a dilute ascent


! Arrays only used for unstable calculations

logical ::                                                                     &
  shmin(nunstable)        & ! Flag for finding min in parcel buoyancy below
                            ! 3km (for shallow Cu)
 ,cumulus_c(nunstable)    & ! compressed cumulus
 ,shallow_c(nunstable)      ! compressed shallow

real(kind=r_bl) ::                                                             &
  t(nunstable, model_levels)            & ! temperature (from theta)
 ,tl(nunstable, model_levels)           & ! Ice/liquid water temperature,
                                          ! but replaced by T in LS_CLD.
 ,qw(nunstable, model_levels)           & ! Total water content
 ,svl(nunstable, model_levels)          & ! Liquid/frozen water virtual
                                          ! static energy over cp.
 ,env_svl(nunstable,model_levels)       & ! Density (virtual) static energy
                                          ! over cp for layer.
 ,par_svl(nunstable,model_levels)       & ! Density (virtual) static energy
                                          ! over cp of parcel for level.
 ,par_svl_dil(nunstable,model_levels)     ! Density (virtual) static energy
                                          ! over cp of parcel for level.

real(kind=r_bl) ::                                                             &
  t_lcl(nunstable)             & ! Temperature at lifting condensation level.
 ,p_lcl(nunstable)             & ! Pressure at lifting condensation level.
 ,sl_plume(nunstable)          & ! Liquid/frozen water static energy
                                 ! over cp for a plume rising without
                                 ! dilution from level 1.
 ,qw_plume(nunstable)          & ! qw for a plume rising without
                                 ! dilution from level 1.
 ,Dt_dens_parc_T(nunstable)    & ! t_dens_parc-t_dens_env at ntpar
 ,Dt_dens_parc_Tmin(nunstable) & ! t_dens_parc-t_dens_env at kshmin
 ,thv_pert(nunstable)          & ! threshold thv of parcel
 ,delthvu_add(nunstable)       & ! adjustment applied for inversion top
! Added for improved parcel top - mainly used for finding an ascent
! capped by an inversion.
 ,Dt_dens_parc_T2(nunstable)   & ! 2nd copy of Dt_dens_parc_T
 ,delthvu2(nunstable)          & !  2nd copy
 ,zh2(nunstable)               & !  2nd copy
 ,max_buoy(nunstable)          & ! max parcel buoyancy
 ,max_buoy_dil(nunstable)      & ! max parcel buoyancy dilute parcel
 ,ql_ad2(nunstable)            & ! ql_ad 2nd copy
 ,pot_en(nunstable)              ! parcel potential energy when overshooting

! parcel calculation

real(kind=r_bl) ::                                                             &
  T_parc(nunstable, model_levels)          & ! Temperature of parcel.
 ,ql_parc_c(nunstable,model_levels)        & ! parcel water
 ,t_dens_env(nunstable, model_levels)      & ! Density potential temperature
                                             ! of environment.
 ,denv_bydz(nunstable, model_levels)       & ! Gradient of density potential
                                             ! temp in the environment.
 ,dpar_bydz(nunstable, model_levels)       & ! Gradient of density potential
                                             ! temperature of the parcel.
 ,buoyancy(nunstable, model_levels)        & ! undilute parcel buoyancy (K)
 ,buoy_use(nunstable, model_levels)        & ! parcel buoyancy actually used (K)
 ,dqsatdz(nunstable, model_levels)           ! dqsat/dz along adiabat

! Arrays added for dilute parcel calculation

real(kind=r_bl) ::                                                             &
  entrain_fraction(nunstable,model_levels) & ! fraction of environmental
                                             ! air to mix with parcel
 ,t_parc_dil(nunstable,model_levels)       & ! dilute parcel temeperature
 ,buoyancy_dil(nunstable,model_levels)       ! dilute parcel buoyancy (K)

real(kind=r_bl) ::                                                             &
  z_surf          &  ! approx height of top of surface layer
 ,z_neut          &  ! interpolated height of neutral buoyancy
 ,dpe             &  ! change in parcel PE between levels
 ,buoy_ave        &  ! average buoyancy between levels
 ,dz              &  ! Layer depth
 ,a_plume_use     &  ! value of a_plume actually used
 ,ent_rate        &  ! entrainment rate
 ,r_over_a           ! radius of a level over radius of earth

! required for average w calculation

real(kind=r_bl) ::                                                             &
  dmass_theta(nunstable,model_levels) & ! r**2rho*dr on theta levels
 ,w_avg(nunstable)                    & ! mean w over layer (m/s)
 ,w_avg2(nunstable)                     ! mean w over layer (m/s)

real(kind=r_bl) :: zcld_thresh   ! threshold depth for diagnosis

real(kind=r_bl) ::                                                             &
  pstar_w_cape_limit(nunstable)   ! scaled critical vertical velocity

character (len=*), parameter ::  RoutineName = 'CONV_DIAG_COMP_6A'
character (len=errormessagelength) :: cmessage        ! error message

!-----------------------------------------------------------------------
! Required for SCM diagnostics
!-----------------------------------------------------------------------

!
! Model constants:
!

real(kind=r_bl) :: zcld                  ! Depth of cloud layer
logical :: Ttest              ! Test on temperature at top of cloud layer

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Mixing ratio, r,  versus specific humidity, q
!
! In most cases the expression to first order are the same
!
!  tl = T - (lc/cp)qcl - [(lc+lf)/cp]qcf
!  tl = T - (lc/cp)rcl - [(lc+lf)/cp]rcf  - equally correct definition
!
! thetav = theta(1+cvq)         accurate
!        = theta(1+r/repsilon)/(1+r) ~ theta(1+cvr) approximate
!
! svl = (tl+gz/cp)*(1+(1/repsilon-1)qt)
!     ~ (tl+gz/cp)*(1+(1/repsilon-1)rt)
!
! dqsat/dT = repsilon*Lc*qsat/(R*T*T)
! drsat/dT = repsilon*Lc*rsat/(R*T*T)  equally approximate
!
! Only altering the expression for vapour pressure
!
!  e = qp/repsilon       - approximation
!  e = rp/(repsilon+r)   - accurate
!
!-----------------------------------------------------------------------
! 1.0 set variables
!-----------------------------------------------------------------------
!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

mbl = bl_levels - 1

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(nunstable,cumulus_c,shallow_c,freeze_lev,ntml_c,                  &
!$OMP        nlcl_c,nlfc_c,model_levels,index_i,index_j,r_theta_levels,        &
!$OMP        planet_radius,r_rho_levels,dmass_theta,rho_theta,zh_c,zh,         &
!$OMP        z_lcl_c,z_half,pstar_c,pstar,t,theta,exner_theta_levels,          &
!$OMP        t_parc,q_c,q,qcl_c,qcl,qcf_c,qcf,z_full_c,z_full,z_half_c,        &
!$OMP        exner_theta_levels_c,exner_rho_c,exner_rho,p_c,p,                 &
!$OMP        p_theta_lev_c,p_theta_lev,cloud_fraction_c,cloud_fraction,        &
!$OMP        qw,tl,lcrcp,lsrcp,svl,gamma_dry,c_virtual,t_dens_env,             &
!$OMP        k_plume,forced_cu,limit_pert_opt,sl_plume,tv1_sd,thv_pert,        &
!$OMP        qw_plume,rows,row_length,ntml)                                    &
!$OMP private(ii,i,j,k,r_over_a,z_surf,a_plume_use)
!-----------------------------------------------------------------------
! 1.1 initialise just unstable points
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do ii=1,nunstable

  cumulus_c(ii)    = .false.
  shallow_c(ii)    = .false.

  freeze_lev(ii) = 1

  ntml_c(ii) = 1
  nlcl_c(ii) = 1
  nlfc_c(ii) = 0  ! implies unset

end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 2.0 Calculate mass of layers for use later
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do k=1, model_levels-1
  do ii=1, nunstable
    i = index_i(ii)
    j = index_j(ii)
    r_over_a = r_theta_levels(i,j,k)/planet_radius
    dmass_theta(ii,k) = rho_theta(i,j,k)*r_over_a*r_over_a                     &
                           *(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 3.0 Calculate various quantities required by the parcel calculations
!-----------------------------------------------------------------------
! compress boundary layer depth
!$OMP do SCHEDULE(STATIC)
do ii=1,nunstable
  i = index_i(ii)
  j = index_j(ii)
  zh_c(ii)  = zh(i,j)
  z_lcl_c(ii) = z_half(i,j,nlcl_c(ii)+1)
  pstar_c(ii) = pstar(i,j)
end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do k = 1, model_levels
  do ii=1, nunstable
    i = index_i(ii)
    j = index_j(ii)
    t(ii,k) = theta(i,j,k) * exner_theta_levels(i,j,k)

    ! initialise t_parc at all points
    ! added for safety of qsat_mix calls later

    t_parc(ii,k) = t(ii,k)
    q_c(ii,k)   = q(i,j,k)
    qcl_c(ii,k) = qcl(i,j,k)
    qcf_c(ii,k) = qcf(i,j,k)
    z_full_c(ii,k) = z_full(i,j,k)
    z_half_c(ii,k) = z_half(i,j,k)
    exner_theta_levels_c(ii,k) = exner_theta_levels(i,j,k)
    exner_rho_c(ii,k)          = exner_rho(i,j,k)
    P_c(ii,k)              = p(i,j,k)
    P_theta_lev_c(ii,k)    = P_theta_lev(i,j,k)
    cloud_fraction_c(ii,k) = cloud_fraction(i,j,k)
  end do
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! 3.1 Calculate total water content, qw and Liquid water temperature, tl
!     Definitions for qw and tl the same whether q, qcl, qcf  are
!     specific humidities  or mixing ratio.
!
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do k=1,model_levels
  do ii=1, nunstable

    ! Total water   - BL doc P243.10

    qw(ii,k) = q_c(ii,k) + qcl_c(ii,k) + qcf_c(ii,k)

    ! Liquid water temperature as defined  BL doc P243.9

    tl(ii,k) = t(ii,k) - lcrcp*qcl_c(ii,k) - lsrcp*qcf_c(ii,k)


    ! Calculate svl: conserved variable  a form of moist static energy /cp
    !       svl = (tl+gz/cp)*(1+(1/repsilon-1)qt)  - specific humidity

    svl(ii,k) = ( tl(ii,k) + gamma_dry * z_full_c(ii,k) )                      &
                                     * ( one + c_virtual*qw(ii,k) )

    ! Density potential temperature of environment (K)

    t_dens_env(ii,k)=t(ii,k)*(one+                                             &
         c_virtual*q_c(ii,k)-qcl_c(ii,k)-qcf_c(ii,k))

  end do
end do
!$OMP end do
!-----------------------------------------------------------------------
! 4.0 Parts of parcel calculation common to all options
!-----------------------------------------------------------------------
! 4.1 Work out initial parcel properties and LCL
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do ii=1, nunstable

  k_plume(ii) = 1
  !-----------------------------------------------------------------------
  ! Only perform parcel ascent if unstable
  ! Start plume ascent from grid-level above top of surface layer, taken
  ! to be at a height, z_surf, given by 0.1*zh
  !-----------------------------------------------------------------------
  z_surf = 0.1_r_bl * zh_c(ii)

  do while( z_full_c(ii,k_plume(ii))  <   z_surf .and.                         &
    !                   ! not reached z_surf
                svl(ii,k_plume(ii)+1)  <   svl(ii,k_plume(ii)) )
    !                   ! not reached inversion

    k_plume(ii) = k_plume(ii) + 1

  end do
end do          ! loop over ii
!$OMP end do NOWAIT

a_plume_use = a_plume
if (forced_cu >= on) a_plume_use = zero  ! no need to add bonus buoyancy

if (limit_pert_opt == 2) then
!$OMP do SCHEDULE(STATIC)
  do ii=1, nunstable
    sl_plume(ii) = tl(ii,k_plume(ii)) + gamma_dry * z_full_c(ii,k_plume(ii))
    thv_pert(ii) = min( max( a_plume_use,                                      &
                       min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) ),       &
                       max_diag_thpert )
    qw_plume(ii) = qw(ii,k_plume(ii))
  end do
!$OMP end do NOWAIT
else if (limit_pert_opt == 0 .or. limit_pert_opt == 1) then
!$OMP do SCHEDULE(STATIC)
  do ii=1, nunstable
    sl_plume(ii) = tl(ii,k_plume(ii)) + gamma_dry * z_full_c(ii,k_plume(ii))
    thv_pert(ii) = max( a_plume_use,                                           &
                       min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
    qw_plume(ii) = qw(ii,k_plume(ii))
  end do
!$OMP end do NOWAIT
end if

!-----------------------------------------------------------------------
! Reset zh  (at this point in the code ntml is initialised as =1)
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do j=1,rows
  do i=1,row_length
    zh(i,j) = z_half(i,j,ntml(i,j)+1)
  end do
end do
!$OMP end do NOWAIT

!$OMP end PARALLEL
!-----------------------------------------------------------------------
! 4.2 Calculate temperature and pressure of lifting condensation level
!       using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
call lift_cond_lev ( nunstable, model_levels, k_plume,                         &
                     pstar_c, q_c, t,                                          &
                     P_theta_lev_c, exner_rho_c, z_half_c,                     &
                     t_lcl, p_lcl, zlcl_c, qsat_lcl_c )

!-----------------------------------------------------------------------
! Convection scheme requires NLCL to be at least 2.
! For model stability over mountains, also require ZLCL > 150m
! (approx nlcl=2 for G3 levels)
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(nunstable,z_half_c,mbl,nlcl_min,                                  &
!$OMP        index_i,index_j,zh_c,zh,zh2,z_lcl,zlcl_c,qsat_lcl,                &
!$OMP        qsat_lcl_c,model_levels,p_lcl,p_c,nlcl_c,z_lcl_c,z_full_c,        &
!$OMP        t,freeze_lev,nlcl,z_lcl_uv)                                       &
!$OMP private(ii,k,i,j)
!$OMP do SCHEDULE(STATIC)
do ii=1, nunstable

  k=3
  do while ( z_half_c(ii,k) < 150.0_r_bl .and. k < mbl )
    k=k+1
  end do
  nlcl_min(ii) = k-1

end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)
  zh_c(ii) = zh(i,j)
  zh2(ii)  = zh_c(ii)
  z_lcl(i,j)    = zlcl_c(ii)         ! expand up accurate z_lcl
  qsat_lcl(i,j) = qsat_lcl_c(ii)     ! expand up qsat at cloud base
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p,T,q       nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv,p,rho    nlcl+1,  z_lcl , p(nlcl+1)     either
!     + + + + + + + +   lcl, Plcl, not a model level            lower part
!    ---------------   p,T,q       nlcl , p_theta(nlcl+1)          of layer
!
!    - - - - - - - -   uv,p,rho    nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p,T,q       nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv,p,rho    nlcl+1,  z_lcl , p(nlcl+1)        or
!                                                                upper part
!    ---------------   p,T,q       nlcl , p_theta_lev(nlcl+1)      of layer
!
!    - - - - - - - -   uv,p,rho    nlcl   p(nlcl)
!
!-----------------------------------------------------------------------

do k = 2,model_levels
!$OMP do SCHEDULE(STATIC)
  do ii=1, nunstable

    ! NLCL level
    if ( p_lcl(ii)  <   P_c(ii,k) ) then

      ! compressed copies
      nlcl_c(ii)  = k-1
      z_lcl_c(ii) = z_half_c(ii,nlcl_c(ii)+1)
    end if     ! test on p_lcl

    ! Freezing level

    if (t(ii,k) <  tm .and. t(ii,k-1) >= tm) then
      if (freeze_lev(ii) == 1) then
        freeze_lev(ii) = k
      end if
    end if

  end do       ! ii loop
!$OMP end do NOWAIT
end do         ! k loop

! expand to full arrays
!$OMP do SCHEDULE(STATIC)
do ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)
  nlcl(i,j) = nlcl_c(ii)
  z_lcl_uv(i,j)   = z_full_c(ii,nlcl_c(ii))
end do       ! ii loop
!$OMP end do NOWAIT
!$OMP end PARALLEL

!=======================================================================
! 5.0 Parcel calculation - various options available
!                          (only on unstable points)
!=======================================================================
!
!  icvdiag option   |      explaination
!  -------------------------------------------------------------------
!     1             |  undilute parcel calculation
!  -  > 2 - - - -   |  - dilute parcel calculation - - - - - - - - -
!     2             | 0.55/z entrainment rates
!     3             |  1/z entrainment rate
!  -------------------------------------------------------------------
!     >3        Experimental options
!     4             |  Diurnal cycle over land entrainment rates,
!                   !  Ocean used 0.55/Z dilution of parcel
!     5             |  Diurnal cycle over land entrainment rates,
!                   !  Ocean undilute ascent of parcel
!=======================================================================
! Original 5A code - undilute parcel
!=======================================================================
if (icvdiag == 1) then
  !=======================================================================

  l_dilute = .false.

  !=======================================================================
  ! Dilute parcel ascent options
  !=======================================================================

else if (icvdiag >= 2 .and. icvdiag < 15 ) then

  l_dilute = .true.       ! dilute ascent required as well as undilute

  ! Set top wet level entrainment to zero
  k=model_levels
  do ii=1, nunstable
    entrain_fraction(ii,k) = zero
  end do

  !-----------------------------------------------------------------------
  !  Option 2 - 0.55/z entrainment rate as Jakob & Siebesma
  !             Entrainment rate not passed to convection scheme.
  !-----------------------------------------------------------------------

  if (icvdiag == 2) then

    do k= 1,model_levels-1
      do ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)
        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = 0.55_r_bl * dz/z_full_c(ii,k)
      end do
    end do

    !-----------------------------------------------------------------------
    !  Option 3 - 1/z entrainment rate, similar to many CRM results for
    !             deep and shallow convection.
    !             Entrainment rate not passed to convection scheme.
    !-----------------------------------------------------------------------
  else if (icvdiag == 3) then

    do k= 1,model_levels-1
      do ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)
        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = one * dz/z_full_c(ii,k)
      end do
    end do

    !-----------------------------------------------------------------------
    !  Option 4 - 0.55/z entrainment rate  for ocean
    !             Diurnal varying rate over land.
    !             Entrainment rate for land passed to convection scheme.
    !-----------------------------------------------------------------------

  else if (icvdiag == 4 ) then

    !--------------------------------------------------------------------
    ! Diurnal cycle changes
    !--------------------------------------------------------------------
    ! Only applied to land points (points with > 80% land if coastal).
    ! x/z entrainment rates - where x depends on BL depth (zlcl)
    ! If zlcl > 1000m then assumes convection uses so called equilbrium
    ! entrainment rates for land, i.e. 0.55/z ?
    ! No constraint on maximum entrainment rate (formula => 10.23,
    ! if zlcl=0.0)
    !--------------------------------------------------------------------

    do k= 1,model_levels-1
      do ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)

        if (frac_land(ii) > 0.8_r_bl) then
          ent_rate = 0.55_r_bl +                                               &
                        8.0_r_bl*(1.2_r_bl-(z_lcl(i,j)/1000.0_r_bl))*          &
                        (1.2_r_bl-(z_lcl(i,j)/1000.0_r_bl))
          entrain_coef(i,j) = ent_rate

          ! While very high entrainment rates are ok in the diagnosis
          ! the shallow convection scheme does not perform well
          ! if passed an entrainment rate ~10. Therefore restricting
          ! entrainment rates passed to convection.
          if (entrain_coef(i,j) > 3.5) then
            entrain_coef(i,j) = 3.5
          end if
          if (z_lcl(i,j) > 1200.0_r_bl) then
            ent_rate = 0.55_r_bl
            entrain_coef(i,j) = 0.55_r_bl
          end if

          !--------------------------------------------------------------------
          ! Suppressing convection over land where there is steep orography
          ! and very high vertical velocities may tend to cause gridpoint
          ! storms. A test has been added to revert to the standard sea
          ! entrainment rates where the w based CAPE closure is being used to
          ! control gridpoint storms.
          !--------------------------------------------------------------------

          if (cldbase_opt_dp == 3 .or. cldbase_opt_dp == 4 .or.                &
              cldbase_opt_dp == 5 .or. cldbase_opt_dp == 6 .or.                &
              cldbase_opt_dp == 7 .or. cldbase_opt_dp == 8 ) then
            ! Pressure weighted w_max test
            pstar_w_cape_limit(ii) = w_cape_limit * pstar_c(ii) / 1.0e5_r_bl
            if (w_max(i,j) > pstar_w_cape_limit(ii)) then
              ent_rate = 0.55_r_bl
              entrain_coef(i,j) = -99.0_r_bl
            end if
          end if

        else               ! Sea points  or coastal points
          !------------------------------------------------------------------
          ! Sea points - using a dilute 0.55/z diagnosis but standard
          !              entrainment rates for shallow and deep convection
          !------------------------------------------------------------------
          ent_rate = 0.55_r_bl
          entrain_coef(i,j) = -99.0_r_bl

        end if        ! land fraction test

        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = ent_rate * dz/z_full_c(ii,k)

      end do          ! unstable point loop
    end do            ! level loop


    !-----------------------------------------------------------------------
    !  Option 5 - undilute ocean diagnosis
    !             Diurnal varying dilute ascent over land.
    !             Entrainment rate for land passed to convection scheme.
    !-----------------------------------------------------------------------

  else if (icvdiag == 5 ) then

    !--------------------------------------------------------------------
    ! Diurnal cycle changes
    !--------------------------------------------------------------------
    ! Only applied to land points (points with > 80% land if coastal).
    ! x/z entrainment rates - where x depends on BL depth (zlcl)
    ! If zlcl > 1000m then assumes convection uses so called equilbrium
    ! entrainment rates for land, i.e. 0.55/z ?
    ! No constraint on maximum entrainment rate (formula => 10.23,
    ! if zlcl=0.0)
    !--------------------------------------------------------------------

    do k= 1,model_levels-1
      do ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)

        if (frac_land(ii) > 0.8_r_bl) then

          ent_rate = 0.55_r_bl +                                               &
                        8.0_r_bl*(1.2_r_bl-(z_lcl(i,j)/1000.0_r_bl))*          &
                        (1.2_r_bl-(z_lcl(i,j)/1000.0_r_bl))
          entrain_coef(i,j) = ent_rate

          ! While very high entrainment rates are ok in the diagnosis
          ! the convection scheme does not perform well
          ! if passed an entrainment rate ~10. Therefore restricting
          ! entrainment rates passed to convection.
          if (entrain_coef(i,j) > 3.5_r_bl) then
            entrain_coef(i,j) = 3.5_r_bl
          end if
          if (z_lcl(i,j) > 1200.0_r_bl) then
            ent_rate = 0.55_r_bl
            entrain_coef(i,j) = 0.55_r_bl   ! uses 0.55/z
          end if

          !--------------------------------------------------------------------
          ! Suppressing convection over land where there is steep orography
          ! and very higher vertical velocities may tend to cause gridpoint
          ! storms. A test has been added to revert to the standard sea
          ! entrainment rates where the w based CAPE closure is being used to
          ! control gridpoint storms.
          !--------------------------------------------------------------------

          if (cldbase_opt_dp == 3 .or. cldbase_opt_dp == 4 .or.                &
              cldbase_opt_dp == 5 .or. cldbase_opt_dp == 6 .or.                &
              cldbase_opt_dp == 7 .or. cldbase_opt_dp == 8 ) then
             ! Pressure weighted w_max test
            pstar_w_cape_limit(ii) = w_cape_limit * pstar_c(ii) / 1.0e5_r_bl
            if (w_max(i,j) > pstar_w_cape_limit(ii)) then
              ent_rate = 0.55_r_bl
              entrain_coef(i,j) = -99.0_r_bl
            end if
          end if

        else               ! Sea points  or coastal points
          !------------------------------------------------------------------
          ! Sea points - no dilute ascent
          !------------------------------------------------------------------
          ent_rate = zero
          entrain_coef(i,j) = -99.0_r_bl

        end if        ! land fraction test

        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = ent_rate * dz/z_full_c(ii,k)

      end do          ! unstable point loop
    end do            ! level loop

    !-----------------------------------------------------------------------
    !  Option 6 - prognostic convection diagnosis with p/pstar^2 entrainment
    ! based on the prognostic at the initiation level
    !-----------------------------------------------------------------------

  else if (icvdiag == 6) then

    do k= 1,model_levels-1
      do ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)
        if (conv_prog_precip(i,j,k) > 1.0e-10_r_bl) then
          ent_rate = prog_ent_grad *                                           &
                     log10(conv_prog_precip(i,j,nlcl_c(ii)) *                  &
                     refqsat/qsat_lcl_c(ii)) + prog_ent_int
          !Limit scaling
          ent_rate = min(max(ent_rate, prog_ent_min), prog_ent_max)
        else
          ent_rate = prog_ent_max
        end if

        ent_rate = entcoef * ae2 * ent_rate * ent_fac_dp *                     &
                   p_theta_lev(i,j,k) / (pstar(i,j)**2)  *                     &
                   rho_theta(i,j,k) * g

        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = ent_rate * dz
        entrain_coef(i,j) = -99.0_r_bl
      end do          ! unstable point loop
    end do            ! level loop

    !-----------------------------------------------------------------------
    !  Option 7 - prognostic convection diagnosis with p/pstar^2 entrainment
    ! based on the prognostic at the current level
    !-----------------------------------------------------------------------

  else if (icvdiag == 7) then

    do k= 1,model_levels-1
      do ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)
        if (conv_prog_precip(i,j,k) > 1.0e-10_r_bl) then
          ent_rate = prog_ent_grad *                                           &
                     log10(conv_prog_precip(i,j,k) *                           &
                     refqsat/qsat_lcl_c(ii)) + prog_ent_int
          !Limit scaling
          ent_rate = min(max(ent_rate, prog_ent_min), prog_ent_max)
        else
          ent_rate = prog_ent_max
        end if

        ent_rate = entcoef * ae2 * ent_rate * ent_fac_dp *                     &
                   p_theta_lev(i,j,k) / (pstar(i,j)**2)  *                     &
                   rho_theta(i,j,k) * g

        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = ent_rate * dz
        entrain_coef(i,j) = -99.0_r_bl
      end do          ! unstable point loop
    end do            ! level loop

    !-----------------------------------------------------------------------
    !  Option 8 - convection diagnosis with p/pstar^2 entrainment as in the
    !  main deep ascent under ent_opt_dp=0
    !-----------------------------------------------------------------------

  else if (icvdiag == 8) then

    do k= 1,model_levels-1
      do ii=1, nunstable
        i = index_i(ii)
        j = index_j(ii)

        ent_rate = entcoef * ae2 * ent_fac_dp *                                &
                   p_theta_lev(i,j,k) / (pstar(i,j)**2)  *                     &
                   rho_theta(i,j,k) * g

        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = ent_rate * dz
        entrain_coef(i,j) = -99.0_r_bl
      end do          ! unstable point loop
    end do            ! level loop

  else if (icvdiag == 10 ) then ! Ben's diagnosis

    call calc_column_rh(nunstable, nlcl_c                                      &
       , row_length,rows                                                       &
       , index_i, index_j                                                      &
       , p_c, t, q_c, rho_only                                                 &
       , z_half_c, 20000.0_r_bl                                                &
       , qsat_c, column_rh, column_rh_bl                                       &
       , column_q )

    do ii=1, nunstable
      i=index_i(ii)
      j=index_j(ii)
      if ( column_rh(ii) >0.8_r_bl) then
        ent_rate = zero
        entrain_coef(i,j) = -99.0_r_bl
      else if ( column_rh(ii) >.7 ) then
        ent_rate = 0.55_r_bl
        entrain_coef(i,j) = -99.0_r_bl
      else
        ent_rate = 1.55_r_bl
        entrain_coef(i,j) = -99.0_r_bl
      end if

      do k= 1,model_levels-1
        dz = r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k)
        entrain_fraction(ii,k) = ent_rate * dz/z_full_c(ii,k)
      end do
    end do




  else

    icode = 3     ! fatal return code
    ! Statement assumes icvdiag is in the range that format statement covers.
    write(cmessage,'(a40,i6)')                                                 &
       'Convection diagnosis option not allowed ',icvdiag
    call ereport(routinename,icode,cmessage)
  end if    !   diagnosis option

  !=======================================================================
  ! Option ?. New diagnosis - yet to be written
  !=======================================================================

else

  icode = 2     ! fatal return code
  ! Statement assumes icvdiag is in the range that format statement covers.
  write(cmessage,'(a40,i6)')                                                   &
     'Convection diagnosis option not allowed ',icvdiag
  call ereport(routinename,icode,cmessage)

  !=======================================================================
  !    End of choice of diagnosis code
  !=======================================================================
end if

!-----------------------------------------------------------------------
! Parcel ascent
!-----------------------------------------------------------------------
! Parcel starts from level k_plume and is lifted up.
! If a dilute parcel ascent - mix in environmental air above lifting
! condensation level
!----------------------------------------------------------------------
call parcel_ascent ( nunstable, nSCMDpkgs,                                     &
                       nlcl_c, k_plume,                                        &
                       l_dilute, L_SCMDiags,                                   &
                       sl_plume, qw_plume,                                     &
                       t, tl, qcl_c, qcf_c, qw, t_dens_env,                    &
                       p_theta_lev_c, exner_theta_levels_c,                    &
                       z_full_c, entrain_fraction,                             &
                       T_parc, t_parc_dil, ql_parc_c,                          &
                       buoyancy,buoyancy_dil,                                  &
                       env_svl, par_svl, par_svl_dil,                          &
                       denv_bydz, dpar_bydz, dqsatdz )


!-----------------------------------------------------------------------
! Tests on parcel ascent - look for an inversion or not
!-----------------------------------------------------------------------

select case(cvdiag_inv)

case (0)    ! no inversion testing

  !-----------------------------------------------------------------------
  !   Now compare plume s_VL with each model layer s_VL in turn to
  !     find the first time that plume has negative buoyancy.
  !-----------------------------------------------------------------------
  call cv_parcel_neutral_dil(nunstable,                                        &
             nlcl_c,k_plume,l_dilute,                                          &
             z_lcl_c, thv_pert, z_full_c, z_half_c, exner_theta_levels_c,      &
             buoyancy, buoyancy_dil, t_dens_env, dqsatdz,                      &
             zh_c,                                                             &
             k_max,k_max_dil,k_neutral,k_neutral_dil,                          &
             max_buoy,max_buoy_dil,ql_ad_c,delthvu_c,                          &
             cape_c,cin_c)

  !-----------------------------------------------------------------------
  ! Default parcel top properties are assumed to be those when the
  ! ascent reaches the level of neutral buoyancy.
  !-----------------------------------------------------------------------

  if (l_dilute) then

    do ii=1,nunstable
      ntml_c(ii) = k_neutral_dil(ii)
      ntinv_c(ii)   = ntml_c(ii)
      zh_itop_c(ii) = -one   ! to test when set sensibly
    end do
    do k=1,model_levels
      do ii=1,nunstable
        buoy_use(ii,k)=buoyancy_dil(ii,k)+thv_pert(ii)
      end do
    end do

  else     ! undilute ascent

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(nunstable,ntml_c,k_neutral,ntinv_c,zh_itop_c,model_levels,        &
!$OMP        buoy_use,buoyancy,thv_pert)                                       &
!$OMP private(ii,k)
!$OMP do SCHEDULE(STATIC)
    do ii=1,nunstable
      ntml_c(ii) = k_neutral(ii)
      ntinv_c(ii)   = ntml_c(ii)
      zh_itop_c(ii) = -one   ! to test when set sensibly
    end do    ! ii loop
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
    do k=1,model_levels
      do ii=1,nunstable
        buoy_use(ii,k)=buoyancy(ii,k)+thv_pert(ii)
      end do ! ii loop
    end do
!$OMP end do
!$OMP end PARALLEL
  end if ! test on l_dilute
  !-----------------------------------------------------------------------
  ! Estimate inversion thickness
  !-----------------------------------------------------------------------
  if (forced_cu >= on .or. bl_res_inv /= off) then
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(model_levels,nunstable,nlfc_c,nlcl_c,buoy_use,delthvu_add,        &
!$OMP        z_full_c,pot_en,bl_vscale2,zh_itop_c,zh_c,t_dens_env,g,           &
!$OMP        z_half_c,ntinv_c,delthvu_c,ntml_c,exner_theta_levels_c,           &
!$OMP        ntml_save,zh_save,denv_bydz,dpar_bydz,cumulus_c)                  &
!$OMP private(k,ii,z_neut,dz,buoy_ave,dpe)
    ! find the level of free convection
    !  - nlfc is the first grid level above the LCL that is buoyant
    do k=2,model_levels-1
!$OMP do SCHEDULE(STATIC)
      do ii=1,nunstable
        if (nlfc_c(ii) == 0 .and. k > nlcl_c(ii) .and.                         &
            ! LFC unset but above LCL
            buoy_use(ii,k) > zero) then
            ! become positively buoyant
          nlfc_c(ii) = k
        end if
      end do
!$OMP end do NOWAIT
    end do

!$OMP do SCHEDULE(STATIC)
    do ii=1,nunstable
      delthvu_add(ii)=zero   ! initiailise
      k = ntml_c(ii)
      ! Find height where buoy_use=0, linearly interpolating
      if ( abs( buoy_use(ii,k) - buoy_use(ii,k+1) ) > zero ) then
        z_neut = z_full_c(ii,k) - buoy_use(ii,k) *                             &
                                (z_full_c(ii,k+1)-z_full_c(ii,k))/             &
                                     (buoy_use(ii,k+1)-buoy_use(ii,k))
      else
        ! Avoid divide-by-zero caused by occasional unexpected values
        z_neut = z_full_c(ii,k)
      end if
      dz = z_full_c(ii,k+1) - z_neut
      ! Note buoy_use=0 at z_neut, hence integrating to z_full(k+1) gives:
      pot_en(ii) = -one_half*buoy_use(ii,k+1)*dz*g/t_dens_env(ii,k+1)
      if (pot_en(ii) > bl_vscale2(ii) ) then
        zh_itop_c(ii) = max( zh_c(ii), z_neut -                                &
            2.0_r_bl*bl_vscale2(ii)*t_dens_env(ii,k+1)/(buoy_use(ii,k+1)*g) )
      end if
    end do
!$OMP end do NOWAIT

    do k=2,model_levels-1
!$OMP do SCHEDULE(STATIC)
      do ii=1,nunstable
        if (zh_itop_c(ii) < zero .and. k >= ntml_c(ii)+1) then
          ! look for inversion top
          dz = z_full_c(ii,k+1) - z_full_c(ii,k)
          buoy_ave = one_half*(buoy_use(ii,k)+buoy_use(ii,k+1))
          dpe = - buoy_ave * dz * g / t_dens_env(ii,k)
          if ( pot_en(ii)+dpe > bl_vscale2(ii) ) then
            zh_itop_c(ii) = max( zh_c(ii), z_full_c(ii,k) +                    &
                         (pot_en(ii)-bl_vscale2(ii))*t_dens_env(ii,k)/         &
                         ( buoy_ave*g ) )
            if (zh_itop_c(ii) >= z_half_c(ii,k+1)) then
              ntinv_c(ii) = k  ! marks top level within inversion
            else
              ntinv_c(ii) = k-1  ! marks top level within inversion
            end if
            if (delthvu_add(ii) > zero) then
              ! if parcel has overcome inhibition to become positively
              ! buoyant then adjust delthvu to allow convection triggering
              delthvu_c(ii)=delthvu_c(ii)+delthvu_add(ii)
            end if
          else
            pot_en(ii) = pot_en(ii) + dpe
            dz = z_half_c(ii,k+1) - z_half_c(ii,k)
            delthvu_add(ii) = delthvu_add(ii) + buoy_use(ii,k)*dz              &
                                          /exner_theta_levels_c(ii,k)
          end if
        end if
      end do  ! over ii
!$OMP end do NOWAIT
    end do  ! over k
!$OMP do SCHEDULE(STATIC)
    do ii=1,nunstable
      ntml_save(ii) = ntml_c(ii)
      zh_save(ii)   = zh_c(ii)
      if (ntinv_c(ii) > nlfc_c(ii)) then
        ! Inversion top is above LFC so implies free rather than forced cu
        !  - use inversion top as top of parcel ascent in cumulus diagnosis
        ntml_c(ii) = ntinv_c(ii)
        zh_c(ii)   = zh_itop_c(ii)
        ! Test whether cloud within inversion is capped by a significant
        ! inversion (lapse rate greater than moist adiabatic, given
        ! by dpar_bydz) - if not, trigger cumulus
        k = ntinv_c(ii)
        if ( ntinv_c(ii) > nlcl_c(ii) .and.                                    &
             denv_bydz(ii,k+1) < 1.25_r_bl*dpar_bydz(ii,k+1) ) then
          cumulus_c(ii) = .true.
        end if
      end if   ! ntinv above lfc
    end do
!$OMP end do NOWAIT
!$OMP end PARALLEL
  end if  ! test on forced_cu >= on or bl_res_inv /= off

  !-----------------------------------------------------------------------
  ! Average vertical velocity over a layer  - required for shallow
  !   convection test.
  !-----------------------------------------------------------------------
  ! Layer from top in cloud w to value 1500km above cloud top
  !-----------------------------------------------------------------------

  call mean_w_layer(nunstable,row_length,rows,model_levels,                    &
             ntml_c, index_i, index_j,                                         &
             1500.0_r_bl, z_full_c, z_half_c, w_copy, dmass_theta,             &
             w_avg)


case (1)    ! Original inversion test for shallow convection
           ! Only available for undilute ascent

  !-----------------------------------------------------------------------
  !   Now compare plume s_VL with each model layer s_VL in turn to
  !     find the first time that plume has negative buoyancy.
  !-----------------------------------------------------------------------

  call cv_parcel_neutral_inv(nunstable,                                        &
             nlcl_c,k_plume,                                                   &
             z_lcl_c, thv_pert, z_full_c, z_half_c, exner_theta_levels_c,      &
             buoyancy, t_dens_env, dqsatdz, denv_bydz, dpar_bydz,              &
             zh_c, zh2,                                                        &
             k_max,k_neutral,k_inv, kshmin, shmin,                             &
             max_buoy,dt_dens_parc_t,ql_ad_c,delthvu_c,cape_c,cin_c,           &
             dt_dens_parc_t2,dt_dens_parc_tmin,ql_ad2,                         &
             delthvu2)

  !-----------------------------------------------------------------------
  ! Average vertical velocity over a layer  - required for shallow
  !   convection test.
  !-----------------------------------------------------------------------
  ! Layer from top in cloud w to value 1500km above cloud top
  !-----------------------------------------------------------------------

  call mean_w_layer(nunstable,row_length,rows,model_levels,                    &
             k_neutral, index_i, index_j,                                      &
             1500.0_r_bl, z_full_c, z_half_c, w_copy, dmass_theta,             &
             w_avg)

  call mean_w_layer(nunstable,row_length,rows,model_levels,                    &
             k_inv, index_i, index_j,                                          &
             1500.0_r_bl, z_full_c, z_half_c, w_copy, dmass_theta,             &
             w_avg2)

  !-----------------------------------------------------------------------
  ! Default parcel top properties are assumed to be those when the
  ! ascent reaches the level of neutral buoyancy. These may not be those
  ! required in the case of shallow convection.
  ! Shallow convection requires the possible identification of an inversion
  ! at the top of the ascent. This may not be detected by the LNB test.
  ! The gradient tests are designed to detect the shallow top.
  !-----------------------------------------------------------------------
  ! Modify top if ascent is likely to be shallow

  do ii=1,nunstable

    if (shmin(ii) ) then    ! found an inversion
      ! points where k_inv not the same as k_neutral and level below freezing
      ! may be shallow or congestus or deep

      if (k_inv(ii) == k_neutral(ii)) then
        !  Both methods give same answer for top level leave shmin set
        ntml_c(ii) = k_neutral(ii)

        ! Inversion top lower than level of neutral buoyancy.
        ! Check also, either below freezing level or less than 2500m for
        ! shallow convection.

      else if ((k_inv(ii) <  freeze_lev(ii) .or.                               &
                         z_full_C(ii,k_inv(ii)+1)  <=  2500.0_r_bl )           &
                 .and. k_inv(ii) <  k_neutral(ii) ) then


        if ( (z_full_c(ii,kshmin(ii)) - z_full_c(ii,k_inv(ii)))                &
             <=  1.25_r_bl*(z_half_c(ii,k_inv(ii)+1) - z_lcl_c(ii)) .and.      &
             (dt_dens_parc_tmin(ii)  <=  0.55_r_bl*dt_dens_parc_t2(ii))        &
             .and.     (w_avg2(ii)  <   zero)  ) then

          ! May be shallow or congestus
          ! set values to those found from inversion testing
          ntml_c(ii)  = k_inv(ii)
          delthvu_c(ii) = delthvu2(ii)
          zh_c(ii)    = zh2(ii)
          w_avg(ii)   = w_avg2(ii)
          ql_ad_c(ii) = ql_ad2(ii)
          dt_dens_parc_t(ii) = dt_dens_parc_t2(ii)

        else   ! Assume not shallow or congestus

          ntml_c(ii) = k_neutral(ii)
          shmin(ii) = .false.  ! inversion top found not good
                                 ! don't do shallow tests
        end if

      else   ! Assume deep  and therefore top LNB
        ntml_c(ii) = k_neutral(ii)
        shmin(ii) = .false.  ! inversion top found not good
                                 ! don't do shallow tests
      end if        ! tests on k_inv

    else    !  No inversion found  i.e. shmin=false

      ntml_c(ii) = k_neutral(ii)

    end if     ! shmin test

  end do    ! ii loop

end select  ! test on type of testing inversion of not.

! Set convective type to none unless set otherwise later.
conv_type(:,:)=0

!-----------------------------------------------------------------------
! Test which points are cumulus
!-----------------------------------------------------------------------
call cumulus_test (nunstable,mbl,                                              &
                      ntml_c, k_plume,                                         &
                      zh_c, frac_land, qw, cloud_fraction_c, z_full_c,         &
                      z_lcl_c, nlcl_c, cumulus_c )

!==============================================================================
! Original shallow/deep diagnosis, no congestus diagnosed
! Or no parametrized convection scheme being called so expect iconv_congestus
! to be set to imdi
if (iconv_congestus == 0 .or. iconv_congestus == imdi) then


  !==============================================================================

  select case(cvdiag_inv)

  case (0)    ! no inversion testing
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED(nunstable,index_i,index_j,cumulus_c,conv_type,ntml_c,             &
!$OMP        z_full_c,t,l_dilute,iconv_deep,frac_land,shallow_c,               &
!$OMP        k_neutral,w_avg,cvdiag_sh_wtest)                                  &
!$OMP private(ii,i,j,ntpar_l)
    !CDIR NODEP
    do ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)
      if ( cumulus_c(ii) ) then
        ! Set convective type to deep unless set otherwise later.
        conv_type(i,j)=8

        !---------------------------------------------------------------------
        ! If cumulus has been diagnosed, determine whether it is shallow
        ! or deep convection
        !---------------------------------------------------------------------
        ! Conditions for shallow convection
        !  (1) top of parcel ascent < 2500. or T (top of parcel) > TM
        !  (2) Additional condition   w_avg < cvdiag_sh_wtest
        !---------------------------------------------------------------------
        ntpar_l = ntml_c(ii)

        if ( z_full_c(ii,ntpar_l)  <=  2500.0_r_bl .or.                        &
                 t(ii,ntpar_l)  >=  tm ) then

          if (l_dilute .or. iconv_deep == 2) then
            if (frac_land(ii) > 0.8_r_bl) then      ! land points no w test

              shallow_c(ii) = .true.
              conv_type(i,j)=1

              ! If undilute parcel would have been deep then reset to go
              ! through deep scheme but using the new entrain rate.

              if (t(ii,k_neutral(ii)) < tm) then
                shallow_c(ii) = .false.
                conv_type(i,j)=8
              end if

            else                       ! ocean points do w test

              if (w_avg(ii)  <  cvdiag_sh_wtest) then
                shallow_c(ii) = .true.
                conv_type(i,j)=1
              end if
            end if

          else   ! undilute ascent (icvdiag=1)

            if ( w_avg(ii)  <  cvdiag_sh_wtest ) then
              shallow_c(ii) = .true.
              conv_type(i,j)=1
            end if

          end if  ! l_dilute


        end if  !  height and temp test

      end if        ! test on cumulus

    end do          ! ii loop
!$OMP end PARALLEL do

  case (1)    ! Orignal inversion test
    !CDIR NODEP
    do ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)
      if ( cumulus_c(ii) ) then
        ! Set convective type to deep unless set otherwise later.
        conv_type(i,j)=8

        !---------------------------------------------------------------------
        ! If cumulus has been diagnosed, determine whether it is shallow
        ! or deep convection
        !---------------------------------------------------------------------
        ! Conditions for shallow convection
        ! (1) w_avg < cvdiag_sh_wtest     (descending air or weak ascent)
        ! (2) top of parcel ascent < 2500. or T (top of parcel) > TM
        ! (3) height of min buoyancy (above Bl) - height of parcel top T level
        !      <1.25(height parcel top - z lifting condensation level)
        ! (4) t_dens_parc -t_dens_env at kshmin <0.55t_dens_parc -t_dens_env
        !     at ntpar
        !
        ! The last 2 conditions are looking for a strong inversion at the top
        ! of the shallow cumulus.
        !---------------------------------------------------------------------
        if ( shmin(ii) ) then
          ntpar_l = ntml_c(ii)
          ! New deep turbulence scheme no w test over land for shallow
          if (iconv_deep == 2 .and. frac_land(ii) > 0.8_r_bl) then
            if ( (z_full_c(ii,ntpar_l)  <=  2500.0_r_bl .or.                   &
                     t(ii,ntpar_l)  >=  tm)                                    &
                .and. (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar_l))         &
                   <=  1.25_r_bl*(zh_c(ii) - z_lcl_c(ii)) .and.                &
                 Dt_dens_parc_Tmin(ii)  <=  0.55_r_bl*Dt_dens_parc_T(ii) ) then

              shallow_c(ii) = .true.
              conv_type(i,j)=1

              ! may be problem with ntpar diagnosis for deep if wadv test sets
              ! L_shallow  false
            end if

          else
            if ( w_avg(ii)  <  cvdiag_sh_wtest .and.                           &
             (z_full_c(ii,ntpar_l)  <=  2500.0_r_bl .or.                       &
                 t(ii,ntpar_l)  >=  tm)                                        &
            .and. (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar_l))             &
               <=  1.25_r_bl*(zh_c(ii) - z_lcl_c(ii)) .and.                    &
              Dt_dens_parc_Tmin(ii)  <=  0.55_r_bl*Dt_dens_parc_T(ii) ) then

              shallow_c(ii) = .true.
              conv_type(i,j)=1

              ! may be problem with ntpar diagnosis for deep if wadv test sets
              ! L_shallow  false

            end if
          end if  ! iconv_deep

        end if       ! test on shmin

      end if        ! test on cumulus

    end do          ! ii loop

  end select  ! test on type of testing inversion of not.

  !==============================================================================

else if (iconv_congestus == 1) then   ! Mass flux congestus option

  !==============================================================================

  select case(cvdiag_inv)

  case (0)    ! no inversion testing
    !CDIR NODEP
    do ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)

      if ( cumulus_c(ii) ) then
        ! Set convective type to deep unless set otherwise later.
        conv_type(i,j)=8
        !
                !-----------------------------------------------------------------------
                ! Conditions for congestus ;
                !    T(at top) >= freezing level ?
                !    height of top > 2.5km
                !    No need for air to be descending above ?
                !
                ! Conditions for shallow if congestus diagnosis are;
                !
                !    height of top < 2.5km
                !    Air to be descending above
                !
                ! Net result for tests mean shallow plus congestus points > shallow
                ! from old method.
                !-----------------------------------------------------------------------

        ntpar_l = ntml_c(ii)

        if ( (t(ii,ntpar_l)  >=  tm) .and.                                     &
             (z_full_c(ii,ntpar_l) >  2500.0_r_bl)) then

          l_congestus(i,j) = .true.

        else if (z_full_c(ii,ntpar_l) <= 2500.0_r_bl) then

          ! Top of convection below 2500km
          if (l_dilute) then   ! no w test over land
            if (frac_land(ii) > 0.8_r_bl) then      ! land points
              shallow_c(ii) = .true.

              if (t(ii,k_neutral(ii)) < tm) then
                shallow_c(ii) = .false.
              end if

            else         ! ocean/coastal points do w test

              if (w_avg(ii)  <  cvdiag_sh_wtest) then
                shallow_c(ii) = .true.
              else
                l_congestus(i,j) = .true.
              end if
            end if

          else    ! undilute ascents (icvdiag=1)

            ! Air above shallow convection less than a set value
            if (w_avg(ii)  <  cvdiag_sh_wtest) then
              shallow_c(ii) = .true.
            else   ! Air rising therefore class as congestus rather than deep
              l_congestus(i,j) = .true.
            end if    ! w_avg test

          end if      ! (l_dilute)
        end if        ! test on TM etc

        !  l_congestus - all congestus points
        !  l_congestus2 - points with descending air

        if (l_congestus(i,j)) then
          if (w_avg(ii)  <  cvdiag_sh_wtest ) then
            l_congestus2(i,j) = .true.
          end if     ! descending air
        end if       ! on congestus

      end if     ! test on cumulus

    end do      ! ii loop

  case (1)    ! original inversion test
    !CDIR NODEP
    do ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)

      if ( cumulus_c(ii) ) then
        ! Set convective type to deep unless set otherwise later.
        conv_type(i,j)=8

        !-----------------------------------------------------------------------
        ! Conditions for congestus ;
        !   As for shallow ie shmin = .true.
        !    T(at top) >= freezing level ?
        !    height of top > 2.5km  < 3.5km?
        !    No need for air to be descending above ?
        !
        ! Conditions for shallow if congestus diagnosis are;
        !
        !    height of top < 2.5km
        !    Air to be descending above
        !
        ! Net result for tests mean shallow plus congestus points > shallow
        ! from old method.
        !-----------------------------------------------------------------------
        if ( shmin(ii) ) then
          ntpar_l = ntml_c(ii)
          if ( (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar_l))                &
                    <=  1.25_r_bl*(zh_c(ii) - z_lcl_c(ii)) .and.               &
                 Dt_dens_parc_Tmin(ii)  <=  0.55_r_bl*Dt_dens_parc_T(ii) )     &
            then

            ! top of convection above 2500km and still below freezing level
            ! No check on descending air

            if ( (t(ii,ntpar_l)  >=  tm) .and.                                 &
               (z_full_c(ii,ntpar_l) >  2500.0_r_bl)) then
              l_congestus(i,j) = .true.
            else

              ! Top of convection below 2500km

              if (z_full_c(ii,ntpar_l) <= 2500.0_r_bl) then
                ! Air above shallow convection less than a set value
                if (w_avg(ii)  <  cvdiag_sh_wtest) then
                  shallow_c(ii) = .true.
                  conv_type(i,j)=2
                else   ! Air rising therefore class as congestus
                       ! rather than deep
                  if ((t(ii,ntpar_l)  <  tm)) then
                    conv_type(i,j)=4
                  else
                    conv_type(i,j)=3
                  end if
                  l_congestus(i,j) = .true.
                end if    ! w_avg test
              end if      ! test on 2500m

            end if        ! test on TM etc

          end if         ! shmin

          !  l_congestus - all congestus points
          !  l_congestus2 - points with descending air

          if (l_congestus(i,j)) then
            if (w_avg(ii)  <  cvdiag_sh_wtest ) then
              l_congestus2(i,j) = .true.
            end if     ! descending air
          end if       ! on congestus

        end if       ! test on shmin

      end if     ! test on cumulus

    end do      ! ii loop

  end select  ! test on type of testing inversion of not.

  do ii=1,nunstable
    i = index_i(ii)
    j = index_j(ii)
    select case (conv_type(i,j))
    case (1:4)
      shallow_c(ii) = .true. ! types 1-4 go through shallow logical
    case DEFAULT
      shallow_c(ii) = .false.
    end select
  end do      ! ii loop

  !=======================================================================

else if (iconv_congestus > 1) then  ! Warm TCS options

  !=======================================================================

  select case(iwtcs_diag1)
  case (1)
    zcld_thresh = 6000.0_r_bl
  case DEFAULT
    zcld_thresh = 4000.0_r_bl
  end select

  select case(cvdiag_inv)

  case (0)    ! no inversion testing

    ! Do we want this option ?

  case (1)    ! Standard Inversion testing

    !CDIR NODEP
    do ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)
      zcld=zh_c(ii) - z_lcl_c(ii)     ! depth of cloud layer
      Ttest = (t(ii,ntml_c(ii))  >=  tm) ! cloud top temperature > freezing
      if ( cumulus_c(ii) ) then
        ! Set convective type to deep unless set otherwise later.
        conv_type(i,j)=8
        !-----------------------------------------------------------------------
        ! Conditions for setting conv_type ;
        !       conv_type = 8:  depth of cloud layer > 4km
        !       conv_type = 4:  depth of cloud layer > 2km and < 4km
        !                       and T(at top) > freezing
        !       conv_type = 3:  depth of cloud layer > 2km and < 4km
        !                       and T(at top) > freezing
        !       conv_type = 2:  depth of cloud layer < 2km and > 500m
        !       conv_type = 1:  depth of cloud layer < 500m
        !-----------------------------------------------------------------------
        if ( shmin(ii) .and. zcld < zcld_thresh ) then

          ntpar_l = ntml_c(ii)
          if ( (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar_l))                &
                    <=  1.25_r_bl*zcld .and.                                   &
                 Dt_dens_parc_Tmin(ii)  <=  0.55_r_bl*Dt_dens_parc_T(ii) )     &
             then

            shallow_c(ii)=.true. ! types 1-4 go through shallow logical

            if ( zcld >  2000.0_r_bl) then
              ! depth of convection >  2km and < 4km
              if (Ttest) then
                conv_type(i,j)=3
              else
                conv_type(i,j)=4
              end if

            else if (zcld <=  2000.0_r_bl .and. zcld >  500.0_r_bl) then
              ! depth of convection < 2km and > 500m
              conv_type(i,j)=2
            else
              ! depth of convection < 500m
              conv_type(i,j)=1
            end if

          end if        ! test on TM etc

        end if       ! test on shmin

      end if     ! test on cumulus

    end do      ! ii loop

  end select

  !=======================================================================

end if            ! test on congestus diagnosis

!=======================================================================


!=======================================================================
! Manipulate diagnosis depending on switches
!=======================================================================
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(nunstable,index_i,index_j,iwtcs_diag2,conv_type,                  &
!$OMP        shallow_c,cumulus_c,nlcl_c,p_lcl,p_theta_lev_c,ntml_c,            &
!$OMP        mbl,zh_c,z_half_c,nlcl_min,l_congestus,delthvu,delthvu_c,         &
!$OMP        cape,cape_c,cin,cin_c,ql_ad,ql_ad_c,ntpar,zh,zhpar,ntml,          &
!$OMP        nlcl,z_lcl_uv,z_full_c,cldbase_opt_sh,w_copy,                     &
!$OMP        cvdiag_sh_wtest,cumulus,l_shallow,entrain_coef,                   &
!$OMP        l_reset_neg_delthvu)                                              &
!$OMP private(ii,i,j,k)

!CDIR NODEP
!$OMP do SCHEDULE(STATIC)
do ii=1,nunstable
  i = index_i(ii)
  j = index_j(ii)
  select case (iwtcs_diag2)
  case (1)
    !-------------------------------------------
    ! congestus points through deep (MF) scheme
    !-------------------------------------------
    if (conv_type(i,j) >= 3) then
      conv_type(i,j) = 8
    end if
  case (2)
    !-------------------------------------------
    ! conv_type=4 through conv_type=3
    !-------------------------------------------
    if (conv_type(i,j) == 4) then
      conv_type(i,j) = 3
    end if
  case (3)
    !-------------------------------------------
    ! conv_type=1,3,4 through conv_type=2
    !-------------------------------------------
    if (conv_type(i,j) == 1 .or.                                               &
        conv_type(i,j) == 3 .or.                                               &
        conv_type(i,j) == 4) then
      conv_type(i,j) = 2
    end if
  end select

  !-------------------------------------------
  ! Make sure shallow_c set appropriately
  !-------------------------------------------
  select case (conv_type(i,j))
  case (1:4)
    shallow_c(ii) = .true. ! types 1-4 go through shallow logical
  case DEFAULT
    shallow_c(ii) = .false.
  end select

end do
!$OMP end do NOWAIT
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)
  if (cumulus_c(ii)) then
    !-------------------------------------------------------------------
    ! Set mixed layer depth to z_lcl
    !-------------------------------------------------------------------
    if (p_lcl(ii)  <   (P_theta_lev_c(ii,nlcl_c(ii)+1))) then
      !-------------------------------------------------------------------
      ! If LCL is diagnosed in the upper half of the layer set z_lcl to
      ! the height of the upper layer interface
      ! (in code above LCL is always set to the lower interface).
      !-------------------------------------------------------------------
      nlcl_c(ii) = nlcl_c(ii)+1
    end if

  else      ! not cumulus

    !---------------------------------------------------------------------
    !      If not cumulus, reset parameters to within bl_levels
    !---------------------------------------------------------------------
    if (ntml_c(ii)  >   mbl) then
      ntml_c(ii)  = mbl
      zh_c(ii)    = z_half_c(ii,mbl+1)
    end if

  end if        ! test on cumulus

  !---------------------------------------------------------
  ! nlcl is not permitted to be less than nlcl_min
  !---------------------------------------------------------

  nlcl_c(ii) = max(nlcl_min(ii), nlcl_c(ii))

  if ( ntml_c(ii)-nlcl_c(ii) <=2 ) then
    ! Cloud layer now too shallow so rediagnose as well-mixed
    cumulus_c(ii) = .false.
    shallow_c(ii) = .false.
    l_congestus(i,j) = .false.
  end if

end do       ! ii loop
!$OMP end do NOWAIT
!----------------------------------------------------------------------
! Expand back up to full arrays
!----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)

  delthvu(i,j) = delthvu_c(ii)
  cape(i,j)    = CAPE_c(ii)
  cin(i,j)     = CIN_c(ii)
  ql_ad(i,j)   = ql_ad_c(ii)

  ntpar(i,j)   = ntml_c(ii)   ! holds parcel top level
  zh(i,j)      = zh_c(ii)     ! holds parcel top at this point
  zhpar(i,j)   = zh_c(ii)

  if (cumulus_c(ii)) then
    ntml(i,j)  = nlcl_c(ii)   ! now sets to LCL
    zh(i,j)    = z_half_c(ii,nlcl_c(ii)+1)    ! reset to zlcl
  else
    ntml(i,j)  = ntml_c(ii)
  end if
  nlcl(i,j)  = nlcl_c(ii)

  z_lcl_uv(i,j)   = z_full_c(ii,nlcl_c(ii))    ! LCL height on uv

  ! Cumulus points always have nlcl <= MBL from previous checks
  if (nlcl_c(ii)  >   mbl) then    ! only applied if not cumulus
    nlcl(i,j)    = mbl
    z_lcl_uv(i,j)   = z_full_c(ii,mbl-1)
  end if

end do       ! ii loop
!$OMP end do NOWAIT
!      If cumulus has been diagnosed but delthvu is negative, reset
!      cumulus and L_shallow to false but leave zh and ntml at LCL

if ( l_reset_neg_delthvu ) then
!$OMP do SCHEDULE(STATIC)
  do ii=1, nunstable
    if (cumulus_c(ii) .and. delthvu_c(ii)  <=  zero) then
      i = index_i(ii)
      j = index_j(ii)
      conv_type(i,j)   = 0
      cumulus_c(ii)    = .false.
      shallow_c(ii)    = .false.
    end if
  end do
!$OMP end do NOWAIT
end if


! Additional requirement for grey zone (so high res) option, that there
! is low level convergence (ie significant positive w) before convection
! scheme is triggered (note deep and mid-level are disabled)
if (cldbase_opt_sh == sh_grey_closure) then
!$OMP do SCHEDULE(STATIC)
  do ii=1, nunstable
    ! Start by making all convection shallow
    shallow_c(ii) = cumulus_c(ii)
    i = index_i(ii)
    j = index_j(ii)
    k = nlcl_c(ii)+1
    if ( cumulus_c(ii) .and. w_copy(i,j,k) < cvdiag_sh_wtest) then
      shallow_c(ii) = .false.   ! disables convection scheme
    end if
  end do
!$OMP end do NOWAIT
end if

! Expand back shallow and cumulus arrays
!$OMP do SCHEDULE(STATIC)
do ii=1, nunstable
  i = index_i(ii)
  j = index_j(ii)
  cumulus(i,j)   = cumulus_c(ii)
  l_shallow(i,j) = shallow_c(ii)

  ! No points with altered entrainment go through shallow scheme.

  if (l_shallow(i,j) .and. entrain_coef(i,j) > zero) then
    entrain_coef(i,j) = -99.0_r_bl
  end if

end do
!$OMP end do NOWAIT
!$OMP end PARALLEL
!  If cumulus has been diagnosed but delthvu is negative, reset
!      congestus to false but leave zh and ntml at LCL

if (iconv_congestus == 1) then     ! mass flux option only
  do ii=1, nunstable
    i = index_i(ii)
    j = index_j(ii)
    if (cumulus_c(ii) .and. delthvu_c(ii)  <=  zero) then
      l_congestus(i,j)  = .false.
      l_congestus2(i,j) = .false.
    end if
  end do       ! ii loop
end if            ! test on congestus diagnosis

!-----------------------------------------------------------------------
! Finalise inversion thickness calculation
!-----------------------------------------------------------------------
if (forced_cu >= on .or. bl_res_inv /= off) then
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED(nunstable,index_i,index_j,cumulus,dzh,ntml,ntml_save,zh,          &
!$OMP        zh_save,mbl,z_half_c,zh_itop_c,model_levels,ntinv_c,              &
!$OMP        qcl_inv_top,z_full_c,ql_parc_c)                                   &
!$OMP private(ii,i,j,k,dz)
  do ii=1, nunstable
    i = index_i(ii)
    j = index_j(ii)
    if (cumulus(i,j)) then
      dzh(i,j)  = rmdi
    else
      ! need to reset parcel top to inversion base
      ! (ie without inversion included)
      ntml(i,j) = ntml_save(ii)
      zh(i,j)   = zh_save(ii)
      if (ntml(i,j) > mbl) then
        ntml(i,j) = mbl
        zh(i,j)   = z_half_c(ii,mbl+1)
      end if
      ! restrict inversion thickness to be less than the BL depth
      dzh(i,j) = min( zh(i,j), zh_itop_c(ii) - zh(i,j) )
    end if
    ! forced_cu options need parcel water content at the top of the inversion
    if ( dzh(i,j) > zero ) then
      zh_itop_c(ii) = zh(i,j)+dzh(i,j)
    else
      ! used as a lower limit on in-cloud qcl near cu cloud base
      zh_itop_c(ii) = zh(i,j)+300.0_r_bl
    end if
    k = ntml(i,j)+1
    ntinv_c(ii) = k
    do while (z_half_c(ii,k+1) < zh_itop_c(ii) .and. k < model_levels-1)
      ntinv_c(ii) = k+1  ! marks top level within inversion
      k=k+1
    end do
    k=ntinv_c(ii)
    dz = z_full_c(ii,k+1)-z_full_c(ii,k)
    qcl_inv_top(i,j) = ql_parc_c(ii,k) + (zh_itop_c(ii)-z_full_c(ii,k)) *      &
                                     (ql_parc_c(ii,k+1)-ql_parc_c(ii,k))/dz
  end do
!$OMP end PARALLEL do
end if


!-----------------------------------------------------------------------
!     SCM  Diagnostics all versions
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine conv_diag_comp_6a

end module conv_diag_comp_6a_mod
