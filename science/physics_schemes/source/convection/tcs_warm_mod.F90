! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module containing tcs warm rain subroutine
!
module tcs_warm_mod

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

use ereport_mod, only: ereport
use umPrintMgr,  only: umPrint, umMessage
use um_types, only: real_umphys

implicit none

!
! Description:
!   This routine calculates convective tendencies for warm
!   precipitating convection.
!
! Method:
!   Currently a development version, this uses a turbulence-based
!   approach rather than a traditional mass flux approach.
!   <reference to documentation to go here, once available>
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!

character(len=*), parameter, private :: ModuleName='TCS_WARM_MOD'

contains

subroutine tcs_warm(                                                           &
     nlev,n_xx,ntra,trlev ,n_cca_lev                                           &
   , ntml, ntpar, conv_type                                                    &
   , L_calc_dxek,L_tracer                                                      &
   , timestep                                                                  &
   , delthvu,ql_ad,uw0,vw0                                                     &
   , wstar_dn,wthvs,zlcl_uv,ztop_uv                                            &
   , p_layer_centres,p_layer_boundaries                                        &
   , exner_layer_centres                                                       &
   , exner_layer_boundaries                                                    &
   , z_theta,z_rho, r_rho                                                      &
   , rho, rho_theta, r2rho, r2rho_theta                                        &
   , dr_across_rh,dr_across_th                                                 &
   , q_mix, theta, u,v                                                         &
   , qse                                                                       &
                              ! Inout fields (PC2 ones not used)
   , cf_liquid                                                                 &
   , qcl, tracer                                                               &
                              ! Output fields
   , cape_out,cclwp,ccw,cca                                                    &
   , dbcfbydt,dcffbydt,dcflbydt                                                &
   , dqbydt,dqcfbydt,dqclbydt                                                  &
   , dthbydt,dubydt,dvbydt,dtrabydt                                            &
   , iccb,icct,lcca,lcbase,lctop                                               &
   , rain,snow, up_flux                                                        &
   , uw_shall,vw_shall                                                         &
   , wqt,wthetal,wthetavl,wthetav, wh, wql                                     &
   , wstar_up,mb1,mb2                                                          &
   , tcw, cca_2d)

  ! Modules Used:

use tcs_classes, only:                                                         &
  tcs_allocate, tcs_deallocate, & ! overloaded allocation routines
  similarity, cloud_input,      & ! derived types
  assign_interface_input, assign_cloud_input
use tcs_constants,    only:                                                    &
   g, lc, r, ra2, c_virtual, rv, gamma_dry,                                    &
   lc_o_cp, recip_epsilon, repsilon
use tcs_calc_scales,  only:                                                    &
    calc_scales_warm
use tcs_cloudbase,    only:                                                    &
    calc_cloudbase
use tcs_precip,       only:                                                    &
    calc_precip
use tcs_inversion,    only:                                                    &
    calc_inversion_warm
use tcs_similarity,   only:                                                    &
    calc_similarity
use tcs_qlup,         only:                                                    &
    calc_qlup_warm
use tcs_wql,          only:                                                    &
    calc_wql_warm
use tcs_grad_h,       only:                                                    &
    calc_grad_h
use tcs_turb_fluxes,  only:                                                    &
    calc_turb_fluxes
use tcs_grad_flux,    only:                                                    &
    calc_grad_flux
use tcs_grad_stress,  only:                                                    &
    calc_grad_stress
use tcs_base_stress,  only:                                                    &
    calc_base_stress
use tcs_cmt_incr_6a,  only:                                                    &
    calc_cmt_incr
use tcs_pc2,          only:                                                    &
    calc_pc2
use tcs_parameters_warm, only:                                                 &
    asmooth_inv,  bsmooth_inv


use tcs_common_warm,  only:                                                    &
   scales, cb, cb_p1, cb_m1, inv, inv_p1, inv_m1

use cv_run_mod,       only:                                                    &
    l_mom

use cv_stash_flg_mod, only:                                                    &
   flg_uw_shall                                                                &
   ,flg_vw_shall


use model_domain_mod,       only: model_type, mt_single_column
use s_scmop_mod,            only: default_streams,                             &
                                 t_inst,d_wet,d_point

implicit none


!------------------------------------------------------------------
! Subroutine Arguments
!------------------------------------------------------------------
!
! Arguments with intent in:
!

integer, intent(in), target ::                                                 &
   nlev                                                                        &
                            ! No. of model layers
   , n_xx                                                                      &
                            ! No. of warm convection points
   , ntra                                                                      &
                            ! No. of tracer fields
   , trlev
                            ! No. of model levels on
                            ! which tracers are included

integer, intent(in) ::                                                         &
     n_cca_lev              ! No. of convective cloud
                            ! amount levels (1 for 2D, nlevs for 3D)

integer, intent(in) ::                                                         &
   ntml(n_xx)                                                                  &
                            ! Top level of surface mixed layer defined
                            ! relative to theta,q grid
   , ntpar(n_xx)
                            ! Top level of initial parcel ascent in BL
                            ! scheme defined relative to theta,q grid

integer, intent(in) ::                                                         &
   conv_type(n_xx)
                            ! Integer index describing convective type:
                            !    1=non-precipitating shallow
                            !    2=drizzling shallow
                            !    3=warm congestus

logical, intent(in) ::                                                         &
   L_calc_dxek                                                                 &
                            ! Switch for calculation of condensate increment
                            ! (PC2)
   , L_tracer
                            ! Switch for inclusion of tracers

real(kind=real_umphys), intent(in) ::                                          &
   timestep           ! Model timestep (s)

real(kind=real_umphys), intent(in) ::                                          &
   delthvu(n_xx)                                                               &
                            ! Integral of undilute parcel buoyancy over
                            ! convective cloud layer (Km)
   , ql_ad(n_xx)                                                               &
                            ! adiabatic liquid water content at inversion
                            ! (kg/kg)
   , uw0(n_xx)                                                                 &
                            ! U-comp of surface stress (N/m2)
   , vw0(n_xx)                                                                 &
                            ! V-comp of surface stress (N/m2)
   , wstar_dn(n_xx)                                                            &
                            ! Mixed layer convective velocity scale
   , wthvs(n_xx)                                                               &
                            ! Surface flux of THV (Km/s)
   , zlcl_uv(n_xx)                                                             &
                            ! Lifting condensation level defined for
                            ! the uv grid (m)
   , ztop_uv(n_xx)
! Top of cloud layer defined for the uv
! grid (m)

real(kind=real_umphys), intent(in)    ::                                       &
     p_layer_centres(n_xx,0:nlev)                                              &
                            ! Pressure (Pa)
   , p_layer_boundaries(n_xx,0:nlev)                                           &
                            ! Pressure at half level above
                            ! p_layer_centres (Pa)
   , exner_layer_centres(n_xx,0:nlev)                                          &
                            !Exner
   , exner_layer_boundaries(n_xx,0:nlev)
! Exner at half level above
! exner_layer_centres
real(kind=real_umphys), intent(in)    ::                                       &
   z_theta(n_xx,nlev)                                                          &
                            ! height of theta levels(m)
   , z_rho(n_xx,nlev)                                                          &
                            ! height of rho levels (m)
   , r_rho(n_xx,nlev)                                                          &
                            ! radius of rho levels (m)
   , rho(n_xx,nlev)                                                            &
                            ! density on rho levels kg/m3
   , rho_theta(n_xx,nlev)                                                      &
                            ! density on theta levels kg/m3
   , r2rho(n_xx,nlev)                                                          &
                            ! r2*rho rho levels (kg/m)
   , r2rho_theta(n_xx,nlev)                                                    &
                            ! r2*rho theta levels (kg/m)
   , dr_across_rh(n_xx,nlev)                                                   &
                            ! dr across rho levels  (m)
   , dr_across_th(n_xx,nlev)
! dr across theta levels (m)

real(kind=real_umphys), intent(in), target    ::                               &
   q_mix(n_xx,nlev)                                                            &
                            ! q as mixing ratio(kg/kg)
   , theta(n_xx,nlev)                                                          &
                            ! Model potential temperature (K)
   , u(n_xx,nlev)                                                              &
                            ! Model U field (m/s)
   , v(n_xx,nlev)                                                              &
                            ! Model V field (m/s)
   , qse(n_xx,nlev)
! Saturation mixing ratio of cloud
! environment (kg/kg)

!
! Arguments with intent INOUT:
!
! Variables used only by PC2

real(kind=real_umphys), intent(in out) ::                                      &
     cf_liquid(n_xx,nlev)                                                      &
                            ! Liq water cloud volume
   , qcl(n_xx,nlev)
                            ! Liq condensate mix ratio (kg/kg)

! Variable only used if LTRACER

real(kind=real_umphys), intent(in out)    :: tracer(n_xx,trlev,ntra)
!Model tracer fields (kg/kg)

!
! Arguments with intent out:
!

real(kind=real_umphys), intent(out) ::                                         &
   cape_out(n_xx)                                                              &
                            ! Saved convective available potential energy
                            ! for diagnostic output (J/kg)
                            ! not used in 5A
   , cclwp(n_xx)                                                               &
                            ! Condensed water path (kg/m2)
   , ccw(n_xx,nlev)                                                            &
                            ! Convective cloud liquid water on model
                            ! levels (g/kg)
   , cca(n_xx,n_cca_lev)
                            ! Convective cloud amount on model levels
                            ! (fraction)

real(kind=real_umphys), intent(out) ::                                         &
   dbcfbydt(n_xx,nlev)                                                         &
                            ! Increments to total cld volume due to
                            ! convection(/s)
   , dcffbydt(n_xx,nlev)                                                       &
                            ! Increments to ice cloud volume due to
                            !  convection(/s)
   , dcflbydt(n_xx,nlev)                                                       &
                            ! Increments to liq cloud volume due to
                            ! convection(/s)
   , dqbydt(n_xx,nlev)                                                         &
                            ! Increments to q (water vap mixing ratio)
                            ! due to convection (kg/kg/s)
   , dqcfbydt(n_xx,nlev)                                                       &
                            ! Increments to ice condensate due to
                            ! convection (kg/kg/s)
   , dqclbydt(n_xx,nlev)                                                       &
                            ! Increments to liq condensate due to
                            ! convection (kg/kg/s)
   , dthbydt(n_xx,nlev)                                                        &
                            ! Increments to potential temp. due to
                            ! convection (K/s)
   , dubydt(n_xx,nlev+1)                                                       &
                            ! Increments to U due to CMT (m/s2)
   , dvbydt(n_xx,nlev+1)                                                       &
                            ! Increments to V due to CMT (m/s2)
   , dtrabydt(n_xx,nlev,ntra)
!Increment to tracer due to
! convection (kg/kg/s)

integer, intent(out) ::                                                        &
   iccb(n_xx)                                                                  &
                            ! Convective cloud base level (m)
   , icct(n_xx)
! Convective cloud top level (m)

real(kind=real_umphys), intent(out) ::                                         &
   lcca(n_xx)   ! Lowest conv. cloud amt. (%)

integer, intent(out) ::                                                        &
   lcbase(n_xx)                                                                &
                            ! Lowest conv. cloud base level (m)
   , lctop(n_xx)
! Lowest conv. cloud top level (m)

real(kind=real_umphys), intent(out) ::                                         &
   rain(n_xx)                                                                  &
                            ! Surface convective rainfall (kg/m2/s)
   , snow(n_xx)
! Surface convective snowfall (kg/m2/s)

real(kind=real_umphys), intent(out) ::                                         &
   up_flux(n_xx,nlev)  ! mass flux (Pa/s)

real(kind=real_umphys), intent(out) ::                                         &
   uw_shall(n_xx,nlev)                                                         &
                            ! X-comp. of stress from congestus convection
                            !(kg/m/s2)
   , vw_shall(n_xx,nlev)                                                       &
                            ! Y-comp. of stress from congestus convection
                            !(kg/m/s2)
   , wqt(n_xx,nlev)                                                            &
                            ! w'qt' flux  (m/s kg/kg)
   , wthetal(n_xx,nlev)                                                        &
                            ! w'thetal' flux  (Km/s)
   , wthetavl(n_xx,nlev)                                                       &
                            ! w'thetavl' flux  (Km/s)
   , wthetav(n_xx,nlev)                                                        &
                            ! w'thetav' flux  (Km/s)
   , wh(n_xx,nlev)                                                             &
                            ! w'h' flux  (K m/s)
   , wql(n_xx,nlev)
                            ! w'ql'   (kg m/s/kg)

real(kind=real_umphys), intent(out) ::                                         &
   wstar_up(n_xx)                                                              &
                            ! convective layer velocity scale
   , mb1(n_xx)                                                                 &
                            ! cloud base massflux using parametrization
                            !  method 1
   , mb2(n_xx)
                  ! cloud base massflux using parametrization method 2

real(kind=real_umphys), intent(out) ::                                         &
   tcw(n_xx)                                                                   &
                            ! Total condensed water(kg/m2/s)
   , cca_2d(n_xx)
                            ! 2D convective cloud amount (%)

!--------------------------------
! Local variables
!--------------------------------

character(len=*), parameter ::  RoutineName = 'TCS_WARM'
! Commented out variables used in commented out code
!character(len=errormessagelength) :: Message
!integer :: ErrorStatus ! Return code:
!   0 = Normal exit
! +ve = Fatal Error
! -ve = Warning

real(kind=real_umphys) ::                                                      &
   wup(n_xx,nlev)                                                              &
                            ! ensemble vertical velocity (m/s)
   ,  wtheta(n_xx,nlev)                                                        &
                            ! wtheta flux
   ,  wq(n_xx,nlev)                                                            &
                            ! wq flux
   ,  qlup(n_xx,nlev)
                            ! ensemble qlup   (kg/kg)

real(kind=real_umphys) ::                                                      &
   wtracer(n_xx,trlev,ntra) ! wtracer flux on uv levels (kgm/kg/s)

!
! Limit nlev loop to those levels actually required using ntpar
! diagnosed in conv_diag
!
integer :: ntpar_max         ! max ntpar value
integer :: ntml_max          ! max ntml value

!-------------------------------------
! Local arrays - In cloud arrays
!-------------------------------------

type(cloud_input) :: cld_in

!
! values evaluated for cloud levels
!
real(kind=real_umphys) ::                                                      &
   wup_cld(n_xx,nlev)                                                          &
                            ! cumulus ensemble vertical velocity
   , wup_h_cld(n_xx,nlev)                                                      &
                            ! cumulus ensemble vertical velocity
   , ql_up_cld(n_xx,nlev)                                                      &
                            ! cumulus ensemble ql
   , dwql_dz(n_xx,nlev)                                                        &
                            ! dwql/dz all level
   , h_cld(n_xx,nlev)                                                          &
                            ! moist static energy / cp in cloud
   , t_cld(n_xx,nlev)                                                          &
                            ! temperature in degrees Celcius in cloud
   , dqsatdt_cld(n_xx,nlev)                                                    &
                            ! dqsat/dt
   , mf_cld(n_xx,nlev)                                                         &
                            ! mass flux in cld (theta levels)
   , mf_h_cld(n_xx,nlev)                                                       &
                            ! mass flux in cld (rho levels)
   , kdthetavdz(n_xx,nlev)                                                     &
                            ! Kdthetav/dz
   , kdhdz(n_xx,nlev)                                                          &
                            ! Kdh/dz
   , kdtracerdz(n_xx,nlev,ntra)                                                &
                            ! K*dtracer/dz in cloud
   , wq_cld(n_xx,nlev)                                                         &
                            ! wq in cloud
   , wqt_cld(n_xx,nlev)                                                        &
                            ! wqt in cloud
   , wql_cld(n_xx,nlev)                                                        &
                            ! wql in cloud with cb & inv values
   , wtheta_cld(n_xx,nlev)                                                     &
                            ! wtheta in cloud
   , wthetal_cld(n_xx,nlev)                                                    &
                            ! wthetal in cloud
   , wh_cld(n_xx,nlev)                                                         &
                            ! wh in cloud
   , wthetavl_cld(n_xx,nlev)                                                   &
                            ! wthetavl in cloud
   , wthetav_cld(n_xx,nlev)                                                    &
                            ! wthetavl in cloud
   , ztheta_cld(n_xx,nlev)                                                     &
                            ! height of model theta levels in cloud
   , thetav_cld(n_xx,nlev)                                                     &
                            ! thetav on cloud levels
   , r2rho_cld(n_xx,nlev)                                                      &
                            ! rho on uv cloud levels
   , r2rho_theta_cld(n_xx,nlev)                                                &
                            ! rho on theta cloud levels
   , u_cld(n_xx,nlev)                                                          &
                            ! u on cloud levels
   , v_cld(n_xx,nlev)                                                          &
                            ! v on cloud levels
   , uw_cld(n_xx,nlev)                                                         &
                            ! uw on cloud levels
   , vw_cld(n_xx,nlev)                                                         &
                            ! vw on cloud levels
   , dr_across_rh_cld(n_xx,nlev)                                               &
                            ! thickness of rho cloud levels
   , dr_across_th_cld(n_xx,nlev)                                               &
                            ! thickness of theta cloud levels
   , wqr_cld(n_xx,nlev)                                                        &
                            ! flux of rain water (kg m /kg/s)
   , precip_product_hcld(n_xx,nlev)                                            &
                            !precipitation production
                            ! (kg/m3/s) uv levels
   , precip_product_cld(n_xx,nlev)                                             &
                            ! precipitation production
                            ! (kg/m3/s) theta levels
   , wtracer_cld(n_xx,trlev,ntra)                                              &
                            ! wtracer in cloud values(kgm/kg/s)
   , coeff_a, coeff_b
                            ! Dummy coefficients to simplify code


!
! Flux gradients on all model levels
!
real(kind=real_umphys) ::                                                      &
   dwthetal_dz(n_xx,nlev)                                                      &
                            ! dwthetal/dz
   , dwqt_dz(n_xx,nlev)                                                        &
                            ! dwqt/dz
   , dwthetavl_dz(n_xx,nlev)                                                   &
                            ! dwthetavl/dz
   , dwtracer_dz(n_xx,trlev,ntra)
                            ! dwtracer/dz

!
! Cloud base fluxes
!
real(kind=real_umphys) ::                                                      &
   wthetav_cb(n_xx)                                                            &
                            ! wthetav at cloud base
   , wthetav_cb_m2(n_xx)                                                       &
                            ! wthetav at cloud base using method 2
   , wtheta_plus_cb(n_xx)                                                      &
                            ! wtheta at cloud base (in cloud value)
   , wthetal_cb(n_xx)                                                          &
                            ! wthetal at cloud base
   , wq_plus_cb(n_xx)                                                          &
                            ! wq at cloud base (in cloud value)
   , wql_cb(n_xx)                                                              &
                            ! wql at cloud base
   , wqt_cb(n_xx)                                                              &
                            ! wqt at cloud base
   , wtracer_cb(n_xx,ntra)
                            ! wtracer at cloud base

!
! Arrays for cloud base calculations
!
real(kind=real_umphys) ::                                                      &
   dthetal_cb(n_xx)                                                            &
                            ! dthetal across cloud base
   , dthetav_cb(n_xx)                                                          &
                            ! dthetav across cloud base
   , dqt_cb(n_xx)                                                              &
                            ! dqt across cloud base
   , dz_cb(n_xx)                                                               &
                            ! dz across cloud base
   , dtracer_cb(n_xx,ntra)
                            ! dtracer across cloud base

real(kind=real_umphys) ::                                                      &
   du_cb(n_xx)                                                                 &
                            ! du across cloud base
   , dv_cb(n_xx)
                            ! dv across cloud base
!
! Arrays for across cloud
!
integer ::                                                                     &
   ncld_thlev(n_xx)       ! number of theta levels in cloud

real(kind=real_umphys) ::                                                      &
   dqt_cld(n_xx)                                                               &
                            ! qt across cloud
   , dthetal_cld(n_xx)                                                         &
                            ! thetal across cloud
   , dp_cld(n_xx)                                                              &
                            ! p across cloud
   , dthetavdz_ad_lcl(n_xx)                                                    &
                            ! dthetav/dz for a moist adiabat at lcl
   , dthetav_cld(n_xx)
                            ! needed for one option of inv

!
! Base of inversion fluxes
!
real(kind=real_umphys) ::                                                      &
   wtheta_inv(n_xx)                                                            &
                            ! wtheta  across inversion
   ,wthetal_inv(n_xx)                                                          &
                            ! wthetal across inversion
   ,wq_inv(n_xx)                                                               &
                            ! flux of wq across inversion
   ,wqt_inv(n_xx)                                                              &
                            ! flux of wqt across inversion
   ,wql_inv(n_xx)                                                              &
                            ! flux of wql across inversion
   ,wthetavl_inv(n_xx)                                                         &
                            ! wthetavl across inversion
   ,wthetav_inv(n_xx)                                                          &
                            ! wthetav across inversion
   ,wh_inv(n_xx)                                                               &
                            ! wh across inversion
   ,wqr_inv(n_xx)                                                              &
                            ! wqr flux at the base of inversion
   ,precip_product_inv(n_xx)                                                   &
                            ! integral of precip production for
                            ! inversion
   ,wtracer_inv(n_xx,ntra)
                            ! wtracer flux across inversion

!
! Arrays for across inversion calculation
!
real(kind=real_umphys) ::                                                      &
   dz_inv(n_xx)
                            ! depth of inversion (m)

!
! Arrays for similarity functions
!
type(similarity) :: sim

!
! Arrays for CMT calculation
!
integer ::                                                                     &
   nlcl_uv(n_xx)                                                               &
                            ! Level index for LCL
   , ntop_uv(n_xx)
                            ! Level index for top of layer

real(kind=real_umphys) ::                                                      &
   uw(n_xx,nlev)                                                               &
                            ! U- comp stress profile (m2s-2)
   , vw(n_xx,nlev)
                            ! V-comp stress profile (m2s-2)

!
! Arrays for CMT calculation - possible other calculations
!
real(kind=real_umphys) ::                                                      &
   dr_across_th1(n_xx,nlev)       ! dr across theta levels

!
! Arrays for precipitation calculations
!
real(kind=real_umphys) ::                                                      &
   wqr(n_xx,nlev)                                                              &
                            ! flux of rain water (kg m /kg/s)
   , precip_product(n_xx,nlev)
                            ! precipitation production/density (/s)

!
!   required by check on -ve q
!
real(kind=real_umphys) :: qMinInColumn(n_xx)     ! Minimum value for q in column
                               ! (kg/kg)
real(kind=real_umphys) :: tempx(n_xx)            ! work array

real(kind=real_umphys), parameter :: qmin = 1.0e-8 ! Global minimum allowed Q
!
! temporary storeage arrays
!
real(kind=real_umphys) ::                                                      &
   wup2_h, wup2, dmass                                                         &
   , factor                                                                    &
   , t_lcl_temp, temp1, temp2

!
! Loop counters
!
integer :: i,k,ktra,klev,j
integer :: max_cldlev, ilev, ibase, kinv, max_cldtrlev

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle


!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!                                                                      !
!        <----------- End of variable declarations ----------->        !
!                                                                      !
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

!-----------------------------------------------------------------------
!
!  Congestus scheme assumes levels ntml+1 to ntpar are in cloud
!  The parametrisation scheme assumes the atmosphere is hydrostatic and
!  Boussinesq. (This is contrary to the assumptions made in the model
!  dynamics.)
!
!  Diagram of congestus cloud levels
!
!
!       ----------------------------  ntpar + 1
!
!       ++++++++++++++++++++++++++++   Inversion on uv level (ntpar+1)
!
!       ----------------------------  ntpar
!
!       + + + + + + + + + + + + + +
!
!       ----------------------------
!          "          "
!          "          "               Centre of cloud
!          "          "
!       + + + + + + + + + + + + + +
!
!       ----------------------------  ntml + 1   (plus level)
!
!       ++++++++++++++++++++++++++++  cloud base on uv level (ntml+1)
!
!       ----------------------------  ntml       (minus level)
!
!
!  Key
!  ----
!
!  -------   theta levels (a.k.a half levels)
!
!  + + + +   uv levels (a.k.a. rho/flux/full levels)
!
!

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1.0 Allocate and initialize local storage
!-----------------------------------------------------------------------
!-----------------------------------------------------------------
! work out maximum level values - used to reduce allocation
! and loops
!-----------------------------------------------------------------
ncld_thlev(:) = ntpar(:) - ntml(:) + 1
                                ! number of theta levels in cloud
ntpar_max  = maxval(ntpar)+1    ! max top of parcel ascents
ntml_max   = maxval(ntml)       ! max top of subcloud layer
max_cldlev = maxval(ncld_thlev)        ! max number of cloud levels
max_cldtrlev = min(max_cldlev,trlev-1) ! max number of tracer levels

!-----------------------------------------------------------------
! Call allocation routines
!-----------------------------------------------------------------
call tcs_allocate(cb, n_xx=n_xx, ntra=ntra)     ! interface_input
call tcs_allocate(cb_p1, n_xx=n_xx, ntra=ntra)  ! interface_input
call tcs_allocate(cb_m1, n_xx=n_xx, ntra=ntra)  ! interface_input

call tcs_allocate(inv, n_xx=n_xx, ntra=ntra)    ! interface_input
call tcs_allocate(inv_p1, n_xx=n_xx, ntra=ntra) ! interface_input
call tcs_allocate(inv_m1, n_xx=n_xx, ntra=ntra) ! interface_input

call tcs_allocate(cld_in,                     & ! cloud_input
   n_xx=n_xx, nlev=max_cldlev+1, ntra=ntra, initval=0.0)
call tcs_allocate(sim,                        & ! similarity
   n_xx=n_xx, nlev=cld_in%nlev, initval=0.0)
call tcs_allocate(scales, n_xx=n_xx)            ! scales_conv

!-----------------------------------------------------------------
! Initialise some PC2 variables for safety.
!-----------------------------------------------------------------
dqclbydt(:,:) = 0.0
dqcfbydt(:,:) = 0.0
dbcfbydt(:,:) = 0.0
dcflbydt(:,:) = 0.0
dcffbydt(:,:) = 0.0

!-----------------------------------------------------------------
! Initialize in-cloud diagnosed fields
!-----------------------------------------------------------------
thetav_cld(:,:)      =  0.0
h_cld(:,:)           =  0.0
t_cld(:,:)           =  0.0
ztheta_cld(:,:)      =  0.0
r2rho_theta_cld(:,:) =  0.0

!-----------------------------------------------------------------
! Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
! increment arrays and cca
!-----------------------------------------------------------------
dthbydt(:,:) = 0.0
dqbydt(:,:)  = 0.0
if (L_mom) then
  dubydt(:,:) = 0.0
  dvbydt(:,:) = 0.0
end if  ! L_mom

!-----------------------------------------------------------------
! Initialise tracer increments if tracers present.
!-----------------------------------------------------------------
if (L_tracer) then
  dtrabydt(:,:,:) = 0.0
end if  ! L_tracer

!-----------------------------------------------------------------
! Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------
ccw(:,:)    = 0.0
rain(:)     = 0.0
snow(:)     = 0.0
cclwp(:)    = 0.0
lcca(:)     = 0.0
cape_out(:) = 0.0
tcw(:)      = 0.0
if (L_mom) then
  if (flg_uw_shall) then
    uw_shall(:,:) = 0.0
  end if
  if (flg_vw_shall) then
    vw_shall(:,:) = 0.0
  end if
end if  ! L_mom

!-----------------------------------------------------------------
! diagnostics on full levels (i.e. all model levels)
! calculations on just cloud levels  - need 2 arrays for each ?
!-----------------------------------------------------------------
uw(:,:)             = 0.0
vw(:,:)             = 0.0
wtheta(:,:)         = 0.0
wq(:,:)             = 0.0
wqt(:,:)            = 0.0
wqr(:,:)            = 0.0
wql(:,:)            = 0.0
wthetavl(:,:)       = 0.0
wthetav(:,:)        = 0.0
wh(:,:)             = 0.0
wthetal(:,:)        = 0.0
dwql_dz(:,:)        = 0.0
up_flux(:,:)        = 0.0
wup(:,:)            = 0.0
qlup(:,:)           = 0.0
precip_product(:,:) = 0.0
if (L_tracer) then
  wtracer(:,:,:) = 0.0
end if  ! L_tracer

!-----------------------------------------------------------------
! Distance between rho levels (not thickness of levels in the case
! of k=1)
!-----------------------------------------------------------------
do i=1,n_xx
  do k=1,nlev-1
    dr_across_th1(i,k) = r_rho(i,k+1) - r_rho(i,k)
  end do
end do

!-----------------------------------------------------------------------
!
! 2.0  Interface level and jump calculations
!
!-----------------------------------------------------------------------

! Working note: it would be useful to give a detailed description
! of what these layers represent, with particular attention to p
! and exner on theta and rho levels

! Working note: Some further jumps are calculated within the
! subroutines, e.g. tcs_cloudbase - these should be pulled out
! here

!-----------------------------------------------------------------
! Set arrays around cloud base
!-----------------------------------------------------------------
call assign_interface_input(cb_m1, ntml, theta,                                &
   q_mix, qse, p_layer_centres, p_layer_boundaries,                            &
   exner_layer_centres, exner_layer_boundaries,                                &
   z_theta, z_rho, rho,tracer )

call assign_interface_input(cb, ntml+1, theta,                                 &
   q_mix, qse, p_layer_centres, p_layer_boundaries,                            &
   exner_layer_centres, exner_layer_boundaries,                                &
   z_theta, z_rho, rho,tracer )

call assign_interface_input(cb_p1, ntml+2, theta,                              &
   q_mix, qse, p_layer_centres, p_layer_boundaries,                            &
   exner_layer_centres, exner_layer_boundaries,                                &
   z_theta, z_rho, rho,tracer )

!-----------------------------------------------------------------
! Set arrays around the base of the inversion
!-----------------------------------------------------------------
call assign_interface_input(inv_m1, ntpar, theta,                              &
   q_mix, qse, p_layer_centres, p_layer_boundaries,                            &
   exner_layer_centres, exner_layer_boundaries,                                &
   z_theta, z_rho, rho,tracer )

call assign_interface_input(inv, ntpar+1, theta,                               &
   q_mix, qse, p_layer_centres, p_layer_boundaries,                            &
   exner_layer_centres, exner_layer_boundaries,                                &
   z_theta, z_rho, rho,tracer )

call assign_interface_input(inv_p1, ntpar+2, theta,                            &
   q_mix, qse, p_layer_centres, p_layer_boundaries,                            &
   exner_layer_centres, exner_layer_boundaries,                                &
   z_theta, z_rho, rho,tracer )

!-----------------------------------------------------------------
! calculate delta values across cloud
!-----------------------------------------------------------------
dthetal_cld(:) = inv%theta(:) - cb_m1%theta(:)
dqt_cld(:)     = inv%q_mix(:) - cb_m1%q_mix(:)
dp_cld(:)      = inv_m1%p_rho(:) - cb_m1%p_rho(:)
dthetav_cld(:) = inv%theta(:)*                                                 &
   (1.0+inv%q_mix(:)/repsilon)/(1.0+inv%q_mix(:))                              &
   - cb_m1%theta(:)                                                            &
   *(1.0+cb_m1%q_mix(:)/repsilon)/(1.0+cb_m1%q_mix(:))

!-----------------------------------------------------------------
! calculate delta values across cloud base
!-----------------------------------------------------------------
dz_cb(:)      = cb%z_rho(:) - cb_m1%z_rho(:)
dthetav_cb(:) = cb%theta(:)*                                                   &
   (1.0+cb%q_mix(:)/repsilon)/(1.0+cb%q_mix(:))                                &
   - cb_m1%theta(:)*                                                           &
   (1.0+cb_m1%q_mix(:)/repsilon)/(1.0+cb_m1%q_mix(:))

!-----------------------------------------------------------------
! dthetav/dz along a moist adiabat at LCL
! ntml level below cloud base ?
!-----------------------------------------------------------------
do i=1,n_xx
  t_lcl_temp = cb_m1%theta(i)*cb_m1%exner_theta(i)
  temp1      = cb_m1%qse(i)*lc/(rv*t_lcl_temp*t_lcl_temp)
  temp2      = -gamma_dry*(1.0+t_lcl_temp*temp1*recip_epsilon)                 &
     /(1.0+lc_o_cp*temp1)
  dthetavdz_ad_lcl(i) = temp2*(1.0+0.61*cb_m1%theta(i)*temp1)+                 &
     gamma_dry +0.61*g*cb_m1%qse(i)*cb_m1%theta(i)                             &
     / (r*t_lcl_temp)
end do
if (L_tracer) then
  dtracer_cb(:,:) = cb%tracer(:,:) - cb_m1%tracer(:,:)
end if  ! L_tracer

!-----------------------------------------------------------------------
! 3.0 Calculate scalings
!-----------------------------------------------------------------------
! wstar_dn (from BL) sub cloud layer velocity scale
! scales%wstar_up  cumulus layer convective velocity scale
! scales%wstar_up = (scales%mb CAPE)**1/3
! scales%mb = 0.04*wstar_dn - cloud base mass flux
! zlcl_uv  - convection cloud base  from conv_diag
! ztop_uv  - base of inversion  from conv_diag
! scales%zcld - depth of congestus cloud
! delthuv  - from conv_diag some form of dz
!-----------------------------------------------------------------------
call calc_scales_warm(scales, conv_type, ntml, ntpar, wstar_dn, delthvu,       &
   theta, q_mix,                                                               &
   ztop_uv, zlcl_uv, z_rho)

!-----------------------------------------------------------------------
!
! 4.0 Calculate in-cloud height, corresponding compressed arrays and
!     similarity functions
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------
! Set in-cloud fields
!-----------------------------------------------------------------
!Working note: do we really need max_cldlev+1?  See other
! instances in the code.

call assign_cloud_input(cld_in, ntml+1, max_cldlev+1,                          &
   max_cldtrlev+1, theta,                                                      &
   q_mix, qse, p_layer_centres, p_layer_boundaries,                            &
   exner_layer_centres, exner_layer_boundaries,                                &
   z_theta, z_rho, r2rho, r2rho_theta, rho, tracer )

!-----------------------------------------------------------------
! Set derived in-cloud fields
!-----------------------------------------------------------------
do i=1,n_xx
  do k=1,max_cldlev+1
    h_cld(i,k)      = cld_in%theta(i,k) + lc_o_cp*cld_in%q_mix(i,k)
    thetav_cld(i,k) = cld_in%theta(i,k)                                        &
       *(1.0+cld_in%q_mix(i,k)/repsilon)                                       &
       /(1+cld_in%q_mix(i,k))
    t_cld(i,k)      = cld_in%theta(i,k)*cld_in%exner_theta(i,k)
    ztheta_cld(i,k) = cld_in%z_theta(i,k)
    r2rho_theta_cld(i,k) = cld_in%r2rho_theta(i,k)
    r2rho_cld(i,k)  = cld_in%r2rho(i,k)
  end do
end do

!-----------------------------------------------------------------
! Calculate LES similarity functions of non-dimensional height
!-----------------------------------------------------------------

call calc_similarity(cld_in%eta_theta, cld_in%eta_rho,                         &
   conv_type, sim )

!-----------------------------------------------------------------------
!
! 5.0 Calculation of cloud base fluxes, precipitation and inversion fluxes
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------
! Calculate cloud base fluxes
!-----------------------------------------------------------------
call calc_cloudbase( n_xx, ntra, l_tracer, wthvs, dtracer_cb,                  &
   dthetav_cb, dthetal_cb, dqt_cb,                                             &
   wtheta_plus_cb, wthetal_cb,wq_plus_cb,                                      &
   wqt_cb,wql_cb, wthetav_cb,wthetav_cb_m2,wtracer_cb)

!-----------------------------------------------------------------
! Calculate precip
!-----------------------------------------------------------------
! Different methods can be controlled by icong_options - for
! now there is just one method.
!    if (icong_options >= 1) then
dz_inv(:) = inv%z_theta(:) - inv_m1%z_theta(:)
call calc_precip( n_xx, nlev, conv_type,                                       &
   ql_ad, dz_inv, inv%rho,                                                     &
   cld_in%eta_rho, cld_in%eta_theta,                                           &
   rain, wqr_inv, precip_product_inv,                                          &
   wqr_cld, precip_product_cld,precip_product_hcld)

!    end if

    !-----------------------------------------------------------------
    ! Calculate inversion fluxes
    !-----------------------------------------------------------------
call calc_inversion_warm( n_xx, ntra, l_tracer,wthvs,                          &
   dthetav_cld, dp_cld,wqr_inv, precip_product_inv,                            &
   wtheta_inv,wq_inv,wql_inv,wthetavl_inv,wthetav_inv,                         &
   wthetal_inv,wqt_inv,wh_inv,wtracer_inv)

!-----------------------------------------------------------------------
! 6.0 Calculate ensemble vertical velocity and mass flux
!     On cld_in%eta_rho levels or full levels ?
!     At present doing both.
!-----------------------------------------------------------------------

! Working note: would like to put all of updraught calculations
! into a single plume subroutine (and possibly precip?)

do k=1,max_cldlev
  do i = 1,n_xx
    wup2 = scales%wup2_cb(i) +                                                 &
       scales%wstar_up(i)*scales%wstar_up(i)*sim%fw_func(i,k)
    wup2_h = scales%wup2_cb(i) +                                               &
       scales%wstar_up(i)*scales%wstar_up(i)*sim%fw_func_rho(i,k)

    ! Removed alternate options here
    mf_cld(i,k) =   sim%cmask(i,k) * scales%mb_new(i)                          &
                  * (  scales%wup2_cb(i)  + scales%wstar_up(i)                 &
                     * scales%wstar_up(i) * sim%g_func(i,k) )                  &
                  / wup2

    mf_h_cld(i,k) = sim%cmask_rho(i,k) * scales%mb_new(i)                      &
                  * (  scales%wup2_cb(i)  + scales%wstar_up(i)                 &
                     * scales%wstar_up(i) * sim%g_func_rho(i,k) )              &
                  / wup2_h

    wup_cld(i,k) = sqrt(sim%cmask(i,k)*wup2)
    wup_h_cld(i,k) = sqrt(sim%cmask_rho(i,k)*wup2_h)
  end do
end do

!-----------------------------------------------------------------------
! 7.0 Calculate qlup (ie CCW) (all convective cloud levels)
!-----------------------------------------------------------------------
! This routine is purely diagnostic, calculating qlup (CCW)

call calc_qlup_warm( n_xx, max_cldlev, wthetal_cb, wqt_cb,                     &
   dthetal_cb, dqt_cb, dthetal_cld, dqt_cld, cld_in, sim,                      &
   ql_up_cld)

!-----------------------------------------------------------------------
! 8.0 Calculate liquid water flux and dqsat/dT in cloud
!-----------------------------------------------------------------------
call calc_wql_warm(n_xx, nlev, max_cldlev, ncld_thlev, wql_cb,                 &
   wql_inv, cld_in, t_cld, sim, wql_cld, dqsatdt_cld)

!----------------------------------------------------------------
! add on wqr to wql calculated so far to take account of precip
!----------------------------------------------------------------

do k=1,max_cldlev+1
  do i=1,n_xx
    wql_cld(i,k) = wql_cld(i,k) + wqr_cld(i,k)
  end do
end do

!-----------------------------------------------------------------------
! 9.0 turbulent transports  - implicit Fluxes on half levels
!                             therefore need cloud base and cloud top
!                             values.
!
!  Requires mass_flux, w_up, sim%k_func, sim%fng_func, sim%b_func
!           wthetal_cb, wqt_cb
!
! Note
!  w'h' = cp* w'theta' +(L/exner)*w'q'
!       = cp*w'thetal' + (L/exner)*w'qt'
!-----------------------------------------------------------------------
! Problems here because of levels required for this calculation

! 9.1 Full K and B functions on cld_in%eta_rho

do k =1,max_cldlev
  do i=1,n_xx
    sim%k_func_rho(i,k) = mf_h_cld(i,k) *                                      &
       (wup_h_cld(i,k)/scales%wstar_up(i)) *                                   &
       scales%zcld(i)*sim%k_func_rho(i,k)
  end do
end do

! 9.2 Evaluation of thetav at t+1

call calc_grad_h(n_xx, max_cldlev, nlev                                        &
   ,                      timestep, ncld_thlev                                 &
   ,                      cld_in%z_rho, ztheta_cld, r2rho_cld                  &
   ,                      r2rho_theta_cld                                      &
   ,                      thetav_cld, sim%k_func_rho                           &
   ,                      wthetavl_inv                                         &
   ,                      Kdthetavdz )

! 9.2 Evaluation of moist static energy at t+1

call calc_grad_h(n_xx, max_cldlev, nlev                                        &
   ,                      timestep, ncld_thlev                                 &
   ,                      cld_in%z_rho, ztheta_cld, r2rho_cld                  &
   ,                      r2rho_theta_cld                                      &
   ,                      h_cld, sim%k_func_rho                                &
   ,                      wh_inv                                               &
   ,                      Kdhdz )

! Tracer at t+1
! loop over tracers calculating dtracer/dz

if (L_tracer) then
  do ktra = 1,ntra

    ! Note kdtracerdz dimensioned on nlev so subroutine works
    call calc_grad_h(n_xx, max_cldtrlev, nlev                                  &
       ,                      timestep, ncld_thlev                             &
       ,                      cld_in%z_rho, ztheta_cld, r2rho_cld              &
       ,                      r2rho_theta_cld                                  &
       ,                      cld_in%tracer(:,:,ktra), sim%k_func_rho          &
       ,                      wtracer_inv(:,ktra)                              &
       ,                      kdtracerdz(:,:,ktra) )
  end do

end if  ! L_tracer
!--------------------------------------------------------------------
! 9.3 Evaluation wtheta and wq, wthetavl, wqt, wthetal fluxes in cloud

call calc_turb_fluxes(n_xx, ntra, max_cldlev, nlev                             &
   ,                      max_cldtrlev,trlev                                   &
   ,                      l_tracer                                             &
   ,                      Kdthetavdz, kdtracerdz                               &
   ,                      sim                                                  &
   ,                      wql_cld, mf_cld, cld_in                              &
   ,                      dqsatdt_cld                                          &
   ,                      precip_product_hcld,wup_h_cld                        &
   ,                      scales, wthetav_cb_m2, wtracer_cb                    &
   ,               wtheta_cld, wq_cld , wthetavl_cld, wthetav_cld              &
   ,               wthetal_cld,wqt_cld, wh_cld, wtracer_cld)



!--------------------------------------------------------------------
! 9.4 store various fluxes back on full model levels for diagnostic
!     output and for calculating gradients

! Cloud base and inversion fluxes

do i = 1,n_xx
  ibase              = ntml(i)+1
  kinv               = ntpar(i) + 1
  wq(i,ibase)        = wq_plus_cb(i)
  wq(i,kinv)         = wq_inv(i)
  wqt(i,ibase)       = wqt_cb(i)
  wqt(i,kinv)        = wqt_inv(i)
  wtheta(i,ibase)    = wtheta_plus_cb(i)
  wtheta(i,kinv)     = wtheta_inv(i)
  wthetavl(i,ibase)  = wthetav_cb(i)
  wthetavl(i,kinv)   = wthetavl_inv(i)
  wthetav(i,ibase)   = wthetav_cb(i)
  wthetav(i,kinv)    = wthetav_inv(i)
  wthetal(i,ibase)   = wthetal_cb(i)
  wthetal(i,kinv)    = wthetal_inv(i)
  wql(i,kinv)        = wql_inv(i)
  !
  ! Leave flux of ql as zero through cloud base for now
  !
  !      wql(i,ibase)       = wql_cb(i)
  wqr(i,kinv)        = wqr_inv(i)
end do

if (L_tracer) then
  do ktra = 1,ntra
    do i = 1,n_xx
      ibase                 = ntml(i)+1
      kinv                  = ntpar(i) + 1
      wtracer(i,ibase,ktra) = wtracer_cb(i,ktra)
      wtracer(i,kinv,ktra)  = wtracer_inv(i,ktra)
    end do
  end do
end if

!
! Expand in cloud values on uv levels
!
! Also added in smoothing of/to inversion fluxes
!
do k = 2,max_cldlev
  do i = 1,n_xx
    kinv = ntpar(i)+1
    ilev = ntml(i)+k
    if (ilev <  kinv) then

      coeff_a =  bsmooth_inv*cld_in%eta_theta(i,k)**asmooth_inv
      coeff_b =  1.0 - coeff_a

      wq(i,ilev)       = coeff_a * wq(i,kinv)       + coeff_b * wq_cld(i,k)
      wqt(i,ilev)      = coeff_a * wqt(i,kinv)      + coeff_b * wqt_cld(i,k)
      wtheta(i,ilev)   = coeff_a * wtheta(i,kinv)   + coeff_b * wtheta_cld(i,k)
      wthetal(i,ilev)  = coeff_a * wthetal(i,kinv)  + coeff_b * wthetal_cld(i,k)
      wthetavl(i,ilev) = coeff_a * wthetavl(i,kinv) + coeff_b * wthetavl_cld(i,k)
      wthetav(i,ilev)  = coeff_a * wthetav(i,kinv)  + coeff_b * wthetav_cld(i,k)
      wh(i,ilev)       = coeff_a * wh(i,kinv)       + coeff_b * wh_cld(i,k)
      wql(i,ilev)      = coeff_a * wql(i,kinv)      + coeff_b * wql_cld(i,k)
      wqr(i,ilev)      = coeff_a * wqr(i,kinv)      + coeff_b * wqr_cld(i,k)

    end if
  end do
end do

if (L_tracer) then
  do ktra = 1,ntra
    do k = 2,max_cldtrlev
      do i = 1,n_xx
        kinv = ntpar(i)+1
        ilev = ntml(i)+k
        if (ilev <  kinv) then
          wtracer(i,ilev,ktra) = wtracer_cld(i,k,ktra)
        end if
      end do
    end do
  end do
end if
!
! theta levels
!
do k = 1,max_cldlev
  do i = 1,n_xx
    kinv = ntpar(i)+1
    ilev = ntml(i)+k
    if (ilev <= kinv) then
      wup(i,ilev)            = wup_cld(i,k)
      ccw(i,ilev)            = ql_up_cld(i,k)
      up_flux(i,ilev)        = mf_cld(i,k)
      precip_product(i,ilev) = precip_product_cld(i,k)/rho_theta(i,ilev)
    end if
  end do

end do

! Fluxes below cloud base
! Know fluxes at cloud base for w'q' and w'theta' assume they
! go linearly to zero at surface below cloud base
! At bottom of cloud base transition region
!  w'q' = w'qt'_cb    & w'theta' = w'thetal'_cb

do k=1,ntml_max
  do i=1,n_xx
    if (k <= ntml(i)) then

      factor       = z_rho(i,k)/z_rho(i,ntml(i)+1)
      wqt(i,k)     = factor*wqt_cb(i)
      wthetal(i,k) = factor*wthetal_cb(i)

    end if
  end do
end do

if (L_tracer) then
  do ktra = 1,ntra
    do k = 1,ntml_max
      do i = 1,n_xx
        if (k <= ntml(i)) then
          factor            = z_rho(i,k)/z_rho(i,ntml(i)+1)
          wtracer(i,k,ktra) = factor*wtracer_cb(i,ktra)
        end if
      end do
    end do
  end do

end if

!--------------------------------------------------------------------
! 9.5 calculate gradient of fluxes in cloud
call calc_grad_flux (n_xx,ntra, nlev,trlev,ntpar_max,l_tracer                  &
   ,                  r2rho,r2rho_theta,dr_across_th1                          &
   ,                  wthetavl, wthetal, wqt, wtracer                          &
   ,                  dwthetavl_dz, dwthetal_dz, dwqt_dz, dwtracer_dz)
!
! Increments to model
! method for calculating increments to theta and rv
!
do k = 1,nlev
  do i = 1,n_xx
    dthbydt(i,k) = -dwthetal_dz(i,k)                                           &
       + lc_o_cp*precip_product(i,k)/exner_layer_centres(i,k)
    dqbydt(i,k)  = -dwqt_dz(i,k) - precip_product(i,k)
  end do
end do

! tracer increments

if (L_tracer) then
  do ktra = 1,ntra
    do k = 1,trlev
      do i = 1,n_xx
        dtrabydt(i,k,ktra) = -dwtracer_dz(i,k,ktra)
      end do
    end do
  end do

end if

!
! Convective cloud liquid water path
!
do i = 1,n_xx
  cclwp(i) = 0.0
end do

do k=1,nlev-1
  do i=1,n_xx

    dmass    = r2rho_theta(i,k)*dr_across_th(i,k)*ra2
    cclwp(i) = cclwp(i) + ccw(i,k)*dmass

    ! Is this still correct given up_flux nolonger upward ?

    tcw(i) = tcw(i) - ccw(i,k)*up_flux(i,k)*dmass

  end do
end do

!-----------------------------------------------------------------------
! 10.0 Convective Momentum Transport as existing scheme (if L_mom = .T.)
!-----------------------------------------------------------------------
! But recoded to use wup and mass_flux from above plus height directly
!-----------------------------------------------------------------------
!       Because of the definition of nlcl, the pressure of the top of
!       the mixed layer is phalf_uv(nlcl,*)
!
! cld_in%eta_theta levels for this calculation different from above as cloud
!  base and top taken to be theta levels for this calculation
!-----------------------------------------------------------------------
! Initialize arrays required for Convective Momentum Transport(CMT)
!-----------------------------------------------------------------------

if (L_mom) then

  do i = 1,n_xx
    nlcl_uv(i)    = ntml(i) + 1
    ntop_uv(i)    = ntpar(i) + 1
  end do

  !-----------------------------------------------------------------------
  ! 10.1 compression to convective cloud levels - U & V
  !-----------------------------------------------------------------------
  ! Note fill arrays even if above cloud top

  do k=1,max_cldlev+1
    do i = 1,n_xx
      klev = ntml(i)+k
      u_cld(i,k)            = u(i,klev)
      v_cld(i,k)            = v(i,klev)
      dr_across_rh_cld(i,k) = dr_across_rh(i,klev)
      dr_across_th_cld(i,k) = dr_across_th(i,klev)
      r2rho_cld(i,k)        = r2rho(i,klev)
      r2rho_theta_cld(i,k)  = r2rho_theta(i,klev)
    end do
  end do

  !-----------------------------------------------------------------------
  ! 10.2  Convective Momentum Transport (if L_mom = .T.)
  !-----------------------------------------------------------------------

  ! Altered value of scales%wsc_o_mb ? here as well need original ?
  ! what values of mf_cld etc do I pass in ? those
  ! on theta, uv levels relative to zcld_th or zcld_cmt?
  ! need to try different ones to see which work best.

  call calc_grad_stress(n_xx, nlev,max_cldlev                                  &
     ,                        ncld_thlev                                       &
     ,                        timestep,scales                                  &
     ,                        mf_cld,wup_cld,sim                               &
     ,                        u_cld,v_cld                                      &
     ,    r2rho_cld,r2rho_theta_cld,dr_across_rh_cld,dr_across_th_cld          &
                            ! out
     ,                        uw_cld,vw_cld)


  !
  ! Expand in cloud values on uv levels
  !
  do k = 1,max_cldlev
    do i = 1,n_xx
      kinv = ntpar(i)+1
      ilev = ntml(i)+k
      if (ilev <  kinv) then
        uw(i,ilev)      = uw_cld(i,k)
        vw(i,ilev)      = vw_cld(i,k)
      end if
    end do
  end do

  !
  !  change in winds across cloud base
  !
  do i = 1,n_xx
    du_cb(i) = u(i,ntml(i)+1) - u(i,ntml(i))
    dv_cb(i) = v(i,ntml(i)+1) - v(i,ntml(i))
  end do

  call calc_base_stress(n_xx,nlev,ntml,ntpar,ntpar_max                         &
     ,                            timestep, scales                             &
     ,                            uw0,vw0, du_cb,dv_cb                         &
     ,                            rho_theta,z_rho,z_theta                      &
     ,                            flg_uw_shall,flg_vw_shall                    &
                            ! INOUT
     ,                            uw,vw                                        &
                            ! out
     ,                            uw_shall,vw_shall)


  call calc_cmt_incr(n_xx,nlev, ntpar_max                                      &
     ,                          r2rho, r2rho_theta                             &
     ,                          dr_across_rh                                   &
     ,                          uw ,vw                                         &
                            !out
     ,                          dubydt,dvbydt)


end if ! L_mom


!-----------------------------------------------------------------------
! 11.0 Calculate increments for PC2
!-----------------------------------------------------------------------

if (L_calc_dxek) then
  call calc_pc2( n_xx, max_cldlev, timestep, ntml, ntpar,                      &
 mf_h_cld, ql_up_cld, wql_cld, qcl, cf_liquid,                                 &
 cld_in, sim, dqclbydt, dqcfbydt, dbcfbydt, dcflbydt, dcffbydt)
end if

!-----------------------------------------------------------------------
! 12.0  Check for negative or very small water vapour after convection
!-----------------------------------------------------------------------
! Trys to prevent negative q by adjusting the vertical integral of dq/dt
! in a conservative way.
!-----------------------------------------------------------------------
! Find minimum value of q in column

do i = 1,n_xx
  qMinInColumn(i) = q_mix(i,nlev)
end do
do k = 1,nlev-1
  do i = 1,n_xx
    if (q_mix(i,k)  <   qMinInColumn(i)) then
      qMinInColumn(i) = q_mix(i,k)
    end if
  end do
end do

!
! Ensure Q does not go below global allowed minimum (QMIN)
!
do i = 1,n_xx
  qMinInColumn(i)=max(qmin,qMinInColumn(i))
end do

!
! Apply an artificial upwards flux from k-1 level to ensure Q
! remains above minimum value in the column.
! Only safe if increment below is positive (definitely a problem if do
! below ntml  where all increments are negative).
!

do k = nlev,2,-1
  do i = 1,n_xx

    tempx(i)=q_mix(i,k) + dqbydt(i,k) * timestep

    if (tempx(i)  <   qMinInColumn(i)) then

      ! safe to correct level below has a positive increment > negative
      if (dqbydt(i,k-1) >= 0.0                                                 &
         .and. dqbydt(i,k-1) >  (-1.0*dqbydt(i,k))) then
        dqbydt(i,k-1) = dqbydt(i,k-1) -                                        &
           ((qMinInColumn(i) - q_mix(i,k)) / timestep-dqbydt(i,k))             &
           * (r2rho_theta(i,k)*dr_across_th(i,k))                              &
           / (r2rho_theta(i,k-1)*dr_across_th(i,k-1))

        dqbydt(i,k) = (qMinInColumn(i) - q_mix(i,k)) / timestep

      else
        ! unsafe to correct ?
        ! try checking on level 2 below
        if (dqbydt(i,k-2) >= 0.0) then

          dqbydt(i,k-1) = dqbydt(i,k-1) -                                      &
             ((qMinInColumn(i) - q_mix(i,k)) / timestep-dqbydt(i,k))           &
             * (r2rho_theta(i,k)*dr_across_th(i,k))                            &
             / (r2rho_theta(i,k-1)*dr_across_th(i,k-1))

          dqbydt(i,k) = (qMinInColumn(i) - q_mix(i,k)) / timestep

        else
          write(umMessage,*) 'TCS Warm problem negative q ',i,k,               &
             tempx(i),dqbydt(i,k), q_mix(i,k),                                 &
             (dqbydt(i,j),j=1,ntpar(i)+1)
          call umPrint(umMessage,src='tcs_warm_mod')

        end if
      end if

    end if   ! resultant less than qmin
  end do     ! n_xx loop
end do       ! nlev loop

! check q on bottom level

k=1
do i = 1,n_xx

  tempx(i)=q_mix(i,k) + dqbydt(i,k) * timestep

  if (tempx(i)  <   qMinInColumn(i)) then
    write(umMessage,*) 'TCS Warm Problem negative q k=1 ',i,                   &
       tempx(i),dqbydt(i,k),q_mix(i,k),                                        &
       (dqbydt(i,j),j=1,ntpar(i)+1)
    call umPrint(umMessage,src='tcs_warm_mod')

  end if
end do ! n_xx loop


!-----------------------------------------------------------------------
! 13.0 Tidy up and Deallocate local storage
!-----------------------------------------------------------------------

!----------------------------------------
! Some diagnostic output
! Congestus convective cloud amounts etc
! Taken from mods for HadGEM1
!----------------------------------------
iccb(:)   = ntml(:)+1
lcbase(:) = iccb(:)
icct(:)   = ntpar(:)+1
lctop(:)  = icct(:)
cca_2d(:) = 2.0*scales%mb_o_wsc(:)
lcca(:)   = cca_2d(:)
! Working note: Transfer some output variables (This is messy, but
! is necessary to maintain the interface with the calling procedure)
wstar_up(:) = scales%wstar_up(:)
mb1(:)      = scales%mb(:)
mb2(:)      = scales%mb_new(:)

!-----------------------------------------------------------------------
! Convective cloud amount 3d - no anvil
!-----------------------------------------------------------------------
! Initialise output array

do k = 1,n_cca_lev
  do i = 1,n_xx
    cca(i,k) = 0.0
  end do
end do

! Whether PC2 or not set 3D cca values for convective cloud
! The new shallow turbulence scheme provides diagnostic CCA and CCW
! to PC2

do k = 1,n_cca_lev
  do i = 1,n_xx
    if (k >= iccb(i) .and. k <= icct(i)) then
      cca(i,k) = cca_2d(i)
    end if
  end do
end do

! Deallocation/disassociation
call tcs_deallocate(scales) ! scales_conv
call tcs_deallocate(cld_in) ! cloud_input
call tcs_deallocate(sim)    ! similarity

call tcs_deallocate(inv_m1) ! interface_input
call tcs_deallocate(inv_p1) ! interface_input
call tcs_deallocate(inv)    ! interface_input

call tcs_deallocate(cb_m1)  ! interface_input
call tcs_deallocate(cb_p1)  ! interface_input
call tcs_deallocate(cb)     ! interface_input

!     ErrorStatus=1
!     Message='tcs_warm scheme not yet fully functional.'
!
!     call Ereport(routinename, ErrorStatus, Message)



if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine tcs_warm

end module tcs_warm_mod
