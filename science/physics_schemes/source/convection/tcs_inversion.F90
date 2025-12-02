! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calculate values of theta, q etc in inversion. These
!   values are then used to calculate the fluxes at the inversion.
!
module tcs_inversion


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
!   Module to calculate values of theta, q etc inversion. These
!   values are then used to calculate the fluxes at the inversion.
!
! Method:
!   <reference to documentation to go here, once available>
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!

character(len=*), parameter, private :: ModuleName='TCS_INVERSION'

contains

subroutine calc_inversion_warm( n_xx, ntra, l_tracer,                          &
   wthetav_surf, dthetav_cld, dp_cld, wqr_inv,                                 &
   precip_product_inv_base, wtheta_inv_base,wq_inv,wql_inv,                    &
   wthetavl_inv, wthetav_inv,wthetal_inv, wqt_inv_base,wh_inv,                 &
   wtracer_inv)

  !-------------------------------------------------------------------
  !
  ! Description:
  !   Calculate values of theta, q etc at the inversion level. These
  !   values are then used to calculate the fluxes at the inversion.
  !   The inversion base level is defined as a uv level in the scheme.
  !
  !-------------------------------------------------------------------

use tcs_parameters_warm,  only:                                                &
   beta_cld, jump_inv_a1, jump_inv_a2, jump_inv_a3                             &
   , mb_inv_a1 ,sat_inv_a1, q_min
use tcs_constants,        only:                                                &
   repsilon, lc_o_cp, a_bolton, b_bolton, c_bolton                             &
   ,d_bolton, pref, c_virtual, lc, rv, kappa, g                                &
   ,recip_kappa
use tcs_common_warm ,     only:                                                &
   scales, inv_m1, inv, inv_p1

use qsat_mod, only: qsat_mix

implicit none
!------------------------------------------------------------------
! Subroutine Arguments
!------------------------------------------------------------------
integer, intent(in) ::                                                         &
   n_xx                                                                        &
                            ! No. of congestus convection points
   ,  ntra
! Number of tracer fields

logical, intent(in) ::                                                         &
   L_tracer         ! true - tracers present.

real(kind=real_umphys), intent(in) ::                                          &
   wthetav_surf(n_xx)                                                          &
                            ! wthetav flux at the surface
   , dthetav_cld(n_xx)                                                         &
                            !  dthetav across cloud
   , dp_cld(n_xx)                                                              &
                            !  dp across cloud
   , wqr_inv(n_xx)                                                             &
                            ! rain water flux (kg/kg m/s)
   , precip_product_inv_base(n_xx)
! precip production from inversion
! (kg/m2/s) / density ?

real(kind=real_umphys), intent(out) ::                                         &
   wthetal_inv(n_xx)                                                           &
                            !  w'thetal' at inversion (Km/s)
   , wtheta_inv_base(n_xx)                                                     &
                            !  w'theta' at inversion (Km/s)
   , wq_inv(n_xx)                                                              &
                            !  w'q'  at inversion (kg/kg m/s)
   , wqt_inv_base(n_xx)                                                        &
                            !  w'qt'  at inversion (kg/kg m/s)
   , wql_inv(n_xx)                                                             &
                            !  w'ql at inversion (kg/kg m/s)
   , wthetavl_inv(n_xx)                                                        &
                            !  w'thetavl' at inversion (Km/s)
   , wthetav_inv(n_xx)                                                         &
                            !  w'thetav' at inversion (Km/s)
   , wh_inv(n_xx)                                                              &
                            !  w'h' at inversion (Km/s)
   , wtracer_inv(n_xx,ntra)
                            !  w'tracer' at inversion (kgm/kg/s)

!------------------------------------------------------------------
! Variables defined locally
!------------------------------------------------------------------

real(kind=real_umphys) ::                                                      &
   dthetav_inv(n_xx)                                                           &
                            ! thetav across inversion (k)
   , dh_inv(n_xx)                                                              &
                            ! h across inversion (k)
   , thetav_mid(n_xx)                                                          &
   , thetav_below(n_xx)                                                        &
   , thetav_above(n_xx)                                                        &
   , h_mid(n_xx)                                                               &
   , h_below(n_xx)                                                             &
   , t_mid(n_xx)                                                               &
   , t_above(n_xx)                                                             &
   , t_below(n_xx)                                                             &
   , r_inv_base(n_xx)                                                          &
                            ! water vapour at inversion (mixing
                            !  ratio, kg/kg)
   , t_inv_base(n_xx)                                                          &
                            ! T at inversion (K)
   , rsat_inv_base(n_xx)                                                       &
                            ! rsat at inversion (kg/kg)
   , drsatdt_inv_base(n_xx)                                                    &
                            ! drsat/dT at inversion (kg/kg/K)
   , t_lcl_below(n_xx)                                                         &
                            ! T at LCL for level below inversion base
   , t_lcl_mid(n_xx)                                                           &
                            ! T at LCL for level above inversion base
   , t_lcl_above(n_xx)                                                         &
                            ! T at LCL for level above inversion top
   , t_lcl_inv(n_xx)                                                           &
                            ! T at LCL for  inversion base
   , pstar_below(n_xx)                                                         &
                            ! p at LCL for level below inversion base
   , pstar_mid(n_xx)                                                           &
                            ! p at LCL for level above inversion base
   , pstar_above(n_xx)                                                         &
                            ! p at LCL for level above inversion top
   , dthetav_inv2(n_xx)                                                        &
                            !  dthetav across inversion
   , dthetav_invnew(n_xx)                                                      &
                            !  dthetav across inversion
   , r_rsat_inv_base(n_xx)                                                     &
                            ! r-rsat at inversion
   , dpstar_inv(n_xx)                                                          &
                            ! dpstar across inversion
   , dp_inv(n_xx)                                                              &
                            ! dp across inversion
   , rsat_lcl(n_xx)                                                            &
                            ! rsat at lcl
   , dthetav_dpstar_top(n_xx)                                                  &
                            ! dthetav/dp* at top of inversion
   , dthetav_dpstar_base(n_xx)                                                 &
                            ! dthetav/dp* at base of inversion
   , dthetav_dpstar_layer(n_xx)                                                &
                            ! dthetav/dp* across inversion
   , beta_inv(n_xx)                                                            &
                            ! beta for the inversion layer
   , dp_lower(n_xx)                                                            &
                            ! dp across lower part of inversion
   , dp_upper(n_xx)                                                            &
                            ! dp across upper part of inversion
   , thetav_inv_base(n_xx)                                                     &
                            ! thetav and base of inversion (k)
   , thetav_inv_top(n_xx)                                                      &
                            ! thetav and top of inversion (k)
   , theta_inv_base(n_xx)                                                      &
                            ! thetav and base of inversion (k)
   , pstar_inv_base(n_xx)
                            ! pstar at inversion base

! temporary variables

real(kind=real_umphys) ::                                                      &
   term_a                                                                      &
   , term_b                                                                    &
   , term_c                                                                    &
   , vap_press                                                                 &
   , div_term2, divisor_term


logical :: l_3lev_inv     ! true for 3 level inversion cal
! current default is false.

!-------------------------
! Loop counters
!-------------------------
integer :: i, ktra

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_INVERSION_WARM'


!----------------------------------------------------------------
! 1.0 Initialise arrays
!----------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

l_3lev_inv = .false. ! Flag to use a 3 level inversion scheme
! Working note: this doesn't yet work, but
! if/when it does, then this switch should
! be put at a higher level

wtheta_inv_base(:) = 0.0
wthetavl_inv(:)    = 0.0
wthetav_inv(:)     = 0.0
wh_inv(:)          = 0.0
wthetal_inv(:)     = 0.0
wq_inv(:)          = 0.0
wqt_inv_base(:)    = 0.0
wql_inv(:)         = 0.0
wtracer_inv(:,:)   = 0.0

!----------------------------------------------------------------
! derived values at levels above and below inversion base and top
!----------------------------------------------------------------
thetav_below(:) = inv_m1%theta(:)*                                             &
   (1.0+inv_m1%q_mix(:)/repsilon)/(1.0+inv_m1%q_mix(:))
thetav_mid(:)   = inv%theta(:)*                                                &
   (1.0+inv%q_mix(:)/repsilon)/(1.0+inv%q_mix(:))
thetav_above(:) = inv_p1%theta(:)*                                             &
   (1.0+inv_p1%q_mix(:)/repsilon)/(1.0+inv_p1%q_mix(:))

h_below(:) = inv_m1%theta(:) + lc_o_cp*inv_m1%q_mix(:)
h_mid(:)   = inv%theta(:)   + lc_o_cp*inv%q_mix(:)

t_below(:) = inv_m1%theta(:)*inv_m1%exner_theta(:)
t_mid(:)   = inv%theta(:)  *inv%exner_theta(:)
t_above(:) = inv_p1%theta(:)*inv_p1%exner_theta(:)


!----------------------------------------------------------------
! Calculate of r and theta at inversion
! Require dp*/dp to do this  dtheta/dp* and dr/dp* are ~ constant
! through inversion
!
! Need p* above and below inversion base and top
!
! Calculate temperature and pressure of lifting condensation level
!  using approximations from Bolton (1980)
!----------------------------------------------------------------

if (l_3lev_inv) then
  ! Working note: This 3 level inversion method does not
  ! currently work.
  do i = 1,n_xx
    if (inv%q_mix(i) >  0.0) then
      vap_press = 0.01*inv%q_mix(i)*inv%p_theta(i)/                            &
                      (repsilon+inv%q_mix(i))
    else
      vap_press = 0.01*q_min*inv%p_theta(i)/(repsilon+q_min)
    end if
    t_lcl_mid(i) =  a_bolton + b_bolton / ( c_bolton*log(t_mid(i))             &
       - log(vap_press) - d_bolton )
    pstar_mid(i) = inv%p_theta(i) *                                            &
       ( t_lcl_mid(i) / t_mid(i) )**recip_kappa

    if (inv_m1%q_mix(i) >  0.0) then
      vap_press = 0.01*inv_m1%q_mix(i)*inv_m1%p_theta(i)                       &
         /(repsilon+inv_m1%q_mix(i))
    else
      vap_press = 0.01*q_min*inv_m1%p_theta(i)/(repsilon+q_min)
    end if
    t_lcl_below(i) = a_bolton + b_bolton/( c_bolton*log(t_below(i))            &
       - log(vap_press) - d_bolton )
    pstar_below(i) = inv_m1%p_theta(i) *                                       &
       ( t_lcl_below(i) / t_below(i) )**recip_kappa

    if (inv_p1%q_mix(i) >  0.0) then
      vap_press = 0.01*inv_p1%q_mix(i)*inv_p1%p_theta(i)                       &
         /(repsilon+inv_p1%q_mix(i))
    else
      vap_press = 0.01*q_min*inv_p1%p_theta(i)/(repsilon+q_min)
    end if
    t_lcl_above(i) = a_bolton + b_bolton/( c_bolton*log(t_above(i))            &
       - log(vap_press) - d_bolton )
    pstar_above(i) = inv_p1%p_theta(i) *                                       &
       ( t_lcl_above(i) / t_above(i) )**recip_kappa
  end do

  ! gradients of dthetav/dp* across inversion base and inversion top
  !  and across the whole layer

  do i = 1,n_xx
    dthetav_dpstar_top(i) = (thetav_above(i)-thetav_mid(i))/                   &
       (pstar_above(i)-pstar_mid(i))

    dthetav_dpstar_base(i) = (thetav_mid(i)-thetav_below(i))/                  &
       (pstar_mid(i)-pstar_below(i))

    dthetav_dpstar_layer(i) = (thetav_above(i)-thetav_below(i))/               &
       (pstar_above(i)-pstar_below(i))

    ! pressure across inversion layer
    dp_inv(i) = inv%p_rho(i)-inv_m1%p_rho(i)


    ! dthetav across the inversion layer

    dthetav_inv(i) = (0.16+0.55*scales%zcld(i)/scales%zlcl(i))                 &
       *wthetav_surf(i)/(0.2*scales%mb(i))


    ! dp* across inversion layer

    dpstar_inv(i) = dthetav_inv(i)/dthetav_dpstar_layer(i)

    beta_inv(i)  = dpstar_inv(i)/dp_inv(i)

    dp_lower(i) = inv%p_theta(i) - inv_m1%p_rho(i)
    dp_upper(i) = inv%p_rho(i) - inv%p_theta(i)
  end do

  ! calculate thetav at inversion top and inversion base using
  ! interpolation from value at inversion middle.

  do i = 1,n_xx

    thetav_inv_base(i) = thetav_mid(i)                                         &
       - beta_inv(i)*dthetav_dpstar_base(i)*dp_lower(i)

    thetav_inv_top(i) = thetav_mid(i)                                          &
       + beta_inv(i)*dthetav_dpstar_top(i)*dp_upper(i)

    pstar_inv_base(i) = pstar_mid(i) - beta_inv(i)*dp_lower(i)

    dthetav_inv2(i) = thetav_inv_top(i) - thetav_inv_base(i)

    ! There are points with problems
    !  dthetav_inv2   -ve
    !  dthetav_inv2 >> dthetav_inv
    ! For these points need to use another method ?


  end do

else   ! 2 level inversion

  do i = 1,n_xx

    if (inv%q_mix(i) >  0.0) then
      vap_press = 0.01*inv%q_mix(i)*inv%p_theta(i)/                            &
                            (repsilon+inv%q_mix(i))
    else
      vap_press = 0.01*q_min*inv%p_theta(i)/(repsilon+q_min)
    end if
    t_lcl_mid(i) = a_bolton + b_bolton /( c_bolton*log(t_mid(i))               &
       - log(vap_press) - d_bolton )
    pstar_mid(i) = inv%p_theta(i) *                                            &
       ( t_lcl_mid(i) / t_mid(i) )**recip_kappa

  end do

  do i = 1,n_xx
    dp_inv(i) = inv%p_theta(i)-inv_m1%p_theta(i)
    dthetav_invnew(i) = (0.16+0.55*scales%zcld(i)/scales%zlcl(i))              &
       *wthetav_surf(i)/(0.2*scales%mb(i))

    ! restriction on dthetav_cld to avoid problems with strange
    ! profiles

    if (dthetav_cld(i) >  0.01) then
      dpstar_inv(i)     = beta_cld*dthetav_invnew(i)*dp_cld(i)                 &
         /dthetav_cld(i)

    else
      dpstar_inv(i)     = beta_cld*dthetav_invnew(i)*dp_cld(i)                 &
         /0.01
    end if

    pstar_inv_base(i)  = pstar_mid(i) - dpstar_inv(i)
    thetav_inv_base(i) = thetav_mid(i) - dthetav_invnew(i)
    !
    ! value of thetav across inversion
    !
    dthetav_inv2(i) = thetav_mid(i) - thetav_below(i)
    dh_inv(i) = h_mid(i) - h_below(i)

  end do
end if  ! inversion cal type

!---------------------------------------------------------------------
! Calculate of r-rsat at inversion
!---------------------------------------------------------------------

do i = 1,n_xx
  t_lcl_inv(i) = thetav_inv_base(i)*                                           &
     ((pstar_inv_base(i)/pref)**kappa)
end do

call qsat_mix(rsat_lcl,t_lcl_inv,pstar_inv_base,n_xx)

do i = 1,n_xx
  r_inv_base(i) = rsat_lcl(i)/                                                 &
     (1.0+c_virtual*(lc/(rv*t_lcl_inv(i)))*rsat_lcl(i))

  theta_inv_base(i) = thetav_inv_base(i)                                       &
     /(1.0+c_virtual*r_inv_base(i))
  t_inv_base(i) = theta_inv_base(i)*inv_m1%exner_rho(i)
end do

call qsat_mix(rsat_inv_base,t_inv_base,inv_m1%p_rho,n_xx)

do i = 1,n_xx

  drsatdt_inv_base(i) = Lc *rsat_inv_base(i)/                                  &
     (Rv * t_inv_base(i)*t_inv_base(i))
  r_rsat_inv_base(i) = r_inv_base(i) - rsat_inv_base(i)

end do


!-----------------------------------------------------------------------
!  Calculate fluxes at inversion uses values calculated above
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
! Solution for fluxes at inversion
!----------------------------------------------------------------------


do i = 1,n_xx

  ! precipitation adds extra terms to A and B

  term_a = -jump_inv_a1*scales%mb_new(i)*dthetav_inv2(i) +                     &
     (lc_o_cp/inv_m1%exner_rho(i)                                              &
     -0.61*theta_inv_base(i))*precip_product_inv_base(i)

  term_b = thetav_inv_base(i)*scales%root_mb_new_o_wsc(i)                      &
     *(scales%wstar_up(i)**3)/                                                 &
     (g*scales%zcld(i))                                                        &
     +theta_inv_base(i)*wqr_inv(i)

  term_c = mb_inv_a1*scales%mb_new(i)*r_rsat_inv_base(i)

  ! Corrected version
  div_term2    = lc_o_cp/inv_m1%exner_rho(i)                                   &
     -(1.0+c_virtual)*theta_inv_base(i)
  divisor_term = (1.0+c_virtual*t_inv_base(i)*drsatdt_inv_base(i))

  !
  ! viscosity based value
  !
  wthetavl_inv(i) = term_a
  !
  ! Liquid water flux correct for mixing ratio
  !
  wql_inv(i) =(term_b-term_a)/div_term2

  wqt_inv_base(i) = ((1.0+lc_o_cp*drsatdt_inv_base(i))*                        &
     wql_inv(i)/sat_inv_a1-term_c+                                             &
     term_a*inv_m1%exner_rho(i)*drsatdt_inv_base(i))/divisor_term




  wthetal_inv(i) = term_a                                                      &
     -c_virtual*theta_inv_base(i)*wqt_inv_base(i)

  !
  ! above implies
  !
  wq_inv(i)=wqt_inv_base(i) - wql_inv(i)

  wtheta_inv_base(i) = wthetal_inv(i)                                          &
     +lc_o_cp*wql_inv(i)/inv_m1%exner_rho(i)

  !
  ! LEM derived wh_inv
  !
  wh_inv(i) = jump_inv_a2*scales%mb_new(i)*dh_inv(i)

  !
  ! LEM derived wthetav_inv
  !
  wthetav_inv(i) = jump_inv_a3*scales%mb_new(i)*dthetav_inv2(i)

end do

! Note at base of inversion expect wtheta & wthetal to be negative and
! wq, wqt to be positive

!-----------------------------------------------------------------------
! Tracer flux at inversion - is this accurate (Do something slightly
! different for dtheta and dq across inversion ?)
!
if (l_tracer) then
  do ktra = 1,ntra
    do i = 1,n_xx
      wtracer_inv(i,ktra) = -jump_inv_a1*scales%mb_new(i)*                     &
         (inv%tracer(i,ktra) - inv_m1%tracer(i,ktra))
    end do
  end do
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_inversion_warm

end module tcs_inversion
