! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calculate precipation for warm rain tcs convection
!
module tcs_precip


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
! Module to calculate precipation for warm rain tcs convection
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

character(len=*), parameter, private :: ModuleName='TCS_PRECIP'

contains

subroutine calc_precip( n_xx, nlev, conv_type                                  &
   ,                     ql_ad, dz_inv, rho_inv                                &
   ,                     eta_rho,eta_theta                                     &
   ,                     rainfall, wqr_inv, precip_product_inv                 &
   ,                     wqr, precip_product_th, precip_product_uv )

  !-----------------------------------------------------------------
  !
  ! Description:
  !   Calculate the warm rain precipitation
  !
  !-----------------------------------------------------------------

use tcs_parameters_warm,     only:                                             &
   ql_t, epsilon_rain, beta_rain, gamma_rain, kappa_rain
use tcs_constants,           only:                                             &
   g
use tcs_common_warm ,        only:                                             &
   scales

implicit none
!-----------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------
integer, intent(in) ::                                                         &
   n_xx                                                                        &
                            ! No. of congestus convection points
   ,  nlev
! No. model levels

integer, intent(in) ::                                                         &
   conv_type(n_xx)
! Integer index describing convective type:
!    1=non-precipitating shallow
!    2=drizzling shallow
!    3=warm congestus

real(kind=real_umphys), intent(in) ::                                          &
   ql_ad(n_xx)                                                                 &
                            ! adiabatic liquid water content (kg/kg)
   , dz_inv(n_xx)                                                              &
                            ! dz across inversion (m)
   , rho_inv(n_xx)                                                             &
                            ! density at inversion (kg/m3)
   , eta_rho(:,:)                                                              &
                            ! eta rho levels
   , eta_theta(:,:)
! eta theta levels

real(kind=real_umphys), intent(out) ::                                         &
   rainfall(n_xx)                                                              &
                            ! surface rainfall rate
   , wqr_inv(n_xx)                                                             &
                            ! flux of w'qr' at inversion
   , precip_product_inv(n_xx)                                                  &
                            ! Integral of precip_product over inv
   , wqr(n_xx,nlev)                                                            &
                            !  w'qr'
   , precip_product_th(n_xx,nlev)                                              &
                            ! precipitation production rate
   , precip_product_uv(n_xx,nlev) ! precipitation production rate


!-----------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------
real(kind=real_umphys)  ::                                                     &
   frac_of_ens(n_xx)                                                           &
                  ! fraction of ensemble with lwc above threshold
   , ql_up_inv(n_xx)                                                           &
                  ! mean liquid water content of ensemble at base
                  ! of inversion
   , p0(n_xx)
                  ! precip productio p0 value
real(kind=real_umphys)  ::                                                     &
   fr(n_xx,nlev)  ! gravitional flux of precipitation


integer :: neta   ! number of levels in eta_rho/eta_theta
integer :: i,k    ! loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_PRECIP'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
neta=size(eta_rho(1,:))

!-----------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------
precip_product_th(:,:) = 0.0
precip_product_uv(:,:) = 0.0
fr(:,:) = 0.0
wqr(:,:) = 0.0
frac_of_ens(:) = 0.0
rainfall(:) = 0.0
p0(:) = 0.0
wqr_inv(:) = 0.0
precip_product_inv(:) = 0.0

!-----------------------------------------------------------------
! Calculate liquid water content of buoyant updraughts at the
! base of the inversion
!-----------------------------------------------------------------
ql_up_inv(:) = (1.5/g)*((scales%wstar_up(:)/scales%mb(:))**0.5)*               &
   scales%wstar_up(:)*scales%wstar_up(:)/scales%zcld(:)

! Only produce precip if adiabatic liquid water content is above
! threshold.  The fraction of the ensemble producing rain will
! further depend on whether the mean ql is above the threshold.
!$OMP PARALLEL private(i) DEFAULT(none)                                        &
!$OMP SHARED(n_xx, ql_ad, ql_t, conv_type, rainfall, epsilon_rain,             &
!$OMP scales,ql_up_inv, p0, beta_rain, dz_inv, wqr_inv, gamma_rain,            &
!$OMP precip_product_inv, kappa_rain, rho_inv, frac_of_ens)
!$OMP do SCHEDULE(STATIC)
do i=1,n_xx
  if (ql_ad(i) > ql_t .and. ql_up_inv(i) > ql_t) then
    frac_of_ens(i) = 1.0 - ql_t*ql_t/(ql_up_inv(i)*ql_ad(i))
  else if (ql_ad(i) > ql_t .and. ql_up_inv(i) <= ql_t) then
    frac_of_ens(i) = (1.0-ql_t/ql_ad(i))*(1.0-ql_t/ql_ad(i))/                  &
                    (1.0-ql_up_inv(i)/ql_ad(i))
  end if
end do
!$OMP end do NOWAIT

!$OMP do SCHEDULE(STATIC)
do i=1,n_xx
  if (ql_ad(i) > ql_t .and. conv_type(i) > 1) then
    ! Only produce precip if adiabatic
    ! liquid water content is above
    ! threshold
    ! and is of precipitating type
    !-------------------------------------------------------------
    ! surface precipitation rate
    !-------------------------------------------------------------
    rainfall(i) = epsilon_rain*scales%mb(i)*ql_up_inv(i)*frac_of_ens(i)

    p0(i) = rainfall(i)/(beta_rain*(1.0-0.5*beta_rain)                         &
       *scales%zcld(i)*scales%zcld(i) +                                        &
       0.5*beta_rain*scales%zcld(i)*dz_inv(i))

    !-------------------------------------------------------------
    ! wqr at inversion
    !-------------------------------------------------------------
    wqr_inv(i) = gamma_rain*rainfall(i)                                        &
       - 0.5*p0(i)*beta_rain*scales%zcld(i)*dz_inv(i)

    precip_product_inv(i) = kappa_rain*0.5*p0(i)*beta_rain                     &
       *scales%zcld(i)*dz_inv(i)/rho_inv(i)
  end if
end do
!$OMP end do
!$OMP end PARALLEL


!----------------------------------------------------------------
! calculate in cloud profiles of precipitation
!----------------------------------------------------------------
!--------------------------
! on rho levels
!--------------------------
do k=1,neta
  do i=1,n_xx
    if (eta_rho(i,k) <= beta_rain) then
      ! Lower portion of cloud layer
      precip_product_uv(i,k) = eta_rho(i,k)*scales%zcld(i)                     &
         *p0(i)
      fr(i,k) = rainfall(i)
      wqr(i,k) = fr(i,k) - rainfall(i)                                         &
         + precip_product_uv(i,k)*eta_rho(i,k)*scales%zcld(i)/2.0

    else if ( eta_rho(i,k) <= 1.0 + spacing(eta_rho(i,k)) ) then
      ! Upper portion of cloud layer
      precip_product_uv(i,k) = p0(i)*beta_rain                                 &
         *scales%zcld(i)
      fr(i,k) = (((eta_rho(i,k)-beta_rain)/(1.0 - beta_rain))                  &
         *(gamma_rain - 1.0) + 1.0 )*rainfall(i)
      wqr(i,k) = fr(i,k) - rainfall(i)                                         &
         + p0(i)*beta_rain*scales%zcld(i)                                      &
         * (eta_rho(i,k) - beta_rain/2.0)                                      &
         *scales%zcld(i)
    end if

    !--------------------------
    ! on theta levels
    !--------------------------

    if ( eta_theta(i,k) <= beta_rain ) then
      precip_product_th(i,k) = eta_theta(i,k)*scales%zcld(i)                   &
         *p0(i)
    else if ( eta_theta(i,k) <= 1.0 + spacing(eta_rho(i,k)) ) then
      precip_product_th(i,k) = p0(i)*beta_rain                                 &
         *scales%zcld(i)
    end if
  end do
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_precip

end module tcs_precip
