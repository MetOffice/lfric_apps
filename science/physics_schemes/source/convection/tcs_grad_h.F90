! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calculate dX/dz at time t+1 where X could be h or thetav
! or some other field. Required as part of flux calculations eg w'qt'
! and w'thetal'
!
module tcs_grad_h


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
! Module to calculate the gradient of the turbulent fluxes
! w'theta' and w'q' on in-cloud levels
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


character(len=*), parameter, private :: ModuleName='TCS_GRAD_H'

contains

subroutine calc_grad_h(n_xx, max_cldlev, nlev, timestep                        &
   ,                      nclev                                                &
   ,                      z_rho, z_theta                                       &
   ,                      r2rho, r2rho_theta                                   &
   ,                      x_t, k_f                                             &
   ,                      wx_inv                                               &
   ,                      kdXdz )

use tridiag_all_mod, only: tridiag_all
implicit none
!------------------------------------------------------------------
! Subroutine Arguments
!------------------------------------------------------------------
!
! Arguments with intent in:
!
integer, intent(in) ::                                                         &
   n_xx                                                                        &
                     ! No. of congestus convection points
   , max_cldlev                                                                &
                     ! Maximum number of convective cloud levels
   , nlev                                                                      &
                     ! Maximum number of convective cloud levels
   , nclev(n_xx)
                     ! cloud levels for point

real(kind=real_umphys), intent(in) ::                                          &
   timestep                 ! timestep for convection
real(kind=real_umphys), intent(in) ::                                          &
   x_t(n_xx,nlev)                                                              &
                            ! field X at time t
   , z_rho(n_xx,nlev)                                                          &
                            ! height of rho levels above surface
   , z_theta(n_xx,nlev)                                                        &
                            ! height of theta levels
   , r2rho(n_xx,nlev)                                                          &
                            ! r2rho  on rho levels
   , r2rho_theta(n_xx,nlev)                                                    &
                            ! r2rho on theta levels
   , k_f(n_xx,nlev)                                                            &
                            ! K on rho levels ?
   , wx_inv(n_xx)
                            ! wx at inversion
!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent out:
!
real(kind=real_umphys), intent(out) ::                                         &
   kdXdz(n_xx,nlev)  ! gradient of field X at t+1

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
real(kind=real_umphys) ::                                                      &
   a(n_xx,max_cldlev)                                                          &
   , b(n_xx,max_cldlev)                                                        &
   , c(n_xx,max_cldlev)                                                        &
   , x_tp1(n_xx,max_cldlev)                                                    &
                            ! X at T+1
   , h_t(n_xx,max_cldlev)                                                      &
   , dz_lower, dz_upper, dz_rho

!-------------------------
! Loop counters
!-------------------------
integer :: i,k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_GRAD_H'

!-----------------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------------
! Initialise kdXdz output array

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 1,nlev
  do i = 1,n_xx
    kdXdz(i,k) = 0.0
  end do
end do

!-----------------------------------------------------------------------
! 2.0 Implict calculation of X and dX/dz
!-----------------------------------------------------------------------
! Note -assume fluxes w'X' are zero at cloud base and inversion
! This assumption implies values of H at levels above and below
! cloud base are not required
!-----------------------------------------------------------------------
! Problem w'X' not zero at inversion
! At present assuming -KdX/dz part is zero - may not be true
!
! w'X'(inv) = -K(inv)[x_t(above)-x_t(below)]/[z(above)-z(below)]
!                          +w'X'cb FNG(inv)
!-----------------------------------------------------------------------

k=1
do i=1,n_xx
  dz_rho =( z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)

  dz_upper = z_theta(i,k+1) - z_theta(i,k)

  c(i,k) = -k_f(i,k+1)*r2rho(i,k+1)                                            &
     *timestep/(dz_rho*dz_upper)

  a(i,k) = 0.0            ! w'h' flux=0 at cloud base
  b(i,k) = 1.0 -a(i,k) - c(i,k)
  h_t(i,k) = x_t(i,k)
end do

do k=2,max_cldlev
  do i=1,n_xx
    if ((nclev(i)-1) >= k) then
      dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)

      dz_upper = z_theta(i,k+1) - z_theta(i,k)
      dz_lower = z_theta(i,k) - z_theta(i,k-1)

      c(i,k) = -k_f(i,k+1)*r2rho(i,k+1)                                        &
         *timestep/(dz_rho*dz_upper)
      a(i,k) = -k_f(i,k)*r2rho(i,k)                                            &
         *timestep/(dz_rho*dz_lower)
      b(i,k) = 1.0 -a(i,k) - c(i,k)

      h_t(i,k) = x_t(i,k)          !
    else if (nclev(i) == k) then

      dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)
      dz_lower = z_theta(i,k) - z_theta(i,k-1)
      c(i,k) = 0.0
      a(i,k) = -k_f(i,k)*r2rho(i,k)                                            &
         *timestep/(dz_rho*dz_lower)
      b(i,k) = 1.0 -a(i,k) - c(i,k)
      h_t(i,k) = x_t(i,k)
    else
      ! elements not required in calculation (zero)
      c(i,k) = 0.0
      a(i,k) = 0.0
      b(i,k) = 0.0
      h_t(i,k) = 0.0

    end if

  end do
end do

! CALCULATE NEW TIMESTEP X (moist static energy) using tridiag
!
call TRIDIAG_all(max_cldlev,n_xx,nclev,a,b,c,h_t,x_tp1)


!
! Take account of flux at inversion
!

do i=1,n_xx
  k=nclev(i)
  dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)
  x_tp1(i,k)=x_tp1(i,k)                                                        &
     - wx_inv(i)*r2rho(i,k+1)*timestep/dz_rho

end do

!
! Calculate gradient of X from t+1 values
! Only applies to cloud interior.
!
!  kdXdz= KdX/dz   ie Kterm
! Holding values from cloud base to last in cloud level
! (ie not inversion)

do i = 1,n_xx
  kdXdz(i,1) = 0.0     ! cloud base k_f=0.0
end do


do k = 2,max_cldlev
  do i = 1,n_xx
    if (k <= nclev(i)) then
      dz_lower = z_theta(i,k) - z_theta(i,k-1)
      kdXdz(i,k) =k_f(i,k)*(x_tp1(i,k) - x_tp1(i,k-1))                         &
         /dz_lower
    end if
  end do
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_grad_h

end module tcs_grad_h
