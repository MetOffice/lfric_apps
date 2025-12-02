! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calculate the gradient of the turbulent fluxes
! w'theta' and w'q' on in-cloud levels
!
module tcs_wql

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
!   Module to calculate the liquid water flux
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

character(len=*), parameter, private :: ModuleName='TCS_WQL'

contains

subroutine calc_wql_warm(n_xx, nlev, max_cldlev, ncld_thlev, wql_cb,           &
   wql_inv, cld_in, t_cld, sim, wql_cld, dqsatdt_cld)

use tcs_constants,          only:                                              &
    Lc, Rv, repsilon, g
use tcs_class_similarity,   only:                                              &
    similarity
use tcs_class_cloud,        only:                                              &
    cloud_input
use tcs_common_warm,        only:                                              &
    scales

implicit none
!----------------------------------------------------------------
! Subroutine Arguments
!----------------------------------------------------------------
integer, intent(in) ::                                                         &
   n_xx                                                                        &
                      ! No. of congestus convection points
   , nlev                                                                      &
                      ! No. of model layers
   , max_cldlev
                      ! maximum noumber of in cloud levels

integer, intent(in) ::                                                         &
   ncld_thlev(n_xx)   ! number of theta cloud levels

real(kind=real_umphys), intent(in)    ::                                       &
   wql_cb(n_xx)                                                                &
                      ! flux of wql across cloud base
   , wql_inv(n_xx)
                      ! flux of wql across inversion


type(cloud_input), intent(in) :: cld_in


real(kind=real_umphys), intent(in) :: t_cld(n_xx,nlev) ! temperature in cloud

type(similarity), intent(in) :: sim  ! Similarity functions

real(kind=real_umphys), intent(in out) ::                                      &
   wql_cld(n_xx,nlev)     ! wql in cloud (includes cloud base)

real(kind=real_umphys), intent(out) ::                                         &
   dqsatdt_cld(n_xx,max_cldlev+1)    ! dqsat/dt

!----------------------------------------------------------------
! Variables defined locally
!----------------------------------------------------------------
real(kind=real_umphys) ::                                                      &
   thetav_cld(n_xx,max_cldlev+1)                                               &
   , temp(n_xx)

!-------------------------
! Loop counters
!-------------------------
integer :: i,k,itop

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_WQL_WARM'

!----------------------------------------------------------------
! 1.0  Calculate dqsat/dT  in cloud
!----------------------------------------------------------------

!----------------------------------------------------------------
!   dqsatdt on cloud levels  (expression for mixing ratio)
!----------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 1,max_cldlev+1
  do i = 1,n_xx
    dqsatdt_cld(i,k) = cld_in%qse(i,k)*lc                                      &
       /(Rv*t_cld(i,k)*t_cld(i,k))
    thetav_cld(i,k)=cld_in%theta(i,k)*(1.0+cld_in%q_mix(i,k)/repsilon)         &
       /(1.0+cld_in%q_mix(i,k))
  end do
end do

!----------------------------------------------------------------
! 2.0 Calculate  in cloud values of wql on uv levels using
!      similarity expression
!----------------------------------------------------------------
! values for new water flux parametrisation

k=1
do i=1,n_xx
  wql_cld(i,k) = wql_cb(i)
  temp(i) =scales%root_mb_new_o_wsc(i)*(scales%wstar_up(i)**3)/(scales%zcld(i)*g)
end do

! taking thetav/theta=1.

do k=2,max_cldlev+1
  do i=1,n_xx
    wql_cld(i,k) = 0.5*temp(i)*sim%fql_func_rho(i,k)*                          &
       cld_in%theta(i,k)/thetav_cld(i,k)
  end do
end do

! value at inversion

do i=1,n_xx
  itop=ncld_thlev(i)+1
  wql_cld(i,itop) = wql_inv(i)
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_wql_warm

end module tcs_wql
