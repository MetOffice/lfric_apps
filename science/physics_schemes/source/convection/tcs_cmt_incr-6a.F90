! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to calculate CMT increments to U and V.
!
module tcs_cmt_incr_6a


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
! This module calculates the cmt increments to U and V winds for
! the tcs warm rain convection schem
!
! Method:
!   <reference to documentation to go here, once available>
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.


character(len=*), parameter, private :: ModuleName='TCS_CMT_INCR_6A'

contains

subroutine calc_cmt_incr(n_xx, nlevs,ntpar_max                                 &
   ,                            r2rho, r2rho_th                                &
   ,                            dr_across_rh                                   &
   ,                            uw,vw                                          &
                              ! output arguements
   ,                            dubydt,dvbydt)

implicit none
!-------------------------------------------------------------------
! Subroutine Arguments
!-------------------------------------------------------------------
!
! Arguments with intent in:
!
integer, intent(in) ::                                                         &
   n_xx                                                                        &
   ! total number of congestus convective points
   , nlevs                                                                     &
   ! number of model levels
   , ntpar_max
   ! maximum cloud top level

real(kind=real_umphys), intent(in) ::                                          &
   r2rho(n_xx,nlevs)                                                           &
                            ! r2 rho on rho levels
   , r2rho_th(n_xx,nlevs)                                                      &
                            ! r2 rho on theta levels
   , dr_across_rh(n_xx,nlevs)                                                  &
                            ! dr across rho levels
   , uw(n_xx,nlevs)                                                            &
                            ! U-Component of stress profile (M2)
   , vw(n_xx,nlevs)
! V-Component of stress profile (M2)
!
! Arguments with intent out:
!
real(kind=real_umphys), intent(out) ::                                         &
   dubydt(n_xx,nlevs+1)                                                        &
                            ! Tendency in U (MS-2)
   , dvbydt(n_xx,nlevs+1)
                            ! Tendency in V (MS-2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
!
real(kind=real_umphys) ::                                                      &
   rhodz                    ! r2 rho * dz

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_CMT_INCR'

!-------------------------
! Loop counters
!-------------------------
integer :: i,k

!
!-----------------------------------------------------------------------
! Input to this routine contains stress profile of uw & vw on levels
! from the surface to the top of the inversion.
!-----------------------------------------------------------------------
!
! Lowest uv level surface stress zero
!
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

k=1
do i=1,n_xx
  rhodz  = r2rho(i,k)*dr_across_rh(i,k)
  dubydt(i,k) = -uw(i,k)*r2rho_th(i,k)/rhodz
  dvbydt(i,k) = -vw(i,k)*r2rho_th(i,k)/rhodz
end do

!
! Mixed layer and Cloud layer increments
!
do k=2,ntpar_max+1
  do i=1,n_xx

    rhodz  = r2rho(i,k)*dr_across_rh(i,k)
    dubydt(i,k) =                                                              &
       -(r2rho_th(i,k)*uw(i,k)-r2rho_th(i,k-1)*uw(i,k-1))/rhodz
    dvbydt(i,k) =                                                              &
       -(r2rho_th(i,k)*vw(i,k)-r2rho_th(i,k-1)*vw(i,k-1))/rhodz
  end do
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_cmt_incr

end module tcs_cmt_incr_6a
