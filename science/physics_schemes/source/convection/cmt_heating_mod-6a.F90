! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the heating due to the dissipation of kinetic energy in
! the convective momentum transport
!
module cmt_heating_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Calculates the heating due to the dissipation of kinetic energy in
!   the convective momentum transport
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


character(len=*), parameter, private :: ModuleName='CMT_HEATING_MOD'

contains

subroutine cmt_heating(npnts, nlev,                                            &
                       z_theta, z_rho, exner_layer_centres,                    &
                       u, v, dubydt, dvbydt,                                   &
                       ! Out
                       dthbydt)

use planet_constants_mod, only: cp

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

integer,intent(in) :: npnts               ! Number of points
integer,intent(in) :: nlev                ! Number of model layers

real(kind=real_umphys),intent(in) :: z_theta(npnts,nlev)
                                          ! height of theta levels (m)
real(kind=real_umphys),intent(in) :: z_rho(npnts,nlev)
                                          ! height of rho levels (m)
real(kind=real_umphys),intent(in) :: exner_layer_centres(npnts,0:nlev)
                                          ! exner pressure at layer centres
real(kind=real_umphys),intent(in) :: u(npnts,nlev)
                                          ! Model U field (m/s)
real(kind=real_umphys),intent(in) :: v(npnts,nlev)
                                          ! Model V field (m/s)
real(kind=real_umphys),intent(in) :: dubydt(npnts,nlev)
                                          ! Increments to U due to CMT (m/s2)
real(kind=real_umphys),intent(in) :: dvbydt(npnts,nlev)
                                          ! Increments to V due to CMT (m/s2)

!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------

real(kind=real_umphys),intent(in out) :: dthbydt(npnts,nlev)
                                           ! Increments to potential temp.
                                          ! due to convection (K/s)

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i, k ! loop counters

real(kind=real_umphys) :: dzkm14        ! Thickness of level k-1/4 (m)
real(kind=real_umphys) :: dzkp14        ! Thickness of level k+1/4 (m)
real(kind=real_umphys) :: dzk           ! Thickness of level k (m)
real(kind=real_umphys) :: uhat
                      ! U-wind interpolated to theta level (m/s)
real(kind=real_umphys) :: vhat
                      ! V-wind interpolated to theta level (m/s)
real(kind=real_umphys) :: ududt         ! u*du/dt on theta level (m2/s3)
real(kind=real_umphys) :: vdvdt         ! v*dv/dt on theta level (m2/s3)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CMT_HEATING'

!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 1, nlev-1
  do i = 1, npnts

    dzkm14 = z_theta(i,k)   -  z_rho(i,k)
    dzkp14 = z_rho(i,k+1)   -  z_theta(i,k)
    dzk    = z_rho(i,k+1)   -  z_rho(i,k)

    uhat   = ( dzkp14 * u(i,k) + dzkm14 * u(i,k+1) ) / dzk
    vhat   = ( dzkp14 * v(i,k) + dzkm14 * v(i,k+1) ) / dzk

    ududt  = uhat * ( dzkp14 * dubydt(i,k) + dzkm14 * dubydt(i,k+1) ) / dzk
    vdvdt  = vhat * ( dzkp14 * dvbydt(i,k) + dzkm14 * dvbydt(i,k+1) ) / dzk

    !      dT/dt on theta_level(k)
    dthbydt(i,k)  = dthbydt(i,k) - (ududt + vdvdt) /                           &
                   (cp * exner_layer_centres(i,k))
  end do !Loop over points
end do !Loop over levels

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine cmt_heating
end module cmt_heating_mod
