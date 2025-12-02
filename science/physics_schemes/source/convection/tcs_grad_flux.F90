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
module tcs_grad_flux


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


character(len=*), parameter, private :: ModuleName='TCS_GRAD_FLUX'

contains

subroutine calc_grad_flux(n_xx, ntra, nlev, trlev, maxlev, l_tracer            &
   ,                      r2rho,r2rho_th,dr_across_th                          &
   ,                      wthetavl, wthetal, wqt,wtracer                       &
   ,                      dwthetavl_dz, dwthetal_dz, dwqt_dz                   &
   ,                      dwtracer_dz)

implicit none
!--------------------------------------------------------------------
! Subroutine Arguments
!--------------------------------------------------------------------
!
! Arguments with intent in:
!
integer, intent(in) ::                                                         &
   n_xx                                                                        &
                            ! No. of congestus convection points
   , ntra                                                                      &
                            ! No. of tracer
   , maxlev                                                                    &
                            ! Maximum number of levels where
                            ! gradient is non zero
   , nlev                                                                      &
                            ! Maximum number of convective cloud levels
   , trlev
                            ! Maximum number of tracer levels

logical, intent(in) ::                                                         &
   l_tracer                 ! true - tracers present


real(kind=real_umphys), intent(in) ::                                          &
   r2rho(n_xx,nlev)                                                            &
                            ! radius**2 density on rho lev (kg/m)
   , r2rho_th(n_xx,nlev)                                                       &
                            ! radius**2 density on theta lev (kg/m)
   , dr_across_th(n_xx,nlev)
                            !  thickness on theta levels (m)

!
! fluxes  all held on rho levels
!
real(kind=real_umphys), intent(in) ::                                          &
   wthetavl(n_xx,nlev)                                                         &
                            ! w'thetavl'
   , wthetal(n_xx,nlev)                                                        &
                            ! w'theta'
   , wqt(n_xx,nlev)                                                            &
                            ! w'qt'
   , wtracer(n_xx,trlev,ntra)
                            ! w'tracer' (kgm/kg/s)
!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent out:
!
real(kind=real_umphys), intent(out) ::                                         &
   dwthetavl_dz(n_xx,nlev)                                                     &
                            ! dwthetavl/dz  on theta levels
   ,dwthetal_dz(n_xx,nlev)                                                     &
                            ! dwthetal/dz  on theta levels
   ,dwqt_dz(n_xx,nlev)                                                         &
                            ! dwqt/dz      on theta levels
   ,dwtracer_dz(n_xx,trlev,ntra)
                            ! dwtracer/dz   on theta levels (kg/kg/s)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

integer :: max_trlev

real(kind=real_umphys) ::                                                      &
   rdz                                                                         &
                            ! 1/(dz*r2*rho)
   , r2_kp1                                                                    &
                            ! r**2 rho at k+1 levels
   , r2_k
                            ! r**2 rho at k levels


!-------------------------
! Loop counters
!-------------------------
integer :: i,k,ktra

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_GRAD_FLUX'

!-----------------------------------------------------------------------
! 1.0 Initialise arrays - set all functions to zero
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 1,nlev
  do i = 1,n_xx
    dwqt_dz(i,k) = 0.0
    dwthetal_dz(i,k) = 0.0
    dwthetavl_dz(i,k) = 0.0
  end do
end do

!-----------------------------------------------------------------------
! 2.0 Level loop to calculate gradients of fluxes
!-----------------------------------------------------------------------

do k = 1,maxlev    ! problem if nlev (no rho theta)
  do i = 1,n_xx

    rdz  =1.0/(dr_across_th(i,k)*r2rho_th(i,k))


    r2_kp1 = r2rho(i,k+1)
    r2_k   = r2rho(i,k)

    dwqt_dz(i,k)     = (r2_kp1*wqt(i,k+1)-r2_k*wqt(i,k))*rdz

    dwthetal_dz(i,k) = (r2_kp1*wthetal(i,k+1)                                  &
       -r2_k*wthetal(i,k))*rdz

    dwthetavl_dz(i,k)= (r2_kp1*wthetavl(i,k+1)                                 &
       -r2_k*wthetavl(i,k))*rdz
  end do
end do


!-----------------------------------------------------------------------
! 3.0 Tracers
!-----------------------------------------------------------------------


if (l_tracer) then

  ! Check max levels ?
  ! May be problems if tracers on less levels < maxlev ?
  max_trlev = maxlev
  if (trlev  <=  maxlev) then
    max_trlev = trlev-1
  end if

  do ktra = 1,ntra
    do k = 1,trlev
      do i = 1,n_xx
        dwtracer_dz(i,k,ktra) = 0.0
      end do
    end do

    do k = 1,max_trlev  ! problem if nlev (no rho theta)
      do i = 1,n_xx

        rdz  =1.0/(dr_across_th(i,k)*r2rho_th(i,k))
        r2_kp1 = r2rho(i,k+1)
        r2_k   = r2rho(i,k)
        dwtracer_dz(i,k,ktra) = (r2_kp1*wtracer(i,k+1,ktra)                    &
           -r2_k*wtracer(i,k,ktra))*rdz

      end do
    end do
  end do
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_grad_flux

end module tcs_grad_flux
