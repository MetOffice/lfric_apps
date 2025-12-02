! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calculate the gradient component of the stress
! for the tcs warm rain convection calculations
!
module tcs_grad_stress

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use umPrintMgr, only:                                                          &
    umPrint,                                                                   &
    umMessage
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

character(len=*), parameter, private :: ModuleName='TCS_GRAD_STRESS'

contains

subroutine calc_grad_stress(n_xx, nlev,max_cldlev                              &
   ,                       nclev                                               &
   ,                       timestep, scales                                    &
   ,                       mf_cld, w_up, sim                                   &
   ,                       u,v                                                 &
   ,                      r2rho,r2rho_th,dr_across_rh,dr_across_th             &
                              ! output arguements
   ,                       uw_cld,vw_cld)

use tcs_class_scales,         only:                                            &
   scales_conv
use tcs_class_similarity,     only:                                            &
   similarity

use tridiag_all_mod, only: tridiag_all
implicit none
!------------------------------------------------------------------
! Subroutine Arguments
!------------------------------------------------------------------
!
! Arguments with intent in:
!

integer, intent(in) ::                                                         &
   nlev                                                                        &
                            ! No. of model layers
   , n_xx                                                                      &
                            ! No. of congestus convection points
   , nclev(n_xx)                                                               &
                            ! number of cloud levels
   , max_cldlev
                            ! Maximum number of cloud levels

real(kind=real_umphys), intent(in) :: timestep   ! model timestep (s)

type(scales_conv), intent(in) :: scales

real(kind=real_umphys), intent(in) ::                                          &
   w_up(n_xx,nlev)                                                             &
                            ! ensemble vertical velocity (m/s)
   , mf_cld(n_xx,nlev)
                            ! mass flux in (th lev) ensemble (m/s)

type(similarity), intent(in) :: sim ! similarity functions

real(kind=real_umphys), intent(in) ::                                          &
   u(n_xx,nlev)                                                                &
                            ! U-component of mean wind (MS-1)
   , v(n_xx,nlev)
                            ! V-component of mean wind (MS-1)

real(kind=real_umphys), intent(in) ::                                          &
   r2rho(n_xx,nlev)                                                            &
   , r2rho_th(n_xx,nlev)                                                       &
   , dr_across_th(n_xx,nlev)                                                   &
   , dr_across_rh(n_xx,nlev)

!
! Arguments with intent inout:
!
!
! Arguments with intent out:
!
real(kind=real_umphys), intent(out) ::                                         &
   uw_cld(n_xx,nlev)                                                           &
                            ! U-component of stress
   , vw_cld(n_xx,nlev)      ! V-component of stress

!------------------------------------------------------------------
! Variables defined locally
!------------------------------------------------------------------

real(kind=real_umphys)      ::                                                 &
   visc(n_xx,max_cldlev+1)                                                     &
                            ! viscosity profile (M2S-1)
   , a(n_xx,max_cldlev+1)                                                      &
                            ! tridiagonal matrix elements
   , b(n_xx,max_cldlev+1)                                                      &
   , c(n_xx,max_cldlev+1)                                                      &
   , ue_tp1(n_xx,max_cldlev+1)                                                 &
                            ! after timestep velocity vectors
   , ve_tp1(n_xx,max_cldlev+1)                                                 &
                            ! for subsequent use
   , dz,dz12
integer ::                                                                     &
   nclevp1(n_xx)                                                               &
                            ! cloud level for uv cal
   , max_cldlev1
                            ! max cloud levels plus 1

!-------------------------
! Loop counters
!-------------------------
integer :: i,k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_GRAD_STRESS'

!
!------------------------------------------------------------------
! Note 4A code assumed a 2 level inversion now go to 1 level
! ie drop top condition.
!------------------------------------------------------------------
!
!  wup and mass_flux input on required levels ?
!
! Calculate the eddy viscosity profile K(z/scales%zcld)
!

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
max_cldlev1 = max_cldlev+1

do k=1,max_cldlev1  ! from level above cloud base to just
  do i=1,n_xx

    if (k <  nclev(i)) then         ! in cloud

      visc(i,k)=sim%k_func(i,k)*mf_cld(i,k)                                    &
         *W_up(i,k)*scales%zcld_uv(i)/scales%wstar_up(i)
    else if (k == nclev(i)) then     ! inversion ?

      dz = dr_across_th(i,k)
      visc(i,k)=0.09*scales%mb_new(i)*dz

    else
      visc(i,k)=0.0
    end if
  end do
end do

!
! Calculate GRADIENT component of stress
!
!
! Use implicit timestepping
!

k=1
do i=1,n_xx
  dz12 = dr_across_th(i,k)
  dz   = dr_across_rh(i,k)*r2rho(i,k)
  a(i,k)=0.0
  c(i,k)=-visc(i,k)*r2rho_th(i,k)                                              &
     *timestep/(dz*dz12)
  b(i,k)=1.0-a(i,k)-c(i,k)

  nclevp1(i) = nclev(i) + 1
end do
do k=2,max_cldlev1
  do i=1,n_xx
    dz = dr_across_rh(i,k)*r2rho(i,k)
    if (dz == 0.0) then
      write(umMessage,*) ' dz = 0',dr_across_rh(i,k),r2rho(i,k),i,k
      call umPrint(umMessage,src='tcs_grad_stress')
    end if
    if (k <= (nclev(i))) then
      dz12 = dr_across_th(i,k-1)
      if (dz12 == 0.0) then
        write(umMessage,*) ' dz12 = 0',dr_across_th(i,k-1),i,k
        call umPrint(umMessage,src='tcs_grad_stress')
      end if
      a(i,k)=-visc(i,k-1)*r2rho_th(i,k-1)                                      &
         *timestep/(dz*dz12)
      dz12 = dr_across_th(i,k)
      c(i,k)=-visc(i,k)*r2rho_th(i,k)                                          &
         *timestep/(dz*dz12)
    else if (k == (nclev(i)+1)) then
      dz12 = dr_across_th(i,k-1)
      a(i,k)=-visc(i,k-1)*r2rho_th(i,k-1)                                      &
         *timestep/(dz*dz12)
      c(i,k)=0.0
    else
      ! elements not required in calculation (zero)
      c(i,k) = 0.0
      a(i,k) = 0.0

    end if
    b(i,k)=1.0-a(i,k)-c(i,k)

  end do
end do

!
! Calculate NEW timestep wind conponents using tridiag
!
call TRIDIAG_all(max_cldlev1,n_xx,nclevp1,a,b,c,u,ue_tp1)
call TRIDIAG_all(max_cldlev1,n_xx,nclevp1,a,b,c,v,ve_tp1)
!
! Calculate stress profile -Kdu/dz from latest u/v values
!

! Initialise uw,vw
do k=1,max_cldlev1
  do i = 1,n_xx
    uw_cld(i,k) = 0.0
    vw_cld(i,k) = 0.0
  end do
end do

do k=1,max_cldlev1
  do i=1,n_xx
    if (k <= nclev(i)) then
      dz = dr_across_th(i,k)
      uw_cld(i,k)=-visc(i,k)*(ue_tp1(i,k+1)-ue_tp1(i,k))/dz
      vw_cld(i,k)=-visc(i,k)*(ve_tp1(i,k+1)-ve_tp1(i,k))/dz
    end if
  end do
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_grad_stress

end module tcs_grad_stress
