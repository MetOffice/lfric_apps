! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing tcs warm rain subroutine for calculating
! cloud base stress
!
module tcs_base_stress


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
!   This routine calculates cloud base stress
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
!

character(len=*), parameter, private :: ModuleName='TCS_BASE_STRESS'

contains

subroutine calc_base_stress(n_xx, nlevs, ntml, ntpar,ntpar_max                 &
   ,                              timestep                                     &
   ,                              scales                                       &
   ,                              uw0,vw0,du_cb,dv_cb                          &
   ,                              rho_theta,zrho,ztheta                        &
   ,                              flg_uw_cong,flg_vw_cong                      &
                              ! in/out ARGUMENTS
   ,                              uw,vw                                        &
                              ! OUTPUT ARGUMENTS
   ,                              uw_cong,vw_cong)


use tcs_parameters_warm,   only:                                               &
   beta_cmt, delta_cmt, gamma_cmt
use tcs_class_scales,      only:                                               &
   scales_conv

implicit none

!--------------------------------------------------------------------
! Subroutine Arguments
!--------------------------------------------------------------------
!
! Arguments with intent in:
!
integer, intent(in) ::                                                         &
   n_xx                                                                        &
                            ! Total number of congestus points
   , nlevs                                                                     &
                            ! Number of model levels
   , ntml(n_xx)                                                                &
                            ! levels of LCL
   , ntpar(n_xx)                                                               &
                            ! levels of TOP ofcloud layer
   , ntpar_max
! max value of ntpar +1
!
logical, intent(in) ::                                                         &
   flg_uw_cong                                                                 &
                            ! STASH FLAGS FOR CONGESTUS
   ,flg_vw_cong
! CONVECTION STRESS DIAGNOSTIC
!
real(kind=real_umphys), intent(in) ::                                          &
   timestep         ! MODEL timestep (S)

type(scales_conv), intent(in) :: scales

real(kind=real_umphys), intent(in) ::                                          &
   uw0(n_xx)                                                                   &
                            ! U-component of surface stress (M2S-2)
   ,  vw0(n_xx)                                                                &
                            ! V-component of surface stress (M2S-2)
   ,  du_cb(n_xx)                                                              &
                            ! dU across cloud base (m/s)
   ,  dv_cb(n_xx)                                                              &
                            ! dV across cloud base (m/s)
   ,  rho_theta(n_xx,nlevs)                                                    &
                            ! Density model th levels (kgm-3)
   ,  zrho(n_xx,nlevs)                                                         &
                            ! height of model rho levels (m)
   ,  ztheta(n_xx,nlevs)
! height of model theta levels (m)

!
! Arguments with intent INOUT:
!
real(kind=real_umphys), intent(in out) ::                                      &
   uw(n_xx,nlevs)                                                              &
                            ! U-component of STRESS PROFILE (M2S-2)
   , vw(n_xx,nlevs)
! V-component of STRESS PROFILE (M2S-2)
!
! Arguments with intent out:
!
real(kind=real_umphys), intent(out) ::                                         &
   uw_cong(n_xx,nlevs)                                                         &
                            ! STASH DIAGNOSTIC FOR U-COMP STREss
   , vw_cong(n_xx,nlevs)
! STASH DIAGNOSTIC FOR V-COMP STREss
!
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

real(kind=real_umphys) ::                                                      &
   omg2_jump(n_xx)                                                             &
                            ! jump in Y component of vorticity
   , omg1_jump(n_xx)                                                           &
                            ! jump in X-component of vorticity
   , zlcl_cmt(n_xx)                                                            &
                            ! lcl for CMT    (m)
   , fcmt,dz                                                                   &
   , a,b,c                                                                     &
                            ! coefficients
   , t,dz1                                                                     &
   , z_depth(n_xx)                                                             &
   , expadt                  ! exp (Adt)

!-------------------------
! Loop counters
!-------------------------
integer :: i,k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_BASE_STRESS'

!-------------------------------------------------------------------
!
! Calculate jumps in vorticity across cloud base. (This is done by
! assuming that during the time step du and dv vary as exp(-T/TAU).
! Needs to be done to avoid instability around cloud base.)
!
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i=1,n_xx

  zlcl_cmt(i) = ztheta(i,ntml(i))   ! cloud base for CMT
  !
  ! dz (cloud base - theta level below cloud base (zlcl_cmt) )
  !
  dz = zrho(i,ntml(i)+1)-zlcl_cmt(i)
  !
  ! dz across cloud base
  !
  dz1 = ztheta(i,ntml(i)+1) - zlcl_cmt(i)
  !
  ! depth of cloud as defined for uv calculations (different from th &q)
  ! calculations ?
  !
  z_depth(i) = ztheta(i,ntpar(i)) - zlcl_cmt(i)

  !  f(z/scales%zcld) = exp(-fcmt)
  fcmt=beta_cmt*scales%wsc_o_mb(i)*(dz1/z_depth(i))
  b=(1.0/zlcl_cmt(i)-(exp(-fcmt)-1.0)/dz1)
  ! alpha in doc (but extra /dz factor)
  a=zlcl_cmt(i)*scales%mb(i)*b/(delta_cmt*dz)
  expadt = exp(-a*timestep)

  ! beta in documentation

  c = ( b*(1.0-gamma_cmt/delta_cmt) - 1.0/zlcl_cmt(i) )*uw0(i)

  if (c == 0.0 .and. du_cb(i) == 0.0) then
    omg2_jump=0.0
  else
    t=-log( (c*(1.0-expadt)/a+du_cb(i)*(expadt-1.0))/                          &
                             ((c-a*du_cb(i))*timestep))/a
    omg2_jump(i)=( c*(1.0-exp(-a*t))/a + du_cb(i)*exp(-a*t) )/dz
  end if

  ! beta in documentation
  c = ( b*(1.0-gamma_cmt/delta_cmt) - 1.0/zlcl_cmt(i) )*vw0(i)

  if (c == 0.0 .and. dv_cb(i) == 0.0) then
    omg1_jump=0.0
  else
    t=-log((c*(1.0-expadt)/a+dv_cb(i)*(expadt-1.0))/                           &
                             ((c-a*dv_cb(i))*timestep))/a
    omg1_jump(i)=-(c*(1.0-exp(-a*t))/a+dv_cb(i)*exp(-a*t))/dz
  end if
end do
!
! Calculate the cloud-base stress components
! Equations 11 & 12 section 5.1
!

do i=1,n_xx
  uw(i,ntml(i))=zlcl_cmt(i)*(-scales%mb(i)*omg2_jump(i)-                       &
                       gamma_cmt*uw0(i)/zlcl_cmt(i))/delta_cmt+uw0(i)
  vw(i,ntml(i))=zlcl_cmt(i)*(scales%mb(i)*omg1_jump(i)-                        &
                       gamma_cmt*vw0(i)/zlcl_cmt(i))/delta_cmt+vw0(i)
end do

!
! Calculate non-gradient stress profile in cloud
! Altered numbering of uw arrays
!
! New form of Fcmt   where alpha >1?  (like fng term for thermo?)
! Fcmt(z/zlcd)  = exp(-alpha*(z/scales%zcld))

do k=1,ntpar_max

  do i=1,n_xx


    if (k >= (ntml(i)+1) .and. k <= (ntpar(i)-1)) then
      ! F function

      fcmt=exp(-1.1*(ztheta(i,k)-zlcl_cmt(i))/z_depth(i))

      ! all cloud levels add non-gradient term to eddy viscosity term

      uw(i,k)=uw(i,k)+uw(i,ntml(i))*fcmt
      vw(i,k)=vw(i,k)+vw(i,ntml(i))*fcmt

    else if (k <= (ntml(i)-1)) then

      ! This if assuming uw on rho levels
      uw(i,k) = uw(i,ntml(i))*zrho(i,k)/zlcl_cmt(i)
      vw(i,k) = vw(i,ntml(i))*zrho(i,k)/zlcl_cmt(i)

    end if

  end do

end do

!
! Copy stress to output arrays multiplying by density on theta levels
!
if (flg_uw_cong) then
  do k=1,ntpar_max+1
    do i=1,n_xx
      uw_cong(i,k)=uw(i,k)*rho_theta(i,k)
    end do
  end do
end if
if (flg_vw_cong) then
  do k=1,ntpar_max+1
    do i=1,n_xx
      vw_cong(i,k)=vw(i,k)*rho_theta(i,k)
    end do
  end do
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_base_stress

end module tcs_base_stress
