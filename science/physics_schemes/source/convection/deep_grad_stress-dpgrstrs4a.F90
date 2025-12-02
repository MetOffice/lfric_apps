! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the gradient component of the stress profile for deep CMT
!
module deep_grad_stress_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'DEEP_GRAD_STRESS_MOD'
contains

subroutine deep_grad_stress(nconv,nlevs,nlcl,ntop,                             &
                            nterm,cu_term,                                     &
                            ue,ve,visc,phalf,p,timestep,                       &
                            ! Output
                            uw,vw)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use tridiag_mod, only: tridiag
implicit none

!-----------------------------------------------------------------------
! Description :
!   To calculate the gradient component of the stress profile for deep
!   convection. Calculation an be done explicitly or implicitly.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
!------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  nconv                & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nlcl(nconv)          & ! Lifting condensation level
 ,ntop(nconv)          & ! Top level of convection
 ,nterm                & ! Number of points terminating
 ,cu_term(nterm)         ! Indices for terminating points

real(kind=real_umphys), intent(in)    ::                                       &
  ue(nlevs,nconv)      & ! Environment U-wind component (m/s)
 ,ve(nlevs,nconv)      & ! Environment V-wind component (m/s)
 ,visc(nlevs,nconv)    & ! Viscosity
 ,phalf(nlevs,nconv)   & ! Pressure on model half levels (hPa)
 ,p(nlevs,nconv)       & ! Pressure on model levels (hPa)
 ,timestep               ! Model timestep (s)

real(kind=real_umphys), intent(out) ::                                         &
  uw(nlevs,nconv)      & ! U-component of viscous stress
 ,vw(nlevs,nconv)        ! V-component of viscous stress

integer       ::                                                               &
  i,j,k,m        &  ! loop counters
  ,nlev             ! Number of levels

real(kind=real_umphys) ::                                                      &
  a(nlevs)        & ! Implicit solver variables
 ,b(nlevs)        & !
 ,c(nlevs)        & !
 ,u_t(nlevs)      & ! Current U wind
 ,u_tp1(nlevs)    & ! U Wind at T+1
 ,v_t(nlevs)      & ! Current V wind
 ,v_tp1(nlevs)      ! V wind at T+1

real(kind=real_umphys) ::                                                      &
  ue_tp1(nlevs,nconv)   & ! Implicitly updated U wind
 ,ve_tp1(nlevs,nconv)     ! Implicitly updated V wind


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEEP_GRAD_STRESS'
!------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------

! Use implicit timestepping

do i=1,nterm
  j=cu_term(i)
  nlev=0

  ! Calculate components of tridiagonal matirx and construct a vector of
  ! current timestep wind components

  do k=nlcl(j),ntop(j)+1
    nlev=nlev+1
    if (k == nlcl(j)) then
      a(nlev)=0.0
      c(nlev)=-visc(k+1,j)*timestep/                                           &
              ((p(k+1,j)-p(k,j))*(phalf(k+1,j)-phalf(k,j)))
    else if (k <= ntop(j)) then
      a(nlev)=-visc(k,j)*timestep/                                             &
              ((p(k,j)-p(k-1,j))*(phalf(k+1,j)-phalf(k,j)))
      c(nlev)=-visc(k+1,j)*timestep/                                           &
              ((p(k+1,j)-p(k,j))*(phalf(k+1,j)-phalf(k,j)))
    else if (k == ntop(j)+1) then
      a(nlev)=-visc(k,j)*timestep/                                             &
              ((p(k,j)-p(k-1,j))*(phalf(k+1,j)-phalf(k,j)))
      c(nlev)=0.0
    end if
    b(nlev)=1.0-a(nlev)-c(nlev)
    u_t(nlev)=ue(k,j)
    v_t(nlev)=ve(k,j)
  end do

  ! Calculate new timestep wind components using tridiag

  call tridiag(a,b,c,u_t,u_tp1,nlev)
  call tridiag(a,b,c,v_t,v_tp1,nlev)

  ! Store updated wind components for later

  nlev=0
  do k=nlcl(j),ntop(j)+1
    nlev=nlev+1
    ue_tp1(k,j)=u_tp1(nlev)
    ve_tp1(k,j)=v_tp1(nlev)
  end do
end do      ! nterm

! Calculate stress profiles

do i=1,nterm
  m=cu_term(i)
  j=nlcl(m)
  uw(j,m)=0.0
  vw(j,m)=0.0
  do k=j+1,ntop(m)+1
    uw(k,m)=-visc(k,m)*(ue_tp1(k,m)-ue_tp1(k-1,m))/(p(k-1,m)-p(k,m))
    vw(k,m)=-visc(k,m)*(ve_tp1(k,m)-ve_tp1(k-1,m))/(p(k-1,m)-p(k,m))
  end do
  uw(ntop(m)+2,m)=0.0
  vw(ntop(m)+2,m)=0.0
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine deep_grad_stress
end module deep_grad_stress_mod
