! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the gradient component of the stress due to shallow convection
!

module shallow_grad_stress_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'SHALLOW_GRAD_STRESS_MOD'
contains

subroutine shallow_grad_stress(npnts,n_cumulus,nterm,nlevs,                    &
                       cu_ind,nlcl,ntop,mb,wsc,wstr,zcld,plcl,                 &
                       ptop,p,phalf,rho,ue,ve,timestep,weight_param,           &
                      ! Outputs
                       uw,vw)

use planet_constants_mod, only: g

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use tridiag_mod, only: tridiag
implicit none

!------------------------------------------------------------------------
! Description:
!
! Calculates the gradient component of the stress due to shallow
! convection.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts                & ! Total number of points in segment
 ,n_cumulus            & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,cu_ind(nterm)        & ! Indices for terminating points
 ,nlcl(n_cumulus)      & ! Lifting condensation level
 ,ntop(n_cumulus)        ! Top level of convection

real(kind=real_umphys), intent(in)    ::                                       &
  mb(n_cumulus)            & ! Cloud base mass flux for shallow cu (m/s)
 ,wsc(n_cumulus)           & ! Cloud-layer velocity scale (m/s)
 ,wstr(n_cumulus)          & ! Mixed layer velocity scale (m/s)
 ,zcld(npnts)              & ! cloud layer depth (m)
 ,plcl(n_cumulus)          & ! Pressure at lifting condensation level (Pa)
 ,ptop(n_cumulus)          & ! Pressure at top of cloud layer (Pa)
 ,p(nlevs,n_cumulus)       & ! Pressure on model levels (Pa)
 ,phalf(nlevs,n_cumulus)   & ! Pressure on model half levels (Pa)
 ,rho(nlevs,n_cumulus)     & ! Density, model uv levels (kg/m3/s)
 ,ue(nlevs,n_cumulus)      & ! U-component of mean wind (m/s)
 ,ve(nlevs,n_cumulus)      & ! V-component of mean wind (m/s)
 ,weight_param(n_cumulus)  & ! Grey zone weighting factor
 ,timestep                   ! Model timestep (s)

real(kind=real_umphys), intent(out) ::                                         &
  uw(nlevs,n_cumulus)      & ! U-component of stress (m2/s2)
 ,vw(nlevs,n_cumulus)        ! V-component of stress (m2/s2)

! Local variables

integer ::                                                                     &
  i,j,k,m,nlev         ! Loop counters


real(kind=real_umphys) ::                                                      &
  w(nlevs,n_cumulus)        & ! Non-dimensional plume vertical velocity
 ,mass(nlevs,n_cumulus)     & ! Mass flux profile (m/s)
 ,visc(nlevs,nterm)         & ! Viscosity profile (m2/s)
 ,a(nlevs)                  & ! Tridiagonal matrix elements
 ,b(nlevs)                  & !
 ,c(nlevs)                  & !
 ,u_t(nlevs)                & ! Current velocity vectors (m/s)
 ,v_t(nlevs)                & !
 ,u_tp1(nlevs)              & ! After timestep velocity vectors
 ,v_tp1(nlevs)              & !
 ,ue_tp1(nlevs,n_cumulus)   & ! After timestep velocity vectors
 ,ve_tp1(nlevs,n_cumulus)   & ! After timestep velocity vectors
 ,w02                       & !
 ,p_depth                   & !
 ,zeta                      & !
 ,entr_sc                   & !
 ,exp_k                     & !
 ,exp_kp1                   & !
 ,dz                        & !
 ,dz12                        !

! Parameters (in future consider putting in a module).

real(kind=real_umphys), parameter ::                                           &
  a_stress=0.3      & !
 ,a_w02=10.24         !

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SHALLOW_GRAD_STRESS'

!------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------
! Calculate vertical velocity profile in updraughts

do i=1,nterm
  j=cu_ind(i)
  w02=a_w02*(mb(j)*wstr(j)**2)**0.6667
  p_depth=ptop(j)-plcl(j)
  do k=nlcl(j),ntop(j)
    zeta=(phalf(k,j)-plcl(j))/p_depth
    dz = min(6.0*zeta,6.0)
    w(k,i)=sqrt(w02+wsc(j)**2*dz)/wsc(j)
  end do
end do

! Calculates mass flux profile
! Uses fractional detrainment=1.3  fractional entrainment

do i=1,nterm
  j=cu_ind(i)
  entr_sc=0.04*wsc(j)/(mb(j)*zcld(j))
  p_depth=ptop(j)-plcl(j)
  mass(nlcl(j),i)=mb(j)
  exp_k=exp(-(phalf(nlcl(j),j)-plcl(j))/p_depth)
  do k=nlcl(j),ntop(j)-1
    exp_kp1=exp(-(phalf(k+1,j)-plcl(j))/p_depth)
    zeta=(1.0-1.3)*entr_sc*p_depth*(exp_kp1-exp_k)/(g*rho(k+1,j))
    mass(k+1,i)=mass(k,i)*exp(zeta)
    exp_k=exp_kp1
  end do
end do

! Calculate the eddy viscosity profile

do i=1,nterm
  j=cu_ind(i)
  do k=nlcl(j)+1,ntop(j)+1
    if (k <  ntop(j)) then
      visc(k,i)=a_stress*mass(k,i)*w(k,i)*zcld(j)
    else if (k == ntop(j)) then
      visc(k,i)=0.162*mb(j)*zcld(j)
    else if (k == ntop(j)+1) then
      dz=-(p(k,j)-p(k-1,j))/(g*(rho(k,j)+rho(k-1,j))/2.0)
      visc(k,i)=0.09*mb(j)*dz
    end if
    visc(k,i)=weight_param(j)*visc(k,i)
  end do
end do

! Calculate gradient component of stress

! Use implicit timestepping

do i=1,nterm
  j=cu_ind(i)
  nlev=0

  ! Calculate components of tridiagnol matrix and construct vector of
  ! current timestep wind components

  do k=nlcl(j),ntop(j)+1
    nlev=nlev+1
    if (k == nlcl(j)) then
      dz=-(phalf(k+1,j)-phalf(k,j))/(g*rho(k,j))
      dz12=-(p(k+1,j)-p(k,j))/(g*(rho(k+1,j)+rho(k,j))/2.0)
      a(nlev)=0.0
      c(nlev)=-visc(k+1,i)*timestep/(dz*dz12)
    else if (k <= ntop(j)) then
      dz=-(phalf(k+1,j)-phalf(k,j))/(g*rho(k,j))
      dz12=-(p(k,j)-p(k-1,j))/(g*(rho(k,j)+rho(k-1,j))/2.0)
      a(nlev)=-visc(k,i)*timestep/(dz*dz12)
      dz12=-(p(k+1,j)-p(k,j))/(g*(rho(k+1,j)+rho(k,j))/2.0)
      c(nlev)=-visc(k+1,i)*timestep/(dz*dz12)
    else if (k == ntop(j)+1) then
      dz=-(phalf(k+1,j)-phalf(k,j))/(g*rho(k,j))
      dz12=-(p(k,j)-p(k-1,j))/(g*(rho(k,j)+rho(k-1,j))/2.0)
      a(nlev)=-visc(k,i)*timestep/(dz*dz12)
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
end do   ! nterm

! Calculate stress profiles

#if defined(NEC)
! Initialise UW,VW
do j = 1,n_cumulus+1
  do m = 1,nlevs
    uw(m,j) = 0.0
    vw(m,j) = 0.0
  end do
end do
#endif
do i=1,nterm
  j=cu_ind(i)
  m=nlcl(j)
#if !defined(NEC)
  uw(m,j)=0.0
  vw(m,j)=0.0
#endif
  !CDIR NODEP
  do k=m+1,ntop(j)+1
    dz=-(p(k,j)-p(k-1,j))/(g*(rho(k,j)+rho(k-1,j))/2.0)
    uw(k,j)=-visc(k,i)*(ue_tp1(k,j)-ue_tp1(k-1,j))/dz
    vw(k,j)=-visc(k,i)*(ve_tp1(k,j)-ve_tp1(k-1,j))/dz
  end do
#if !defined(NEC)
  uw(ntop(j)+2,j)=0.0
  vw(ntop(j)+2,j)=0.0
#endif
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine shallow_grad_stress

end module shallow_grad_stress_mod
