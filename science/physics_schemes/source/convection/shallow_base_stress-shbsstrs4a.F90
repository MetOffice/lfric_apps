! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates cloud base stress for shallow convection
!

module shallow_base_stress_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'SHALLOW_BASE_STRESS_MOD'
contains

subroutine shallow_base_stress(np_field,npnts,n_cumulus,nlevs,                 &
                               nterm,                                          &
                               cu_ind,cu_full,nlcl,ntop,mb,wsc,                &
                               zlcl,uw0,vw0,plcl,ptop,ue,                      &
                               ve,phalf,p,rho,timestep,weight_param,           &
                                     ! in/out ARGUMENTS
                               uw,vw,                                          &
                                     ! OUTPUT ARGUMENTS
                               uw_shall,vw_shall)

use planet_constants_mod, only: g
use cv_stash_flg_mod,    only: flg_uw_shall, flg_vw_shall
use model_domain_mod,    only: model_type, mt_single_column, mt_lfric

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!------------------------------------------------------------------------
! Description:
!   Calculates cloud base stress for shallow convection
!   (also completes caluclation of stress profile).
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
  np_field             & ! Full field length
 ,npnts                & ! Total number of points in segment
 ,n_cumulus            & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,cu_ind(nterm)        & ! Indices for terminating points
 ,cu_full(nterm)       & ! Indices of points in output array
 ,nlcl(n_cumulus)      & ! Lifting condensation level
 ,ntop(n_cumulus)        ! Top level of convection

real(kind=real_umphys), intent(in)    ::                                       &
  mb(n_cumulus)            & ! Cloud base mass flux for shallow cu (m/s)
 ,wsc(n_cumulus)           & ! Cloud-layer velocity scale (m/s)
 ,zlcl(npnts)              & ! Height of LCL (m)
 ,uw0(npnts)               & ! U component of surface stress (m2/s2)
 ,vw0(npnts)               & ! V component of surface stress (m2/s2)
 ,plcl(n_cumulus)          & ! Pressure at lifting condensation level (Pa)
 ,ptop(n_cumulus)          & ! Pressure at top of cloud layer (Pa)
 ,ue(nlevs,n_cumulus)      & ! U-component of mean wind (m/s)
 ,ve(nlevs,n_cumulus)      & ! V-component of mean wind (m/s)
 ,phalf(nlevs,n_cumulus)   & ! Pressure on model half levels (Pa)
 ,p(nlevs,n_cumulus)       & ! Pressure on model levels (Pa)
 ,rho(nlevs,n_cumulus)     & ! Density, model uv levels (kg/m3/s)
 ,weight_param(n_cumulus)  & ! Grey zone weighting factor
 ,timestep                   ! Model timestep (s)

real(kind=real_umphys), intent(in out) ::                                      &
  uw(nlevs,n_cumulus)      & ! U-component of stress (m2/s2)
 ,vw(nlevs,n_cumulus)        ! V-component of stress (m2/s2)

real(kind=real_umphys), intent(out) ::                                         &
  uw_shall(np_field,nlevs) & ! Stash diagnostic for U comp stress
 ,vw_shall(np_field,nlevs)   ! Stash diagnostic for V comp stress

! Local variables

integer ::                                                                     &
  i,j,k,n          ! Loop counters


real(kind=real_umphys) ::                                                      &
  omg2_jump(nterm)     & ! Jump in Y component of vorticity
 ,omg1_jump(nterm)     & ! Jump in X component of vorticity
 ,zeta                 & !
 ,dz                   & ! height difference
 ,p_depth              & !
 ,a                    & !
 ,b                    & !
 ,c                    & !
 ,t                    & !
 ,du                   & !
 ,dv                   & !
 ,dz1                  & ! height difference
 ,rho_h                  !


! Parameters (in future consider putting in a module).

real(kind=real_umphys), parameter ::                                           &
  beta=0.04         & !
 ,delta=2.3         & !
 ,gamma_c=1.63

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SHALLOW_BASE_STRESS'

!------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------
! Calculate jumps in vorticity across cloud base (this is done by assuming
! that during the timestep du and dv vary as exp(-T/TAU). Needs to be done
! to avoid instability around cloud base).

if ((model_type == mt_single_column) .or. (model_type == mt_lfric)) then
  do i=1, nterm
    j       = cu_ind(i)
    n       = cu_full(i)
    dz      = -(p(nlcl(j),j)-plcl(j))/(g*rho(nlcl(j),j))
    du      = (ue(nlcl(j),j)-ue(nlcl(j)-1,j))
    dv      = (ve(nlcl(j),j)-ve(nlcl(j)-1,j))
    dz1     = -(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/(g*rho(nlcl(j)+1,j))
    p_depth = (ptop(j)-plcl(j))
    zeta    = beta*wsc(j)*(phalf(nlcl(j)+1,j)-plcl(j))/(mb(j)*p_depth)
    b       = (1.0/zlcl(n)-(exp(-zeta)-1.0)/dz1)
    a       = zlcl(n)*mb(j)*b/(delta*dz)
    c       = (b*(1.0-gamma_c/delta)-1.0/zlcl(n))*uw0(n)

    if (c == 0.0 .and. (dv == 0.0 .or. du == 0.0)) then
      omg2_jump=0.0
    else
      t = -log((c*(1.0-exp(-a*timestep))/a+                                    &
           du*(exp(-a*timestep)-1.0))/                                         &
           ((c-a*du)*timestep))/a
      omg2_jump(i) = (c*(1.0-exp(-a*t))/a+du*exp(-a*t))/dz
    end if
    c = (b*(1.0-gamma_c/delta)-1.0/zlcl(n))*vw0(n)
    if (c == 0.0 .and. (dv == 0.0 .or. du == 0.0)) then
      omg1_jump=0.0
    else
      t = -log((c*(1.0-exp(-a*timestep))/a+                                    &
           dv*(exp(-a*timestep)-1.0))/                                         &
           ((c-a*dv)*timestep))/a
      omg1_jump(i) = -(c*(1.0-exp(-a*t))/a+dv*exp(-a*t))/dz
    end if
  end do
else !not scm or LFRic
  do i=1, nterm
    j       = cu_ind(i)
    n       = cu_full(i)
    dz      = -(p(nlcl(j),j)-plcl(j))/(g*rho(nlcl(j),j))
    du      = (ue(nlcl(j),j)-ue(nlcl(j)-1,j))
    dv      = (ve(nlcl(j),j)-ve(nlcl(j)-1,j))
    dz1     = -(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/(g*rho(nlcl(j)+1,j))
    p_depth = (ptop(j)-plcl(j))
    zeta    = beta*wsc(j)*(phalf(nlcl(j)+1,j)-plcl(j))/(mb(j)*p_depth)
    b       = (1.0/zlcl(n)-(exp(-zeta)-1.0)/dz1)
    a       = zlcl(n)*mb(j)*b/(delta*dz)
    c       = (b*(1.0-gamma_c/delta)-1.0/zlcl(n))*uw0(n)

    t = -log((c*(1.0-exp(-a*timestep))/a+                                      &
          du*(exp(-a*timestep)-1.0))/                                          &
          ((c-a*du)*timestep))/a
    omg2_jump(i) = (c*(1.0-exp(-a*t))/a+du*exp(-a*t))/dz

    c = (b*(1.0-gamma_c/delta)-1.0/zlcl(n))*vw0(n)

    t = -log((c*(1.0-exp(-a*timestep))/a+                                      &
          dv*(exp(-a*timestep)-1.0))/                                          &
          ((c-a*dv)*timestep))/a
    omg1_jump(i) = -(c*(1.0-exp(-a*t))/a+dv*exp(-a*t))/dz
  end do
end if !model_type

! Calculate the cloud-base stress components
do i=1,nterm
  j=cu_ind(i)
  n=cu_full(i)
  uw(nlcl(j),j) = zlcl(n)*( -mb(j)*omg2_jump(i) -                              &
                             gamma_c*uw0(n)/zlcl(n) )/delta + uw0(n)
  vw(nlcl(j),j) = zlcl(n)*( mb(j)*omg1_jump(i) -                               &
                            gamma_c*vw0(n)/zlcl(n)  )/delta + vw0(n)
  uw(nlcl(j),j) = weight_param(j)*uw(nlcl(j),j)
  vw(nlcl(j),j) = weight_param(j)*vw(nlcl(j),j)
end do

! Calculate non-gradient stress profile

do i=1,nterm
  j=cu_ind(i)
  n=cu_full(i)
  p_depth=(ptop(j)-plcl(j))
  do k=nlcl(j)+1,ntop(j)+1
    zeta=beta*wsc(j)*(phalf(k,j)-plcl(j))/(mb(j)*p_depth)
    rho_h=rho(k-1,j)+(rho(k,j)-rho(k-1,j))/(p(k,j)-p(k-1,j))*                  &
                                     (phalf(k,j)-p(k-1,j))
    if (k <  ntop(j)) then
      uw(k,j)=rho_h*(uw(k,j)+uw(nlcl(j),j)*exp(-zeta))
      vw(k,j)=rho_h*(vw(k,j)+vw(nlcl(j),j)*exp(-zeta))
    else if (k == ntop(j)) then
      uw(k,j)=rho_h*(uw(k,j)+uw(nlcl(j),j)*exp(-beta*wsc(j)/mb(j)))
      vw(k,j)=rho_h*(vw(k,j)+vw(nlcl(j),j)*exp(-beta*wsc(j)/mb(j)))
    else if (k == ntop(j)+1) then
      uw(k,j)=rho_h*(uw(k,j)+uw(nlcl(j),j)*                                    &
                          exp(-1.75*beta*wsc(j)/mb(j)))
      vw(k,j)=rho_h*(vw(k,j)+vw(nlcl(j),j)*                                    &
                          exp(-1.75*beta*wsc(j)/mb(j)))
    end if
  end do

  ! Weight cloud base stress by rho (omitted from above level loop)
  ! Needs to be done after level loop.

  k=nlcl(j)
  rho_h=rho(k-1,j)+(rho(k,j)-rho(k-1,j))/(p(k,j)-p(k-1,j))*                    &
                            (phalf(k,j)-p(k-1,j))
  uw(nlcl(j),j) = rho_h*uw(nlcl(j),j)
  vw(nlcl(j),j) = rho_h*vw(nlcl(j),j)


  if (flg_uw_shall) then
    do k=nlcl(j),ntop(j)+1
      uw_shall(n,k)=uw(k,j)
    end do
  end if
  if (flg_vw_shall) then
    do k=nlcl(j),ntop(j)+1
      vw_shall(n,k)=vw(k,j)
    end do
  end if
end do   ! nterm

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine shallow_base_stress
end module shallow_base_stress_mod
