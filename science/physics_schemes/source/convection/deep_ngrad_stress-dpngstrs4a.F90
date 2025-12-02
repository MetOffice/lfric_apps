! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the cloud base stress for deep CMT
!

module deep_ngrad_stress_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'DEEP_NGRAD_STRESS_MOD'
contains

subroutine deep_ngrad_stress(np_field,npnts,nconv,nterm,nlevs,                 &
                             nlcl,ntop,cu_term,cu_comp,cu_tend,                &
                             pstar,uw0,vw0,zlcl,ue,ve,                         &
                             mass,p,phalf,rho,timestep,                        &
                             ! Input/output
                             uw,vw,                                            &
                             ! Output
                             uw_base,vw_base,uw_dp,vw_dp)

use cv_stash_flg_mod, only:                                                    &
    flg_uw_dp, flg_vw_dp

use cv_param_mod, only:                                                        &
    dp_cmt_gamma, dp_cmt_delta

use model_domain_mod, only: model_type, mt_single_column, mt_lfric
use planet_constants_mod, only: g

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!-----------------------------------------------------------------------
! Description :
!   To calculate the cloud base stress for deep convection and complete
!   the calculation of the stress profile.
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
  np_field             & ! Full field length
 ,npnts                & ! Total number of points in segment
 ,nconv                & ! Number of convecting points
 ,nterm                & ! Number of points terminating
 ,nlevs                & ! Number of model levels
 ,nlcl(nconv)          & ! Lifting condensation level
 ,ntop(nconv)          & ! Top level of convection
 ,cu_comp(npnts)       & ! Index array for convecting points
 ,cu_term(nterm)       & ! Indices for terminating points
 ,cu_tend(nterm)         ! Index of points in output array

real(kind=real_umphys), intent(in)    ::                                       &
  pstar(npnts)         & ! surface pressure
 ,uw0(npnts)           & ! Surface shear stress x-component (m2/s2)
 ,vw0(npnts)           & ! Surface shear stress x-component (m2/s2)
 ,zlcl(npnts)          & ! Height of LCL (m)
 ,ue(nlevs,nconv)      & ! Environment U-wind component (m/s)
 ,ve(nlevs,nconv)      & ! Environment V-wind component (m/s)
 ,mass(nlevs,nconv)    & ! Updraught mass flux (Pa/s)
 ,p(nlevs,nconv)       & ! Pressure on model levels (hPa)
 ,phalf(nlevs,nconv)   & ! Pressure on model half levels (hPa)
 ,rho(nlevs,nconv)     & ! Density, model uv levels (kg/m3/s)
 ,timestep               ! Model timestep (s)

real(kind=real_umphys), intent(in out) ::                                      &
  uw(nlevs,nconv)         & ! U-component of stress
 ,vw(nlevs,nconv)           ! V-component of stress

real(kind=real_umphys), intent(out) ::                                         &
  uw_base(nconv)          & ! cloud base U-component of stress
 ,vw_base(nconv)          & ! cloud base V-component of stress
 ,uw_dp(np_field,nlevs)   & ! U-component of stress for stash
 ,vw_dp(np_field,nlevs)     ! V-component of stress for stash


! Local variables

integer ::                                                                     &
  i              & ! local array index
 ,k              & ! Level index
 ,j              & ! Indexes points in compressed input arrays
 ,m              & !
 ,n                ! Indexes points in uncompressed input arrays

real(kind=real_umphys) ::                                                      &
  omg2_jump(nterm) & ! Cloud base jump in Y component of vorticity
 ,omg1_jump(nterm) & ! Cloud base jump in X component of vorticity
 ,mb(nterm)        & ! Cloud base mass flux (m/s)
 ,dz               & !
 ,beta             & !
 ,du               & !
 ,dv               & !
 ,dz1              & !
 ,zeta             & !
 ,a                & !
 ,b                & !
 ,c                & !
 ,t                  !

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEEP_NGRAD_STRESS'
!------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------
! Convert cloud-base mass flux from Pa/s to m/s

do i=1,nterm
  j=cu_term(i)
  n=cu_comp(i)
  k=nlcl(j)
  mb(i)=mass(k,j)/g
end do

! Calculate jumps in vorticity components across cloud-base.
! 'Implicit technique' assumes du, dv vary as exp(-t/tau) through timestep
! needed because explicit calculation can lead to instability in du, dv
! under some circumstances.
if ((model_type == mt_single_column) .or. (model_type == mt_lfric)) then
  do i=1,nterm
    j    = cu_term(i)
    n    = cu_comp(i)
    dz   = -(p(nlcl(j),j)-phalf(nlcl(j),j))/(g*rho(nlcl(j),j))
    du   =  (ue(nlcl(j),j)-ue(nlcl(j)-1,j))
    dv   =  (ve(nlcl(j),j)-ve(nlcl(j)-1,j))
    dz1  = -(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/(g*rho(nlcl(j)+1,j))
    zeta = -(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/25000.0
    b    = (1.0/zlcl(n)-(exp(-zeta)-1.0)/dz1)
    a    = zlcl(n)*mb(i)*b/(dp_cmt_delta*dz)
    c    = (b*(1.0-dp_cmt_gamma/dp_cmt_delta)-1.0/zlcl(n))*uw0(n)

    if (c == 0.0 .and. (dv == 0.0 .or. du == 0.0)) then
      omg2_jump(i)=0.0
    else
      t=-log((c*(1.0-exp(-a*timestep))/a+                                      &
           du*(exp(-a*timestep)-1.0))/                                         &
           ((c-a*du)*timestep))/a
      omg2_jump(i)=(c*(1.0-exp(-a*t))/a+du*exp(-a*t))/dz
    end if

    c=(b*(1.0-dp_cmt_gamma/dp_cmt_delta)-1/zlcl(n))*vw0(n)
    if (c == 0.0 .and. (dv == 0.0 .or. du == 0.0)) then
      omg1_jump(i)=0.0
    else
      t=-log((c*(1.0-exp(-a*timestep))/a+                                      &
           dv*(exp(-a*timestep)-1.0))/                                         &
           ((c-a*dv)*timestep))/a
      omg1_jump(i)=-(c*(1.0-exp(-a*t))/a+dv*exp(-a*t))/dz
    end if

  end do

else !not scm or LFRic
  ! So far there have been no reported failures from full models due to
  ! du & uw0 or dv & vw0 being exactly zero and leading to a (c-a*du) giving
  ! a divide by zero.

  do i=1,nterm
    j    = cu_term(i)
    n    = cu_comp(i)
    dz   = -(p(nlcl(j),j)-phalf(nlcl(j),j))/(g*rho(nlcl(j),j))
    du   = (ue(nlcl(j),j)-ue(nlcl(j)-1,j))
    dv   = (ve(nlcl(j),j)-ve(nlcl(j)-1,j))
    dz1  = -(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/(g*rho(nlcl(j)+1,j))
    zeta = -(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/25000.0
    b    = (1.0/zlcl(n)-(exp(-zeta)-1.0)/dz1)
    a    = zlcl(n)*mb(i)*b/(dp_cmt_delta*dz)
    c    = (b*(1.0-dp_cmt_gamma/dp_cmt_delta)-1.0/zlcl(n))*uw0(n)

    t = -log((c*(1.0-exp(-a*timestep))/a+                                      &
         du*(exp(-a*timestep)-1.0))/                                           &
         ((c-a*du)*timestep))/a

    omg2_jump(i) = (c*(1.0-exp(-a*t))/a+du*exp(-a*t))/dz
    c = (b*(1.0-dp_cmt_gamma/dp_cmt_delta)-1/zlcl(n))*vw0(n)
    t = -log((c*(1.0-exp(-a*timestep))/a+                                      &
         dv*(exp(-a*timestep)-1.0))/                                           &
         ((c-a*dv)*timestep))/a
    omg1_jump(i) = -(c*(1.0-exp(-a*t))/a+dv*exp(-a*t))/dz

  end do
end if !model_type


! Calculate cloud base stress components. Note factor of g to convect
! back to Pa/s.

do i=1,nterm
  j=cu_term(i)
  n=cu_comp(i)
  uw_base(j)=g*(zlcl(n)*(-mb(i)*omg2_jump(i)-                                  &
                 dp_cmt_gamma*uw0(n)/zlcl(n))/dp_cmt_delta+uw0(n))
  vw_base(j)=g*(zlcl(n)*(mb(i)*omg1_jump(i)-                                   &
                 dp_cmt_gamma*vw0(n)/zlcl(n))/dp_cmt_delta+vw0(n))
end do

! Calculate total stress

! Calculate stress profiles the function beta was again tuned to TOGA-COARE
! CRM simulation

do i=1,nterm
  m=cu_term(i)
  n=cu_tend(i)
  j=nlcl(m)
  uw(j,m)=uw_base(m)
  vw(j,m)=vw_base(m)

  ! below cloud base
  do k=1,nlcl(m)-1
    beta=(phalf(k,m)-pstar(m))/(phalf(j,m)-pstar(m))
    uw(k,m)=uw(k,m)+beta*uw_base(m)
    vw(k,m)=vw(k,m)+beta*vw_base(m)
  end do

  do k=j+1,ntop(m)+1
    beta=exp(((phalf(k,m)-phalf(j,m))/25000.0))
    uw(k,m)=uw(k,m)+beta*uw_base(m)
    vw(k,m)=vw(k,m)+beta*vw_base(m)
  end do
  uw(ntop(m)+2,m)=0.0
  vw(ntop(m)+2,m)=0.0

  ! Stash diagnostics (NOTE diagnostics are in m/s for direct comparison
  ! with shallow convection stresses)

  if (flg_uw_dp) then
    do k=j,ntop(m)+2
      uw_dp(n,k)=uw(k,m)/g
    end do
  end if
  if (flg_vw_dp) then
    do k=j,ntop(m)+2
      vw_dp(n,k)=vw(k,m)/g
    end do
  end if
end do ! nterm

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


return
end subroutine deep_ngrad_stress
end module deep_ngrad_stress_mod
