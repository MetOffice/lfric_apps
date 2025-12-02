! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculates increments to U and V due to deep convection.

module deep_cmt_incr_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'DEEP_CMT_INCR_MOD'
contains

subroutine deep_cmt_incr(np_field,npnts,nconv,nlevs,nterm,                     &
                         nlcl,ntop,cu_term,cu_tend,                            &
                         zlcl,phalf,                                           &
                         uw,vw,                                                &
                         ! Output
                         dubydt,dvbydt)

use planet_constants_mod, only:g
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!-----------------------------------------------------------------------
! Description :
!  Calculates increments to U and V due to deep convection.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.
!
!------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  np_field             & ! Full field length
 ,npnts                & ! Total number of points in segment
 ,nconv                & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,nlcl(nconv)          & ! Lifting condensation level
 ,ntop(nconv)          & ! Top level of convection
 ,cu_term(nterm)       & ! Indices for terminating points
 ,cu_tend(nterm)         ! Indices of points in output array

real(kind=real_umphys), intent(in)    ::                                       &
  zlcl(npnts)          & ! Height of LCL (m)
 ,phalf(nlevs,nconv)   & ! Pressure on model half levels (Pa)
 ,uw(nlevs,nconv)      & ! U-component of stress (Pam/s2)
 ,vw(nlevs,nconv)        ! V-component of stress (Pam/s2)


real(kind=real_umphys), intent(out) ::                                         &
  dubydt(np_field,nlevs) & ! U increment (m/s2)
 ,dvbydt(np_field,nlevs)   ! V increment (m/s2)

! Local variables

integer ::                                                                     &
  i,j,k,m,n         ! Loop counters

real(kind=real_umphys) ::                                                      &
  dudt_bl        & ! U CMT tendency in subcloud-layer (m/s2)
 ,dvdt_bl        & ! V CMT tendency in subcloud-layer (m/s2)
 ,dp               ! Pressure difference between adjacent half levels (Pa)


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEEP_CMT_INCR'

!------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------
! Calculate U and V wind increments by differentiating stress profile

do i=1,nterm
  m=cu_term(i)
  n=cu_tend(i)
  j=nlcl(m)

  ! CMT tendencies in the subcloud layer are constant with height

  dudt_bl=-uw(j,m)/(g*zlcl(n))
  dvdt_bl=-vw(j,m)/(g*zlcl(n))
  do k=1,j-1
    dubydt(n,k)=dudt_bl
    dvbydt(n,k)=dvdt_bl
  end do
  do k=j,ntop(m)+1
    dp=-(phalf(k+1,m)-phalf(k,m))
    dubydt(n,k)=-(uw(k+1,m)-uw(k,m))/dp
    dvbydt(n,k)=-(vw(k+1,m)-vw(k,m))/dp
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine deep_cmt_incr
end module deep_cmt_incr_mod
