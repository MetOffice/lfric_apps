! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

module shallow_cmt_incr_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'SHALLOW_CMT_INCR_MOD'
contains

subroutine shallow_cmt_incr(np_field,npnts,n_cumulus,nlevs,nterm,              &
                            cu_ind,cu_full,nlcl,ntop,uw,vw,phalf,              &
                            rho,zlcl,                                          &
                            ! Output
                            dubydt,dvbydt)

use planet_constants_mod, only: g

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!------------------------------------------------------------------------
! Description:
!   Calculates increments to U and V due to shallow convection.
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
  uw(nlevs,n_cumulus)      & ! U-component of stress (m2/s2)
 ,vw(nlevs,n_cumulus)      & ! V-component of stress (m2/s2)
 ,phalf(nlevs,n_cumulus)   & ! Pressure on model half levels (Pa)
 ,rho(nlevs,n_cumulus)     & ! Density, model uv levels (kg/m3/s)
 ,zlcl(npnts)                ! Height of LCL (m)


real(kind=real_umphys), intent(out) ::                                         &
  dubydt(np_field,nlevs)   & ! U increment (m/s2)
 ,dvbydt(np_field,nlevs)     ! V increment (m/s2)

! Local variables

integer ::                                                                     &
  i,j,k,m        ! Loop counters


real(kind=real_umphys) ::                                                      &
  dudt_bl        & ! U CMT tendency in subcloud-layer (m/s2)
 ,dvdt_bl        & ! V CMT tendency in subcloud-layer (m/s2)
 ,rhodz            ! Density times height difference between half level

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SHALLOW_CMT_INCR'
!------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------

do i=1,nterm
  j=cu_ind(i)
  m=cu_full(i)
  dudt_bl=-uw(nlcl(j),j)/zlcl(m)
  dvdt_bl=-vw(nlcl(j),j)/zlcl(m)

  ! Mixed layer increments (assumed constant in mixed layer)

  do k=1,nlcl(j)-1
    dubydt(m,k)=dudt_bl/rho(k,j)
    dvbydt(m,k)=dvdt_bl/rho(k,j)
  end do

  ! Cloud layer increments

  do k=nlcl(j),ntop(j)
    rhodz=-(phalf(k+1,j)-phalf(k,j))/g
    dubydt(m,k)=-(uw(k+1,j)-uw(k,j))/rhodz
    dvbydt(m,k)=-(vw(k+1,j)-vw(k,j))/rhodz
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine shallow_cmt_incr
end module shallow_cmt_incr_mod
