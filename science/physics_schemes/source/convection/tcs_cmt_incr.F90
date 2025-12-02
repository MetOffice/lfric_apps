! *****************************COPYRIGHT*******************************
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

module tcs_cmt_incr_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'TCS_CMT_INCR_MOD'
contains

subroutine tcs_cmt_incr(n_npts, nlevs, ntop_max, kterm                         &
                        ,r2rho, r2rho_th                                       &
                        ,dr_across_rh                                          &
                        ,uw,vw                                                 &
                                       ! output arguements
                                       !
                        ,dubydt,dvbydt)

! Description:
!
!  To calculate increments to U and V due to shallow convection.
!


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent in:
!
integer, intent(in) ::                                                         &
  n_npts               & ! total number of convective columns
, nlevs                & ! number of model levels
, ntop_max             & ! maximum cloud top level
, kterm(n_npts)          ! top of convection theta levels

real(kind=real_umphys), intent(in) ::                                          &
  r2rho(n_npts,nlevs)        & ! r2 rho on rho levels (kg/m)
, r2rho_th(n_npts,nlevs)     & ! r2 rho on theta levels (kg/m)
, dr_across_rh(n_npts,nlevs) & ! dr across rho levels (m)
, uw(n_npts,nlevs)           & ! U-Component of stress profile (m2/s2)
, vw(n_npts,nlevs)             ! V-Component of stress profile (m2/s2)

!
! Arguments with intent out:
!
real(kind=real_umphys), intent(out) ::                                         &
  dubydt(n_npts,nlevs)       & ! Tendency in U (ms-2)
, dvbydt(n_npts,nlevs)         ! Tendency in V (ms-2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

integer :: i,k      ! loop counters

real(kind=real_umphys) ::   rhodz     ! r2 rho * dz

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='TCS_CMT_INCR'


!
!-----------------------------------------------------------------------
! Input to this routine contains stress profile of uw & vw on levels
! from the surface to the top of the cloud.
!-----------------------------------------------------------------------
!
! Lowest uv level surface stress zero
!
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
k=1
do i=1,n_npts
  if (0 < kterm(i)) then    ! Column has convection
    rhodz  = r2rho(i,k)*dr_across_rh(i,k)
    dubydt(i,k) = -uw(i,k)*r2rho_th(i,k)/rhodz
    dvbydt(i,k) = -vw(i,k)*r2rho_th(i,k)/rhodz
  end if
end do

!
! Mixed layer and Cloud layer increments
! Note restricted to levels where increments apply as otherwise get
! non-zero increments even though uw & vw = 0.0 due to numerics.
!
!CDIR NOUNROLL
do k=2,ntop_max+1    ! Should reduce cost where large stratosphere
  !   Do k=2,nlevs   !
  do i=1,n_npts
    if (k <= kterm(i)+1) then   ! level where increments apply
      rhodz  = r2rho(i,k)*dr_across_rh(i,k)
      dubydt(i,k) =                                                            &
           -(r2rho_th(i,k)*uw(i,k)-r2rho_th(i,k-1)*uw(i,k-1))/rhodz
      dvbydt(i,k) =                                                            &
           -(r2rho_th(i,k)*vw(i,k)-r2rho_th(i,k-1)*vw(i,k-1))/rhodz
    end if
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine tcs_cmt_incr
end module tcs_cmt_incr_mod
