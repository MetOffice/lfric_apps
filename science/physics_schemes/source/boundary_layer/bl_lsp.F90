! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Convert temperature from liquid ice to liquid, and convert
!           the vapour+liquid+ice variable (Q) to vapour+liquid. This
!           subroutine is used if the mixed phase precipitation scheme
!           is selected and a full boundary layer treatment is not
!           performed.

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

module bl_lsp_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'BL_LSP_MOD'
contains

subroutine bl_lsp( bl_levels,qcf,q,t )

use atm_fields_bounds_mod, only: tdims
use planet_constants_mod, only: lsrcp
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

integer, intent(in) ::                                                         &
  bl_levels             ! in   Number of boundary layer levels

real(kind=real_umphys), intent(in out) ::                                      &
  qcf(tdims%i_start:tdims%i_end,                                               &
      bl_levels),                                                              &
                                 ! INOUT Ice water content
  q(tdims%i_start:tdims%i_end,                                                 &
    bl_levels),                                                                &
                                 ! INOUT
!                                  in    Vapour+liquid+ice content
!                                  out   Vapour+liquid content
    t(tdims%i_start:tdims%i_end,                                               &
      bl_levels)                   ! INOUT
!                                  in    Liquid ice temperature
!                                  out   Liquid temperature
! Temporary Space
integer ::                                                                     &
        i,                                                                     &
                               ! Counter over points
        k                ! Counter over boundary layer levels
real(kind=real_umphys) :: newqcf              ! Temporary variable for QCF

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BL_LSP'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP          private(i,k,newqcf)                                             &
!$OMP          SHARED(bl_levels,tdims,q,qcf,t,lsrcp)
do k = 1, bl_levels
  do i = tdims%i_start, tdims%i_end
    ! Convert Q (vapour+liquid+ice) to (vapour+liquid)
    q(i,k)=q(i,k)-qcf(i,k)
    ! Check that Q is not negative
    if (q(i,k)  <   0.0) then
      ! Evaporate ice to keep Q positive, but don't let ice go negative
      ! itself
      newqcf=max(qcf(i,k)+q(i,k),0.0)
      q(i,k)=q(i,k)+(qcf(i,k)-newqcf)
      qcf(i,k)=newqcf
    end if
    ! Adjust T from T liquid ice to T liquid
    t(i,k)=t(i,k)+lsrcp*qcf(i,k)
  end do
end do
!$OMP end PARALLEL do
! End the subroutine
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bl_lsp
end module bl_lsp_mod
