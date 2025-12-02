! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing tcs warm rain subroutine for calculating
! scalings for tcs scheme
!
module tcs_calc_scales


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
!   This routine calculates scalings for tcs warm rain calculations
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

character(len=*), parameter, private :: ModuleName='TCS_CALC_SCALES'

contains

subroutine calc_scales_warm(scales, conv_type, ntml, ntpar,                    &
   wstar_dn, delthvu, theta, q_mix,                                            &
   ztop_uv,                                                                    &
   zlcl_uv, z_rho)
!-----------------------------------------------------------------------
! Description:
!   Calculate scales for convection scheme - warm rain version.
!
!   Allocation for the scales array is done in the calling procedure
!-----------------------------------------------------------------------

use tcs_parameters_warm, only:                                                 &
   mb_a1, wup_a1
  !, mb_a2

use tcs_constants, only:                                                       &
   g, c_virtual

use tcs_class_scales, only:                                                    &
   scales_conv

implicit none
type(scales_conv), intent(in out) :: scales

integer, intent(in) :: conv_type(:)
              ! Indicator of type of convection
              !    1=non-precipitating shallow
              !    2=drizzling shallow
              !    3=warm congestus

integer, intent(in) ::  ntml(:)

integer, intent(in) ::  ntpar(:)

real(kind=real_umphys), intent(in) :: wstar_dn(:)

real(kind=real_umphys), intent(in) :: delthvu(:)

real(kind=real_umphys), intent(in) :: theta(:,:)

real(kind=real_umphys), intent(in) :: q_mix(:,:)

real(kind=real_umphys), intent(in) :: ztop_uv(:)

real(kind=real_umphys), intent(in) :: zlcl_uv(:)

real(kind=real_umphys), intent(in) :: z_rho(:,:)

!-------------------------------------------------------------------
! Loop counters
!-------------------------------------------------------------------
integer :: i

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_SCALES_WARM'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------
! Put wstar_dn (from the BL) into the scales structure.
!-------------------------------------------------------------------
scales%wstar_dn(:)=wstar_dn(:)

!-------------------------------------------------------------------
! Cloud base mass flux closure
!-------------------------------------------------------------------

scales%mb(:)  = mb_a1 * wstar_dn(:)
!
! The offset for congestus needs to be more carefully
! considered and related to evaporation of precip/Ri
!

do i=1,size(conv_type)
  if (conv_type(i)>=3) then
    scales%mb(i)=max(scales%mb(i) - 0.016, .005)
  end if
end do

!-------------------------------------------------------------------
! Alternative cloud base mass flux
!-------------------------------------------------------------------
scales%mb_new(:)=scales%mb(:)

!  Not currently doing this...
!     where (abs(wthvs) < tiny(wthvs(1))) ! Since wthvs can be zero from
!                                         ! BL in full model
!        scales%mb_new(:) = mb_a2*wstar_dn(:)
!     ELSEWHERE
!        ! added condition on dthetav_cb to avoid problems
!        scales%mb_new(:) = mb_a2*wthvs(:)/               &
!             merge(dthetav_cb(:), wthvs(:)/wstar_dn(:),  &
!             dthetav_cb(:) > wthvs(:)/wstar_dn(:))
!     end where

    !-------------------------------------------------------------------
    ! Ensemble vertical velocity at cloud base
    !-------------------------------------------------------------------
scales%wup2_cb(:) = wup_a1 * wstar_dn(:) * wstar_dn(:)

!-------------------------------------------------------------------
! Working note: Check this may wish to use rho and dz rather than
! delthuv
!-------------------------------------------------------------------
do i=1,size(ntml)
  scales%wstar_up(i) = (delthvu(i) * scales%mb(i) * g                          &
       / (theta(i,ntml(i))                                                     &
       * (1.0 + c_virtual * q_mix(i,ntml(i)))))**0.3333
  scales%zlcl(i) = z_rho(i,ntml(i)+1)
  scales%zcld(i) = z_rho(i,ntpar(i)+1) - z_rho(i,ntml(i)+1)
end do

scales%zcld_uv(:) = ztop_uv(:) - zlcl_uv(:)

!-------------------------------------------------------------------
! pre-calculate frequently used expressions
!-------------------------------------------------------------------
scales%wsc_o_mb(:) = scales%wstar_up(:)/scales%mb(:)
scales%mb_o_wsc(:) = scales%mb(:)/scales%wstar_up(:)
scales%root_mb_o_wsc(:) = sqrt(scales%mb_o_wsc(:))
scales%mb_new_o_wsc(:) = scales%mb_new(:)/scales%wstar_up(:)
scales%root_mb_new_o_wsc(:) = sqrt(scales%mb_new_o_wsc(:))
scales%wstar_up3(:)=scales%root_mb_new_o_wsc(:) *(scales%wstar_up(:)**3)
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_scales_warm

end module tcs_calc_scales
