! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module um2lfric_apply_masked_field_adjustments_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64, real64

implicit none

private

public :: um2lfric_apply_masked_field_adjustments


contains

subroutine um2lfric_apply_masked_field_adjustments(stashcode, src, dst)

use lfricinp_stashmaster_mod,              only: land_compressed,      &
                                                 grid,                 &
                                                 get_stashmaster_item, &
                                                 p_points_values_over_sea
use um2lfric_masked_field_adjustments_mod, only: land_field_adjustments, &
                                                 maritime_field_adjustments

implicit none

! Arguments
integer(kind=int64), intent(in) :: stashcode
real(kind=real64), intent(in) :: src(:,:)
real(kind=real64), intent(in out) :: dst(:)

! Identify fields which need to be processed as sea-point only fields.
integer(kind=int64), parameter :: sea_pt_fld(3)= [195,440,441]

! We can't identify sea-point based fields purely by STASH grid code
! because things like multi-category seaice fields are just defined
! with a standard STASH grid code of 1. So we have to identify them by
! individual grid code!

if (any(stashcode == sea_pt_fld)) then

  call maritime_field_adjustments%apply_masked_adjustment_src_2d_dst_1d(src, &
                                                                     dst)
else

  ! Get grid code from stashmaster
  select case (get_stashmaster_item(stashcode, grid))

    case(land_compressed)

      call land_field_adjustments%apply_masked_adjustment_src_2d_dst_1d(src,     &
                                                                      dst)
    case(p_points_values_over_sea)

      call maritime_field_adjustments%apply_masked_adjustment_src_2d_dst_1d(src, &
                                                                     dst)
   end select
end if

end subroutine um2lfric_apply_masked_field_adjustments

end module um2lfric_apply_masked_field_adjustments_mod
