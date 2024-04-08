! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_apply_masked_field_adjustments_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64, real64

IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_apply_masked_field_adjustments


CONTAINS

SUBROUTINE um2lfric_apply_masked_field_adjustments(stashcode, src, dst)

USE lfricinp_stashmaster_mod, ONLY: land_compressed, grid,                     &
  get_stashmaster_item, p_points_values_over_sea
USE um2lfric_masked_field_adjustments_mod, ONLY: land_field_adjustments,       &
                                                 maritime_field_adjustments

IMPLICIT NONE

! Arguments
INTEGER(KIND=int64), INTENT(IN) :: stashcode
REAL(KIND=real64), INTENT(IN) :: src(:,:)
REAL(KIND=real64), INTENT(IN OUT) :: dst(:)

! Get grid code from stashmaster
SELECT CASE (get_stashmaster_item(stashcode, grid))

  CASE(land_compressed)
    CALL land_field_adjustments%apply_masked_adjustment_src_2d_dst_1d(src,     &
                                                                      dst)
  CASE(p_points_values_over_sea)
    CALL maritime_field_adjustments%apply_masked_adjustment_src_2d_dst_1d(src, &
                                                                      dst)

END SELECT

END SUBROUTINE um2lfric_apply_masked_field_adjustments

END MODULE um2lfric_apply_masked_field_adjustments_mod

