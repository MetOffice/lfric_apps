! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_masked_field_adjust_type_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int32, int64, real64

! lfricinputs modules
USE lfricinp_regrid_weights_type_mod, ONLY: lfricinp_regrid_weights_type

! Shumlib modules
USE f_shum_field_mod, ONLY: shum_field_type

! LFRic modules
USE log_mod, ONLY: log_event, log_scratch_space, LOG_LEVEL_ERROR, &
  LOG_LEVEL_INFO
USE constants_mod, ONLY: imdi, rmdi

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_masked_field_adjust_type

! Type to contain information relating to masked fields adjustments/corrections

TYPE :: lfricinp_masked_field_adjust_type

  ! Number of destination points that do not have a full weight contribution
  ! from the source points of the same field mask type. ie some contributions
  ! come from points with different field mask type
  INTEGER(KIND=int32) :: num_adjusted_points

  ! Indices of adjusted points on the destination grid/mesh
  INTEGER(KIND=int32), ALLOCATABLE :: adjusted_dst_indices_1D(:)
  ! Map that links the which single source point has been selected to replace
  ! the adjusted destination data point
  INTEGER(KIND=int32), ALLOCATABLE :: adjusted_dst_to_src_map_2D(:,:)

  ! Flag to check whether the adjustment type has been initialised
  LOGICAL :: initialised = .FALSE.

  ! Destination mask - logical is true when point is VALID, logical false
  ! point should be ignored/masked out
  LOGICAL, ALLOCATABLE :: dst_mask_1D(:)

CONTAINS
  PROCEDURE :: find_adjusted_points_src_2d_dst_1d
  PROCEDURE :: apply_masked_adjustment_src_2d_dst_1d

END TYPE lfricinp_masked_field_adjust_type

CONTAINS
!---------------------------------------------------------
! Start of type bound procedures
!---------------------------------------------------------

SUBROUTINE find_adjusted_points_src_2d_dst_1d(self, src_mask, dst_mask, weights)
!
! This routine finds all masked field destination points that will require
! post regridding adjustment
!
! Argument(s)
!
CLASS(lfricinp_masked_field_adjust_type)       :: self
LOGICAL,                            INTENT(IN) :: src_mask(:,:)
LOGICAL,                            INTENT(IN) :: dst_mask(:)
TYPE(lfricinp_regrid_weights_type), INTENT(IN) :: weights

!
! Local variables
!
INTEGER(KIND=int32)              :: i, j, l, w
INTEGER(KIND=int32)              :: src_index1, src_index2, dst_index
INTEGER(KIND=int32), ALLOCATABLE :: dst_point_contrb_record(:)
REAL(KIND=real64)                :: weight_value
INTEGER(KIND=int32), PARAMETER   :: unchecked = 0, src_mask_contrb_only = 1,   &
                                    off_src_mask_contrb = 2
LOGICAL                          :: l_on_src_mask, l_on_dst_mask

ALLOCATE(self%dst_mask_1D(SIZE(dst_mask)))
self%dst_mask_1D(:) = dst_mask(:)

! Initialise arrays that records whether dst points had any or no
! contribution from on mask src points and the src point data
! to replace the dst point data
ALLOCATE(dst_point_contrb_record(SIZE(dst_mask)))
dst_point_contrb_record = unchecked

! Loop over remap matrix, considering only non-zero weight elements, to
! determine whether a dst point has contribution from any off mask src points
DO w = 1, weights%num_wgts
  DO l = 1, weights%num_links

    dst_index = weights%dst_address(l)
    l_on_dst_mask = dst_mask(dst_index)

    weight_value = ABS(weights%remap_matrix(w,l))
    IF (weight_value > 0.0 .AND. l_on_dst_mask) THEN

      src_index1 = weights%src_address_2d(l,1)
      src_index2 = weights%src_address_2d(l,2)
      l_on_src_mask= src_mask(src_index1, src_index2)

      ! Update records on whether the dst point has contributions from only
      ! masked src points, or some off mask source points.
      SELECT CASE (dst_point_contrb_record(dst_index))

        CASE (unchecked)
          IF (l_on_src_mask) THEN
            dst_point_contrb_record(dst_index) = src_mask_contrb_only
          ELSE
            dst_point_contrb_record(dst_index) = off_src_mask_contrb
          END IF

        CASE (src_mask_contrb_only)
          IF ( .NOT. l_on_src_mask) THEN
            dst_point_contrb_record(dst_index) = off_src_mask_contrb
          END IF

      END SELECT

    END IF
  END DO
END DO

! Set number of part resolved dst points.
self%num_adjusted_points = COUNT((dst_point_contrb_record ==                   &
                                       off_src_mask_contrb))
WRITE (log_scratch_space, '(A,I0)') "Number of adjusted points = ",            &
     self%num_adjusted_points
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

! Generate array of dst indices that requires post regridding masked adjustment
ALLOCATE(self%adjusted_dst_indices_1D(self%num_adjusted_points))
j = 0
DO i = 1, SIZE(dst_mask)
  IF (dst_point_contrb_record(i) == off_src_mask_contrb) THEN
    j = j + 1
    self%adjusted_dst_indices_1D(j) = i
  END IF
END DO

DEALLOCATE(dst_point_contrb_record)

END SUBROUTINE find_adjusted_points_src_2d_dst_1d

!---------------------------------------------------------

SUBROUTINE apply_masked_adjustment_src_2d_dst_1d(self, src, dst)
!
! Apply post regridding adjustments to masked field destination points
! Uses real arrays, could be overloaded for different types,
! precision and shapes
!
USE lfricinp_um_parameters_mod, ONLY: um_rmdi

!
! Argument(s)
!
REAL(KIND=real64), INTENT(IN) :: src(:,:)
REAL(KIND=real64), INTENT(IN OUT) :: dst(:)
CLASS(lfricinp_masked_field_adjust_type) :: self
!
! Local variables
!
INTEGER(KIND=int32) :: i

! Check if masked field adjust type has been initialised. If not
! report a warning
IF (self%initialised) THEN

  DO i = 1, self%num_adjusted_points
    dst(self%adjusted_dst_indices_1D(i)) =                                     &
                                     src(self%adjusted_dst_to_src_map_2D(i,1), &
                                         self%adjusted_dst_to_src_map_2D(i,2))
  END DO

  DO i = 1, SIZE(dst)
    IF (.NOT. self%dst_mask_1D(i)) THEN
      dst(i) = um_rmdi
    END IF
  END DO

ELSE

  log_scratch_space = 'Masked field adjustment type not initialised.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

END IF

END SUBROUTINE apply_masked_adjustment_src_2d_dst_1d

!---------------------------------------------------------

END MODULE lfricinp_masked_field_adjust_type_mod
