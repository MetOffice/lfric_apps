! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!> @brief     Module for masked field adjustments TYPE
!> @details   Holds data and procedures for making adjustments to points
!!            that are contributed to by multiple field mask types
module lfricinp_masked_field_adjust_type_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int32, int64, real64

! lfricinputs modules
use lfricinp_regrid_weights_type_mod, only: lfricinp_regrid_weights_type

! Shumlib modules
use f_shum_field_mod, only: shum_field_type

! LFRic modules
use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_ERROR, &
                   LOG_LEVEL_INFO
use constants_mod, only: imdi, rmdi

implicit none

private

public :: lfricinp_masked_field_adjust_type

! Type to contain information relating to masked fields adjustments/corrections

type :: lfricinp_masked_field_adjust_type

  ! Number of destination points that do not have a full weight contribution
  ! from the source points of the same field mask type. ie some contributions
  ! come from points with different field mask type
  integer(kind=int32) :: num_adjusted_points

  ! Indices of adjusted points on the destination grid/mesh
  integer(kind=int32), allocatable :: adjusted_dst_indices_1D(:)
  ! Map that links the which single source point has been selected to replace
  ! the adjusted destination data point
  integer(kind=int32), allocatable :: adjusted_dst_to_src_map_2D(:,:)

  ! Flag to check whether the adjustment type has been initialised
  logical :: initialised = .false.

  ! Destination mask - logical is true when point is VALID, logical false
  ! point should be ignored/masked out
  logical, allocatable :: dst_mask_1D(:)

contains
  procedure :: find_adjusted_points_src_2d_dst_1d
  procedure :: apply_masked_adjustment_src_2d_dst_1d

end type lfricinp_masked_field_adjust_type

contains
!---------------------------------------------------------
! Start of type bound procedures
!---------------------------------------------------------
!> @brief   Finds all masked field destination points that will require
!!          post regridding adjustment
!> @param[in] src_mask  Mask that applies to the source grid
!> @param[in] dst_mask  Mask that applies to the destination grid
!> @param[in] weights   Regridding weights for masked data
subroutine find_adjusted_points_src_2d_dst_1d(self, src_mask, dst_mask, weights)
!
! Argument(s)
!
class(lfricinp_masked_field_adjust_type)       :: self
logical,                            intent(in) :: src_mask(:,:)
logical,                            intent(in) :: dst_mask(:)
type(lfricinp_regrid_weights_type), intent(in) :: weights

!
! Local variables
!
integer(kind=int32)              :: i, j, l, w
integer(kind=int32)              :: src_index1, src_index2, dst_index
integer(kind=int32), allocatable :: dst_point_contrb_record(:)
real(kind=real64)                :: weight_value
integer(kind=int32), parameter   :: unchecked = 0, src_mask_contrb_only = 1,   &
                                    off_src_mask_contrb = 2
logical                          :: l_on_src_mask, l_on_dst_mask

allocate(self%dst_mask_1D(size(dst_mask)))
self%dst_mask_1D(:) = dst_mask(:)

! Initialise arrays that records whether dst points had any or no
! contribution from on mask src points and the src point data
! to replace the dst point data
allocate(dst_point_contrb_record(size(dst_mask)))
dst_point_contrb_record = unchecked

! Loop over remap matrix, considering only non-zero weight elements, to
! determine whether a dst point has contribution from any off mask src points
do w = 1, weights%num_wgts
  do l = 1, weights%num_links

    dst_index = weights%dst_address(l)
    l_on_dst_mask = dst_mask(dst_index)

    weight_value = abs(weights%remap_matrix(w,l))
    if (weight_value > 0.0_real64 .and. l_on_dst_mask) then

      src_index1 = weights%src_address_2d(l,1)
      src_index2 = weights%src_address_2d(l,2)
      l_on_src_mask= src_mask(src_index1, src_index2)

      ! Update records on whether the dst point has contributions from only
      ! masked src points, or some off mask source points.
      select case (dst_point_contrb_record(dst_index))

        case (unchecked)
          if (l_on_src_mask) then
            dst_point_contrb_record(dst_index) = src_mask_contrb_only
          else
            dst_point_contrb_record(dst_index) = off_src_mask_contrb
          end if

        case (src_mask_contrb_only)
          if ( .not. l_on_src_mask) then
            dst_point_contrb_record(dst_index) = off_src_mask_contrb
          end if

      end select

    end if
  end do
end do

! Set number of part resolved dst points.
self%num_adjusted_points = count((dst_point_contrb_record ==                   &
                                       off_src_mask_contrb))
write (log_scratch_space, '(A,I0)') "Number of adjusted points = ",            &
     self%num_adjusted_points
call log_event(log_scratch_space, LOG_LEVEL_INFO)

! Generate array of dst indices that requires post regridding masked adjustment
allocate(self%adjusted_dst_indices_1D(self%num_adjusted_points))
j = 0
do i = 1, size(dst_mask)
  if (dst_point_contrb_record(i) == off_src_mask_contrb) then
    j = j + 1
    self%adjusted_dst_indices_1D(j) = i
  end if
end do

deallocate(dst_point_contrb_record)

end subroutine find_adjusted_points_src_2d_dst_1d

!---------------------------------------------------------
!> @brief Applies post regridding adjustments to masked field destination points
!> @param[in]     src   2d source data array
!> @param[inout]  dst   1d destination array
subroutine apply_masked_adjustment_src_2d_dst_1d(self, src, dst)
!
! Uses real arrays, could be overloaded for different types,
! precision and shapes
!
use lfricinp_um_parameters_mod, only: um_rmdi

!
! Argument(s)
!
real(kind=real64), intent(in) :: src(:,:)
real(kind=real64), intent(in out) :: dst(:)
class(lfricinp_masked_field_adjust_type) :: self
!
! Local variables
!
integer(kind=int32) :: i

! Check if masked field adjust type has been initialised. If not
! report a warning
if (self%initialised) then

  do i = 1, self%num_adjusted_points
    dst(self%adjusted_dst_indices_1D(i)) =                                     &
                                     src(self%adjusted_dst_to_src_map_2D(i,1), &
                                         self%adjusted_dst_to_src_map_2D(i,2))
  end do

  do i = 1, size(dst)
    if (.not. self%dst_mask_1D(i)) then
      dst(i) = um_rmdi
    end if
  end do

else

  log_scratch_space = 'Masked field adjustment type not initialised.'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)

end if

end subroutine apply_masked_adjustment_src_2d_dst_1d

!---------------------------------------------------------

end module lfricinp_masked_field_adjust_type_mod
