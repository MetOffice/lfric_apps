!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Module with functions to help initialise arbitrary/random masks
!! for use in limited area models.

module set_mask_by_threshold_kernel_mod

  use argument_mod,                          only : ANY_SPACE_1,           &
                                                    arg_type, CELL_COLUMN, &
                                                    GH_FIELD, GH_SCALAR,   &
                                                    GH_REAL, GH_READINC,   &
                                                    GH_READ
  use kernel_mod,                            only : kernel_type
  use constants_mod,                         only : i_def, r_single, r_double

  implicit none

  type, extends(kernel_type) :: set_mask_by_threshold_kernel_type
     type(arg_type), dimension(2) :: meta_args = (/                      &
       arg_type(GH_FIELD,  GH_REAL, GH_READINC, ANY_SPACE_1),            & ! mask
       arg_type(GH_SCALAR, GH_REAL, GH_READ)                             & ! threshold
     /)
     integer :: operates_on = CELL_COLUMN
  end type set_mask_by_threshold_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: set_mask_by_threshold_kernel_code

  interface set_mask_by_threshold_kernel_code
    module procedure &
      set_mask_by_threshold_kernel_code_r_single, &
      set_mask_by_threshold_kernel_code_r_double
  end interface

contains
  !> @brief  Generates r_single mask using threshold.
  !>         For values of the input mask above the threshold, assign 1,
  !>         else assign 0.
  !> @param[in]    nlayers     No. of layers.
  !> @param[inout] mask        Mask to be assigned.
  !> @param[in]    threshold   Threshold to determine mask assignment.
  !> @param[in]    ndf_mask    No. degrees of freedom for mask.
  !> @param[in]    undf_mask   Unique no. degrees of freedom for mask.
  !> @param[in]    map_mask    Map for mask.
  subroutine set_mask_by_threshold_kernel_code_r_single( nlayers,       &
                                                         mask,          &
                                                         threshold,     &
                                                         ndf_mask,      &
                                                         undf_mask,     &
                                                         map_mask )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_mask
    integer(kind=i_def), intent(in) :: undf_mask
    real(kind=r_single), intent(in) :: threshold
    real(kind=r_single), intent(inout), dimension(undf_mask) :: mask
    integer(kind=i_def), intent(in),    dimension(ndf_mask)  :: map_mask

    ! Internal variables
    integer(kind=i_def) :: k, df

    do k = 0, nlayers-1
      do df = 1, ndf_mask
        if (mask(map_mask(df) + k) <= threshold) then
          mask(map_mask(df) + k) = 0.0_r_single
        else
          mask(map_mask(df) + k) = 1.0_r_single
        end if
      end do
    end do

  end subroutine set_mask_by_threshold_kernel_code_r_single

  !> @brief  Generates r_double mask using threshold.
  !>         For values of the input mask above the threshold, assign 1,
  !>         else assign 0.
  !> @param[in]    nlayers     No. of layers.
  !> @param[inout] mask        Mask to be assigned.
  !> @param[in]    threshold   Threshold to determine mask assignment.
  !> @param[in]    ndf_mask    No. degrees of freedom for mask.
  !> @param[in]    undf_mask   Unique no. degrees of freedom for mask.
  !> @param[in]    map_mask    Map for mask.
  subroutine set_mask_by_threshold_kernel_code_r_double( nlayers,       &
                                                         mask,          &
                                                         threshold,     &
                                                         ndf_mask,      &
                                                         undf_mask,     &
                                                         map_mask )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_mask
    integer(kind=i_def), intent(in) :: undf_mask
    real(kind=r_double), intent(in) :: threshold
    real(kind=r_double), intent(inout), dimension(undf_mask) :: mask
    integer(kind=i_def), intent(in),    dimension(ndf_mask)  :: map_mask

    ! Internal variables
    integer(kind=i_def) :: k, df

    do k = 0, nlayers-1
      do df = 1, ndf_mask
        if (mask(map_mask(df) + k) <= threshold) then
          mask(map_mask(df) + k) = 0.0_r_double
        else
          mask(map_mask(df) + k) = 1.0_r_double
        end if
      end do
    end do

  end subroutine set_mask_by_threshold_kernel_code_r_double

end module set_mask_by_threshold_kernel_mod
