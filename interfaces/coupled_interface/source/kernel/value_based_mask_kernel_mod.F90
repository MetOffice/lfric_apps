!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the
! terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Sets mask to 1 in a particular value range
!>
module value_based_mask_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_REAL, GH_READ, &
                               CELL_COLUMN, GH_SCALAR, GH_WRITE,     &
                               ANY_DISCONTINUOUS_SPACE_1,            &
                               ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,     only: r_def, i_def, degrees_to_radians,     &
                               radians_to_degrees
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed
  !by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: value_based_mask_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/&
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),          &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)           &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: value_based_mask
  end type value_based_mask_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: value_based_mask

contains

  !> @details       Sets mask to 1.0 where values field is in specified range
  !> @param[in]     nlayers        number of levels
  !> @param[in,out] mask           1.0 within specified range
  !> @param[in]     values         field containing values to check
  !> @param[in]     min_value      minimum value where mask=1.0
  !> @param[in]     max_value        maximum value where mask=1.0
  !> @param[in]     ndf_2d         The number of dofs per cell for 2d field
  !> @param[in]     undf_2d        The number of unique dofs for 2d field
  !> @param[in]     map_2d         array holding the dofmap for 2d field
  subroutine value_based_mask(nlayers, mask, values,     &
                              min_value, max_value,     &
                              ndf_2d, undf_2d, map_2d         &
                              )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers

    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d

    real(kind=r_def), dimension(undf_2d), intent(inout) :: mask
    real(kind=r_def), dimension(undf_2d), intent(in) :: values
    integer(kind=i_def), dimension(ndf_2d),  intent(in) :: map_2d
    real(kind=r_def), intent(in) :: min_value, max_value


    if ( (values(map_2d(1)) .ge. min_value) .and.                    &
         (values(map_2d(1)) .le. max_value)) then
      mask(map_2d(1)) = 1.0_r_def
    else
      mask(map_2d(1)) = 0.0_r_def
    end if

  end subroutine value_based_mask

end module value_based_mask_kernel_mod

