!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the advective increment for the FFSL scheme
!> @details Compute the advective increment based upon field value
!!          and the transported unity field
module ffsl_advective_increment_kernel_mod

  use argument_mod,  only : arg_type,              &
                            GH_FIELD, GH_REAL,     &
                            GH_READWRITE, GH_READ, &
                            GH_SCALAR,             &
                            DOF, ANY_SPACE_1
  use constants_mod, only : r_tran
  use kernel_mod,    only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: ffsl_advective_increment_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                           &
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_SPACE_1), & ! increment
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_SPACE_1), & ! field
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                   & ! one_over_dt
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_SPACE_1)  & ! adv_one
         /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: ffsl_advective_increment_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: ffsl_advective_increment_code

contains

  !> @brief Compute the advective increment for the FFSL scheme
  !> @param[inout] Increment   Updated increment to compute
  !> @param[in]    field       Transported field
  !> @param[in]    one_over_dt Inverse timestep
  !> @param[in]    adv_one     Transported unity field
  subroutine ffsl_advective_increment_code(increment,   &
                                           field,       &
                                           one_over_dt, &
                                           adv_one)

    implicit none

    real(r_tran), intent(inout) :: increment
    real(r_tran), intent(in)    :: adv_one,      &
                                   field,        &
                                   one_over_dt

    increment = (increment + (adv_one - 1.0_r_tran)*field*one_over_dt)/adv_one

  end subroutine ffsl_advective_increment_code

end module ffsl_advective_increment_kernel_mod



