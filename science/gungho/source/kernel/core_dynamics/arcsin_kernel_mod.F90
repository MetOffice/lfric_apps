!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Applies the arcsin function to a field

module arcsin_kernel_mod

  use argument_mod,  only : arg_type,          &
                            GH_FIELD, GH_REAL, &
                            GH_WRITE, GH_READ, &
                            DOF, ANY_SPACE_1
  use constants_mod, only : r_def
  use kernel_mod,    only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: arcsin_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                       &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_SPACE_1), & ! field_out
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1)  & ! field_in
         /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: arcsin_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: arcsin_code

contains
  !> @brief Takes arcsin function
  !> @param[in,out] field_out   Output field
  !> @param[in]     field_in    Input field
  subroutine arcsin_code(field_out, field_in)

    implicit none

    real(kind=r_def), intent(inout) :: field_out
    real(kind=r_def), intent(in)    :: field_in

    field_out = asin(field_in)

  end subroutine arcsin_code

end module arcsin_kernel_mod



