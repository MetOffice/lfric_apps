!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Inner update for the SWIFT scheme
!> @details Computes the value of the new field based upon the
!!          transported field and dry mass in either the x- or y-direction
!!
module swift_inner_update_kernel_mod

  use argument_mod,  only : arg_type,          &
                            GH_FIELD, GH_REAL, &
                            GH_WRITE, GH_READ, &
                            GH_SCALAR,         &
                            DOF, ANY_SPACE_1
  use constants_mod, only : r_tran, i_def
  use kernel_mod,    only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: swift_inner_update_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                       &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_SPACE_1), & ! field_np1
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! field_n
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! dry_mass_np1
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! dry_mass_n
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! increment_x
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                & ! dt
         /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: swift_inner_update_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: swift_inner_update_code

contains
  !> @brief Inner update for the SWIFT scheme
  !> @param[inout] field_np1    Updated field to compute
  !> @param[in]    field_n      Initial field
  !> @param[in]    dry_mass_np1 Final value of the transported dry mass
  !> @param[in]    dry_mass_n   Initial mass
  !> @param[in]    increment    Increment from the sweep
  !> @param[in]    dt           Timestep
  subroutine swift_inner_update_code(field_np1,    &
                                     field_n,      &
                                     dry_mass_np1, &
                                     dry_mass_n,   &
                                     increment,    &
                                     dt)

    implicit none

    real(r_tran), intent(inout) :: field_np1
    real(r_tran), intent(in)    :: field_n,      &
                                   dry_mass_np1, &
                                   dry_mass_n,   &
                                   increment,    &
                                   dt

    field_np1 = (field_n*dry_mass_n - dt*increment)/dry_mass_np1

  end subroutine swift_inner_update_code

end module swift_inner_update_kernel_mod



