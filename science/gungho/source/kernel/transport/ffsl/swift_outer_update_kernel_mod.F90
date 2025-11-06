!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Outer update for the SWIFT scheme
!> @details Computes the value of the new field based upon the
!!          transported field and dry mass in x- and y-sweeps
!!          and their increments.
!!

module swift_outer_update_kernel_mod

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
  type, public, extends(kernel_type) :: swift_outer_update_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                       &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_SPACE_1), & ! field_np1
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! field_x
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! dry_mass_x
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! field_y
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! dry_mass_y
         arg_type(GH_SCALAR, GH_REAL, GH_READ),               & ! dt
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! increment_x
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), & ! increment_y
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1)  & ! dry_mass_np1
         /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: swift_outer_update_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: swift_outer_update_code

contains
  !> @brief Outer update for the SWIFT scheme
  !> @param[inout] field_np1    Updated field to compute
  !> @param[in]    field_x      Field from the x-sweep
  !> @param[in]    dry_mass_x   Dry mass from the x-sweep
  !> @param[in]    field_y      Field from the y-sweep
  !> @param[in]    dry_mass_y   Dry mass from the y-sweep
  !> @param[in]    dt           Timestep
  !> @param[in]    increment_x  Increment from the x-sweep
  !> @param[in]    increment_y  Increment from the y-sweep
  !> @param[in]    dry_mass_np1 Final value of the transported dry mass
  subroutine swift_outer_update_code( field_np1,    &
                                      field_x,      &
                                      dry_mass_x,   &
                                      field_y,      &
                                      dry_mass_y,   &
                                      dt,           &
                                      increment_x,  &
                                      increment_y,  &
                                      dry_mass_np1 )

    implicit none

    real(r_tran), intent(inout) :: field_np1
    real(r_tran), intent(in)    :: field_x,     &
                                   dry_mass_x,  &
                                   field_y,     &
                                   dry_mass_y,  &
                                   dt,          &
                                   increment_x, &
                                   increment_y, &
                                   dry_mass_np1


    field_np1 = ( - 0.5_r_tran*dt*increment_x       &
                  - 0.5_r_tran*dt*increment_y       &
                  + 0.5_r_tran*field_x*dry_mass_x   &
                  + 0.5_r_tran*field_y*dry_mass_y ) &
              /dry_mass_np1

  end subroutine swift_outer_update_code

end module swift_outer_update_kernel_mod



