!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Filters the drag increments in W2 to prevent drag from accelerating
!>        the wind
module drag_incs_kernel_mod

  use argument_mod,  only : arg_type,                   &
                            GH_FIELD, GH_REAL, GH_READ, &
                            GH_WRITE, CELL_COLUMN
  use constants_mod, only : i_def, r_def
  use fs_continuity_mod, only: W2
  use kernel_mod,    only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  ! The type declaration for the kernel. Contains the metadata needed by the
  ! Psy layer.
  !
  type, public, extends(kernel_type) :: drag_incs_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                     &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2),        &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),        &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2)         &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: drag_incs_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: drag_incs_code

contains

!> @param[in]  nlayers Number of layers
!> @param[out] du_out  Filtered output increment
!> @param[in]  du_in   Unfiltered input increment
!> @param[in]  u_in    Wind to be incremented
!> @param[in]  ndf     Number of degrees of freedom per cell
!> @param[in]  undf    Total number of degrees of freedom
!> @param[in]  map     Dofmap for the cell at the base of the column
subroutine drag_incs_code(nlayers, du_out, du_in, u_in, &
                          ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, undf
  integer(kind=i_def), dimension(ndf), intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout) :: du_out
  real(kind=r_def), dimension(undf), intent(in) :: du_in, u_in

  ! Internal variables
  integer(kind=i_def) :: df, k

  ! Only loop over horizontal DoFs
  do df = 1, 4
    do k = 0, nlayers-1

      ! Only act on points not already done
      if (du_out(map(df)+k) == 0.0_r_def) then

        ! Ensure drag is always acting in opposite sense to the wind,
        ! i.e. u and du aren't the same sign
        if (sign(1.0_r_def,u_in(map(df)+k)) == sign(1.0_r_def,du_in(map(df)+k))) then
          ! Set output increment to 0
          du_out(map(df)+k) = 0.0_r_def
        ! Ensure drag does not decelerate the wind beyond 0,
        ! i.e. u+du is the same sign as u
        else if ( sign(1.0_r_def,u_in(map(df)+k)+du_in(map(df)+k)) &
             /= sign(1.0_r_def,u_in(map(df)+k)) ) then
          ! Set output increment to negative of input wind
          du_out(map(df)+k) = -u_in(map(df)+k)
        else
          ! Output increment is input increment
          du_out(map(df)+k) = du_in(map(df)+k)
        end if
      end if

    end do
  end do

end subroutine drag_incs_code

end module drag_incs_kernel_mod
