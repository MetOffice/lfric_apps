!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates Wbig: set 3D field to 1 where vertical velocity exceeds
!!        1 m s-1 and output.Produces the wbig field on each level.
!> @details Wbig is given a value of 1 if the value of w (vertical velocity)
!!          at that location is larger than 1 m s-1, otherwise it is set to
!!          zero. The diagnostic's primary use is in long, climate runs where
!!          it is averaged in time to help indicate whether the model has any
!!          stability problems.
module calc_wbig_kernel_mod

  use argument_mod,      only: arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_READ,           &
                               GH_READWRITE,      &
                               CELL_COLUMN
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: WTHETA
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_wbig_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                     &
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),      & !w_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA) & !wbig
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_wbig_code
  end type calc_wbig_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: calc_wbig_code

contains
!> @brief The identifies values of w_in_wth greater than 1 and sets them to 1
!>        in the wbig field.
!! @param[in] nlayers Integer the number of layers
!! @param[in] w_in_wth Real array, w component of u_physics
!! @param[in,out] wbig - the wbig field
!! @param[in] ndf_wth The number of degrees of freedom per cell for wth
!! @param[in] undf_wth The number of unique degrees of freedom for wth
!! @param[in] map_wth Integer array holding the dofmap for the cell at the
!>            base of the column for wth
subroutine calc_wbig_code(nlayers,                    &
                          w_in_wth,                  &
                          wbig,                       &
                          ndf_wth, undf_wth, map_wth  &
                          )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth

  real(kind=r_def), dimension(undf_wth), intent(in)    :: w_in_wth
  real(kind=r_def), dimension(undf_wth), intent(inout) :: wbig
  integer(kind=i_def), dimension(ndf_wth),  intent(in) :: map_wth

  ! Internal variables
  integer(kind=i_def) :: k

  do k = 0, nlayers

    if ( (w_in_wth(map_wth(1) + k)) > 1.0_r_def ) then

      wbig(map_wth(1) + k) = 1.0_r_def

    else

      wbig(map_wth(1) + k) = 0.0_r_def

    end if

  end do

end subroutine calc_wbig_code

end module calc_wbig_kernel_mod
