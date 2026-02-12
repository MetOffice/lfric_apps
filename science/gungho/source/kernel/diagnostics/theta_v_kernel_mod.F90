!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Computes theta_v, the virtual potential temperature.
module theta_v_kernel_mod

  use argument_mod,               only : arg_type, GH_SCALAR, &
                                         GH_FIELD, GH_REAL,   &
                                         GH_WRITE, GH_READ,   &
                                         CELL_COLUMN
  use constants_mod,              only : r_def, i_def
  use fs_continuity_mod,          only : Wtheta
  use kernel_mod,                 only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: theta_v_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                  &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)           &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: theta_v_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: theta_v_code

contains

!> @brief   Computes theta_v, the virtual potential temperature.
subroutine theta_v_code( nlayers,       &
                         theta_v,       &
                         theta,         &
                         mr_v,          &
                         recip_epsilon, &
                         ndf_wt,        &
                         undf_wt,       &
                         map_wt         )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: undf_wt, ndf_wt

  integer(kind=i_def), dimension(ndf_wt),  intent(in)    :: map_wt

  real(kind=r_def),    dimension(undf_wt), intent(inout) :: theta_v
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: theta
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: mr_v
  real(kind=r_def),                        intent(in)    :: recip_epsilon

  ! Internal variables
  integer(kind=i_def) :: df, min_col_index, max_col_index

  ! Find minimum and maximum DoF numberings for this column
  min_col_index = minval(map_wt)
  max_col_index = maxval(map_wt) + nlayers - 1

  ! Directly loop over all DoFs in the column
  do df = min_col_index, max_col_index
    theta_v(df) = theta(df) * (1.0_r_def + mr_v(df) * recip_epsilon) / (1.0_r_def + mr_v(df))
  end do

end subroutine theta_v_code

end module theta_v_kernel_mod
