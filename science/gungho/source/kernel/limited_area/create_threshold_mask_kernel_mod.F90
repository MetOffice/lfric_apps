!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Sets output field to flag_value where input field is greater than threshold
!>        at the base of the column.
!> @details Uses an input scalar to determine threshold.
!>        Note that this kernel can be used with continuous spaces and will
!>        set dofs that have ANY neighbouring cell meeting the threshold criterion.
module create_threshold_mask_kernel_mod

  use argument_mod,              only : arg_type, func_type, &
                                        mesh_data_type,      &
                                        GH_SCALAR, GH_FIELD, &
                                        GH_READ, GH_WRITE,   &
                                        GH_REAL, GH_BASIS,   &
                                        CELL_COLUMN,         &
                                        ANY_SPACE_1
  use fs_continuity_mod,         only : W3
  use constants_mod,             only : r_def, i_def, l_def
  use kernel_mod,                only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: create_threshold_mask_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_SPACE_1),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3),                         &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: create_threshold_mask_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: create_threshold_mask_code

contains

!> @param[in] nlayers    Number of layers
!> @param[in,out] threshold_mask The output field
!> @param[in] reference_field    The input field
!> @param[in] threshold  Threshold value determining if it's on or off
!> @param[in] flag_value real value (nominally 0 or 1) which sets things on or off
!> @param[in] ndf_out    Number of degrees of freedom for threshold_mask
!> @param[in] undf_out   Total number of degrees of freedom for threshold_mask
!> @param[in] map_out    Dofmap for the cell at the base of the column for threshold_mask
!> @param[in] ndf_in     Number of degrees of freedom for reference_field
!> @param[in] undf_in    Total number of degrees of freedom for reference_field
!> @param[in] map_in     Dofmap for the cell at the base of the column for reference_field
subroutine create_threshold_mask_code( nlayers,         &
                                       threshold_mask,  &
                                       reference_field, &
                                       threshold,       &
                                       flag_value,      &
                                       ndf_out,         &
                                       undf_out,        &
                                       map_out,         &
                                       ndf_in,          &
                                       undf_in,         &
                                       map_in)

  implicit none

  ! Arguments
  integer(kind=i_def),       intent(in) :: nlayers
  integer(kind=i_def),       intent(in) :: ndf_out, undf_out
  integer(kind=i_def),       intent(in) :: ndf_in, undf_in
  real(kind=r_def), dimension(undf_out), intent(inout)  :: threshold_mask
  real(kind=r_def), dimension(undf_in),  intent(in)     :: reference_field
  real(kind=r_def),          intent(in) :: threshold
  real(kind=r_def),          intent(in) :: flag_value
  integer(kind=i_def), dimension(ndf_out), intent(in) :: map_out
  integer(kind=i_def), dimension(ndf_in),  intent(in) :: map_in

  ! Internal variables
  integer(kind=i_def) :: k, df

  ! Note that we only use the base of the column for the input field, so
  ! a 2D field can be used as input if required.
  if (reference_field(map_in(1)) > threshold)then
    do k=0,nlayers-1
      do df=1,ndf_out
        threshold_mask(map_out(df)+k) = flag_value
      end do
    end do
  end if


end subroutine create_threshold_mask_code

end module create_threshold_mask_kernel_mod
