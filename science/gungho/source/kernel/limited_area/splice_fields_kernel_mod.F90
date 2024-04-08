!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Where the input mask above a threshold, replace the values
!>        of the output field with the input field.
!> @details Uses an input scalar to determine threshold
module splice_fields_kernel_mod

  use argument_mod,              only : arg_type, func_type, &
                                        mesh_data_type,      &
                                        GH_SCALAR, GH_FIELD, &
                                        GH_READ, GH_WRITE,   &
                                        GH_REAL, GH_BASIS,   &
                                        CELL_COLUMN,         &
                                        ANY_SPACE_1
  use fs_continuity_mod,         only : W3
  use constants_mod,             only : r_def, i_def, l_def, r_single, r_double
  use kernel_mod,                only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: splice_fields_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                       &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_SPACE_1),& ! field_A
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_SPACE_1), & ! field_B
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3),          & ! onion_layers
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)               & ! threshold
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: splice_fields_code

  ! Generic interface for real32 and real64 types
  interface splice_fields_code
    module procedure  &
      splice_fields_code_r_single, &
      splice_fields_code_r_double
  end interface

contains

!> @param[in] nlayers      Number of layers
!> @param[in,out] field_A  Field to be overwritten where threshold exceded
!> @param[in] field_B      Field to use in the overwritten region
!> @param[in] threshold    Threshold value
!> @param[in] ndf          Number of degrees of freedom for field_A and field_B
!> @param[in] undf         Total number of degrees of freedom for field_A and field_B
!> @param[in] map          Dofmap for the cell at the base of the column for field_A and field_B
!> @param[in] ndf_w3       Number of degrees of freedom for onion_layers
!> @param[in] undf_w3      Total number of degrees of freedom for onion_layers
!> @param[in] map_w3       Dofmap for the cell at the base of the column for onion_layers


! R_SINGLE PRECISION
! ==================
subroutine splice_fields_code_r_single( nlayers,      &
                                        field_A,      &
                                        field_B,      &
                                        onion_layers, &
                                        threshold,    &
                                        ndf,          &
                                        undf,         &
                                        map,          &
                                        ndf_w3,       &
                                        undf_w3,      &
                                        map_w3)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: nlayers
  integer(kind=i_def),                    intent(in) :: ndf, undf
  integer(kind=i_def),                    intent(in) :: ndf_w3, undf_w3
  real(kind=r_single), dimension(undf),   intent(inout) :: field_A
  real(kind=r_single), dimension(undf),   intent(in) :: field_B
  real(kind=r_def),    dimension(undf),   intent(in) :: onion_layers
  real(kind=r_def),                       intent(in) :: threshold
  integer(kind=i_def), dimension(ndf),    intent(in) :: map
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  ! Internal variables
  integer(kind=i_def)                    :: k, df

  if (onion_layers(map_w3(1)) > threshold) then
    do k=0,nlayers-1
      do df=1,ndf
        field_A(map(df)+k) = field_B(map(df)+k)
      end do
    end do
  end if

end subroutine splice_fields_code_r_single

! R_DOUBLE PRECISION
! ==================
subroutine splice_fields_code_r_double( nlayers,      &
                                        field_A,      &
                                        field_B,      &
                                        onion_layers, &
                                        threshold,    &
                                        ndf,          &
                                        undf,         &
                                        map,          &
                                        ndf_w3,       &
                                        undf_w3,      &
                                        map_w3)

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: nlayers
  integer(kind=i_def),                     intent(in) :: ndf, undf
  integer(kind=i_def),                     intent(in) :: ndf_w3, undf_w3
  real(kind=r_double), dimension(undf), intent(inout) :: field_A
  real(kind=r_double), dimension(undf),    intent(in) :: field_B
  real(kind=r_def),    dimension(undf),    intent(in) :: onion_layers
  real(kind=r_def),                        intent(in) :: threshold
  integer(kind=i_def), dimension(ndf),     intent(in) :: map
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

  ! Internal variables
  integer(kind=i_def)                    :: k, df

  if (onion_layers(map_w3(1)) > threshold) then
    do k=0,nlayers-1
      do df=1,ndf
        field_A(map(df)+k) = field_B(map(df)+k)
      end do
    end do
  end if

end subroutine splice_fields_code_r_double

end module splice_fields_kernel_mod
