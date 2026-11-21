!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of the combining of two multidata fields with n1 & n2 data points
!!        into a single field with n = n1 + n2 data points

module adj_sci_combine_multidata_field_kernel_mod

  use, intrinsic :: iso_fortran_env, only: real32, real64

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_REAL,         &
                                GH_SCALAR, GH_INTEGER,     &
                                GH_LOGICAL,                &
                                GH_INC, GH_READ,           &
                                GH_READINC,                &
                                CELL_COLUMN,               &
                                ANY_SPACE_1,               &
                                ANY_SPACE_2,               &
                                ANY_SPACE_3
  use constants_mod,     only : i_def, l_def
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: adj_combine_multidata_field_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                             &
         arg_type(GH_FIELD,  GH_REAL,    GH_READINC, ANY_SPACE_1),  &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                  &
         arg_type(GH_FIELD,  GH_REAL,    GH_INC,     ANY_SPACE_2),  &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                  &
         arg_type(GH_FIELD,  GH_REAL,    GH_INC,     ANY_SPACE_3),  &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                  &
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                   &
         /)
    integer :: operates_on = CELL_COLUMN
  end type
  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_combine_multidata_field_code

  ! Generic interface for real32 and real64 types
  interface adj_combine_multidata_field_code
    module procedure  &
      adj_combine_multidata_field_code_real32, &
      adj_combine_multidata_field_code_real64
  end interface

contains

!> @brief Adjoint of combination of two multidata fields into a third field
!! @param[in]     nlayers     Number of layers
!! @param[in,out] field_out   Output multidata field
!! @param[in]     n           Number of data points in field_out
!! @param[in]     field1_in   First input multidata field
!! @param[in]     n1          Number of data points in field1_in
!! @param[in]     field2_in   Second input multidata field
!! @param[in]     n2          Number of data points in field2_in
!! @param[in]     ndata_first Data layout (data- or layer-first)
!! @param[in]     ndf         Number of dofs per cell for field_out
!! @param[in]     undf        Total number of dofs per cell for field_out
!! @param[in]     map         Cell dofmap for field_out
!! @param[in]     ndf_1       Number of dofs per cell for field1_in
!! @param[in]     undf_1      Total number of dofs per cell for field1_in
!! @param[in]     map_1       Cell dofmap for field1_in
!! @param[in]     ndf_2       Number of dofs per cell for field2_in
!! @param[in]     undf_2      Total number of dofs per cell for field2_in
!! @param[in]     map_2       Cell dofmap for field2_in

! REAL32 PRECISION
! ==================
subroutine adj_combine_multidata_field_code_real32( nlayers,              &
                                                    field_out, n,         &
                                                    field1_in, n1,        &
                                                    field2_in, n2,        &
                                                    ndata_first,          &
                                                    ndf, undf, map,       &
                                                    ndf_1, undf_1, map_1, &
                                                    ndf_2, undf_2, map_2 )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: ndf, ndf_1, ndf_2
  integer(kind=i_def),                     intent(in)    :: undf, undf_1, undf_2
  integer(kind=i_def),                     intent(in)    :: n, n1, n2
  integer(kind=i_def), dimension(ndf),     intent(in)    :: map
  integer(kind=i_def), dimension(ndf_1),   intent(in)    :: map_1
  integer(kind=i_def), dimension(ndf_2),   intent(in)    :: map_2
  real(kind=real32),   dimension(undf),    intent(inout) :: field_out
  real(kind=real32),   dimension(undf_1),  intent(inout) :: field1_in
  real(kind=real32),   dimension(undf_2),  intent(inout) :: field2_in
  logical(kind=l_def),                     intent(in)    :: ndata_first


  ! Internal variables
  integer(kind=i_def) :: df, k, ij

  if ( ndata_first ) then
    do k = nlayers-1, 0, -1
      ij = map(1) + k*n + n1
      do df = n2-1, 0, -1
        field2_in(map_2(1) + k*n2 + df) = field2_in(map_2(1) + k*n2 + df) &
                                        + field_out(ij + df)
        field_out(ij + df) = 0.0_real32
      end do

      ij = map(1) + k*n
      do df = n1-1, 0, -1
        field1_in(map_1(1) + k*n1 + df) = field1_in(map_1(1) + k*n1 + df) &
                                        + field_out(ij + df)
        field_out(ij + df) = 0.0_real32
      end do
    end do
  else
    do df = n2, 1, -1
      ij = map(1) + n1*nlayers + (df-1)*nlayers
      do k = nlayers-1, 0, -1
        field2_in(map_2(1) + (df-1)*nlayers + k) = field2_in(map_2(1) + (df-1)*nlayers + k) &
                                                 + field_out(ij + k)
        field_out(ij + k) = 0.0_real32
      end do
    end do

    do df = n1, 1, -1
      ij = map(1) + (df-1)*nlayers
      do k = nlayers-1, 0, -1
        field1_in(map_1(1) + (df-1)*nlayers + k) = field1_in(map_1(1) + (df-1)*nlayers + k) &
                                                 + field_out(ij + k)
        field_out(ij + k) = 0.0_real32
      end do
    end do
  end if

end subroutine adj_combine_multidata_field_code_real32

! REAL64 PRECISION
! ==================
subroutine adj_combine_multidata_field_code_real64( nlayers,              &
                                                    field_out, n,         &
                                                    field1_in, n1,        &
                                                    field2_in, n2,        &
                                                    ndata_first,          &
                                                    ndf, undf, map,       &
                                                    ndf_1, undf_1, map_1, &
                                                    ndf_2, undf_2, map_2 )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: ndf, ndf_1, ndf_2
  integer(kind=i_def),                     intent(in)    :: undf, undf_1, undf_2
  integer(kind=i_def),                     intent(in)    :: n, n1, n2
  integer(kind=i_def), dimension(ndf),     intent(in)    :: map
  integer(kind=i_def), dimension(ndf_1),   intent(in)    :: map_1
  integer(kind=i_def), dimension(ndf_2),   intent(in)    :: map_2
  real(kind=real64),   dimension(undf),    intent(inout) :: field_out
  real(kind=real64),   dimension(undf_1),  intent(inout) :: field1_in
  real(kind=real64),   dimension(undf_2),  intent(inout) :: field2_in
  logical(kind=l_def),                     intent(in)    :: ndata_first


  ! Internal variables
  integer(kind=i_def) :: df, k, ij

  if ( ndata_first ) then
    do k = nlayers-1, 0, -1
      ij = map(1) + k*n + n1
      do df = n2-1, 0, -1
        field2_in(map_2(1) + k*n2 + df) = field2_in(map_2(1) + k*n2 + df) &
                                        + field_out(ij + df)
        field_out(ij + df) = 0.0_real64
      end do

      ij = map(1) + k*n
      do df = n1-1, 0, -1
        field1_in(map_1(1) + k*n1 + df) = field1_in(map_1(1) + k*n1 + df) &
                                        + field_out(ij + df)
        field_out(ij + df) = 0.0_real64
      end do
    end do
  else
    do df = n2, 1, -1
      ij = map(1) + n1*nlayers + (df-1)*nlayers
      do k = nlayers-1, 0, -1
        field2_in(map_2(1) + (df-1)*nlayers + k) = field2_in(map_2(1) + (df-1)*nlayers + k) &
                                                 + field_out(ij + k)
        field_out(ij + k) = 0.0_real64
      end do
    end do

    do df = n1, 1, -1
      ij = map(1) + (df-1)*nlayers
      do k = nlayers-1, 0, -1
        field1_in(map_1(1) + (df-1)*nlayers + k) = field1_in(map_1(1) + (df-1)*nlayers + k) &
                                                 + field_out(ij + k)
        field_out(ij + k) = 0.0_real64
      end do
    end do
  end if

end subroutine adj_combine_multidata_field_code_real64

end module adj_sci_combine_multidata_field_kernel_mod
