!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Module with functions to initialise operators with average value of two fields.

module setop_average_two_fields_kernel_mod

  use argument_mod,                          only : ANY_SPACE_1, ANY_SPACE_2, &
                                                    arg_type, cell_column,    &
                                                    GH_OPERATOR, GH_FIELD,    &
                                                    GH_REAL, GH_WRITE,        &
                                                    GH_READ
  use kernel_mod,                            only : kernel_type
  use constants_mod,                         only : i_def, r_single, r_double

  implicit none

  type, extends(kernel_type) :: setop_average_two_fields_kernel_type
     type(arg_type), dimension(3) :: meta_args = (/                        &
       arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, ANY_SPACE_1, ANY_SPACE_2), & ! operator_data
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_SPACE_1),              & ! field_aspc1
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_SPACE_2)               & ! field_aspc2
     /)
     integer :: operates_on = CELL_COLUMN
  end type setop_average_two_fields_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: setop_average_two_fields_kernel_code

  interface setop_average_two_fields_kernel_code
    module procedure &
      setop_average_two_fields_kernel_code_r_single, &
      setop_average_two_fields_kernel_code_r_double
  end interface

contains
  !> @brief  Assigns r_single operator as average of two fields.
  !> @param[in] cell             Index of cell to be randomised.
  !> @param[in] nlayers          No. of layers.
  !> @param[in] ncell_3d         Total no. of cells.
  !> @param[inout] operator_data Collection array of the 2D operator data at each cell.
  !> @param[in] ndf_aspc1        No. degrees of freedom for 2D operator data, 1st index.
  !> @param[in] ndf_aspc2        No. degrees of freedom for 2D operator data, 2nd index.
  subroutine setop_average_two_fields_kernel_code_r_single( cell,          &
                                                            nlayers,       &
                                                            ncell_3d,      &
                                                            operator_data, &
                                                            field_aspc1,   &
                                                            field_aspc2,   &
                                                            ndf_aspc1,     &
                                                            undf_aspc1,    &
                                                            map_aspc1,     &
                                                            ndf_aspc2,     &
                                                            undf_aspc2,    &
                                                            map_aspc2 )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: cell
    integer(kind=i_def), intent(in) :: ncell_3d
    integer(kind=i_def), intent(in) :: ndf_aspc1
    integer(kind=i_def), intent(in) :: undf_aspc1
    integer(kind=i_def), intent(in) :: ndf_aspc2
    integer(kind=i_def), intent(in) :: undf_aspc2
    real(kind=r_single), intent(inout), dimension(ncell_3d, ndf_aspc1, ndf_aspc2) :: operator_data
    real(kind=r_single), intent(in), dimension(undf_aspc1) :: field_aspc1
    real(kind=r_single), intent(in), dimension(undf_aspc2) :: field_aspc2
    integer(kind=i_def), intent(in), dimension(ndf_aspc1) :: map_aspc1
    integer(kind=i_def), intent(in), dimension(ndf_aspc2) :: map_aspc2

    ! Internal variables
    integer(kind=i_def) :: k, ik, df1, df2

    do k = 0, nlayers-1
      ik = (cell-1)*nlayers + k + 1
      do df1 = 1, ndf_aspc1
        do df2 = 1, ndf_aspc2
          operator_data(ik,df1,df2) = ( field_aspc1(map_aspc1(df1) + k) + field_aspc2(map_aspc2(df2) + k) ) / 2.0_r_single
        end do
      end do
    end do

  end subroutine setop_average_two_fields_kernel_code_r_single

  !> @brief  Assigns r_doube operator as average of two fields.
  !> @param[in] cell             Index of cell to be randomised.
  !> @param[in] nlayers          No. of layers.
  !> @param[in] ncell_3d         Total no. of cells.
  !> @param[inout] operator_data Collection array of the 2D operator data at each cell.
  !> @param[in] ndf_aspc1        No. degrees of freedom for 2D operator data, 1st index.
  !> @param[in] ndf_aspc2        No. degrees of freedom for 2D operator data, 2nd index.
  subroutine setop_average_two_fields_kernel_code_r_double( cell,          &
                                                            nlayers,       &
                                                            ncell_3d,      &
                                                            operator_data, &
                                                            field_aspc1,   &
                                                            field_aspc2,   &
                                                            ndf_aspc1,     &
                                                            undf_aspc1,    &
                                                            map_aspc1,     &
                                                            ndf_aspc2,     &
                                                            undf_aspc2,    &
                                                            map_aspc2 )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: cell
    integer(kind=i_def), intent(in) :: ncell_3d
    integer(kind=i_def), intent(in) :: ndf_aspc1
    integer(kind=i_def), intent(in) :: undf_aspc1
    integer(kind=i_def), intent(in) :: ndf_aspc2
    integer(kind=i_def), intent(in) :: undf_aspc2
    real(kind=r_double), intent(inout), dimension(ncell_3d, ndf_aspc1, ndf_aspc2) :: operator_data
    real(kind=r_double), intent(in), dimension(undf_aspc1) :: field_aspc1
    real(kind=r_double), intent(in), dimension(undf_aspc2) :: field_aspc2
    integer(kind=i_def), intent(in), dimension(ndf_aspc1) :: map_aspc1
    integer(kind=i_def), intent(in), dimension(ndf_aspc2) :: map_aspc2

    ! Internal variables
    integer(kind=i_def) :: k, ik, df1, df2

    do k = 0, nlayers-1
      ik = (cell-1)*nlayers + k + 1
      do df1 = 1, ndf_aspc1
        do df2 = 1, ndf_aspc2
          operator_data(ik,df1,df2) = ( field_aspc1(map_aspc1(df1) + k) + field_aspc2(map_aspc2(df2) + k) ) / 2.0_r_double
        end do
      end do
    end do

  end subroutine setop_average_two_fields_kernel_code_r_double

end module setop_average_two_fields_kernel_mod
