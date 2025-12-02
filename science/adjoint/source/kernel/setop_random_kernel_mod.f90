!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Module with functions to initialise operators with random values.
!> All psyad generated adjoint tests for LFRic kernels use this
!> module, even if the kernel is never called.

module setop_random_kernel_mod

  use argument_mod,                          only : any_space_1, any_space_2, &
                                                    arg_type, cell_column,    &
                                                    gh_operator, gh_real,     &
                                                    gh_write
  use kernel_mod,                            only : kernel_type
  use constants_mod,                         only : i_def, r_single, r_double

  implicit none

  type, extends(kernel_type) :: setop_random_kernel_type
     type(arg_type), dimension(1) :: meta_args = (/                       &
       arg_type(gh_operator, gh_real, gh_write, any_space_1, any_space_2) &
     /)
     integer :: operates_on = CELL_COLUMN
  end type setop_random_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: setop_random_kernel_code

  interface setop_random_kernel_code
    module procedure &
      setop_random_kernel_code_r_single, &
      setop_random_kernel_code_r_double
  end interface

contains
  !> @brief     Generates randomly assigned r_single operator for use in psyad adjoint tests.
  !> @details   Calls random_number on the r_single 2D operator data at each cell.
  !> @param[in] cell             Index of cell to be randomised.
  !> @param[in] nlayers          No. of layers.
  !> @param[in] ncell_3d         Total no. of cells.
  !> @param[inout] operator_data Collection array of the 2D operator data at each cell.
  !> @param[in] ndf_aspc1        No. degrees of freedom for 2D operator data, 1st index.
  !> @param[in] ndf_aspc2        No. degrees of freedom for 2D operator data, 2nd index.
  subroutine setop_random_kernel_code_r_single(cell, nlayers, ncell_3d, &
                                            operator_data, ndf_aspc1, ndf_aspc2)

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: cell
    integer(kind=i_def), intent(in) :: ncell_3d
    integer(kind=i_def), intent(in) :: ndf_aspc1, ndf_aspc2
    real(kind=r_single), intent(inout), dimension(ncell_3d, ndf_aspc1, ndf_aspc2) :: operator_data
    integer(kind=i_def) :: k, ik

    do k = 0, nlayers-1
      ik = (cell-1)*nlayers + k + 1
      call random_number(operator_data(ik,:,:))
    end do

  end subroutine setop_random_kernel_code_r_single

  !> @brief     Generates randomly assigned r_double operator for use in psyad adjoint tests.
  !> @details   Calls random_number on the r_double 2D operator data at each cell.
  !> @param[in] cell             Index of cell to be randomised.
  !> @param[in] nlayers          No. of layers.
  !> @param[in] ncell_3d         Total no. of cells.
  !> @param[inout] operator_data Collection array of the 2D operator data at each cell.
  !> @param[in] ndf_aspc1        No. degrees of freedom for 2D operator data, 1st index.
  !> @param[in] ndf_aspc2        No. degrees of freedom for 2D operator data, 2nd index.
  subroutine setop_random_kernel_code_r_double(cell, nlayers, ncell_3d, &
                                            operator_data, ndf_aspc1, ndf_aspc2)

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: cell
    integer(kind=i_def), intent(in) :: ncell_3d
    integer(kind=i_def), intent(in) :: ndf_aspc1, ndf_aspc2
    real(kind=r_double), intent(inout), dimension(ncell_3d, ndf_aspc1, ndf_aspc2) :: operator_data
    integer(kind=i_def) :: k, ik

    do k = 0, nlayers-1
      ik = (cell-1)*nlayers + k + 1
      call random_number(operator_data(ik,:,:))
    end do

  end subroutine setop_random_kernel_code_r_double

end module setop_random_kernel_mod
