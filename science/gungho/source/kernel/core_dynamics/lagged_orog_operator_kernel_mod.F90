!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Zero the horizontal-vertical correlations in a given W2 mass matrix
!!        operator.
!> @details This operator is used by the LHS when the sloped orography terms
!!        are lagged.
module lagged_orog_operator_kernel_mod

use constants_mod,           only: r_def, i_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type,             &
                                   GH_OPERATOR, GH_REAL, &
                                   GH_READ, GH_WRITE,    &
                                   CELL_COLUMN
use fs_continuity_mod,       only: W2

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: lagged_orog_operator_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                    &
       arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, W2), &
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W2, W2)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: lagged_orog_operator_kernel_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: lagged_orog_operator_kernel_code
contains

!> @param[in] cell Cell number
!> @param[in] nlayers Number of layers.
!> @param[in] ncell_3d_1 Ncell*ndf
!> @param[in,out] lagged_orog_operator to create
!> @param[in] ncell_3d_2 Ncell*ndf
!> @param[in] mass_matrix1 Non-lagged W2 mass matrix operator
!> @param[in] ndf1 Number of dofs per cell for space 1
subroutine lagged_orog_operator_kernel_code(cell,                 &
                                            nlayers,              &
                                            ncell_3d_1,           &
                                            lagged_orog_operator, &
                                            ncell_3d_2,           &
                                            mass_matrix1,         &
                                            ndf1)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ncell_3d_1
  integer(kind=i_def), intent(in) :: ncell_3d_2
  integer(kind=i_def), intent(in) :: ndf1

  real(kind=r_def), dimension(ncell_3d_1,ndf1,ndf1), intent(inout) :: lagged_orog_operator
  real(kind=r_def), dimension(ncell_3d_2,ndf1,ndf1), intent(in)    :: mass_matrix1

  ! Internal variables
  integer(kind=i_def) :: k, ik, df1, df2

  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers

    lagged_orog_operator(ik,:,:) = 0.0_r_def

    do df1 = 1,4
      do df2 = 1,4
        lagged_orog_operator(ik,df1,df2) = mass_matrix1(ik,df1,df2)
      end do
    end do
    do df1 = 5,6
      do df2 = 5,6
        lagged_orog_operator(ik,df1,df2) = mass_matrix1(ik,df1,df2)
      end do
    end do

  end do

end subroutine lagged_orog_operator_kernel_code

end module lagged_orog_operator_kernel_mod
