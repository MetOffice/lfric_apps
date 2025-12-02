!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module scaled_matrix_vector_kernel_mod

  use argument_mod,      only : arg_type,              &
                                GH_FIELD, GH_OPERATOR, &
                                GH_READ, GH_INC,       &
                                GH_REAL, CELL_COLUMN
  use constants_mod,     only : r_solver, i_def
  use fs_continuity_mod, only : W2, W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: scaled_matrix_vector_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                   &
         arg_type(GH_FIELD,    GH_REAL, GH_INC,  W2),     &
         arg_type(GH_FIELD,    GH_REAL, GH_READ, W3),     &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ, W2, W3), &
         arg_type(GH_FIELD,    GH_REAL, GH_READ, W2),     &
         arg_type(GH_FIELD,    GH_REAL, GH_READ, W2)      &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: scaled_matrix_vector_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: scaled_matrix_vector_code

contains

!> @brief Computes lhs = y*z*matrix*x where matrix maps from x space to lhs space
!>        and y and z are fields in the same space as lhs
!> @param[in] cell Horizontal cell index
!! @param[in] nlayers Number of layers
!> @param[in,out] lhs Output lhs (A*x)
!! @param[in] x Input data
!! @param[in] ncell_3d Total number of cells
!! @param[in] matrix Local matrix assembly form of the operator A
!! @param[in] y Field to scale output by
!! @param[in] z Second field to scale output by
!! @param[in] ndf1 Number of degrees of freedom per cell for the output field
!! @param[in] undf1 Unique number of degrees of freedom  for the output field
!! @param[in] map1 Dofmap for the cell at the base of the column for the output field
!! @param[in] ndf2 Number of degrees of freedom per cell for the input field
!! @param[in] undf2 Unique number of degrees of freedom for the input field
!! @param[in] map2 Dofmap for the cell at the base of the column for the input field
subroutine scaled_matrix_vector_code(cell,              &
                                     nlayers,           &
                                     lhs, x,            &
                                     ncell_3d,          &
                                     matrix,            &
                                     y,                 &
                                     z,                 &
                                     ndf1, undf1, map1, &
                                     ndf2, undf2, map2)

  implicit none

  ! Arguments
  integer(kind=i_def),                   intent(in) :: cell, nlayers, ncell_3d
  integer(kind=i_def),                   intent(in) :: undf1, ndf1
  integer(kind=i_def),                   intent(in) :: undf2, ndf2
  integer(kind=i_def), dimension(ndf1),  intent(in) :: map1
  integer(kind=i_def), dimension(ndf2),  intent(in) :: map2

  real(kind=r_solver), dimension(undf2),              intent(in)    :: x
  real(kind=r_solver), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_solver), dimension(ncell_3d,ndf1,ndf2), intent(in)    :: matrix
  real(kind=r_solver), dimension(undf1),              intent(in)    :: y
  real(kind=r_solver), dimension(undf1),              intent(in)    :: z

  ! Internal variables
  integer(kind=i_def) :: df, df2, ij, nl, i1, i2

  ij = (cell-1)*nlayers + 1
  nl = nlayers - 1
  do df = 1, ndf1
    i1 = map1(df)
    do df2 = 1, ndf2
      i2 = map2(df2)
      lhs(i1:i1+nl) = lhs(i1:i1+nl) &
                    + matrix(ij:ij+nl, df, df2)*x(i2:i2+nl)*y(i1:i1+nl)*z(i1:i1+nl)
    end do
  end do

end subroutine scaled_matrix_vector_code

end module scaled_matrix_vector_kernel_mod
