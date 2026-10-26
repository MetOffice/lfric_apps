!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of dg_matrix_vector_kernel. This version is for use on
!>        discontinous spaces and over-writes the output field entirely,
!>        hence can only be used for W3 spaces
module adj_dg_matrix_vector_kernel_mod

  use constants_mod, only : i_def, r_single, r_double
  use kernel_mod,    only : kernel_type
  use argument_mod,  only : arg_type,                  &
                            GH_FIELD, GH_OPERATOR,     &
                            GH_READ,                   &
                            GH_INC, GH_READWRITE,      &
                            GH_REAL, ANY_SPACE_1,      &
                            ANY_DISCONTINUOUS_SPACE_1, &
                            CELL_COLUMN

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: adj_dg_matrix_vector_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                           &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lhs
         arg_type(GH_FIELD,    GH_REAL, GH_INC,       ANY_SPACE_1),               & ! x
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1,  &
                                                      ANY_SPACE_1) & ! matrix
         /)
    integer :: operates_on = CELL_COLUMN
  end type adj_dg_matrix_vector_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_dg_matrix_vector_code

  ! Generic interface for real32 and real64 types
  interface adj_dg_matrix_vector_code
    module procedure  &
      adj_dg_matrix_vector_code_r_single, &
      adj_dg_matrix_vector_code_r_double
  end interface

  contains

  !> @brief Computes adjoint of lhs = matrix*x for discontinuous function spaces
  !> @param[in]      cell      Horizontal cell index
  !> @param[in]      nlayers   Number of layers
  !> @param[in,out]  lhs       Left hand side field of lhs = matrix*x
  !> @param[in,out]  x         Right hand side field of lhs = matrix*x
  !> @param[in]      ncell_3d  Total number of cells
  !> @param[in]      matrix    Matrix values in LMA form
  !> @param[in]      ndf1      Number of degrees of freedom per cell for the output field
  !> @param[in]      undf1     Unique number of degrees of freedom  for the output field
  !> @param[in]      map1      Dofmap for the cell at the base of the column for the output field
  !> @param[in]      ndf2      Number of degrees of freedom per cell for the input field
  !> @param[in]      undf2     Unique number of degrees of freedom for the input field
  !> @param[in]      map2      Dofmap for the cell at the base of the column for the input field

  ! R_SINGLE PRECISION
  ! ==================
  subroutine adj_dg_matrix_vector_code_r_single(cell,              &
                                                nlayers,           &
                                                lhs, x,            &
                                                ncell_3d,          &
                                                matrix,            &
                                                ndf1, undf1, map1, &
                                                ndf2, undf2, map2)

    implicit none

    ! Arguments
    integer(kind=i_def),                   intent(in) :: cell, nlayers, ncell_3d
    integer(kind=i_def),                   intent(in) :: undf1, ndf1
    integer(kind=i_def),                   intent(in) :: undf2, ndf2
    integer(kind=i_def), dimension(ndf1),  intent(in) :: map1
    integer(kind=i_def), dimension(ndf2),  intent(in) :: map2

    real(kind=r_single), dimension(undf2),              intent(inout) :: x
    real(kind=r_single), dimension(undf1),              intent(inout) :: lhs
    real(kind=r_single), dimension(ncell_3d,ndf1,ndf2), intent(in)    :: matrix

    ! Internal variables
    integer(kind=i_def) :: df1, df2, ik, i1, i2, nl

    nl = nlayers - 1

    ik = cell * nlayers - nlayers + 1
    do df2 = ndf2, 1, -1
      i2 = map2(df2)
      do df1 = ndf1, 1, -1
        i1 = map1(df1)
        x(i2:i2+nl) = x(i2:i2+nl) + matrix(ik:ik+nl,df1,df2)*lhs(i1:i1+nl)
      end do
    end do

    do df1 = ndf1, 1, -1
      i1 = map1(df1)
      lhs(i1:i1+nl) = 0.0_r_single
    end do

  end subroutine adj_dg_matrix_vector_code_r_single

  ! R_DOUBLE PRECISION
  ! ==================
  subroutine adj_dg_matrix_vector_code_r_double(cell,              &
                                                nlayers,           &
                                                lhs, x,            &
                                                ncell_3d,          &
                                                matrix,            &
                                                ndf1, undf1, map1, &
                                                ndf2, undf2, map2)

    implicit none

    ! Arguments
    integer(kind=i_def),                   intent(in) :: cell, nlayers, ncell_3d
    integer(kind=i_def),                   intent(in) :: undf1, ndf1
    integer(kind=i_def),                   intent(in) :: undf2, ndf2
    integer(kind=i_def), dimension(ndf1),  intent(in) :: map1
    integer(kind=i_def), dimension(ndf2),  intent(in) :: map2

    real(kind=r_double), dimension(undf2),              intent(inout) :: x
    real(kind=r_double), dimension(undf1),              intent(inout) :: lhs
    real(kind=r_double), dimension(ncell_3d,ndf1,ndf2), intent(in)    :: matrix

    ! Internal variables
    integer(kind=i_def) :: df1, df2, ik, i1, i2, nl

    nl = nlayers - 1

    ik = cell * nlayers - nlayers + 1
    do df2 = ndf2, 1, -1
      i2 = map2(df2)
      do df1 = ndf1, 1, -1
        i1 = map1(df1)
        x(i2:i2+nl) = x(i2:i2+nl) + matrix(ik:ik+nl,df1,df2)*lhs(i1:i1+nl)
      end do
    end do

    do df1 = ndf1, 1, -1
      i1 = map1(df1)
      lhs(i1:i1+nl) = 0.0_r_double
    end do

  end subroutine adj_dg_matrix_vector_code_r_double

end module adj_dg_matrix_vector_kernel_mod
