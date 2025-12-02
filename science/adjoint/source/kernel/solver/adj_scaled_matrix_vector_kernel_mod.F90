!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of lhs = y*z*matrix*x

module adj_scaled_matrix_vector_kernel_mod

  use argument_mod,      only: arg_type, &
                               GH_FIELD, GH_OPERATOR, GH_REAL, &
                               GH_READ, GH_READWRITE, GH_INC,  &
                               CELL_COLUMN
  use constants_mod,     only: i_def, r_solver
  use fs_continuity_mod, only: W2, W3
  use kernel_mod,        only: kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: adj_scaled_matrix_vector_kernel_type
    private
    type(ARG_TYPE) :: META_ARGS(5) = (/                     &
      arg_type(GH_FIELD,    GH_REAL, GH_READ,      W2),     &
      arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3),     &
      arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W2, W3), &
      arg_type(GH_FIELD,    GH_REAL, GH_READ,      W2),     &
      arg_type(GH_FIELD,    GH_REAL, GH_READ,      W2)      &
      /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_scaled_matrix_vector_code
  end type adj_scaled_matrix_vector_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  public :: adj_scaled_matrix_vector_code

contains

!> @brief Compute x = (matrix^T)*(lhs*y*z)
!> @details matrix maps x space to lhs space, y and z are fields in the same
!!          space as lhs
!> @param[in]     cell     Horizontal cell index
!> @param[in]     nlayers  Number of layers
!> @param[in]     lhs      Input data
!> @param[in,out] x        Output data
!> @param[in]     ncell_3d Total number of cells
!> @param[in]     matrix   Local matrix assembly form of the operator A
!> @param[in]     y        Field to scale output by
!> @param[in]     z        Second field to scale output by
!> @param[in]     ndf1     No. degrees of freedom per cell for y/z/lhs
!> @param[in]     undf1    Unique no. degrees of freedom  for y/z/lhs
!> @param[in]     map1     Dofmap for cell at base of column for y/z/lhs
!> @param[in]     ndf2     Number of degrees of freedom per cell for x
!> @param[in]     undf2    Unique number of degrees of freedom for x
!> @param[in]     map2     Dofmap for cell at base of column for x
subroutine adj_scaled_matrix_vector_code(cell, nlayers, lhs, &
                                         x,                  &
                                         ncell_3d,           &
                                         matrix, y, z,       &
                                         ndf1, undf1, map1,  &
                                         ndf2, undf2, map2)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell, nlayers
  integer(kind=i_def), intent(in) :: undf1, undf2
  real(kind=r_solver), intent(inout) :: lhs(undf1), x(undf2)
  integer(kind=i_def), intent(in) :: ncell_3d, ndf1, ndf2
  real(kind=r_solver), intent(in) :: matrix(ncell_3d,ndf1,ndf2)
  real(kind=r_solver), intent(in) :: y(undf1), z(undf1)
  integer(kind=i_def), dimension(ndf1), intent(in) :: map1(ndf1), map2(ndf2)

  ! Local variables
  integer(kind=i_def) :: df
  integer(kind=i_def) :: df2
  integer(kind=i_def) :: k
  integer(kind=i_def) :: ik

  do df = ndf1, 1, -1
    do df2 = ndf2, 1, -1
      do k = nlayers - 1, 0, -1
        ik = cell * nlayers + k - nlayers + 1
        x(k + map2(df2)) = x(k + map2(df2)) + &
          matrix( ik, df, df2 ) * lhs(map1(df) + k) &
          * y(k + map1(df)) * z(k + map1(df))
      end do
    end do
  end do

end subroutine adj_scaled_matrix_vector_code

end module adj_scaled_matrix_vector_kernel_mod
