!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Adjoint of application of the strong curl
!>
module adj_strong_curl_kernel_mod

  use argument_mod,      only : arg_type, func_type,     &
                                GH_FIELD, GH_REAL,       &
                                GH_READINC, GH_INC,      &
                                GH_DIFF_BASIS, GH_BASIS, &
                                CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W1, W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: adj_strong_curl_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/               &
         arg_type(GH_FIELD, GH_REAL, GH_READINC, W2), &
         arg_type(GH_FIELD, GH_REAL, GH_INC,     W1)  &
         /)
    type(func_type) :: meta_funcs(2) = (/             &
         func_type(W2, GH_BASIS),                     &
         func_type(W1, GH_DIFF_BASIS)                 &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: adj_strong_curl_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_strong_curl_code

contains

!> @brief Adjoint of kernel to compute the strong curl
!! @param[in] nlayers Number of layers
!! @param[in,out] xi Field to contain curl of u
!! @param[in,out] u Wind field
!! @param[in] ndf2 Number of degrees of freedom per cell for W2
!! @param[in] undf2 Number of unique degrees of freedom for W2
!! @param[in] map2 Dofmap for the cell at the base of the column for W2
!! @param[in] basis_w2 W2 basis functions evaluated at nodal points of W2
!! @param[in] ndf1 Number of degrees of freedom per cell for W1
!! @param[in] undf1 Number of unique degrees of freedom for W1
!! @param[in] map1 Dofmap for the cell at the base of the column for the field to be advected
!! @param[in] diff_basis_w1 Differential W1 basis functions evaluated at nodal points of W2
subroutine adj_strong_curl_code(nlayers,                          &
                                xi, u,                            &
                                ndf2, undf2, map2, basis_w2,      &
                                basis_w2_on_w1,                   &
                                ndf1, undf1, map1, diff_basis_w1, &
                                diff_basis_w1_on_w1               &
                               )

  implicit none

  ! Arguments
  integer(kind=i_def),                      intent(in)    :: nlayers
  integer(kind=i_def),                      intent(in)    :: ndf1, undf1, &
                                                             ndf2, undf2
  integer(kind=i_def), dimension(ndf1),     intent(in)    :: map1
  integer(kind=i_def), dimension(ndf2),     intent(in)    :: map2
  real(kind=r_def), dimension(3,ndf2,ndf2), intent(in)    :: basis_w2
  real(kind=r_def), dimension(3,ndf2,ndf1), intent(in)    :: basis_w2_on_w1
  real(kind=r_def), dimension(3,ndf1,ndf2), intent(in)    :: diff_basis_w1
  real(kind=r_def), dimension(3,ndf1,ndf1), intent(in)    :: diff_basis_w1_on_w1
  real(kind=r_def), dimension(undf2),       intent(inout) :: xi
  real(kind=r_def), dimension(undf1),       intent(inout) :: u

  ! Internal variables
  integer(kind=i_def)            :: df1, df2, k
  real(kind=r_def), dimension(3) :: curl_u

  do k = nlayers-1, 0, -1
    do df2 = ndf2, 1, -1
      curl_u = basis_w2(:,df2,df2)*xi( map2(df2)+k )
      xi( map2(df2)+k ) = 0.0_r_def
      do df1 = ndf1, 1, -1
        u(map1(df1)+k) = u(map1(df1)+k) + dot_product(diff_basis_w1(:,df1,df2), curl_u)
      end do
    end do
  end do

end subroutine adj_strong_curl_code

end module adj_strong_curl_kernel_mod
