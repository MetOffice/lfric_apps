!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of map from W2 broken to W2.
module adj_sci_average_w2b_to_w2_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_INC,           &
                                    GH_READWRITE,              &
                                    CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2, W2broken
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : N, E, S, W, T, B

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: adj_average_w2b_to_w2_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_INC,       W2),                   &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W2broken),             &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      W2)                    &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_average_w2b_to_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_average_w2b_to_w2_code

contains

!> @brief Adjoint of conversion of a broken W2 field into a continuous W2 field
!>
!> @param[in] nlayers Number of layers in the mesh
!> @param[in,out] field_w2 Field in the W2 space to be returned.
!> @param[in,out] field_w2_broken Original field in W2 broken to be used.
!> @param[in] rmultiplicity_w2 Reciprocal of nodal multiplicity field for W2
!> @param[in] ndf_w2 Number of degrees of freedom per cell for W2
!> @param[in] undf_w2 Number of (local) unique degrees of freedom for W2
!> @param[in] map_w2 Dofmap for the cell at the base of the column for W2
!> @param[in] ndf_w2_broken Number of degrees of freedom per cell for W2 broken
!> @param[in] undf_w2_broken Number of (local) unique degrees of freedom for W2 broken
!> @param[in] map_w2_broken Dofmap for the cell at the base of the column for W2 broken
subroutine adj_average_w2b_to_w2_code( nlayers,          &
                                       field_w2,         &
                                       field_w2_broken,  &
                                       rmultiplicity_w2, &
                                       ndf_w2,           &
                                       undf_w2,          &
                                       map_w2,           &
                                       ndf_w2_broken,    &
                                       undf_w2_broken,   &
                                       map_w2_broken     &
                                     )

  implicit none

  ! Arguments
  integer(kind=i_def),                               intent(in) :: nlayers
  integer(kind=i_def),                               intent(in) :: ndf_w2_broken, ndf_w2
  integer(kind=i_def),                               intent(in) :: undf_w2_broken, undf_w2
  integer(kind=i_def), dimension(ndf_w2_broken),     intent(in) :: map_w2_broken
  integer(kind=i_def), dimension(ndf_w2),            intent(in) :: map_w2

  real(kind=r_def),    dimension(undf_w2),        intent(inout) :: field_w2
  real(kind=r_def),    dimension(undf_w2_broken), intent(inout) :: field_w2_broken
  real(kind=r_def),    dimension(undf_w2),           intent(in) :: rmultiplicity_w2

  ! Internal variables
  integer(kind=i_def) :: df, k, j
  integer(kind=i_def) :: horizontal_dofs(4)

  horizontal_dofs = (/ N, E, S, W /)

  field_w2_broken(map_w2_broken(T)+nlayers-1) = field_w2_broken(map_w2_broken(T)+nlayers-1) &
                                              + field_w2(map_w2(T)+nlayers-1)
  field_w2(map_w2(T)+nlayers-1) = 0.0_r_def

  do k = nlayers - 1, 1, -1
    field_w2_broken(map_w2_broken(B)+k) = field_w2_broken(map_w2_broken(B)+k) &
                                        + 0.5_r_def*field_w2(map_w2(B)+k)
    field_w2_broken(map_w2_broken(T)+k-1) = field_w2_broken(map_w2_broken(T)+k-1) &
                                          + 0.5_r_def*field_w2(map_w2(B)+k)
    field_w2(map_w2(B)+k) = 0.0_r_def
  end do

  field_w2_broken(map_w2_broken(B)) = field_w2_broken(map_w2_broken(B)) &
                                    + field_w2(map_w2(B))
  field_w2(map_w2(B)) = 0.0_r_def

  do j = size(horizontal_dofs), 1, -1
    df = horizontal_dofs(j)
    do k = nlayers - 1, 0, -1
      field_w2_broken(map_w2_broken(df)+k) = field_w2_broken(map_w2_broken(df)+k) &
        + rmultiplicity_w2(map_w2(df)+k) * field_w2(map_w2(df)+k)
    end do
  end do

end subroutine adj_average_w2b_to_w2_code

end module adj_sci_average_w2b_to_w2_kernel_mod
