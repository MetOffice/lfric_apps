!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Maps a field from W2 broken to W2.
!> @details "Unbreaks" a W2 field by averaging values on either side of
!>          broken facets.
!>          This kernel only works for the lowest-order elements.
module assemble_w2h_from_w2hb_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_INC,           &
                                    CELL_COLUMN
  use constants_mod,         only : r_solver, i_def
  use fs_continuity_mod,     only : W2h, W2broken
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
  type, public, extends(kernel_type) :: assemble_w2h_from_w2hb_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                  &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2h),      &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2broken) &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: assemble_w2h_from_w2hb_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: assemble_w2h_from_w2hb_code

contains

!> @brief Converts a broken W2 field into a continuous W2 field
!>
!> @param[in] nlayers Number of layers in the mesh
!> @param[in,out] field_w2 Field in the W2 space to be returned.
!> @param[in] field_w2_broken Original field in W2 broken to be used.
!> @param[in] rmultiplicity_w2 Reciprocal of nodal multiplicity field for W2
!> @param[in] ndf_w2 Number of degrees of freedom per cell for W2
!> @param[in] undf_w2 Number of (local) unique degrees of freedom for W2
!> @param[in] map_w2 Dofmap for the cell at the base of the column for W2
!> @param[in] ndf_w2_broken Number of degrees of freedom per cell for W2 broken
!> @param[in] undf_w2_broken Number of (local) unique degrees of freedom for W2 broken
!> @param[in] map_w2_broken Dofmap for the cell at the base of the column for W2 broken
subroutine assemble_w2h_from_w2hb_code( nlayers,          &
                                   field_w2,         &
                                   field_w2_broken,  &
                                   ndf_w2,           &
                                   undf_w2,          &
                                   map_w2,           &
                                   ndf_w2_broken,    &
                                   undf_w2_broken,   &
                                   map_w2_broken     &
                                 )

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in) :: nlayers
  integer(kind=i_def),                            intent(in) :: ndf_w2_broken, ndf_w2
  integer(kind=i_def),                            intent(in) :: undf_w2_broken, undf_w2
  integer(kind=i_def), dimension(ndf_w2_broken),  intent(in) :: map_w2_broken
  integer(kind=i_def), dimension(ndf_w2),         intent(in) :: map_w2

  real(kind=r_solver),    dimension(undf_w2),     intent(inout) :: field_w2
  real(kind=r_solver),    dimension(undf_w2_broken), intent(in) :: field_w2_broken

  ! Internal variables
  integer(kind=i_def) :: df, k

  ! Loop over horizontal W2 DoFs
  ! Contribution from broken field is based on multiplicity
  do df = 1, ndf_w2
    ! Loop over layers of mesh
    do k = 0, nlayers - 1
      field_w2(map_w2(df)+k) = field_w2(map_w2(df)+k)              &
        + field_w2_broken(map_w2_broken(df)+k)

    end do
  end do

end subroutine assemble_w2h_from_w2hb_code

end module assemble_w2h_from_w2hb_kernel_mod
