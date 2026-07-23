!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of map of field from W2h broken to W2h.
module adj_assemble_w2h_from_w2hb_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_READWRITE,     &
                                    CELL_COLUMN
  use constants_mod,         only : r_solver, i_def
  use fs_continuity_mod,     only : W2h, W2broken
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: adj_assemble_w2h_from_w2hb_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                      &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      W2h),     & ! field_w2h
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W2broken) & ! field_w2h_broken
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_assemble_w2h_from_w2hb_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_assemble_w2h_from_w2hb_code

contains

!> @brief Converts a broken W2h field into a continuous W2h field
!>
!> @param[in] nlayers Number of layers in the mesh
!> @param[in] field_w2h Field in the W2h space to be returned.
!> @param[in,out] field_w2h_broken Original field in W2h broken to be used.
!> @param[in] ndf_w2h Number of degrees of freedom per cell for W2h
!> @param[in] undf_w2h Number of (local) unique degrees of freedom for W2h
!> @param[in] map_w2h Dofmap for the cell at the base of the column for W2h
!> @param[in] ndf_w2h_broken Number of degrees of freedom per cell for W2h broken
!> @param[in] undf_w2h_broken Number of (local) unique degrees of freedom for W2h broken
!> @param[in] map_w2h_broken Dofmap for the cell at the base of the column for W2h broken
subroutine adj_assemble_w2h_from_w2hb_code( nlayers,          &
                                            field_w2h,        &
                                            field_w2h_broken, &
                                            ndf_w2h,          &
                                            undf_w2h,         &
                                            map_w2h,          &
                                            ndf_w2h_broken,   &
                                            undf_w2h_broken,  &
                                            map_w2h_broken )

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in) :: nlayers
  integer(kind=i_def),                            intent(in) :: ndf_w2h_broken, ndf_w2h
  integer(kind=i_def),                            intent(in) :: undf_w2h_broken, undf_w2h
  integer(kind=i_def), dimension(ndf_w2h_broken), intent(in) :: map_w2h_broken
  integer(kind=i_def), dimension(ndf_w2h),        intent(in) :: map_w2h

  real(kind=r_solver), dimension(undf_w2h),        intent(in)    :: field_w2h
  real(kind=r_solver), dimension(undf_w2h_broken), intent(inout) :: field_w2h_broken

  ! Internal variables
  integer(kind=i_def) :: df, k

  ! Loop over horizontal W2h DoFs
  do df = 1, ndf_w2h
    ! Loop over layers of mesh
    do k = 0, nlayers - 1
      field_w2h_broken(map_w2h_broken(df)+k) = field_w2h_broken(map_w2h_broken(df)+k) &
                                             + field_w2h(map_w2h(df)+k)
    end do
  end do

end subroutine adj_assemble_w2h_from_w2hb_code

end module adj_assemble_w2h_from_w2hb_kernel_mod
