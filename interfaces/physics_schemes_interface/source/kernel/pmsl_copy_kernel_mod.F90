!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Holds code to copy PMSL in halo

module pmsl_copy_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_SCALAR,       &
                                  GH_READ, GH_WRITE,         &
                                  GH_REAL, OWNED_AND_HALO_CELL_COLUMN,      &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  STENCIL, CROSS2D
  use fs_continuity_mod,    only: WTHETA, W2
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: pmsl_copy_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! exner_in
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &
         /) ! exner_out
    integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
  contains
    procedure, nopass :: pmsl_copy_code
  end type pmsl_copy_kernel_type

  public :: pmsl_copy_code

contains

  !> @brief Copy field including a depth in the halo
  !> @param[in]     nlayers       The number of layers
  !> @param[in]     exner_in      Initial guess for exner at MSL
  !> @param[in,out] exner_out     Next guess for exner at MSL
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell for 2d fields
  !> @param[in]     undf_2d       Number of total degrees of freedom for 2d fields
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column for 2d fields
  subroutine pmsl_copy_code(nlayers,                         &
                            exner_in,                        &
                            exner_out,                       &
                            ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    ! Arguments passed explicitly from algorithm
    real(kind=r_def),    intent(in), dimension(undf_2d) :: exner_in
    real(kind=r_def),    intent(inout), dimension(undf_2d) :: exner_out

    exner_out(map_2d(1)) = exner_in(map_2d(1))

  end subroutine pmsl_copy_code

end module pmsl_copy_kernel_mod
