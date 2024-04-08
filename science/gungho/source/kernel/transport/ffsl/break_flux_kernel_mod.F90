!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Maps a field from W2 broken to W2.
!> @details "Unbreaks" a W2 field by averaging values on either side of
!>          broken facets.
!>          This kernel only works for the lowest-order elements.
module break_flux_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_WRITE,         &
                                    CELL_COLUMN
  use constants_mod,         only : r_single, r_double, i_def
  use fs_continuity_mod,     only : W2H, W2broken
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: break_flux_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                   &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2broken), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2H),      &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2H)       &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: break_flux_code

  interface break_flux_code
    module procedure &
      break_flux_code_r_single, &
      break_flux_code_r_double
    end interface


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
subroutine break_flux_code_r_single( nlayers,          &
                                     broken_flux,      &
                                     flux_x,           &
                                     flux_y,           &
                                     ndf_w2h_broken,   &
                                     undf_w2h_broken,  &
                                     map_w2h_broken,   &
                                     ndf_w2h,          &
                                     undf_w2h,         &
                                     map_w2h           &
                                   )

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in) :: nlayers
  integer(kind=i_def),                            intent(in) :: ndf_w2h_broken, ndf_w2h
  integer(kind=i_def),                            intent(in) :: undf_w2h_broken, undf_w2h
  integer(kind=i_def), dimension(ndf_w2h_broken), intent(in) :: map_w2h_broken
  integer(kind=i_def), dimension(ndf_w2h),        intent(in) :: map_w2h

  real(kind=r_single),    dimension(undf_w2h_broken), intent(inout) :: broken_flux
  real(kind=r_single),    dimension(undf_w2h),        intent(in)    :: flux_x, flux_y

  ! Internal variables
  integer(kind=i_def) :: k

  ! Loop over layers of mesh
  do k = 0, nlayers - 1
    broken_flux(map_w2h_broken(1)+k) = flux_x(map_w2h(1)+k)
    broken_flux(map_w2h_broken(2)+k) = flux_y(map_w2h(2)+k)
    broken_flux(map_w2h_broken(3)+k) = flux_x(map_w2h(3)+k)
    broken_flux(map_w2h_broken(4)+k) = flux_y(map_w2h(4)+k)
  end do

end subroutine break_flux_code_r_single

subroutine break_flux_code_r_double( nlayers,          &
                                     broken_flux,      &
                                     flux_x,           &
                                     flux_y,           &
                                     ndf_w2h_broken,   &
                                     undf_w2h_broken,  &
                                     map_w2h_broken,   &
                                     ndf_w2h,          &
                                     undf_w2h,         &
                                     map_w2h           &
                                   )

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in) :: nlayers
  integer(kind=i_def),                            intent(in) :: ndf_w2h_broken, ndf_w2h
  integer(kind=i_def),                            intent(in) :: undf_w2h_broken, undf_w2h
  integer(kind=i_def), dimension(ndf_w2h_broken), intent(in) :: map_w2h_broken
  integer(kind=i_def), dimension(ndf_w2h),        intent(in) :: map_w2h

  real(kind=r_double),    dimension(undf_w2h_broken), intent(inout) :: broken_flux
  real(kind=r_double),    dimension(undf_w2h),        intent(in)    :: flux_x, flux_y

  ! Internal variables
  integer(kind=i_def) :: k

  ! Loop over layers of mesh
  do k = 0, nlayers - 1
    broken_flux(map_w2h_broken(1)+k) = flux_x(map_w2h(1)+k)
    broken_flux(map_w2h_broken(2)+k) = flux_y(map_w2h(2)+k)
    broken_flux(map_w2h_broken(3)+k) = flux_x(map_w2h(3)+k)
    broken_flux(map_w2h_broken(4)+k) = flux_y(map_w2h(4)+k)
  end do


end subroutine break_flux_code_r_double

end module break_flux_kernel_mod
