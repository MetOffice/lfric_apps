!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes adjoint of the vertical mass flux using PCM for the adjoint model.

module atl_ffsl_flux_z_constant_kernel_mod

use argument_mod,                   only : arg_type,              &
                                           GH_FIELD, GH_REAL,     &
                                           GH_READ,               &
                                           GH_READWRITE,          &
                                           GH_SCALAR, CELL_COLUMN
use fs_continuity_mod,              only : W3, W2v
use constants_mod,                  only : r_tran, i_def, l_def, EPS_R_TRAN
use kernel_mod,                     only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: atl_ffsl_flux_z_constant_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                      &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W2v), & ! flux_pert
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W2v), & ! flux_ls
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2v), & ! ls_frac_wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2v), & ! dep_dist
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W2v), & ! frac_wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),  & ! field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),  & ! ls_field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),  & ! detj
       arg_type(GH_SCALAR, GH_REAL,    GH_READ)            & ! dt
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: atl_ffsl_flux_z_constant_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: atl_ffsl_flux_z_constant_code

contains

!> @brief Computes the adjoint of the mass flux for FFSL using PCM in the z direction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] flux_pert      The pert field ls wind flux to be computed
!> @param[in,out] flux_ls        The ls field pert wind flux to be computed
!> @param[in]     ls_frac_wind   The fractional vertical ls wind
!> @param[in]     dep_dist       The vertical ls departure points
!> @param[in,out] frac_wind      The fractional vertical perturbation wind
!> @param[in,out] field          The field to construct the first flux
!> @param[in]     ls_field       The ls_field to construct the second flux
!> @param[in]     detj           Volume of cells
!> @param[in]     dt             Time step
!> @param[in]     ndf_w2v        Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v       Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v        The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3         Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3        Number of unique degrees of freedom for W3
!> @param[in]     map_w3         The dofmap for the cell at the base of the column
subroutine atl_ffsl_flux_z_constant_code( nlayers,        &
                                          flux_pert,      &
                                          flux_ls,        &
                                          ls_frac_wind,   &
                                          dep_dist,       &
                                          frac_wind,      &
                                          field,          &
                                          ls_field,       &
                                          detj,           &
                                          dt,             &
                                          ndf_w2v,        &
                                          undf_w2v,       &
                                          map_w2v,        &
                                          ndf_w3,         &
                                          undf_w3,        &
                                          map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: ndf_w3
  real(kind=r_tran),   intent(inout) :: flux_pert(undf_w2v)
  real(kind=r_tran),   intent(inout) :: flux_ls(undf_w2v)
  real(kind=r_tran),   intent(inout) :: field(undf_w3)
  real(kind=r_tran),   intent(in)    :: ls_field(undf_w3)
  real(kind=r_tran),   intent(in)    :: ls_frac_wind(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
  real(kind=r_tran),   intent(inout) :: frac_wind(undf_w2v)
  real(kind=r_tran),   intent(in)    :: detj(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  real(kind=r_tran),   intent(in)    :: dt

  ! Internal variables
  integer(kind=i_def) :: k, i, w2v_idx, w3_idx
  integer(kind=i_def) :: int_displacement, sign_displacement
  integer(kind=i_def) :: dep_cell_idx, sign_offset
  integer(kind=i_def) :: lowest_whole_cell, highest_whole_cell

  real(kind=r_tran)   :: displacement, frac_dist
  real(kind=r_tran)   :: mass_whole_cells

  w3_idx = map_w3(1)
  w2v_idx = map_w2v(1)

  flux_ls(w2v_idx + nlayers) = 0.0_r_tran
  flux_pert(w2v_idx + nlayers) = 0.0_r_tran

  ! Loop through faces
  do k = nlayers - 1, 1, -1

    ! Pull out departure point, and separate into integer / frac parts
    displacement = dep_dist(w2v_idx + k)
    int_displacement = int(displacement, i_def)
    frac_dist = displacement - real(int_displacement, r_tran)
    sign_displacement = int(sign(1.0_r_tran, displacement))

    ! Set an offset for the stencil index, based on dep point sign
    sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

    ! Determine departure cell
    dep_cell_idx = k - int_displacement + sign_offset - 1

    lowest_whole_cell = min(dep_cell_idx + 1, k)
    highest_whole_cell = max(dep_cell_idx, k) - 1

    frac_wind(w2v_idx + k) = frac_wind(w2v_idx + k)                       &
               + (flux_ls(w2v_idx + k) * ls_field(w3_idx + dep_cell_idx)  &
               ) / dt
    flux_ls(w2v_idx + k) = 0.0_r_tran

    field(w3_idx + dep_cell_idx) = field(w3_idx + dep_cell_idx)           &
                + (ls_frac_wind(w2v_idx + k) * flux_pert(w2v_idx + k)) / dt
    mass_whole_cells = (sign(1.0_r_tran, displacement) * flux_pert(w2v_idx + k)) / dt
    flux_pert(w2v_idx + k) = 0.0_r_tran

    do i = highest_whole_cell, lowest_whole_cell, -1
      field(w3_idx + i) = field(w3_idx + i) + mass_whole_cells * detj(w3_idx + i)
    end do

  end do

  flux_ls(w2v_idx) = 0.0_r_tran
  flux_pert(w2v_idx) = 0.0_r_tran

end subroutine atl_ffsl_flux_z_constant_code

end module atl_ffsl_flux_z_constant_kernel_mod
