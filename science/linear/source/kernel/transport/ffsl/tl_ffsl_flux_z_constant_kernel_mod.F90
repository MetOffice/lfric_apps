!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes the vertical mass flux using PCM for the tangent linear model.
!> @details This kernel computes both parts of the vertical flux using the PCM
!!          scheme for tangent linear model transport. The first part of the flux is
!!          the reconstruction of the perturbation field with the ls wind. The second
!!          part uses a constant reconstruction of the ls field in the ls wind
!!          departure cell, but is multiplied by the fractional perturbation wind.
!!          The PCM reconstruction is equivalent to a first order upwind scheme.
!!          This kernel is designed to work in the vertical direction only and takes
!!          into account the vertical boundaries and grid spacing.
!!
!!          Note that this kernel only works when field is a W3 field at lowest
!!          order, since it is assumed that ndf_w3 = 1.


module tl_ffsl_flux_z_constant_kernel_mod

use argument_mod,                   only : arg_type,              &
                                           GH_FIELD, GH_REAL,     &
                                           GH_READ, GH_WRITE,     &
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
type, public, extends(kernel_type) :: tl_ffsl_flux_z_constant_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                  &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! flux_pert
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! flux_ls
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! frac_wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! frac_wind_pert
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! ls_field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! detj
       arg_type(GH_SCALAR, GH_REAL,    GH_READ)        & ! dt
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: tl_ffsl_flux_z_constant_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: tl_ffsl_flux_z_constant_code

contains

!> @brief Computes the mass flux for FFSL using PCM in the z direction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] flux_pert      The pert field ls wind flux to be computed
!> @param[in,out] flux_ls        The ls field pert wind flux to be computed
!> @param[in]     frac_wind      The fractional vertical ls wind
!> @param[in]     dep_dist       The vertical ls departure points
!> @param[in]     frac_wind_pert The fractional vertical perturbation wind
!> @param[in]     field          The field to construct the first flux
!> @param[in]     ls_field       The ls_field to construct the second flux
!> @param[in]     detj           Volume of cells
!> @param[in]     dt             Time step
!> @param[in]     ndf_w2v        Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v       Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v        The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3         Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3        Number of unique degrees of freedom for W3
!> @param[in]     map_w3         The dofmap for the cell at the base of the column
subroutine tl_ffsl_flux_z_constant_code( nlayers,        &
                                         flux_pert,      &
                                         flux_ls,        &
                                         frac_wind,      &
                                         dep_dist,       &
                                         frac_wind_pert, &
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
  real(kind=r_tran),   intent(in)    :: field(undf_w3)
  real(kind=r_tran),   intent(in)    :: ls_field(undf_w3)
  real(kind=r_tran),   intent(in)    :: frac_wind(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
  real(kind=r_tran),   intent(in)    :: frac_wind_pert(undf_w2v)
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

  ! ========================================================================== !
  ! Build fluxes
  ! ========================================================================== !

  ! Force bottom flux to be zero
  flux_pert(w2v_idx) = 0.0_r_tran
  flux_ls(w2v_idx) = 0.0_r_tran

  ! Loop through faces
  do k = 1, nlayers - 1

    ! Pull out departure point, and separate into integer / frac parts
    displacement = dep_dist(w2v_idx + k)
    int_displacement = INT(displacement, i_def)
    frac_dist = displacement - REAL(int_displacement, r_tran)
    sign_displacement = INT(SIGN(1.0_r_tran, displacement))

    ! Set an offset for the stencil index, based on dep point sign
    sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

    ! Determine departure cell
    dep_cell_idx = k - int_displacement + sign_offset - 1

    ! ======================================================================== !
    ! Integer sum
    mass_whole_cells = 0.0_r_tran
    lowest_whole_cell = MIN(dep_cell_idx + 1, k)
    highest_whole_cell = MAX(dep_cell_idx, k) - 1
    do i = lowest_whole_cell, highest_whole_cell
      mass_whole_cells = mass_whole_cells + field(w3_idx + i) * detj(w3_idx + i)
    end do

    ! ======================================================================== !
    ! Compute flux
    ! Fractional part of reconstruction is just frac_wind multiplying the field
    ! in the departure cell
    flux_pert(w2v_idx + k) =                                                   &
                (frac_wind(w2v_idx + k) * field(w3_idx + dep_cell_idx)         &
                + sign(1.0_r_tran, displacement) * mass_whole_cells) / dt
    flux_ls(w2v_idx + k) =                                                     &
               (frac_wind_pert(w2v_idx + k) * ls_field(w3_idx + dep_cell_idx)  &
               ) / dt

  end do

  ! Force top flux to be zero
  flux_pert(w2v_idx + nlayers) = 0.0_r_tran
  flux_ls(w2v_idx + nlayers) = 0.0_r_tran

end subroutine tl_ffsl_flux_z_constant_code

end module tl_ffsl_flux_z_constant_kernel_mod
