!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes the vertical mass flux using Nirvana.
!> @details This kernel reconstructs a field using the reversible Nirvana scheme
!!          which is equivalent to fitting a quadratic to the cell such that the
!!          integral of the quadratic equals the integral of the field in the
!!          cell, and the quadratic matches the gradient of the field at
!!          cell edges. A limiter can be applied to ensure monotonicity.
!!          This kernel is designed to work in the vertical direction only and
!!          takes into account the vertical boundaries and grid spacing.
!!
!!          Note that this kernel only works when field is a W3 field at lowest
!!          order, since it is assumed that ndf_w3 = 1.

module ffsl_flux_z_nirvana_kernel_mod

use argument_mod,                   only : arg_type,              &
                                           GH_FIELD, GH_REAL,     &
                                           GH_READ, GH_WRITE,     &
                                           GH_SCALAR, GH_INTEGER, &
                                           GH_LOGICAL, CELL_COLUMN
use fs_continuity_mod,              only : W3, W2v
use constants_mod,                  only : r_tran, i_def, l_def, EPS_R_TRAN
use kernel_mod,                     only : kernel_type
use transport_enumerated_types_mod, only : vertical_monotone_relaxed, &
                                           vertical_monotone_strict,  &
                                           vertical_monotone_positive

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ffsl_flux_z_nirvana_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/                  &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! flux
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! frac_wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! dz
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! detj
       arg_type(GH_SCALAR, GH_REAL,    GH_READ),       & ! dt
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)        & ! monotone
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: ffsl_flux_z_nirvana_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: ffsl_flux_z_nirvana_code

contains

!> @brief Compute the flux using the Nirvana reconstruction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] flux      The flux to be computed
!> @param[in]     frac_wind The fractional vertical wind
!> @param[in]     dep_dist  The vertical departure points
!> @param[in]     field     The field to construct the flux
!> @param[in]     dz        Vertical length of the W3 cell
!> @param[in]     detj      Volume of cells
!> @param[in]     dt        Time step
!> @param[in]     monotone  Monotonicity option to use
!> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine ffsl_flux_z_nirvana_code( nlayers,   &
                                     flux,      &
                                     frac_wind, &
                                     dep_dist,  &
                                     field,     &
                                     dz,        &
                                     detj,      &
                                     dt,        &
                                     monotone,  &
                                     ndf_w2v,   &
                                     undf_w2v,  &
                                     map_w2v,   &
                                     ndf_w3,    &
                                     undf_w3,   &
                                     map_w3 )

  use subgrid_vertical_support_mod, only: vertical_nirvana_recon,              &
                                          vertical_nirvana_mono_relax,         &
                                          vertical_nirvana_mono_strict,        &
                                          vertical_nirvana_positive,           &
                                          second_order_vertical_gradient

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: ndf_w3
  real(kind=r_tran),   intent(inout) :: flux(undf_w2v)
  real(kind=r_tran),   intent(in)    :: field(undf_w3)
  real(kind=r_tran),   intent(in)    :: frac_wind(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dz(undf_w3)
  real(kind=r_tran),   intent(in)    :: detj(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  real(kind=r_tran),   intent(in)    :: dt
  integer(kind=i_def), intent(in)    :: monotone

  ! Internal variables
  integer(kind=i_def) :: k, i, w2v_idx, w3_idx
  integer(kind=i_def) :: int_displacement, sign_displacement
  integer(kind=i_def) :: dep_cell_idx, sign_offset
  integer(kind=i_def) :: lowest_whole_cell, highest_whole_cell
  integer(kind=i_def) :: lower_cell, higher_cell

  real(kind=r_tran)   :: displacement, frac_dist
  real(kind=r_tran)   :: reconstruction, mass_whole_cells
  real(kind=r_tran)   :: edge_gradient(0:nlayers)
  real(kind=r_tran)   :: field_mono_stencil(3)

  w3_idx = map_w3(1)
  w2v_idx = map_w2v(1)

  ! ========================================================================== !
  ! Reconstruct edge values
  ! ========================================================================== !

  ! Gradient is zero at the bottom
  edge_gradient(0) = 0.0_r_tran

  ! Build edge gradients on internal layers
  do k = 1, nlayers - 1
    call second_order_vertical_gradient(field(w3_idx+k-1 : w3_idx+k),          &
                                        dz(w3_idx+k-1 : w3_idx+k),             &
                                        edge_gradient(k))
  end do

  ! Gradient at the top is zero
  edge_gradient(nlayers) = 0.0_r_tran

  ! ========================================================================== !
  ! Build fluxes
  ! ========================================================================== !

  ! Force bottom flux to be zero
  flux(w2v_idx) = 0.0_r_tran

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
    ! Perform reversible PPM reconstruction for fractional part
    call vertical_nirvana_recon(reconstruction, frac_dist,       &
                                field(w3_idx + dep_cell_idx),    & ! field in dep cell
                                edge_gradient(dep_cell_idx),     & ! edge below dep cell
                                edge_gradient(dep_cell_idx + 1), & ! edge above dep cell
                                dz(w3_idx + dep_cell_idx))

    ! ======================================================================== !
    ! Apply monotonicity for fractional part
    select case ( monotone )
    case ( vertical_monotone_strict )
      lower_cell = MAX(0, dep_cell_idx - 1)
      higher_cell = MIN(nlayers - 1, dep_cell_idx + 1)
      field_mono_stencil(1) = field(w3_idx + lower_cell)
      field_mono_stencil(2) = field(w3_idx + dep_cell_idx)
      field_mono_stencil(3) = field(w3_idx + higher_cell)
      call vertical_nirvana_mono_strict(reconstruction,                  &
                                        field_mono_stencil,              &
                                        edge_gradient(dep_cell_idx),     &
                                        edge_gradient(dep_cell_idx + 1), &
                                        dz(w3_idx + dep_cell_idx))

    case ( vertical_monotone_relaxed )
      lower_cell = MAX(0, dep_cell_idx - 1)
      higher_cell = MIN(nlayers - 1, dep_cell_idx + 1)
      field_mono_stencil(1) = field(w3_idx + lower_cell)
      field_mono_stencil(2) = field(w3_idx + dep_cell_idx)
      field_mono_stencil(3) = field(w3_idx + higher_cell)
      call vertical_nirvana_mono_relax(reconstruction, frac_dist,       &
                                       field_mono_stencil,              &
                                       edge_gradient(dep_cell_idx),     &
                                       edge_gradient(dep_cell_idx + 1), &
                                       dz(w3_idx + dep_cell_idx))

    case ( vertical_monotone_positive )
      call vertical_nirvana_positive(reconstruction, frac_dist,       &
                                     field(w3_idx + dep_cell_idx),    &
                                     edge_gradient(dep_cell_idx),     &
                                     edge_gradient(dep_cell_idx + 1), &
                                     dz(w3_idx + dep_cell_idx))

    end select

    ! ======================================================================== !
    ! Compute flux
    flux(w2v_idx + k) = (frac_wind(w2v_idx + k) * reconstruction               &
                         + sign(1.0_r_tran, displacement) * mass_whole_cells) / dt

  end do

  ! Force top flux to be zero
  flux(w2v_idx + nlayers) = 0.0_r_tran

end subroutine ffsl_flux_z_nirvana_code

end module ffsl_flux_z_nirvana_kernel_mod
