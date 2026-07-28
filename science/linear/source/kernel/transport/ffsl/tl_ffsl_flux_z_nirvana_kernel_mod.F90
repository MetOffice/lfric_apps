!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes the vertical mass flux using Nirvana for the tangent linear model.
!> @details This kernel computes both parts of the vertical flux using the Nirvana
!!          scheme for tangent linear model transport. The first part of the flux is
!!          the reconstruction of the perturbation field with the ls wind. The second
!!          part uses a constant reconstruction of the ls field in the ls wind
!!          departure cell, but is multiplied by the fractional perturbation wind.
!!          The Nirvana scheme is equivalent to fitting a quadratic to the cell such
!!          that the integral of the quadratic equals the integral of the field in
!!          the cell, and the quadratic matches the gradient of the field at cell edges.
!!          This kernel is designed to work in the vertical direction only and takes
!!          into account the vertical boundaries and grid spacing.
!!
!!          Note that this kernel only works when field is a W3 field at lowest
!!          order, since it is assumed that ndf_w3 = 1.

module tl_ffsl_flux_z_nirvana_kernel_mod

use argument_mod,                   only : arg_type,              &
                                           GH_FIELD, GH_REAL,     &
                                           GH_READ, GH_WRITE,     &
                                           GH_SCALAR, GH_LOGICAL, &
                                           CELL_COLUMN
use fs_continuity_mod,              only : W3, W2v
use constants_mod,                  only : r_tran, i_def, l_def, EPS_R_TRAN
use kernel_mod,                     only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: tl_ffsl_flux_z_nirvana_kernel_type
  private
  type(arg_type) :: meta_args(13) = (/                  &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2v), & ! flux_pert
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2v), & ! flux_ls
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2v), & ! frac_wind
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2v), & ! frac_wind_pert
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),  & ! field
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),  & ! ls_field
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  W3),  & ! dla_dz
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  W3),  & ! dlb_dz
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),  & ! dz
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),  & ! detj
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       & ! dt
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)        & ! log_space
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: tl_ffsl_flux_z_nirvana_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: tl_ffsl_flux_z_nirvana_code

contains

!> @brief Compute the flux using the Nirvana reconstruction.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] flux_pert      The pert field ls wind flux to be computed
!> @param[in,out] flux_ls        The ls field pert wind flux to be computed
!> @param[in]     frac_wind      The fractional vertical ls wind
!> @param[in]     dep_dist       The vertical ls departure points
!> @param[in]     frac_wind_pert The fractional vertical perturbation wind
!> @param[in]     field          The perturbation field to construct the first flux
!> @param[in]     ls_field       The ls_field to construct the second flux
!> @param[in]     dla_dz_1       First set of coefficients for the reconstruction of
!!                               the field at the edge above the cell
!> @param[in]     dla_dz_2       Second set of coefficients for the reconstruction of
!!                               the field at the edge above the cell
!> @param[in]     dla_dz_3       Third set of coefficients for the reconstruction of
!!                               the field at the edge above the cell
!> @param[in]     dlb_dz_1       First set of coefficients for the reconstruction of
!!                               the field at the edge below the cell
!> @param[in]     dlb_dz_2       Second set of coefficients for the reconstruction of
!!                               the field at the edge below the cell
!> @param[in]     dlb_dz_3       Third set of coefficients for the reconstruction of
!!                               the field at the edge below the cell
!> @param[in]     dz             Vertical length of the W3 cell
!> @param[in]     detj           Volume of cells
!> @param[in]     dt             Time step
!> @param[in]     log_space      Switch to use natural logarithmic space
!!                               for edge interpolation
!> @param[in]     ndf_w2v        Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v       Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v        The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3         Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3        Number of unique degrees of freedom for W3
!> @param[in]     map_w3         The dofmap for the cell at the base of the column
subroutine tl_ffsl_flux_z_nirvana_code( nlayers,         &
                                        flux_pert,       &
                                        flux_ls,         &
                                        frac_wind,       &
                                        dep_dist,        &
                                        frac_wind_pert,  &
                                        field,           &
                                        ls_field,        &
                                        dla_dz_1,        &
                                        dla_dz_2,        &
                                        dla_dz_3,        &
                                        dlb_dz_1,        &
                                        dlb_dz_2,        &
                                        dlb_dz_3,        &
                                        dz,              &
                                        detj,            &
                                        dt,              &
                                        log_space,       &
                                        ndf_w2v,         &
                                        undf_w2v,        &
                                        map_w2v,         &
                                        ndf_w3,          &
                                        undf_w3,         &
                                        map_w3 )

  use subgrid_vertical_support_mod,   only: third_order_vertical_edge
  use subgrid_common_support_mod,     only: subgrid_quadratic_recon
  use transport_enumerated_types_mod, only: monotone_none


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
  real(kind=r_tran),   intent(in)    :: frac_wind_pert(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dla_dz_1(undf_w3)
  real(kind=r_tran),   intent(in)    :: dla_dz_2(undf_w3)
  real(kind=r_tran),   intent(in)    :: dla_dz_3(undf_w3)
  real(kind=r_tran),   intent(in)    :: dlb_dz_1(undf_w3)
  real(kind=r_tran),   intent(in)    :: dlb_dz_2(undf_w3)
  real(kind=r_tran),   intent(in)    :: dlb_dz_3(undf_w3)
  real(kind=r_tran),   intent(in)    :: dz(undf_w3)
  real(kind=r_tran),   intent(in)    :: detj(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  real(kind=r_tran),   intent(in)    :: dt
  logical(kind=l_def), intent(in)    :: log_space

  ! Local arrays
  integer(kind=i_def) :: sign_displacement(nlayers-1)
  integer(kind=i_def) :: sign_offset(nlayers-1)
  integer(kind=i_def) :: dep_cell_idx(nlayers-1)
  real(kind=r_tran)   :: reconstruction(nlayers-1)
  real(kind=r_tran)   :: edge_above(0:nlayers-1)
  real(kind=r_tran)   :: edge_below(0:nlayers-1)
  real(kind=r_tran)   :: edge_left(nlayers-1)
  real(kind=r_tran)   :: edge_right(nlayers-1)
  real(kind=r_tran)   :: field_local(nlayers-1)
  real(kind=r_tran)   :: ls_reconstruction(nlayers-1)
  real(kind=r_tran)   :: displacement(nlayers-1)
  real(kind=r_tran)   :: frac_dist(nlayers-1)
  real(kind=r_tran)   :: mass(nlayers)

  ! Local scalars
  integer(kind=i_def) :: k, w2v_idx, w3_idx
  integer(kind=i_def) :: b_idx, t_idx
  integer(kind=i_def) :: lowest_whole_cell, highest_whole_cell
  real(kind=r_tran)   :: inv_dt

  w3_idx = map_w3(1)
  w2v_idx = map_w2v(1)
  inv_dt = 1.0_r_tran / dt

  ! ========================================================================== !
  ! EDGE RECONSTRUCTION
  ! ========================================================================== !
  ! Use Nirvana reconstruction for interpolated edge values, computing
  ! edge_above and edge_below from the perturbation field
  call third_order_vertical_edge(                                              &
          field(w3_idx : w3_idx+nlayers-1),                                    &
          dla_dz_1(w3_idx : w3_idx+nlayers-1),                                 &
          dla_dz_2(w3_idx : w3_idx+nlayers-1),                                 &
          dla_dz_3(w3_idx : w3_idx+nlayers-1),                                 &
          dlb_dz_1(w3_idx : w3_idx+nlayers-1),                                 &
          dlb_dz_2(w3_idx : w3_idx+nlayers-1),                                 &
          dlb_dz_3(w3_idx : w3_idx+nlayers-1),                                 &
          dz(w3_idx : w3_idx+nlayers-1),                                       &
          edge_above, edge_below, log_space, nlayers                           &
  )

  ! ========================================================================== !
  ! Extract departure info
  ! ========================================================================== !
  b_idx = w2v_idx + 1
  t_idx = w2v_idx + nlayers - 1

  ! Pull out departure point, and separate into integer / frac parts
  displacement(:) = dep_dist(b_idx : t_idx)
  frac_dist(:) = ABS(displacement(:) - REAL(INT(displacement(:), i_def), r_tran))
  sign_displacement(:) = INT(SIGN(1.0_r_tran, displacement(:)))

  ! Set an offset for the stencil index, based on dep point sign
  sign_offset(:) = (1 - sign_displacement(:)) / 2   ! 0 if sign == 1, 1 if sign == -1

  ! Determine departure cell
  do k = 1, nlayers - 1
    dep_cell_idx(k) = k - INT(displacement(k), i_def) + sign_offset(k) - 1
  end do

  ! ========================================================================== !
  ! Populate local arrays for fractional flux calculations
  ! ========================================================================== !

  do k = 1, nlayers - 1
    field_local(k) = field(w3_idx + dep_cell_idx(k))
    ! Constant ls_field reconstruction is just the upwind cell
    ls_reconstruction(k) = ls_field(w3_idx + dep_cell_idx(k))

    ! Edge values on left/right of cell needed for subgrid reconstruction
    ! Left value is the *edge* upwind of the departure cell
    ! Right value is the *edge* downwind of the departure cell
    edge_left(k) = (                                                           &
      sign_offset(k) * edge_above(dep_cell_idx(k))                             &
      + (1 - sign_offset(k)) * edge_below(dep_cell_idx(k))                     &
    )
    edge_right(k) = (                                                          &
      (1 - sign_offset(k)) * edge_above(dep_cell_idx(k))                       &
      + sign_offset(k) * edge_below(dep_cell_idx(k))                           &
    )
  end do

  ! ========================================================================== !
  ! FRACTIONAL FIRST FLUX RECONSTRUCTION
  ! ========================================================================== !

  ! Compute reconstruction using field edge values
  ! and quadratic subgrid reconstruction (note, set order to 1 to get correct
  ! size of field_local in subgrid_quadratic_recon)

  call subgrid_quadratic_recon(                                                &
          reconstruction, frac_dist, field_local,                              &
          edge_left, edge_right, monotone_none, 0, nlayers-1                   &
  )

  ! ========================================================================== !
  ! INTEGER FIRST FLUX
  ! ========================================================================== !
  ! Set fluxes to be zero initially
  flux_pert(w2v_idx : w2v_idx + nlayers) = 0.0_r_tran
  flux_ls(w2v_idx : w2v_idx + nlayers)   = 0.0_r_tran

  ! Pre-multiply field and volume to get mass, for use in integer flux
  mass(:) = field(w3_idx : w3_idx+nlayers-1) * detj(w3_idx : w3_idx+nlayers-1)

  ! Integer sum
  do k = 1, nlayers - 1
    lowest_whole_cell = MIN(dep_cell_idx(k) + 1, k) + 1
    highest_whole_cell = MAX(dep_cell_idx(k), k)
    flux_pert(w2v_idx + k) = SUM(mass(lowest_whole_cell : highest_whole_cell))
  end do

  ! ========================================================================== !
  ! Assign fluxes
  ! ========================================================================== !
  b_idx = w2v_idx + 1
  t_idx = w2v_idx + nlayers - 1

  flux_pert(b_idx : t_idx) = inv_dt * (                                        &
      frac_wind(b_idx : t_idx) * reconstruction(:)                             &
      + sign_displacement(:) * flux_pert(b_idx : t_idx)                        &
  )

  flux_ls(b_idx : t_idx) = inv_dt * (                                          &
      frac_wind_pert(b_idx : t_idx) * ls_reconstruction(:)                     &
  )

end subroutine tl_ffsl_flux_z_nirvana_code

end module tl_ffsl_flux_z_nirvana_kernel_mod
