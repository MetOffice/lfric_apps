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

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ffsl_flux_z_nirvana_kernel_type
  private
  type(arg_type) :: meta_args(13) = (/                  &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2v), & ! flux
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2v), & ! frac_wind
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),  & ! field
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  W3),  & ! dla_dz
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  W3),  & ! dlb_dz
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),  & ! dz
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W3),  & ! detj
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       & ! dt
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),       & ! monotone
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ),       & ! min_val
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),       & ! log_space
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)        & ! monotone_above
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
!> @param[in]     dla_dz_1  First set of coefficients for the reconstruction of
!!                          the field at the edge above the cell
!> @param[in]     dla_dz_2  Second set of coefficients for the reconstruction of
!!                          the field at the edge above the cell
!> @param[in]     dla_dz_3  Third set of coefficients for the reconstruction of
!!                          the field at the edge above the cell
!> @param[in]     dlb_dz_1  First set of coefficients for the reconstruction of
!!                          the field at the edge below the cell
!> @param[in]     dlb_dz_2  Second set of coefficients for the reconstruction of
!!                          the field at the edge below the cell
!> @param[in]     dlb_dz_3  Third set of coefficients for the reconstruction of
!!                          the field at the edge below the cell
!> @param[in]     dz        Vertical length of the W3 cell
!> @param[in]     detj      Volume of cells
!> @param[in]     dt        Time step
!> @param[in]     monotone  Monotonicity option to use
!> @param[in]     min_val   Minimum value to enforce when using
!!                          quasi-monotone limiter for PPM
!> @param[in]     log_space Switch to use natural logarithmic space
!!                          for edge interpolation
!> @param[in]     monotone_above
!!                          Levels at and above which to apply the monotonicity
!!                          limiter. Below this, apply no monotonicity
!> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine ffsl_flux_z_nirvana_code( nlayers,         &
                                     flux,            &
                                     frac_wind,       &
                                     dep_dist,        &
                                     field,           &
                                     dla_dz_1,        &
                                     dla_dz_2,        &
                                     dla_dz_3,        &
                                     dlb_dz_1,        &
                                     dlb_dz_2,        &
                                     dlb_dz_3,        &
                                     dz,              &
                                     detj,            &
                                     dt,              &
                                     monotone,        &
                                     min_val,         &
                                     log_space,       &
                                     monotone_above,  &
                                     ndf_w2v,         &
                                     undf_w2v,        &
                                     map_w2v,         &
                                     ndf_w3,          &
                                     undf_w3,         &
                                     map_w3 )

  use subgrid_vertical_support_mod,   only: third_order_vertical_edge
  use subgrid_common_support_mod,     only: monotonic_edge,                    &
                                            subgrid_quadratic_recon
  use transport_enumerated_types_mod, only: monotone_none,                     &
                                            monotone_positive

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
  integer(kind=i_def), intent(in)    :: monotone
  real(kind=r_tran),   intent(in)    :: min_val
  logical(kind=l_def), intent(in)    :: log_space
  integer(kind=i_def), intent(in)    :: monotone_above

  ! Local arrays
  integer(kind=i_def) :: sign_displacement(nlayers-1)
  integer(kind=i_def) :: sign_offset(nlayers-1)
  integer(kind=i_def) :: dep_cell_idx(nlayers-1)
  real(kind=r_tran)   :: reconstruction(nlayers-1)
  real(kind=r_tran)   :: edge_above(0:nlayers-1)
  real(kind=r_tran)   :: edge_below(0:nlayers-1)
  real(kind=r_tran)   :: edge_left(nlayers-1)
  real(kind=r_tran)   :: edge_right(nlayers-1)
  real(kind=r_tran)   :: field_local(nlayers-1, 3)
  real(kind=r_tran)   :: displacement(nlayers-1)
  real(kind=r_tran)   :: frac_dist(nlayers-1)
  real(kind=r_tran)   :: mass(nlayers)
  real(kind=r_tran)   :: field_local_lower(MAX(monotone_above-1,1), 3)
  real(kind=r_tran)   :: field_local_upper(MAX(nlayers-monotone_above,1), 3)

  ! Local scalars
  integer(kind=i_def) :: k, w2v_idx, w3_idx
  integer(kind=i_def) :: b_idx, t_idx, l_idx, r_idx, array_length
  integer(kind=i_def) :: lowest_whole_cell, highest_whole_cell
  real(kind=r_tran)   :: inv_dt

  w3_idx = map_w3(1)
  w2v_idx = map_w2v(1)
  inv_dt = 1.0_r_tran / dt

  ! ========================================================================== !
  ! EDGE RECONSTRUCTION
  ! ========================================================================== !
  ! Use Nirvana reconstruction for interpolated edge values
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
  if (monotone == monotone_none .or. monotone == monotone_positive) then
    ! Don't need the neighbouring cell values
    do k = 1, nlayers - 1
      field_local(k,2) = field(w3_idx + dep_cell_idx(k))

      ! Edge values on left/right of cell needed for subgrid reconstruction
      ! Left value is the *edge* upwind of the departure cell
      ! Right value is the *edge* downwind of the departure cell
      edge_left(k) = (                                                         &
        sign_offset(k) * edge_above(dep_cell_idx(k))                           &
        + (1 - sign_offset(k)) * edge_below(dep_cell_idx(k))                   &
      )
      edge_right(k) = (                                                        &
        (1 - sign_offset(k)) * edge_above(dep_cell_idx(k))                     &
        + sign_offset(k) * edge_below(dep_cell_idx(k))                         &
      )
    end do

  else
    ! When populating arrays, include neighbouring cell values
    do k = 1, nlayers - 1
      ! Fields on left/right of cell needed for bounding edge values
      ! Left value is the *cell* upwind of the departure cell
      ! Right value is the *cell* downwind of the departure cell
      ! At the top/bottom of the domain, use the same cell as the neighbour
      l_idx = MIN(MAX(dep_cell_idx(k) - sign_displacement(k), 0), nlayers-1)
      r_idx = MIN(MAX(dep_cell_idx(k) + sign_displacement(k), 0), nlayers-1)
      field_local(k,1) = field(w3_idx + l_idx)
      field_local(k,2) = field(w3_idx + dep_cell_idx(k))
      field_local(k,3) = field(w3_idx + r_idx)

      ! Edge values on left/right of cell needed for subgrid reconstruction
      ! Left value is the *edge* upwind of the departure cell
      ! Right value is the *edge* downwind of the departure cell
      edge_left(k) = (                                                         &
        sign_offset(k) * edge_above(dep_cell_idx(k))                           &
        + (1 - sign_offset(k)) * edge_below(dep_cell_idx(k))                   &
      )
      edge_right(k) = (                                                        &
        (1 - sign_offset(k)) * edge_above(dep_cell_idx(k))                     &
        + sign_offset(k) * edge_below(dep_cell_idx(k))                         &
      )
    end do
  end if

  ! ========================================================================== !
  ! FRACTIONAL FLUX RECONSTRUCTION
  ! ========================================================================== !
  ! Apply monotonicity to edges if required
  call monotonic_edge(                                                         &
          field_local, monotone, min_val, edge_left, edge_right,               &
          1, monotone_above, nlayers-1                                         &
  )

  ! Compute reconstruction using field edge values
  ! and quadratic subgrid reconstruction

  if (monotone /= monotone_none .and. monotone_above > 1) then
    ! Different monotonicity for different levels, so call this function twice
    ! Application without monotonicity
    b_idx = 1
    t_idx = MIN(monotone_above-1, nlayers-1)
    array_length = t_idx - b_idx + 1
    field_local_lower(:,:) = field_local(b_idx:t_idx, :)
    call subgrid_quadratic_recon(                                              &
            reconstruction(b_idx:t_idx), frac_dist(b_idx:t_idx),               &
            field_local_lower, edge_left(b_idx:t_idx),                         &
            edge_right(b_idx:t_idx), monotone_none, 1, array_length            &
    )

    ! Application with monotonicity
    if (monotone_above < nlayers) then
      b_idx = monotone_above
      t_idx = nlayers-1
      array_length = t_idx - b_idx + 1
      field_local_upper(:,:) = field_local(b_idx:t_idx, :)
      call subgrid_quadratic_recon(                                            &
              reconstruction(b_idx:t_idx), frac_dist(b_idx:t_idx),             &
              field_local_upper, edge_left(b_idx:t_idx),                       &
              edge_right(b_idx:t_idx), monotone, 1, array_length               &
      )
    end if
  else
    call subgrid_quadratic_recon(                                              &
            reconstruction, frac_dist, field_local,                            &
            edge_left, edge_right, monotone, 1, nlayers-1                      &
    )
  end if

  ! ========================================================================== !
  ! INTEGER FLUX
  ! ========================================================================== !
  ! Set flux to be zero initially
  flux(w2v_idx : w2v_idx + nlayers) = 0.0_r_tran

  ! Pre-multiply field and volume to get mass, for use in integer flux
  mass(:) = field(w3_idx : w3_idx+nlayers-1) * detj(w3_idx : w3_idx+nlayers-1)

  ! Integer sum
  do k = 1, nlayers - 1
    lowest_whole_cell = MIN(dep_cell_idx(k) + 1, k) + 1
    highest_whole_cell = MAX(dep_cell_idx(k), k)
    flux(w2v_idx + k) = SUM(mass(lowest_whole_cell : highest_whole_cell))
  end do

  ! ========================================================================== !
  ! Assign flux
  ! ========================================================================== !
  b_idx = w2v_idx + 1
  t_idx = w2v_idx + nlayers - 1

  flux(b_idx : t_idx) = inv_dt * (                                             &
      frac_wind(b_idx : t_idx) * reconstruction(:)                             &
      + sign_displacement(:) * flux(b_idx : t_idx)                             &
  )

end subroutine ffsl_flux_z_nirvana_code

end module ffsl_flux_z_nirvana_kernel_mod
