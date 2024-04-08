!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to compute the horizontal departure distances for FFSL
!> @details This code calculates the distance which is swept through a cell
!!          in one dimension during one time step of consistent transport. The
!!          arrival point is the cell face and the departure point is calculated
!!          and stored as a field.
!!          The departure points are a dimensionless displacement, corresponding
!!          to the number of cells moved by a fluid parcel. The fractional wind
!!          or dry flux is also computed.

module hori_dep_dist_ffsl_kernel_mod

use argument_mod,                only : arg_type,                  &
                                        GH_FIELD, GH_REAL,         &
                                        GH_WRITE, GH_READ,         &
                                        GH_SCALAR, GH_INTEGER,     &
                                        STENCIL, X1D, Y1D,         &
                                        GH_LOGICAL, CELL_COLUMN,   &
                                        ANY_DISCONTINUOUS_SPACE_2, &
                                        ANY_DISCONTINUOUS_SPACE_3
use fs_continuity_mod,           only : W3, W2h
use constants_mod,               only : r_tran, i_def, l_def
use kernel_mod,                  only : kernel_type
use log_mod,                     only : log_event,         &
                                        log_scratch_space, &
                                        LOG_LEVEL_ERROR
use reference_element_mod,       only : E, W, N, S

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: hori_dep_dist_ffsl_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                       &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_2), & ! frac_dry_flux
       arg_type(GH_FIELD,  GH_REAL, GH_READ,    ANY_DISCONTINUOUS_SPACE_2), & ! dry_flux
       arg_type(GH_FIELD,  GH_REAL, GH_READ,    W3, STENCIL(X1D)),          & ! dry_mass_for_x
       arg_type(GH_FIELD,  GH_REAL, GH_READ,    W3, STENCIL(Y1D)),          & ! dry_mass_for_y
       arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_3), & ! face_selector ew
       arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_3), & ! face_selector ns
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                            & ! stencil_extent
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                             & ! cap_dep_points
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: hori_dep_dist_ffsl_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: hori_dep_dist_ffsl_code

contains

!> @brief Kernel to compute the horizontal departure distances for FFSL.
!> @param[in]     nlayers           The number of layers in the mesh
!> @param[in,out] dep_dist          Field with dep distances
!> @param[in,out] frac_dry_flux     Fractional part of the dry flux in the local
!!                                  x direction, to be computed here
!> @param[in]     dry_flux          The mass flux used in a one-dimensional
!!                                  transport step for the dry density
!> @param[in]     dry_mass_for_x    The dry mass field before the y-direction
!!                                  transport step
!> @param[in]     stencil_size_x    Size of the X1D stencil for the mass field
!> @param[in]     stencil_map_x     Map of DoFs in the stencil for mass field
!> @param[in]     dry_mass_for_y    The dry mass field before the x-direction
!!                                  transport step
!> @param[in]     stencil_size_y    Size of the Y1D stencil for the mass field
!> @param[in]     stencil_map_y     Map of DoFs in the stencil for mass field
!> @param[in]     stencil_size      Size of the stencil for the mass field
!> @param[in]     stencil_map       Map of DoFs in the stencil for mass field
!> @param[in]     face_selector_ew  2D field indicating which faces to loop over in x
!> @param[in]     face_selector_ns  2D field indicating which faces to loop over in y
!> @param[in]     stencil_extent    Max number of stencil cells in one direction
!> @param[in]     cap_dep_points    Flag for whether departure points should be
!!                                  capped if they exceed the stencil depth
!> @param[in]     ndf_w2h           Number of DoFs per cell for W2H
!> @param[in]     undf_w2h          Num of W2H DoFs in memory for this partition
!> @param[in]     map_w2h           Map of lowest-cell W2H DoFs
!> @param[in]     ndf_w3            Number of DoFs per cell for W3
!> @param[in]     undf_w3           Num of W3 DoFs in memory for this partition
!> @param[in]     map_w3            Map of lowest-cell W3 DoFs
!> @param[in]     ndf_w3_2d         Number of DoFs for 2D W3 per cell
!> @param[in]     undf_w3_2d        Number of DoFs for this partition for 2D W3
!> @param[in]     map_w3_2d         Map for 2D W3
subroutine hori_dep_dist_ffsl_code( nlayers,             &
                                    dep_dist,            &
                                    frac_dry_flux,       &
                                    dry_flux,            &
                                    dry_mass_for_x,      &
                                    stencil_size_x,      &
                                    stencil_map_x,       &
                                    dry_mass_for_y,      &
                                    stencil_size_y,      &
                                    stencil_map_y,       &
                                    face_selector_ew,    &
                                    face_selector_ns,    &
                                    stencil_extent,      &
                                    cap_dep_points,      &
                                    ndf_w2h,             &
                                    undf_w2h,            &
                                    map_w2h,             &
                                    ndf_w3,              &
                                    undf_w3,             &
                                    map_w3,              &
                                    ndf_w3_2d,           &
                                    undf_w3_2d,          &
                                    map_w3_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  logical(kind=l_def), intent(in)    :: cap_dep_points
  integer(kind=i_def), intent(in)    :: ndf_w2h, ndf_w3, ndf_w3_2d
  integer(kind=i_def), intent(in)    :: undf_w2h, undf_w3, undf_w3_2d
  integer(kind=i_def), intent(in)    :: stencil_size_x
  integer(kind=i_def), intent(in)    :: stencil_size_y
  integer(kind=i_def), intent(in)    :: stencil_extent
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w3_2d(ndf_w3_2d)
  integer(kind=i_def), intent(in)    :: stencil_map_x(ndf_w3, stencil_size_x)
  integer(kind=i_def), intent(in)    :: stencil_map_y(ndf_w3, stencil_size_y)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)
  real(kind=r_tran),   intent(in)    :: dry_mass_for_x(undf_w3)
  real(kind=r_tran),   intent(in)    :: dry_mass_for_y(undf_w3)
  real(kind=r_tran),   intent(in)    :: dry_flux(undf_w2h)
  real(kind=r_tran),   intent(inout) :: dep_dist(undf_w2h)
  real(kind=r_tran),   intent(inout) :: frac_dry_flux(undf_w2h)

  integer(kind=i_def) :: df, df_idx, idx_w2h, j, k, num_int_cells
  integer(kind=i_def) :: sign_flux, sign_offset, dof_offset
  integer(kind=i_def) :: stencil_idx, col_idx, rel_idx
  real(kind=r_tran)   :: flux_face_i, running_flux_i, flux_cell_j
  real(kind=r_tran)   :: frac_dep_dist, frac_flux_face_i
  integer(kind=i_def) :: stencil_size, stencil_half, lam_edge_size
  integer(kind=i_def) :: local_dofs_x(2), local_dofs_y(2)

  ! x-direction
  local_dofs_x = (/ W, E /)
  local_dofs_y = (/ S, N /)

  ! Set stencil info -----------------------------------------------------------
  stencil_size = stencil_size_x
  stencil_half = (stencil_size + 1_i_def) / 2_i_def
  lam_edge_size = 2_i_def*stencil_extent+1_i_def

  ! Calculation near LAM boundaries --------------------------------------------
  if (lam_edge_size > stencil_size) then
    ! Set output to zero
    do df_idx = 1, face_selector_ew(map_w3_2d(1))
      df = local_dofs_x(df_idx)
      do k = 0, nlayers - 1
        frac_dry_flux(map_w2h(df) + k) = 0.0_r_tran
        dep_dist(map_w2h(df) + k) = 0.0_r_tran
      end do
    end do
    do df_idx = 1, face_selector_ns(map_w3_2d(1))
      df = local_dofs_y(df_idx)
      do k = 0, nlayers - 1
        frac_dry_flux(map_w2h(df) + k) = 0.0_r_tran
        dep_dist(map_w2h(df) + k) = 0.0_r_tran
      end do
    end do

  ! ========================================================================== !
  ! Calculation in domain interior
  ! ========================================================================== !
  else
    ! Not at edge of LAM so compute fluxes

    ! Loop through horizontal x DoFs
    do df_idx = 1, face_selector_ew(map_w3_2d(1))
      df = local_dofs_x(df_idx)
      dof_offset = (df+1)/2 - 1 ! 0 for W, 1 for E

      ! Pull out index to avoid multiple indirections
      idx_w2h = map_w2h(df)

      ! Loop through layers
      do k = 0, nlayers - 1

        flux_face_i = dry_flux(idx_w2h+k)

        sign_flux = INT(SIGN(1.0_r_tran, flux_face_i))

        ! Set an offset for the stencil index, based on dep point sign
        sign_offset = (1 - sign_flux) / 2   ! 0 if sign == 1, 1 if sign == -1

        ! -------------------------------------------------------------------- !
        ! Find departure point
        ! -------------------------------------------------------------------- !
        num_int_cells = 0
        running_flux_i = 0.0_r_tran

        ! Step backwards through flux to find the dimensionless departure point:
        ! First, the integer part of departure point (num_int_cells)
        do j = 1, stencil_extent

          ! If this column has idx 0, find relative index along column of
          ! the departure cell, between - stencil_half and stencil_half
          rel_idx = dof_offset - (j - 1) * sign_flux + sign_offset - 1

          ! Determine the index in the stencil from rel_idx
          ! e.g. for extent 4:
          ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
          ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
          stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2
          col_idx = stencil_map_x(1,stencil_idx)

          ! Get mass for the next cell
          flux_cell_j = dry_mass_for_x(col_idx + k)

          ! Check if the mass swept up to this cell now exceeds the flux:
          ! if so, the we have found integer number of cells and the departure
          ! point will lie in cell (num_int_cells + 1), so exit do-loop
          if (running_flux_i + flux_cell_j > abs(flux_face_i)) EXIT

          ! Increment running values, if we have not exited loop
          running_flux_i = running_flux_i + flux_cell_j
          num_int_cells = num_int_cells + 1
        end do

        ! Second, the fractional part of departure point
        if (num_int_cells < stencil_extent) then
          ! Departure point is within stencil: determine fractional departure
          ! distance and fractional flux

          ! Set fractional flux. running_flux_i is now the integer flux
          frac_flux_face_i = sign_flux*(abs(flux_face_i) - running_flux_i)

          ! If this column has idx 0, find relative index along column of
          ! the departure cell, between - stencil_half and stencil_half
          rel_idx = dof_offset - num_int_cells * sign_flux + sign_offset - 1

          ! Determine the index in the stencil from rel_idx
          ! e.g. for extent 4:
          ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
          ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
          stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2
          col_idx = stencil_map_x(1,stencil_idx)

          frac_dep_dist = frac_flux_face_i / dry_mass_for_x(col_idx + k)

        else if (cap_dep_points) then
          ! Departure point has exceeded stencil, but the user has specified
          ! to cap the departure point so set fractional parts to be zero
          frac_flux_face_i = 0.0_r_tran
          frac_dep_dist = 0.0_r_tran

        else
          ! Departure point has exceeded stencil depth so throw an error
          write(log_scratch_space, '(A,I8)') 'hori_x_dep_dist: departure ' // &
            'points have exceeded stencil depth of ', stencil_extent
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        ! Set the values of the output fields
        frac_dry_flux(idx_w2h+k) = frac_flux_face_i
        dep_dist(idx_w2h+k) = real(sign_flux*num_int_cells, r_tran) + frac_dep_dist
      end do
    end do

    ! ----------------------------------------------------------------------- !

    ! Loop through horizontal y DoFs
    do df_idx = 1, face_selector_ns(map_w3_2d(1))
      df = local_dofs_y(df_idx)
      dof_offset = (df+1)/2 - 1 ! 0 for S, 1 for N

      ! Pull out index to avoid multiple indirections
      idx_w2h = map_w2h(df)

      ! Loop through layers
      do k = 0, nlayers - 1

        flux_face_i = dry_flux(idx_w2h+k)

        ! NB: the y-direction for stencils runs in the opposite direction to the
        ! direction used for the winds. Compensate for this with a minus sign
        sign_flux = -INT(SIGN(1.0_r_tran, flux_face_i))

        ! Set an offset for the stencil index, based on dep point sign
        sign_offset = (1 - sign_flux) / 2   ! 0 if sign == 1, 1 if sign == -1

        ! -------------------------------------------------------------------- !
        ! Find departure point
        ! -------------------------------------------------------------------- !
        num_int_cells = 0
        running_flux_i = 0.0_r_tran

        ! Step backwards through flux to find the dimensionless departure point:
        ! First, the integer part of departure point (num_int_cells)
        do j = 1, stencil_extent

          ! If this column has idx 0, find relative index along column of
          ! the departure cell, between - stencil_half and stencil_half
          rel_idx = dof_offset - (j - 1) * sign_flux + sign_offset - 1

          ! Determine the index in the stencil from rel_idx
          ! e.g. for extent 4:
          ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
          ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
          stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2
          col_idx = stencil_map_y(1,stencil_idx)

          ! Get mass for the next cell
          flux_cell_j = dry_mass_for_y(col_idx + k)

          ! Check if the mass swept up to this cell now exceeds the flux:
          ! if so, the we have found integer number of cells and the departure
          ! point will lie in cell (num_int_cells + 1), so exit do-loop
          if (running_flux_i + flux_cell_j > abs(flux_face_i)) EXIT

          ! Increment running values, if we have not exited loop
          running_flux_i = running_flux_i + flux_cell_j
          num_int_cells = num_int_cells + 1
        end do

        ! Second, the fractional part of departure point
        if (num_int_cells < stencil_extent) then
          ! Departure point is within stencil: determine fractional departure
          ! distance and fractional flux

          ! Set fractional flux. running_flux_i is now the integer flux
          frac_flux_face_i = sign_flux*(abs(flux_face_i) - running_flux_i)

          ! If this column has idx 0, find relative index along column of
          ! the departure cell, between - stencil_half and stencil_half
          rel_idx = dof_offset - num_int_cells * sign_flux + sign_offset - 1

          ! Determine the index in the stencil from rel_idx
          ! e.g. for extent 4:
          ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
          ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
          stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2
          col_idx = stencil_map_y(1,stencil_idx)

          frac_dep_dist = frac_flux_face_i / dry_mass_for_y(col_idx + k)

        else if (cap_dep_points) then
          ! Departure point has exceeded stencil, but the user has specified
          ! to cap the departure point so set fractional parts to be zero
          frac_flux_face_i = 0.0_r_tran
          frac_dep_dist = 0.0_r_tran

        else
          ! Departure point has exceeded stencil depth so throw an error
          write(log_scratch_space, '(A,I8)') 'hori_y_dep_dist: departure ' // &
            'points have exceeded stencil depth of ', stencil_extent
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        ! Set the values of the output fields
        ! NB: minus signs to return to conventional y-direction
        frac_dry_flux(idx_w2h+k) = -frac_flux_face_i
        dep_dist(idx_w2h+k) = - (real(sign_flux*num_int_cells, r_tran) + frac_dep_dist)

      end do
    end do

  end if

end subroutine hori_dep_dist_ffsl_code

end module hori_dep_dist_ffsl_kernel_mod
