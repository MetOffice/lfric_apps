!-------------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to compute the horizontal departure distances for FFSL on the
!!        cubed sphere
!> @details This code calculates the distance which is swept through a cell
!!          in one dimension during one time step of consistent transport. The
!!          arrival point is the cell face and the departure point is calculated
!!          and stored as a field.
!!          The departure points are a dimensionless displacement, corresponding
!!          to the number of cells moved by a fluid parcel. The fractional flux
!!          is also computed.
!!          This differs from the other kernel as it takes into account the
!!          swapping of directions at the boundaries of cubed sphere panels.

module hori_dep_dist_ffsl_sphere_kernel_mod

use argument_mod,                only : arg_type,                  &
                                        GH_FIELD, GH_REAL,         &
                                        GH_WRITE, GH_READ,         &
                                        GH_SCALAR, GH_INTEGER,     &
                                        STENCIL, CROSS2D,          &
                                        GH_LOGICAL, CELL_COLUMN,   &
                                        ANY_DISCONTINUOUS_SPACE_2, &
                                        ANY_DISCONTINUOUS_SPACE_3, &
                                        ANY_DISCONTINUOUS_SPACE_1, &
                                        ANY_DISCONTINUOUS_SPACE_4, &
                                        GH_READWRITE
use fs_continuity_mod,           only : W3, W2h
use constants_mod,               only : r_tran, i_def, l_def, r_def
use kernel_mod,                  only : kernel_type
use reference_element_mod,       only : E, W, N, S

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: hori_dep_dist_ffsl_sphere_kernel_type
  private
  type(arg_type) :: meta_args(14) = (/                                         &
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE,   ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
      arg_type(GH_FIELD,   GH_REAL,    GH_WRITE,   ANY_DISCONTINUOUS_SPACE_2), & ! frac_ref_flux
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,    ANY_DISCONTINUOUS_SPACE_2), & ! ref_flux
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,    W3, STENCIL(CROSS2D)),      & ! ref_mass_for_x
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,    W3, STENCIL(CROSS2D)),      & ! ref_mass_for_y
      arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_4), & ! dep_lowest_k
      arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE,   ANY_DISCONTINUOUS_SPACE_4), & ! dep_highest_k
      arg_type(GH_FIELD,   GH_REAL,    GH_READ,    ANY_DISCONTINUOUS_SPACE_1), & ! panel_id
      arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,    ANY_DISCONTINUOUS_SPACE_1), & ! panel_edge_dist
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,    ANY_DISCONTINUOUS_SPACE_3), & ! face_selector ew
      arg_type(GH_FIELD,   GH_INTEGER, GH_READ,    ANY_DISCONTINUOUS_SPACE_3), & ! face_selector ns
      arg_type(GH_FIELD,   GH_INTEGER, GH_READWRITE,                           &
                                                   ANY_DISCONTINUOUS_SPACE_3), & ! error_flag
      arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                               & ! ndep
      arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)                                & ! cap_dep_points
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: hori_dep_dist_ffsl_sphere_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: hori_dep_dist_ffsl_sphere_code

contains

!> @brief Kernel to compute the horizontal departure distances for FFSL on the
!!        cubed sphere, where calculations may cross a panel edge at which the
!!        x and y directions are rotated.
!> @param[in]     nlayers            The number of layers in the mesh
!> @param[in,out] dep_dist           Field with dep distances
!> @param[in,out] frac_ref_flux      Fractional part of the reference flux
!> @param[in]     ref_flux           The mass flux used in a one-dimensional
!!                                   transport step for the reference density
!> @param[in]     ref_mass_for_x     The reference mass field to be used for
!!                                   the x-direction transport step
!> @param[in]     stencil_sizes_x    Sizes of the branches of the cross stencil
!> @param[in]     stencil_max_x      Maximum size of a cross stencil branch
!> @param[in]     stencil_map_x      Map of DoFs in the stencil for mass field
!> @param[in]     ref_mass_for_y     The reference mass field to be used for
!!                                   the y-direction transport step
!> @param[in]     stencil_sizes_y    Sizes of the branches of the cross stencil
!> @param[in]     stencil_max_y      Maximum size of a cross stencil branch
!> @param[in]     stencil_map_y      Map of DoFs in the stencil for mass field
!> @param[in,out] dep_lowest_k       2D integer multidata W2H field, storing
!!                                   the lowest model level to use in each
!!                                   integer flux sum, for each column
!> @param[in,out] dep_highest_k      2D integer multidata W2H field, storing
!!                                   the highest model level to use in each
!!                                   integer flux sum, for each column
!> @param[in]     panel_id           Field containing IDs of mesh panels
!> @param[in]     panel_edge_dist_W  2D field containing the distance of each
!!                                   column from the panel edge to the West
!> @param[in]     panel_edge_dist_E  2D field containing the distance of each
!!                                   column from the panel edge to the East
!> @param[in]     panel_edge_dist_S  2D field containing the distance of each
!!                                   column from the panel edge to the South
!> @param[in]     panel_edge_dist_N  2D field containing the distance of each
!!                                   column from the panel edge to the North
!> @param[in]     face_selector_ew   2D field indicating faces to loop over in x
!> @param[in]     face_selector_ns   2D field indicating faces to loop over in y
!> @param[in,out] error_flag         2D field for storing error flag, only used
!!                                   optionally if departure points not capped
!> @param[in]     cap_dep_points     Flag for whether departure points should be
!!                                   capped if they exceed the stencil depth
!> @param[in]     ndep               Number of multidata points for departure
!!                                   index fields
!> @param[in]     ndf_w2h            Number of DoFs per cell for W2H
!> @param[in]     undf_w2h           Num of W2H DoFs for this partition
!> @param[in]     map_w2h            Map of lowest-cell W2H DoFs
!> @param[in]     ndf_w3             Number of DoFs per cell for W3
!> @param[in]     undf_w3            Num of W3 DoFs for this partition
!> @param[in]     map_w3             Map of lowest-cell W3 DoFs
!> @param[in]     ndf_depk           Num of DoFs for dep idx fields per cell
!> @param[in]     undf_depk          Num of DoFs for this partition for dep
!!                                   idx fields
!> @param[in]     map_depk           Map for departure index fields
!> @param[in]     ndf_pid            Num of DoFs for panel ID per cell
!> @param[in]     undf_pid           Num of DoFs for this partition for panel ID
!> @param[in]     map_pid            Map for panel ID
!> @param[in]     ndf_w3_2d          Number of DoFs for 2D W3 per cell
!> @param[in]     undf_w3_2d         Number of DoFs for this partition for 2D W3
!> @param[in]     map_w3_2d          Map for 2D W3
subroutine hori_dep_dist_ffsl_sphere_code( nlayers,             &
                                           dep_dist,            &
                                           frac_ref_flux,       &
                                           ref_flux,            &
                                           ref_mass_for_x,      &
                                           stencil_sizes_x,     &
                                           stencil_max_x,       &
                                           stencil_map_x,       &
                                           ref_mass_for_y,      &
                                           stencil_sizes_y,     &
                                           stencil_max_y,       &
                                           stencil_map_y,       &
                                           dep_lowest_k,        &
                                           dep_highest_k,       &
                                           panel_id,            &
                                           panel_edge_dist_W,   &
                                           panel_edge_dist_E,   &
                                           panel_edge_dist_S,   &
                                           panel_edge_dist_N,   &
                                           face_selector_ew,    &
                                           face_selector_ns,    &
                                           error_flag,          &
                                           ndep,                &
                                           cap_dep_points,      &
                                           ndf_w2h,             &
                                           undf_w2h,            &
                                           map_w2h,             &
                                           ndf_w3,              &
                                           undf_w3,             &
                                           map_w3,              &
                                           ndf_depk,            &
                                           undf_depk,           &
                                           map_depk,            &
                                           ndf_pid,             &
                                           undf_pid,            &
                                           map_pid,             &
                                           ndf_w3_2d,           &
                                           undf_w3_2d,          &
                                           map_w3_2d )

  use hori_dep_dist_ffsl_kernel_mod,  only: hori_dep_dist_ffsl_1d
  use panel_edge_support_mod,         only: crosses_rotated_panel_edge

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers, ndep
  logical(kind=l_def), intent(in)    :: cap_dep_points
  integer(kind=i_def), intent(in)    :: ndf_w2h, ndf_w3, ndf_w3_2d, ndf_depk
  integer(kind=i_def), intent(in)    :: undf_w2h, undf_w3, undf_w3_2d, undf_depk
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: stencil_sizes_x(4)
  integer(kind=i_def), intent(in)    :: stencil_sizes_y(4)
  integer(kind=i_def), intent(in)    :: stencil_max_x, stencil_max_y
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_depk(ndf_depk)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: map_w3_2d(ndf_w3_2d)
  integer(kind=i_def), intent(in)    :: stencil_map_x(ndf_w3, stencil_max_x, 4)
  integer(kind=i_def), intent(in)    :: stencil_map_y(ndf_w3, stencil_max_y, 4)
  integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
  integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)
  integer(kind=i_def), intent(inout) :: error_flag(undf_w3_2d)
  integer(kind=i_def), intent(inout) :: dep_lowest_k(undf_depk)
  integer(kind=i_def), intent(inout) :: dep_highest_k(undf_depk)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  real(kind=r_tran),   intent(in)    :: ref_mass_for_x(undf_w3)
  real(kind=r_tran),   intent(in)    :: ref_mass_for_y(undf_w3)
  real(kind=r_tran),   intent(in)    :: ref_flux(undf_w2h)
  real(kind=r_tran),   intent(inout) :: dep_dist(undf_w2h)
  real(kind=r_tran),   intent(inout) :: frac_ref_flux(undf_w2h)

  ! Internal arguments
  integer(kind=i_def) :: i, ipanel
  integer(kind=i_def) :: edge_dist_E, edge_dist_W, edge_dist_S, edge_dist_N
  logical(kind=l_def) :: use_spherical_treatment
  integer(kind=i_def) :: stencil_extent_xl, stencil_extent_xr
  integer(kind=i_def) :: stencil_extent_yl, stencil_extent_yr
  integer(kind=i_def) :: stencil_map_x_1d(-stencil_max_x:stencil_max_x)
  integer(kind=i_def) :: stencil_map_y_1d(-stencil_max_y:stencil_max_y)

  ipanel = INT(panel_id(map_pid(1)), i_def)
  edge_dist_E = panel_edge_dist_E(map_pid(1))
  edge_dist_W = panel_edge_dist_W(map_pid(1))
  edge_dist_S = panel_edge_dist_S(map_pid(1))
  edge_dist_N = panel_edge_dist_N(map_pid(1))

  ! Form X and Y 1D stencils
  stencil_extent_xl = stencil_sizes_x(1) - 1
  stencil_extent_xr = stencil_sizes_x(3) - 1
  stencil_extent_yl = stencil_sizes_y(2) - 1
  stencil_extent_yr = stencil_sizes_y(4) - 1

  do i = -stencil_extent_xl, 0
    stencil_map_x_1d(i) = stencil_map_x(1, 1-i, 1)
  end do
  do i = 1, stencil_extent_xr
    stencil_map_x_1d(i) = stencil_map_x(1, i+1, 3)
  end do

  do i = -stencil_extent_yl, 0
    stencil_map_y_1d(i) = stencil_map_y(1, 1-i, 2)
  end do
  do i = 1, stencil_extent_yr
    stencil_map_y_1d(i) = stencil_map_y(1, i+1, 4)
  end do

  ! X direction ================================================================
  use_spherical_treatment = crosses_rotated_panel_edge(                        &
      edge_dist_W, edge_dist_E, MAX(stencil_extent_xl, stencil_extent_xr),     &
      ipanel, 1                                                                &
  )

  if (use_spherical_treatment) then
    call hori_dep_dist_ffsl_sphere_1d( nlayers,             &
                                       .true.,              &
                                       dep_dist,            &
                                       frac_ref_flux,       &
                                       ref_flux,            &
                                       ref_mass_for_x,      &
                                       ref_mass_for_y,      &
                                       stencil_extent_xl,   &
                                       stencil_extent_xr,   &
                                       stencil_max_x,       &
                                       stencil_map_x_1d,    &
                                       dep_lowest_k,        &
                                       dep_highest_k,       &
                                       ipanel,              &
                                       edge_dist_W,         &
                                       edge_dist_E,         &
                                       face_selector_ew,    &
                                       error_flag,          &
                                       ndep,                &
                                       cap_dep_points,      &
                                       ndf_w2h,             &
                                       undf_w2h,            &
                                       map_w2h,             &
                                       ndf_w3,              &
                                       undf_w3,             &
                                       map_w3,              &
                                       ndf_depk,            &
                                       undf_depk,           &
                                       map_depk,            &
                                       ndf_w3_2d,           &
                                       undf_w3_2d,          &
                                       map_w3_2d )
  else
    call hori_dep_dist_ffsl_1d( nlayers,             &
                                .true.,              &
                                dep_dist,            &
                                frac_ref_flux,       &
                                ref_flux,            &
                                ref_mass_for_x,      &
                                stencil_extent_xl,   &
                                stencil_extent_xr,   &
                                stencil_max_x,       &
                                stencil_map_x_1d,    &
                                dep_lowest_k,        &
                                dep_highest_k,       &
                                face_selector_ew,    &
                                error_flag,          &
                                ndep,                &
                                cap_dep_points,      &
                                ndf_w2h,             &
                                undf_w2h,            &
                                map_w2h,             &
                                ndf_w3,              &
                                undf_w3,             &
                                map_w3,              &
                                ndf_depk,            &
                                undf_depk,           &
                                map_depk,            &
                                ndf_w3_2d,           &
                                undf_w3_2d,          &
                                map_w3_2d )
  end if

  ! Y direction ================================================================
  use_spherical_treatment = crosses_rotated_panel_edge(                        &
      edge_dist_S, edge_dist_N, MAX(stencil_extent_yl, stencil_extent_yr),     &
      ipanel, 2                                                                &
  )

  if (use_spherical_treatment) then
    call hori_dep_dist_ffsl_sphere_1d( nlayers,             &
                                       .false.,             &
                                       dep_dist,            &
                                       frac_ref_flux,       &
                                       ref_flux,            &
                                       ref_mass_for_y,      &
                                       ref_mass_for_x,      &
                                       stencil_extent_yl,   &
                                       stencil_extent_yr,   &
                                       stencil_max_y,       &
                                       stencil_map_y_1d,    &
                                       dep_lowest_k,        &
                                       dep_highest_k,       &
                                       ipanel,              &
                                       edge_dist_S,         &
                                       edge_dist_N,         &
                                       face_selector_ns,    &
                                       error_flag,          &
                                       ndep,                &
                                       cap_dep_points,      &
                                       ndf_w2h,             &
                                       undf_w2h,            &
                                       map_w2h,             &
                                       ndf_w3,              &
                                       undf_w3,             &
                                       map_w3,              &
                                       ndf_depk,            &
                                       undf_depk,           &
                                       map_depk,            &
                                       ndf_w3_2d,           &
                                       undf_w3_2d,          &
                                       map_w3_2d )
  else
    call hori_dep_dist_ffsl_1d( nlayers,             &
                                .false.,             &
                                dep_dist,            &
                                frac_ref_flux,       &
                                ref_flux,            &
                                ref_mass_for_y,      &
                                stencil_extent_yl,   &
                                stencil_extent_yr,   &
                                stencil_max_y,       &
                                stencil_map_y_1d,    &
                                dep_lowest_k,        &
                                dep_highest_k,       &
                                face_selector_ns,    &
                                error_flag,          &
                                ndep,                &
                                cap_dep_points,      &
                                ndf_w2h,             &
                                undf_w2h,            &
                                map_w2h,             &
                                ndf_w3,              &
                                undf_w3,             &
                                map_w3,              &
                                ndf_depk,            &
                                undf_depk,           &
                                map_depk,            &
                                ndf_w3_2d,           &
                                undf_w3_2d,          &
                                map_w3_2d )
  end if

end subroutine hori_dep_dist_ffsl_sphere_code

! ============================================================================ !
! SINGLE UNDERLYING 1D ROUTINE
! ============================================================================ !

!> @brief Computes the horizontal departure distances for FFSL in one direction
subroutine hori_dep_dist_ffsl_sphere_1d( nlayers,             &
                                         x_direction,         &
                                         dep_dist,            &
                                         frac_ref_flux,       &
                                         ref_flux,            &
                                         ref_mass,            &
                                         ref_mass_swapped,    &
                                         stencil_extent_l,    &
                                         stencil_extent_r,    &
                                         stencil_array_max,   &
                                         stencil_map,         &
                                         dep_lowest_k,        &
                                         dep_highest_k,       &
                                         ipanel,              &
                                         panel_edge_dist_l,   &
                                         panel_edge_dist_r,   &
                                         face_selector,       &
                                         error_flag,          &
                                         ndep,                &
                                         cap_dep_points,      &
                                         ndf_w2h,             &
                                         undf_w2h,            &
                                         map_w2h,             &
                                         ndf_w3,              &
                                         undf_w3,             &
                                         map_w3,              &
                                         ndf_depk,            &
                                         undf_depk,           &
                                         map_depk,            &
                                         ndf_w3_2d,           &
                                         undf_w3_2d,          &
                                         map_w3_2d )

  use panel_edge_support_mod, only: rotated_panel_neighbour, FAR_AWAY

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers, ndep
  logical(kind=l_def), intent(in)    :: cap_dep_points
  integer(kind=i_def), intent(in)    :: ndf_w2h, ndf_w3, ndf_w3_2d, ndf_depk
  integer(kind=i_def), intent(in)    :: undf_w2h, undf_w3, undf_w3_2d, undf_depk
  integer(kind=i_def), intent(in)    :: stencil_extent_l, stencil_extent_r
  integer(kind=i_def), intent(in)    :: stencil_array_max
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_depk(ndf_depk)
  integer(kind=i_def), intent(in)    :: map_w3_2d(ndf_w3_2d)
  integer(kind=i_def), intent(in)    :: stencil_map(-stencil_array_max:stencil_array_max)
  integer(kind=i_def), intent(in)    :: face_selector(undf_w3_2d)
  integer(kind=i_def), intent(inout) :: error_flag(undf_w3_2d)
  integer(kind=i_def), intent(inout) :: dep_lowest_k(undf_depk)
  integer(kind=i_def), intent(inout) :: dep_highest_k(undf_depk)
  integer(kind=i_def), intent(in)    :: ipanel
  integer(kind=i_def), intent(in)    :: panel_edge_dist_l
  integer(kind=i_def), intent(in)    :: panel_edge_dist_r
  real(kind=r_tran),   intent(in)    :: ref_flux(undf_w2h)
  real(kind=r_tran),   intent(inout) :: dep_dist(undf_w2h)
  real(kind=r_tran),   intent(inout) :: frac_ref_flux(undf_w2h)
  logical(kind=l_def), intent(in)    :: x_direction

  real(kind=r_tran), target, intent(in) :: ref_mass(undf_w3)
  real(kind=r_tran), target, intent(in) :: ref_mass_swapped(undf_w3)

  ! Local arrays
  integer(kind=i_def) :: local_dofs(2)
  integer(kind=i_def) :: sign_flux(nlayers)
  integer(kind=i_def) :: sign_switch(nlayers)
  integer(kind=i_def) :: running_num_cells(nlayers)
  integer(kind=i_def) :: rel_idx(nlayers)
  integer(kind=i_def) :: sign_offset(nlayers)
  integer(kind=i_def) :: stencil_max(nlayers)
  real(kind=r_tran)   :: frac_flux_tmp(nlayers)
  real(kind=r_tran)   :: running_int_flux(nlayers)
  real(kind=r_tran)   :: swept_mass(nlayers)
  real(kind=r_tran)   :: frac_dep_dist(nlayers)
  integer(kind=i_def) :: swap_switch(nlayers)
  integer(kind=i_def) :: rounding_factor(nlayers)
  integer(kind=i_def) :: int_cell_switch(nlayers)
  integer(kind=i_def) :: k_low, k_high, k_low_prev, k_high_prev
  integer(kind=i_def) :: stencil_sign, sign_offset_k

  real(kind=r_tran), pointer :: ref_mass_ptr(:)

  ! Local scalars
  integer(kind=i_def) :: j, k
  integer(kind=i_def) :: df, df_idx, w2h_idx, depk_idx, depk_idx_l
  integer(kind=i_def) :: dof_offset
  integer(kind=i_def) :: col_idx, rel_idx_k
  integer(kind=i_def) :: stencil_max_k
  integer(kind=i_def) :: rot_edge_dist_l, rot_edge_dist_r
  integer(kind=i_def) :: rotated_panel_l, rotated_panel_r
  integer(kind=i_def) :: ndep_half
  real(kind=r_tran)   :: direction

  if (x_direction) then
    local_dofs = (/ W, E /)
    direction = 1.0_r_tran
  else
    ! y-direction
    local_dofs = (/ S, N /)
    direction = -1.0_r_tran
  end if

  ! Set stencil info -----------------------------------------------------------
  ndep_half = (ndep - 1_i_def) / 2_i_def

  ! Pre-determine aspects of crossing panel boundary ---------------------------
  rotated_panel_l = rotated_panel_neighbour(ipanel, local_dofs(1))
  rotated_panel_r = rotated_panel_neighbour(ipanel, local_dofs(2))

  if (ABS(rotated_panel_l) > 0) then
    rot_edge_dist_l = -ABS(panel_edge_dist_l)
  else
    rot_edge_dist_l = -FAR_AWAY
  end if
  if (ABS(rotated_panel_r) > 0) then
    rot_edge_dist_r = ABS(panel_edge_dist_r)
  else
    rot_edge_dist_r = FAR_AWAY
  end if

  ! Loop through horizontal DoFs -----------------------------------------------
  do df_idx = 1, face_selector(map_w3_2d(1))
    df = local_dofs(df_idx)

    ! Set a local offset, dependent on the face we are looping over
    select case (df)
    case (W, S)
      dof_offset = 0
    case (E, N)
      dof_offset = 1
    end select

    ! Pull out index to avoid multiple indirections
    w2h_idx = map_w2h(df)
    depk_idx = map_depk(df) + ndep_half

    ! NB: the y-direction for stencils runs in the opposite direction to the
    ! direction used for the winds. Compensate for this with a minus sign
    sign_flux(:) = INT(SIGN(                                                   &
        1.0_r_tran, direction * ref_flux(w2h_idx : w2h_idx + nlayers - 1)      &
    ), i_def)

    ! Set an offset for the stencil index, based on dep point sign
    sign_offset = (1 - sign_flux) / 2   ! 0 if sign == 1, 1 if sign == -1

    running_num_cells(:) = 0
    running_int_flux(:) = 0.0_r_tran

    ! Loop over both potential signs of the wind (so both sides of stencil)
    do stencil_sign = -1, 1, 2

      ! Set an offset for the stencil index, based on dep point sign
      sign_offset_k = (1 - stencil_sign) / 2   ! 0 if sign == 1, 1 if sign == -1

      ! Determine extent to loop to for each direction
      stencil_max_k = (                                                        &
        sign_offset_k*(stencil_extent_r + 1 - dof_offset)                      &
        + (1 - sign_offset_k)*(stencil_extent_l + dof_offset)                  &
      )

      swept_mass(:) = 0.0_r_tran

      ! Initialise indices for lowest and highest levels
      k_low = 0
      k_high = nlayers - 1

      ! Step backwards through flux to find the dimensionless departure point:
      do j = 1, stencil_max_k

        ! Get index for this column ------------------------------------------
        ! If this column has idx 0, find relative index for the column of the
        ! departure cell, between -stencil_extent_l and stencil_extent_r, e.g.
        ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
        rel_idx_k = MIN(stencil_extent_r, MAX(-stencil_extent_l,               &
            dof_offset - (j - 1) * stencil_sign + sign_offset_k - 1            &
        ))

        ! Determine if we have crossed into a rotated panel
        if (rel_idx_k < rot_edge_dist_l + 1 .or. rel_idx_k > rot_edge_dist_r - 1) then
          ref_mass_ptr => ref_mass_swapped
        else
          ref_mass_ptr => ref_mass
        end if

        col_idx = stencil_map(rel_idx_k)

        ! Calculate mask for whether to include each cell in the calculation
        ! This is 1 if:
        ! (a) the wind sign is correct, and
        ! (b) the departure cell has not yet been found before considering
        !     this column, and
        ! (c) when including the volume/mass of the cell for this column in
        !     the running flux total, the total is still less than the flux
        sign_switch(k_low+1 : k_high+1) = (                                    &
          ABS((sign_flux(k_low+1 : k_high+1) + stencil_sign) / 2)              &
        )
        int_cell_switch(k_low+1 : k_high+1) = (                                &
          ! 1 if the wind sign is correct, 0 otherwise
          sign_switch(k_low+1 : k_high+1)                                      &
          ! 1 if dep cell has not been found prior to this column, 0 otherwise
          * ABS(                                                               &
            (SIGN(1, running_num_cells(k_low+1 : k_high+1) + 2 - j) + 1) / 2   &
          )                                                                    &
          ! 1 if this cell is not yet the dep cell, 0 otherwise
          * ABS(                                                               &
            (INT(SIGN(1.0_r_tran,                                              &
              ABS(ref_flux(w2h_idx + k_low : w2h_idx + k_high))                &
              - swept_mass(k_low+1 : k_high+1)                                 &
              - ref_mass_ptr(col_idx + k_low : col_idx + k_high)               &
          ), i_def) + 1) / 2)                                                  &
        )
        ! We need to keep track of the total swept mass, to determine if the
        ! total mass has been swept through or not
        swept_mass(k_low+1 : k_high+1) = (                                     &
          swept_mass(k_low+1 : k_high+1)                                       &
          + sign_switch(k_low+1 : k_high+1)                                    &
          * ref_mass_ptr(col_idx + k_low : col_idx + k_high)                   &
        )
        ! Increment running counts of integer flux and num cells
        running_int_flux(k_low+1 : k_high+1) = (                               &
          running_int_flux(k_low+1 : k_high+1)                                 &
          + int_cell_switch(k_low+1 : k_high+1)                                &
          * ref_mass_ptr(col_idx + k_low : col_idx + k_high)                   &
        )
        running_num_cells(k_low+1 : k_high+1) = (                              &
          running_num_cells(k_low+1 : k_high+1)                                &
          + int_cell_switch(k_low+1 : k_high+1)                                &
        )

        ! Set new lowest and highest levels ------------------------------------
        ! Convert rel_idx_k into multidata index for storing these values
        depk_idx_l = depk_idx + rel_idx_k

        ! These indices are the bounds of the array to include in the next
        ! calculation (these will also be saved for use in flux calculations)
        if (MAXVAL(int_cell_switch(k_low+1 : k_high+1)) > 0_i_def) then
          k_low_prev = k_low
          k_high_prev = k_high

          ! Find lowest index
          do k = k_low_prev, k_high_prev
            if (int_cell_switch(k+1) == 1) then
              k_low = k
              EXIT
            end if
          end do

          ! Find highest index
          do k = k_high_prev, k_low_prev, -1
            if (int_cell_switch(k+1) == 1) then
              k_high = k
              EXIT
            end if
          end do

          dep_lowest_k(depk_idx_l) = k_low
          dep_highest_k(depk_idx_l) = k_high

        else
          ! All of the departure cells have now been found for this sweep
          k_low = 0
          k_high = -1
          dep_lowest_k(depk_idx_l) = k_low
          dep_highest_k(depk_idx_l) = k_high
          ! We don't need to continue find dep distances in this direction,
          ! so we can exit this sweep's loop
          EXIT
        end if
      end do  ! Loop through columns (on one side of stencil)
    end do  ! Loop through possible wind directions

    ! ======================================================================== !
    ! FRACTIONAL PART OF DEPARTURE DISTANCE
    ! ======================================================================== !

    frac_flux_tmp(:) = sign_flux(:)*(                                          &
      ABS(ref_flux(w2h_idx : w2h_idx+nlayers-1)) - running_int_flux(:)         &
    )

    ! If this column has idx 0, find relative index for the column of the
    ! departure cell, between -stencil_extent_l and stencil_extent_r, e.g.
    ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
    rel_idx(:) = MIN(stencil_extent_r, MAX(-stencil_extent_l,                  &
      dof_offset - running_num_cells(:) * sign_flux(:) + sign_offset - 1       &
    ))

    stencil_max = (                                                            &
      sign_offset*(stencil_extent_r + 1 - dof_offset)                          &
      + (1 - sign_offset)*(stencil_extent_l + dof_offset)                      &
    )

    ! Determine whether this column has crossed a panel edge
    ! This factor is 1 if the column has crossed a panel edge, and 0 if it
    ! has not. It is made from (1 - switch_l*switch_r) where:
    ! switch_l is 0 when on a rotated left panel, 1 otherwise
    ! switch_r is 0 when on a rotated right panel, 1 otherwise
    swap_switch(:) = (                                                         &
        1 - (1 + SIGN(1, rel_idx - rot_edge_dist_l - 1))                       &
        ! 0 if crossed to rotated panel, 1 if not
        * (1 + SIGN(1, rot_edge_dist_r - rel_idx - 1)) / 4                     &
    )

    do k = 1, nlayers
      col_idx = stencil_map(rel_idx(k))

      if (running_num_cells(k) < stencil_max(k) .and. swap_switch(k) == 0_i_def) then
        ! Safe to set the fractional dep dist. Not crossed to rotated panel.
        frac_dep_dist(k) = frac_flux_tmp(k) / ref_mass(col_idx + k - 1)

        ! NB: here we could use the swap_switch array to combine the cases of
        ! having crossed / not crossed to a rotated panel. However, this changes
        ! answers at the bit level compared with the simple non-spherical
        ! version of this kernel. Whether we perform this kernel or the simple
        ! version depends on whether there is a panel edge somewhere in the
        ! stencil -- thus changing the size of the stencil could change the
        ! answer at the bit level. Therefore we do not combine these cases to
        ! ensure that answers are not changed by changing the stencil size

      else if (running_num_cells(k) < stencil_max(k)) then
        ! Safe to set the fractional dep dist. Crossed to rotated panel.
        frac_dep_dist(k) = frac_flux_tmp(k) / ref_mass_swapped(col_idx + k - 1)

      else if (cap_dep_points) then
        ! Departure point has exceeded stencil, but the user has specified
        ! to cap the departure point. Do this by setting the fractional
        ! part of the departure distance to be 0
        frac_dep_dist(k) = 0.0_r_tran
        ! Do not set the fractional flux to zero, as this will still be used

      else
        ! Departure point has exceeded stencil depth so indicate that an
        ! error needs to be thrown, rather than either seg-faulting or
        ! having the model fail ungracefully at a later point
        frac_dep_dist(k) = 0.0_r_tran
        error_flag(map_w3_2d(1)) = error_flag(map_w3_2d(1)) + 1_i_def
      end if
    end do

    ! NB: Y calculation has a minus sign to return to conventional y-direction
    dep_dist(w2h_idx : w2h_idx+nlayers-1) = (                                  &
      direction                                                                &
      * (REAL(sign_flux(:)*running_num_cells(:), r_tran) + frac_dep_dist(:))   &
    )

    ! ======================================================================== !
    ! FRACTIONAL FLUX
    ! ======================================================================== !

    ! If floating point rounding puts the dep_dist into the next integer value
    ! then set the fractional part back to zero with rounding factor
    rounding_factor(:) = (                                                     &
      1_i_def + ABS(running_num_cells(:))                                      &
      - ABS(INT(dep_dist(w2h_idx : w2h_idx+nlayers-1), i_def))                 &
    )

    ! If the rounding factor is 0, correct the dep_low/highest_k fields, as
    ! these may no longer be right as the dep point crosses into a new cell
    if (MINVAL(rounding_factor) == 0_i_def) then
      do k = 1, nlayers
        if (rounding_factor(k) == 0_i_def) then
          ! Use the Courant number to get the index of dep_low/highest field
          rel_idx_k = MIN(ndep_half, MAX(-ndep_half,                           &
            dof_offset - sign_flux(k)                                          &
            * (running_num_cells(k) + 1_i_def - sign_offset(k))                &
          ))
          depk_idx_l = depk_idx + rel_idx_k
          ! If this level is now the low/highest level, update the fields
          if (k - 1 < dep_lowest_k(depk_idx_l)) then
            dep_lowest_k(depk_idx_l) = k - 1
          end if
          if (k - 1 > dep_highest_k(depk_idx_l)) then
            dep_highest_k(depk_idx_l) = k - 1
          end if
        end if
      end do
    end if

    ! Fractional flux is the difference between the total and integer fluxes
    ! NB: Y calculation has a minus sign to return to conventional y-direction
    frac_ref_flux(w2h_idx : w2h_idx+nlayers-1) = (                             &
      direction * rounding_factor(:) * frac_flux_tmp(:)                        &
    )
  end do  ! DoF

end subroutine hori_dep_dist_ffsl_sphere_1d

end module hori_dep_dist_ffsl_sphere_kernel_mod
