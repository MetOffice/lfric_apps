!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute the interpolation weights and indices to remap a W3 or Wtheta
!!        scalar field from the standard cubed-sphere mesh to another panel at
!!        the edges of the mesh panels.
module panel_edge_weights_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type, func_type,                          &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_LOGICAL, GH_INTEGER,              &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 ANY_DISCONTINUOUS_SPACE_9,                    &
                                 GH_BASIS, GH_EVALUATOR,                       &
                                 CELL_COLUMN, STENCIL, CROSS2D
use constants_mod,         only: r_tran, r_def, i_def, l_def, LARGE_REAL_POSITIVE
use fs_continuity_mod,     only: Wchi
use reference_element_mod, only: W, S, N, E

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_edge_weights_kernel_type
  private
  type(arg_type) :: meta_args(11) = (/                                         &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  Wchi, STENCIL(CROSS2D)),     &
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),  &
       arg_type(GH_FIELD*2, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,   &
                                                        STENCIL(CROSS2D)),     &
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                              &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                               &
  /)
  type(func_type) :: meta_funcs(2) = (/                                        &
      func_type(Wchi, GH_BASIS),                                               &
      func_type(ANY_DISCONTINUOUS_SPACE_9, GH_BASIS)                           &
  /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: panel_edge_weights_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_edge_weights_code

contains

!> @brief Compute the interpolation weights and indices to remap a W3 or Wtheta
!!        scalar field from the standard cubed-sphere mesh to another panel at
!!        the edges of the mesh panels.
!> @param[in]     nlayers               Number of layers
!> @param[in,out] panel_edge_weights_x  Remapping interpolation weights for x
!> @param[in,out] panel_edge_weights_y  Remapping interpolation weights for y
!> @param[in,out] panel_edge_indices_x  Remapping interpolation indices for x
!> @param[in,out] panel_edge_indices_y  Remapping interpolation indices for y
!> @param[in]     alpha                 Alpha coordinate
!> @param[in]     beta                  Beta coordinate
!> @param[in]     height                Height coordinate
!> @param[in]     wx_stencil_size       Size of the Wchi stencil (num cells)
!> @param[in]     wx_stencil            DoF map for the Wchi-stencil
!> @param[in]     wx_max_length         Maximum stencil branch length
!> @param[in]     alpha_x               2D alpha coord for mesh extended in x
!> @param[in]     beta_x                2D beta coord for mesh extended in x
!> @param[in]     alpha_y               2D alpha coord for mesh extended in y
!> @param[in]     beta_y                2D beta coord for mesh extended in y
!> @param[in]     panel_id              ID of cubed sphere panel for each cell
!> @param[in]     pid_stencil_size      Size of panel ID stencil (num cells)
!> @param[in]     pid_stencil           DoF map for the panelid-stencil
!> @param[in]     pid_max_length        Maximum stencil branch length
!> @param[in]     panel_id              Field containing IDs of mesh panels
!> @param[in]     panel_edge_dist_W     2D field containing the distance of each
!!                                      column from the panel edge to the West
!> @param[in]     panel_edge_dist_E     2D field containing the distance of each
!!                                      column from the panel edge to the East
!> @param[in]     panel_edge_dist_S     2D field containing the distance of each
!!                                      column from the panel edge to the South
!> @param[in]     panel_edge_dist_N     2D field containing the distance of each
!!                                      column from the panel edge to the North
!> @param[in]     cubic_remap           Whether to use cubic remapping. If false
!!                                      linear remapping is used
!> @param[in]     remap_depth           Depth to which the remapping may happen
!> @param[in]     ndf_ww                Num DoFs per cell for remapping weights
!> @param[in]     undf_ww               Num DoFs for this partition for weights
!> @param[in]     map_ww                DoF map for remapping weights
!> @param[in]     ndf_wx                Num DoFs per cell for coord fields
!> @param[in]     undf_wx               Num DoFs for this partition for coords
!> @param[in]     map_wx                DoF map for coord fields
!> @param[in]     basis_wx              Basis functions of coordinate fields,
!!                                      evaluated at nodes of scalar fields
!> @param[in]     ndf_wx_2d             Num DoFs per cell for 2D coord fields
!> @param[in]     undf_wx_2d            Num DoFs in this partition for 2D coords
!> @param[in]     map_wx_2d             DoF map for 2D coord fields
!> @param[in]     basis_wx_2d           Basis functions of 2D coordinate fields,
!!                                      evaluated at nodes of scalar fields
!> @param[in]     ndf_pid               Num DoFs per cell for panel ID field
!> @param[in]     undf_pid              Num DoFs for this partition for panel ID
!> @param[in]     map_pid               DoF map for panel ID field
subroutine panel_edge_weights_code( nlayers,                                   &
                                    panel_edge_weights_x,                      &
                                    panel_edge_weights_y,                      &
                                    panel_edge_indices_x,                      &
                                    panel_edge_indices_y,                      &
                                    alpha, beta, height,                       &
                                    wx_stencil_size,                           &
                                    wx_max_length,                             &
                                    wx_stencil_map,                            &
                                    alpha_x, beta_x,                           &
                                    alpha_y, beta_y,                           &
                                    panel_id,                                  &
                                    pid_stencil_size,                          &
                                    pid_max_length,                            &
                                    pid_stencil_map,                           &
                                    panel_edge_dist_W,                         &
                                    panel_edge_dist_E,                         &
                                    panel_edge_dist_S,                         &
                                    panel_edge_dist_N,                         &
                                    cubic_remap,                               &
                                    remap_depth,                               &
                                    ndf_ww, undf_ww, map_ww,                   &
                                    ndf_wx, undf_wx, map_wx, basis_wx,         &
                                    ndf_wx_2d, undf_wx_2d, map_wx_2d,          &
                                    basis_wx_2d,                               &
                                    ndf_pid, undf_pid, map_pid                 &
  )

  use panel_edge_support_mod, only: panel_neighbour

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: wx_stencil_size(4)
  integer(kind=i_def), intent(in)    :: pid_stencil_size(4)
  integer(kind=i_def), intent(in)    :: wx_max_length
  integer(kind=i_def), intent(in)    :: pid_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: remap_depth
  logical(kind=l_def), intent(in)    :: cubic_remap

  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wx_stencil_map(ndf_wx, wx_max_length, 4)
  integer(kind=i_def), intent(in)    :: pid_stencil_map(ndf_pid, pid_max_length, 4)
  real(kind=r_def),    intent(in)    :: basis_wx(1, ndf_wx, ndf_ww)
  real(kind=r_def),    intent(in)    :: basis_wx_2d(1, ndf_wx_2d, ndf_ww)

  real(kind=r_tran),   intent(inout) :: panel_edge_weights_x(undf_ww)
  real(kind=r_tran),   intent(inout) :: panel_edge_weights_y(undf_ww)
  integer(kind=i_def), intent(inout) :: panel_edge_indices_x(undf_ww)
  integer(kind=i_def), intent(inout) :: panel_edge_indices_y(undf_ww)
  real(kind=r_def),    intent(in)    :: alpha(undf_wx)
  real(kind=r_def),    intent(in)    :: beta(undf_wx)
  real(kind=r_def),    intent(in)    :: height(undf_wx)
  real(kind=r_def),    intent(in)    :: alpha_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_x(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: alpha_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)

  ! Internal variables
  integer(kind=i_def) :: pid_idx, owned_panel, swapped_panel_x, swapped_panel_y

  ! Panel id for this column
  pid_idx = map_pid(1)

  swapped_panel_x = 0
  swapped_panel_y = 0

  ! Compute mapping weights ----------------------------------------------------
  if (ABS(panel_edge_dist_W(pid_idx)) <= remap_depth                           &
      .or. ABS(panel_edge_dist_E(pid_idx)) <= remap_depth) then

    owned_panel = int(panel_id(pid_idx), i_def)

    if (ABS(panel_edge_dist_W(pid_idx)) < ABS(panel_edge_dist_E(pid_idx))) then
      ! If the W edge is closer, remap to the W panel
      swapped_panel_x = panel_neighbour(owned_panel, W)
    else
      ! If the E edge is closer, remap to the E panel
      swapped_panel_x = panel_neighbour(owned_panel, E)
    end if

    ! Compute remapping weights and indices for x direction
    call panel_edge_weights_1d(                                                &
            panel_edge_weights_x, panel_edge_indices_x,                        &
            owned_panel, swapped_panel_x,                                      &
            alpha, beta, height, wx_max_length,                                &
            ! Indices S and N, because crossing the x edge of a panel, and
            ! performing the remapping parallel to the edge
            wx_stencil_size(S), wx_stencil_map(:,:,S),                         &
            wx_stencil_size(N), wx_stencil_map(:,:,N),                         &
            alpha_x, beta_x,                                                   &
            panel_id, pid_max_length,                                          &
            pid_stencil_size(S), pid_stencil_map(:,:,S),                       &
            pid_stencil_size(N), pid_stencil_map(:,:,N),                       &
            cubic_remap,                                                       &
            ndf_ww, undf_ww, map_ww,                                           &
            ndf_wx, undf_wx, map_wx, basis_wx,                                 &
            ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,                     &
            ndf_pid, undf_pid, map_pid                                         &
    )
  end if

  if (ABS(panel_edge_dist_S(pid_idx)) <= remap_depth                           &
      .or. ABS(panel_edge_dist_N(pid_idx)) <= remap_depth) then

    owned_panel = int(panel_id(pid_idx), i_def)

    if (ABS(panel_edge_dist_S(pid_idx)) < ABS(panel_edge_dist_N(pid_idx))) then
      ! If the S edge is closer, remap to the S panel
      swapped_panel_y = panel_neighbour(owned_panel, S)
    else
      ! If the N edge is closer, remap to the N panel
      swapped_panel_y = panel_neighbour(owned_panel, N)
    end if

    ! Compute remapping weights and indices for y direction
    call panel_edge_weights_1d(                                                &
            panel_edge_weights_y, panel_edge_indices_y,                        &
            owned_panel, swapped_panel_y,                                      &
            alpha, beta, height, wx_max_length,                                &
            ! Indices W and E, because crossing the y edge of a panel, and
            ! performing the remapping parallel to the edge
            wx_stencil_size(W), wx_stencil_map(:,:,W),                         &
            wx_stencil_size(E), wx_stencil_map(:,:,E),                         &
            alpha_y, beta_y, panel_id, pid_max_length,                         &
            pid_stencil_size(W), pid_stencil_map(:,:,W),                       &
            pid_stencil_size(E), pid_stencil_map(:,:,E),                       &
            cubic_remap,                                                       &
            ndf_ww, undf_ww, map_ww,                                           &
            ndf_wx, undf_wx, map_wx, basis_wx,                                 &
            ndf_wx_2d, undf_wx_2d, map_wx_2d, basis_wx_2d,                     &
            ndf_pid, undf_pid, map_pid                                         &
    )
  end if

end subroutine panel_edge_weights_code

!> @brief Private routine for calculation of edge remapping weights and
!!        indices for *a single direction* (x or y)
subroutine panel_edge_weights_1d( remap_weights, remap_indices,                &
                                  owned_panel, swapped_panel,                  &
                                  alpha, beta, height, wx_max_length,          &
                                  wx_stencil_size_l, wx_stencil_l,             &
                                  wx_stencil_size_r, wx_stencil_r,             &
                                  alpha_ext, beta_ext,                         &
                                  panel_id, pid_max_length,                    &
                                  pid_stencil_size_l, pid_stencil_l,           &
                                  pid_stencil_size_r, pid_stencil_r,           &
                                  cubic_remap,                                 &
                                  ndf_ww, undf_ww, map_ww,                     &
                                  ndf_wx, undf_wx, map_wx,                     &
                                  basis_wx,                                    &
                                  ndf_wx_2d, undf_wx_2d, map_wx_2d,            &
                                  basis_wx_2d,                                 &
                                  ndf_pid, undf_pid, map_pid                   &
)

  use sci_chi_transform_mod, only: chi2xyz, get_to_rotate, get_to_stretch,     &
                                   get_inverse_mesh_rotation_matrix,           &
                                   get_stretch_factor
  use coord_transform_mod,   only: alphabetar2xyz, xyz2alphabetar,             &
                                   inverse_schmidt_transform_xyz

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: wx_stencil_size_l
  integer(kind=i_def), intent(in)    :: wx_stencil_size_r
  integer(kind=i_def), intent(in)    :: wx_max_length
  integer(kind=i_def), intent(in)    :: pid_stencil_size_l
  integer(kind=i_def), intent(in)    :: pid_stencil_size_r
  integer(kind=i_def), intent(in)    :: pid_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in)    :: ndf_wx_2d, undf_wx_2d
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: owned_panel
  integer(kind=i_def), intent(in)    :: swapped_panel
  logical(kind=l_def), intent(in)    :: cubic_remap
  ! Maps
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wx_stencil_l(ndf_wx, wx_stencil_size_l)
  integer(kind=i_def), intent(in)    :: wx_stencil_r(ndf_wx, wx_stencil_size_r)
  integer(kind=i_def), intent(in)    :: pid_stencil_l(ndf_pid, pid_stencil_size_l)
  integer(kind=i_def), intent(in)    :: pid_stencil_r(ndf_pid, pid_stencil_size_r)
  ! Fields
  real(kind=r_def),    intent(in)    :: basis_wx(1, ndf_wx, ndf_ww)
  real(kind=r_def),    intent(in)    :: basis_wx_2d(1, ndf_wx_2d, ndf_ww)
  real(kind=r_def),    intent(in)    :: alpha(undf_wx)
  real(kind=r_def),    intent(in)    :: beta(undf_wx)
  real(kind=r_def),    intent(in)    :: height(undf_wx)
  real(kind=r_def),    intent(in)    :: alpha_ext(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: beta_ext(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_tran),   intent(inout) :: remap_weights(undf_ww)
  integer(kind=i_def), intent(inout) :: remap_indices(undf_ww)

  ! Internal variables
  integer(kind=i_def), parameter   :: interp_dir_alpha = 1
  integer(kind=i_def), parameter   :: interp_dir_beta  = 2
  real(kind=r_def),    parameter   :: unit_radius = 1.0_r_def
  real(kind=r_def),    allocatable :: x1(:), dx(:)

  integer(kind=i_def) :: panel, panel_edge
  integer(kind=i_def) :: df, id1, id2, id3, id4, n
  integer(kind=i_def) :: interp_dir
  integer(kind=i_def) :: ncells_in_stencil
  integer(kind=i_def) :: wx_stencil_1d(ndf_wx, 2*wx_max_length-1)
  integer(kind=i_def) :: pid_stencil_1d(2*pid_max_length-1)
  real(kind=r_def)    :: abh(3)
  real(kind=r_def)    :: x0, w1, w2, w3, w4
  real(kind=r_def)    :: xyz(3), h_dummy
  real(kind=r_def)    :: alpha_owned, beta_owned
  real(kind=r_def)    :: stretch_factor
  real(kind=r_def)    :: inverse_rot_matrix(3,3)
  logical(kind=l_def) :: to_rotate, to_stretch

  to_rotate = get_to_rotate()
  to_stretch = get_to_stretch()
  if (to_rotate) then
    inverse_rot_matrix = get_inverse_mesh_rotation_matrix()
  end if
  if (to_stretch) then
    stretch_factor = get_stretch_factor()
  end if

  ! -------------------------------------------------------------------------- !
  ! Determine direction of interpolation (alpha or beta)
  ! -------------------------------------------------------------------------- !
  panel_edge = 10*owned_panel + swapped_panel
  select case (panel_edge)
  case (15, 26, 34, 43, 51, 62)
    ! alpha direction on both panels
    interp_dir = interp_dir_alpha
  case (12, 21, 36, 45, 54, 63)
    ! beta direction on both panels
    interp_dir = interp_dir_beta
  case (16, 41, 32, 25, 64, 53)
    ! alpha on owned panel, beta on swapped panel
    interp_dir = interp_dir_alpha
  case (61, 14, 23, 52, 46, 35)
    ! beta on owned panel, alpha on swapped panel
    interp_dir = interp_dir_beta
  end select

  ! -------------------------------------------------------------------------- !
  ! Create a 1D stencil
  ! -------------------------------------------------------------------------- !
  ! We have the x (or y) parts of the 2D cross stencil. Positive and negative
  ! directions are stored in separate parts of the array, so here we unify them
  ! e.g. we may have the x-part of the 2D cross stencil of the form:
  !  | 2 | 1 | 4 | 6 |
  ! Stored as stencil =[
  ! 1, 2;
  ! 1, 4, 6]
  ! We now extract this to be = [1, 2, 4, 6] to match 1D stencils

  ncells_in_stencil = wx_stencil_size_l + wx_stencil_size_r - 1

  wx_stencil_1d = 0
  do n = 1, wx_stencil_size_l
    wx_stencil_1d(:,n) = wx_stencil_l(:,n)
  end do
  ! Omit 1 from the second part, as the central cell is already included
  do n = 1, wx_stencil_size_r - 1
    wx_stencil_1d(:,n+wx_stencil_size_l) = wx_stencil_r(:,n+1)
  end do

  pid_stencil_1d = 0
  do n = 1, pid_stencil_size_l
    pid_stencil_1d(n) = pid_stencil_l(1,n)
  end do
  do n = 1, pid_stencil_size_r - 1
    pid_stencil_1d(n+pid_stencil_size_l) = pid_stencil_r(1,n+1)
  end do

  ! -------------------------------------------------------------------------- !
  ! Compute 1D locations of points used in interpolation
  ! -------------------------------------------------------------------------- !
  ! Compute interpolation point (x0) in centre of remapped cell
  abh = 0.0_r_def
  do df = 1, ndf_wx_2d
    abh(1) = abh(1) + alpha_ext(map_wx_2d(df))*basis_wx_2d(1,df,1)
    abh(2) = abh(2) + beta_ext(map_wx_2d(df))*basis_wx_2d(1,df,1)
  end do

  ! Convert (alpha,beta,r) on panel swapped_panel to xyz
  call alphabetar2xyz(                                                         &
          abh(1), abh(2), unit_radius, swapped_panel, xyz(1), xyz(2), xyz(3)   &
  )
  ! Now convert back to (alpha, beta, r) on owned panel
  ! As (alpha,beta) are already native coordinates, no need to consider
  ! rotation or stretching here, as this would be undone in the next step
  call xyz2alphabetar(                                                         &
          xyz(1), xyz(2), xyz(3), owned_panel, abh(1), abh(2), h_dummy         &
  )

  x0 = abh(interp_dir)

  ! Now compute all points on original mesh in coordinates of extended mesh
  allocate( x1(ncells_in_stencil), dx(ncells_in_stencil) )

  do n = 1, ncells_in_stencil

    panel = int(panel_id(pid_stencil_1d(n)), i_def)

    if ( owned_panel /= panel ) then
      ! If point is on a different panel to the swapped panel (i.e we have gone
      ! 'around' the cubed sphere corner) then discard that point, by setting
      ! the distance to some large number
      x1(n) = LARGE_REAL_POSITIVE
    else
      ! Otherwise compute the coordinates of the point in terms of the
      ! owned_panel
      abh = 0.0_r_def
      do df = 1, ndf_wx
        ! Convert chi fields to native cubed-sphere coordinates, on owned panel
        call chi2xyz(                                                          &
                alpha(wx_stencil_1d(df, n)), beta(wx_stencil_1d(df, n)),       &
                unit_radius, owned_panel, xyz(1), xyz(2), xyz(3)               &
        )
        ! Transform to the Cartesian coordinates in the *native* coordinate
        ! system by applying the inverse of any mesh rotation and stretching:
        if (to_rotate) then
          xyz = matmul(inverse_rot_matrix, xyz)
        end if
        if (to_stretch) then
          xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
        end if
        call xyz2alphabetar(                                                   &
                xyz(1), xyz(2), xyz(3), owned_panel,                           &
                alpha_owned, beta_owned, h_dummy                               &
        )
        abh(1) = abh(1) + alpha_owned * basis_wx(1,df,1)
        abh(2) = abh(2) + beta_owned * basis_wx(1,df,1)
      end do

      ! Fill array with the points to use for the interpolation
      x1(n) = abh(interp_dir)
    end if
  end do

  ! -------------------------------------------------------------------------- !
  ! Calculate interpolation weights
  ! -------------------------------------------------------------------------- !
  ! Distance between interpolation point x0 and data points x1
  dx = ABS(x1 - x0)
  ! Now find two points (for linear interpolation) in x1 that are closest to x0
  id1 = MINLOC(dx, 1)
  dx(id1) = LARGE_REAL_POSITIVE
  id2 = MINLOC(dx, 1)

  ! Compute weights for cubic interpolation ------------------------------------
  if (cubic_remap) then
    ! First find the next two closest points in stencil
    dx(id2) = LARGE_REAL_POSITIVE
    id3 = MINLOC(dx, 1)
    dx(id3) = LARGE_REAL_POSITIVE
    id4 = MINLOC(dx, 1)

    ! Now compute the weights
    w1 = (x0 - x1(id2))/(x1(id1) - x1(id2))                                    &
          * (x0 - x1(id3))/(x1(id1) - x1(id3))                                 &
          * (x0 - x1(id4))/(x1(id1) - x1(id4))
    w2 = (x0 - x1(id1))/(x1(id2) - x1(id1))                                    &
          * (x0 - x1(id3))/(x1(id2) - x1(id3))                                 &
          * (x0 - x1(id4))/(x1(id2) - x1(id4))
    w3 = (x0 - x1(id1))/(x1(id3) - x1(id1))                                    &
          * (x0 - x1(id2))/(x1(id3) - x1(id2))                                 &
          * (x0 - x1(id4))/(x1(id3) - x1(id4))
    w4 = (x0 - x1(id1))/(x1(id4) - x1(id1))                                    &
          * (x0 - x1(id2))/(x1(id4) - x1(id2))                                 &
          * (x0 - x1(id3))/(x1(id4) - x1(id3))

    remap_weights(map_ww(1)) = REAL(w1, r_tran)
    remap_weights(map_ww(1)+1) = REAL(w2, r_tran)
    remap_weights(map_ww(1)+2) = REAL(w3, r_tran)
    remap_weights(map_ww(1)+3) = REAL(w4, r_tran)

    remap_indices(map_ww(1)) = id1
    remap_indices(map_ww(1)+1) = id2
    remap_indices(map_ww(1)+2) = id3
    remap_indices(map_ww(1)+3) = id4

  else
    ! Compute weights for linear interpolation ---------------------------------
    w1 = (x0 - x1(id2))/(x1(id1) - x1(id2))
    w2 = (x0 - x1(id1))/(x1(id2) - x1(id1))

    remap_weights(map_ww(1)) = REAL(w1, r_tran)
    remap_weights(map_ww(1)+1) = REAL(w2, r_tran)

    remap_indices(map_ww(1)) = id1
    remap_indices(map_ww(1)+1) = id2
  end if

  deallocate( x1, dx )

end subroutine panel_edge_weights_1d

end module panel_edge_weights_kernel_mod
