!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Computes the coordinates corresponding to extending the neighbouring
!!        cubed-sphere panel into this cell
module panel_edge_coords_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_INTEGER,                          &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 ANY_DISCONTINUOUS_SPACE_7,                    &
                                 OWNED_AND_HALO_CELL_COLUMN
use constants_mod,         only: r_def, i_def, l_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_edge_coords_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                          &
       arg_type(GH_FIELD*2, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_7),  &
       arg_type(GH_FIELD*2, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_7),  &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                               &
  /)
  integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
contains
  procedure, nopass :: panel_edge_coords_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_edge_coords_code

contains

!> @brief Extend the equiangular coordinate fields from panel boundaries into
!!        the haloes
!> @param[in]     nlayers            Number of layers
!> @param[in,out] alpha_x            Alpha coordinate field extended in x
!> @param[in,out] beta_x             Beta coordinate field extended in x
!> @param[in,out] alpha_y            Alpha coordinate field extended in y
!> @param[in,out] beta_y             Beta coordinate field extended in y
!> @param[in]     chi_1              First coodinate field
!> @param[in]     chi_2              Second coodinate field
!> @param[in]     chi_3              Third coodinate field
!> @param[in]     panel_id           ID of the panel for each cell column
!> @param[in]     panel_edge_dist_W  2D field containing the distance of each
!!                                   column from the panel edge to the West
!> @param[in]     panel_edge_dist_E  2D field containing the distance of each
!!                                   column from the panel edge to the East
!> @param[in]     panel_edge_dist_S  2D field containing the distance of each
!!                                   column from the panel edge to the South
!> @param[in]     panel_edge_dist_N  2D field containing the distance of each
!!                                   column from the panel edge to the North
!> @param[in]     stencil_extent     Max stencil extent to be used
!> @param[in]     ndf_wx_2d          Num DoFs per cell for 2D coord fields
!> @param[in]     undf_wx_2d         Num DoFs for this partition for 2D coords
!> @param[in]     map_wx_2d          DoFmap for 2D coord fields
!> @param[in]     ndf_wx             Num DoFs per cell for 3D coord fields
!> @param[in]     undf_wx            Num DoFs for this partition for 3D coords
!> @param[in]     map_wx             DoFmap for 3D coord fields
!> @param[in]     ndf_pid            Num DoFs per cell for panel ID
!> @param[in]     undf_pid           Num DoFs for this partition for panel ID
!> @param[in]     map_pid            DoFmap for panel ID
subroutine panel_edge_coords_code( nlayers,                                    &
                                   alpha_x, beta_x, alpha_y, beta_y,           &
                                   chi_1, chi_2, chi_3, panel_id,              &
                                   panel_edge_dist_W, panel_edge_dist_E,       &
                                   panel_edge_dist_S, panel_edge_dist_N,       &
                                   stencil_extent,                             &
                                   ndf_wx_2d, undf_wx_2d, map_wx_2d,           &
                                   ndf_wx, undf_wx, map_wx,                    &
                                   ndf_pid, undf_pid, map_pid                  &
                                  )

  use panel_edge_support_mod, only: panel_neighbour
  use reference_element_mod,  only: W, S, N, E

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_wx, ndf_pid, ndf_wx_2d
  integer(kind=i_def), intent(in)    :: undf_wx, undf_pid, undf_wx_2d
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: stencil_extent
  real(kind=r_def),    intent(inout) :: alpha_x(undf_wx_2d)
  real(kind=r_def),    intent(inout) :: beta_x(undf_wx_2d)
  real(kind=r_def),    intent(inout) :: alpha_y(undf_wx_2d)
  real(kind=r_def),    intent(inout) :: beta_y(undf_wx_2d)
  real(kind=r_def),    intent(in)    :: chi_1(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_2(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_3(undf_wx)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)

  integer(kind=i_def) :: owned_panel, swapped_panel_x, swapped_panel_y
  integer(kind=i_def) :: panel_W, panel_E, panel_S, panel_N

  ! Panel id for this column
  owned_panel = int(panel_id(map_pid(1)), i_def)

  panel_W = panel_neighbour(owned_panel, W)
  panel_E = panel_neighbour(owned_panel, E)
  panel_S = panel_neighbour(owned_panel, S)
  panel_N = panel_neighbour(owned_panel, N)

  ! Determine if we are near the edge of a panel -------------------------------
  ! Initialise swapped panels to zero
  swapped_panel_x = 0
  swapped_panel_y = 0

  ! At the moment we can only set one neighbour for the x/y directions, so use
  ! whichever edge is closest
  if (ABS(panel_edge_dist_W(map_pid(1))) <= stencil_extent) then
    swapped_panel_x = panel_W
  end if
  if (ABS(panel_edge_dist_E(map_pid(1))) <= stencil_extent .and.               &
      ! Take panel to the E if this edge is closer
      ABS(panel_edge_dist_E(map_pid(1))) < ABS(panel_edge_dist_W(map_pid(1)))) then
    swapped_panel_x = panel_E
  end if

  if (ABS(panel_edge_dist_S(map_pid(1))) <= stencil_extent) then
    swapped_panel_y = panel_S
  end if
  if (ABS(panel_edge_dist_N(map_pid(1))) <= stencil_extent .and.               &
      ! Take panel to the N if this edge is closer
      ABS(panel_edge_dist_N(map_pid(1))) < ABS(panel_edge_dist_S(map_pid(1)))) then
    swapped_panel_y = panel_N
  end if

  ! Now extend the coordinates -------------------------------------------------
  if ( swapped_panel_x /= 0_i_def ) then
    call panel_edge_coords_1d( alpha_x, beta_x, chi_1, chi_2, chi_3,           &
                               owned_panel, swapped_panel_x,                   &
                               ndf_wx_2d, undf_wx_2d, map_wx_2d,               &
                               ndf_wx, undf_wx, map_wx                         &
    )
  end if

  if ( swapped_panel_y /= 0_i_def ) then
    call panel_edge_coords_1d( alpha_y, beta_y, chi_1, chi_2, chi_3,           &
                               owned_panel, swapped_panel_y,                   &
                               ndf_wx_2d, undf_wx_2d, map_wx_2d,               &
                               ndf_wx, undf_wx, map_wx                         &
    )
  end if

end subroutine panel_edge_coords_code

!> @brief Private routine to perform the extension of the cubed-sphere coords
!!        for a single direction. This reduces code duplication.
subroutine panel_edge_coords_1d( alpha, beta, chi_1, chi_2, chi_3,             &
                                 owned_panel, swapped_panel,                   &
                                 ndf_wx_2d, undf_wx_2d, map_wx_2d,             &
                                 ndf_wx, undf_wx, map_wx                       &
                               )

  use sci_chi_transform_mod, only: chi2xyz, get_to_rotate, get_to_stretch,     &
                                   get_inverse_mesh_rotation_matrix,           &
                                   get_stretch_factor
  use coord_transform_mod,   only: xyz2alphabetar, inverse_schmidt_transform_xyz

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: ndf_wx, ndf_wx_2d
  integer(kind=i_def), intent(in)    :: undf_wx, undf_wx_2d
  integer(kind=i_def), intent(in)    :: map_wx(ndf_wx)
  integer(kind=i_def), intent(in)    :: map_wx_2d(ndf_wx_2d)
  real(kind=r_def),    intent(in)    :: chi_1(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_2(undf_wx)
  real(kind=r_def),    intent(in)    :: chi_3(undf_wx)
  real(kind=r_def),    intent(inout) :: alpha(undf_wx_2d)
  real(kind=r_def),    intent(inout) :: beta(undf_wx_2d)
  integer(kind=i_def), intent(in)    :: owned_panel, swapped_panel

  ! Local variables
  integer(kind=i_def) :: df
  integer(kind=i_def) :: panel_edge
  real(kind=r_def)    :: h_dummy, xyz(3)
  real(kind=r_def)    :: alpha_swapped(ndf_wx), beta_swapped(ndf_wx)
  real(kind=r_def)    :: alpha_owned(ndf_wx), beta_owned(ndf_wx)
  real(kind=r_def)    :: alpha_extended(ndf_wx), beta_extended(ndf_wx)
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

  ! Calculate equiangular coordinates for both the owned panel of this column,
  ! and the panel of its neighbour
  ! Fill small local arrays with the coordinates
  do df = 1, ndf_wx
    ! Ignore height coordinate as this is not needed
    call chi2xyz(                                                              &
            chi_1(map_wx(df)), chi_2(map_wx(df)), chi_3(map_wx(df)),           &
            owned_panel, xyz(1), xyz(2), xyz(3)                                &
    )
    ! Transform to the Cartesian coordinates in the *native* coordinate system
    ! by applying the inverse of any mesh rotation and stretching:
    if (to_rotate) then
      xyz = matmul(inverse_rot_matrix, xyz)
    end if
    if (to_stretch) then
      xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
    end if
    ! Convert back to equiangular coordinates, both for the system on this panel
    ! but also the system on the neighbouring panel
    call xyz2alphabetar(                                                       &
            xyz(1), xyz(2), xyz(3), owned_panel,                               &
            alpha_owned(df), beta_owned(df), h_dummy                           &
    )
    call xyz2alphabetar(                                                       &
            xyz(1), xyz(2), xyz(3), swapped_panel,                             &
            alpha_swapped(df), beta_swapped(df), h_dummy                       &
    )
  end do

  ! alpha_swapped now contains the (alpha,beta) coords of the neighbouring panel
  ! alpha_owned contains the same physical points but in coords of owned panel

  ! To obtain the extended coordinates we take one component of alpha/beta_owned
  ! and one component of alpha/beta_swapped
  ! e.g between panels 1 & 2 we have
  ! ------------------------------
  ! |             |              |
  ! | beta        | beta         |
  ! | ^     1     | ^     2      |
  ! | |           | |            |
  ! | --> alpha   | --> alpha    |
  ! |             |              |
  ! ------------------------------
  ! and so when extending panel 1 we set alpha_extended to be the alpha
  ! coordinate from panel 1 (alpha_extended = alpha_owned)
  ! and the beta coordinate is from the beta coordinate on panel 2
  ! (beta_extended = beta_swapped)

  ! Now extend coordinates in the halo regions
  panel_edge = 10*owned_panel + swapped_panel
  select case ( panel_edge )
  case(12, 21, 36, 45, 54, 63)
    ! 12: EAST edge of panel 1 to 2
    ! 21: WEST edge of panel 2 to 1
    ! 36: SOUTH edge of panel 3 to 6
    ! 45: NORTH edge of panel 4 to 5
    ! 54: WEST edge of panel 5 to 4
    ! 63: EAST edge of panel 6 to 3
    alpha_extended = alpha_swapped
    beta_extended = beta_owned
  case(15, 26, 34, 43, 51, 62)
    ! 15: NORTH edge of panel 1 to 5
    ! 26: SOUTH edge of panel 2 to 6
    ! 34: WEST edge of panel 3 to 4
    ! 43: WEST edge of panel 4 to 3
    ! 51: SOUTH edge of panel 5 to 1
    ! 62: NORTH edge of panel 6 to 2
    alpha_extended = alpha_owned
    beta_extended = beta_swapped
  case(14, 23, 35, 46, 52, 61)
    ! 14: WEST edge of panel 1 to 4
    ! 23: EAST edge of panel 2 to 3
    ! 35: NORTH edge of panel 3 to 5
    ! 46: SOUTH edge of panel 4 to 6
    ! 52: EAST edge of panel 5 to 2
    ! 61: WEST edge of panel 6 to 1
    alpha_extended = beta_owned
    beta_extended = beta_swapped
  case(16, 25, 32, 41, 53, 64)
    ! 16: SOUTH edge of panel 1 to 6
    ! 25: NORTH edge of panel 2 to 5
    ! 32: EAST edge of panel 3 to 2
    ! 41: EAST edge of panel 4 to 1
    ! 53: NORTH edge of panel 5 to 3
    ! 64: SOUTH edge of panel 6 to 4
    beta_extended = alpha_owned
    alpha_extended = alpha_swapped
  end select

  ! Write back to coordinate field
  do df = 1, ndf_wx
    alpha(map_wx_2d(df)) = alpha_extended(df)
    beta(map_wx_2d(df)) = beta_extended(df)
  end do

end subroutine panel_edge_coords_1d

end module panel_edge_coords_kernel_mod
