!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains routines useful for identifying whether a calculation will
!!        cross a panel edge
module panel_edge_support_mod


use constants_mod,         only : i_def, l_def
use reference_element_mod, only : W, S, N, E

implicit none

private

public :: FAR_AWAY
public :: crosses_panel_edge
public :: crosses_rotated_panel_edge
public :: rotated_panel_neighbour
public :: panel_neighbour

integer(kind=i_def), parameter :: FAR_AWAY = 10000

contains

!> @brief Determines whether a calculation will cross the edge of a cubed sphere
!!        panel. This is done by checking the distances to the panel edges, and
!!        determining the maximum stencil extent needed for this calculation
!> @param[in]  edge_dist_left       Distance to the left panel edge
!> @param[in]  edge_dist_right      Distance to the right panel edge
!> @param[in]  extent_size          Size of the stencil extent
!> @param[in]  order                Polynomial order used in calculation
!> @param[in]  num_faces_to_check   Number of faces to check calculation for.
!!                                  This corresponds to the values in the
!!                                  face_selector fields. If 1, only do the left
!!                                  face, if 2 do both left and right faces.
!> @param[in]  only_rotated         If true, only check for rotated panel edges.
!!                                  If false, check for all panel edges.
!> @param[in]  panel_id             The ID of the panel for this column
!> @param[in]  direction            The direction this calculation is performed
!!                                  in. If 1, then the stencil is in x, if 2
!!                                  then it is in the y direction.
!> @param[in]  dof_L                Local DoF index for the left face
!> @param[in]  dof_R                Local DoF index for the right face
!> @param[in]  dep_highest_k        Highest k index (layer) for each column to
!!                                  be used in the integer part of a flux. It
!!                                  is used to determine the maximum departure
!!                                  distance needed for this calculation
!> @param[in]  ndf_depk             Number of DoFs per cell for dep_highest_k
!> @param[in]  undf_depk            Length of the dep_highest_k array
!> @param[in]  map_depk             DoFmap for dep_highest_k
!> @return                          Flag to indicate whether calculation will
!!                                  cross the panel edge
function crosses_panel_edge(edge_dist_left, edge_dist_right,                   &
                            extent_size, order, num_faces_to_check,            &
                            only_rotated, panel_id, direction,                 &
                            dof_L, dof_R, dep_highest_k, ndep,                 &
                            ndf_depk, undf_depk, map_depk)                     &
                            result(cross_edge_flag)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: edge_dist_left
  integer(kind=i_def), intent(in) :: edge_dist_right
  integer(kind=i_def), intent(in) :: extent_size
  integer(kind=i_def), intent(in) :: order
  integer(kind=i_def), intent(in) :: num_faces_to_check
  logical(kind=l_def), intent(in) :: only_rotated
  integer(kind=i_def), intent(in) :: panel_id
  integer(kind=i_def), intent(in) :: direction
  integer(kind=i_def), intent(in) :: dof_L
  integer(kind=i_def), intent(in) :: dof_R
  integer(kind=i_def), intent(in) :: ndep
  integer(kind=i_def), intent(in) :: ndf_depk
  integer(kind=i_def), intent(in) :: undf_depk
  integer(kind=i_def), intent(in) :: dep_highest_k(undf_depk)
  integer(kind=i_def), intent(in) :: map_depk(ndf_depk)

  ! Local variables
  integer(kind=i_def) :: closest_edge, max_int_dep_dist, dep_idx, face_dir
  logical(kind=l_def) :: cross_edge_flag, consider_left, consider_right

  cross_edge_flag = .false.

  ! 1) Check if there is an edge somewhere within the stencil
  consider_left = edge_dist_left > 0 .and. edge_dist_left < extent_size + 1
  consider_right = edge_dist_right > 0 .and. edge_dist_right < extent_size + 1

  ! 2) Take into account whether rotation needs considering for this panel edge
  if (only_rotated) then
    consider_left = consider_left                                              &
      .and. ABS(rotated_panel_neighbour(panel_id, direction)) > 0

    face_dir = 2_i_def + direction
    consider_right = consider_right                                            &
      .and. ABS(rotated_panel_neighbour(panel_id, face_dir)) > 0
  end if

  ! 3) There is a panel edge somewhere in the stencil. Now departure if the
  !    calculation will involve crossing the panel edge
  if (consider_left .or. consider_right) then

    ! a) If we are very close to the edge then some part of the calculation will
    !    definitely cross the panel edge
    if (consider_left .and. consider_right) then
      closest_edge = MIN(ABS(edge_dist_left), ABS(edge_dist_right))
    else if (consider_left) then
      closest_edge = edge_dist_left
    else
      closest_edge = edge_dist_right
    end if

    if (closest_edge < order + 2) then
      cross_edge_flag = .true.

    ! b) The closest edge is somewhere within the stencil. Check the maximum
    !    departure points to left and right to determine whether there is a
    !    calculation that will cross the panel edge
    else
      ! i) Check left side of stencil, left face
      if (consider_left) then
        max_int_dep_dist = -edge_dist_left + order + 1
        dep_idx = map_depk(dof_L) + (ndep - 1) / 2 + max_int_dep_dist
        if (dep_highest_k(dep_idx) > -1) cross_edge_flag = .true.

        ! ii) Check left side of stencil, right face
        if (.not. cross_edge_flag .and. num_faces_to_check == 2) then
          max_int_dep_dist = -edge_dist_left + order + 1
          dep_idx = map_depk(dof_R) + (ndep - 1) / 2 + max_int_dep_dist
          if (dep_highest_k(dep_idx) > -1) cross_edge_flag = .true.
        end if

      end if
      ! iii) Check right side of stencil, left face
      if (.not. cross_edge_flag .and. consider_right) then
        max_int_dep_dist = edge_dist_right - order - 1
        dep_idx = map_depk(dof_L) + (ndep - 1) / 2 + max_int_dep_dist
        if (dep_highest_k(dep_idx) > -1) cross_edge_flag = .true.

        ! iv) Check right side of stencil, right face
        if (.not. cross_edge_flag .and. num_faces_to_check == 2) then
          max_int_dep_dist = edge_dist_right - order - 1
          dep_idx = map_depk(dof_R) + (ndep - 1) / 2 + max_int_dep_dist
          if (dep_highest_k(dep_idx) > -1) cross_edge_flag = .true.
        end if
      end if
    end if
  end if

end function crosses_panel_edge


!> @brief Determines whether a calculation will cross the edge of a cubed sphere
!!        panel, by checking the distances to the panel edges. No consideration
!!        of the calculation is used other than the stencil size.
!> @param[in]  edge_dist_left       Distance to the left panel edge
!> @param[in]  edge_dist_right      Distance to the right panel edge
!> @param[in]  extent_size          Size of the stencil extent
!> @param[in]  panel_id             The ID of the panel for this column
!> @param[in]  direction            The direction this calculation is performed
!> @return                          Flag to indicate whether calculation will
!!                                  cross the panel edge
function crosses_rotated_panel_edge(edge_dist_left, edge_dist_right,           &
                                    extent_size, panel_id, direction)          &
                                    result(cross_edge_flag)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: edge_dist_left
  integer(kind=i_def), intent(in) :: edge_dist_right
  integer(kind=i_def), intent(in) :: extent_size
  integer(kind=i_def), intent(in) :: panel_id
  integer(kind=i_def), intent(in) :: direction

  ! Local variables
  integer(kind=i_def) :: face_dir
  logical(kind=l_def) :: cross_edge_flag, consider_left, consider_right

  cross_edge_flag = .false.

  consider_left = edge_dist_left > 0 .and. edge_dist_left < extent_size + 1
  consider_right = edge_dist_right > 0 .and. edge_dist_right < extent_size + 1

  consider_left = (                                                            &
    consider_left .and. ABS(rotated_panel_neighbour(panel_id, direction)) > 0  &
  )
  face_dir = 2_i_def + direction
  consider_right = (                                                           &
    consider_right .and. ABS(rotated_panel_neighbour(panel_id, face_dir)) > 0  &
  )

  cross_edge_flag = (consider_left .or. consider_right)

end function crosses_rotated_panel_edge


!> @brief Determines the neighbour of a panel in a given direction
!> @param[in]  panel_id    The ID of the current panel
!> @param[in]  direction   The direction of the edge to check
!> @return                 The ID of the neighbouring panel
function panel_neighbour(panel_id, direction) result(neighbour)

  implicit none

  integer(kind=i_def), intent(in) :: panel_id
  integer(kind=i_def), intent(in) :: direction

  integer(kind=i_def) :: neighbour
  integer(kind=i_def) :: panel_W, panel_E, panel_S, panel_N

  ! Encode each of the neighbouring panels to this panel
  select case (panel_id)
  case (1)
    panel_W = 4
    panel_E = 2
    panel_S = 6
    panel_N = 5
  case (2)
    panel_W = 1
    panel_E = 3
    panel_S = 6
    panel_N = 5
  case (3)
    panel_W = 6
    panel_E = 5
    panel_S = 4
    panel_N = 2
  case (4)
    panel_W = 6
    panel_E = 5
    panel_S = 1
    panel_N = 3
  case (5)
    panel_W = 4
    panel_E = 2
    panel_S = 1
    panel_N = 3
  case (6)
    panel_W = 1
    panel_E = 3
    panel_S = 4
    panel_N = 2
  end select

  select case (direction)
  case (W)
    neighbour = panel_W
  case (E)
    neighbour = panel_E
  case (S)
    neighbour = panel_S
  case (N)
    neighbour = panel_N
  end select

end function panel_neighbour


!> @details Determines whether crossing an edge of a cubed sphere panel edge
!!          involves a rotation of the the local directions on the panel.
!> @param[in] panel_id    The ID of the current panel
!> @param[in] direction   The direction of the edge to check
!> @return                Flag to indicate whether the panel edge is rotated
function rotated_panel_neighbour(panel_id, direction) result(rotation_flag)

  implicit none

  integer(kind=i_def), intent(in) :: panel_id
  integer(kind=i_def), intent(in) :: direction

  integer(kind=i_def) :: neighbour
  integer(kind=i_def) :: panel_edge
  integer(kind=i_def) :: rotation_flag

  neighbour = panel_neighbour(panel_id, direction)
  panel_edge = 10*panel_id + neighbour

  ! This is a list of panel edges will involve a rotation in directions
  select case (panel_edge)
  ! Clockwise rotations
  case (41, 32, 16, 25, 64, 53)
    rotation_flag = 1

  ! Anti-clockwise rotations
  case (14, 23, 61, 52, 46, 35)
    rotation_flag = -1

  ! No rotation
  case default
    rotation_flag = 0

  end select

end function rotated_panel_neighbour

end module panel_edge_support_mod
