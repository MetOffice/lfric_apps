!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the distance of columns from the edges of cubed sphere panels
module panel_edge_dist_kernel_mod

  use argument_mod,          only : arg_type, func_type,         &
                                    GH_FIELD, GH_REAL, GH_WRITE, &
                                    GH_READ, GH_INTEGER,         &
                                    ANY_DISCONTINUOUS_SPACE_3,   &
                                    STENCIL, CROSS2D,            &
                                    CELL_COLUMN

  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W3
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : W, S, N, E

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: panel_edge_dist_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                        &
         arg_type(GH_FIELD*4, GH_INTEGER, GH_WRITE, W3),                       &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3, &
                                                            STENCIL(CROSS2D))  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: panel_edge_dist_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: panel_edge_dist_code

contains

!> @brief Computes the distance of columns from the edges of cubed sphere panels
!> @param[in]     nlayers            Number of layers in mesh
!> @param[in,out] panel_edge_dist_W  2D field containing the distance of this
!!                                   column from the panel edge to the West
!> @param[in,out] panel_edge_dist_E  2D field containing the distance of this
!!                                   column from the panel edge to the East
!> @param[in,out] panel_edge_dist_S  2D field containing the distance of this
!!                                   column from the panel edge to the South
!> @param[in,out] panel_edge_dist_N  2D field containing the distance of this
!!                                   column from the panel edge to the North
!> @param[in]     panel_id           Field giving the ID for mesh panels.
!> @param[in]     stencil_sizes      Size of each branch of the cross stencil
!> @param[in]     stencil_map        Map for W3 cross stencil
!> @param[in]     max_stencil_size   Maximum size of a cross stencil branch
!> @param[in]     ndf_w3             Num DoFs per cell for W3
!> @param[in]     undf_w3            Num DoFs for this partition for W3
!> @param[in]     map_w3             DoFmap for W3
!> @param[in]     ndf_pid            Num DoFs per cell for panel ID
!> @param[in]     undf_pid           Num DoFs for this partition for panel ID
!> @param[in]     map_pid            DoFmap for panel ID
subroutine panel_edge_dist_code( nlayers,                                      &
                                 panel_edge_dist_W,                            &
                                 panel_edge_dist_E,                            &
                                 panel_edge_dist_S,                            &
                                 panel_edge_dist_N,                            &
                                 panel_id,                                     &
                                 stencil_sizes,                                &
                                 max_stencil_size,                             &
                                 stencil_map,                                  &
                                 ndf_w3, undf_w3, map_w3,                      &
                                 ndf_pid, undf_pid, map_pid                    &
                               )

  use panel_edge_support_mod, only: FAR_AWAY

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w3
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: ndf_pid
  integer(kind=i_def), intent(in)    :: undf_pid
  integer(kind=i_def), intent(in)    :: max_stencil_size
  integer(kind=i_def), intent(in)    :: stencil_sizes(4)

  integer(kind=i_def), intent(inout) :: panel_edge_dist_W(undf_w3)
  integer(kind=i_def), intent(inout) :: panel_edge_dist_E(undf_w3)
  integer(kind=i_def), intent(inout) :: panel_edge_dist_N(undf_w3)
  integer(kind=i_def), intent(inout) :: panel_edge_dist_S(undf_w3)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: stencil_map(ndf_w3,max_stencil_size,4)

  ! Local variables
  integer(kind=i_def) :: j
  integer(kind=i_def) :: ipanel, ipanel_W, ipanel_S, ipanel_E, ipanel_N

  ! Get panel IDs at the end of each part of the cross stencil
  ipanel = int(panel_id(map_pid(1)), i_def)
  ipanel_W = int(panel_id(stencil_map(1, stencil_sizes(W), W)), i_def)
  ipanel_S = int(panel_id(stencil_map(1, stencil_sizes(S), S)), i_def)
  ipanel_E = int(panel_id(stencil_map(1, stencil_sizes(E), E)), i_def)
  ipanel_N = int(panel_id(stencil_map(1, stencil_sizes(N), N)), i_def)

  ! Initialise dists to large negative number
  panel_edge_dist_W(map_w3(1)) = -FAR_AWAY
  panel_edge_dist_E(map_w3(1)) = -FAR_AWAY
  panel_edge_dist_N(map_w3(1)) = -FAR_AWAY
  panel_edge_dist_S(map_w3(1)) = -FAR_AWAY

  ! If the panel ID changes along the the stencil, we are near a panel edge
  ! Indicate in which direction the panel is changing
  ! NB: this won't work for incredibly coarse meshes (e.g. C1 or C2)
  ! Change of panel to W -------------------------------------------------------
  if (ipanel /= ipanel_W) then
    ! Determine how close this cell is to the panel edge
    do j = 2, stencil_sizes(W)
      if (INT(panel_id(stencil_map(1, j, W)), i_def) == ipanel_W) then
        panel_edge_dist_W(map_w3(1)) = j - 1
        EXIT
      end if
    end do
  end if

  ! Change of panel to E -------------------------------------------------------
  if (ipanel /= ipanel_E) then
    ! Determine how close this cell is to the panel edge
    do j = 2, stencil_sizes(E)
      if (INT(panel_id(stencil_map(1, j, E)), i_def) == ipanel_E) then
        panel_edge_dist_E(map_w3(1)) = j - 1
        EXIT
      end if
    end do
  end if

  ! Change of panel to S -------------------------------------------------------
  if (ipanel /= ipanel_S) then
    ! Determine how close this cell is to the panel edge
    do j = 2, stencil_sizes(S)
      if (INT(panel_id(stencil_map(1, j, S)), i_def) == ipanel_S) then
        panel_edge_dist_S(map_w3(1)) = j - 1
        EXIT
      end if
    end do
  end if

  ! Change of panel to N -------------------------------------------------------
  if (ipanel /= ipanel_N) then
    ! Determine how close this cell is to the panel edge
    do j = 2, stencil_sizes(N)
      if (INT(panel_id(stencil_map(1, j, N)), i_def) == ipanel_N) then
        panel_edge_dist_N(map_w3(1)) = j - 1
        EXIT
      end if
    end do
  end if

end subroutine panel_edge_dist_code

end module panel_edge_dist_kernel_mod
