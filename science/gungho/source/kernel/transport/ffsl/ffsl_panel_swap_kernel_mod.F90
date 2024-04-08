!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Swaps the values of FFSL fields from x- and y-sweeps by rotated panels
!> @details FFSL uses 1D sweeps to calculate horizontal mass fluxes. In the
!!          outer stage, a y-sweep is performed on a field that has undergone a
!!          sweep in the x-direction, and vice versa.
!!          However when this sweep involves a stencil that crosses a cubed
!!          sphere panel edge, the x- and y-directions may be rotated. It is
!!          therefore appropriate to swap the values in the halo regions of
!!          these fields.
!!          This kernel performs that operation.

module ffsl_panel_swap_kernel_mod

  use argument_mod,       only : arg_type,             &
                                 GH_FIELD, GH_REAL,    &
                                 CELL_COLUMN, GH_READ, &
                                 GH_READWRITE,         &
                                 ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,      only : r_tran, i_def, r_def
  use fs_continuity_mod,  only : W3
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: ffsl_panel_swap_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                            &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                        & ! field_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                        & ! field_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3 ) & ! panel_id
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_panel_swap_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: ffsl_panel_swap_code

contains

  !> @brief Swaps the values of FFSL fields from x- and y-sweeps, to account for
  !!        the rotation of cubed sphere panel directions
  !> @param[in]     nlayers       Number of layers
  !> @param[in,out] field_x       Transported field from a step in x
  !> @param[in,out] field_y       Transported field from a step in y
  !> @param[in]     panel_id      Field containing identifiers of mesh panels
  !> @param[in]     ndf_w3        Num of DoFs for W3 per cell
  !> @param[in]     undf_w3       Num of DoFs for this partition for W3
  !> @param[in]     map_w3        Map for W3
  !> @param[in]     ndf_wp        Num of DoFs for panel ID
  !> @param[in]     undf_wp       Num of DoFs for this partition for panel ID
  !> @param[in]     map_wp        Map for panel ID function space
  subroutine ffsl_panel_swap_code( nlayers,        &
                                   field_x,        &
                                   field_y,        &
                                   panel_id,       &
                                   ndf_w3,         &
                                   undf_w3,        &
                                   map_w3,         &
                                   ndf_wp,         &
                                   undf_wp,        &
                                   map_wp )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_wp
    integer(kind=i_def), intent(in) :: ndf_wp

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_wp),  intent(in) :: map_wp

    ! Arguments: Fields
    real(kind=r_tran),   dimension(undf_w3),  intent(inout) :: field_x
    real(kind=r_tran),   dimension(undf_w3),  intent(inout) :: field_y
    real(kind=r_def),   dimension(undf_wp),   intent(in)    :: panel_id

    ! Indices
    real(kind=r_tran)   :: tmp_x, tmp_y
    integer(kind=i_def) :: k, cell_panel, owned_panel, panel_boundary, col_idx

    owned_panel = INT(panel_id(1), i_def)
    cell_panel = INT(panel_id(map_wp(1)), i_def)
    panel_boundary = 10*owned_panel + cell_panel

    ! Check if this cell has crossed a panel-rotation boundary
    select case (panel_boundary)
    case (14, 41, 25, 52, 23, 32, 16, 61, 46, 64, 35, 53)

      col_idx = map_w3(1)

      ! Loop up the column and swap field_x and field_y values
      do k = 0, nlayers - 1
        tmp_x = field_x(col_idx + k)
        tmp_y = field_y(col_idx + k)
        field_x(col_idx + k) = tmp_y
        field_y(col_idx + k) = tmp_x
      end do

    ! Do nothing otherwise
    end select

  end subroutine ffsl_panel_swap_code

end module ffsl_panel_swap_kernel_mod
