!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Remap a scalar field (W3 or Wtheta) from the standard cubed sphere mesh
!!        to the extended cubed sphere
module remap_on_extended_mesh_kernel_mod

use kernel_mod,        only: kernel_type
use argument_mod,      only: arg_type, func_type,        &
                             GH_FIELD, GH_SCALAR,        &
                             GH_REAL, GH_INTEGER,        &
                             GH_LOGICAL,                 &
                             GH_READ, GH_WRITE,          &
                             ANY_DISCONTINUOUS_SPACE_1,  &
                             ANY_DISCONTINUOUS_SPACE_2,  &
                             ANY_DISCONTINUOUS_SPACE_3,  &
                             GH_BASIS, HALO_CELL_COLUMN, &
                             STENCIL, CROSS2D,           &
                             GH_EVALUATOR
use constants_mod,     only: r_def, r_tran, i_def, l_def, LARGE_REAL_POSITIVE
use fs_continuity_mod, only: Wchi

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: remap_on_extended_mesh_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                                           &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                   &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2),                   &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),                   &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3),                   &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                                               &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                                               &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                                               &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                                                &
       /)
  integer :: operates_on = HALO_CELL_COLUMN
contains
  procedure, nopass :: remap_on_extended_mesh_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: remap_on_extended_mesh_code
contains

!> @brief Remap a W3 or Wtheta scalar field from the standard cubed sphere mesh
!!        to the extended cubed sphere mesh.
!> @param[in]     nlayers          Number of layers
!> @param[in,out] remap_field      Field to compute remap values in
!> @param[in]     field            Field to map from (on original mesh)
!> @param[in]     stencil_size     Size of the x-stencil (number of cells)
!> @param[in]     stencil          Dofmaps for the x-stencil
!> @param[in]     max_length       Maximum stencil branch length
!> @param[in]     remap_weights    Interpolation weights for remapping
!> @param[in]     remap_indices    Interpolation indices for remapping
!> @param[in]     panel_id         Id of cubed sphere panel
!> @param[in]     ndata            Number of data points for a linear (=2) or cubic (=4) remapping
!> @param[in]     monotone         Flag to enforce a monotone interpolation
!> @param[in]     enforce_minvalue Flag to enforce a minimum value on the interpolation
!> @param[in]     minvalue         Minimum value for the interpolation
!> @param[in]     ndf_ws           Number of degrees of freedom per cell for scalar fields
!> @param[in]     undf_ws          Number of unique degrees of freedom for scalar fields
!> @param[in]     map_ws           Dofmap for the cell at the base of the column for scalar fields
!> @param[in]     ndf_wr           Number of degrees of freedom per cell for remap fields
!> @param[in]     undf_wr          Number of unique degrees of freedom for remap fields
!> @param[in]     map_wr           Dofmap for the cell at the base of the column for remap fields
!> @param[in]     ndf_pid          Number of degrees of freedom per cell for the panel id field
!> @param[in]     undf_pid         Number of unique degrees of freedom for the panel id field
!> @param[in]     map_pid          Dofmap for the cell at the base of the column for the panel id field
subroutine remap_on_extended_mesh_code(nlayers,                              &
                                       remap_field,                          &
                                       field,                                &
                                       stencil_size, stencil, max_length,    &
                                       remap_weights, remap_indices,         &
                                       panel_id,                             &
                                       ndata,                                &
                                       monotone, enforce_minvalue, minvalue, &
                                       ndf_ws, undf_ws, map_ws,              &
                                       ndf_wr, undf_wr, map_wr,              &
                                       ndf_pid, undf_pid, map_pid            &
                                      )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: nlayers
  integer(kind=i_def), dimension(4),       intent(in) :: stencil_size
  integer(kind=i_def),                     intent(in) :: max_length
  integer(kind=i_def),                     intent(in) :: ndf_ws, undf_ws
  integer(kind=i_def),                     intent(in) :: ndf_wr, undf_wr
  integer(kind=i_def),                     intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf_ws),  intent(in) :: map_ws
  integer(kind=i_def), dimension(ndf_wr),  intent(in) :: map_wr
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  integer(kind=i_def), intent(in) :: ndata
  logical(kind=l_def), intent(in) :: monotone
  logical(kind=l_def), intent(in) :: enforce_minvalue
  real(kind=r_tran),   intent(in) :: minvalue

  integer(kind=i_def), dimension(ndf_ws,  max_length,    4), intent(in) :: stencil

  real(kind=r_tran),   dimension(undf_ws),  intent(inout) :: remap_field
  real(kind=r_tran),   dimension(undf_ws),  intent(in)    :: field
  real(kind=r_tran),   dimension(undf_wr),  intent(in)    :: remap_weights
  integer(kind=i_def), dimension(undf_wr),  intent(in)    :: remap_indices
  real(kind=r_def),    dimension(undf_pid), intent(in)    :: panel_id

  integer(kind=i_def)            :: owned_panel, halo_panel, &
                                    panel_edge, k, n, nl
  integer(kind=i_def)            :: interp_dir, dir_id
  integer(kind=i_def), parameter :: interp_dir_alpha = 1
  integer(kind=i_def), parameter :: interp_dir_beta  = 2
  integer(kind=i_def)            :: ncells_in_stencil
  integer(kind=i_def)            :: alpha_dir_id, beta_dir_id

  integer(kind=i_def), dimension(2*max_length-1) :: stencil_1d

  real(kind=r_tran), dimension(0:3) :: field_in_stencil
  real(kind=r_tran)                 :: min_value_in_stencil
  real(kind=r_tran)                 :: max_value_in_stencil

  ! Assume the first entry in panel id corresponds to an owned (not halo) cell
  owned_panel = int(panel_id(1),i_def)
  ! Panel id for this column
  halo_panel = int(panel_id(map_pid(1)),i_def)

  ! nl = nlayers for Wtheta fields (ndf_ws=2)
  ! nl = nlayers-1 for W3 fields (ndf_ws=1)
  nl = (nlayers - 1) + (ndf_ws - 1)

  ! Only need to remap if the halo is on a different panel to the
  ! owned cell (otherwise the remaped field is identical to the original
  ! field and should have been picked up by a copy_field call)
  if ( halo_panel /= owned_panel ) then

    ! Compute interpolation direction (alpha or beta)
    ! the first digit of panel_edge is the panel we are interpolating to and the
    ! second is the one we are interpolating from, i.e. panel_edge = 36 means we
    ! are interpolating from panel 6 to panel 3
    panel_edge = 10*owned_panel + halo_panel
    select case (panel_edge)
      case (15, 16, 25, 26, 32, 34, 41, 43, 51, 53, 62, 64)
        interp_dir = interp_dir_alpha
      case (12, 14, 21, 23, 35, 36, 45, 46, 52, 54, 61, 63)
        interp_dir = interp_dir_beta
    end select

    ! The stencils are ordered (W,S,E,N) on their local panel with (W,E) in the
    ! alpha direction and (S,N) in the beta direction. However, this doesn't
    ! necessarily correspond to the same direction on the owned panel so we
    ! need to potentially rotate the directions.
    ! i.e on the halo of panel 1 corresponding to panel 4 we want to
    ! interpolate in the beta direction of panel 1 but this corresponds to
    ! alpha direction on panel 4 so we set alpha_dir_id = 2 here
    select case( panel_edge )
      case (14, 16, 23, 25, 32, 35, 41, 46, 52, 53, 61, 64)
        ! Change in orientation of alpha and beta directions
        alpha_dir_id = 2
        beta_dir_id  = 1
      case default
        ! Default values for when there is no change in orientation
        alpha_dir_id = 1
        beta_dir_id  = 2
    end select

    if ( interp_dir == interp_dir_alpha ) then
      dir_id = alpha_dir_id
    else
      dir_id = beta_dir_id
    end if

    ! Compute combined 1D stencils
    ! for a Cross stencil of the form:
    !      | 7 |
    !      | 5 |
    !  | 2 | 1 | 4 | 6 |
    !      | 3 |
    ! Stored as stencil =[
    ! 1, 2;
    ! 1, 3;
    ! 1, 4, 6;
    ! 1, 5, 7]
    ! Then the new stencil_1d = [
    ! 1, 2, 4, 6;
    ! 1, 3, 5, 7]

    ncells_in_stencil = stencil_size(dir_id) + stencil_size(dir_id+2) - 1

    stencil_1d = 0
    do n = 1, stencil_size(dir_id)
      stencil_1d(n) = stencil(1,n,dir_id)
    end do
    do n = 1, stencil_size(dir_id+2)-1
      stencil_1d(n+stencil_size(dir_id)) = stencil(1,n+1,dir_id+2)
    end do

    if ( ndata == 2 ) then
      do k = 0, nl
        remap_field(map_ws(1)+k) = remap_weights(map_wr(1)+0)*field(stencil_1d(remap_indices(map_wr(1)+0))+k) &
                                 + remap_weights(map_wr(1)+1)*field(stencil_1d(remap_indices(map_wr(1)+1))+k)

      end do
    else
      if ( monotone ) then
        do k = 0, nl
          do n = 0, 3
            field_in_stencil(n) = field(stencil_1d(remap_indices(map_wr(1)+n))+k)
          end do
          min_value_in_stencil = minval( field_in_stencil )
          max_value_in_stencil = maxval( field_in_stencil )

          remap_field(map_ws(1)+k) = remap_weights(map_wr(1)+0)*field_in_stencil(0) &
                                   + remap_weights(map_wr(1)+1)*field_in_stencil(1) &
                                   + remap_weights(map_wr(1)+2)*field_in_stencil(2) &
                                   + remap_weights(map_wr(1)+3)*field_in_stencil(3)
          remap_field(map_ws(1)+k) = max( min_value_in_stencil,      &
                                          min( max_value_in_stencil, &
                                               remap_field(map_ws(1)+k) ) )

        end do

      else
        do k = 0, nl
          remap_field(map_ws(1)+k) = remap_weights(map_wr(1)+0)*field(stencil_1d(remap_indices(map_wr(1)+0))+k) &
                                   + remap_weights(map_wr(1)+1)*field(stencil_1d(remap_indices(map_wr(1)+1))+k) &
                                   + remap_weights(map_wr(1)+2)*field(stencil_1d(remap_indices(map_wr(1)+2))+k) &
                                   + remap_weights(map_wr(1)+3)*field(stencil_1d(remap_indices(map_wr(1)+3))+k)
        end do
      end if
    end if
    if ( enforce_minvalue ) then
      do k = 0, nl
        remap_field(map_ws(1)+k) = max(remap_field(map_ws(1)+k), minvalue)
      end do
    end if
  else
    ! Halo value is on same panel as owned value so just copy the field
    do k = 0, nl
      remap_field(map_ws(1)+k) = field(map_ws(1)+k)
    end do

  end if

end subroutine remap_on_extended_mesh_code

end module remap_on_extended_mesh_kernel_mod
