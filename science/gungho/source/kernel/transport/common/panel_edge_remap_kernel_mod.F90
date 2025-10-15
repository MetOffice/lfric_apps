!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Remap a W3 or Wtheta scalar field at the edges of cubed-sphere panels
!!        to the corresponding neighbouring panels.
!!        Field values that are further away from the panel edge than the
!!        specified depth are set to zero.
module panel_edge_remap_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_INTEGER, GH_LOGICAL,              &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_1,                    &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 ANY_DISCONTINUOUS_SPACE_7,                    &
                                 OWNED_AND_HALO_CELL_COLUMN,                   &
                                 STENCIL, CROSS2D
use constants_mod,         only: r_def, r_tran, i_def, l_def
use reference_element_mod, only: W, S, N, E

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: panel_edge_remap_kernel_type
  private
  type(arg_type) :: meta_args(17) = (/                                         &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_7,   &
                                                            STENCIL(CROSS2D)), &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_7,   &
                                                            STENCIL(CROSS2D)), &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD*4, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_FIELD,   GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                              &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                              &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                              &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                              &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),                              &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                               &
  /)
  integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
contains
  procedure, nopass :: panel_edge_remap_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: panel_edge_remap_code

contains

!> @brief Remap a W3 or Wtheta scalar field at the edges of cubed-sphere panels
!!        to the corresponding neighbouring panels.
!> @param[in]     nlayers               Number of layers
!> @param[in,out] remapped_in_x         Remapped field in local x direction
!> @param[in,out] remapped_in_y         Remapped field in local y direction
!> @param[in]     field_for_x           Field to be remapped in x direction
!> @param[in]     wsx_stencil_size      Size of field stencil (num cells)
!> @param[in]     wsx_stencil_map       DoF map for the field stencil
!> @param[in]     wsx_max_length        Maximum stencil branch length
!> @param[in]     field_for_t           Field to be remapped in t direction
!> @param[in]     wsy_stencil_size      Size of field stencil (num cells)
!> @param[in]     wsy_stencil_map       DoF map for the field stencil
!> @param[in]     wsy_max_length        Maximum stencil branch length
!> @param[in]     panel_edge_weights_x  Remapping interpolation weights for x
!> @param[in]     panel_edge_weights_y  Remapping interpolation weights for y
!> @param[in]     panel_edge_indices_x  Remapping interpolation indices for x
!> @param[in]     panel_edge_indices_y  Remapping interpolation indices for y
!> @param[in]     panel_edge_dist_W     2D field containing the distance of each
!!                                      column from the panel edge to the West
!> @param[in]     panel_edge_dist_E     2D field containing the distance of each
!!                                      column from the panel edge to the East
!> @param[in]     panel_edge_dist_S     2D field containing the distance of each
!!                                      column from the panel edge to the South
!> @param[in]     panel_edge_dist_N     2D field containing the distance of each
!!                                      column from the panel edge to the North
!> @param[in]     halo_mask_x           Mask for halo cells in the x direction
!> @param[in]     halo_mask_y           Mask for halo cells in the y direction
!> @param[in]     compute_all_halo      Flag indicating whether to perform
!!                                      remapping in both directions in halo
!!                                      cells, or in the direction of the halo
!> @param[in]     depth                 Maximum halo depth to consider
!> @param[in]     ndata                 Number of data points in the remapping
!> @param[in]     monotone              Flag to enforce a monotone interpolation
!> @param[in]     enforce_minvalue      Flag to enforce a bounded interpolation
!> @param[in]     minvalue              Minimum value for the interpolation
!> @param[in]     ndf_wr                Num DoFs per cell for remapped fields
!> @param[in]     undf_wr               Num DoFs for this partition for remapped
!!                                      fields
!> @param[in]     map_wr                DoF map for remapped fields
!> @param[in]     ndf_ws                Num DoFs per cell for input scalar field
!> @param[in]     undf_ws               Num DoFs for this partition for input
!!                                      scalar field
!> @param[in]     map_ws                DoF map for input scalar field
!> @param[in]     ndf_ww                Num DoFs per cell for remapping weights
!> @param[in]     undf_ww               Num DoFs for this partition for weights
!> @param[in]     map_ww                DoF map for remapping weights
!> @param[in]     ndf_pid               Num DoFs per cell for panel ID field
!> @param[in]     undf_pid              Num DoFs for this partition for panel ID
!> @param[in]     map_pid               DoF map for panel ID field
subroutine panel_edge_remap_code( nlayers,                                     &
                                  remapped_in_x,                               &
                                  remapped_in_y,                               &
                                  field_for_x,                                 &
                                  wsx_stencil_size,                            &
                                  wsx_max_length,                              &
                                  wsx_stencil_map,                             &
                                  field_for_y,                                 &
                                  wsy_stencil_size,                            &
                                  wsy_max_length,                              &
                                  wsy_stencil_map,                             &
                                  panel_edge_weights_x,                        &
                                  panel_edge_weights_y,                        &
                                  panel_edge_indices_x,                        &
                                  panel_edge_indices_y,                        &
                                  panel_edge_dist_W,                           &
                                  panel_edge_dist_E,                           &
                                  panel_edge_dist_S,                           &
                                  panel_edge_dist_N,                           &
                                  halo_mask_x,                                 &
                                  halo_mask_y,                                 &
                                  compute_all_halo,                            &
                                  depth,                                       &
                                  ndata,                                       &
                                  monotone, enforce_minvalue, minvalue,        &
                                  ndf_wr, undf_wr, map_wr,                     &
                                  ndf_ws, undf_ws, map_ws,                     &
                                  ndf_ww, undf_ww, map_ww,                     &
                                  ndf_pid, undf_pid, map_pid                   &
  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: wsx_stencil_size(4), wsy_stencil_size(4)
  integer(kind=i_def), intent(in)    :: wsx_max_length, wsy_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_ws, undf_ws
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: ndata
  integer(kind=i_def), intent(in)    :: depth
  logical(kind=l_def), intent(in)    :: monotone
  logical(kind=l_def), intent(in)    :: enforce_minvalue
  logical(kind=l_def), intent(in)    :: compute_all_halo
  real(kind=r_tran),   intent(in)    :: minvalue

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_ws(ndf_ws)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: wsx_stencil_map(ndf_ws, wsx_max_length, 4)
  integer(kind=i_def), intent(in)    :: wsy_stencil_map(ndf_ws, wsy_max_length, 4)

  real(kind=r_tran),   intent(inout) :: remapped_in_x(undf_wr)
  real(kind=r_tran),   intent(inout) :: remapped_in_y(undf_wr)
  real(kind=r_tran),   intent(in)    :: field_for_x(undf_ws)
  real(kind=r_tran),   intent(in)    :: field_for_y(undf_ws)
  real(kind=r_tran),   intent(in)    :: panel_edge_weights_x(undf_ww)
  real(kind=r_tran),   intent(in)    :: panel_edge_weights_y(undf_ww)
  integer(kind=i_def), intent(in)    :: panel_edge_indices_x(undf_ww)
  integer(kind=i_def), intent(in)    :: panel_edge_indices_y(undf_ww)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_W(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_E(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_S(undf_pid)
  integer(kind=i_def), intent(in)    :: panel_edge_dist_N(undf_pid)
  integer(kind=i_def), intent(in)    :: halo_mask_x(undf_pid)
  integer(kind=i_def), intent(in)    :: halo_mask_y(undf_pid)

  ! Internal variables
  integer(kind=i_def) :: r_idx, pid_idx, nvert
  logical(kind=l_def) :: compute_x, compute_y

  ! Number of vertical levels to loop over, to handle both W3 and Wtheta
  ! nvert = nlayers for Wtheta fields (ndf_ws=2)
  ! nvert = nlayers-1 for W3 fields (ndf_ws=1)
  nvert = (nlayers - 1) + (ndf_ws - 1)

  r_idx = map_wr(1)
  pid_idx = map_pid(1)

  ! Determine whether to perform remapping --------------------------------------
  compute_x = (                                                                &
    (compute_all_halo .or. halo_mask_y(pid_idx) == 0) .and. (                  &
      ABS(panel_edge_dist_W(pid_idx)) <= depth .or.                            &
      ABS(panel_edge_dist_E(pid_idx)) <= depth                                 &
    )                                                                          &
  )

  compute_y = (                                                                &
    (compute_all_halo .or. halo_mask_x(pid_idx) == 0) .and. (                  &
      ABS(panel_edge_dist_N(pid_idx)) <= depth .or.                            &
      ABS(panel_edge_dist_S(pid_idx)) <= depth                                 &
    )                                                                          &
  )

  ! Remap at x edge of panel ---------------------------------------------------
  if (compute_x) then
    call panel_edge_remap_1d(                                                  &
            nlayers,                                                           &
            remapped_in_x,                                                     &
            field_for_x, wsx_max_length,                                       &
            ! Indices S and N, because crossing the x edge of a panel, and
            ! performing the remapping parallel to the edge
            wsx_stencil_size(S), wsx_stencil_map(:,:,S),                       &
            wsx_stencil_size(N), wsx_stencil_map(:,:,N),                       &
            panel_edge_weights_x, panel_edge_indices_x,                        &
            ndata, monotone, enforce_minvalue, minvalue,                       &
            ndf_wr, undf_wr, map_wr,                                           &
            ndf_ws, undf_ws, map_ws,                                           &
            ndf_ww, undf_ww, map_ww                                            &
    )
  else
    ! These field values should not ever be used properly, but the data may
    ! be accessed so should be filled with some value. Set it to a bad value
    ! so it can be detected if it is used
    remapped_in_x(r_idx : r_idx+nvert) = -1900.0_r_tran
  end if

  ! Remap at y edge of panel ---------------------------------------------------
  if (compute_y) then
    call panel_edge_remap_1d(                                                  &
            nlayers,                                                           &
            remapped_in_y,                                                     &
            field_for_y, wsy_max_length,                                       &
            ! Indices W and E, because crossing the y edge of a panel, and
            ! performing the remapping parallel to the edge
            wsy_stencil_size(W), wsy_stencil_map(:,:,W),                       &
            wsy_stencil_size(E), wsy_stencil_map(:,:,E),                       &
            panel_edge_weights_y, panel_edge_indices_y,                        &
            ndata, monotone, enforce_minvalue, minvalue,                       &
            ndf_wr, undf_wr, map_wr,                                           &
            ndf_ws, undf_ws, map_ws,                                           &
            ndf_ww, undf_ww, map_ww                                            &
    )
  else
    ! These field values should not ever be used properly, but the data may
    ! be accessed so should be filled with some value. Set it to a bad value
    ! so it can be detected if it is used
    remapped_in_y(r_idx : r_idx+nvert) = -1900.0_r_tran
  end if

end subroutine panel_edge_remap_code

!> @brief Private routine for remapping a field over a cubed-sphere panel edge
!!        for *a single direction* (x or y)
subroutine panel_edge_remap_1d( nlayers,                                       &
                                remapped_field,                                &
                                field, ws_max_length,                          &
                                ws_stencil_size_l, ws_stencil_l,               &
                                ws_stencil_size_r, ws_stencil_r,               &
                                panel_edge_weights, panel_edge_indices,        &
                                ndata, monotone, enforce_minvalue, minvalue,   &
                                ndf_wr, undf_wr, map_wr,                       &
                                ndf_ws, undf_ws, map_ws,                       &
                                ndf_ww, undf_ww, map_ww                        &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ws_stencil_size_l
  integer(kind=i_def), intent(in)    :: ws_stencil_size_r
  integer(kind=i_def), intent(in)    :: ws_max_length
  integer(kind=i_def), intent(in)    :: ndf_ww, undf_ww
  integer(kind=i_def), intent(in)    :: ndf_wr, undf_wr
  integer(kind=i_def), intent(in)    :: ndf_ws, undf_ws
  integer(kind=i_def), intent(in)    :: ndata
  logical(kind=l_def), intent(in)    :: monotone
  logical(kind=l_def), intent(in)    :: enforce_minvalue
  real(kind=r_tran),   intent(in)    :: minvalue

  integer(kind=i_def), intent(in)    :: map_wr(ndf_wr)
  integer(kind=i_def), intent(in)    :: map_ww(ndf_ww)
  integer(kind=i_def), intent(in)    :: map_ws(ndf_ws)
  integer(kind=i_def), intent(in)    :: ws_stencil_l(ndf_ws, ws_stencil_size_l)
  integer(kind=i_def), intent(in)    :: ws_stencil_r(ndf_ws, ws_stencil_size_r)

  real(kind=r_tran),   intent(inout) :: remapped_field(undf_wr)
  real(kind=r_tran),   intent(in)    :: field(undf_ws)
  real(kind=r_tran),   intent(in)    :: panel_edge_weights(undf_ww)
  integer(kind=i_def), intent(in)    :: panel_edge_indices(undf_ww)

  ! Internal variables
  integer(kind=i_def) :: n
  integer(kind=i_def) :: w_idx, r_idx, f_idx
  integer(kind=i_def) :: nvert
  integer(kind=i_def) :: ncells_in_stencil
  integer(kind=i_def) :: stencil_1d(2*ws_max_length-1)

  ! Stencils will have size of "nvert" variable declared below
  real(kind=r_tran)   :: min_stencil((nlayers - 1) + (ndf_ws - 1) + 1)
  real(kind=r_tran)   :: max_stencil((nlayers - 1) + (ndf_ws - 1) + 1)

  ! Number of vertical levels to loop over, to handle both W3 and Wtheta
  ! nvert = nlayers for Wtheta fields (ndf_ws=2)
  ! nvert = nlayers-1 for W3 fields (ndf_ws=1)
  nvert = (nlayers - 1) + (ndf_ws - 1)

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

  ncells_in_stencil = ws_stencil_size_l + ws_stencil_size_r - 1

  stencil_1d(:) = 0
  do n = 1, ws_stencil_size_l
    stencil_1d(n) = ws_stencil_l(1,n)
  end do
  ! Omit 1 from the second part, as the central cell is already included
  do n = 1, ws_stencil_size_r - 1
    stencil_1d(n+ws_stencil_size_l) = ws_stencil_r(1,n+1)
  end do

  if (monotone) then
    ! Initialise to extreme numbers
    min_stencil(:) = HUGE(0.0_r_tran)
    max_stencil(:) = -HUGE(0.0_r_tran)
  end if

  ! -------------------------------------------------------------------------- !
  ! Remap
  ! -------------------------------------------------------------------------- !
  ! Indices at base of column
  w_idx = map_ww(1)
  r_idx = map_wr(1)

  ! Set the output field for this column to be zero
  remapped_field(r_idx : r_idx+nvert) = 0.0_r_tran

  ! Loop through (multidata) interpolation points
  do n = 0, ndata - 1
    ! Base index of field to be remapped
    f_idx = stencil_1d(panel_edge_indices(w_idx+n))

    ! Increment the remapped field by the weighted field value
    remapped_field(r_idx : r_idx+nvert) = remapped_field(r_idx : r_idx+nvert)  &
      + panel_edge_weights(w_idx+n) * field(f_idx : f_idx+nvert)

    ! If applying monotonicity, store min/max values in an array
    if (monotone) then
      min_stencil(:) = MIN(min_stencil(:), field(f_idx : f_idx+nvert))
      max_stencil(:) = MAX(max_stencil(:), field(f_idx : f_idx+nvert))
    end if
  end do

  ! If specified, apply monotonicity by looping through column for final time
  if (monotone) then
    remapped_field(r_idx : r_idx+nvert) =                                      &
        MIN(max_stencil(:),                                                    &
        MAX(min_stencil(:), remapped_field(r_idx : r_idx+nvert))               &
    )

  ! No need to enforce minimum value if already enforcing monotonicity
  else if (enforce_minvalue) then
    remapped_field(r_idx : r_idx+nvert) =                                      &
        MAX(remapped_field(r_idx : r_idx+nvert), minvalue)
  end if

end subroutine panel_edge_remap_1d

end module panel_edge_remap_kernel_mod
