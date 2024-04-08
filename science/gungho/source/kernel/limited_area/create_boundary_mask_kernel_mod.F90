!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes a mask with 1s on the solver boundary and zero elsewhere
!> @details Makes use of the onion_layers to determine location
!>          if dofs are between the specified distance from the edge and
!>          that distance+1.
!>          Requires stencil with depth 1.
module create_boundary_mask_kernel_mod

  use argument_mod,              only : arg_type,            &
                                        GH_SCALAR, GH_FIELD, &
                                        GH_READ, GH_WRITE,   &
                                        GH_REAL,             &
                                        CELL_COLUMN,         &
                                        STENCIL, CROSS
  use constants_mod,             only : r_def, i_def
  use kernel_mod,                only : kernel_type
  use fs_continuity_mod,         only : W3, W2
  use reference_element_mod,     only : N,S,E,W

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: create_boundary_mask_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                              &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W2),                & ! mask
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3, STENCIL(CROSS)), & ! onion_layers
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                      & ! boundary_inner_layer
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: create_boundary_mask_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: create_boundary_mask_code

contains

!> @param[in] nlayers              Number of layers
!> @param[in,out] mask             The output mask
!> @param[in] onion_layers         The mask defining 'onion layers'
!> @param[in] stencil_size         Local stencil size (will be reduced near a boundary)
!> @param[in] stencil_map          Stencil map
!> @param[in] boundary_inner_layer Onion layer adjacent to the solver boundary
!> @param[in] ndf                  Number of degrees of freedom per cell
!>                                 for W2 mask
!> @param[in] undf                 Total number of degrees of freedom
!>                                 for W2 mask
!> @param[in] map                  Dofmap for the cell at the base of the column
!>                                 for W2 mask
!> @param[in] ndf_w3               Number of degrees of freedom per cell
!>                                 for onion_layers
!> @param[in] undf_w3              Total number of degrees of freedom
!>                                 for onion_layers
!> @param[in] map_w3               Dofmap for the cell at the base of the column
!>                                 for onion_layers
subroutine create_boundary_mask_code( nlayers,              &
                                      mask,                 &
                                      onion_layers,         &
                                      stencil_size,         &
                                      stencil_map,          &
                                      boundary_inner_layer, &
                                      ndf,                  &
                                      undf,                 &
                                      map,                  &
                                      ndf_w3,               &
                                      undf_w3,              &
                                      map_w3)

  implicit none

  ! Arguments
  integer(kind=i_def),                  intent(in) :: nlayers,     &
                                                      ndf,         &
                                                      undf,        &
                                                      ndf_w3,      &
                                                      undf_w3
  real(kind=r_def), dimension(undf),    intent(inout) :: mask
  real(kind=r_def), dimension(undf_w3), intent(in)    :: onion_layers
  integer(kind=i_def), intent(in)                     :: stencil_size
  real(kind=r_def), intent(in)                        :: boundary_inner_layer
  integer(kind=i_def), dimension(ndf_w3,stencil_size), intent(in) :: stencil_map
  integer(kind=i_def), dimension(ndf),    intent(in)  :: map
  integer(kind=i_def), dimension(ndf_w3), intent(in)  :: map_w3

  ! Internal variables
  integer(kind=i_def) :: k, dir
  integer(kind=i_def) :: onion_layer, cell_next_layer

  onion_layer = int(onion_layers(stencil_map(1,1)))

  ! The solver boundary must be in the blending zone, so only operate
  ! where onion_layers are positive for cost saving
  ! Don't operate at the edge where cells do not have all 4 directions
  if (onion_layer > 0 .and. stencil_size > 4)then
    !loop over the directions of the stencil that coincide with the horizontal W2 dofs
    !and identify the horizontal W2 DoFs common to those cells that straddle the
    !solver boundary.
    do dir=1,4
      cell_next_layer = int(onion_layers(stencil_map(1,dir+1)))
      if (cell_next_layer  == boundary_inner_layer + 1 &
         .and. onion_layer == boundary_inner_layer )then
        do k=0,nlayers-1
          mask(map(dir)+k) = 1.0_r_def
        end do
      end if
    end do
  end if


end subroutine create_boundary_mask_code

end module create_boundary_mask_kernel_mod
