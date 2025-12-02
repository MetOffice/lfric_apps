!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Routine to calculate surface albedo (grid-box mean)

module calc_albedo_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_REAL,         &
                              GH_READ, GH_WRITE,         &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

public :: calc_albedo_kernel_type
public :: calc_albedo_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: calc_albedo_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                     &
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! surf_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf    
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2)  & ! sw_up_tile    
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: calc_albedo_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                Number of layers
!> @param[in,out] surf_albedo            Surface albedo (gridbox)
!> @param[in]     sw_down_surf           Total down shortwave flux (gridbox)
!> @param[in]     tile_fraction          Surface tile fractions
!> @param[in]     sw_up_tile             Total upwelling sw flux (tile)
!> @param[in]     ndf_2d                 DOFs per cell for 2-D / surface field
!> @param[in]     undf_2d                Total DOFs for for 2-D / surface field
!> @param[in]     map_2d                 Dofmap for cell at the base of the column
!> @param[in]     ndf_tile               Number of DOFs per cell for tiles
!> @param[in]     undf_tile              Number of total DOFs for tiles
!> @param[in]     map_tile               Dofmap for cell at the base of the column
subroutine calc_albedo_code(nlayers,                                &
                            surf_albedo,                            &
                            sw_down_surf,                           &
                            tile_fraction,                          &
                            sw_up_tile,                             &
                            ndf_2d, undf_2d, map_2d,                &
                            ndf_tile, undf_tile, map_tile)

  use jules_control_init_mod, only: n_surf_tile  

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: map_2d(ndf_2d)
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: map_tile(ndf_tile)

  real(r_def), intent(inout) :: surf_albedo(undf_2d)
  real(r_def), intent(in)    :: sw_down_surf(undf_2d)
  
  real(r_def), intent(in) :: tile_fraction(undf_tile)
  real(r_def), intent(in) :: sw_up_tile(undf_tile)

  ! Local variables for the kernel
  integer(i_def) :: i  
  real(r_def) :: sw_up_gb
  
  ! Calculate gridcell mean upward SW flux
  sw_up_gb = 0.0_r_def
  do i = 1, n_surf_tile
    sw_up_gb = sw_up_gb + ( tile_fraction( map_tile(1) + i - 1 ) *         &
                            sw_up_tile( map_tile(1) + i - 1 ) )        
  end do
  if ( sw_down_surf(map_2d(1)) > 0.0 ) then
    surf_albedo(map_2d(1)) = sw_up_gb / sw_down_surf(map_2d(1))
  end if  
  
  surf_albedo(map_2d(1)) = max(surf_albedo(map_2d(1)), 0.0_r_def)
  surf_albedo(map_2d(1)) = min(surf_albedo(map_2d(1)), 0.99_r_def)

  return
end subroutine calc_albedo_code

end module calc_albedo_kernel_mod
