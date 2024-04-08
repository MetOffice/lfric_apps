!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Kernel which converts sea ice fractions from a fraction of whole grid
!> @brief box to fraction of the marine portion of the grid box

module convert_to_marine_fraction_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,                  &
                                    GH_REAL,                   &
                                    GH_INTEGER,                &
                                    GH_FIELD,                  &
                                    GH_SCALAR,                 &
                                    GH_WRITE, GH_READ,         &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: convert_to_marine_fraction_type
  private
  type(arg_type) :: meta_args(6) = (/                                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: convert_to_marine_fraction_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: convert_to_marine_fraction_code
contains

!> @brief Kernel which converts sea ice fractions from a fraction whole grid
!> @brief box to fraction of the marine portion of the grid box
!! @param[in] nlayers Number of layers
!! @param[in,out] sea_ice_marine_fraction Output sea ice as a marine fraction (on sea ice categories)
!! @param[in,out] sea_marine_fraction Output sea area as a marine fraction
!! @param[in] tile_fraction   Input tile fractions (on tiles)
!! @param[in] first_sea_ice_tile Index of the first sea ice tile
!! @param[in] n_sea_ice_tile Number of sea ice tiles
!! @param[in] first_sea_tile Index of the first sea tile
!! @param[in] ndf_sea_ice Number of degrees of freedom for sea ice fraction
!! @param[in] undf_sea_ice Number of unique degrees of freedom for sea ice fraction
!! @param[in] map_sea_ice Dofmap for sea ice fraction
!! @param[in] ndf_sea Number of degrees of freedom for sea ice fraction
!! @param[in] undf_sea Number of unique degrees of freedom for sea ice fraction
!! @param[in] map_sea Dofmap for sea ice fraction
!! @param[in] ndf_tile Number of degrees of freedom for tile fraction
!! @param[in] undf_tile Number of unique degrees of freedom for tile fraction
!! @param[in] map_tile Dofmap for tile fraction

subroutine convert_to_marine_fraction_code(nlayers,                  &
                      sea_ice_marine_fraction, sea_marine_fraction,  &
                      tile_fraction,                                 &
                      first_sea_ice_tile,                            &
                      n_sea_ice_tile,                                &
                      first_sea_tile,                                &
                      ndf_sea_ice, undf_sea_ice, map_sea_ice,        &
                      ndf_sea, undf_sea, map_sea,                    &
                      ndf_tile, undf_tile, map_tile)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_sea_ice
  integer(kind=i_def), intent(in) :: undf_sea_ice
  integer(kind=i_def), intent(in) :: map_sea_ice(ndf_sea_ice)
  integer(kind=i_def), intent(in) :: ndf_sea
  integer(kind=i_def), intent(in) :: undf_sea
  integer(kind=i_def), intent(in) :: map_sea(ndf_sea)
  integer(kind=i_def), intent(in) :: ndf_tile
  integer(kind=i_def), intent(in) :: undf_tile
  integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
  real(kind=r_def), intent(inout) :: sea_ice_marine_fraction(undf_sea_ice)
  real(kind=r_def), intent(inout) :: sea_marine_fraction(undf_sea)
  real(kind=r_def), intent(in)    :: tile_fraction(undf_tile)
  integer(kind=i_def), intent(in) :: first_sea_ice_tile
  integer(kind=i_def), intent(in) :: n_sea_ice_tile
  integer(kind=i_def), intent(in) :: first_sea_tile


  ! Local variables
  integer(kind=i_def) :: i
  real(kind=r_def)    :: total_marine_fraction

  ! Calculate the total marine fraction by adding together
  ! the ocean and sea ice fractions.
  total_marine_fraction = tile_fraction(map_tile(1) + first_sea_tile - 1)
  do i = first_sea_ice_tile, first_sea_ice_tile+n_sea_ice_tile-1
    total_marine_fraction = total_marine_fraction + tile_fraction(map_tile(1) + i-1)
  end do

  ! Divide by the total marine fraction to get the sea ice fractions
  ! as a fraction of the marine portion of the grid box
  do i = first_sea_ice_tile, first_sea_ice_tile+n_sea_ice_tile-1
    if (total_marine_fraction > 0.0_r_def) then
      sea_ice_marine_fraction(map_sea_ice(1) + i-first_sea_ice_tile) =         &
                       tile_fraction(map_tile(1) + i-1) / total_marine_fraction
    else
      sea_ice_marine_fraction(map_sea_ice(1) + i-first_sea_ice_tile) = 0.0_r_def
    end if
  end do

  ! Divide by the total marine fraction to get the sea fraction
  ! as a fraction of the marine portion of the grid box
  if (total_marine_fraction > 0.0_r_def) then
    sea_marine_fraction(map_sea(1)) =                                          &
       tile_fraction(map_tile(1) + first_sea_tile - 1) / total_marine_fraction
  else
    sea_marine_fraction(map_sea(1)) = 0.0_r_def
  end if

end subroutine convert_to_marine_fraction_code

end module convert_to_marine_fraction_mod
