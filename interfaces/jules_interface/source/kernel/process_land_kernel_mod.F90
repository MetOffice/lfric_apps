!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Process Jules land ancillaries
!> @details Kernel used without proper PSyclone support of multi-data fields
!>          see https://github.com/stfc/PSyclone/issues/868
module process_land_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_REAL,         &
                           GH_READ,                   &
                           GH_WRITE, CELL_COLUMN,     &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2, &
                           ANY_DISCONTINUOUS_SPACE_3
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  use jules_control_init_mod, only: n_land_tile, first_sea_tile, &
                                    first_sea_ice_tile, n_sea_ice_tile

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: process_land_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: process_land_code
  end type process_land_kernel_type

  public :: process_land_code

contains

  !> @param[in]     nlayers            The number of layers
  !> @param[in]     land_area_fraction Fraction of land in grid-box
  !> @param[in]     land_tile_fraction Land tile fractions from ancil
  !> @param[in,out] tile_fraction      Surface tile fractions
  !> @param[in]     ndf_2d             Number of DOFs per cell for 2d fields
  !> @param[in]     undf_2d            Number of total DOFs for 2d fields
  !> @param[in]     map_2d             Dofmap for cell for surface 2d fields
  !> @param[in]     ndf_land           Number of DOFs per cell for land
  !> @param[in]     undf_land          Number of total DOFs for land
  !> @param[in]     map_land           Dofmap for cell for surface land
  !> @param[in]     ndf_tile           Number of DOFs per cell for tiles
  !> @param[in]     undf_tile          Number of total DOFs for tiles
  !> @param[in]     map_tile           Dofmap for cell for surface tiles
  subroutine process_land_code(nlayers,                       &
                               land_area_fraction,            &
                               land_tile_fraction,            &
                               tile_fraction,                 &
                               ndf_2d, undf_2d, map_2d,       &
                               ndf_land, undf_land, map_land, &
                               ndf_tile, undf_tile, map_tile)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
    integer(kind=i_def), intent(in) :: ndf_land, undf_land
    integer(kind=i_def), intent(in) :: map_land(ndf_land)
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

    real(kind=r_def), intent(in)     :: land_area_fraction(undf_2d)
    real(kind=r_def), intent(in)     :: land_tile_fraction(undf_land)
    real(kind=r_def), intent(inout)  :: tile_fraction(undf_tile)

    ! Internal variables
    integer(kind=i_def) :: i
    real(kind=r_def) :: tot_land

    ! Set land tile fractions from ancillary
    if (land_area_fraction(map_2d(1)) > 0.0_r_def) then
      tot_land = sum(land_tile_fraction(map_land(1):map_land(1)+n_land_tile-1))
      do i = 1, n_land_tile
        tile_fraction(map_tile(1)+i-1) = land_tile_fraction(map_land(1)+i-1) &
                                       * land_area_fraction(map_2d(1))       &
                                       / tot_land
      end do
    else
      tile_fraction(map_tile(1):map_tile(1)+n_land_tile-1) = 0.0_r_def
    end if

    ! Ensure these contain 0 where land=1
    tile_fraction(map_tile(1)+first_sea_tile-1) = &
         1.0_r_def - land_area_fraction(map_2d(1))
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      tile_fraction(map_tile(1)+i-1) = 1.0_r_def - land_area_fraction(map_2d(1))
    end do

  end subroutine process_land_code

end module process_land_kernel_mod
