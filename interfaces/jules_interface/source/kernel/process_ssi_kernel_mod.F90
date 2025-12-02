!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Process Jules sea and sea-ice ancillaries
!> @details Kernel used without proper PSyclone support of multi-data fields
!>          see https://github.com/stfc/PSyclone/issues/868
module process_ssi_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_REAL,         &
                           GH_READ, GH_READWRITE,     &
                           CELL_COLUMN,               &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2, &
                           ANY_DISCONTINUOUS_SPACE_3, &
                           GH_SCALAR, GH_LOGICAL
  use constants_mod, only: r_def, i_def, l_def
  use kernel_mod,    only: kernel_type

  use jules_control_init_mod, only: n_sea_ice_tile,     &
                                    first_sea_tile,     &
                                    first_sea_ice_tile, &
                                    n_land_tile,        &
                                    n_surf_tile
  use jules_physics_init_mod, only : min_sea_ice_frac

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: process_ssi_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                             &
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: process_ssi_code
  end type process_ssi_kernel_type

  public :: process_ssi_code

contains

  !> @param[in]     nlayers            The number of layers
  !> @param[in]     sea_ice_fraction   Fraction of sea-ice in grid-box
  !> @param[in,out  sea_ice_thickness  Thickness of the sea-ice
  !> @param[in,out] tile_fraction      Surface tile fractions
  !> @param[in]     latitude           Latitude of the grid-cell
  !> @param[in]     amip_ice_thick     Use AMIP ice thickness calculation
  !> @param[in]     l_couple_sea_ice   Model is coupled to a sea-ice model
  !> @param[in]     ndf_sice           Number of DOFs per cell for sea ice
  !> @param[in]     undf_sice          Number of total DOFs for sea ice
  !> @param[in]     map_sice           Dofmap for cell for surface sea ice
  !> @param[in]     ndf_tile           Number of DOFs per cell for tiles
  !> @param[in]     undf_tile          Number of total DOFs for tiles
  !> @param[in]     map_tile           Dofmap for cell for surface tiles
  !> @param[in]     ndf_2d             Number of DOFs per cell for 2d
  !> @param[in]     undf_2d            Number of total DOFs for 2d
  !> @param[in]     map_2d             Dofmap for cell for 2d
  subroutine process_ssi_code(nlayers,                       &
                              sea_ice_fraction,              &
                              sea_ice_thickness,             &
                              tile_fraction,                 &
                              latitude,                      &
                              amip_ice_thick,                &
                              l_couple_sea_ice,              &
                              ndf_sice, undf_sice, map_sice, &
                              ndf_tile, undf_tile, map_tile, &
                              ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
    integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
    integer(kind=i_def), intent(in) :: map_sice(ndf_sice)
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

    real(kind=r_def), intent(in)     :: sea_ice_fraction(undf_sice)
    real(kind=r_def), intent(inout)  :: sea_ice_thickness(undf_sice)
    real(kind=r_def), intent(inout)  :: tile_fraction(undf_tile)
    real(kind=r_def), intent(in)     :: latitude(undf_2d)

    logical(kind=l_def), intent(in) :: amip_ice_thick, l_couple_sea_ice

    ! Internal variables
    integer(kind=i_def) :: i, i_sice
    real(kind=r_def) :: tot_ice, tot_land, ice_thick

    ! Calculate the current ice fraction for use below
    tot_ice = 0.0_r_def
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      tot_ice = tot_ice + tile_fraction(map_tile(1)+i-1)
    end do

    ! Set up the land fraction from the input data
    tot_land = min(sum(tile_fraction(map_tile(1):map_tile(1)+n_land_tile-1)),1.0_r_def)
    if (tot_land < 1.0_r_def .and. &
         tile_fraction(map_tile(1)+first_sea_tile-1) == 0.0_r_def .and. &
         tot_ice == 0.0_r_def ) then
      tot_land = 1.0_r_def
    end if

    ! Set the new sea ice fraction from an ancillary
    tot_ice = 0.0_r_def
    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      ! Only use where field contains valid data
      if (sea_ice_fraction(map_sice(1)+i_sice-1) > min_sea_ice_frac .and. &
          tot_land < 1.0_r_def) then
        tile_fraction(map_tile(1)+i-1) = sea_ice_fraction(map_sice(1)+i_sice-1)&
                                       * (1.0_r_def - tot_land)
      else
        tile_fraction(map_tile(1)+i-1) = 0.0_r_def
      end if
      tot_ice = tot_ice + tile_fraction(map_tile(1)+i-1)
    end do

    ! Now set the sea fraction
    tile_fraction(map_tile(1)+first_sea_tile-1) = max(1.0_r_def &
                                                - tot_land &
                                                - tot_ice, 0.0_r_def)

    ! Calculate sea-ice thickness if required
    if (amip_ice_thick .and. .not. l_couple_sea_ice) then
      if (latitude(map_2d(1)) > 0.0_r_def) then
        ice_thick = 2.0_r_def ! Arctic sea-ice
      else
        ice_thick = 1.0_r_def ! Antarctic sea-ice
      end if

      if (tot_land < 1.0_r_def) then
        i_sice = 0
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = i_sice + 1
          sea_ice_thickness(map_sice(1)+i_sice-1) = ice_thick * &
                                          tot_ice / (1.0_r_def - tot_land)
        end do
      else
        i_sice = 0
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = i_sice + 1
          sea_ice_thickness(map_sice(1)+i_sice-1) = 0.0_r_def
        end do
      end if

    end if

  end subroutine process_ssi_code

end module process_ssi_kernel_mod
