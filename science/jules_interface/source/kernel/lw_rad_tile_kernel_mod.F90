!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Jules for LW surface tile radiative properties

module lw_rad_tile_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_READ, GH_WRITE,         &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def, r_um, i_um
use kernel_mod,        only : kernel_type

implicit none

private

public :: lw_rad_tile_kernel_type
public :: lw_rad_tile_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: lw_rad_tile_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                     &
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_lw_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! tile_lw_grey_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ                            )  & ! n_band
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: lw_rad_tile_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                Number of layers
!> @param[in,out] tile_lw_albedo         LW tile albedos
!> @param[in,out] tile_lw_grey_albedo    LW tile grey albedos
!> @param[in]     tile_fraction          Surface tile fractions
!> @param[in]     tile_temperature       Surface tile temperature
!> @param[in]     n_band                 Number of spectral bands
!> @param[in]     ndf_lw_tile            DOFs per cell for tiles and lw bands
!> @param[in]     undf_lw_tile           Total DOFs for tiles and lw bands
!> @param[in]     map_lw_tile            Dofmap for cell at the base of the column
!> @param[in]     ndf_tile               Number of DOFs per cell for tiles
!> @param[in]     undf_tile              Number of total DOFs for tiles
!> @param[in]     map_tile               Dofmap for cell at the base of the column
subroutine lw_rad_tile_code(nlayers,                                &
                            tile_lw_albedo,                         &
                            tile_lw_grey_albedo,                    &
                            tile_fraction,                          &
                            tile_temperature,                       &
                            n_band,                                 &
                            ndf_lw_tile, undf_lw_tile, map_lw_tile, &
                            ndf_tile, undf_tile, map_tile)

  use socrates_init_mod, only: n_band_exclude, index_exclude, &
    wavelength_short, wavelength_long
  use jules_control_init_mod, only: &
    n_surf_tile, n_land_tile, n_sea_tile, n_sea_ice_tile, &
    first_sea_tile, first_sea_ice_tile
  use jules_surface_types_mod, only: soil
  use specemis_mod, only: specemis
  use surface_config_mod, only: &
    emis_method_sea,emis_method_sea_fixed,emis_method_sea_feldman,emis_method_sea_iremis, &
    emis_method_soil,emis_method_soil_fixed,emis_method_soil_feldman_desert

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, n_band
  integer(i_def), intent(in) :: ndf_lw_tile, undf_lw_tile
  integer(i_def), intent(in) :: map_lw_tile(ndf_lw_tile)
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: map_tile(ndf_tile)

  real(r_def), intent(inout) :: tile_lw_albedo(undf_lw_tile)

  real(r_def), intent(inout) :: tile_lw_grey_albedo(undf_tile)

  real(r_def), intent(in) :: tile_fraction(undf_tile)

  real(r_def), intent(in) :: tile_temperature(undf_tile)

  ! Local variables for the kernel
  integer(i_def) :: i_tile, i_band
  integer(i_def) :: df_rtile

  real(r_def) :: specemis_sea(n_band)
  real(r_def) :: specemis_soil(n_band)

  real(r_def) :: greyemis_sea
  real(r_def) :: greyemis_soil

  ! If spectrally varying albedos/emissivities are used for desert tiles
  ! call specemis and set tile_lw_grey_albedo and tile_lw_albedo in the
  ! later part of this routine accordingly
  if ( (tile_fraction(map_tile(1)+soil-1) > 0.0_r_def) .and. &
       (emis_method_soil == emis_method_soil_feldman_desert) ) then

    call specemis('desert_feldman', n_band, wavelength_short, wavelength_long, &
                  tile_temperature(map_tile(1)+soil-1), &
                  n_band_exclude, index_exclude, specemis_soil, greyemis_soil)

    tile_lw_grey_albedo(map_tile(1)+soil-1) &
      = 1.0_r_def - greyemis_soil

  end if


  ! Perform the same steps if spectrally varying albedos/emissivities
  ! are used for sea tiles (from feldman or iremis)
  if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def) then

    if (emis_method_sea == emis_method_sea_feldman) then
      call specemis('sea_feldman', n_band, wavelength_short, wavelength_long, &
                   tile_temperature(map_tile(1)+first_sea_tile-1), &
                   n_band_exclude, index_exclude, specemis_sea, greyemis_sea)

      tile_lw_grey_albedo(map_tile(1)+first_sea_tile-1) &
        = 1.0_r_def - greyemis_sea

    end if

    if (emis_method_sea == emis_method_sea_iremis) then
      call specemis('sea_iremis', n_band, wavelength_short, wavelength_long, &
                   tile_temperature(map_tile(1)+first_sea_tile-1), &
                   n_band_exclude, index_exclude, specemis_sea, greyemis_sea)

      tile_lw_grey_albedo(map_tile(1)+first_sea_tile-1) &
        = 1.0_r_def - greyemis_sea

    end if

  end if


  ! Now set the tile_lw_albedo for land, sea and sea-ice tiles
  ! If constant albedos/emissivities are used these can are copied from tile_lw_grey_albedo
  ! If spectrally varying emissivities are used these have been calculated
  ! per band in the calls to specemis above and are applied below to set tile_lw_albedo
  do i_band = 1, n_band

    ! Land tile albedos
    df_rtile = n_surf_tile * (i_band-1)
    do i_tile = 1, n_land_tile
      df_rtile = df_rtile + 1
      if ( (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def).and. &
           ((i_tile == soil) .and. (emis_method_soil /= emis_method_soil_fixed)) ) then
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = 1.0_r_def - specemis_soil(i_band)
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = tile_lw_grey_albedo(map_tile(1)+i_tile-1)
      end if
    end do

    ! Sea tile albedos
    df_rtile = first_sea_tile-1 + n_surf_tile * (i_band-1)
    do i_tile = first_sea_tile, first_sea_tile + n_sea_tile - 1
      df_rtile = df_rtile + 1
      if ( (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def).and. &
           (emis_method_sea /= emis_method_sea_fixed)) then
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = 1.0_r_def - specemis_sea(i_band)
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = tile_lw_grey_albedo(map_tile(1)+i_tile-1)
      end if
    end do

    ! Sea-ice tile albedos
    df_rtile = first_sea_ice_tile - 1 + n_surf_tile * (i_band-1)
    do i_tile = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      df_rtile = df_rtile + 1
      tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
        = tile_lw_grey_albedo(map_tile(1)+i_tile-1)
    end do

  end do

end subroutine lw_rad_tile_code

end module lw_rad_tile_kernel_mod
