!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Process Jules snow ancillaries
module process_snow_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_REAL,         &
                           GH_INTEGER,                &
                           GH_READ, GH_READWRITE,     &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2, &
                           ANY_DISCONTINUOUS_SPACE_3, &
                           CELL_COLUMN
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: process_snow_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_INTEGER,GH_READWRITE,ANY_DISCONTINUOUS_SPACE_1),&
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: process_snow_code
  end type process_snow_kernel_type

  public :: process_snow_code

contains

  !> @param[in]     nlayers                The number of layers
  !> @param[in,out] tile_snow_mass         Snow mass on tiles (kg m-2)
  !> @param[in,out] n_snow_layers          Number of snow layers on tiles
  !> @param[in,out] snow_under_canopy      Amount of snow under canopy (kg m-2)
  !> @param[in,out] snowpack_density       Density of snow on ground (kg m-3)
  !> @param[in,out] snow_depth             Snow depth on tiles (m)
  !> @param[in,out] snow_layer_thickness   Thickness of snow layers (m)
  !> @param[in,out] snow_layer_ice_mass    Mass of ice in snow layers (kg m-2)
  !> @param[in,out] snow_layer_liq_mass    Mass of liquid in snow layers (kg m-2)
  !> @param[in,out] snow_layer_temp        Temperature of snow layer (K)
  !> @param[in,out] snow_layer_rgrain      Grain radius of snow layer (microns)
  !> @param[in]     soil_temperature       Soil temperature (K)
  !> @param[in]     ndf_tile               Total DOFs per cell for surface tiles
  !> @param[in]     undf_tile              Unique DOFs per cell for surface tile
  !> @param[in]     map_tile               DOFmap for cells for surface tiles
  !> @param[in]     ndf_snow               Total DOFs per cell for snow layers
  !> @param[in]     undf_snow              Unique DOFs per cell for snow layers
  !> @param[in]     map_snow               DOFmap for cells for snow layers
  !> @param[in]     ndf_soil               Number of DOFs per cell for soil levs
  !> @param[in]     undf_soil              Number of total DOFs for soil levels
  !> @param[in]     map_soil               Dofmap for cell for soil levels
  subroutine process_snow_code(nlayers,                       &
                               tile_snow_mass,                &
                               n_snow_layers,                 &
                               snow_under_canopy,             &
                               snowpack_density,              &
                               snow_depth,                    &
                               snow_layer_thickness,          &
                               snow_layer_ice_mass,           &
                               snow_layer_liq_mass,           &
                               snow_layer_temp,               &
                               snow_layer_rgrain,             &
                               soil_temperature,              &
                               ndf_tile, undf_tile, map_tile, &
                               ndf_snow, undf_snow, map_snow, &
                               ndf_soil, undf_soil, map_soil)

    use jules_control_init_mod, only: n_land_tile
    use jules_snow_mod, only: nsmax, rho_snow_const, dzsnow, cansnowtile
    use jules_surface_types_mod, only: ice
    use jules_vegetation_mod, only: can_model
    use water_constants_mod, only: tm

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
    integer(kind=i_def), intent(in) :: ndf_snow, undf_snow
    integer(kind=i_def), intent(in) :: map_snow(ndf_snow)
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil)

    real(kind=r_def), intent(inout) :: tile_snow_mass(undf_tile)
    integer(kind=i_def), intent(inout) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(inout) :: snow_under_canopy(undf_tile)
    real(kind=r_def), intent(inout) :: snowpack_density(undf_tile)
    real(kind=r_def), intent(inout) :: snow_depth(undf_tile)

    real(kind=r_def), intent(inout) :: snow_layer_thickness(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_ice_mass(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_liq_mass(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_temp(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_rgrain(undf_snow)

    real(kind=r_def), intent(in)    :: soil_temperature(undf_soil)

    integer(kind=i_def) :: i, n, i_snow
    real(kind=r_def)    :: d_snow_mass, d_snow_thick

    real(kind=r_def), parameter :: snow_landice_min = 5.0e4_r_def
    real(kind=r_def), parameter :: snow_icefree_max = 5.0e3_r_def
    real(kind=r_def), parameter :: rgrain_icesheet = 150.0_r_def
    ! A value which is bigger than the UM MDI to ensure we don't modify
    ! sea-points
    real(kind=r_def), parameter :: snow_land_min = -1000.0_r_def

    do i = 1, n_land_tile

      ! Remove any small negative values in dump
      if (tile_snow_mass(map_tile(1)+i-1) < 0.0_r_def .and. &
          tile_snow_mass(map_tile(1)+i-1) > snow_land_min) then
        tile_snow_mass(map_tile(1)+i-1) = 0.0_r_def
      end if

      if (i == ice) then
        ! If this is an ice point, ensure that the mass of snow is large
        ! enough: if it has been reconfigured from an ice-free point this
        ! may not be the case.
        if (tile_snow_mass(map_tile(1)+i-1) >= 0.0_r_def .and. &
            tile_snow_mass(map_tile(1)+i-1) < snow_landice_min ) then

          if (n_snow_layers(map_tile(1)+i-1) < nsmax ) then
            ! The ice is so thin that the snow pack should be reset.
            tile_snow_mass(map_tile(1)+i-1) = snow_landice_min
            if (can_model == 4) snow_under_canopy(map_tile(1)+i-1) = 0.0_r_def
            snowpack_density(map_tile(1)+i-1) = rho_snow_const
            snow_depth(map_tile(1)+i-1) = snow_landice_min / rho_snow_const
            n_snow_layers(map_tile(1)+i-1) = nsmax
            do n = 1, nsmax
              i_snow = (i-1)*nsmax + n-1
              snow_layer_thickness(map_snow(1)+i_snow) = dzsnow(n)
              if (n==nsmax) then
                snow_layer_thickness(map_snow(1)+i_snow) = &
                     snow_depth(map_tile(1)+i-1) - sum(dzsnow(1:nsmax-1))
              end if
              ! Keep it frozen.
              snow_layer_ice_mass(map_snow(1)+i_snow) = &
                   rho_snow_const * snow_layer_thickness(map_snow(1)+i_snow)
              snow_layer_liq_mass(map_snow(1)+i_snow) = 0.0_r_def
              snow_layer_temp(map_snow(1)+i_snow) = &
                   min(soil_temperature(map_soil(1)), tm-1.0_r_def)
              snow_layer_rgrain(map_snow(1)+i_snow) = rgrain_icesheet
            end do
          else
            ! The snow pack is essentially correct. Simply reset the
            ! lowest layer.
            i_snow = (i-1)*nsmax + nsmax-1
            d_snow_mass = snow_landice_min - tile_snow_mass(map_tile(1)+i-1)
            d_snow_thick = d_snow_mass * &
                 snow_layer_thickness(map_snow(1)+i_snow) / &
                 ( snow_layer_ice_mass(map_snow(1)+i_snow) + &
                   snow_layer_liq_mass(map_snow(1)+i_snow) )
            snow_layer_ice_mass(map_snow(1)+i_snow) = &
                 snow_layer_ice_mass(map_snow(1)+i_snow) + d_snow_mass
            snow_layer_thickness(map_snow(1)+i_snow) = &
                 snow_layer_thickness(map_snow(1)+i_snow) + d_snow_thick
            tile_snow_mass(map_tile(1)+i-1) = snow_landice_min
            snow_depth(map_tile(1)+i-1) = snow_depth(map_tile(1)+i-1) + &
                                          d_snow_thick
            snowpack_density(map_tile(1)+i-1) = snow_landice_min / &
                                                snow_depth(map_tile(1)+i-1)
          end if
        end if

      else
        ! If this is an ice-free point, ensure that the mass of snow is not
        ! too large: if it has been reconfigured from a land ice point this
        ! may not be the case.
        if (tile_snow_mass(map_tile(1)+i-1) > snow_icefree_max ) then
          if (can_model == 4 .and. cansnowtile(i) ) then
            tile_snow_mass(map_tile(1)+i-1) = 0.0_r_def
            snow_under_canopy(map_tile(1)+i-1) = snow_icefree_max
          else
            tile_snow_mass(map_tile(1)+i-1) = snow_icefree_max
          end if
          snowpack_density(map_tile(1)+i-1) = rho_snow_const
          snow_depth(map_tile(1)+i-1) = snow_icefree_max / rho_snow_const
          n_snow_layers(map_tile(1)+i-1) = nsmax
          do n = 1, nsmax
            i_snow = (i-1)*nsmax + n-1
            snow_layer_thickness(map_snow(1)+i_snow) = dzsnow(n)
            if (n==nsmax) then
              snow_layer_thickness(map_snow(1)+i_snow) = &
                   snow_depth(map_tile(1)+i-1) - sum(dzsnow(1:nsmax-1))
            end if
            ! Keep it frozen.
            snow_layer_ice_mass(map_snow(1)+i_snow) = &
                 rho_snow_const * snow_layer_thickness(map_snow(1)+i_snow)
            snow_layer_liq_mass(map_snow(1)+i_snow) = 0.0_r_def
            snow_layer_temp(map_snow(1)+i_snow) = &
                 min(soil_temperature(map_soil(1)), tm-1.0_r_def)
            ! Use a typical value for ice sheets. The point is indeed
            ! not land ice, but will be close to land ice, so we assume
            ! that surface snow will have aged in much the same way.
            snow_layer_rgrain(map_snow(1)+i_snow) = rgrain_icesheet
          end do
        end if
      end if
    end do

  end subroutine process_snow_code

end module process_snow_kernel_mod
