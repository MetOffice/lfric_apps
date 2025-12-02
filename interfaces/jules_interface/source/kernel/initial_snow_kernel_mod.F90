!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Initialise Jules snow fields on snow levels
!> @details Non-standard Surface fields (pseudo-levels) aren't as yet not
!>  implemented in LFRic. As an interim measure Higher-order W3 fields have
!>  been used to mimic psuedo-level field behaviour. This code is written
!>  based on this interim measure and will need to be updated when
!>  suitable infrastructure is available (Ticket #2081)
module initial_snow_kernel_mod

  use argument_mod,  only: arg_type,                      &
                           GH_FIELD, GH_REAL, GH_INTEGER, &
                           GH_WRITE, CELL_COLUMN,         &
                           ANY_DISCONTINUOUS_SPACE_1,     &
                           ANY_DISCONTINUOUS_SPACE_2
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  use jules_control_init_mod, only: n_land_tile, n_surf_tile
  use jules_snow_mod,         only: nsmax
  use ideal_surface_config_mod, only: tile_snow_mass_in => tile_snow_mass, &
       n_snow_layers_in => n_snow_layers, snow_depth_in => snow_depth, &
       snow_layer_thickness_in => snow_layer_thickness, &
       snow_layer_ice_mass_in => snow_layer_ice_mass, &
       snow_layer_temp_in => snow_layer_temp

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: initial_snow_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                       &
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: initial_snow_code
  end type initial_snow_kernel_type

  public :: initial_snow_code

contains

  !> @param[in]     nlayers              The number of layers
  !> @param[in,out] tile_snow_mass       Mass of snow on surface tiles (kg)
  !> @param[in,out] n_snow_layers        Number of snow layers on tiles
  !> @param[in,out] snow_depth           Snow depth on surface tiles (m)
  !> @param[in,out] snow_layer_thickness Thickness of snow layers (m)
  !> @param[in,out] snow_layer_ice_mass  Mass of ice in snow layers (kg m-2)
  !> @param[in,out] snow_layer_temp      Temperature of snow layer (K)
  !> @param[in]     ndf_snow             Number of DOFs per cell for snow
  !> @param[in]     undf_snow            Number of total DOFs for snow
  !> @param[in]     map_snow             Dofmap for cell for snow fields
  subroutine initial_snow_code(nlayers,                       &
                               tile_snow_mass,                &
                               n_snow_layers,                 &
                               snow_depth,                    &
                               snow_layer_thickness,          &
                               snow_layer_ice_mass,           &
                               snow_layer_temp,               &
                               ndf_tile, undf_tile, map_tile, &
                               ndf_snow, undf_snow, map_snow)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_snow, undf_snow
    integer(kind=i_def), intent(in) :: map_snow(ndf_snow)
    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)

    real(kind=r_def), intent(inout) :: tile_snow_mass(undf_tile)
    integer(kind=i_def), intent(inout) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(inout) :: snow_depth(undf_tile)
    real(kind=r_def), intent(inout) :: snow_layer_thickness(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_ice_mass(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_temp(undf_snow)

    ! Internal variables
    integer(kind=i_def) :: i, j, i_snow, indexes(nsmax)

    ! Initialise snow prognostics for SCM testing or when no values are
    ! provided by um2lfric start dump
    do i = 1, n_surf_tile
      tile_snow_mass(map_tile(1)+i-1) = tile_snow_mass_in(i)
      n_snow_layers(map_tile(1)+i-1) = n_snow_layers_in(i)
      snow_depth(map_tile(1)+i-1) = snow_depth_in(i)
    end do

    ! Set from a namelist variable when no start dump read
    i_snow = 0
    do i = 1, n_land_tile

      do j = 1, nsmax
        i_snow = i_snow + 1

        snow_layer_thickness(map_snow(1)+i_snow-1) = &
             snow_layer_thickness_in(i+(j-1)*n_land_tile)
        snow_layer_ice_mass(map_snow(1)+i_snow-1) = &
             snow_layer_ice_mass_in(i+(j-1)*n_land_tile)
        snow_layer_temp(map_snow(1)+i_snow-1) = &
             snow_layer_temp_in(i+(j-1)*n_land_tile)

      end do
    end do

  end subroutine initial_snow_code

end module initial_snow_kernel_mod
