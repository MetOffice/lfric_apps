!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Initialise Jules surface fields on tiles
!> @details Non-standard Surface fields (pseudo-levels) aren't as yet not
!>  implemented in LFRic. As an interim measure Higher-order W3 fields have
!>  been used to mimic psuedo-level field behaviour. This code is written
!>  based on this interim measure and will need to be updated when
!>  suitable infrastructure is available (Ticket #2081)
module initial_surface_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_REAL,         &
                           GH_WRITE, CELL_COLUMN,     &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  use jules_control_init_mod, only: n_land_tile, &
                                    n_sea_tile, n_sea_ice_tile, &
       first_sea_tile, first_sea_ice_tile, n_surf_tile
  use ideal_surface_config_mod, only: surf_tile_fracs, surf_tile_temps, &
                                      canopy_height_in => canopy_height, leaf_area_index_in => leaf_area_index
  use jules_surface_types_mod, only: npft
  use nvegparm, only: emis_nvg
  use pftparm, only: emis_pft
  use jules_sea_seaice_mod, only: emis_sea, emis_sice
  use section_choice_config_mod, only: surface, surface_jules

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: initial_surface_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: initial_surface_code
  end type initial_surface_kernel_type

  public :: initial_surface_code

contains

  !> @param[in]     nlayers             The number of layers
  !> @param[in,out] tile_fraction       Surface tile fractions
  !> @param[in,out] leaf_area_index     Leaf Area Index
  !> @param[in,out] canopy_height       Canopy height (m)
  !> @param[in,out] tile_temperature    Surface tile temperatures (K)
  !> @param[in,out] tile_lw_grey_albedo Surface tile grey albedo
  !> @param[in]     ndf_tile            Number of DOFs per cell for tiles
  !> @param[in]     undf_tile           Number of total DOFs for tiles
  !> @param[in]     map_tile            Dofmap for cell for surface tiles
  !> @param[in]     ndf_pft             Number of DOFs per cell for PFTs
  !> @param[in]     undf_pft            Number of total DOFs for PFTs
  !> @param[in]     map_pft             Dofmap for cell for PFTs
  subroutine initial_surface_code(nlayers,                       &
                                  tile_fraction,                 &
                                  leaf_area_index,               &
                                  canopy_height,                 &
                                  tile_temperature,              &
                                  tile_lw_grey_albedo,           &
                                  ndf_tile, undf_tile, map_tile, &
                                  ndf_pft, undf_pft, map_pft)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
    integer(kind=i_def), intent(in) :: ndf_pft, undf_pft
    integer(kind=i_def), intent(in) :: map_pft(ndf_pft)

    real(kind=r_def), intent(inout) :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(inout) :: tile_lw_grey_albedo(undf_tile)

    real(kind=r_def), intent(inout) :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(inout) :: canopy_height(undf_pft)

    ! Internal variables
    integer(kind=i_def) :: i

    do i = 1, n_surf_tile
      tile_fraction(map_tile(1)+i-1) = surf_tile_fracs(i)
      tile_temperature(map_tile(1)+i-1) = surf_tile_temps(i)
    end do

    if (surface == surface_jules) then

      ! Set from a namelist variable when no ancil read
      do i = 1, npft
        leaf_area_index(map_pft(1)+i-1) = leaf_area_index_in(i)
        canopy_height(map_pft(1)+i-1) = canopy_height_in(i)
      end do

      ! Initialise grey_albedos for all tiles
      do i = 1, n_land_tile
        if (i <= npft) then
          tile_lw_grey_albedo(map_tile(1)+i-1) = 1.0_r_def - real(emis_pft(i), r_def)
        else
          tile_lw_grey_albedo(map_tile(1)+i-1) = 1.0_r_def - real(emis_nvg(i-npft), r_def)
        endif
      end do

      do i = first_sea_tile, first_sea_tile + n_sea_tile - 1
        tile_lw_grey_albedo(map_tile(1)+i-1) = 1.0_r_def - real(emis_sea, r_def)
      end do

      do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        tile_lw_grey_albedo(map_tile(1)+i-1) = 1.0_r_def - real(emis_sice, r_def)
      end do

    end if

  end subroutine initial_surface_code

end module initial_surface_kernel_mod
