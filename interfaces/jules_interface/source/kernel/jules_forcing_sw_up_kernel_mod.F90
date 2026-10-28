!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the upward SW flux on tiles for JULES site forcing
!>
module jules_forcing_sw_up_kernel_mod

  use argument_mod,       only : arg_type,                  &
                                 GH_FIELD, GH_REAL,         &
                                 GH_READ, GH_WRITE,         &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 ANY_DISCONTINUOUS_SPACE_2, &
                                 ANY_DISCONTINUOUS_SPACE_3, &
                                 CELL_COLUMN
  use constants_mod,      only : i_def, r_def
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: jules_forcing_sw_up_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),     & ! sw_up_tile
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),     & ! tile_sw_diffuse_albedo
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)      & ! sw_down_surf
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: jules_forcing_sw_up_code
  end type

  public :: jules_forcing_sw_up_code

contains

  !> @brief Calculates the upward SW flux on tiles for JULES site forcing
  !> @details The tile albedos are for a single band covering the whole SW
  !>          spectrum. As in standalone JULES with its default albedo
  !>          weights, the diffuse albedo is used for all of the downward SW.
  !> @param[in]     nlayers                Number of layers
  !> @param[in,out] sw_up_tile             Upward SW flux on tiles (W m-2)
  !> @param[in]     tile_sw_diffuse_albedo Diffuse SW albedo on tiles
  !> @param[in]     sw_down_surf           Downward SW flux at surface (W m-2)
  !> @param[in]     ndf_tile               Number of DOFs per cell for tiles
  !> @param[in]     undf_tile              Number of unique DOFs for tiles
  !> @param[in]     map_tile               Dofmap for the cell at the base of the column for tiles
  !> @param[in]     ndf_sw_tile            Number of DOFs per cell for SW tile albedos
  !> @param[in]     undf_sw_tile           Number of unique DOFs for SW tile albedos
  !> @param[in]     map_sw_tile            Dofmap for the cell at the base of the column for SW tile albedos
  !> @param[in]     ndf_2d                 Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d                Number of unique DOFs for 2D fields
  !> @param[in]     map_2d                 Dofmap for the cell at the base of the column for 2D fields
  subroutine jules_forcing_sw_up_code( nlayers,                            &
                                       sw_up_tile,                         &
                                       tile_sw_diffuse_albedo,             &
                                       sw_down_surf,                       &
                                       ndf_tile, undf_tile, map_tile,      &
                                       ndf_sw_tile, undf_sw_tile,          &
                                       map_sw_tile,                        &
                                       ndf_2d, undf_2d, map_2d )

    use jules_control_init_mod, only: n_surf_tile

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: ndf_sw_tile, undf_sw_tile
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), dimension(ndf_tile),    intent(in) :: map_tile
    integer(kind=i_def), dimension(ndf_sw_tile), intent(in) :: map_sw_tile
    integer(kind=i_def), dimension(ndf_2d),      intent(in) :: map_2d

    real(kind=r_def), intent(inout) :: sw_up_tile(undf_tile)
    real(kind=r_def), intent(in)    :: tile_sw_diffuse_albedo(undf_sw_tile)
    real(kind=r_def), intent(in)    :: sw_down_surf(undf_2d)

    ! Local variables
    integer(kind=i_def) :: n

    do n = 1, n_surf_tile
      sw_up_tile(map_tile(1)+n-1) = tile_sw_diffuse_albedo(map_sw_tile(1)+n-1) &
                                  * sw_down_surf(map_2d(1))
    end do

  end subroutine jules_forcing_sw_up_code

end module jules_forcing_sw_up_kernel_mod
