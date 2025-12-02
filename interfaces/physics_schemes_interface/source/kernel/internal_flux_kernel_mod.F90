!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adds a surface temperature forcing via prescribed internal flux.
!>
!> @details Kernel that changes the surface temperature on tiles according to
!!          a map of prescribed internal flux from the planet's interior using
!!          the blackbody radiation law. The flux is prescribed via an ancillary
!!          file or a namelist value, and is positive when the flux is directed
!!          "up", i.e. from the planet's interior to the atmosphere.

module internal_flux_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_INTEGER, GH_REAL,       &
                              GH_READ, GH_WRITE,         &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

public :: internal_flux_kernel_type
public :: internal_flux_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: internal_flux_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                 &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_temperature
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! lw_down_surf
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! sw_down_surf
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_tile
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! internal_flux
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         & ! n_surf_tile
    arg_type(GH_SCALAR, GH_REAL, GH_READ)                             & ! planet_emissivity
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: internal_flux_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers           Number of layers
!> @param[in,out] tile_temperature  Surface temperature on tiles
!> @param[in]     lw_down_surf      LW downward radiation flux at the surface
!> @param[in]     sw_down_surf      SW downward radiation flux at the surface
!> @param[in]     sw_up_tile        SW upward radiation flux at the surface
!> @param[in]     internal_flux     Prescribed internal flux at the surface
!> @param[in]     n_surf_tile       Number of surface tiles
!> @param[in]     planet_emissivity Planet emissivity
!> @param[in]     ndf_2d            No. of degrees of freedom per cell for 2D space
!> @param[in]     undf_2d           No. unique of degrees of freedom for 2D space
!> @param[in]     map_2d            Dofmap for cell at base of column for 2D space
!> @param[in]     ndf_tile          Number of DOFs per cell for tiles
!> @param[in]     undf_tile         Number of total DOFs for tiles
!> @param[in]     map_tile          Dofmap for cell at the base of the column
subroutine internal_flux_code(nlayers,                       &
                              tile_temperature,              &
                              lw_down_surf,                  &
                              sw_down_surf,                  &
                              sw_up_tile,                    &
                              internal_flux,                 &
                              n_surf_tile,                   &
                              planet_emissivity,             &
                              ndf_tile, undf_tile, map_tile, &
                              ndf_2d, undf_2d, map_2d)

  use science_chemistry_constants_mod, only : stefan_boltzmann

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: n_surf_tile
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: map_2d(ndf_2d)
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: map_tile(ndf_tile)

  real(r_def), intent(inout) :: tile_temperature(undf_tile)
  real(r_def), intent(in)    :: lw_down_surf(undf_2d)
  real(r_def), intent(in)    :: sw_down_surf(undf_2d)
  real(r_def), intent(in)    :: sw_up_tile(undf_tile)
  real(r_def), intent(in)    :: internal_flux(undf_2d)
  real(r_def), intent(in)    :: planet_emissivity

  ! Local variables for the kernel
  integer(i_def) :: i_tile

  real(r_def) :: rad_sum

  do i_tile = 1, n_surf_tile
    ! Net radiation on tiles, exluding the LW up part, which is
    ! effectively diagnosed by this kernel
    rad_sum = sw_down_surf(map_2d(1)) - sw_up_tile(map_tile(1)+i_tile-1) &
      + lw_down_surf(map_2d(1)) * planet_emissivity

    ! Updated surface temperature
    ! T_surf = ( (F_int + F_down_rad) / (eps * sigma) )^1/4
    tile_temperature(map_tile(1)+i_tile-1) = &
      ( (internal_flux(map_2d(1)) + rad_sum) / &
        (planet_emissivity * stefan_boltzmann) ) ** 0.25_r_def
  end do

end subroutine internal_flux_code

end module internal_flux_kernel_mod
