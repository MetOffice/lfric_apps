!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Socrates model timestep interpolation for longwave (thermal) fluxes

module lw_mts_kernel_mod

use argument_mod,      only : arg_type, &
                              GH_FIELD, GH_SCALAR, &
                              GH_REAL, GH_LOGICAL, &
                              GH_READ, GH_WRITE, &
                              DOMAIN, &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: lw_mts_kernel_type
  private
  type(arg_type) :: meta_args(15) = (/ &
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, Wtheta),                    & ! lw_heating_rate
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_up_surf
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_up_toa
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! lw_up_tile
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! lw_heating_rate_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! lw_up_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! lw_up_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! lw_up_tile_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta),                    & ! exner_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_lw_grey_albedo
    arg_type(GH_SCALAR, GH_LOGICAL, GH_READ                            )  & ! rad_this_tstep
    /)
  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: lw_mts_code
end type

public :: lw_mts_code

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                  Number of layers
!> @param[in]     n_profile                Number of columns
!> @param[in,out] lw_heating_rate          LW heating rate
!> @param[in,out] lw_down_surf             LW downward surface flux
!> @param[in,out] lw_up_surf               LW upward surface flux
!> @param[in,out] lw_up_toa                LW upward top-of-atmosphere flux
!> @param[in,out] lw_up_tile               LW upward tiled surface flux
!> @param[in]     lw_heating_rate_rts      LW heating rate
!> @param[in]     lw_down_surf_rts         LW downward surface flux
!> @param[in]     lw_up_surf_rts           LW upward surface flux
!> @param[in]     lw_up_toa_rts            LW upward top-of-atmosphere flux
!> @param[in]     lw_up_tile_rts           LW upward tiled surface flux
!> @param[in]     exner_in_wth             Exner pressure in wth space
!> @param[in]     tile_fraction            Surface tile fractions
!> @param[in]     tile_temperature         Surface tile temperature
!> @param[in]     tile_lw_grey_albedo      LW tile grey albedos
!> @param[in]     rad_this_tstep           Full radiation call this timestep
!> @param[in]     ndf_wth                  No. DOFs per cell for wth space
!> @param[in]     undf_wth                 No. unique of DOFs for wth space
!> @param[in]     map_wth                  Dofmap for wth space column base cell
!> @param[in]     ndf_2d                   No. of DOFs per cell for 2D space
!> @param[in]     undf_2d                  No. unique of DOFs for 2D space
!> @param[in]     map_2d                   Dofmap for 2D space column base cell
!> @param[in]     ndf_tile                 Number of DOFs per cell for tiles
!> @param[in]     undf_tile                Number of total DOFs for tiles
!> @param[in]     map_tile                 Dofmap for tile space column base cell
subroutine lw_mts_code(nlayers, n_profile,                                     &
                   lw_heating_rate, lw_down_surf, lw_up_surf,                  &
                   lw_up_toa, lw_up_tile,                                      &
                   lw_heating_rate_rts, lw_down_surf_rts, lw_up_surf_rts,      &
                   lw_up_toa_rts, lw_up_tile_rts,                              &
                   exner_in_wth,                                               &
                   tile_fraction, tile_temperature,                            &
                   tile_lw_grey_albedo,                                        &
                   rad_this_tstep,                                             &
                   ndf_wth, undf_wth, map_wth,                                 &
                   ndf_2d, undf_2d, map_2d,                                    &
                   ndf_tile, undf_tile, map_tile)

  use jules_control_init_mod, only: n_surf_tile
  use socrates_bones, only: bones

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, n_profile
  integer(i_def), intent(in) :: ndf_wth, ndf_2d, ndf_tile
  integer(i_def), intent(in) :: undf_wth, undf_2d, undf_tile

  integer(i_def), dimension(ndf_wth, n_profile),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d, n_profile),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile, n_profile),  intent(in) :: map_tile

  real(r_def), dimension(undf_wth),  intent(inout) :: lw_heating_rate
  real(r_def), dimension(undf_2d),   intent(inout) :: lw_down_surf, &
    lw_up_surf, lw_up_toa
  real(r_def), dimension(undf_tile), intent(inout) :: lw_up_tile

  real(r_def), dimension(undf_wth),  intent(in) :: lw_heating_rate_rts
  real(r_def), dimension(undf_2d),   intent(in) :: lw_down_surf_rts, &
    lw_up_surf_rts, lw_up_toa_rts
  real(r_def), dimension(undf_tile), intent(in) :: lw_up_tile_rts

  real(r_def), dimension(undf_wth),  intent(in) :: exner_in_wth

  real(r_def), dimension(undf_tile), intent(in) :: tile_fraction
  real(r_def), dimension(undf_tile), intent(in) :: tile_temperature
  real(r_def), dimension(undf_tile), intent(in) :: tile_lw_grey_albedo

  logical(l_def), intent(in) :: rad_this_tstep

  ! Local variables for the kernel
  integer(i_def) :: k
  integer(i_def) :: wth_0, wth_1, wth_last
  integer(i_def) :: tile_1, tile_last
  integer(i_def) :: twod_1, twod_last


  ! Set indexing
  wth_0 = map_wth(1,1)
  wth_1 = map_wth(1,1)+1
  wth_last = map_wth(1,1)+n_profile*(nlayers+1)-1
  tile_1 = map_tile(1,1)
  tile_last = map_tile(1,1)+n_profile*n_surf_tile-1
  twod_1 = map_2d(1,1)
  twod_last = map_2d(1,1)+n_profile-1

  if (rad_this_tstep) then
    ! No corrections needed for model timestep
    lw_heating_rate(wth_0:wth_last) = lw_heating_rate_rts(wth_0:wth_last)
    lw_down_surf(twod_1:twod_last)  = lw_down_surf_rts(twod_1:twod_last)
    lw_up_surf(twod_1:twod_last)    = lw_up_surf_rts(twod_1:twod_last)
    lw_up_toa(twod_1:twod_last)     = lw_up_toa_rts(twod_1:twod_last)
    lw_up_tile(tile_1:tile_last)    = lw_up_tile_rts(tile_1:tile_last)
  else
    ! Not a radiation time-step: apply corrections to the model timestep

    ! The bare "bones" of a simple radiative transfer calculation.
    ! Update the upward fluxes for the change in surface temperature.
    call bones(n_profile, nlayers,                                             &
      n_layer_stride         = nlayers+1,                                      &
      n_tile                 = n_surf_tile,                                    &
      l_grey_emis_correction = .true.,                                         &
      grey_albedo_tile_1d    = tile_lw_grey_albedo(tile_1:tile_last),          &
      frac_tile_1d           = tile_fraction(tile_1:tile_last),                &
      t_tile_1d              = tile_temperature(tile_1:tile_last),             &
      heating_rate_1d_rts    = lw_heating_rate_rts(wth_1:wth_last),            &
      flux_down_surf_rts     = lw_down_surf_rts(twod_1:twod_last),             &
      flux_up_surf_rts       = lw_up_surf_rts(twod_1:twod_last),               &
      flux_up_toa_rts        = lw_up_toa_rts(twod_1:twod_last),                &
      heating_rate_1d_mts    = lw_heating_rate(wth_1:wth_last),                &
      flux_down_surf_mts     = lw_down_surf(twod_1:twod_last),                 &
      flux_up_surf_mts       = lw_up_surf(twod_1:twod_last),                   &
      flux_up_toa_mts        = lw_up_toa(twod_1:twod_last),                    &
      flux_up_tile_1d_mts    = lw_up_tile(tile_1:tile_last),                   &
      l_profile_last         = .true.)
  end if

  ! Set level 0 increment such that theta increment will equal level 1
  do k=wth_0, wth_last, nlayers+1
    lw_heating_rate(k) = lw_heating_rate(k+1) &
                       * exner_in_wth(k) / exner_in_wth(k+1)
  end do

end subroutine lw_mts_code
end module lw_mts_kernel_mod
