!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Socrates model timestep interpolation for shortwave fluxes

module sw_mts_kernel_mod

use argument_mod,      only : arg_type, &
                              GH_FIELD, &
                              GH_REAL, &
                              GH_READ, GH_WRITE, &
                              DOMAIN, &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: sw_mts_kernel_type
  private
  type(arg_type) :: meta_args(26) = (/ &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! sw_heating_rate
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_surf
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_blue_surf
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_blue_surf
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_surf
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_toa
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_toa
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_tile
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_blue_tile
    arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! sw_heating_rate_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_surf_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_blue_surf_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_blue_surf_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_surf_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_toa_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_toa_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_tile_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_blue_tile_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! exner_in_wth
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1)  & ! orographic_correction_rts
    /)
  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: sw_mts_code
end type

public :: sw_mts_code

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                   Number of layers
!> @param[in]     n_profile                 Number of columns
!> @param[in,out] sw_heating_rate           SW heating rate
!> @param[in,out] sw_down_surf              SW downward surface flux
!> @param[in,out] sw_direct_surf            SW unscattered surface flux
!> @param[in,out] sw_down_blue_surf         SW blue downward surface flux
!> @param[in,out] sw_direct_blue_surf       SW blue unscattered surface flux
!> @param[in,out] sw_up_surf                SW upward surface flux
!> @param[in,out] sw_up_toa                 SW upward top-of-atmosphere flux
!> @param[in,out] sw_direct_toa             SW unscattered top-of-atmosphere flux
!> @param[in,out] sw_up_tile                SW upward tiled surface flux
!> @param[in,out] sw_up_blue_tile           SW blue upward tiled surface flux
!> @param[in]     sw_heating_rate_rts       SW heating rate
!> @param[in]     sw_down_surf_rts          SW downward surface flux
!> @param[in]     sw_direct_surf_rts        SW unscattered surface flux
!> @param[in]     sw_down_blue_surf_rts     SW blue downward surface flux
!> @param[in]     sw_direct_blue_surf_rts   SW blue unscattered surface flux
!> @param[in]     sw_up_surf_rts            SW upward surface flux
!> @param[in]     sw_up_toa_rts             SW upward top-of-atmosphere flux
!> @param[in]     sw_direct_toa_rts         SW unscattered top-of-atmosphere flux
!> @param[in]     sw_up_tile_rts            SW upward tiled surface flux
!> @param[in]     sw_up_blue_tile_rts       SW blue upward tiled surface flux
!> @param[in]     exner_in_wth              Exner pressure in wth space
!> @param[in]     cos_zenith_angle          Cosine of the stellar zenith angle
!> @param[in]     lit_fraction              Lit fraction of the timestep
!> @param[in]     cos_zenith_angle_rts      Cosine of the stellar zenith angle
!> @param[in]     lit_fraction_rts          Lit fraction of the timestep
!> @param[in]     orographic_correction_rts Orographic Correction
!> @param[in]     ndf_wth                   No. DOFs per cell for wth space
!> @param[in]     undf_wth                  No. unique of DOFs for wth space
!> @param[in]     map_wth                   Dofmap for wth space column base cell
!> @param[in]     ndf_2d                    No. of DOFs per cell for 2D space
!> @param[in]     undf_2d                   No. unique of DOFs for 2D space
!> @param[in]     map_2d                    Dofmap for 2D space column base cell
!> @param[in]     ndf_tile                  Number of DOFs per cell for tiles
!> @param[in]     undf_tile                 Number of total DOFs for tiles
!> @param[in]     map_tile                  Dofmap for tile space column base cell
subroutine sw_mts_code(nlayers, n_profile,                                     &
                   sw_heating_rate, sw_down_surf, sw_direct_surf,              &
                   sw_down_blue_surf, sw_direct_blue_surf,                     &
                   sw_up_surf, sw_up_toa, sw_direct_toa,                       &
                   sw_up_tile, sw_up_blue_tile,                                &
                   sw_heating_rate_rts, sw_down_surf_rts, sw_direct_surf_rts,  &
                   sw_down_blue_surf_rts, sw_direct_blue_surf_rts,             &
                   sw_up_surf_rts, sw_up_toa_rts, sw_direct_toa_rts,           &
                   sw_up_tile_rts, sw_up_blue_tile_rts,                        &
                   exner_in_wth,                                               &
                   cos_zenith_angle, lit_fraction,                             &
                   cos_zenith_angle_rts, lit_fraction_rts,                     &
                   orographic_correction_rts,                                  &
                   ndf_wth, undf_wth, map_wth,                                 &
                   ndf_2d, undf_2d, map_2d,                                    &
                   ndf_tile, undf_tile, map_tile)

  use radiation_config_mod, only: n_radstep, l_trans_zen_correction
  use jules_control_init_mod, only: n_surf_tile
  use socrates_init_mod, only: l_orog
  use socrates_bones, only: bones

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, n_profile
  integer(i_def), intent(in) :: ndf_wth, ndf_2d, ndf_tile
  integer(i_def), intent(in) :: undf_wth, undf_2d, undf_tile

  integer(i_def), dimension(ndf_wth, n_profile),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d, n_profile),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile, n_profile),  intent(in) :: map_tile

  real(r_def), dimension(undf_wth),  intent(inout) :: sw_heating_rate
  real(r_def), dimension(undf_2d),   intent(inout) :: sw_down_surf, &
    sw_direct_surf, &
    sw_down_blue_surf, sw_direct_blue_surf, &
    sw_up_surf, sw_up_toa, sw_direct_toa
  real(r_def), dimension(undf_tile), intent(inout) :: sw_up_tile, &
    sw_up_blue_tile

  real(r_def), dimension(undf_wth),  intent(in) :: sw_heating_rate_rts
  real(r_def), dimension(undf_2d),   intent(in) :: sw_down_surf_rts, &
    sw_direct_surf_rts, &
    sw_down_blue_surf_rts, sw_direct_blue_surf_rts, &
    sw_up_surf_rts, sw_up_toa_rts, sw_direct_toa_rts
  real(r_def), dimension(undf_tile), intent(in) :: sw_up_tile_rts, &
    sw_up_blue_tile_rts

  real(r_def), dimension(undf_wth), intent(in) :: exner_in_wth

  real(r_def), dimension(undf_2d), intent(in) :: &
    cos_zenith_angle, lit_fraction, &
    cos_zenith_angle_rts, lit_fraction_rts, &
    orographic_correction_rts

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

  if (n_radstep == 1) then
    ! Radiation timestep = model timestep
    sw_heating_rate(wth_0:wth_last)  = sw_heating_rate_rts(wth_0:wth_last)
    sw_down_surf(twod_1:twod_last)   = sw_down_surf_rts(twod_1:twod_last)
    sw_up_surf(twod_1:twod_last)     = sw_up_surf_rts(twod_1:twod_last)
    sw_up_toa(twod_1:twod_last)      = sw_up_toa_rts(twod_1:twod_last)
    sw_up_tile(tile_1:tile_last)     = sw_up_tile_rts(tile_1:tile_last)
    sw_direct_surf(twod_1:twod_last) = sw_direct_surf_rts(twod_1:twod_last)
    sw_direct_toa(twod_1:twod_last)  = sw_direct_toa_rts(twod_1:twod_last)
    sw_down_blue_surf(twod_1:twod_last) &
                                     = sw_down_blue_surf_rts(twod_1:twod_last)
    sw_direct_blue_surf(twod_1:twod_last) &
                                     = sw_direct_blue_surf_rts(twod_1:twod_last)
    sw_up_blue_tile(tile_1:tile_last) &
                                     = sw_up_blue_tile_rts(tile_1:tile_last)
  else
    ! Corrections to model timestep. The radiative fluxes have been calculated
    ! for a mean sun angle over the radiation timestep and must be converted
    ! to fluxes appropriate for the mean sun angle over the shorter model
    ! timestep.

    ! The bare "bones" of a simple radiative transfer calculation.
    ! Apply the solar zenith angle correction.
    call bones(n_profile, nlayers,                                             &
      n_layer_stride            = nlayers+1,                                   &
      n_tile                    = n_surf_tile,                                 &
      l_cos_zen_correction      = .true.,                                      &
      cos_zen_rts               = cos_zenith_angle_rts(twod_1:twod_last),      &
      lit_frac_rts              = lit_fraction_rts(twod_1:twod_last),          &
      cos_zen_mts               = cos_zenith_angle(twod_1:twod_last),          &
      lit_frac_mts              = lit_fraction(twod_1:twod_last),              &
      l_trans_zen_correction    = l_trans_zen_correction,                      &
      l_orog_corr_rts           = l_orog,                                      &
      orog_corr_rts             = orographic_correction_rts(twod_1:twod_last), &
      heating_rate_1d_rts       = sw_heating_rate_rts(wth_1:wth_last),         &
      flux_up_tile_1d_rts       = sw_up_tile_rts(tile_1:tile_last),            &
      flux_up_blue_tile_1d_rts  = sw_up_blue_tile_rts(tile_1:tile_last),       &
      flux_direct_toa_rts       = sw_direct_toa_rts(twod_1:twod_last),         &
      flux_up_toa_rts           = sw_up_toa_rts(twod_1:twod_last),             &
      flux_direct_surf_rts      = sw_direct_surf_rts(twod_1:twod_last),        &
      flux_down_surf_rts        = sw_down_surf_rts(twod_1:twod_last),          &
      flux_up_surf_rts          = sw_up_surf_rts(twod_1:twod_last),            &
      flux_direct_blue_surf_rts = sw_direct_blue_surf_rts(twod_1:twod_last),   &
      flux_down_blue_surf_rts   = sw_down_blue_surf_rts(twod_1:twod_last),     &
      heating_rate_1d_mts       = sw_heating_rate(wth_1:wth_last),             &
      flux_up_tile_1d_mts       = sw_up_tile(tile_1:tile_last),                &
      flux_up_blue_tile_1d_mts  = sw_up_blue_tile(tile_1:tile_last),           &
      flux_direct_toa_mts       = sw_direct_toa(twod_1:twod_last),             &
      flux_up_toa_mts           = sw_up_toa(twod_1:twod_last),                 &
      flux_direct_surf_mts      = sw_direct_surf(twod_1:twod_last),            &
      flux_down_surf_mts        = sw_down_surf(twod_1:twod_last),              &
      flux_up_surf_mts          = sw_up_surf(twod_1:twod_last),                &
      flux_direct_blue_surf_mts = sw_direct_blue_surf(twod_1:twod_last),       &
      flux_down_blue_surf_mts   = sw_down_blue_surf(twod_1:twod_last),         &
      l_profile_last            = .true.)
  end if

  ! Set level 0 increment such that theta increment will equal level 1
  do k=wth_0, wth_last, nlayers+1
    sw_heating_rate(k) = sw_heating_rate(k+1) &
                       * exner_in_wth(k) / exner_in_wth(k+1)
  end do

end subroutine sw_mts_code
end module sw_mts_kernel_mod
