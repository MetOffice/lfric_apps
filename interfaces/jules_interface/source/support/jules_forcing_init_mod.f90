!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Initialisation for running JULES driven by site forcing data
module jules_forcing_init_mod

  use constants_mod, only: r_def
  use log_mod,       only: log_event, LOG_LEVEL_INFO

  implicit none

  private
  public :: jules_forcing_init, visible_fraction, diffuse_fraction

  !> Fraction of the downward SW flux that is visible (photosynthetically
  !> active), as used for both albedos and PAR in standalone JULES
  real(kind=r_def), parameter :: visible_fraction = 0.5_r_def

  !> Fraction of the downward SW flux that is diffuse, consistent with the
  !> value of jules_vars%diff_frac set in the JULES kernels
  real(kind=r_def), parameter :: diffuse_fraction = 0.4_r_def

contains

  !> @brief Sets up the SW bands used for the JULES tile albedos
  !> @details Without a radiation scheme, the SW tile albedos are calculated
  !>          for a single broadband SW band, as in standalone JULES. The band
  !>          albedos are weighted equally between visible and near-infrared,
  !>          matching the default albedo weights (wght_alb) of standalone
  !>          JULES. This sets the band information otherwise read by
  !>          socrates_init.
  subroutine jules_forcing_init()

    use socrates_init_mod, only: n_sw_band, sw_wavelength_short, &
                                 sw_wavelength_long, sw_weight_blue

    implicit none

    call log_event( 'JULES site forcing: using a broadband SW albedo', &
                    LOG_LEVEL_INFO )

    n_sw_band = 1

    if (allocated(sw_wavelength_short)) deallocate(sw_wavelength_short)
    if (allocated(sw_wavelength_long))  deallocate(sw_wavelength_long)
    if (allocated(sw_weight_blue))      deallocate(sw_weight_blue)
    allocate( sw_wavelength_short(n_sw_band), &
              sw_wavelength_long(n_sw_band),  &
              sw_weight_blue(n_sw_band) )

    ! Negative band limits signify broadband to JULES (spec_band_bb in
    ! standalone JULES), and the visible fraction weights the band albedo
    sw_wavelength_short = -1.0_r_def
    sw_wavelength_long  = -1.0_r_def
    sw_weight_blue      = visible_fraction

  end subroutine jules_forcing_init

end module jules_forcing_init_mod
