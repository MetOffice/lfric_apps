!----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Read all the lut namelist files

module um_radaer_lut_init_mod

  use constants_mod,             only : r_um
  use aerosol_config_mod,        only : aclw_file,                 &
                                        acsw_file,                 &
                                        anlw_file,                 &
                                        answ_file,                 &
                                        crlw_file,                 &
                                        crsw_file,                 &
                                        glomap_mode,               &
                                        glomap_mode_climatology,   &
                                        glomap_mode_dust_and_clim, &
                                        glomap_mode_off,           &
                                        glomap_mode_radaer_test,   &
                                        glomap_mode_ukca,          &
                                        l_radaer,                  &
                                        prec_file
  use section_choice_config_mod, only : aerosol, &
                                        aerosol_um
  use socrates_init_mod,         only : sw_wavelength_short, &
                                        sw_wavelength_long,  &
                                        lw_wavelength_short, &
                                        lw_wavelength_long,  &
                                        n_sw_band,           &
                                        n_lw_band
  use ukca_radaer_lut,           only : ip_ukca_lut_accum, ip_ukca_lut_coarse, &
                                        ip_ukca_lut_accnarrow,                 &
                                        ip_ukca_lut_cornarrow,                 &
                                        ip_ukca_lut_supercoarse,               &
                                        ip_ukca_lut_sw, ip_ukca_lut_lw
  use ukca_radaer_read_precalc_mod, &
                                 only : ukca_radaer_read_precalc
  use um_read_radaer_lut_mod,    only : um_read_radaer_lut

  implicit none

  private
  public :: um_radaer_lut_init

contains

  !>@brief    Read all the lut namelist files.
  !>@details  Read all the lut namelist files and populate the UM to ukca_lut
  !>          structure.

  subroutine um_radaer_lut_init()

    use filenamelength_mod, only: filenamelength

    implicit none

    ! Local variables

    character(len=filenamelength) :: dummy_file

    dummy_file = 'unset'

    if ( aerosol == aerosol_um ) then

      select case (glomap_mode)
      case( glomap_mode_climatology,   &
            glomap_mode_dust_and_clim, &
            glomap_mode_radaer_test,   &
            glomap_mode_ukca )

          if ( l_radaer ) then

            call um_read_radaer_lut ( aclw_file, &
                                      ip_ukca_lut_accum, ip_ukca_lut_lw )

            call um_read_radaer_lut ( acsw_file, &
                                      ip_ukca_lut_accum, ip_ukca_lut_sw )

            call um_read_radaer_lut ( anlw_file, &
                                      ip_ukca_lut_accnarrow, ip_ukca_lut_lw )

            call um_read_radaer_lut ( answ_file, &
                                      ip_ukca_lut_accnarrow, ip_ukca_lut_sw )

            call um_read_radaer_lut ( crlw_file, &
                                      ip_ukca_lut_coarse, ip_ukca_lut_lw )

            call um_read_radaer_lut ( crsw_file, &
                                      ip_ukca_lut_coarse, ip_ukca_lut_sw )

            call um_read_radaer_lut ( dummy_file, &
                                      ip_ukca_lut_cornarrow, ip_ukca_lut_lw )

            call um_read_radaer_lut ( dummy_file, &
                                      ip_ukca_lut_cornarrow, ip_ukca_lut_sw )

            call um_read_radaer_lut ( dummy_file, &
                                      ip_ukca_lut_supercoarse, ip_ukca_lut_lw )

            call um_read_radaer_lut ( dummy_file, &
                                      ip_ukca_lut_supercoarse, ip_ukca_lut_sw )

            call ukca_radaer_read_precalc( prec_file,                      &
                                           real(sw_wavelength_short,r_um), &
                                           real(sw_wavelength_long,r_um),  &
                                           real(lw_wavelength_short,r_um), &
                                           real(lw_wavelength_long,r_um),  &
                                           n_sw_band,                      &
                                           n_lw_band )
          end if

      end select

    end if

  end subroutine um_radaer_lut_init

end module um_radaer_lut_init_mod
