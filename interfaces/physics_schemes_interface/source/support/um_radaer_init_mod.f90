!----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief RADAER initialisation subroutine for science configuration

module um_radaer_init_mod

  ! LFRic namelists which have been read
  use aerosol_config_mod,         only: glomap_mode,                           &
                                        glomap_mode_climatology,               &
                                        glomap_mode_dust_and_clim,             &
                                        glomap_mode_off,                       &
                                        glomap_mode_radaer_test,               &
                                        glomap_mode_ukca,                      &
                                        i_mode_setup

  use constants_mod,              only: i_um

  use ukca_radaer_lfric_init_mod, only: ukca_radaer_lfric_init

  implicit none

  private
  public :: um_radaer_init

contains

subroutine um_radaer_init()

  implicit none

  ! Some options have a different value of i_mode_setup
  ! between ukca and radaer.
  ! For example, we use prognostic dust (6) with other aerosol components (8)
  ! in operational NWP.
  integer(i_um) :: i_mode_setup_radaer_local

  integer, parameter :: i_radaer_mode_setup_eight = 8

  if ( glomap_mode == glomap_mode_climatology ) then
    ! i_mode_setup is not set in the namelist for glomap_mode_climatology
    ! this is always fixed to eight.
    i_mode_setup_radaer_local = i_radaer_mode_setup_eight
 
  else if ( glomap_mode == glomap_mode_dust_and_clim ) then
    ! dust_and_clim runs with a diffent i_mode_setup between ukca and radaer
    ! this is always fixed to eight.
    i_mode_setup_radaer_local = i_radaer_mode_setup_eight
 
  else if ( glomap_mode == glomap_mode_radaer_test ) then
    ! This was developed for aqua planet runs and may be redundant
    ! For now fix this to eight.
    i_mode_setup_radaer_local = i_radaer_mode_setup_eight

  else if ( glomap_mode == glomap_mode_ukca ) then
    ! UKCA and RADAER will use the same value for i_mode_setup
    i_mode_setup_radaer_local = i_mode_setup

  end if

  call ukca_radaer_lfric_init( i_mode_setup_radaer_local )

end subroutine um_radaer_init

end module um_radaer_init_mod

