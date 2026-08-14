!----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief RADAER initialisation subroutine for science configuration

module um_radaer_init_mod

  ! LFRic namelists which have been read
  use aerosol_config_mod,            only: glomap_mode,                        &
                                           glomap_mode_climatology,            &
                                           glomap_mode_dust_and_clim,          &
                                           glomap_mode_off,                    &
                                           glomap_mode_radaer_test,            &
                                           glomap_mode_ukca,                   &
                                           mode_setup,                         &
                                           mode_setup_SUBCOCSSDU_7mode,        &
                                           mode_setup_DUonly_2mode,            &
                                           l_dust_mp_ageing,                   &
                                           l_ukca_radaer_sustrat

  use ukca_mode_setup,               only: i_ukca_bc_tuned

  use constants_mod,                 only: i_um

  use ukca_radaer_lfric_init_mod,    only: ukca_radaer_lfric_init

  use ukca_config_specification_mod, only: i_sussbcocdu_7mode,                 &
                                           i_du_2mode

  use log_mod,                       only: log_event,                          &
                                           log_scratch_space,                  &
                                           LOG_LEVEL_ERROR

  implicit none

  private
  public :: um_radaer_init

contains

subroutine um_radaer_init()

  implicit none

  ! Some options have a different value of mode_setup
  ! between ukca and radaer.
  ! For example, we use prognostic dust (6) with other aerosol components (8)
  ! in operational NWP.
  integer(i_um) :: i_mode_setup_radaer_local

  integer(i_um) :: i_ukca_tune_bc_local
  
  logical :: l_dust_mp_ageing_local

  logical :: l_ukca_radaer_sustrat_local

  integer, parameter :: i_radaer_mode_setup_eight = 8

  integer, parameter :: i_ukca_bc_tuned_zero = 0

  if ( glomap_mode == glomap_mode_climatology ) then
    ! mode_setup is not set in the namelist for glomap_mode_climatology
    ! this is always fixed to eight.
    i_mode_setup_radaer_local = i_radaer_mode_setup_eight

    ! Tune BC turned off
    i_ukca_tune_bc_local      = i_ukca_bc_tuned_zero

    ! Dust ageing turned off
    l_dust_mp_ageing_local    = .false.

    ! sustrat turned off
    l_ukca_radaer_sustrat_local = .false.

  else if ( glomap_mode == glomap_mode_dust_and_clim ) then
    ! dust_and_clim runs with a diffent mode_setup between ukca and radaer
    ! this is always fixed to eight.
    i_mode_setup_radaer_local = i_radaer_mode_setup_eight

    ! Tune BC turned off
    i_ukca_tune_bc_local      = i_ukca_bc_tuned_zero

    ! Dust ageing not allowed for dust only ukca
    l_dust_mp_ageing_local    = .false.

    ! sustrat turned off
    l_ukca_radaer_sustrat_local = .false.
 
  else if ( glomap_mode == glomap_mode_radaer_test ) then
    ! This was developed for aqua planet runs and may be redundant
    ! For now fix this to eight.
    i_mode_setup_radaer_local = i_radaer_mode_setup_eight

    ! Tune BC turned off
    i_ukca_tune_bc_local      = i_ukca_bc_tuned_zero

    ! Dust ageing turned off
    l_dust_mp_ageing_local    = .false.

    ! sustrat turned off
    l_ukca_radaer_sustrat_local = .false.

  else if ( glomap_mode == glomap_mode_ukca ) then
    ! UKCA and RADAER will use the same value for mode_setup

    ! Match rose-meta integers with those used in UKCA
    select case ( mode_setup )
    case ( mode_setup_SUBCOCSSDU_7mode )    
      i_mode_setup_radaer_local = i_sussbcocdu_7mode
    case ( mode_setup_DUonly_2mode )
      i_mode_setup_radaer_local = i_du_2mode
    case default
      write( log_scratch_space, '(A,I0)' )                                     &
      'Developers should include additional mode settings here: ', mode_setup
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select

    ! UKCA currently set to i_ukca_bc_tuned
    i_ukca_tune_bc_local      = i_ukca_bc_tuned

    ! Dust ageing set by namelist
    l_dust_mp_ageing_local    = l_dust_mp_ageing

    ! sustrat set by namelist
    l_ukca_radaer_sustrat_local = l_ukca_radaer_sustrat

  end if

  call ukca_radaer_lfric_init( i_mode_setup_radaer_local,                      &
                               i_ukca_tune_bc_local,                           &
                               l_dust_mp_ageing_local,                         &
                               l_ukca_radaer_sustrat_local )

end subroutine um_radaer_init

end module um_radaer_init_mod

