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

  use ukca_api_mod,                  only: ukca_i_sussbcocdu_7mode,            &
                                           ukca_i_du_2mode,                    &
                                           ukca_i_ukca_bc_tuned,               &
                                           ukca_i_ukca_tune_bc_off

  use constants_mod,                 only: i_um

  use ukca_radaer_lfric_init_mod,    only: ukca_radaer_lfric_init

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

  if ( glomap_mode == glomap_mode_climatology ) then
    ! mode_setup is not set in the namelist for glomap_mode_climatology
    ! this is always fixed to ukca_i_sussbcocdu_7mode.
    i_mode_setup_radaer_local = ukca_i_sussbcocdu_7mode

    ! Tune BC turned off
    i_ukca_tune_bc_local      = ukca_i_ukca_tune_bc_off

    ! Dust ageing turned off
    l_dust_mp_ageing_local    = .false.

    ! sustrat turned on
    l_ukca_radaer_sustrat_local = .true.

    call ukca_radaer_lfric_init( i_mode_setup_radaer_local,                    &
                                 i_ukca_tune_bc_local,                         &
                                 l_dust_mp_ageing_local,                       &
                                 l_ukca_radaer_sustrat_local )

  else if ( ( glomap_mode == glomap_mode_dust_and_clim ) .or.                  &
            ( glomap_mode == glomap_mode_radaer_test   ) ) then
    ! dust_and_clim runs with a diffent mode_setup between ukca and radaer
    ! this is always fixed to ukca_i_sussbcocdu_7mode.
    i_mode_setup_radaer_local = ukca_i_sussbcocdu_7mode

    ! Tune BC turned off
    i_ukca_tune_bc_local      = ukca_i_ukca_tune_bc_off

    ! Dust ageing not allowed for dust only ukca
    l_dust_mp_ageing_local    = .false.

    ! sustrat turned on
    l_ukca_radaer_sustrat_local = .true.

    call ukca_radaer_lfric_init( i_mode_setup_radaer_local,                    &
                                 i_ukca_tune_bc_local,                         &
                                 l_dust_mp_ageing_local,                       &
                                 l_ukca_radaer_sustrat_local )

  else if ( glomap_mode == glomap_mode_ukca ) then
    ! UKCA and RADAER will use the same value for mode_setup

    ! Match rose-meta integers with those used in UKCA
    select case ( mode_setup )
    case ( mode_setup_SUBCOCSSDU_7mode )    
      i_mode_setup_radaer_local = ukca_i_sussbcocdu_7mode
    case ( mode_setup_DUonly_2mode )
      i_mode_setup_radaer_local = ukca_i_du_2mode
    case default
      write( log_scratch_space, '(A,I0)' )                                     &
      'Developers should include additional mode settings here: ', mode_setup
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select

    ! UKCA currently set to ukca_i_ukca_bc_tuned
    i_ukca_tune_bc_local      = ukca_i_ukca_bc_tuned

    ! Dust ageing set by namelist
    l_dust_mp_ageing_local    = l_dust_mp_ageing

    ! sustrat set by namelist
    l_ukca_radaer_sustrat_local = l_ukca_radaer_sustrat

    call ukca_radaer_lfric_init( i_mode_setup_radaer_local,                    &
                                 i_ukca_tune_bc_local,                         &
                                 l_dust_mp_ageing_local,                       &
                                 l_ukca_radaer_sustrat_local )

  end if

end subroutine um_radaer_init

end module um_radaer_init_mod

