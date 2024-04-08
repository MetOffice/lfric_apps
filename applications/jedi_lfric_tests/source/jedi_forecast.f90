!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page jedi_forecast program

!> @brief Main program for running fake model forecast with jedi emulator
!>        objects.

!> @details Setup and run a fake model forecast using the jedi emulator
!>          objects. The jedi objects are constructed via an initialiser call
!>          and the forecast is handled by the model object.
!>

! Note: This program file represents generic OOPS code and so it should not be
!       edited. If you need to make changes at the program level then please
!       contact darth@metofice.gov.uk for advice.
program jedi_forecast

  use constants_mod,           only : PRECISION_REAL, i_def, i_timestep
  use log_mod,                 only : log_event, log_scratch_space, &
                                      LOG_LEVEL_ALWAYS
  use field_collection_mod,    only : field_collection_type

  ! Data types and methods to get/store configurations
  use jedi_state_config_mod,    only : jedi_state_config_type
  use jedi_geometry_config_mod, only : jedi_geometry_config_type
  use cli_mod,                  only : get_initial_filename

  ! Jedi emulator objects
  use jedi_checksum_mod,             only : output_checksum
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_run_mod,                  only : jedi_run_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use jedi_state_mod,                only : jedi_state_type
  use jedi_model_mod,                only : jedi_model_type
  use jedi_post_processor_empty_mod, only : jedi_post_processor_empty_type

  implicit none

  ! Emulator objects
  type( jedi_geometry_type )             :: jedi_geometry
  type( jedi_state_type )                :: jedi_state
  type( jedi_model_type )                :: jedi_model
  type( jedi_run_type )                  :: jedi_run
  type( jedi_post_processor_empty_type ) :: jedi_pp_empty

  ! Emulator object configs
  type( jedi_state_config_type )    :: jedi_state_config
  type( jedi_geometry_config_type ) :: jedi_geometry_config
  type( jedi_duration_type )        :: forecast_length
  type( jedi_duration_type )        :: model_time_step

  ! Local
  character(:), allocatable      :: filename
  integer( i_def )               :: model_communicator

  character(*), parameter        :: program_name = "jedi_forecast"

  type( field_collection_type ), pointer :: depository => null()

  call log_event( 'Running ' // program_name // ' ...', LOG_LEVEL_ALWAYS )
  write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  ! Infrastructure config
  call get_initial_filename( filename )

  ! Run object - handles initialization and finalization of required
  ! infrastructure. Initialize external libraries such as XIOS
  call jedi_run%initialise( program_name, model_communicator )

  ! Ensemble applications would split the communicator here

  ! Initialize LFRic infrastructure
  call jedi_run%initialise_infrastructure( filename, model_communicator )

  ! Configs for for the jedi emulator objects
  ! State config
  call jedi_state_config%initialise( use_pseudo_model = .false. )

  ! Geometry config
  call jedi_geometry_config%initialise( filename )

  ! Model config - model time step duration
  call model_time_step%init( 'P0DT1H0M0S' )

  ! Forecast config - forcast duration
  call forecast_length%init( 'P0DT6H0M0S' )

  ! Create geometry
  call jedi_geometry%initialise( model_communicator, jedi_geometry_config )

  ! Create state
  call jedi_state%initialise( jedi_geometry, jedi_state_config )

  ! Create non-linear model
  call jedi_model%initialise( model_time_step )

  ! Run non-linear model forecast
  call jedi_model%forecast( jedi_state, forecast_length, jedi_pp_empty )

  call log_event( 'Finalising ' // program_name // ' ...', LOG_LEVEL_ALWAYS )
  ! To provide KGO
  depository => jedi_state%model_data%get_field_collection("depository")
  call output_checksum( program_name, depository )

end program jedi_forecast
