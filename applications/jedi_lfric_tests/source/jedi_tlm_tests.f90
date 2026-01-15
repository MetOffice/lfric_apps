!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page jedi_tlm_tests program

!> @brief Main program for running tlm adjoint tests with jedi emulator
!>        objects. This version is used for the full linear and adjoint models.
!>
!> @details Setup and run the adjoint tests using the JEDI emulator objects.
!>          The linear state trajectory is provided via the pseudo model
!>          forecast. The jedi objects are constructed via an initialiser call
!>          and the forecasts are handled by the linear model object.
!>
!>          The standard dot product adjoint tests is employed here which
!>          relies on the following identity for the inner product denoted with
!>          angled braces <>:
!>
!>          <Mx,Mx> == <AMx,x>
!>
!>          where M is the linear model forecast and A is the adjoint model forecast.
!>
!>          This is true for any perturbation state vector x and so in
!>          the test, a random one is used. x is called inc in the code.
!>
!>          Note this test includes scaling of the prognostic fields so that individual
!>          fields do not dominate the total inner product. Otherwise, the test is not fair.

! Note: This program file represents generic JEDI code and so it should not be
!       edited. If you need to make changes at the program level then please
!       contact darth@metofice.gov.uk for advice.

program jedi_tlm_tests

  use cli_mod,                      only : parse_command_line
  use constants_mod,                only : PRECISION_REAL, i_def, str_def, r_def
  use field_collection_mod,         only : field_collection_type
  use log_mod,                      only : log_event, log_scratch_space, &
                                           LOG_LEVEL_ALWAYS, LOG_LEVEL_ERROR, &
                                           LOG_LEVEL_INFO
  use namelist_collection_mod,      only : namelist_collection_type
  use namelist_mod,                 only : namelist_type

  ! Jedi emulator objects
  use jedi_lfric_duration_mod,      only : jedi_duration_type
  use jedi_run_mod,                 only : jedi_run_type
  use jedi_geometry_mod,            only : jedi_geometry_type
  use jedi_state_mod,               only : jedi_state_type
  use jedi_increment_mod,           only : jedi_increment_type
  use jedi_pseudo_model_mod,        only : jedi_pseudo_model_type
  use jedi_linear_model_mod,        only : jedi_linear_model_type
  use jedi_post_processor_traj_mod, only : jedi_post_processor_traj_type

  implicit none

  ! Emulator objects
  type( jedi_geometry_type )            :: geometry
  type( jedi_state_type )               :: state
  type( jedi_increment_type )           :: inc
  type( jedi_increment_type )           :: inc_initial
  type( jedi_pseudo_model_type )        :: pseudo_model
  type( jedi_linear_model_type )        :: linear_model
  type( jedi_run_type )                 :: run
  type( jedi_post_processor_traj_type ) :: pp_traj

  ! Local
  type( namelist_collection_type ), pointer :: configuration
  character(:),                 allocatable :: filename
  integer( kind=i_def )                     :: model_communicator
  type( jedi_duration_type )                :: forecast_length
  type( namelist_type ),            pointer :: jedi_lfric_settings_config
  character( str_def )                      :: forecast_length_str
  real( kind=r_def )                        :: dot_product_1
  real( kind=r_def )                        :: dot_product_2
  real( kind=r_def ),             parameter :: absolute_tolerance = 1.0E-4_r_def
  real( kind=r_def )                        :: machine_tolerance
  real( kind=r_def )                        :: absolute_diff
  real( kind=r_def )                        :: relative_diff

  character(*), parameter :: program_name = "jedi_tlm_tests"

  ! Infrastructure config
  call parse_command_line( filename )

  ! Run object - handles initialization and finalization of required
  ! infrastructure. Initialize external libraries such as XIOS
  call run%initialise( program_name, model_communicator )

  ! Ensemble applications would split the communicator here

  ! Initialize LFRic infrastructure
  call run%initialise_infrastructure( filename, model_communicator )

  call log_event( 'Running ' // program_name // ' ...', LOG_LEVEL_ALWAYS )
  write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  ! Get the configuration
  configuration => run%get_configuration()

  ! Get the forecast length
  jedi_lfric_settings_config => configuration%get_namelist('jedi_lfric_settings')
  call jedi_lfric_settings_config%get_value( 'forecast_length', forecast_length_str )
  call forecast_length%init(forecast_length_str)

  ! Create geometry
  call geometry%initialise( model_communicator, configuration )

  ! Create state
  call state%initialise( geometry, configuration )

  ! Create linear model
  call linear_model%initialise( geometry, filename )

  ! Initialise trajectory post processor with instance of linear_model
  call pp_traj%initialise( linear_model )

  ! Create non-linear model
  call pseudo_model%initialise( configuration )

  ! Run non-linear model forecast to populate the trajectory object
  call pseudo_model%forecast( state, forecast_length, pp_traj )

  ! ---- Perform the adjoint test

  ! Create inc_initial and randomise
  call inc_initial%initialise( geometry, configuration )
  call inc_initial%random()

  ! Check the norm is not zero
  if (inc_initial%norm() <= 0.0_r_def) then
    call log_event("inc_initial norm not > 0.0", LOG_LEVEL_ERROR)
  end if

  ! Create inc via copy constructor using inc_initial
  call inc%initialise( inc_initial )

  ! Propagate via TL model
  call linear_model%forecastTL( inc, forecast_length )

  ! Check the norm is not zero
  if (inc%norm() <= 0.0_r_def) then
    call log_event("inc norm not > 0.0", LOG_LEVEL_ERROR)
  end if

  ! Compute <Mx,Mx>
  dot_product_1 = real(inc%scaled_dot_product_with_itself(), r_def)

  ! Propagate via AD model
  call linear_model%forecastAD( inc, forecast_length )
  ! Check the norm is not zero
  if (inc%norm() <= 0.0_r_def) then
    call log_event("inc norm not > 0.0", LOG_LEVEL_ERROR)
  end if

  ! Compute <AMx,x>
  dot_product_2 = real(inc%dot_product_with(inc_initial), r_def)

  ! The two dot products should be nearly identical. The tolerance is included
  ! due to differences in order of operations and solver non-convergence.
  absolute_diff = abs( dot_product_1 - dot_product_2 )
  machine_tolerance = spacing( max( abs( dot_product_1 ), abs( dot_product_2 ) ) )
  relative_diff = absolute_diff / machine_tolerance
  if (absolute_diff > absolute_tolerance ) then
    call run%finalise_timers()  ! We still want timing info even if the test fails
    write( log_scratch_space, * ) "Adjoint test FAILED", &
      dot_product_1, dot_product_2, absolute_diff, relative_diff
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  else
    write( log_scratch_space, * ) "Adjoint test PASSED", &
      dot_product_1, dot_product_2, absolute_diff, relative_diff
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  endif

  call log_event( 'Finalising ' // program_name // ' ...', LOG_LEVEL_ALWAYS )

  call run%finalise()

end program jedi_tlm_tests
