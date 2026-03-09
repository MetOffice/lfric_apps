!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page jedi_id_tlm_tests program

!> @brief Main program for running tlm adjoint tests with jedi emulator
!>        objects. This version is for the identity linear and adjoint models.
!>
!> @details Setup and run the adjoint tests using the JEDI emulator objects.
!>          The linear state trajectory is provided via the pseudo model
!>          forecast. The jedi objects are constructed via an initialiser call
!>          and the forecasts are handled by the linear model object.
!>
!>          The standard dot product adjoint test is employed here which
!>          relies on the following identity for the inner product denoted with
!>          angled braces <>:
!>
!>          <x,My> == <M^{T}x,y>
!>
!>          where M is the linear model forecast and M^{T} is the adjoint
!>          (transpose) model forecast. x and y represent the models
!>          initial perturbation state . Using arbitrary notation for a perturbation
!>          state "I" and defining the following identities:
!>
!>                I22: x
!>                I21: My
!>          and
!>                I21: M^{T}x
!>                I11: y
!>
!>          Then we can test for the correctness of the adjoint by computing
!>          the inner product and checking they are equal within machine
!>          tolerance:
!>
!>            <I11,I21> == <I12,I22>
!>
!>          This is true for any perturbation state vector x and y and so in
!>          the test, different random vector are computed (I11 and I22).
!>

! Note: This program file represents generic OOPS code and so it should not be
!       edited. If you need to make changes at the program level then please
!       contact darth@metofice.gov.uk for advice.
program jedi_id_tlm_tests

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
  use jedi_id_linear_model_mod,     only : jedi_id_linear_model_type
  use jedi_post_processor_traj_mod, only : jedi_post_processor_traj_type

  implicit none

  ! Emulator objects
  type( jedi_geometry_type )            :: geometry
  type( jedi_state_type )               :: state
  type( jedi_increment_type )           :: increment_11
  type( jedi_increment_type )           :: increment_12
  type( jedi_increment_type )           :: increment_21
  type( jedi_increment_type )           :: increment_22
  type( jedi_pseudo_model_type )        :: nonlinear_model
  type( jedi_id_linear_model_type )     :: linear_model
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
  real( kind=r_def),              parameter :: overall_tolerance = 1500.0_r_def
  real( kind=r_def )                        :: machine_tolerance
  real( kind=r_def )                        :: absolute_diff
  real( kind=r_def )                        :: relative_diff

  character(*), parameter :: program_name = "jedi_id_tlm_tests"

  ! Infrastructure configuration
  call parse_command_line( filename )

  ! Run object - handles initialization and finalization of required
  ! infrastructure. Initialize external libraries such as XIOS
  call run%initialise( program_name, model_communicator )

  ! Ensemble applications would split the communicator here

  ! Initialize LFRic infrastructure
  call run%initialise_infrastructure( filename, model_communicator )

  call log_event( 'Running ' // program_name // ' ...', LOG_LEVEL_ALWAYS )
  write(log_scratch_space,'(A)')                          &
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

  ! Initialise trajectory post processor with instance of linear_model. This is
  ! an object which post-processes a nonlinear model forecast to store a
  ! trajectory/linearisation states for the linear model
  call pp_traj%initialise( linear_model )

  ! Initialise non-linear model
  call nonlinear_model%initialise( configuration )

  ! Run non-linear model forecast to populate the trajectory object
  call nonlinear_model%forecast( state, forecast_length, pp_traj )

  ! ---- Perform the adjoint test

  ! I11 (= y). Create I11 and randomise.
  call increment_11%initialise( geometry, configuration )
  call increment_11%random()

  ! Check the norm is not zero
  if (increment_11%norm() <= 0.0_r_def) then
    call log_event("increment_11 norm not > 0.0", LOG_LEVEL_ERROR)
  end if

  ! I12 (= My). Create via copy constructor using I11 (=y). Propagate via
  !                TL forecast.
  call increment_12%initialise( increment_11 )

  ! Propagate via TL model to obtain I12 (=My).
  call linear_model%forecastTL( increment_12, forecast_length )

  ! Check the norm is not zero
  if (increment_12%norm() <= 0.0_r_def) then
    call log_event("increment_12 norm not > 0.0", LOG_LEVEL_ERROR)
  end if

  ! I22 (= x). Create and randomise. Copy construct using I12 to ensure the
  !               valid time = t_start + forecast_length.
  call increment_22%initialise( increment_12 )
  call increment_22%random()

  ! Check the norm is not zero
  if (increment_22%norm() <= 0.0_r_def) then
    call log_event("increment_22 norm not > 0.0", LOG_LEVEL_ERROR)
  end if

  ! I21 (= M^{T}x). Create via copy constructor using I11(=x) and propagate
  !                    via AD forecast.
  call increment_21%initialise( increment_22 )

  ! Propagate via TL model to obtain I21 (=M^{T}x).
  call linear_model%forecastAD( increment_21, forecast_length )

  ! Check the norm is not zero
  if (increment_21%norm() <= 0.0_r_def) then
    call log_event("increment_21 norm not > 0.0", LOG_LEVEL_ERROR)
  end if

  ! Check the initail norms are different
  if ( increment_22%norm() == increment_11%norm() ) then
    call log_event("norms should not be identical", LOG_LEVEL_ERROR)
  end if

  ! Apply the Adjoint dot product test: <I11,I21> == <I12,I22>

  ! Evaluate dot product of increments I11 and I21 <I11,I21>
  dot_product_1 = real(increment_11%dot_product_with(increment_21), r_def)

  ! Evaluate dot product of increments I12 and I22 <I12,I22>
  dot_product_2 = real(increment_12%dot_product_with(increment_22), r_def)

  ! <I11,I21> == <I12,I22>
  !    The two dot products should be nearly identical. The tolerance is
  !    included due to order of operation differences computed as the smallest
  !    delta between two numbers of a given type using the larger dot_product
  absolute_diff = abs( dot_product_1 - dot_product_2 )
  machine_tolerance = spacing( max( abs( dot_product_1 ), abs( dot_product_2 ) ) )
  relative_diff = absolute_diff / machine_tolerance
  if (absolute_diff > overall_tolerance * machine_tolerance ) then
    call run%finalise_timers()  ! We still want timing info even if the test fails
    write( log_scratch_space, * ) "Adjoint test FAILED", &
      dot_product_1, dot_product_2, absolute_diff, relative_diff
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  else
    write( log_scratch_space, * ) "Adjoint test PASSED", &
      dot_product_1, dot_product_2, absolute_diff, relative_diff
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end if

  call log_event( 'Finalising ' // program_name // ' ...', LOG_LEVEL_ALWAYS )

  call run%finalise()

end program jedi_id_tlm_tests
