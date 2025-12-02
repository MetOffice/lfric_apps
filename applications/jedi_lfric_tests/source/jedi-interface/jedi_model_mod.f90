!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a class than handles the non-linear model
!>
!> @details This module includes a class that handles the non-linear model
!>          time stepping. The class includes init, step and final methods.
!>          These are used by the forecast method also included within the
!>          class. In JEDI, the forecast method is defined in the OOPS base
!>          class and the init, step and final are defined in the model
!>          interface (LFRIC-JEDI). An included forecast application uses the
!>          model forecast method to propagate the state.
module jedi_model_mod

  use constants_mod,                 only : str_def
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_state_mod,                only : jedi_state_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use namelist_collection_mod,       only : namelist_collection_type
  use namelist_mod,                  only : namelist_type

  implicit none

  private

type, public :: jedi_model_type
  private

  !> The model time step duration
  type( jedi_duration_type ) :: time_step

contains

  !> Model initialiser.
  procedure, public  :: initialise

  !> Methods
  procedure, private :: model_init
  procedure, private :: model_step
  procedure, private :: model_final

  !> Run a forecast
  procedure, public  :: forecast

  !> Finalizer
  final              :: jedi_model_destructor

end type jedi_model_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_model_type
!>
!> @param [in] configuration Configuration used to setup the model class
subroutine initialise( self, configuration )

  implicit none

  class( jedi_model_type ),      intent(inout) :: self
  type( namelist_collection_type ), intent(in) :: configuration

  ! Local
  type( namelist_type ), pointer :: jedi_model_config
  character( str_def )           :: time_step_str

  ! Get config info and setup
  jedi_model_config => configuration%get_namelist('jedi_model')
  call jedi_model_config%get_value( 'time_step', time_step_str )
  call self%time_step%init(time_step_str)

end subroutine initialise


!> @brief    Initialise the model
!>
!> @param [inout] state State object to be used in the model initialise
subroutine model_init(self, state)

  implicit none

  class( jedi_model_type ), target, intent(in) :: self
  type( jedi_state_type ),  intent(inout)      :: state

end subroutine model_init

!> @brief    Step the model
!>
!> @param [inout] state State object to propagated
subroutine model_step(self, state)

  use jedi_lfric_nl_modeldb_driver_mod, only : step_nl

  implicit none

  class( jedi_model_type ), target, intent(inout) :: self
  type( jedi_state_type ),          intent(inout) :: state

  ! Local
  logical :: clock_running

  ! check the clock
  clock_running = state%modeldb%clock%tick()
  ! If the clock has finished then it will just get the
  ! data at the end of the file - this prevents that
  if ( .not. clock_running ) then
    write ( log_scratch_space, '(A)' ) &
            "Model::model_step::The LFRic clock has stopped."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  ! step gungho
  call step_nl( state%modeldb )

  ! Copy fields from model data
  call state%from_modeldb()

  ! update the state time
  call state%update_time( self%time_step )

end subroutine model_step

!> @brief    Finalise the model
!>
!> @param [inout] state State object to be used in the model finalise
subroutine model_final(self, state)

  implicit none

  class( jedi_model_type ), target, intent(in) :: self
  type( jedi_state_type ),       intent(inout) :: state

end subroutine model_final

!> @brief    Finalize the jedi_pseudo_model_type
!>
subroutine jedi_model_destructor(self)

  implicit none

  type(jedi_model_type), intent(inout) :: self

end subroutine jedi_model_destructor

!------------------------------------------------------------------------------
! OOPS defined forecast method
!------------------------------------------------------------------------------

!> @brief    Run a forecast using the model init, step and final
!>
!> @param [inout] state           The state object to propagate
!> @param [in]    forecast_length The duration of the forecast
!> @param [inout] post_processor  Post processing object
subroutine forecast(self, state, forecast_length, post_processor)

  use jedi_post_processor_mod,    only : jedi_post_processor_type

  implicit none

  class( jedi_model_type ),          intent(inout) :: self
  type( jedi_state_type ),           intent(inout) :: state
  type( jedi_duration_type ),           intent(in) :: forecast_length
  class( jedi_post_processor_type ), intent(inout) :: post_processor

  ! Local
  type( jedi_datetime_type ) :: end_time

  ! Get the end time
  end_time = state%valid_time() + forecast_length

  ! Initialize the model
  call self%model_init( state )
  ! Initialize the post processor and call first process
  call post_processor%pp_init( state )
  call post_processor%process( state )

  ! Loop until date_time_end
  do while ( end_time > state%valid_time() )
    call self%model_step( state )
    call post_processor%process( state )
  end do

  ! Finalize model and post processor
  call post_processor%pp_final( state )
  call self%model_final( state )

end subroutine forecast

end module jedi_model_mod
