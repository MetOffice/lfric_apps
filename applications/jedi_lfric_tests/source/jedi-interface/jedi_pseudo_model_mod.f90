!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the JEDI pseudo model emulator class.
!>
!> @details This module holds a JEDI pseudo model emulator where the
!>          time-stepping uses the read_file state method. An included forecast
!>          application uses the model init, step and final to run a forecast.
!>
module jedi_pseudo_model_mod

  use constants_mod,                 only : i_def, str_def
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_state_mod,                only : jedi_state_type
  use log_mod,                       only : log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_ERROR
  use namelist_collection_mod,       only : namelist_collection_type
  use namelist_mod,                  only : namelist_type

  implicit none

  private

type, public :: jedi_pseudo_model_type
  private

  !> A date_time list to read
  type( jedi_datetime_type ), allocatable  :: state_times(:)
  !> Index for the current state in date_time_states
  integer( kind=i_def )                    :: current_state
  !> The number of states
  integer( kind=i_def )                    :: n_states

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
  final              :: jedi_pseudo_model_destructor

end type jedi_pseudo_model_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_pseudo_model_type
!>
!> @param [in] config Configuration used to setup the model class
subroutine initialise( self, configuration )

  implicit none

  class( jedi_pseudo_model_type ),  intent(inout) :: self
  type( namelist_collection_type ), intent(in)    :: configuration

  ! Local
  type( namelist_type ), pointer :: jedi_model_config
  character( str_def )           :: initial_time
  character( str_def )           :: time_step_str
  integer( i_def )               :: number_of_steps
  integer( i_def )               :: i
  type( jedi_datetime_type )     :: next_datetime
  type( jedi_duration_type )     :: time_step

  ! Setup the pseudo model
  self%current_state = 1_i_def

  ! Get config info and setup
  jedi_model_config => configuration%get_namelist('jedi_pseudo_model')

  call jedi_model_config%get_value( 'number_of_steps', number_of_steps )
  self%n_states = number_of_steps
  allocate( self%state_times(self%n_states) )

  call jedi_model_config%get_value( 'time_step', time_step_str )
  call time_step%init( time_step_str )

  call jedi_model_config%get_value( 'initial_time', initial_time )
  ! Initialise datetime states 1 - number_of_steps time steps after
  ! lfric calendar_start namelist variable time
  call next_datetime%init( initial_time )
  do i = 1, self%n_states
    next_datetime = next_datetime + time_step
    self%state_times(i) = next_datetime
  end do

end subroutine initialise

!> @brief    Initialise the model
!>
!> @param [inout] state State object to be used in the model initialise
subroutine model_init( self, state )

  implicit none

  class( jedi_pseudo_model_type ), intent(in)    :: self
  type( jedi_state_type ),         intent(inout) :: state

end subroutine model_init

!> @brief    Step the model
!>
!> @param [inout] state State object to propagated
subroutine model_step( self, state )

  implicit none

  class( jedi_pseudo_model_type ), intent(inout) :: self
  type( jedi_state_type ),         intent(inout) :: state

  ! Local
  type( jedi_datetime_type ) :: state_time

  if ( self%current_state > self%n_states ) then
    write ( log_scratch_space, '(A)' ) "self%current_state>self%n_states."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  state_time = self%state_times( self%current_state )
  call state%read_file( state_time )

  ! Iterate the current state
  self%current_state = self%current_state + 1_i_def

end subroutine model_step

!> @brief    Finalise the model
!>
!> @param [inout] state State object to be used in the model finalise
subroutine model_final( self, state )

  implicit none

  class( jedi_pseudo_model_type ), intent(in)    :: self
  type( jedi_state_type ),         intent(inout) :: state

end subroutine model_final

!> @brief    Finalize the jedi_pseudo_model_type
!>
subroutine jedi_pseudo_model_destructor( self )

  implicit none

  type( jedi_pseudo_model_type ), intent(inout) :: self

  if ( allocated(self%state_times) ) deallocate(self%state_times)

end subroutine jedi_pseudo_model_destructor

!------------------------------------------------------------------------------
! OOPS defined forecast method
!------------------------------------------------------------------------------

!> @brief    Run a forecast using the model init, step and final
!>
!> @param [inout] state           The state object to propagate
!> @param [in]    forecast_length The duration of the forecast
!> @param [inout] post_processor  Post processing object
subroutine forecast( self, state, forecast_length, post_processor )

  use jedi_post_processor_mod,    only : jedi_post_processor_type

  implicit none

  class( jedi_pseudo_model_type ),   intent(inout) :: self
  type( jedi_state_type ),           intent(inout) :: state
  type( jedi_duration_type ),           intent(in) :: forecast_length
  class( jedi_post_processor_type ), intent(inout) :: post_processor

  ! Local
  type( jedi_datetime_type ) :: end_time

  ! End time
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

end module jedi_pseudo_model_mod
