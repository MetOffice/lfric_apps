!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a class than handles the linear model (tlm)
!>
!> @details This module includes a class that handles the linear models time
!>          stepping. The class includes init, step and final methods for the
!>          Tangent Linear (TL) and Adjoint (AD). These methods are used by the
!>          forecastTL and forecastAD methods also included within the class.
!>          The forecast methods require a linear-state produced by pre-running
!>          the non-linear model and storing the result in a trajectory
!>          object. The set_trajectory method is included to provide the means
!>          to create and populate the linear state fields. In JEDI, the
!>          forecast* methods are defined in the OOPS base class and the init,
!>          step and final are defined in the model interface (LFRIC-JEDI). An
!>          included forecast application (jedi_tlm_forecast_tl) uses the model
!>          forecastTL method to propagate the increment.
!>
module jedi_linear_model_mod

  use constants_mod,                 only : str_def
  use field_collection_mod,          only : field_collection_type
  use field_array_mod,               only : field_array_type
  use gungho_modeldb_mod,            only : modeldb_type
  use jedi_lfric_moist_fields_mod,   only : update_ls_moist_fields, &
                                            init_moist_fields
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_lfric_wind_fields_mod,    only : create_scalar_winds, &
                                            setup_vector_wind
  use jedi_geometry_mod,             only : jedi_geometry_type
  use jedi_state_mod,                only : jedi_state_type
  use jedi_increment_mod,            only : jedi_increment_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use linear_state_trajectory_mod,   only : linear_state_trajectory_type

  implicit none

  private

type, public :: jedi_linear_model_type
  private

  !> The model time step duration
  type ( jedi_duration_type )          :: time_step

  !> Trajectory of linear states obtained by running the non-linear model
  type( linear_state_trajectory_type ) :: linear_state_trajectory

  !> Modeldb that stores the model fields to propagate
  !> @todo: Required public for checksum but need to move to atlas checksum
  !>        so make it private when that work is done.
  type(modeldb_type), public           :: modeldb

contains

  !> Model initialiser.
  procedure, public  :: initialise

  !> Methods
  procedure, public  :: set_trajectory

  procedure, private :: model_initTL
  procedure, private :: model_stepTL
  procedure, private :: model_finalTL

  !> @todo: When Adjoint is available, add:
  !>        model_initAD,
  !>        model_stepAD,
  !>        model_finalAD

  !> Run a TL and AD forecasts
  procedure, public  :: forecastTL
  !> @todo: When Adjoint is available, add:
  !>        forecastAD

  !> Finalizer
  final              :: jedi_linear_model_destructor

end type jedi_linear_model_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_linear_model_type
!>
!> @param [in] jedi_geometry A JEDI geometry object
!> @param [in] config        The linear model configuration
subroutine initialise( self, jedi_geometry, config )

  use jedi_linear_model_config_mod,         only : jedi_linear_model_config_type
  use jedi_lfric_linear_modeldb_driver_mod, only : initialise_modeldb
  use jedi_lfric_timestep_mod,              only : get_configuration_timestep

  implicit none

  class( jedi_linear_model_type ),    intent(inout) :: self
  type( jedi_geometry_type ),            intent(in) :: jedi_geometry
  type( jedi_linear_model_config_type ), intent(in) :: config

  ! 1. Setup modeldb

  ! 1.1 Initialise the modeldb
  call initialise_modeldb( "linear modeldb", config%config_filename, &
                            jedi_geometry%get_mpi_comm(), self%modeldb )

  ! 1.2 Add scalar winds that link the Atlas fields. These are used to
  ! perform interpolation between horizontally cell-centred and
  ! edge based W2 winds.
  call create_scalar_winds( self%modeldb, jedi_geometry%get_mesh() )

  ! 2. Setup time
  self%time_step = get_configuration_timestep( self%modeldb%configuration )

  ! 3. Setup trajactory
  call self%linear_state_trajectory%initialise( config%forecast_length, &
                                                self%time_step )

end subroutine initialise

!> @brief    Set an instance of the trajectory
!>
!> @param [in] jedi_state The state to add to the trajectory
subroutine set_trajectory( self, jedi_state )

  use jedi_lfric_linear_fields_mod,  only : variable_names, &
                                            create_linear_fields

  implicit none

  class( jedi_linear_model_type ),   intent(inout) :: self
  type( jedi_state_type ),           intent(inout) :: jedi_state

  ! Local
  type( field_collection_type ) :: next_linear_state

  ! Create field collection that contains the linear state fields
  ! without "ls_" prepended.
  call create_linear_fields(jedi_state%geometry%get_mesh(), next_linear_state)

  ! Copy data from the input state into next_linear_state
  call jedi_state%get_to_field_collection( variable_names, &
                                           next_linear_state )

  ! Create W2 wind, interpolate from scaler winds (W3/Wtheta) then
  ! remove the scaler winds
  call setup_vector_wind(next_linear_state)

  ! Add the new linear state to the trajectory (prepends fieldnames with "ls_")
  call self%linear_state_trajectory%add_linear_state( &
                                              jedi_state%valid_time(), &
                                              next_linear_state )

end subroutine set_trajectory

!> @brief    Initialise the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model initialise
subroutine model_initTL(self, increment)

  implicit none

  class( jedi_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),     intent(inout) :: increment

  ! Local
  type( field_collection_type),  pointer :: moisture_fields
  type( field_collection_type),  pointer :: prognostic_fields
  type( field_array_type ),      pointer :: mr_array
  type( field_array_type ),      pointer :: moist_dyn_array

  nullify(moisture_fields, mr_array, moist_dyn_array)

  ! Update the LFRic modeldb pertubation fields

  ! Update the prognostic fields: copy from Atlas
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  call increment%get_to_model_prognostics(prognostic_fields)

  ! Update the missing mixing ratio and moist_dynamics fields. These fields are
  ! computed analytically as outlined in jedi_lfric_linear_fields_mod
  moisture_fields => self%modeldb%fields%get_field_collection("moisture_fields")
  call moisture_fields%get_field("mr", mr_array)
  call moisture_fields%get_field("moist_dyn", moist_dyn_array)
  call init_moist_fields( mr_array%bundle, &
                          moist_dyn_array%bundle )

end subroutine model_initTL

!> @brief    Step the Tangent Linear model
!>
!> @param [inout] increment Increment object to be propagated
subroutine model_stepTL(self, increment)

  use jedi_lfric_linear_modeldb_driver_mod, only: step

  implicit none

  class( jedi_linear_model_type ), target, intent(inout) :: self
  type( jedi_increment_type ),             intent(inout) :: increment

  ! Local
  type( field_collection_type ), pointer :: depository
  type( field_collection_type),  pointer :: moisture_fields
  type( field_collection_type),  pointer :: prognostic_fields
  type( field_array_type ),      pointer :: mr_array
  type( field_array_type ),      pointer :: moist_dyn_array

  nullify(depository, moisture_fields, mr_array, moist_dyn_array)

  ! 1. Update the LFRic modeldb linear state fields

  ! 1.1 Copy from the trajectory into the model_data
  depository => self%modeldb%fields%get_field_collection("depository")
  call self%linear_state_trajectory%get_linear_state( increment%valid_time(), &
                                                      depository )

  ! 1.2 Update the missing mixing ratio and moist_dynamics fields
  moisture_fields => self%modeldb%fields%get_field_collection("moisture_fields")
  call moisture_fields%get_field("ls_mr", mr_array)
  call moisture_fields%get_field("ls_moist_dyn", moist_dyn_array)
  call update_ls_moist_fields( mr_array%bundle, &
                               moist_dyn_array%bundle )

  ! 2. Step the linear model
  call step( self%modeldb )

  ! 3. Update the Atlas fields from the LFRic prognostic fields
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  call increment%set_from_model_prognostics( prognostic_fields )

  ! 4. Update the increment time
  call increment%update_time( self%time_step )

end subroutine model_stepTL

!> @brief    Finalise the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model finalise
subroutine model_finalTL(self, increment)

  implicit none

  class( jedi_linear_model_type ), target, intent(inout) :: self
  type( jedi_increment_type ),             intent(inout) :: increment

end subroutine model_finalTL

!> @todo: Will need to add model_initAD, model_stepAD and model_finalAD
!>        which will be used by forecastAD.

!> @brief    Finalize the jedi_linear_model_type
!>
subroutine jedi_linear_model_destructor(self)

  use jedi_lfric_linear_modeldb_driver_mod, only : finalise_modeldb

  implicit none

  type(jedi_linear_model_type), intent(inout) :: self

  call finalise_modeldb( self%modeldb )

end subroutine jedi_linear_model_destructor

!------------------------------------------------------------------------------
! OOPS defined forecastTL and forecastAD methods
!------------------------------------------------------------------------------

!> @brief    Run a Tangent Linear forecast using the model init, step and final
!>
!> @param [inout] increment       The Increment object to propagate
!> @param [in]    forecast_length The duration of the forecastTL
subroutine forecastTL( self, increment, forecast_length )

  implicit none

  class( jedi_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),     intent(inout) :: increment
  type( jedi_duration_type ),         intent(in) :: forecast_length

  ! Local
  type( jedi_datetime_type ) :: end_time

  ! End time
  end_time = increment%valid_time() + forecast_length

  ! Initialize the model
  call self%model_initTL( increment )
  ! Initialize the post processor and call first process
  ! call post_processor%pp_init( increment )
  ! call post_processor%process( increment )

  ! Loop until date_time_end
  do while ( end_time%is_ahead( increment%valid_time() ) )
    call self%model_stepTL( increment )
    ! call post_processor%process( increment )
  end do

  ! Finalize model and post processor
  ! call post_processor%pp_init( increment )
  call self%model_finalTL( increment )

end subroutine forecastTL

!> @todo: Need to add forecastAD that uses the model_initAD, model_stepAD and
!>        model_finalAD methods in a similar way to forecastTL

end module jedi_linear_model_mod
