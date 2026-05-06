!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a class that handles the full linear model
!>
!> @details This module includes a class that handles the full linear model
!>          initialisation, time stepping and finalisation for both the Tangent
!>          Linear (TL) and Adjoint (AD). These are required by the base
!>          interface class that uses them to define the forecastTL and
!>          forecastAD. The linear model forecasts (TL and AD) require a
!>          lineariastion-state (LS) trajectory. The LS trajectory is created
!>          by running the non-linear model and storing the result in an object
!>          of type: linear_state_trajectory_type. The set_trajectory method is
!>          included to provide the means to create and populate the LS fields.
!>
!>          An included forecast application (jedi_tlm_forecast_tl) uses the
!>          model forecastTL method to propagate the increment.
!>
module jedi_linear_model_mod

  use atl_si_timestep_alg_mod,       only : atl_si_timestep_type
  use base_wind_transform_mod,       only : base_wind_transform_type
  use constants_mod,                 only : str_def, i_def, l_def
  use driver_time_mod,               only : init_time, final_time
  use field_collection_mod,          only : field_collection_type
  use driver_modeldb_mod,            only : modeldb_type
  use incremental_wind_transform_mod, &
                                     only : incremental_wind_transform_type
  use jedi_base_linear_model_mod,    only : jedi_base_linear_model_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use jedi_lfric_moist_fields_mod,   only : update_ls_moist_fields,            &
                                            init_moist_fields,                 &
                                            adj_init_moist_fields,             &
                                            copy_moist_fields_from_prognostic, &
                                            copy_moist_fields_to_prognostic,   &
                                            zero_moist_fields
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_lfric_linear_fields_mod,  only : variable_names, &
                                            ls_variable_names, &
                                            create_linear_fields
  use jedi_lfric_wind_fields_mod,    only : create_scalar_winds, &
                                            setup_vector_wind
  use jedi_state_mod,                only : jedi_state_type
  use jedi_increment_mod,            only : jedi_increment_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use linear_state_trajectory_mod,   only : linear_state_trajectory_type
  use namelist_mod,                  only : namelist_type
  use normal_wind_transform_mod,     only : normal_wind_transform_type
  use zero_field_collection_mod,     only : zero_field_collection

  implicit none

  private

type, public, extends(jedi_base_linear_model_type) :: jedi_linear_model_type
  private

  !> The model time step duration
  type ( jedi_duration_type )          :: time_step

  !> Trajectory of linear states obtained by running the non-linear model
  type( linear_state_trajectory_type ) :: linear_state_trajectory

  !> Modeldb that stores the model fields to propagate
  !> @todo: Required public for checksum but need to move to atlas checksum
  !>        so make it private when that work is done.
  type(modeldb_type), public           :: modeldb

  !> Object encapsulating adjoint model timestep routines
  !> @todo: Will be incorporated into modeldb in #620
  type( atl_si_timestep_type )         :: atl_si_timestep

  !> Object encapsulating transformation between JEDI analysis wind variables and LFRic prognostic wind variables
  class(base_wind_transform_type), allocatable :: wind_transform

contains

  !> Model initialiser.
  procedure, public  :: initialise

  !> Methods
  procedure, public  :: set_trajectory

  procedure, public :: model_initTL
  procedure, public :: model_stepTL
  procedure, public :: model_finalTL

  procedure, public :: model_initAD
  procedure, public :: model_stepAD
  procedure, public :: model_finalAD

  !> Finalizer
  final              :: jedi_linear_model_destructor

end type jedi_linear_model_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_linear_model_type
!>
!> @param [in] jedi_geometry   A JEDI geometry object
!> @param [in] config_filename The name of the configuration file
subroutine initialise( self, jedi_geometry, config_filename )

  use jedi_lfric_linear_modeldb_driver_mod, only : initialise_modeldb
  use jedi_lfric_timestep_mod,              only : get_configuration_timestep

  implicit none

  class( jedi_linear_model_type ), target, intent(inout) :: self
  type( jedi_geometry_type ),                 intent(in) :: jedi_geometry
  character(len=*),                           intent(in) :: config_filename

  ! Local
  type(field_collection_type), pointer :: prognostic_fields
  type( namelist_type ),       pointer :: jedi_lfric_settings_config
  type( namelist_type ),       pointer :: jedi_linear_model_config
  character( str_def )                 :: forecast_length_str
  character( str_def )                 :: nl_time_step_str
  type( jedi_duration_type )           :: forecast_length
  type( jedi_duration_type )           :: nl_time_step
  logical(kind=l_def)                  :: incremental_wind_interpolation

  ! 1. Setup modeldb

  ! 1.1 Initialise the modeldb
  call initialise_modeldb( "linear modeldb", config_filename,           &
                            jedi_geometry%get_mpi_comm(), self%modeldb, &
                            self%atl_si_timestep )

  ! 1.2 Add scalar winds that link the Atlas fields. These are used to
  ! perform interpolation between horizontally cell-centred and
  ! edge based W2 winds.
  call create_scalar_winds( self%modeldb, jedi_geometry%get_mesh() )

  ! Set up extra fields for incremental wind interpolation, if required
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  jedi_linear_model_config => self%modeldb%configuration%get_namelist('jedi_linear_model')
  call jedi_linear_model_config%get_value( 'incremental_wind_interpolation', incremental_wind_interpolation )
  if (incremental_wind_interpolation) then
    allocate (incremental_wind_transform_type :: self%wind_transform)
  else
    allocate (normal_wind_transform_type :: self%wind_transform)
  end if
  call self%wind_transform%initialise(prognostic_fields)

  ! 2. Setup time
  self%time_step = get_configuration_timestep( self%modeldb%configuration )

  jedi_lfric_settings_config => self%modeldb%configuration%get_namelist('jedi_lfric_settings')
  call jedi_lfric_settings_config%get_value( 'forecast_length', forecast_length_str )
  call forecast_length%init(forecast_length_str)

  ! 3. Setup trajactory
  call jedi_linear_model_config%get_value( 'nl_time_step', nl_time_step_str )
  call nl_time_step%init(nl_time_step_str)

  call self%linear_state_trajectory%initialise( forecast_length, &
                                                nl_time_step )

end subroutine initialise

!> @brief    Set an instance of the trajectory
!>
!> @param [in] jedi_state The state to add to the trajectory
subroutine set_trajectory( self, jedi_state )

  implicit none

  class( jedi_linear_model_type ),   intent(inout) :: self
  type( jedi_state_type ),           intent(inout) :: jedi_state

  ! Local
  type( field_collection_type ) :: next_linear_state

  ! Create field collection that contains the linear state fields
  ! without "ls_" prepended.
  call create_linear_fields(jedi_state%geometry%get_mesh(), jedi_state%geometry%get_twod_mesh(), next_linear_state)

  ! Copy data from the input state into next_linear_state
  call jedi_state%get_to_field_collection( ls_variable_names, &
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

  ! Update the LFRic modeldb pertubation fields

  ! Update the prognostic fields: copy from Atlas
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  call increment%get_to_field_collection(variable_names, prognostic_fields)
  call self%wind_transform%scalar_to_vector(prognostic_fields)

  ! Update the missing mixing ratio and moist_dynamics fields. These fields are
  ! computed analytically as outlined in jedi_lfric_linear_fields_mod
  moisture_fields => self%modeldb%fields%get_field_collection("moisture_fields")
  call copy_moist_fields_from_prognostic( moisture_fields, prognostic_fields )
  call init_moist_fields( moisture_fields )

  ! Initialise clock and calendar
  call init_time( self%modeldb )

end subroutine model_initTL

!> @brief    Step the Tangent Linear model
!>
!> @param [inout] increment Increment object to be propagated
subroutine model_stepTL(self, increment)

  use jedi_lfric_linear_modeldb_driver_mod, only : step_tl

  implicit none

  class( jedi_linear_model_type ), target, intent(inout) :: self
  type( jedi_increment_type ),             intent(inout) :: increment

  ! Local
  type( field_collection_type ), pointer :: ls_fields
  type( field_collection_type),  pointer :: moisture_fields
  type( field_collection_type),  pointer :: prognostic_fields

  ! 1. Update the LFRic modeldb linear state fields

  ! 1.1 Copy from the trajectory into the model_data
  ls_fields => self%modeldb%fields%get_field_collection("ls_fields")
  call self%linear_state_trajectory%get_linear_state( increment%valid_time(), &
                                                      ls_fields )

  ! 1.2 Update the missing mixing ratio and moist_dynamics fields
  moisture_fields => self%modeldb%fields%get_field_collection("moisture_fields")
  call update_ls_moist_fields( ls_fields, moisture_fields )

  ! Copy pre-step vector winds if doing incremental wind interpolation
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  call self%wind_transform%process(prognostic_fields)

  ! 2. Step the linear model
  call step_tl( self%modeldb )

  ! 3. Update the Atlas fields from the LFRic prognostic fields
  call self%wind_transform%vector_to_scalar(prognostic_fields)
  call copy_moist_fields_to_prognostic( moisture_fields, prognostic_fields )
  call increment%set_from_field_collection( variable_names, prognostic_fields )

  ! 4. Update the increment time
  call increment%update_time( self%time_step )

end subroutine model_stepTL

!> @brief    Finalise the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model finalise
subroutine model_finalTL(self, increment)

  implicit none

  class( jedi_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),     intent(inout) :: increment

  ! Finalise clock and calendar
  call final_time( self%modeldb )

end subroutine model_finalTL

!> @brief    Initialise the Adjoint of the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model initialise
subroutine model_initAD(self, increment)

  implicit none

  class( jedi_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),     intent(inout) :: increment

  !>@todo: in ticket #267
  !>       Replace clock with reversible clock

  ! Local
  type( field_collection_type), pointer :: prognostic_fields
  type( field_collection_type), pointer :: moisture_fields

  ! Initialise clock and calendar
  call init_time( self%modeldb )

  ! Update the prognostic fields: zero LFRic fields
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  call zero_field_collection(prognostic_fields)
  moisture_fields => self%modeldb%fields%get_field_collection("moisture_fields")
  call zero_moist_fields(moisture_fields)
  call self%wind_transform%initialise_for_adjoint()

end subroutine model_initAD

!> @brief    Step the Adjoint of the Tangent Linear model
!>
!> @param [inout] increment Increment object to be propagated
subroutine model_stepAD(self, increment)

  use jedi_lfric_linear_modeldb_driver_mod, only : step_ad

  implicit none

  class( jedi_linear_model_type ), target, intent(inout) :: self
  type( jedi_increment_type ),             intent(inout) :: increment

  ! Local
  type( field_collection_type ), pointer :: ls_fields
  type( field_collection_type),  pointer :: moisture_fields
  type( field_collection_type),  pointer :: prognostic_fields

  ! Adjoint model integrates backwards in time, hence negative timestep
  call increment%update_time( self%time_step * (-1_i_def) )

  ls_fields => self%modeldb%fields%get_field_collection("ls_fields")
  call self%linear_state_trajectory%get_linear_state( increment%valid_time(), &
                                                      ls_fields )

  moisture_fields => self%modeldb%fields%get_field_collection("moisture_fields")
  call update_ls_moist_fields( ls_fields, moisture_fields )
  !>@todo: In ticket #775
  !>       This is always a repeated computation
  !>       Could be removed by storing the moist dynamics fields after the
  !>       equivalent call in the TL step

  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  call increment%set_from_field_collection_ad( variable_names, prognostic_fields )
  call copy_moist_fields_from_prognostic( moisture_fields, prognostic_fields )
  call self%wind_transform%adj_vector_to_scalar(prognostic_fields)

  call step_ad( self%modeldb, self%atl_si_timestep )

  ! Do adjoint of copy of vector winds if doing incremental wind interpolation
  call self%wind_transform%adj_process(prognostic_fields)

end subroutine model_stepAD

!> @brief    Finalise the Adjoint of the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model finalise
subroutine model_finalAD(self, increment)

  implicit none

  class( jedi_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),     intent(inout) :: increment

  !>@todo: in ticket #267
  !>       Replace clock with reversible clock

  ! Local
  type( field_collection_type),  pointer :: moisture_fields
  type( field_collection_type),  pointer :: prognostic_fields

  moisture_fields => self%modeldb%fields%get_field_collection("moisture_fields")
  call adj_init_moist_fields( moisture_fields )
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  call copy_moist_fields_to_prognostic( moisture_fields, prognostic_fields )
  call self%wind_transform%adj_scalar_to_vector( prognostic_fields )
  call increment%get_to_field_collection_ad(variable_names, prognostic_fields)

  ! Finalise clock and calendar
  call final_time( self%modeldb )

end subroutine model_finalAD

!> @brief    Finalize the jedi_linear_model_type
!>
subroutine jedi_linear_model_destructor(self)

  use jedi_lfric_linear_modeldb_driver_mod, only : finalise_modeldb

  implicit none

  type(jedi_linear_model_type), intent(inout) :: self

  call finalise_modeldb( self%modeldb, self%atl_si_timestep )

end subroutine jedi_linear_model_destructor

end module jedi_linear_model_mod
