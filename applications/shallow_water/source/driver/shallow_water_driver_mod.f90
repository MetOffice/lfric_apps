!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Drives the execution of the shallow_water miniapp.
!!
module shallow_water_driver_mod

  use base_mesh_config_mod,          only: prime_mesh_name
  use calendar_mod,                  only: calendar_type
  use constants_mod,                 only: i_def, imdi, r_def
  use driver_modeldb_mod,            only: modeldb_type
  use field_mod,                     only: field_type
  use field_collection_mod,          only: field_collection_type
  use io_config_mod,                 only: write_diag,           &
                                           diagnostic_frequency, &
                                           nodal_output_on_w3
  use log_mod,                       only: log_event,         &
                                           log_scratch_space, &
                                           log_level_error,   &
                                           LOG_LEVEL_INFO
  use mesh_collection_mod,           only: mesh_collection
  use mesh_mod,                      only: mesh_type
  use runtime_constants_mod,         only: final_runtime_constants
  use shallow_water_diagnostics_mod, only: shallow_water_diagnostics
  use shallow_water_model_mod,       only: initialise_infrastructure, &
                                           initialise_model,          &
                                           finalise_infrastructure,   &
                                           finalise_model
  use shallow_water_init_fields_mod, only: create_model_data,     &
                                           initialise_model_data, &
                                           output_model_data,     &
                                           finalise_model_data
  use shallow_water_step_mod,        only: shallow_water_step

  implicit none

  private

  public :: initialise, &
            step,       &
            finalise

  ! Mesh
  type(mesh_type), pointer :: mesh => null()

contains

  !=============================================================================
  !> @brief   Sets up required state in preparation for run.
  !> @details Initialises the infrastructure and the fields stored in
  !!          model_data, then sets the initial conditions for the run.
  !> @param [in,out] modeldb      The structure that holds model state
  !> @param [in]     program_name An identifier given to the model begin run
  !> @param [in]     calendar     Interprets date strings.
  !>
  subroutine initialise( modeldb, program_name, calendar )

    implicit none

    type(modeldb_type),   intent(inout) :: modeldb
    character(*),         intent(in)    :: program_name
    class(calendar_type), intent(in)    :: calendar

    ! Initialise infrastructure (from shallow_water_model_mod.F90) and setup constants
    call initialise_infrastructure( program_name, modeldb )

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    call log_event( 'Creating model data ...', LOG_LEVEL_INFO )
    ! Instantiate the fields stored in model_data
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    call create_model_data( modeldb, mesh )

    call log_event( 'Initialising model data ...', LOG_LEVEL_INFO )
    ! Initialise the fields stored in the model_data
    call initialise_model_data( modeldb, mesh )

    ! Initial output: we only want these once at the beginning of a run
    if (modeldb%clock%is_initialisation() .and. write_diag) then
      call log_event( 'Output of initial diagnostics ...', LOG_LEVEL_INFO )
      ! Calculation and output of diagnostics
      call shallow_water_diagnostics( mesh,               &
                                      modeldb,            &
                                      nodal_output_on_w3, &
                                      1_i_def )
    end if

    ! Model configuration initialisation
    call initialise_model( mesh, modeldb )

  end subroutine initialise

  !=============================================================================
  !> @brief Performs a time step for the shallow water miniapp.
  !> @param [in,out] modeldb The structure that holds model state
  !!
  subroutine step( modeldb )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call shallow_water_step( modeldb  )

    ! Use diagnostic output frequency to determine whether to write
    ! diagnostics on this timestep

    if ( ( mod(modeldb%clock%get_step(), diagnostic_frequency) == 0 ) &
         .and. ( write_diag ) ) then

      call shallow_water_diagnostics(mesh,               &
                                     modeldb,            &
                                     nodal_output_on_w3, &
                                     0_i_def)

    end if

  end subroutine step

  !=============================================================================
  !> @brief Tidies up after a run.
  !> @param [in,out] modeldb The structure that holds model state
  !!
  subroutine finalise( modeldb, program_name )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb
    character(*), intent(in) :: program_name

    ! Output the fields stored in the model_data (checkpoint and dump)
    call output_model_data( modeldb )

    ! Model configuration finalisation
    call finalise_model( modeldb, &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( modeldb )

    !-------------------------------------------------------------------------
    ! Finalise constants
    !-------------------------------------------------------------------------

    call final_runtime_constants()

    ! Finalise infrastructure and constants
    call finalise_infrastructure()

  end subroutine finalise

end module shallow_water_driver_mod
